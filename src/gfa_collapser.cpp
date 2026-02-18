#include "../include/gfa_collapser.hpp"
#include "../include/gfa_name.hpp"
#include "../include/aligner.hpp"
#include "../include/minimizer_dna.hpp"

#include <algorithm>
#include <limits>
#include <cstdint>
#include <unordered_set>


/* -------------------------------------------- collapse unitigs -------------------------------------------- */
bool GfaCollapser::try_build_chain_(
    uint32_t start,
    std::vector<uint32_t>& chain,
    std::vector<size_t>& edge_idx_seq
) {
    chain.clear();
    edge_idx_seq.clear();

    // -------- Bacnktrack to the true "chain head" --------
    // Rule: As long as "indegree == 1 and predecessor's outdegree == 1", keep going back
    uint32_t head = start;
    for (;;) {
        auto in_idxs = getArcsIdxToVertex(head, true);
        if (in_idxs.size() != 1) break;
        const GfaArc& in_e = arcs_[in_idxs[0]];
        uint32_t pred = in_e.get_source_vertex_id();
        if (getArcsIdxFromVertex(pred, true).size() != 1) break;
        head = pred;
    }

    // -------- Extend --------
    chain.push_back(head);
    uint32_t cur = head;

    while (true) {
        // 1. Current vertex must have exactly one outgoing, non-deleted edge.
        auto out_idxs = getArcsIdxFromVertex(cur, true);
        if (out_idxs.size() != 1) break;

        const GfaArc& arc = arcs_[out_idxs[0]];
        uint32_t w = arc.get_target_vertex_id();

        // 2. The next vertex must be linear on the incoming side: indegree == 1.
        if (getArcsIdxToVertex(w, true).size() != 1) break;

        // Record this step
        edge_idx_seq.push_back(out_idxs[0]);
        chain.push_back(w);

        // 3. The next vertex must have exactly one outgoing, non-deleted edge.
        if (getArcsIdxFromVertex(w, true).size() != 1) break;

        cur = w;
    }

    return chain.size() >= 2;
}

void GfaCollapser::collapse_unitigs() {
    log_stream() << "Collapsing linear unitigs (head outdeg=1, tail indeg=1, internal indeg=outdeg=1) ...\n";

    // Mark segments that have been included in a chain to avoid duplicate collection in the same round
    std::vector<uint8_t> used(nodes_.size(), 0);
    size_t nodes_merged_num = 0;

    gfaName namer;
    ProgressTracker prog(nodes_.size());

    // A chain's pending changes
    struct Pending {
        std::string new_name;
        std::string merged_seq;
        std::vector<uint32_t> chain_vtx;
        std::vector<size_t>   edge_idx_seq;
        std::vector<uint32_t> in_src_vs;
        std::vector<size_t>   in_del_idxs;
        std::vector<uint32_t> out_dst_ws;
        std::vector<size_t>   out_del_idxs;
        uint32_t head_forward = UINT32_MAX;
        uint32_t tail_forward = UINT32_MAX;
    };

    std::vector<Pending> pendings;
    pendings.reserve(1024);

    auto outdeg = [&](uint32_t v){ return getArcsIdxFromVertex(v).size(); };
    auto indeg  = [&](uint32_t v){ return getArcsIdxToVertex(v).size(); };

    // 2. Traverse all segments and collect candidate chains by "true chain head"
    for (uint32_t seg_id = 0; seg_id < nodes_.size(); ++seg_id) {
        prog.hit();
        if (used[seg_id] || nodes_[seg_id].deleted) continue;

        uint32_t vf = (seg_id << 1) | 0;
        uint32_t vr = vf ^ 1;

        std::vector<uint32_t> candidates;
        if (outdeg(vf) == 1) candidates.push_back(vf);
        if (outdeg(vr) == 1) candidates.push_back(vr);
        if (candidates.empty()) continue;

        std::vector<uint32_t> chain;
        std::vector<size_t>   edge_idx_seq;
        bool built = false;
        for (uint32_t cand : candidates) {
            if (try_build_chain_(cand, chain, edge_idx_seq)) { built = true; break; }
        }
        if (!built) continue;

        bool conflict = false;
        for (uint32_t vtx : chain) {
            uint32_t sid = NodeHandle::get_segment_id(vtx);
            if (sid < used.size() && (used[sid] || nodes_[sid].deleted)) { conflict = true; break; }
        }
        if (conflict) continue;

        std::string merged = get_oriented_sequence(chain.front());
        for (size_t i = 0; i < edge_idx_seq.size(); ++i) {
            const GfaArc& a = arcs_[edge_idx_seq[i]];
            uint32_t w = a.get_target_vertex_id();
            std::string wseq = get_oriented_sequence(w);
            uint32_t ow = (a.ow == INT32_MAX ? 0u : (uint32_t)std::max(0, a.ow));
            if (ow < wseq.size()) merged.append(wseq.begin() + ow, wseq.end());
        }

        std::vector<std::string> oriented_parts;
        oriented_parts.reserve(chain.size());
        for (uint32_t vtx : chain) {
            uint32_t sid = NodeHandle::get_segment_id(vtx);
            bool  is_rev = NodeHandle::get_is_reverse(vtx);
            std::string base = nodes_[sid].name;
            if (is_rev) base = namer.force_name_dir(base, is_rev);
            oriented_parts.push_back(base);
        }
        std::string new_name = namer.merge_chain_names(oriented_parts);

        Pending P;
        P.new_name   = std::move(new_name);
        P.merged_seq = std::move(merged);
        P.chain_vtx  = chain;
        P.edge_idx_seq = edge_idx_seq;

        // Record head/tail forward vertices
        P.head_forward = chain.front();
        P.tail_forward = chain.back();

        for (size_t ei : getArcsIdxToVertex(P.head_forward)) {
            const auto &e = arcs_[ei];
            if (e.get_del()) continue;
            P.in_src_vs.push_back(e.get_source_vertex_id());
            P.in_del_idxs.push_back(ei);
        }
        for (size_t ei : getArcsIdxFromVertex(P.tail_forward)) {
            const auto &e = arcs_[ei];
            if (e.get_del()) continue;
            P.out_dst_ws.push_back(e.get_target_vertex_id());
            P.out_del_idxs.push_back(ei);
        }

        for (uint32_t vtx : P.chain_vtx) {
            uint32_t sid = NodeHandle::get_segment_id(vtx);
            if (sid < used.size()) used[sid] = 1;
        }
        pendings.emplace_back(std::move(P));
        nodes_merged_num += chain.size();
    }

    // First create all new segments and establish the mapping from "old endpoints -> new vertices"
    std::vector<uint32_t> new_vertices; new_vertices.reserve(pendings.size());
    // Mapping: the "head/tail (directed vertex)" of a chain -> the new_vertex after merging
    std::unordered_map<uint32_t,uint32_t> bound2new; bound2new.reserve(pendings.size()*2);

    // Add all new segments first
    for (auto &P : pendings) {
        uint32_t new_seg_id = add_segment(P.new_name, P.merged_seq);
        uint32_t new_vertex = (new_seg_id << 1) | 0;
        new_vertices.push_back(new_vertex);
    }
    // Build head/tail to new_vertex maps
    for (size_t i = 0; i < pendings.size(); ++i) {
        bound2new[pendings[i].head_forward] = new_vertices[i];
        bound2new[pendings[i].tail_forward] = new_vertices[i];
        auto head_rev = pendings[i].head_forward ^ 1;
        auto tail_rev = pendings[i].tail_forward ^ 1;
        bound2new[head_rev] = new_vertices[i] ^ 1;
        bound2new[tail_rev] = new_vertices[i] ^ 1;
    }

    auto enc = [](uint32_t s, uint32_t t) -> uint64_t {
        return (uint64_t(s) << 32) | t;
    };
    std::unordered_set<uint64_t> edge_seen;
    edge_seen.reserve(pendings.size() * 8);

    // Reconnect
    for (size_t i = 0; i < pendings.size(); ++i) {
        auto &P = pendings[i];
        uint32_t newv = new_vertices[i];

        // v -> newv
        for (uint32_t v : P.in_src_vs) {
            uint32_t src = v;
            auto it = bound2new.find(v);
            if (it != bound2new.end()) src = it->second;

            if (src == newv) continue;  // avoid self-loop
            uint64_t k = enc(src, newv);
            if (edge_seen.insert(k).second) {  // deduplicate
                add_arc(src, newv, 0, 0, -1, /*comp=*/false);
            }
        }

        // newv -> w
        for (uint32_t w : P.out_dst_ws) {
            uint32_t dst = w;
            auto it = bound2new.find(w);
            if (it != bound2new.end()) dst = it->second;

            if (dst == newv) continue;  // avoid self-loop
            uint64_t k = enc(newv, dst);
            if (edge_seen.insert(k).second) {  // deduplicate
                add_arc(newv, dst, 0, 0, -1, /*comp=*/false);
            }
        }

        // Mark old edges for deletion
        for (size_t ei : P.in_del_idxs)   if (ei < arcs_.size()) arcs_[ei].set_del(true);
        for (size_t ei : P.out_del_idxs)  if (ei < arcs_.size()) arcs_[ei].set_del(true);
        for (size_t ei : P.edge_idx_seq)  if (ei < arcs_.size()) arcs_[ei].set_del(true);
    }

    for (auto &P : pendings) {
        if (P.chain_vtx.size() >= 2) {
            for (uint32_t vtx : P.chain_vtx) {
                uint32_t sid = NodeHandle::get_segment_id(vtx);
                delete_segment(sid);
            }
        }
    }

    rebuild_after_edits();
    log_stream() << "   - Total unitigs collapsed (vertices involved): " << nodes_merged_num << "\n";
}

/* -------------------------------------------- collapse bubbles -------------------------------------------- */
std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> GfaCollapser::split_paths_(
    const std::vector<uint32_t> &path1,
    const std::vector<uint32_t> &path2
) {
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> result;
    if (path1.empty() || path2.empty()) return result;

    std::unordered_map<uint32_t, size_t> pos2;
    pos2.reserve(path2.size() * 2);
    for (size_t i = 0; i < path2.size(); ++i) {
        pos2[path2[i]] = i;
    }

    std::vector<std::pair<size_t, size_t>> anchors;
    anchors.reserve(std::min(path1.size(), path2.size()));
    std::ptrdiff_t last_j = -1;
    for (size_t i = 0; i < path1.size(); ++i) {
        auto it = pos2.find(path1[i]);
        if (it == pos2.end()) continue;
        if (static_cast<std::ptrdiff_t>(it->second) <= last_j) continue;
        anchors.emplace_back(i, it->second);
        last_j = static_cast<std::ptrdiff_t>(it->second);
    }

    if (anchors.empty()) {
        result.emplace_back(path1, path2);
        return result;
    }

    size_t prev1 = 0, prev2 = 0;
    for (const auto &ac : anchors) {
        size_t i1 = ac.first;
        size_t i2 = ac.second;
        if (i1 > prev1 && i2 > prev2) {
            result.emplace_back(
                std::vector<uint32_t>(path1.begin() + prev1, path1.begin() + i1),
                std::vector<uint32_t>(path2.begin() + prev2, path2.begin() + i2)
            );
        }
        prev1 = i1 + 1;
        prev2 = i2 + 1;
    }

    if (prev1 < path1.size() && prev2 < path2.size()) {
        result.emplace_back(
            std::vector<uint32_t>(path1.begin() + prev1, path1.end()),
            std::vector<uint32_t>(path2.begin() + prev2, path2.end())
        );
    }

    return result;
}

std::string GfaCollapser::vertex_name_(uint32_t vtx) const {
    uint32_t sid = NodeHandle::get_segment_id(vtx);
    return getNodeName(sid);
}

void GfaCollapser::build_concat_seq_(const std::vector<uint32_t>& path, std::string &seq, std::vector<uint32_t> &offsets)
{
    offsets.clear();
    offsets.reserve(path.size() + 1);
    offsets.push_back(0u);

    seq.clear();
    for (uint32_t vtx : path) {
        std::string s = get_oriented_sequence(vtx);
        seq.append(s);
        offsets.push_back(static_cast<uint32_t>(seq.size()));
    }
}

std::pair<size_t,uint32_t> GfaCollapser::find_node_idx_(uint32_t pos, const std::vector<uint32_t> &offs) const
{
    size_t lo = 0, hi = offs.size() - 1;
    while (lo + 1 < hi) {
        size_t mid = (lo + hi) >> 1;
        if (offs[mid] <= pos) lo = mid;
        else hi = mid;
    }
    size_t idx   = lo;
    uint32_t in  = pos - offs[idx];
    return { idx, in };
}

std::vector<GfaCollapser::BubbleAlignment> GfaCollapser::split_cigar_(
    const std::vector<uint32_t> &pathA,
    const std::vector<uint32_t> &pathB,
    const std::vector<uint32_t> &offsA,
    const std::vector<uint32_t> &offsB,
    const BubbleAlignment &big
) {
    std::vector<BubbleAlignment> outs;
    if (pathA.empty() || pathB.empty()) return outs;
    if (big.ops.empty()) return outs;

    const size_t nA = pathA.size();
    const size_t nB = pathB.size();
    const uint32_t LA = offsA.back();
    const uint32_t LB = offsB.back();
    if (LA == 0 || LB == 0) return outs;

    uint32_t posA = big.beg_a;
    uint32_t posB = big.beg_b;

    auto tmpA = find_node_idx_(posA, offsA);
    auto tmpB = find_node_idx_(posB, offsB);
    size_t   idxA = tmpA.first;
    size_t   idxB = tmpB.first;
    uint32_t offA = tmpA.second;
    uint32_t offB = tmpB.second;

    if (idxA >= nA || idxB >= nB) return outs;

    std::vector<CIGAR::COp> cur_ops;
    cur_ops.reserve(big.ops.size());

    bool     pair_started  = false;
    uint32_t pair_begA_abs = posA;
    uint32_t pair_begB_abs = posB;

    auto flush_pair = [&]() {
        if (!pair_started || cur_ops.empty()) return;
        if (idxA >= nA || idxB >= nB) {
            cur_ops.clear();
            pair_started = false;
            return;
        }

        uint32_t vtxA = pathA[idxA];
        uint32_t vtxB = pathB[idxB];

        uint32_t node_startA = offsA[idxA];
        uint32_t node_endA   = offsA[idxA + 1];
        uint32_t node_startB = offsB[idxB];
        uint32_t node_endB   = offsB[idxB + 1];
        uint32_t lenA_node   = node_endA - node_startA;
        uint32_t lenB_node   = node_endB - node_startB;

        // local (oriented) coords in current node
        uint32_t local_beg_a = (pair_begA_abs > node_startA) ? (pair_begA_abs - node_startA) : 0u;
        uint32_t local_end_a = (posA > node_startA) ? (posA - node_startA) : 0u;
        uint32_t local_beg_b = (pair_begB_abs > node_startB) ? (pair_begB_abs - node_startB) : 0u;
        uint32_t local_end_b = (posB > node_startB) ? (posB - node_startB) : 0u;

        if (local_beg_a > lenA_node) {
            error_stream() << vertex_name_(vtxA) << " vs " << vertex_name_(vtxB) << ": local_beg_a (" << local_beg_a << ") > lenA_node (" << lenA_node << ")\n";
            std::exit(1);
        }
        if (local_end_a > lenA_node) {
            error_stream() << vertex_name_(vtxA) << " vs " << vertex_name_(vtxB) << ": local_end_a (" << local_end_a << ") > lenA_node (" << lenA_node << ")\n";
            std::exit(1);
        }
        if (local_beg_b > lenB_node) {
            error_stream() << vertex_name_(vtxA) << " vs " << vertex_name_(vtxB) << ": local_beg_b (" << local_beg_b << ") > lenB_node (" << lenB_node << ")\n";
            std::exit(1);
        }
        if (local_end_b > lenB_node) {
            error_stream() << vertex_name_(vtxA) << " vs " << vertex_name_(vtxB) << ": local_end_b (" << local_end_b << ") > lenB_node (" << lenB_node << ")\n";
            std::exit(1);
        }

        // convert to forward coordinates if vtx is reverse
        bool revA = NodeHandle::get_is_reverse(vtxA);
        bool revB = NodeHandle::get_is_reverse(vtxB);

        uint32_t beg_a, end_a, beg_b, end_b;
        if (revA) {
            beg_a = lenA_node - local_end_a;
            end_a = lenA_node - local_beg_a;
        } else {
            beg_a = local_beg_a;
            end_a = local_end_a;
        }

        if (revB) {
            beg_b = lenB_node - local_end_b;
            end_b = lenB_node - local_beg_b;
        } else {
            beg_b = local_beg_b;
            end_b = local_end_b;
        }

        std::string nameA = vertex_name_(vtxA);
        std::string nameB = vertex_name_(vtxB);

        outs.emplace_back(
            std::move(nameA),
            std::move(nameB),
            beg_a, end_a,
            beg_b, end_b,
            vtxA, vtxB,
            std::move(cur_ops)
        );
        cur_ops.clear();
        pair_started  = false;
        pair_begA_abs = posA;
        pair_begB_abs = posB;
    };

    for (const auto &cop : big.ops) {
        char     op  = cop.op;
        uint32_t len = cop.len;
        if (len == 0) continue;

        bool consumeA = (op == 'M' || op == '=' || op == 'X' || op == 'D');
        bool consumeB = (op == 'M' || op == '=' || op == 'X' || op == 'I');

        while (len > 0 && idxA < nA && idxB < nB) {
            uint32_t remA = consumeA
                          ? (offsA[idxA + 1] - posA)
                          : std::numeric_limits<uint32_t>::max();
            uint32_t remB = consumeB
                          ? (offsB[idxB + 1] - posB)
                          : std::numeric_limits<uint32_t>::max();

            uint32_t step = len;
            if (consumeA && step > remA) step = remA;
            if (consumeB && step > remB) step = remB;

            if (step == 0) {
                flush_pair();

                bool moved = false;
                if (consumeA && remA == 0 && idxA + 1 < nA) {
                    ++idxA;
                    posA = offsA[idxA];
                    offA = 0;
                    moved = true;
                }
                if (consumeB && remB == 0 && idxB + 1 < nB) {
                    ++idxB;
                    posB = offsB[idxB];
                    offB = 0;
                    moved = true;
                }
                if (!moved) {
                    len = 0;
                    break;
                }
                continue;
            }

            if (!pair_started) {
                pair_started  = true;
                pair_begA_abs = posA;
                pair_begB_abs = posB;
            }

            if (!cur_ops.empty() && cur_ops.back().op == op) {
                cur_ops.back().len += step;
            } else {
                cur_ops.emplace_back(step, op);
            }

            if (consumeA) { posA += step; offA += step; }
            if (consumeB) { posB += step; offB += step; }

            len -= step;

            bool at_endA = consumeA && (posA == offsA[idxA + 1]);
            bool at_endB = consumeB && (posB == offsB[idxB + 1]);

            if (at_endA || at_endB) {
                flush_pair();

                if (at_endA && idxA + 1 < nA) {
                    ++idxA;
                    posA = offsA[idxA];
                    offA = 0;
                }
                if (at_endB && idxB + 1 < nB) {
                    ++idxB;
                    posB = offsB[idxB];
                    offB = 0;
                }
            }
        }
    }

    flush_pair();
    return outs;
}

std::vector<GfaCollapser::BubbleAlignment> GfaCollapser::align_subpaths_(
    const std::vector<uint32_t> &subA, 
    const std::vector<uint32_t> &subB
) {
    std::vector<BubbleAlignment> outs;
    if (subA.empty() || subB.empty()) return outs;

    static thread_local minimizerdna::Options mzopt{
        /*k=*/17, /*w=*/30, /*canonical=*/false,
        /*seed=*/0x8a5cd789635d2dffULL
    };
    static thread_local minimizerdna::MinimizerBuilder mzb(mzopt);

    // ---------- Single node vs single node ----------
    if (subA.size() == 1 && subB.size() == 1) {
        const uint32_t vtx_a = subA[0];
        const uint32_t vtx_b = subB[0];

        const uint32_t sid_a = NodeHandle::get_segment_id(vtx_a);
        const uint32_t sid_b = NodeHandle::get_segment_id(vtx_b);
        if (sid_a == sid_b) return outs;

        std::string seq_a = get_oriented_sequence(vtx_a);
        std::string seq_b = get_oriented_sequence(vtx_b);
        const uint32_t vb = 0u, ve = static_cast<uint32_t>(seq_a.size());
        const uint32_t wb = 0u, we = static_cast<uint32_t>(seq_b.size());

        if (ve <= 2 || we <= 2) return outs;

        // minimizer filtering
        double score = mzb.max_containment(seq_a, seq_b);
        if (score < MIN_JACCARD_FOR_ALIGN_) return outs;

        const std::string v_name = getNodeName(sid_a);
        const std::string w_name = getNodeName(sid_b);

        std::vector<BubbleAlignment> alns = align_and_pack_(
            vtx_a, vtx_b,
            v_name, w_name,
            seq_a, seq_b,
            vb, ve, wb, we
        );

        if (alns.empty()) return outs;

        outs.insert(outs.end(), std::make_move_iterator(alns.begin()), std::make_move_iterator(alns.end()));

        return outs;
    }

    // ---------- Multi-node vs multi-node ----------
    std::string seq_a, seq_b;
    std::vector<uint32_t> offsA, offsB;
    build_concat_seq_(subA, seq_a, offsA);
    build_concat_seq_(subB, seq_b, offsB);

    if (seq_a.size() <= chainOpts_.k || seq_b.size() <= chainOpts_.k) return outs;

    // minimizer filtering
    double score = mzb.max_containment(seq_a, seq_b);
    if (score < MIN_JACCARD_FOR_ALIGN_) return outs;

    const uint32_t vtx_a0 = (NodeHandle::get_segment_id(subA.front()) << 1) | 0u;
    const uint32_t vtx_b0 = (NodeHandle::get_segment_id(subB.front()) << 1) | 0u;

    std::string name_a = ""; std::string name_b = "";

    for (auto vtx : subA) {
        if (!name_a.empty()) name_a += ";";
        uint32_t sid = NodeHandle::get_segment_id(vtx);
        name_a += getNodeName(sid);
    }

    for (auto vtx : subB) {
        if (!name_b.empty()) name_b += ";";
        uint32_t sid = NodeHandle::get_segment_id(vtx);
        name_b += getNodeName(sid);
    }

    std::vector<BubbleAlignment> alns = align_and_pack_(
        vtx_a0, vtx_b0,
        name_a, name_b,
        seq_a, seq_b,
        0u, static_cast<uint32_t>(seq_a.size()),
        0u, static_cast<uint32_t>(seq_b.size())
    );
    if (alns.empty()) return outs;

    for (const auto& aln : alns) {
        std::vector<BubbleAlignment> out = split_cigar_(subA, subB, offsA, offsB, aln);

        for (auto &aln_tmp : out) {
            debug_stream() << "   - " << aln_tmp.name_a << " vs " << aln_tmp.name_b << ": "
                << aln_tmp.beg_a << "-" << aln_tmp.end_a << " | "
                << aln_tmp.beg_b << "-" << aln_tmp.end_b << " | "
                << CIGAR::pack(aln_tmp.ops) << "\n";
        }
        outs.insert(outs.end(), std::make_move_iterator(out.begin()), std::make_move_iterator(out.end()));
    }
    
    return outs;
}

void GfaCollapser::bubbles_align_(const std::vector<GfaBubble::Bubble>& bubbles, const std::vector<GfaBubble::ForkGroup>& forks, const GfaBubble::GfaBubbleFinder& bubble_finder) {
    log_stream() << "Aligning bubbles ...\n";
    bubble_aligns_.clear();

    // ---------------- Thread pool ----------------
    ThreadPool pool(alignOpts_.threads);
    std::vector<std::future<std::vector<BubbleAlignment>>> futs;

    auto submit_subpath_pair =
        [&] (std::vector<uint32_t> subA, std::vector<uint32_t> subB)
        {
            if (subA.empty() || subB.empty()) return;
            futs.emplace_back(
                pool.submit([this, subA = std::move(subA), subB = std::move(subB)]() {
                    return align_subpaths_(subA, subB);
                })
            );
        };

    // Deduplication
    std::unordered_set<uint64_t> seens; 
    seens.reserve(1024);
    auto pair_key_seg = [] (uint32_t s1, uint32_t s2) -> uint64_t {
        if (s1 > s2) std::swap(s1, s2);
        return (static_cast<uint64_t>(s1) << 32) | static_cast<uint64_t>(s2);
    };
    auto should_submit_pair =
        [&] (const std::vector<uint32_t>& subA,
             const std::vector<uint32_t>& subB) -> bool
    {
        if (subA.empty() || subB.empty()) return false;

        uint64_t total_len = 0, new_len = 0;

        for (uint32_t a : subA) {
            for (uint32_t b : subB) {
                uint32_t s1 = NodeHandle::get_segment_id(a);
                uint32_t s2 = NodeHandle::get_segment_id(b);
                uint64_t k = pair_key_seg(s1, s2);
                auto [it, inserted] = seens.insert(k);
                total_len += getNodeLength(s1) + getNodeLength(s2);
                new_len += inserted ? (getNodeLength(s1) + getNodeLength(s2)) : 0;
            }
        }

        return ((double)new_len / (double)total_len) >= MIN_NEW_LEN_FRAC_;
    };

    // Submit all path pairs in bubbles
    for (const auto &bb : bubbles) {
        const auto &paths = bb.get_paths();
        if (paths.size() < 2) continue;

        std::vector<std::vector<uint32_t>> paths_filter = bubble_finder.pick_representative_paths(paths);

        for (size_t i = 0; i + 1 < paths_filter.size(); ++i) {
            for (size_t j = i + 1; j < paths_filter.size(); ++j) {
                const auto &pathA = paths_filter[i];
                const auto &pathB = paths_filter[j];
                if (pathA.empty() || pathB.empty()) continue;

                auto chunks = split_paths_(pathA, pathB);
                for (auto &p : chunks) {
                    auto &subA = p.first;
                    auto &subB = p.second;
                    if (subA.empty() || subB.empty()) continue;

                    // deduplication according to vertex sets
                    if (!should_submit_pair(subA, subB)) continue;

                    submit_subpath_pair(std::move(p.first), std::move(p.second));
                }
            }
        }
    }

    for (const auto& fork : forks) {
        const auto& branches = fork.get_branches();
        if (branches.size() < 2) continue;

        for (size_t j = 0; j + 1 < branches.size(); ++j) {
            for (size_t k = j + 1; k < branches.size(); ++k) {
                const uint32_t bj = branches[j];
                const uint32_t bk = branches[k];

                const uint32_t sid_a = NodeHandle::get_segment_id(bj);
                const uint32_t sid_b = NodeHandle::get_segment_id(bk);
                if (sid_a == sid_b) continue;

                // deduplication according to vertex sets
                if (!should_submit_pair({bj}, {bk})) continue;

                submit_subpath_pair({bj}, {bk});
            }
        }
    }

    // Collect results
    ProgressTracker prog(futs.size());
    for (auto &f : futs) {
        prog.hit();
        auto vec = f.get();
        for (auto &aln : vec) {
            if (!aln.ops.empty()) {
                bubble_aligns_.emplace_back(std::move(aln));
            }
        }
    }

    pool.stop();

    log_stream() << "   - Total alignments produced: " << bubble_aligns_.size() << "\n";
}

void GfaCollapser::expand_and_rewire_edges_(
    const SegReplace::Expander& ex
) {
    log_stream() << "Expanding segments and rewiring edges ...\n";
    
    auto query = [&](SegReplace::Seg s){ return ex.query(s); };

    // make sure each edge is only added once
    std::unordered_set<uint64_t> seen_edges;

    gfaName namer;  // for generating new segment names

    uint64_t arc_num = arcs_.size();
    
    // Progress tracker for merging unitigs
    ProgressTracker prog(arc_num);

    for (size_t ai = 0; ai < arc_num; ++ai) {
        prog.hit();  // Update progress

        const GfaArc& arc = arcs_[ai];
        if (arc.get_del() || arc.get_comp()) continue;

        uint32_t v_seg_id = arc.get_source_segment_id();
        uint32_t w_seg_id = arc.get_target_segment_id();
        bool     v_is_rev = arc.get_source_is_reverse();
        bool     w_is_rev = arc.get_target_is_reverse();
        std::string v_name = nodes_[v_seg_id].name, w_name = nodes_[w_seg_id].name;

        if (arc.ov > 0 || arc.ow > 0) {
            error_stream() << "Non-zero overlap edges should have been removed in the de-overlapping step!\n";
            std::exit(1);
        }

        auto& v_cuts = const_cast<std::vector<uint32_t>&>(cuts_[v_seg_id].v);
        auto& w_cuts = const_cast<std::vector<uint32_t>&>(cuts_[w_seg_id].v);

        // ============= Collect the expansion of source/target/overlap =============
        auto collect_non_overlap = [&](
            uint32_t seg_id, bool seg_rev,
            const std::vector<uint32_t>& cuts_vec
        ) -> std::vector<SegReplace::Expansion> {
            std::vector<SegReplace::Expansion> res;
            for (size_t i = 0; i + 1 < cuts_vec.size(); i++) {
                uint32_t beg = cuts_vec[i], end = cuts_vec[i + 1];
                SegReplace::Seg window_u128 = SegReplace::Interval::pack(seg_id, beg, end, seg_rev);
                const auto leafs_u128 = query(window_u128);
                res.push_back(leafs_u128);
            }
            return res;
        };

        // no-overlap
        std::vector<SegReplace::Expansion> source_expansions = collect_non_overlap(v_seg_id, v_is_rev, v_cuts);
        std::vector<SegReplace::Expansion> target_expansions = collect_non_overlap(w_seg_id, w_is_rev, w_cuts);

        if (v_is_rev) std::reverse(source_expansions.begin(), source_expansions.end());
        if (w_is_rev) std::reverse(target_expansions.begin(), target_expansions.end());

        SegReplace::Expansion node_expansion;
        for (const auto& seq : source_expansions) node_expansion.insert(node_expansion.end(), seq.begin(), seq.end());
        for (const auto& seq : target_expansions) node_expansion.insert(node_expansion.end(), seq.begin(), seq.end());

        emit_node_expansion_edges_(
            node_expansion,
            seen_edges,
            namer,
            v_name, w_name,
            v_is_rev, w_is_rev
        );
    }

    // Mark all non-complementary arcs for deletion
    for (auto& e : arcs_) {
        if (!e.get_comp()) e.set_del(true);
    }
    rebuild_after_edits();

    // Add all new edges
    for (const auto& k : seen_edges) {
        auto [v, w] = decode_edge_u64_(k);
        add_arc(v, w, 0, 0, -1, false);
    }
    rebuild_after_edits();

    return;
}

void GfaCollapser::collapse_bubbles(const std::vector<GfaBubble::Bubble>& bubbles, const std::vector<GfaBubble::ForkGroup>& forks, const std::string& prefix, const GfaBubble::GfaBubbleFinder& bubble_finder) {
    log_stream() << "Collapsing bubbles ...\n";

    // 0.1. Align bubbles to find cut points
    bubbles_align_(bubbles, forks, bubble_finder);

    // 0.2. Deduplicate alignments
    dedup_aligns_();

    // 1. Initialize cuts_: boundaries
    initialize_cuts_();

    // 2. Build cuts_ from alignments
    build_cuts_from_cigar_();

    // 3. Propagate cuts_ through alignments
    propagate_cuts_();
    if (DEBUG_ENABLED) print_cuts(cuts_);

    // 4. Build rulemap_
    build_rulemap_();
    SegReplace::Expander ex = build_SegReplace_(false);
    ex.save_map(prefix + ".collapse.map");  // save rulemap

    // 5. Apply rules: expand and rewire 0M edges
    expand_and_rewire_edges_(ex);

    // 6. Remove unused nodes
    remove_unused_nodes_();

    log_stream() << "Finished collapsing bubbles.\n";
}