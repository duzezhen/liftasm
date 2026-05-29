#include "../include/gfa_collapser.hpp"
#include "../include/gfa_name.hpp"
#include "../include/aligner.hpp"
#include "../include/minimizer_dna.hpp"

#include <algorithm>
#include <limits>
#include <cstdint>
#include <unordered_set>
#include <deque>
#include <chrono>


// ------------------------------------------------ GFA deoverlapper ------------------------------------------------
static inline void append_cop_(std::vector<CIGAR::COp>& out, uint32_t len, char op) {
    if (len == 0) return;
    if (!out.empty() && out.back().op == op) out.back().len += len;
    else out.emplace_back(len, op);
}

static std::vector<uint64_t> build_tandem_repeat_mask_(
    const std::string& s,
    uint32_t min_len,
    uint32_t max_period,
    uint32_t max_mismatch
) {
    const uint32_t n = static_cast<uint32_t>(s.size());
    std::vector<uint64_t> mask((n + 63) >> 6, 0);
    if (n < min_len) return mask;
    if (min_len == 0 || max_period == 0) return mask;

    auto set_range = [&](uint32_t beg, uint32_t end) {
        if (beg >= end) return;
        if (end > n) end = n;

        const uint32_t w0 = beg >> 6;
        const uint32_t w1 = (end - 1) >> 6;

        if (w0 == w1) {
            mask[w0] |= (~uint64_t(0) << (beg & 63)) & (~uint64_t(0) >> (63 - ((end - 1) & 63)));
            return;
        }

        mask[w0] |= (~uint64_t(0) << (beg & 63));
        for (uint32_t w = w0 + 1; w < w1; ++w) mask[w] = ~uint64_t(0);
        mask[w1] |= (~uint64_t(0) >> (63 - ((end - 1) & 63)));
    };

    const uint32_t max_p = std::min(max_period, n >> 1);

    for (uint32_t p = 1; p <= max_p; ++p) {
        uint32_t beg = 0;
        uint32_t mismatches = 0;

        for (uint32_t i = 0; i + p < n; ++i) {
            if (s[i] != s[i + p]) {
                ++mismatches;
            }

            while (mismatches > max_mismatch && beg <= i) {
                if (s[beg] != s[beg + p]) {
                    --mismatches;
                }
                ++beg;
            }

            const uint32_t end = i + p + 1;

            if (end > beg && end - beg >= min_len) {
                set_range(beg, end);
            }
        }
    }

    return mask;
}

static void mask_repetitive_match_ops_(
    std::vector<CIGAR::COp>& ops,
    const std::vector<uint64_t>& mask_a,
    const std::vector<uint64_t>& mask_b,
    uint32_t beg_a,
    uint32_t beg_b, 
    uint32_t min_len,
    uint32_t max_period
) {
    auto test_bit = [](const std::vector<uint64_t>& mask, uint32_t pos) -> bool {
        const uint32_t w = pos >> 6;
        return w < mask.size() && (mask[w] & (uint64_t(1) << (pos & 63)));
    };

    if (min_len == 0 || max_period == 0) return;

    std::vector<CIGAR::COp> out;
    out.reserve(ops.size() + 8);

    uint32_t pa = beg_a;
    uint32_t pb = beg_b;

    for (const auto& op : ops) {
        const bool ca = path_pair_split::consumes_a(op.op);
        const bool cb = path_pair_split::consumes_b(op.op);

        if ((op.op == '=' || op.op == 'M') && ca && cb) {
            uint32_t left = op.len;

            while (left > 0) {
                const bool masked = test_bit(mask_a, pa) || test_bit(mask_b, pb);

                uint32_t step = 1;
                while (step < left && (test_bit(mask_a, pa + step) || test_bit(mask_b, pb + step)) == masked) {
                    ++step;
                }

                append_cop_(out, step, masked ? 'X' : op.op);

                pa += step;
                pb += step;
                left -= step;
            }
        } else {
            append_cop_(out, op.len, op.op);
            if (ca) pa += op.len;
            if (cb) pb += op.len;
        }
    }

    ops.swap(out);
}


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

void GfaCollapser::merge_linear_chains() {
    log_stream() << "Merging linear chains (internal nodes: indeg=1, outdeg=1) ...\n";

    bool has_overlap = false;
    for (const auto& a : arcs_) {
        if (a.get_del()) continue;
        if ((a.ov != INT32_MAX && a.ov > 0) || (a.ow != INT32_MAX && a.ow > 0)) {
            has_overlap = true;
            break;
        }
    }

    if (has_overlap) {
        warning_stream() << "  ! Overlap graph detected; skip linear-chain merge.\n";
        return;
    }

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

        std::string merged = get_oriented_sequence(Vertex(chain.front()));
        for (size_t i = 0; i < edge_idx_seq.size(); ++i) {
            const GfaArc& a = arcs_[edge_idx_seq[i]];
            uint32_t w = a.get_target_vertex_id();
            std::string wseq = get_oriented_sequence(Vertex(w));
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
    log_stream() << "  - Total merged nodes: " << nodes_merged_num << "\n" << "\n";
}


/* -------------------------------------------- Collapse Homologous Sequences -------------------------------------------- */
std::vector<GfaCollapser::BubbleAlignment> GfaCollapser::align_subpaths_(
    const std::vector<uint32_t>& subA,
    const std::vector<uint32_t>& subB,
    const path_pair_split::Span& spanA,
    const path_pair_split::Span& spanB
) {
    // ------------------------------------------------ Structure ------------------------------------------------
    struct Piece {
        uint32_t vtx = 0;
        uint32_t cb = 0, ce = 0; // concat coords in sliced sequence
        uint32_t ob = 0, oe = 0; // oriented node coords
    };

    // ------------------------------------------------ Helper functions ------------------------------------------------
    auto len_of = [this](uint32_t v) {
        return getNodeLength(NodeHandle::get_segment_id(v));
    };
    auto rev_of = [](uint32_t v) {
        return NodeHandle::get_is_reverse(v);
    };

    auto build = [&](const std::vector<uint32_t>& path, const path_pair_split::Span& span, std::string& seq, std::vector<Piece>& pieces) {
        seq.clear();
        pieces.clear();
        if (path.empty()) return;

        auto overlap_of = [this](uint32_t v, uint32_t w) {
            return get_edge_ow(v, w);
        };

        path_pair_split::PathLayout<uint32_t, decltype(len_of), decltype(rev_of), decltype(overlap_of)> layout(path, len_of, rev_of, overlap_of);

        uint32_t cur = layout.abs(span.beg);
        const uint32_t end = layout.abs(span.end);
        uint32_t out_pos = 0;

        while (cur < end) {
            auto p = layout.pos(cur);
            if (p.node >= path.size()) break;

            uint32_t vtx = path[p.node];
            uint32_t sid = NodeHandle::get_segment_id(vtx);
            uint32_t L = getNodeLength(sid);

            uint32_t ob = p.off;
            uint32_t oe = std::min<uint32_t>(L, ob + (end - cur));

            if (ob < oe) {
                std::string s = get_oriented_sequence(Vertex(vtx));
                seq.append(s.data() + ob, oe - ob);
                pieces.push_back({vtx, out_pos, out_pos + (oe - ob), ob, oe});
                out_pos += oe - ob;
            }

            cur = layout.abs({p.node + 1, 0});
        }
    };

    auto split = [&](const std::vector<Piece>& A, const std::vector<Piece>& B, const BubbleAlignment& big) {
        std::vector<BubbleAlignment> out;
        if (A.empty() || B.empty() || big.ops.empty()) return out;

        size_t ia = 0, ib = 0;
        uint32_t pa = big.beg_a, pb = big.beg_b;

        size_t seg_ia = 0, seg_ib = 0;
        uint32_t ba = pa, bb = pb;
        bool started = false;
        std::vector<CIGAR::COp> ops;

        auto advance_if_needed = [&]() {
            while (ia + 1 < A.size() && pa >= A[ia].ce) ++ia;
            while (ib + 1 < B.size() && pb >= B[ib].ce) ++ib;
        };

        auto flush = [&]() {
            if (!started || ops.empty()) {
                ops.clear();
                started = false;
                return;
            }

            if (seg_ia >= A.size() || seg_ib >= B.size()) {
                ops.clear();
                started = false;
                return;
            }

            const Piece& x = A[seg_ia];
            const Piece& y = B[seg_ib];

            if (ba < x.cb || pa > x.ce || bb < y.cb || pb > y.ce) {
                ops.clear();
                started = false;
                return;
            }

            uint32_t sid_a = NodeHandle::get_segment_id(x.vtx);
            uint32_t sid_b = NodeHandle::get_segment_id(y.vtx);
            uint32_t len_a = getNodeLength(sid_a);
            uint32_t len_b = getNodeLength(sid_b);

            uint32_t oa0 = x.ob + (ba - x.cb);
            uint32_t oa1 = x.ob + (pa - x.cb);
            uint32_t ob0 = y.ob + (bb - y.cb);
            uint32_t ob1 = y.ob + (pb - y.cb);

            if (oa0 >= oa1 || ob0 >= ob1) {
                ops.clear();
                started = false;
                ba = pa;
                bb = pb;
                return;
            }

            uint32_t fa0, fa1, fb0, fb1;
            if (NodeHandle::get_is_reverse(x.vtx)) {
                fa0 = len_a - oa1;
                fa1 = len_a - oa0;
            } else {
                fa0 = oa0;
                fa1 = oa1;
            }

            if (NodeHandle::get_is_reverse(y.vtx)) {
                fb0 = len_b - ob1;
                fb1 = len_b - ob0;
            } else {
                fb0 = ob0;
                fb1 = ob1;
            }

            out.emplace_back(
                getNodeName(sid_a),
                getNodeName(sid_b),
                fa0, fa1,
                fb0, fb1,
                x.vtx, y.vtx,
                std::move(ops)
            );

            ops.clear();
            started = false;
            ba = pa;
            bb = pb;
        };

        for (const auto& cop : big.ops) {
            uint32_t left = cop.len;
            const bool ca = path_pair_split::consumes_a(cop.op);
            const bool cb = path_pair_split::consumes_b(cop.op);

            while (left && ia < A.size() && ib < B.size()) {
                if (!started) advance_if_needed();
                if (ia >= A.size() || ib >= B.size()) break;

                uint32_t step = left;
                if (ca) step = std::min(step, A[ia].ce - pa);
                if (cb) step = std::min(step, B[ib].ce - pb);

                if (step == 0) {
                    flush();
                    advance_if_needed();
                    ba = pa;
                    bb = pb;
                    continue;
                }

                if (!started) {
                    started = true;
                    seg_ia = ia;
                    seg_ib = ib;
                    ba = pa;
                    bb = pb;
                }

                if (!ops.empty() && ops.back().op == cop.op) {
                    ops.back().len += step;
                } else {
                    ops.emplace_back(step, cop.op);
                }

                if (ca) pa += step;
                if (cb) pb += step;
                left -= step;

                const bool hit_a = ca && pa == A[ia].ce;
                const bool hit_b = cb && pb == B[ib].ce;

                if (hit_a || hit_b) {
                    flush();
                    advance_if_needed();
                    ba = pa;
                    bb = pb;
                }
            }
        }

        flush();
        return out;
    };

    // ------------------------------------------------ Alignment ------------------------------------------------
    std::string seq_a, seq_b;
    std::vector<Piece> piecesA, piecesB;
    build(subA, spanA, seq_a, piecesA);
    build(subB, spanB, seq_b, piecesB);

    if (seq_a.empty() || seq_b.empty()) return {};

    const minimizerdna::MinimizerBuilder mzb(mm_opt_);
    uint32_t k_w_max = std::max({
        static_cast<uint32_t>(chainOpts_.k),
        static_cast<uint32_t>(chainOpts_.w),
        static_cast<uint32_t>(mm_opt_.k),
        static_cast<uint32_t>(mm_opt_.w)
    });

    if (seq_a.size() <= k_w_max || seq_b.size() <= k_w_max) {
        if (seq_a != seq_b && std::min(seq_a.size(), seq_b.size()) < static_cast<size_t>(MIN_EQ_FOR_CUT_)) {
            return {};
        }
    } else if (mzb.max_containment(seq_a, seq_b) < MIN_JACCARD_FOR_ALIGN_) {
        return {};
    }

    std::string name_a, name_b;
    for (const auto& p : piecesA) {
        if (!name_a.empty()) name_a += ";";
        name_a += getNodeName(NodeHandle::get_segment_id(p.vtx));
    }
    for (const auto& p : piecesB) {
        if (!name_b.empty()) name_b += ";";
        name_b += getNodeName(NodeHandle::get_segment_id(p.vtx));
    }

    auto big = align_and_pack_(
        (NodeHandle::get_segment_id(piecesA.front().vtx) << 1) | 0u,
        (NodeHandle::get_segment_id(piecesB.front().vtx) << 1) | 0u,
        name_a, name_b,
        seq_a, seq_b,
        0, static_cast<uint32_t>(seq_a.size()),
        0, static_cast<uint32_t>(seq_b.size())
    );

    const std::vector<uint64_t> repeat_mask_a = build_tandem_repeat_mask_(seq_a, REPEAT_MASK_MIN_LEN_, REPEAT_MASK_MAX_PERIOD_, REPEAT_MASK_MAX_MISMATCH_);
    const std::vector<uint64_t> repeat_mask_b = build_tandem_repeat_mask_(seq_b, REPEAT_MASK_MIN_LEN_, REPEAT_MASK_MAX_PERIOD_, REPEAT_MASK_MAX_MISMATCH_);

    std::vector<BubbleAlignment> outs;
    for (auto& aln : big) {
        if (aln.ops.size() > 1) mask_repetitive_match_ops_(aln.ops, repeat_mask_a, repeat_mask_b, aln.beg_a, aln.beg_b, REPEAT_MASK_MIN_LEN_, REPEAT_MASK_MAX_PERIOD_);

        auto v = split(piecesA, piecesB, aln);
        outs.insert(outs.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
    }

    // Filter out self-alignments
    outs.erase(
        std::remove_if(
            outs.begin(), outs.end(),
            [](const BubbleAlignment& aln) {
                return NodeHandle::get_segment_id(aln.v_a) == NodeHandle::get_segment_id(aln.v_b);
            }
        ),
        outs.end()
    );

    if (DEBUG_ENABLED) {
        for (auto& out : outs) {
            uint32_t da_beg = out.beg_a;
            uint32_t da_end = out.end_a;
            uint32_t db_beg = out.beg_b;
            uint32_t db_end = out.end_b;

            if (NodeHandle::get_is_reverse(out.v_a)) { std::swap(da_beg, da_end); }
            if (NodeHandle::get_is_reverse(out.v_b)) { std::swap(db_beg, db_end); }

            debug_stream() << "        - " << out.name_a << " vs " << out.name_b << " | " << da_beg << "-" << da_end << " vs " << db_beg << "-" << db_end << " | " << CIGAR::pack(out.ops) << "\n";
        }
    }

    return outs;
}


std::vector<GfaCollapser::BubbleAlignment> GfaCollapser::align_paths_(
    const std::vector<std::vector<uint32_t>>& paths, 
    const bool skip_same_start,
    std::shared_ptr<path_pair_split::ComparedIndex> shared_cmp,
    std::shared_ptr<std::shared_mutex> shared_cmp_mutex, 
    uint32_t idx
) {
    std::vector<BubbleAlignment> outs;
    path_pair_split::ComparedIndex* cmp = shared_cmp.get();

    // ------------------------------------------------ Helper functions ------------------------------------------------
    auto name_of = [this](uint32_t v) {
        return getNodeName(NodeHandle::get_segment_id(v));
    };
    auto len_of = [this](uint32_t v) {
        return getNodeLength(NodeHandle::get_segment_id(v));
    };
    auto rev_of = [](uint32_t v) {
        return NodeHandle::get_is_reverse(v);
    };
    auto overlap_of = [this](uint32_t v, uint32_t w) {
        return get_edge_ow(v, w);
    };

    auto remember_node_alignment = [&](const BubbleAlignment& aln) {
        if (aln.ops.empty()) return;

        std::vector<uint32_t> a{aln.v_a};
        std::vector<uint32_t> b{aln.v_b};

        std::vector<path_pair_split::COp> ops;
        ops.reserve(aln.ops.size());
        for (const auto& op : aln.ops) ops.emplace_back(op.len, op.op);

        auto to_oriented_beg = [&](uint32_t vtx, uint32_t fwd_beg, uint32_t fwd_end) {
            uint32_t L = getNodeLength(NodeHandle::get_segment_id(vtx));
            return NodeHandle::get_is_reverse(vtx) ? L - fwd_end : fwd_beg;
        };

        std::unique_lock<std::shared_mutex> lock(*shared_cmp_mutex);
        path_pair_split::add_cigar(
            cmp,
            a, b,
            name_of, len_of, rev_of, overlap_of,
            ops,
            to_oriented_beg(aln.v_a, aln.beg_a, aln.end_a),
            to_oriented_beg(aln.v_b, aln.beg_b, aln.end_b)
        );
    };

    auto append = [&](std::vector<BubbleAlignment>&& alns) {
        if (alns.empty()) return;
        for (const auto& aln : alns) remember_node_alignment(aln);
        for (auto& aln : alns) {
            aln.idx = idx;
            if (!aln.ops.empty()) outs.emplace_back(std::move(aln));
        }
    };

    // ------------------------------------------------ Alignment ------------------------------------------------
    for (size_t i = 0; i + 1 < paths.size(); ++i) {
        for (size_t j = i + 1; j < paths.size(); ++j) {
            const auto& pathA = paths[i];
            const auto& pathB = paths[j];
            if (pathA.empty() || pathB.empty()) continue;
            if (skip_same_start && pathA[0] == pathB[0]) continue;  // Skip pairs of paths that start at the same node for homologous path alignment

            std::vector<path_pair_split::SplitTask<uint32_t>> tasks;

            {
                std::shared_lock<std::shared_mutex> lock(*shared_cmp_mutex);
                tasks = path_pair_split::split_path_pair(
                    pathA, pathB,
                    name_of, len_of, rev_of, overlap_of,
                    cmp
                );
            }

            if (DEBUG_ENABLED) {
                const std::string path_a_name = path_pair_split::join_path_names(pathA, name_of);
                const std::string path_b_name = path_pair_split::join_path_names(pathB, name_of);

                debug_stream() << "\n";
                debug_stream() << "Tasks=" << tasks.size() << "\n";
                debug_stream() << "  - ref: " << path_a_name << "\n";
                debug_stream() << "  - qry: " << path_b_name << "\n";

                for (size_t ti = 0; ti < tasks.size(); ++ti) {
                    const auto& t = tasks[ti];

                    const std::string sub_a = path_pair_split::join_span_names(pathA, t.a_span, name_of, len_of, rev_of, overlap_of);
                    const std::string sub_b = path_pair_split::join_span_names(pathB, t.b_span, name_of, len_of, rev_of, overlap_of);

                    debug_stream() << "    - task[" << ti << "] " << path_pair_split::task_type_name(t) << " | " << sub_a << " vs " << sub_b << "\n";
                }
            }

            for (const auto& t : tasks) {
                if (t.type == path_pair_split::SplitTask<uint32_t>::Type::ComparedAnchor) {
                    continue;
                }

                append(align_subpaths_(pathA, pathB, t.a_span, t.b_span));
            }
        }
    }

    return outs;
}

void GfaCollapser::homologous_align_(
    const std::vector<GfaBubble::Bubble>& bubbles, 
    const std::vector<GfaBubble::ForkGroup>& forks, 
    const std::vector<GfaBubble::HomologousPath>& homologous_paths,
    const GfaBubble::GfaBubbleFinder& bubble_finder
) {
    log_stream() << "Aligning forks, bubbles and homologous paths ...\n";

    // Thread pool
    ThreadPool pool(alignOpts_.threads);
    std::deque<std::future<std::vector<BubbleAlignment>>> futs;

    // To avoid excessive memory usage
    const size_t max_inflight = std::max<size_t>(1, alignOpts_.threads * 4);

    // Shared index for path pair comparison to speed up alignment of multiple path pairs with shared subpaths
    auto shared_cmp = std::make_shared<path_pair_split::ComparedIndex>();
    auto shared_cmp_mutex = std::make_shared<std::shared_mutex>();

    // Used for deduplication
    uint32_t submit_idx = 0;

    // ------------------------------------------------ Helper functions ------------------------------------------------
    auto consume_alignments = [&](std::vector<BubbleAlignment>&& vec) {
        for (auto& aln : vec) {
            if (!aln.ops.empty()) {
                bubble_aligns_.emplace_back(std::move(aln));
            }
        }
    };

    auto collect_ready = [&]() {
        for (auto it = futs.begin(); it != futs.end(); ) {
            if (it->wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
                consume_alignments(it->get());
                it = futs.erase(it);
            } else {
                ++it;
            }
        }
    };

    auto submit_alignment = [&](auto&& task) {
        futs.emplace_back(pool.submit(std::forward<decltype(task)>(task)));

        if (futs.size() >= max_inflight) {
            collect_ready();

            if (futs.size() >= max_inflight) {
                consume_alignments(futs.front().get());
                futs.pop_front();
            }
        }
    };

    // ------------------------------------------------ Alignment ------------------------------------------------
    ProgressTracker prog(bubbles.size() + homologous_paths.size() + forks.size());

    // Bubbles
    for (const auto& bb : bubbles) {
        prog.hit();

        const auto paths = bb.get_paths();
        if (paths.size() < 2) continue;

        submit_alignment([this, paths, shared_cmp, shared_cmp_mutex, idx = submit_idx++]() {
            return align_paths_(paths, false, shared_cmp, shared_cmp_mutex, idx);
        });
    }

    // Homologous paths
    for (const auto& hp : homologous_paths) {
        prog.hit();

        const auto paths = hp.get_paths_u32();
        if (paths.size() < 2) continue;

        // Reduce cycle in the graph. (2026-05-25, v0.1.3-r6)
        bool skip_same_start = (hp.type == GfaBubble::Type::DiffSource);

        submit_alignment([this, paths, shared_cmp, shared_cmp_mutex, skip_same_start, idx = submit_idx++]() {
            return align_paths_(paths, skip_same_start, shared_cmp, shared_cmp_mutex, idx);  // Skip pairs of paths that start at the same node for homologous path alignment
        });
    }

    // Forks
    for (const auto& fork : forks) {
        prog.hit();

        const auto branches = fork.get_branches();
        if (branches.size() < 2) continue;

        std::vector<std::vector<uint32_t>> branch_vecs;
        for (uint32_t b : branches) {
            branch_vecs.push_back({b});
        }

        submit_alignment([this, branch_vecs = std::move(branch_vecs), shared_cmp, shared_cmp_mutex, idx = submit_idx++]() {
            return align_paths_(branch_vecs, false, shared_cmp, shared_cmp_mutex, idx);
        });
    }

    // ------------------------------------------------ Collect results ------------------------------------------------
    while (!futs.empty()) {
        collect_ready();

        if (!futs.empty()) {
            consume_alignments(futs.front().get());
            futs.pop_front();
        }
    }

    pool.stop();

    log_stream() << "  - Total alignments produced: " << bubble_aligns_.size() << "\n" << "\n";
}

void GfaCollapser::collapse_homologous_seq(
    const std::vector<GfaBubble::Bubble>& bubbles,
    const std::vector<GfaBubble::ForkGroup>& forks,
    const std::vector<GfaBubble::HomologousPath>& homologous_paths, 
    const std::string& prefix,
    const GfaBubble::GfaBubbleFinder& bubble_finder
) {
    log_stream() << "Transforming graph to de-overlap form and collapsing homologous sequences ...\n" << "\n";

    // 0. clear old state
    bubble_aligns_.clear();
    cuts_.clear();
    rulemap_.clear();

    // 1. normalize overlap metadata first
    finalize_();
    prune_overlaps_();

    // 2. initialize all segment cuts: boundaries
    initialize_cuts_();

    // 3. add overlap-based alignments
    overlaps_align_();

    // 4. add bubble-based alignments
    homologous_align_(bubbles, forks, homologous_paths, bubble_finder);

    // 5. group alignments by connected segments
    const auto align_groups = build_align_groups_();

    // 6. build cuts, propagate cuts, and prune abnormal cuts
    build_propagate_prune_cuts_(align_groups);
    if (DEBUG_ENABLED) print_cuts(cuts_);

    // 7. build rulemap from final cuts
    build_rulemap_(align_groups);
    SegReplace::Expander ex = build_SegReplace_(false);
    ex.save_map(prefix + ".collapse.map");

    // 8. apply once
    expand_and_rewire_edges_(ex);

    // 9. cleanup
    remove_unused_nodes_();

    // 10. merge linear chains
    merge_linear_chains();

    log_stream() << "Finished transforming graph to de-overlap form and collapsing homologous sequences.\n" << "\n";
}
