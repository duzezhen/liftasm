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


uint64_t LongNodeSplitter::split(GfaGraph& graph, uint32_t min_length, uint32_t chunk_length) {
    log_stream() << "Splitting long segments" << " (min_length=" << min_length << ", chunk_length=" << chunk_length << ") ..." << std::endl;

    if (min_length == 0 || chunk_length == 0) return 0;

    const size_t old_node_count = graph.nodes_.size();
    const size_t old_arc_count = graph.arcs_.size();
    std::vector<uint32_t> left_overlap(old_node_count, 0);
    std::vector<uint32_t> right_overlap(old_node_count, 0);

    // Calculate the maximum overlap length for each segment
    for (size_t i = 0; i < old_arc_count; ++i) {
        const GfaArc& a = graph.arcs_[i];
        if (a.get_del()) continue;
        const uint32_t vs = a.get_source_segment_id();
        const uint32_t ws = a.get_target_segment_id();
        const uint32_t ov = a.ov > 0 ? static_cast<uint32_t>(a.ov) : 0;
        const uint32_t ow = a.ow > 0 ? static_cast<uint32_t>(a.ow) : 0;
        if (vs < old_node_count) {
            uint32_t& end = a.get_source_is_reverse() ? left_overlap[vs] : right_overlap[vs];
            end = std::max(end, ov);
        }
        if (ws < old_node_count) {
            uint32_t& end = a.get_target_is_reverse() ? right_overlap[ws] : left_overlap[ws];
            end = std::max(end, ow);
        }
    }

    // S-lines
    std::vector<std::vector<uint32_t>> pieces(old_node_count);
    const gfaName namer;
    uint64_t split_count = 0;
    for (uint32_t sid = 0; sid < old_node_count; ++sid) {
        const GfaNode original = graph.nodes_[sid];
        if (original.deleted || original.length <= min_length || original.is_complex) continue;

        const uint32_t left = std::min(left_overlap[sid], original.length);
        const uint32_t right = std::min(right_overlap[sid], original.length);
        if (uint64_t(left) + right >= original.length) continue;  // terminal overlap regions intersect

        std::vector<uint32_t> cuts{0};
        if (left > 0) cuts.push_back(left);
        uint32_t pos = left;
        const uint32_t middle_end = original.length - right;
        while (uint64_t(pos) + chunk_length < middle_end) {
            pos += chunk_length;
            cuts.push_back(pos);
        }
        if (middle_end > cuts.back()) cuts.push_back(middle_end);
        if (cuts.back() < original.length) cuts.push_back(original.length);
        if (cuts.size() <= 2) continue;

        auto& ids = pieces[sid];
        ids.reserve(cuts.size() - 1);
        for (size_t p = 1; p < cuts.size(); ++p) {
            const uint32_t beg = cuts[p - 1], end = cuts[p];
            std::string seq = "*";
            if (original.sequence != "*" && end <= original.sequence.size()) {
                seq = original.sequence.substr(beg, end - beg);
            }
            const std::string name = namer.format_interval_name(original.name, beg, end, false);
            const uint32_t id = graph.add_segment(name, seq, true, original.sample_ids);
            GfaNode& node = graph.nodes_[id];
            if (seq == "*") graph.total_segment_length_ += end - beg;
            node.length = end - beg;
            node.is_complex = original.is_complex;
            node.aux = original.aux;
            ids.push_back(id);
        }
        ++split_count;
    }
    if (split_count == 0) {
        log_stream() << "  - Split " << split_count << " segment(s), creating " << 0 << " new segment(s)\n\n";
        return 0;
    }

    auto endpoint = [&](uint32_t vertex, bool source) {
        const uint32_t sid = vertex >> 1;
        if (sid >= pieces.size() || pieces[sid].empty()) return vertex;
        const bool rev = vertex & 1u;
        const auto& ps = pieces[sid];
        const size_t index = source ? (rev ? 0 : ps.size() - 1) : (rev ? ps.size() - 1 : 0);
        return Vertex::make_vertex(ps[index], rev);
    };

    size_t internal_arc_count = 0;
    for (const auto& ps : pieces) {
        if (ps.size() > 1) internal_arc_count += ps.size() - 1;
    }
    const size_t arc_capacity = graph.arcs_.size() + old_arc_count + internal_arc_count;
    graph.arcs_.reserve(arc_capacity);
    graph.link_aux_.reserve(arc_capacity);

    // L-lines
    for (size_t i = 0; i < old_arc_count; ++i) {
        const GfaArc a = graph.arcs_[i];
        if (a.get_del()) continue;
        const uint32_t v = endpoint(a.get_source_vertex_id(), true);
        const uint32_t w = endpoint(a.get_target_vertex_id(), false);
        if (v == a.get_source_vertex_id() && w == a.get_target_vertex_id()) continue;
        GfaArc* copy = graph.add_arc(v, w, a.ov, a.ow, -1, a.get_comp());
        copy->rank = a.rank;
        copy->set_strong(a.get_strong());
    }
    for (uint32_t sid = 0; sid < old_node_count; ++sid) {
        const auto& ps = pieces[sid];
        if (ps.empty()) continue;
        for (size_t i = 1; i < ps.size(); ++i) {
            graph.add_arc(Vertex::make_vertex(ps[i - 1], false), Vertex::make_vertex(ps[i], false), 0, 0, -1, false);
        }
        graph.delete_segment(sid);
    }

    // P-lines
    for (GfaPath& path : graph.paths_) {
        std::vector<PathSegment> expanded;
        expanded.reserve(path.segments.size());
        for (const PathSegment& s : path.segments) {
            if (s.node_id >= pieces.size() || pieces[s.node_id].empty()) {
                expanded.push_back(s);
                continue;
            }
            const auto& ps = pieces[s.node_id];
            if (s.is_reverse) {
                for (auto it = ps.rbegin(); it != ps.rend(); ++it) expanded.push_back({*it, true});
            } else {
                for (uint32_t id : ps) expanded.push_back({id, false});
            }
        }
        path.segments.swap(expanded);
    }

    graph.rebuild_after_edits();

    const size_t new_node_count = graph.nodes_.size();

    log_stream() << "  - Split " << split_count << " segment(s), creating " << (new_node_count - old_node_count) << " new segment(s)\n\n";

    return split_count;
}


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


/* -------------------------------------------- Normalize Homopolymer Bubbles -------------------------------------------- */
static void print_alignment_with_mask_(
    const std::string& ref_seq,
    const std::string& qry_seq,
    const std::vector<uint64_t>& ref_mask,
    const std::vector<uint64_t>& qry_mask,
    const std::string& cigar
) {
    auto test_bit = [](const std::vector<uint64_t>& mask, uint32_t pos) -> bool {
        const uint32_t w = pos >> 6;
        return w < mask.size() && (mask[w] & (uint64_t(1) << (pos & 63)));
    };

    std::string r_aln;
    std::string q_aln;
    std::string mid;
    std::string r_msk;
    std::string q_msk;

    const std::vector<CIGAR::COp> ops = CIGAR::parse(cigar);

    uint32_t pr = 0;
    uint32_t pq = 0;

    for (const CIGAR::COp& op : ops) {
        for (uint32_t k = 0; k < op.len; ++k) {
            switch (op.op) {
                case '=':
                case 'M': {
                    const char rc = pr < ref_seq.size() ? ref_seq[pr] : '?';
                    const char qc = pq < qry_seq.size() ? qry_seq[pq] : '?';

                    r_aln += rc;
                    q_aln += qc;
                    mid += (rc == qc ? '|' : '.');

                    r_msk += test_bit(ref_mask, pr) ? '.' : ' ';
                    q_msk += test_bit(qry_mask, pq) ? '.' : ' ';

                    ++pr;
                    ++pq;
                } break;

                case 'X': {
                    const char rc = pr < ref_seq.size() ? ref_seq[pr] : '?';
                    const char qc = pq < qry_seq.size() ? qry_seq[pq] : '?';

                    r_aln += rc;
                    q_aln += qc;
                    mid += '*';

                    r_msk += test_bit(ref_mask, pr) ? '.' : ' ';
                    q_msk += test_bit(qry_mask, pq) ? '.' : ' ';

                    ++pr;
                    ++pq;
                } break;

                case 'D': {
                    const char rc = pr < ref_seq.size() ? ref_seq[pr] : '?';

                    r_aln += rc;
                    q_aln += '-';
                    mid += ' ';

                    r_msk += test_bit(ref_mask, pr) ? '.' : ' ';
                    q_msk += ' ';

                    ++pr;
                } break;

                case 'I': {
                    const char qc = pq < qry_seq.size() ? qry_seq[pq] : '?';

                    r_aln += '-';
                    q_aln += qc;
                    mid += ' ';

                    r_msk += ' ';
                    q_msk += test_bit(qry_mask, pq) ? '.' : ' ';

                    ++pq;
                } break;

                case 'S': {
                    const char qc = pq < qry_seq.size() ? qry_seq[pq] : '?';

                    r_aln += '-';
                    q_aln += qc;
                    mid += ' ';

                    r_msk += ' ';
                    q_msk += test_bit(qry_mask, pq) ? '.' : ' ';

                    ++pq;
                } break;

                case 'H':
                    break;

                default:
                    r_aln += '?';
                    q_aln += '?';
                    mid += '?';
                    r_msk += ' ';
                    q_msk += ' ';
                    break;
            }
        }
    }

    debug_stream() << "  - " << " mask" << ": " << r_msk << "\n";
    debug_stream() << "  - " << "r_seq" << ": " << r_aln << "\n";
    debug_stream() << "  - " << "cigar" << ": " << mid << "\n";
    debug_stream() << "  - " << "q_seq" << ": " << q_aln << "\n";
    debug_stream() << "  - " << " mask" << ": " << q_msk << "\n";
}

void GfaCollapser::normalize_homopolymer_bubbles(
    const std::vector<GfaBubble::Bubble>& bubbles,
    uint32_t max_path_len
) {
    log_stream() << "Normalizing homopolymer bubbles ...\n";

    struct Branch {
        uint32_t vertex{UINT32_MAX};
        uint32_t segment{UINT32_MAX};
        std::string seq;
        std::vector<uint64_t> mask;
    };

    struct Insert {
        uint32_t pos{0};
        std::string seq;
    };

    auto trim_path = [](const std::vector<uint32_t>& path, uint32_t source, uint32_t sink) {
        size_t beg = 0;
        size_t end = path.size();

        if (beg < end && path[beg] == source) ++beg;
        if (beg < end && path[end - 1] == sink) --end;

        return std::vector<uint32_t>(path.begin() + beg, path.begin() + end);
    };

    auto mask_all = [](const std::vector<uint64_t>& mask, uint32_t beg, uint32_t end) {
        if (beg >= end) return false;

        for (uint32_t p = beg; p < end; ++p) {
            const uint32_t w = p >> 6;

            if (w >= mask.size()) return false;
            if ((mask[w] & (uint64_t(1) << (p & 63))) == 0) return false;
        }

        return true;
    };

    auto mask_any = [](const std::vector<uint64_t>& mask, uint32_t beg, uint32_t end) {
        if (beg >= end) return false;

        for (uint32_t p = beg; p < end; ++p) {
            const uint32_t w = p >> 6;

            if (w < mask.size() && (mask[w] & (uint64_t(1) << (p & 63)))) {
                return true;
            }
        }

        return false;
    };

    auto has_mask = [&](const std::vector<uint64_t>& mask) {
        for (uint64_t x : mask) {
            if (x != 0) return true;
        }

        return false;
    };

    auto coord_like_name = [](const std::string& name) {
        const size_t c = name.rfind(':');
        const size_t d = name.rfind('-');

        return c != std::string::npos && d != std::string::npos && c < d && d + 1 < name.size();
    };

    auto next_clean_id = [&]() -> uint64_t {
        uint64_t next = 0;

        for (const GfaNode& n : nodes_) {
            if (n.deleted) continue;

            const std::string& name = n.name;

            if (name.size() <= 5) continue;
            if (name.compare(0, 5, "clean") != 0) continue;

            bool numeric = true;

            for (size_t i = 5; i < name.size(); ++i) {
                if (!std::isdigit(static_cast<unsigned char>(name[i]))) {
                    numeric = false;
                    break;
                }
            }

            if (!numeric) continue;

            const uint64_t id = std::stoull(name.substr(5));
            next = std::max<uint64_t>(next, id + 1);
        }

        return next;
    };

    uint64_t clean_id = next_clean_id();

    auto rename_if_needed = [&](uint32_t sid, uint32_t len) {
        if (sid >= nodes_.size()) return;

        if (!coord_like_name(nodes_[sid].name)) return;

        name_to_id_map_.erase(nodes_[sid].name);

        while (true) {
            std::string new_name = "clean" + std::to_string(clean_id++) + ":0-" + std::to_string(len);

            if (name_to_id_map_.find(new_name) == name_to_id_map_.end()) {
                nodes_[sid].name = std::move(new_name);
                name_to_id_map_[nodes_[sid].name] = sid;
                return;
            }
        }
    };

    auto build_branch = [&](const std::vector<uint32_t>& raw_path, uint32_t source, uint32_t sink, Branch& b) {
        const std::vector<uint32_t> inner = trim_path(raw_path, source, sink);

        if (inner.size() != 1) return false;

        b.vertex = inner.front();
        b.segment = NodeHandle::get_segment_id(b.vertex);

        if (b.segment >= nodes_.size()) return false;
        if (nodes_[b.segment].deleted || nodes_[b.segment].is_complex) return false;
        if (nodes_[b.segment].sequence.empty() || nodes_[b.segment].sequence == "*") return false;

        b.seq = get_oriented_sequence(Vertex(b.vertex));

        if (b.seq.empty() || b.seq == "*") return false;
        if (b.seq.size() > max_path_len) return false;

        b.mask = build_tandem_repeat_mask_(
            b.seq,
            REPEAT_MASK_MIN_LEN_,
            REPEAT_MASK_MAX_PERIOD_,
            /*max_mismatch=*/0
        );

        return true;
    };

    auto repeat_interval_ok = [&](const Branch& b, uint32_t beg, uint32_t end) {
        if (mask_all(b.mask, beg, end)) return true;

        const uint32_t n = static_cast<uint32_t>(b.seq.size());
        const uint32_t pad = std::max<uint32_t>(REPEAT_MASK_MAX_PERIOD_, end - beg);

        const uint32_t ctx_beg = beg > pad ? beg - pad : 0;
        const uint32_t ctx_end = std::min<uint32_t>(n, end + pad);

        return mask_any(b.mask, ctx_beg, ctx_end);
    };

    auto rebuild_with_inserts = [](
        const std::string& seq,
        std::vector<Insert>& inserts,
        std::string& out
    ) {
        if (inserts.empty()) return false;

        std::stable_sort(
            inserts.begin(),
            inserts.end(),
            [](const Insert& a, const Insert& b) {
                return a.pos < b.pos;
            }
        );

        size_t extra_len = 0;

        for (const Insert& ins : inserts) {
            extra_len += ins.seq.size();
        }

        out.clear();
        out.reserve(seq.size() + extra_len);

        uint32_t cur = 0;

        for (const Insert& ins : inserts) {
            const uint32_t p = std::min<uint32_t>(ins.pos, seq.size());

            if (cur < p) {
                out.append(seq.data() + cur, p - cur);
            }

            out += ins.seq;
            cur = p;
        }

        if (cur < seq.size()) {
            out.append(seq.data() + cur, seq.size() - cur);
        }

        return out.size() > seq.size();
    };

    auto make_new_oriented_seq = [&](
        const Branch& ref,
        const Branch& query,
        std::string& new_ref_seq,
        std::string& new_query_seq,
        std::string& cigar_out
    ) {
        const auto hits = align_short_mm2_(
        // const auto hits = align_short_wfa_(
            getNodeName(ref.segment),
            getNodeName(query.segment),
            ref.seq,
            query.seq
        );

        const auto hit = std::find_if(hits.begin(), hits.end(), [&](const MmWfaHit& h) {
            return h.r_beg == 0 && h.q_beg == 0 && h.r_end == ref.seq.size() && h.q_end == query.seq.size();
        });
        if (hit == hits.end()) return false;

        const std::string& cigar = hit->cigar;
        cigar_out = cigar;

        if (cigar.empty() || cigar == "*") return false;

        const std::vector<CIGAR::COp> ops = CIGAR::parse(cigar);

        std::vector<Insert> ref_inserts;
        std::vector<Insert> query_inserts;

        uint32_t pa = 0;
        uint32_t pb = 0;

        auto choose_insert_pos_for_repeat_base = [](
            const std::string& seq,
            uint32_t x,
            char base,
            uint32_t max_period
        ) -> uint32_t {
            const uint32_t n = static_cast<uint32_t>(seq.size());

            if (x >= n || max_period == 0) {
                return std::min<uint32_t>(x, n);
            }

            struct Score {
                uint32_t support{0};
                uint32_t period{1};
            };

            auto evaluate = [&](uint32_t pos) {
                const uint32_t m = n + 1;

                auto at = [&](uint32_t i) -> char {
                    if (i < pos) return seq[i];
                    if (i == pos) return base;
                    return seq[i - 1];
                };

                Score best;

                for (uint32_t p = 1; p <= max_period && p < m; ++p) {
                    uint32_t left = 0;
                    uint32_t right = 0;

                    for (uint32_t i = pos; i >= p; --i) {
                        if (at(i) != at(i - p)) break;
                        ++left;
                        if (i == p) break;
                    }

                    for (uint32_t i = pos + p; i < m; ++i) {
                        if (at(i) != at(i - p)) break;
                        ++right;
                    }

                    const uint32_t support = left + right;

                    const uint64_t lhs = uint64_t(support) * best.period;
                    const uint64_t rhs = uint64_t(best.support) * p;

                    if (lhs > rhs || (lhs == rhs && support > best.support) || (lhs == rhs && support == best.support && p < best.period)) {
                        best.support = support;
                        best.period = p;
                    }
                }

                return best;
            };

            const Score before = evaluate(x);
            const Score after = evaluate(x + 1);

            const uint64_t before_units = uint64_t(before.support) * after.period;
            const uint64_t after_units = uint64_t(after.support) * before.period;

            if (after_units > before_units) {
                return x + 1;
            }

            if (before_units > after_units) {
                return x;
            }

            if (after.support > before.support) {
                return x + 1;
            }

            if (before.support > after.support) {
                return x;
            }

            if (after.period < before.period) {
                return x + 1;
            }

            return x;
        };

        for (const CIGAR::COp& op : ops) {
            switch (op.op) {
                case '=':
                case 'M': {
                    pa += op.len;
                    pb += op.len;
                } break;

                case 'X': {
                    for (uint32_t k = 0; k < op.len; ++k) {
                        const uint32_t xa = pa + k;
                        const uint32_t xb = pb + k;

                        const bool ref_repeat = mask_all(ref.mask, xa, xa + 1);
                        const bool query_repeat = mask_all(query.mask, xb, xb + 1);

                        if (ref_repeat && query_repeat) {
                            const char ref_base = ref.seq[xa];
                            const char query_base = query.seq[xb];

                            const uint32_t query_pos = choose_insert_pos_for_repeat_base(query.seq, xb, ref_base, REPEAT_MASK_MAX_PERIOD_);
                            const uint32_t ref_pos = choose_insert_pos_for_repeat_base(ref.seq, xa, query_base, REPEAT_MASK_MAX_PERIOD_);

                            query_inserts.push_back(Insert{query_pos, std::string(1, ref_base)});
                            ref_inserts.push_back(Insert{ref_pos, std::string(1, query_base)});
                        }
                    }

                    pa += op.len;
                    pb += op.len;
                } break;

                case 'D': {
                    if (repeat_interval_ok(ref, pa, pa + op.len)) {
                        query_inserts.push_back(Insert{pb, ref.seq.substr(pa, op.len)});
                    }

                    pa += op.len;
                } break;

                case 'I': {
                    if (repeat_interval_ok(query, pb, pb + op.len)) {
                        ref_inserts.push_back(Insert{pa, query.seq.substr(pb, op.len)});
                    }

                    pb += op.len;
                } break;

                case 'S':
                case 'H':
                    break;

                default:
                    return false;
            }
        }

        const bool changed_ref = rebuild_with_inserts(ref.seq, ref_inserts, new_ref_seq);
        const bool changed_query = rebuild_with_inserts(query.seq, query_inserts, new_query_seq);

        return changed_ref || changed_query;
    };

    auto apply_oriented_seq = [&](const Branch& b, const std::string& oriented_seq) {
        if (oriented_seq.empty()) return false;
        if (b.segment >= nodes_.size()) return false;
        if (nodes_[b.segment].deleted || nodes_[b.segment].sequence == "*") return false;

        const uint32_t old_len = nodes_[b.segment].length;

        if (NodeHandle::get_is_reverse(b.vertex)) {
            nodes_[b.segment].sequence = seqUtils::revcomp(oriented_seq);
        } else {
            nodes_[b.segment].sequence = oriented_seq;
        }

        nodes_[b.segment].length = static_cast<uint32_t>(nodes_[b.segment].sequence.size());

        if (nodes_[b.segment].length <= old_len) return false;

        total_segment_length_ += nodes_[b.segment].length - old_len;

        rename_if_needed(b.segment, nodes_[b.segment].length);

        return true;
    };

    size_t normalized_bubbles = 0;
    size_t extended_nodes = 0;
    uint64_t added_bp = 0;

    ProgressTracker prog(bubbles.size());

    for (const GfaBubble::Bubble& bubble : bubbles) {
        prog.hit();

        const uint32_t source = bubble.get_source();
        const uint32_t sink = bubble.get_sink();

        if (source == UINT32_MAX || sink == UINT32_MAX) continue;

        const auto& paths = bubble.get_paths();

        if (paths.size() < 2) continue;

        std::vector<Branch> branches;
        branches.reserve(paths.size());

        bool ok = true;
        bool has_repeat = false;

        for (const auto& path : paths) {
            Branch b;

            if (!build_branch(path, source, sink, b)) {
                ok = false;
                break;
            }

            has_repeat = has_repeat || has_mask(b.mask);
            branches.emplace_back(std::move(b));
        }

        if (!ok || branches.size() < 2) continue;
        if (!has_repeat) continue;

        size_t template_i = 0;

        for (size_t i = 1; i < branches.size(); ++i) {
            if (branches[i].seq.size() > branches[template_i].seq.size()) {
                template_i = i;
            }
        }

        const Branch& tmpl = branches[template_i];

        bool changed_this_bubble = false;

        if (DEBUG_ENABLED) {
            debug_stream() << "Branches: " << branches.size() << "\n";
        }

        for (size_t i = 0; i < branches.size(); ++i) {
            if (i == template_i) continue;

            Branch& target = branches[i];

            std::string new_ref_seq;
            std::string new_query_seq;
            std::string cigar;

            if (!make_new_oriented_seq(tmpl, target, new_ref_seq, new_query_seq, cigar)) { continue; }

            if (!new_ref_seq.empty() && new_ref_seq.size() <= max_path_len) {
                const uint32_t old_len = static_cast<uint32_t>(tmpl.seq.size());
                const std::string old_name = getNodeName(tmpl.segment);
                const std::string old_seq = tmpl.seq;

                if (apply_oriented_seq(tmpl, new_ref_seq)) {
                    const uint32_t delta = static_cast<uint32_t>(new_ref_seq.size()) - old_len;

                    added_bp += delta;
                    ++extended_nodes;
                    changed_this_bubble = true;

                    if (DEBUG_ENABLED) {
                        debug_stream() << "  - change: " << old_name << " -> " << getNodeName(tmpl.segment) << "\n";
                        debug_stream() << "  - cigar: " << cigar << "\n";
                        print_alignment_with_mask_(old_seq, target.seq, tmpl.mask, target.mask, cigar);
                        debug_stream() << "  - n_seq: " << new_ref_seq << "\n";
                    }
                }
            }

            if (!new_query_seq.empty() && new_query_seq.size() <= max_path_len) {
                const uint32_t old_len = static_cast<uint32_t>(target.seq.size());
                const std::string old_name = getNodeName(target.segment);
                const std::string old_seq = target.seq;

                if (apply_oriented_seq(target, new_query_seq)) {
                    const uint32_t delta = static_cast<uint32_t>(new_query_seq.size()) - old_len;

                    added_bp += delta;
                    ++extended_nodes;
                    changed_this_bubble = true;

                    if (DEBUG_ENABLED) {
                        debug_stream() << "  - change: " << old_name << " -> " << getNodeName(target.segment) << "\n";
                        debug_stream() << "  - cigar: " << cigar << "\n";
                        print_alignment_with_mask_(tmpl.seq, old_seq, tmpl.mask, target.mask, cigar);
                        debug_stream() << "  - n_seq: " << new_query_seq << "\n";
                    }
                }
            }

            if (DEBUG_ENABLED) {
                debug_stream() << "\n";
            }
        }

        if (changed_this_bubble) {
            ++normalized_bubbles;
        }
    }

    if (extended_nodes > 0) {
        rebuild_after_edits();
    }

    log_stream() << "  - Normalized bubbles: " << normalized_bubbles << "\n";
    log_stream() << "  - Extended nodes: " << extended_nodes << "\n";
    log_stream() << "  - Added bp: " << added_bp << "\n\n";
}


/* -------------------------------------------- collapse unitigs -------------------------------------------- */
bool GfaCollapser::try_build_chain_(
    uint32_t start,
    std::vector<uint32_t>& chain,
    std::vector<size_t>& edge_idx_seq
) {
    chain.clear();
    edge_idx_seq.clear();

    uint32_t head = start;

    for (;;) {
        auto in_idxs = getArcsIdxToVertex(head, true);
        if (in_idxs.size() != 1) break;

        const GfaArc& in_e = arcs_[in_idxs[0]];
        const uint32_t pred = in_e.get_source_vertex_id();

        if (getArcsIdxFromVertex(pred, true).size() != 1) break;
        if (!same_vertex_sample(pred, head)) break;

        head = pred;
    }

    chain.push_back(head);
    uint32_t cur = head;

    while (true) {
        auto out_idxs = getArcsIdxFromVertex(cur, true);
        if (out_idxs.size() != 1) break;

        const GfaArc& arc = arcs_[out_idxs[0]];
        const uint32_t w = arc.get_target_vertex_id();

        if (getArcsIdxToVertex(w, true).size() != 1) break;
        if (!same_vertex_sample(cur, w)) break;

        edge_idx_seq.push_back(out_idxs[0]);
        chain.push_back(w);

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
        std::vector<uint32_t> sample_ids;
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
            if (sid < used.size() && (used[sid] || nodes_[sid].deleted) || nodes_[sid].is_complex) { conflict = true; break; }
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

        // sample IDs
        std::vector<uint32_t> merged_sample_ids;
        for (uint32_t vtx : chain) {
            uint32_t sid = NodeHandle::get_segment_id(vtx);
            if (sid < nodes_.size()) {
                merge_sample_ids_into(merged_sample_ids, nodes_[sid].sample_ids);
            }
        }

        Pending P;
        P.new_name   = std::move(new_name);
        P.merged_seq = std::move(merged);
        P.chain_vtx  = chain;
        P.edge_idx_seq = edge_idx_seq;
        P.sample_ids = std::move(merged_sample_ids);

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

    // Add all new segments
    for (auto &P : pendings) {
        uint32_t new_seg_id = add_segment(P.new_name, P.merged_seq, true, P.sample_ids);
        for (uint32_t vtx : P.chain_vtx) {
            const uint32_t sid = NodeHandle::get_segment_id(vtx);
            merge_variant_tags_into(nodes_[new_seg_id].aux, nodes_[sid].aux);
        }
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

    auto build = [&](const std::vector<uint32_t>& path, const path_pair_split::Span& span, std::string& seq, std::vector<Piece>& pieces) -> bool {
        seq.clear();
        pieces.clear();
        if (path.empty()) return false;

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
            const uint32_t take = ob < L ? std::min<uint32_t>(L - ob, end - cur) : 0;
            const uint32_t oe = ob + take;

            if (ob < oe) {
                std::string s = get_oriented_sequence(Vertex(vtx));
                if (s.empty() || s == "*" || oe > s.size()) {
                    seq.clear();
                    pieces.clear();
                    return false;
                }
                seq.append(s.data() + ob, oe - ob);
                pieces.push_back({vtx, out_pos, out_pos + (oe - ob), ob, oe});
                out_pos += oe - ob;
            }

            cur = layout.abs({p.node + 1, 0});
        }

        return cur >= end && !pieces.empty();
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

            // If either of the two segments is complex, ignore the alignment result (2026-06-21, 0.1.3-r12)
            if (getNodeComplex(sid_a) || getNodeComplex(sid_b)) {
                ops.clear();
                started = false;
                ba = pa;
                bb = pb;
                return;
            }

            // If the two segments share any root name, ignore the alignment result (2026-06-02, 0.1.3-r9)
            {
                std::vector<std::string> roots_a;
                std::vector<std::string> roots_b;

                gfaName::collect_roots_from_name(getNodeName(sid_a), roots_a);
                gfaName::collect_roots_from_name(getNodeName(sid_b), roots_b);

                if (gfaName::roots_intersect(roots_a, roots_b)) {
                    ops.clear();
                    started = false;
                    ba = pa;
                    bb = pb;
                    return;
                }
            }

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
                big.mapq,
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
    if (!build(subA, spanA, seq_a, piecesA) || !build(subB, spanB, seq_b, piecesB)) return {};

    const minimizerdna::MinimizerBuilder mzb(mm_opt_);
    uint32_t k_w_max = (std::max(static_cast<uint32_t>(chainOpts_.k), static_cast<uint32_t>(chainOpts_.w)) + 1);
    uint32_t hk_hw_max = (std::max(static_cast<uint32_t>(mm_opt_.k), static_cast<uint32_t>(mm_opt_.w)) + 1) * 5;

    std::string name_a, name_b;
    for (const auto& p : piecesA) {
        if (!name_a.empty()) name_a += ";";
        name_a += getNodeName(NodeHandle::get_segment_id(p.vtx));
    }
    for (const auto& p : piecesB) {
        if (!name_b.empty()) name_b += ";";
        name_b += getNodeName(NodeHandle::get_segment_id(p.vtx));
    }

    if (seq_a.size() <= k_w_max || seq_b.size() <= k_w_max) {
        if (seq_a != seq_b && std::min(seq_a.size(), seq_b.size()) < static_cast<size_t>(MIN_EQ_FOR_CUT_)) {
            debug_stream() << "    - skip: " << name_a << " vs " << name_b << " | too short and not equal" << "\n";
            return {};
        }
    }
    
    if (seq_a.size() >= hk_hw_max && seq_b.size() >= k_w_max) {
        if (mzb.max_containment(seq_a, seq_b) < MIN_JACCARD_FOR_ALIGN_) {
            debug_stream() << "    - skip: " << name_a << " vs " << name_b << " | jaccard: " << mzb.max_containment(seq_a, seq_b) << "\n";
            return {};
        }
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
    const std::vector<uint32_t>& clusters,
    GfaBubble::Type path_type,
    const bool skip_same_start,
    std::shared_ptr<path_pair_split::ComparedIndex> shared_cmp,
    std::shared_ptr<std::shared_mutex> shared_cmp_mutex, 
    uint32_t idx
) {
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

    auto path_len = [&](const std::vector<uint32_t>& p) -> uint64_t {
        uint64_t L = 0;
        for (size_t i = 0; i < p.size(); ++i) {
            uint32_t x = len_of(p[i]);
            if (i > 0) {
                const uint32_t ow = overlap_of(p[i - 1], p[i]);
                x -= std::min(x, ow);
            }
            L += x;
        }
        return L;
    };

    auto pair_key = [](uint32_t a, uint32_t b) -> uint64_t {
        return (uint64_t(a) << 32) | uint64_t(b);
    };

    auto has_match = [](const std::vector<CIGAR::COp>& ops) -> bool {
        return !ops.empty() && CIGAR::match_ratio(CIGAR::pack(ops)) > 0.0;
    };

    auto pick_longest = [&](const std::vector<size_t>& xs, const std::vector<uint64_t>& lens) -> size_t {
        size_t best = xs.front();
        for (size_t x : xs) {
            if (lens[x] > lens[best]) best = x;
        }
        return best;
    };

    // ------------------------------------------------ Structure ------------------------------------------------
    struct ComparedSpan {
        uint32_t v_a = 0;
        uint32_t v_b = 0;
        uint32_t beg_a = UINT32_MAX;
        uint32_t end_a = 0;
        uint32_t beg_b = UINT32_MAX;
        uint32_t end_b = 0;
    };

    // ------------------------------------------------ Initialization ------------------------------------------------
    std::vector<BubbleAlignment> outs;
    path_pair_split::ComparedIndex* cmp = shared_cmp.get();

    if (paths.size() < 2) return outs;

    std::vector<uint64_t> path_lens(paths.size(), 0);
    std::vector<size_t> candidates;
    candidates.reserve(paths.size());

    for (size_t i = 0; i < paths.size(); ++i) {
        if (paths[i].empty()) continue;

        const std::vector<uint32_t> comparable = GfaBubble::GfaBubbleFinder::trim_path_cluster_anchors(paths[i], path_type);
        if (comparable.empty()) continue;

        const bool all_complex = std::all_of(
            comparable.begin(), comparable.end(),
            [this](uint32_t v) { return getNodeComplex(NodeHandle::get_segment_id(v)); }
        );
        if (all_complex) continue;

        path_lens[i] = path_len(comparable);
        candidates.push_back(i);
    }

    if (candidates.size() < 2) return outs;

    // ------------------------------------------------ Alignment ------------------------------------------------
    auto align_pair = [&](size_t ai, size_t bi) -> bool {
        const auto& pathA = paths[ai];
        const auto& pathB = paths[bi];

        if (pathA.empty() || pathB.empty()) return false;
        if (skip_same_start && pathA[0] == pathB[0]) return false;

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

        std::unordered_map<uint64_t, ComparedSpan> compared_spans;
        compared_spans.reserve(tasks.size() * 4 + 8);

        bool produced = false;

        for (const auto& t : tasks) {
            if (t.type == path_pair_split::SplitTask<uint32_t>::Type::ComparedAnchor) {
                continue;
            }

            auto alns = align_subpaths_(pathA, pathB, t.a_span, t.b_span);
            if (alns.empty()) continue;

            for (auto& aln : alns) {
                if (!has_match(aln.ops)) continue;

                const uint64_t key = pair_key(aln.v_a, aln.v_b);
                auto& s = compared_spans[key];

                if (s.beg_a == UINT32_MAX) {
                    s.v_a = aln.v_a;
                    s.v_b = aln.v_b;
                }

                s.beg_a = std::min(s.beg_a, aln.beg_a);
                s.end_a = std::max(s.end_a, aln.end_a);
                s.beg_b = std::min(s.beg_b, aln.beg_b);
                s.end_b = std::max(s.end_b, aln.end_b);

                aln.idx = idx;
                outs.emplace_back(std::move(aln));
                produced = true;
            }
        }

        if (!compared_spans.empty()) {
            std::unique_lock<std::shared_mutex> lock(*shared_cmp_mutex);

            for (const auto& kv : compared_spans) {
                const ComparedSpan& s = kv.second;
                if (s.beg_a >= s.end_a || s.beg_b >= s.end_b) continue;

                const std::string& name_a = getNodeName(NodeHandle::get_segment_id(s.v_a));
                const std::string& name_b = getNodeName(NodeHandle::get_segment_id(s.v_b));
                if (name_a == name_b) continue;

                cmp->add(name_a, s.beg_a, s.end_a, name_b, s.beg_b, s.end_b);
            }
        }

        return produced;
    };

    bool all_paths_short = true;
    for (size_t i : candidates) {
        if (path_lens[i] >= ALL_PAIR_LEN_) {
            all_paths_short = false;
            break;
        }
    }

    // If all paths are short, perform all-vs-all alignment without clustering
    if (all_paths_short) {
        for (size_t a = 0; a + 1 < candidates.size(); ++a) {
            for (size_t b = a + 1; b < candidates.size(); ++b) {
                align_pair(candidates[a], candidates[b]);
            }
        }

        return outs;
    }

    std::unordered_map<uint32_t, std::vector<size_t>> cluster_members;
    std::vector<uint32_t> cluster_order;
    cluster_members.reserve(candidates.size());
    cluster_order.reserve(candidates.size());

    for (size_t i : candidates) {
        const uint32_t cluster = i < clusters.size() ? clusters[i] : static_cast<uint32_t>(i);
        std::pair<std::unordered_map<uint32_t, std::vector<size_t>>::iterator, bool> inserted = cluster_members.emplace(cluster, std::vector<size_t>{});

        if (inserted.second) {
            cluster_order.push_back(cluster);
        }

        inserted.first->second.push_back(i);
    }

    std::vector<size_t> representatives;
    representatives.reserve(cluster_order.size());

    // Align all paths in each cluster to the longest path in the cluster
    for (uint32_t cluster : cluster_order) {
        const std::vector<size_t>& members = cluster_members.at(cluster);
        const size_t ref = pick_longest(members, path_lens);
        representatives.push_back(ref);

        for (size_t member : members) {
            if (member != ref) {
                align_pair(ref, member);
            }
        }
    }

    // Align representatives of different clusters to each other
    for (size_t a = 0; a + 1 < representatives.size(); ++a) {
        for (size_t b = a + 1; b < representatives.size(); ++b) {
            align_pair(representatives[a], representatives[b]);
        }
    }

    return outs;
}

void GfaCollapser::homologous_align_(
    const std::vector<GfaBubble::Bubble>& bubbles, 
    const std::vector<GfaBubble::HomologousPath>& homologous_paths,
    const GfaBubble::GfaBubbleFinder& bubble_finder
) {
    log_stream() << "Aligning bubbles and homologous paths ...\n";

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
    auto consume_alignments = [&](std::vector<BubbleAlignment>&& alignments) {
        for (BubbleAlignment& alignment : alignments) {
            const uint32_t sid_a = NodeHandle::get_segment_id(alignment.v_a);
            const uint32_t sid_b = NodeHandle::get_segment_id(alignment.v_b);

            // If either of the two segments is complex, ignore the alignment result (2026-06-21, 0.1.3-r12)
            if (getNodeComplex(sid_a) || getNodeComplex(sid_b)) {
                continue;
            }

            if (alignment.ops.empty()) continue;

            bubble_aligns_.emplace_back(
                std::move(alignment)
            );
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
    ProgressTracker prog(bubbles.size() + homologous_paths.size());

    // Bubbles
    for (const auto& bb : bubbles) {
        prog.hit();

        const auto paths = bb.get_paths();
        if (paths.size() < 2) continue;

        const auto clusters = bb.get_path_clusters();

        submit_alignment([this, paths, clusters, shared_cmp, shared_cmp_mutex, idx = submit_idx++]() {
            return align_paths_(paths, clusters, GfaBubble::Type::Normal, false, shared_cmp, shared_cmp_mutex, idx);
        });
    }

    // Homologous paths
    for (const auto& hp : homologous_paths) {
        prog.hit();

        const auto paths = hp.get_paths_u32();
        if (paths.size() < 2) continue;

        // Reduce cycle in the graph. (2026-05-25, v0.1.3-r6)
        bool skip_same_start = (hp.type == GfaBubble::Type::DiffSource);

        const auto clusters = hp.get_path_clusters();

        submit_alignment([this, paths, clusters, shared_cmp, shared_cmp_mutex, skip_same_start, type = hp.type, idx = submit_idx++]() {
            return align_paths_(paths, clusters, type, skip_same_start, shared_cmp, shared_cmp_mutex, idx);  // Skip pairs of paths that start at the same node for homologous path alignment
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
    const std::vector<GfaBubble::HomologousPath>& homologous_paths, 
    const std::string& prefix,
    const GfaBubble::GfaBubbleFinder& bubble_finder
) {
    log_stream() << "Transforming graph to de-overlap form and collapsing homologous sequences ...\n" << "\n";

    bubble_aligns_.clear();
    cuts_.clear();
    rulemap_.clear();

    finalize_();
    prune_overlaps_();
    initialize_cuts_();
    overlaps_align_();
    homologous_align_(bubbles, homologous_paths, bubble_finder);
    dedup_aligns_(bubble_aligns_.size());

    const auto align_groups = build_align_groups_();
    build_propagate_prune_cuts_(align_groups);
    if (DEBUG_ENABLED) print_cuts(cuts_);
    build_rulemap_(align_groups);

    SegReplace::Expander ex = build_SegReplace_();
    ex.save_map(prefix + ".collapse.map");
    expand_and_rewire_edges_(ex);
    remove_unused_nodes_();
    merge_linear_chains();

    log_stream() << "Finished transforming graph to de-overlap form and collapsing homologous sequences.\n" << "\n";
}
