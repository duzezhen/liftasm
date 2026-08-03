#include "../include/gfa_collapser.hpp"
#include "../include/gfa_name.hpp"
#include "../include/aligner.hpp"
#include "../include/minimizer_dna.hpp"
#include "../include/ThreadPool.hpp"
#include "../include/progress_tracker.hpp"

#include <algorithm>
#include <future>
#include <limits>
#include <cstdint>
#include <numeric>
#include <unordered_set>
#include <deque>

// ------------------------------------------------ GFA collapser ------------------------------------------------
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
        // const auto hits = align_short_mm2_(
        const auto hits = align_short_wfa_(
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

    std::unordered_map<uint32_t, uint32_t> path_replacements;
    path_replacements.reserve(nodes_merged_num * 2);
    for (size_t i = 0; i < pendings.size(); ++i) {
        for (uint32_t vertex : pendings[i].chain_vtx) {
            path_replacements[vertex] = new_vertices[i];
            path_replacements[vertex ^ 1u] = new_vertices[i] ^ 1u;
        }
    }
    for (GfaPath& path : paths_) {
        std::vector<PathSegment> rewritten;
        rewritten.reserve(path.segments.size());
        for (const PathSegment& segment : path.segments) {
            uint32_t vertex = Vertex::make_vertex(static_cast<uint32_t>(segment.node_id), segment.is_reverse);
            const auto found = path_replacements.find(vertex);
            if (found != path_replacements.end()) vertex = found->second;
            const PathSegment mapped{Vertex::get_segment_id(vertex), Vertex::get_is_reverse(vertex)};
            if (!rewritten.empty() && rewritten.back().node_id == mapped.node_id && rewritten.back().is_reverse == mapped.is_reverse) continue;
            rewritten.push_back(mapped);
        }
        path.segments.swap(rewritten);
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
    path_pair_split::ComparedIndex& compared,
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

        tasks = path_pair_split::split_path_pair(
            pathA, pathB,
            name_of, len_of, rev_of, overlap_of,
            &compared
        );

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
            for (const auto& kv : compared_spans) {
                const ComparedSpan& s = kv.second;
                if (s.beg_a >= s.end_a || s.beg_b >= s.end_b) continue;

                const std::string& name_a = getNodeName(NodeHandle::get_segment_id(s.v_a));
                const std::string& name_b = getNodeName(NodeHandle::get_segment_id(s.v_b));
                if (name_a == name_b) continue;

                compared.add(name_a, s.beg_a, s.end_a, name_b, s.beg_b, s.end_b);
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

std::vector<uint32_t> GfaCollapser::collect_candidate_region_(
    const std::vector<std::vector<uint32_t>>& paths
) const {
    std::unordered_set<uint32_t> region;
    size_t path_vertices = 0;
    for (const auto& path : paths) path_vertices += path.size();
    region.reserve(path_vertices * 2 + 1);

    std::unordered_map<uint32_t, CandidateTopoSpan_> spans;
    spans.reserve(paths.size() + 1);

    for (const auto& path : paths) {
        for (uint32_t vertex : path) {
            const uint32_t sid = NodeHandle::get_segment_id(vertex);
            if (sid >= nodes_.size() || nodes_[sid].deleted) continue;
            region.insert(vertex);

            if (sid >= connectivity_index_.size() || vertex >= topo_index_.size()) continue;
            const uint32_t component = connectivity_index_[sid];
            const uint32_t rank = topo_index_[vertex].topo_rank;
            if (component == UINT32_MAX || rank == UINT32_MAX) continue;

            CandidateTopoSpan_& span = spans[component];
            span.component = component;
            span.begin = std::min(span.begin, rank);
            span.end = std::max(span.end, rank);
        }
    }

    if (spans.empty()) {
        std::vector<uint32_t> result(region.begin(), region.end());
        std::sort(result.begin(), result.end());
        return result;
    }

    std::deque<uint32_t> queue(region.begin(), region.end());
    while (!queue.empty()) {
        const uint32_t vertex = queue.front();
        queue.pop_front();
        if (vertex >= arc_indexs_.size()) continue;

        const uint64_t packed = arc_indexs_[vertex];
        const size_t arc_begin = packed >> 32;
        const uint32_t arc_count = static_cast<uint32_t>(packed);
        for (uint32_t i = 0; i < arc_count; ++i) {
            const GfaArc& arc = arcs_[arc_begin + i];
            if (arc.get_del()) continue;

            const uint32_t next = arc.get_target_vertex_id();
            const uint32_t sid = NodeHandle::get_segment_id(next);
            if (sid >= nodes_.size() || nodes_[sid].deleted) continue;
            if (sid >= connectivity_index_.size() || next >= topo_index_.size()) continue;

            const uint32_t component = connectivity_index_[sid];
            const auto span_it = spans.find(component);
            if (span_it == spans.end()) continue;

            const uint32_t rank = topo_index_[next].topo_rank;
            const CandidateTopoSpan_& span = span_it->second;
            if (rank == UINT32_MAX || rank < span.begin || rank > span.end) continue;

            if (region.insert(next).second) queue.push_back(next);
        }
    }

    std::vector<uint32_t> result(region.begin(), region.end());
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<GfaCollapser::AlignmentCandidate_> GfaCollapser::build_alignment_candidates_(
    const std::vector<GfaBubble::Bubble>& bubbles,
    const std::vector<GfaBubble::HomologousPath>& homologous_paths
) const {
    std::vector<AlignmentCandidate_> candidates;
    candidates.reserve(bubbles.size() + homologous_paths.size());

    uint32_t idx = 0;
    for (size_t source = 0; source < bubbles.size(); ++source) {
        const auto& paths = bubbles[source].get_paths();
        if (paths.size() < 2) continue;

        AlignmentCandidate_ candidate;
        candidate.source_index = source;
        candidate.type = GfaBubble::Type::Normal;
        candidate.idx = idx++;
        candidate.input_paths = static_cast<uint32_t>(paths.size());

        std::unordered_map<uint32_t, CandidateTopoSpan_> spans;
        std::unordered_set<uint32_t> segments;
        for (const auto& path : paths) {
            for (size_t pos = 0; pos < path.size(); ++pos) {
                const uint32_t vertex = path[pos];
                const uint32_t sid = NodeHandle::get_segment_id(vertex);
                if (sid >= nodes_.size() || nodes_[sid].deleted) continue;
                if (segments.insert(sid).second) candidate.priority += nodes_[sid].length;

                // Shared source/sink boundaries must not join adjacent bubbles.
                if (pos == 0 || pos + 1 == path.size()) continue;

                uint32_t component = sid < connectivity_index_.size() ? connectivity_index_[sid] : sid;
                if (component == UINT32_MAX) component = sid;
                const uint32_t rank = vertex < topo_index_.size() && topo_index_[vertex].topo_rank != UINT32_MAX ? topo_index_[vertex].topo_rank : sid;
                auto& span = spans[component];
                span.component = component;
                span.begin = std::min(span.begin, rank);
                span.end = std::max(span.end, rank);
            }
        }
        for (const auto& item : spans) candidate.spans.push_back(item.second);
        candidate.region_vertices = collect_candidate_region_(paths);
        candidates.push_back(std::move(candidate));
    }

    for (size_t source = 0; source < homologous_paths.size(); ++source) {
        const auto& paths = homologous_paths[source].get_paths();
        if (paths.size() < 2) continue;

        AlignmentCandidate_ candidate;
        candidate.homologous = true;
        candidate.source_index = source;
        candidate.type = homologous_paths[source].type;
        candidate.skip_same_start = candidate.type == GfaBubble::Type::DiffSource;
        candidate.idx = idx++;
        candidate.input_paths = static_cast<uint32_t>(paths.size());

        std::unordered_map<uint32_t, CandidateTopoSpan_> spans;
        std::unordered_set<uint32_t> segments;
        for (const auto& path : paths) {
            for (size_t pos = 0; pos < path.size(); ++pos) {
                const Vertex vertex = path[pos];
                const uint32_t sid = vertex.segment_id();
                if (sid >= nodes_.size() || nodes_[sid].deleted) continue;
                if (segments.insert(sid).second) candidate.priority += nodes_[sid].length;

                // SameSource paths share the first node. DiffSource paths have
                // no shared boundary, so all their nodes define conflicts.
                if (candidate.type == GfaBubble::Type::SameSource && pos == 0) continue;

                uint32_t component = sid < connectivity_index_.size() ? connectivity_index_[sid] : sid;
                if (component == UINT32_MAX) component = sid;
                const uint32_t vid = vertex.vertex_id();
                const uint32_t rank = vid < topo_index_.size() && topo_index_[vid].topo_rank != UINT32_MAX ? topo_index_[vid].topo_rank : sid;
                auto& span = spans[component];
                span.component = component;
                span.begin = std::min(span.begin, rank);
                span.end = std::max(span.end, rank);
            }
        }
        for (const auto& item : spans) candidate.spans.push_back(item.second);
        candidate.region_vertices = collect_candidate_region_(
            homologous_paths[source].get_paths_u32()
        );
        candidates.push_back(std::move(candidate));
    }

    return candidates;
}

std::vector<std::vector<size_t>> GfaCollapser::group_alignment_candidates_(
    const std::vector<AlignmentCandidate_>& candidates
) const {
    struct SpanRef {
        uint32_t component;
        uint32_t begin;
        uint32_t end;
        size_t candidate;
    };

    std::vector<size_t> parent(candidates.size());
    std::vector<uint8_t> rank(candidates.size(), 0);
    std::iota(parent.begin(), parent.end(), 0);

    auto root = [&parent](size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto join = [&parent, &rank, &root](size_t a, size_t b) {
        a = root(a);
        b = root(b);
        if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
    };

    std::vector<SpanRef> spans;
    for (size_t i = 0; i < candidates.size(); ++i) {
        for (const CandidateTopoSpan_& span : candidates[i].spans) {
            spans.push_back({span.component, span.begin, span.end, i});
        }
    }
    std::sort(spans.begin(), spans.end(), [](const SpanRef& a, const SpanRef& b) {
        if (a.component != b.component) return a.component < b.component;
        if (a.begin != b.begin) return a.begin < b.begin;
        if (a.end != b.end) return a.end > b.end;
        return a.candidate < b.candidate;
    });

    size_t i = 0;
    while (i < spans.size()) {
        size_t representative = spans[i].candidate;
        uint32_t component = spans[i].component;
        uint32_t covered_end = spans[i].end;
        ++i;

        while (i < spans.size() && spans[i].component == component) {
            if (spans[i].begin > covered_end) {
                representative = spans[i].candidate;
                covered_end = spans[i].end;
            } else {
                join(representative, spans[i].candidate);
                covered_end = std::max(covered_end, spans[i].end);
            }
            ++i;
        }
    }

    std::unordered_map<size_t, size_t> group_index;
    std::vector<std::vector<size_t>> groups;
    for (size_t candidate = 0; candidate < candidates.size(); ++candidate) {
        const size_t r = root(candidate);
        auto inserted = group_index.emplace(r, groups.size());
        if (inserted.second) groups.emplace_back();
        groups[inserted.first->second].push_back(candidate);
    }

    for (auto& group : groups) {
        std::sort(group.begin(), group.end(), [&candidates](size_t a, size_t b) {
            if (candidates[a].priority != candidates[b].priority) {
                return candidates[a].priority > candidates[b].priority;
            }
            if (candidates[a].homologous != candidates[b].homologous) {
                return !candidates[a].homologous;
            }
            return candidates[a].idx < candidates[b].idx;
        });
    }
    std::sort(groups.begin(), groups.end(), [&candidates](const auto& a, const auto& b) {
        uint32_t ai = UINT32_MAX;
        uint32_t bi = UINT32_MAX;
        for (size_t x : a) ai = std::min(ai, candidates[x].idx);
        for (size_t x : b) bi = std::min(bi, candidates[x].idx);
        return ai < bi;
    });

    return groups;
}

std::vector<GfaCollapser::BubbleAlignment> GfaCollapser::align_candidate_group_(
    const std::vector<size_t>& group,
    const std::vector<AlignmentCandidate_>& candidates,
    const std::vector<GfaBubble::Bubble>& bubbles,
    const std::vector<GfaBubble::HomologousPath>& homologous_paths
) {
    path_pair_split::ComparedIndex compared;
    std::vector<BubbleAlignment> result;

    for (size_t candidate_id : group) {
        const AlignmentCandidate_& candidate = candidates[candidate_id];
        std::vector<BubbleAlignment> alignments;

        if (candidate.homologous) {
            const auto& source = homologous_paths[candidate.source_index];
            alignments = align_paths_(
                source.get_paths_u32(), source.get_path_clusters(), candidate.type,
                candidate.skip_same_start, compared, candidate.idx
            );
        } else {
            const auto& source = bubbles[candidate.source_index];
            alignments = align_paths_(
                source.get_paths(), source.get_path_clusters(), candidate.type,
                false, compared, candidate.idx
            );
        }

        result.insert(
            result.end(),
            std::make_move_iterator(alignments.begin()),
            std::make_move_iterator(alignments.end())
        );
    }
    return result;
}

void GfaCollapser::reset_collapse_cuts_() {
    cuts_.assign(nodes_.size(), Cuts{});
    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        cuts_[sid].name = nodes_[sid].name;
        cuts_[sid].v = {0, nodes_[sid].length};
    }

    for (const GfaArc& arc : arcs_) {
        if (arc.get_del() || arc.get_comp()) continue;

        const uint32_t v_sid = arc.get_source_segment_id();
        const uint32_t w_sid = arc.get_target_segment_id();
        if (v_sid >= nodes_.size() || w_sid >= nodes_.size()) continue;

        const uint32_t v_len = nodes_[v_sid].length;
        const uint32_t w_len = nodes_[w_sid].length;
        const uint32_t ov = arc.ov == INT32_MAX ? 0 : static_cast<uint32_t>(std::clamp<int64_t>(arc.ov, 0, v_len));
        const uint32_t ow = arc.ow == INT32_MAX ? 0 : static_cast<uint32_t>(std::clamp<int64_t>(arc.ow, 0, w_len));
        if (std::min(ov, ow) == 0) continue;

        if (arc.get_source_is_reverse()) {
            cuts_[v_sid].v.push_back(0);
            cuts_[v_sid].v.push_back(ov);
        } else {
            cuts_[v_sid].v.push_back(v_len - ov);
            cuts_[v_sid].v.push_back(v_len);
        }

        if (arc.get_target_is_reverse()) {
            cuts_[w_sid].v.push_back(w_len - ow);
            cuts_[w_sid].v.push_back(w_len);
        } else {
            cuts_[w_sid].v.push_back(0);
            cuts_[w_sid].v.push_back(ow);
        }
    }
}

SegReplace::Expander GfaCollapser::build_safe_expander_(size_t alignment_begin) {
    log_stream() << "Removing bubbles and homologous paths that increase local graph complexity ...\n\n";

    size_t round = 0;
    while (true) {
        ++round;

        reset_collapse_cuts_();
        const auto align_groups = build_align_groups_();
        build_propagate_prune_cuts_(align_groups);
        if (DEBUG_ENABLED) print_cuts(cuts_);
        build_rulemap_(align_groups);

        SegReplace::Expander expander = build_SegReplace_();
        LocalGraphComplexity checker(*this, alignment_begin);
        const std::unordered_set<int32_t> bad = checker.find_increasing_candidates();
        if (bad.empty()) {
            log_stream() << "Stable after round " << round << ": no complexity-increasing candidates\n\n";
            return expander;
        }

        const size_t removed = checker.reject_candidates(bad);

        log_stream() << "Round " << round << ": removed " << bad.size() << " complexity-increasing candidates (" << removed << " alignments)\n\n";
        if (removed == 0) {
            log_stream() << "Stopped after round " << round << ": no matching alignment could be removed\n\n";
            return expander;
        }
    }
}

void GfaCollapser::homologous_align_(
    const std::vector<GfaBubble::Bubble>& bubbles,
    const std::vector<GfaBubble::HomologousPath>& homologous_paths
) {
    log_stream() << "Aligning bubbles and homologous paths ...\n";

    alignment_candidates_ = build_alignment_candidates_(bubbles, homologous_paths);
    const std::vector<AlignmentCandidate_>& candidates = alignment_candidates_;
    const std::vector<std::vector<size_t>> groups = group_alignment_candidates_(candidates);

    log_stream() << "  - Alignment candidates: " << candidates.size() << "\n";
    log_stream() << "  - Topological conflict groups: " << groups.size() << "\n";

    ThreadPool pool(alignOpts_.threads);
    const size_t max_inflight = std::max<size_t>(1, alignOpts_.threads * 2);
    std::deque<std::future<std::vector<BubbleAlignment>>> futures;
    ProgressTracker progress(candidates.size());
    size_t submitted = 0;

    while (submitted < groups.size() || !futures.empty()) {
        while (submitted < groups.size() && futures.size() < max_inflight) {
            const size_t group_id = submitted++;
            futures.emplace_back(pool.submit(
                [this, group_id, &groups, &candidates, &bubbles, &homologous_paths, &progress]() {
                    std::vector<BubbleAlignment> result = align_candidate_group_(
                        groups[group_id], candidates, bubbles, homologous_paths
                    );
                    progress.add(groups[group_id].size());
                    return result;
                }
            ));
        }

        std::vector<BubbleAlignment> alignments = futures.front().get();
        futures.pop_front();

        for (BubbleAlignment& alignment : alignments) {
            const uint32_t sid_a = NodeHandle::get_segment_id(alignment.v_a);
            const uint32_t sid_b = NodeHandle::get_segment_id(alignment.v_b);
            if (getNodeComplex(sid_a) || getNodeComplex(sid_b) || alignment.ops.empty()) continue;
            bubble_aligns_.emplace_back(std::move(alignment));
        }
    }
    progress.finish();
    pool.stop();

    log_stream() << "  - Total alignments produced: " << bubble_aligns_.size() << "\n\n";
}

void GfaCollapser::collapse_homologous_seq(
    const std::vector<GfaBubble::Bubble>& bubbles,
    const std::vector<GfaBubble::HomologousPath>& homologous_paths, 
    const std::string& prefix
) {
    log_stream() << "Transforming graph to de-overlap form and collapsing homologous sequences ...\n" << "\n";

    bubble_aligns_.clear();
    cuts_.clear();
    rulemap_.clear();

    finalize_();
    prune_overlaps_();
    initialize_cuts_();
    overlaps_align_();
    const size_t homologous_align_begin = bubble_aligns_.size();
    homologous_align_(bubbles, homologous_paths);
    dedup_aligns_(homologous_align_begin);

    SegReplace::Expander ex = build_safe_expander_(homologous_align_begin);
    GfaPathCycleGuard(*this).filter_path_rules(ex);
    ex.save_map(prefix + ".collapse.map");
    expand_and_rewire_edges_(ex);
    rewrite_paths_(ex);
    remove_unused_nodes_();
    merge_linear_chains();

    log_stream() << "Finished transforming graph to de-overlap form and collapsing homologous sequences.\n" << "\n";
}


/* ================================================================================================================
 *                                          LOCAL GRAPH COMPLEXITY START
 * ================================================================================================================ */
LocalGraphComplexity::LocalGraphComplexity(
    GfaCollapser& collapser,
    size_t alignment_begin
)
    : collapser_(collapser),
      alignment_begin_(alignment_begin) {}

size_t LocalGraphComplexity::reject_candidates(
    const std::unordered_set<int32_t>& candidates
) {
    if (collapser_.mark_rejected_complex_) {
        std::vector<uint32_t> segments;
        for (size_t i = alignment_begin_; i < collapser_.bubble_aligns_.size(); ++i) {
            const auto& alignment = collapser_.bubble_aligns_[i];
            if (!candidates.contains(alignment.idx)) continue;
            segments.push_back(NodeHandle::get_segment_id(alignment.v_a));
            segments.push_back(NodeHandle::get_segment_id(alignment.v_b));
        }
        std::sort(segments.begin(), segments.end());
        segments.erase(std::unique(segments.begin(), segments.end()), segments.end());
        collapser_.mark_segments_complex(segments);
    }

    const size_t old_size = collapser_.bubble_aligns_.size();
    collapser_.bubble_aligns_.erase(
        std::remove_if(
            collapser_.bubble_aligns_.begin() + alignment_begin_,
            collapser_.bubble_aligns_.end(),
            [this, &candidates](const auto& alignment) {
                if (candidates.contains(alignment.idx)) return true;
                const uint32_t sid_a = NodeHandle::get_segment_id(alignment.v_a);
                const uint32_t sid_b = NodeHandle::get_segment_id(alignment.v_b);
                return collapser_.getNodeComplex(sid_a) || collapser_.getNodeComplex(sid_b);
            }
        ),
        collapser_.bubble_aligns_.end()
    );
    return old_size - collapser_.bubble_aligns_.size();
}

std::unordered_set<int32_t> LocalGraphComplexity::find_increasing_candidates() const {
    std::unordered_set<int32_t> active;
    for (size_t i = alignment_begin_; i < collapser_.bubble_aligns_.size(); ++i) {
        active.insert(collapser_.bubble_aligns_[i].idx);
    }

    std::vector<size_t> work;
    work.reserve(active.size());
    for (size_t i = 0; i < collapser_.alignment_candidates_.size(); ++i) {
        if (active.contains(static_cast<int32_t>(collapser_.alignment_candidates_[i].idx))) {
            work.push_back(i);
        }
    }

    std::unordered_set<int32_t> bad;
    if (work.empty()) return bad;

    std::vector<std::vector<uint32_t>> affected_segments(collapser_.alignment_candidates_.size());
    for (size_t i = alignment_begin_; i < collapser_.bubble_aligns_.size(); ++i) {
        const auto& alignment = collapser_.bubble_aligns_[i];
        if (alignment.idx < 0 || static_cast<size_t>(alignment.idx) >= affected_segments.size()) {
            continue;
        }
        affected_segments[alignment.idx].push_back(NodeHandle::get_segment_id(alignment.v_a));
        affected_segments[alignment.idx].push_back(NodeHandle::get_segment_id(alignment.v_b));
    }
    for (std::vector<uint32_t>& segments : affected_segments) {
        std::sort(segments.begin(), segments.end());
        segments.erase(std::unique(segments.begin(), segments.end()), segments.end());
    }

    using Rule = std::pair<SegReplace::Seg, SegReplace::Expansion>;
    std::unordered_map<uint32_t, std::vector<Rule>> rules_by_segment;
    rules_by_segment.reserve(std::min(collapser_.rulemap_.size(), collapser_.nodes_.size()));
    for (const auto& rule : collapser_.rulemap_) {
        const uint32_t sid = static_cast<uint32_t>(SegReplace::Interval::seg_id(rule.first));
        rules_by_segment[sid].push_back(rule);
    }

    ThreadPool pool(collapser_.alignOpts_.threads);
    const size_t max_inflight = std::max<size_t>(1, collapser_.alignOpts_.threads * 2);
    std::deque<std::future<std::pair<int32_t, bool>>> futures;
    size_t submitted = 0;

    while (submitted < work.size() || !futures.empty()) {
        while (submitted < work.size() && futures.size() < max_inflight) {
            const size_t candidate_id = work[submitted++];
            futures.emplace_back(pool.submit(
                [this, candidate_id, &affected_segments, &rules_by_segment]() {
                    const auto& candidate = collapser_.alignment_candidates_[candidate_id];
                    SegReplace::RuleMap local_rules;
                    for (uint32_t sid : affected_segments[candidate.idx]) {
                        const auto found = rules_by_segment.find(sid);
                        if (found == rules_by_segment.end()) continue;
                        for (const auto& rule : found->second) {
                            local_rules.emplace(rule.first, rule.second);
                        }
                    }

                    SegReplace::Expander expander(
                        local_rules,
                        {},
                        collapser_.MIN_TRANS_LEN_
                    );
                    expander.build_index();
                    expander.filter_nonmonotonic_index();
                    return std::make_pair(
                        static_cast<int32_t>(candidate.idx),
                        candidate_increases_(candidate_id, expander)
                    );
                }
            ));
        }

        const auto result = futures.front().get();
        futures.pop_front();
        if (result.second) bad.insert(result.first);
    }

    pool.stop();
    return bad;
}

LocalGraphComplexity::LocalGraph LocalGraphComplexity::build_local_graph_(
    const std::vector<uint32_t>& vertices,
    const SegReplace::Expander* expander
) const {
    using Seg = SegReplace::Seg;

    std::unordered_map<Seg, uint32_t, SegReplace::U128Hash, SegReplace::U128Eq> node_ids;
    std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> endpoints;
    std::unordered_set<uint32_t> region(vertices.begin(), vertices.end());
    std::unordered_set<uint64_t> directed_edges;
    node_ids.reserve(vertices.size() * 2 + 8);
    endpoints.reserve(vertices.size());
    directed_edges.reserve(vertices.size() * 4 + 8);

    for (uint32_t vertex : vertices) {
        const uint32_t sid = NodeHandle::get_segment_id(vertex);
        const bool reverse = NodeHandle::get_is_reverse(vertex);
        if (sid >= collapser_.nodes_.size() || collapser_.nodes_[sid].deleted) continue;

        std::vector<uint32_t> fallback;
        const std::vector<uint32_t>* boundaries = nullptr;
        if (sid < collapser_.cuts_.size() && collapser_.cuts_[sid].v.size() >= 2) {
            boundaries = &collapser_.cuts_[sid].v;
        } else {
            fallback = {0, collapser_.nodes_[sid].length};
            boundaries = &fallback;
        }

        uint32_t first = UINT32_MAX;
        uint32_t previous = UINT32_MAX;
        for (size_t step = 0; step + 1 < boundaries->size(); ++step) {
            const size_t piece_index = reverse ? boundaries->size() - step - 2 : step;
            const Seg piece = SegReplace::Interval::pack(
                sid,
                (*boundaries)[piece_index],
                (*boundaries)[piece_index + 1],
                reverse
            );
            const SegReplace::Expansion expansion = expander ? expander->query(piece) : SegReplace::Expansion{piece};

            for (Seg mapped : expansion) {
                mapped = SegReplace::Interval::canonical(mapped);
                const auto inserted = node_ids.emplace(
                    mapped,
                    static_cast<uint32_t>(node_ids.size())
                );
                const uint32_t current = inserted.first->second;
                if (current == previous) continue;
                if (first == UINT32_MAX) first = current;
                if (previous != UINT32_MAX) {
                    directed_edges.insert((uint64_t(previous) << 32) | current);
                }
                previous = current;
            }
        }

        if (first != UINT32_MAX) {
            endpoints.emplace(vertex, std::make_pair(first, previous));
        }
    }

    for (uint32_t vertex : vertices) {
        const auto from = endpoints.find(vertex);
        if (from == endpoints.end() || vertex >= collapser_.arc_indexs_.size()) continue;

        const uint64_t packed = collapser_.arc_indexs_[vertex];
        const size_t arc_begin = packed >> 32;
        const uint32_t arc_count = static_cast<uint32_t>(packed);
        for (uint32_t i = 0; i < arc_count; ++i) {
            const GfaArc& arc = collapser_.arcs_[arc_begin + i];
            if (arc.get_del()) continue;

            const uint32_t target = arc.get_target_vertex_id();
            if (!region.contains(target)) continue;

            const auto to = endpoints.find(target);
            if (to == endpoints.end()) continue;
            directed_edges.insert(
                (uint64_t(from->second.second) << 32) | to->second.first
            );
        }
    }

    return {
        static_cast<uint32_t>(node_ids.size()),
        std::move(directed_edges)
    };
}

bool LocalGraphComplexity::candidate_increases_(
    size_t candidate_id,
    const SegReplace::Expander& expander
) const {
    const auto& candidate = collapser_.alignment_candidates_[candidate_id];
    if (candidate.region_vertices.empty() || candidate.spans.size() > 1) return false;

    const LocalGraph before = build_local_graph_(candidate.region_vertices, nullptr);
    const LocalGraph after = build_local_graph_(candidate.region_vertices, &expander);
    const bool increases = increases_(candidate.input_paths, before, after);
    if (increases && DEBUG_ENABLED) {
        const Metrics old_metrics = measure_(
            before.node_count,
            before.directed_edges
        );
        const Metrics new_metrics = measure_(
            after.node_count,
            after.directed_edges
        );
        debug_stream()
            << "  - Complexity candidate " << candidate.idx
            << ": paths=" << candidate.input_paths
            << ", nodes=" << before.node_count << "->" << after.node_count
            << ", block_cycles=" << old_metrics.max_block_cycle_excess
            << "->" << new_metrics.max_block_cycle_excess
            << ", block_branches=" << old_metrics.max_block_branch_excess
            << "->" << new_metrics.max_block_branch_excess << "\n";
    }
    return increases;
}

bool LocalGraphComplexity::increases_(
    uint32_t input_paths,
    const LocalGraph& before_graph,
    const LocalGraph& after_graph
) const {
    const Metrics before = measure_(
        before_graph.node_count,
        before_graph.directed_edges
    );
    const Metrics after = measure_(
        after_graph.node_count,
        after_graph.directed_edges
    );

    const uint64_t path_budget = uint64_t(std::max<uint32_t>(1, input_paths)) * MAX_COMPLEXITY_GROWTH_;
    const uint64_t allowed_block_cycles = std::max<uint64_t>(before.max_block_cycle_excess * MAX_COMPLEXITY_GROWTH_, path_budget);
    const uint64_t allowed_block_branches = std::max(before.max_block_branch_excess * MAX_COMPLEXITY_GROWTH_, path_budget);

    return after.max_block_cycle_excess > allowed_block_cycles && after.max_block_branch_excess > allowed_block_branches;
}

LocalGraphComplexity::Metrics LocalGraphComplexity::measure_(
    uint32_t node_count,
    const std::unordered_set<uint64_t>& directed_edges
) const {
    Metrics result;
    if (node_count == 0) return result;

    std::vector<std::pair<uint32_t, uint32_t>> edges;
    edges.reserve(directed_edges.size());
    for (uint64_t edge : directed_edges) {
        uint32_t from = static_cast<uint32_t>(edge >> 32);
        uint32_t to = static_cast<uint32_t>(edge);
        if (from == to) continue;
        if (from > to) std::swap(from, to);
        edges.emplace_back(from, to);
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    std::vector<std::vector<uint32_t>> adjacency(node_count);
    for (uint32_t edge_id = 0; edge_id < edges.size(); ++edge_id) {
        const auto [from, to] = edges[edge_id];
        adjacency[from].push_back(edge_id);
        adjacency[to].push_back(edge_id);
    }

    struct DfsFrame {
        uint32_t node;
        size_t next_edge;
    };

    std::vector<uint32_t> discovery(node_count, 0);
    std::vector<uint32_t> low(node_count, 0);
    std::vector<uint32_t> parent_edge(node_count, UINT32_MAX);
    std::vector<uint32_t> edge_stack;
    std::vector<uint32_t> block_edges;
    std::vector<DfsFrame> frames;
    edge_stack.reserve(edges.size());
    block_edges.reserve(64);
    frames.reserve(node_count);
    uint32_t timer = 0;

    for (uint32_t start = 0; start < node_count; ++start) {
        if (discovery[start] != 0) continue;
        discovery[start] = low[start] = ++timer;
        frames.push_back({start, 0});

        while (!frames.empty()) {
            DfsFrame& frame = frames.back();
            const uint32_t node = frame.node;
            if (frame.next_edge < adjacency[node].size()) {
                const uint32_t edge_id = adjacency[node][frame.next_edge++];
                if (edge_id == parent_edge[node]) continue;

                const auto [left, right] = edges[edge_id];
                const uint32_t target = left == node ? right : left;
                if (discovery[target] == 0) {
                    parent_edge[target] = edge_id;
                    discovery[target] = low[target] = ++timer;
                    edge_stack.push_back(edge_id);
                    frames.push_back({target, 0});
                } else if (discovery[target] < discovery[node]) {
                    low[node] = std::min(low[node], discovery[target]);
                    edge_stack.push_back(edge_id);
                }
                continue;
            }

            frames.pop_back();
            const uint32_t tree_edge = parent_edge[node];
            if (tree_edge == UINT32_MAX) continue;

            const auto [left, right] = edges[tree_edge];
            const uint32_t parent = left == node ? right : left;
            low[parent] = std::min(low[parent], low[node]);
            if (low[node] < discovery[parent]) continue;

            block_edges.clear();
            while (!edge_stack.empty()) {
                const uint32_t edge_id = edge_stack.back();
                edge_stack.pop_back();
                block_edges.push_back(edge_id);
                if (edge_id == tree_edge) break;
            }
            measure_block_(edges, block_edges, result);
        }
    }

    return result;
}

void LocalGraphComplexity::measure_block_(
    const std::vector<std::pair<uint32_t, uint32_t>>& edges,
    const std::vector<uint32_t>& block_edges,
    Metrics& metrics
) {
    if (block_edges.empty()) return;

    std::vector<uint32_t> vertices;
    vertices.reserve(block_edges.size() * 2);
    for (uint32_t edge_id : block_edges) {
        vertices.push_back(edges[edge_id].first);
        vertices.push_back(edges[edge_id].second);
    }
    std::sort(vertices.begin(), vertices.end());

    uint64_t branch_excess = 0;
    uint64_t vertex_count = 0;
    for (size_t begin = 0; begin < vertices.size();) {
        size_t end = begin + 1;
        while (end < vertices.size() && vertices[end] == vertices[begin]) ++end;
        const uint64_t degree = end - begin;
        if (degree > 2) branch_excess += degree - 2;
        ++vertex_count;
        begin = end;
    }

    const uint64_t cycle_excess = block_edges.size() >= vertex_count ? block_edges.size() - vertex_count + 1 : 0;
    metrics.max_block_cycle_excess = std::max(metrics.max_block_cycle_excess, cycle_excess);
    metrics.max_block_branch_excess = std::max(metrics.max_block_branch_excess, branch_excess);
}
/* ================================================================================================================
 *                                          LOCAL GRAPH COMPLEXITY END
 * ================================================================================================================ */

/* ================================================================================================================
 *                                         PATH CYCLE GUARD START
 * ================================================================================================================ */
GfaPathCycleGuard::GfaPathCycleGuard(GfaCollapser& graph)
    : graph_(graph) {}

size_t GfaPathCycleGuard::filter_rules(
    const std::vector<GfaPath>& source_paths,
    SegReplace::Expander& expander
) const {
    std::vector<const GfaPath*> paths;
    paths.reserve(source_paths.size());
    for (const GfaPath& path : source_paths) {
        if (!path.segments.empty()) paths.push_back(&path);
    }
    return filter_(paths, expander);
}

size_t GfaPathCycleGuard::filter_path_rules(SegReplace::Expander& expander) const {
    return filter_rules(graph_.paths_, expander);
}

size_t GfaPathCycleGuard::filter_(
    const std::vector<const GfaPath*>& paths,
    SegReplace::Expander& expander
) const {
    log_stream() << "Checking paths for replacement cycles ...\n";
    if (paths.empty()) {
        log_stream() << "  - Paths checked: 0\n";
        log_stream() << "  - Replacement rules reverted: 0\n\n";
        return 0;
    }

    size_t total_removed = 0;
    size_t round = 0;
    while (true) {
        ++round;
        const auto bad = find_bad_rules_(paths, expander);
        if (bad.empty()) {
            log_stream() << "  - Round " << round << " reverted rules: 0\n";
            break;
        }

        for (SegReplace::Seg rule : bad) expander.remove_from_index_by_key(rule);
        total_removed += bad.size();
        log_stream() << "  - Round " << round << " reverted rules: " << bad.size() << "\n";
    }

    log_stream() << "  - Paths checked: " << paths.size() << "\n";
    log_stream() << "  - Replacement rules reverted: " << total_removed << "\n\n";
    return total_removed;
}

std::unordered_set<SegReplace::Seg, SegReplace::U128Hash, SegReplace::U128Eq> GfaPathCycleGuard::find_bad_rules_(
    const std::vector<const GfaPath*>& paths,
    const SegReplace::Expander& expander
) const {
    ThreadPool pool(graph_.alignOpts_.threads);
    const size_t max_inflight = std::max<size_t>(1, graph_.alignOpts_.threads * 2);
    std::deque<std::future<PathScan>> futures;
    std::unordered_set<SegReplace::Seg, SegReplace::U128Hash, SegReplace::U128Eq> bad;
    std::unordered_map<SegReplace::Seg, uint32_t, SegReplace::U128Hash, SegReplace::U128Eq> node_ids;
    std::unordered_set<uint64_t> edge_set;
    std::unordered_map<uint32_t, std::vector<SegReplace::Seg>> node_rules;
    ProgressTracker progress(paths.size());
    size_t submitted = 0;

    while (submitted < paths.size() || !futures.empty()) {
        while (submitted < paths.size() && futures.size() < max_inflight) {
            const GfaPath* path = paths[submitted++];
            futures.emplace_back(pool.submit([this, path, &expander, &progress]() {
                PathScan result = scan_path_(*path, expander);
                progress.hit();
                return result;
            }));
        }

        PathScan scan = futures.front().get();
        futures.pop_front();
        bad.insert(scan.bad_rules.begin(), scan.bad_rules.end());

        std::vector<uint32_t> path_nodes;
        path_nodes.reserve(scan.steps.size());
        for (const ExpandedStep& step : scan.steps) {
            const auto inserted = node_ids.emplace(step.interval, static_cast<uint32_t>(node_ids.size()));
            path_nodes.push_back(inserted.first->second);
            if (step.replaced) node_rules[inserted.first->second].push_back(step.owner);
        }
        for (size_t i = 1; i < path_nodes.size(); ++i) {
            const uint32_t from = path_nodes[i - 1];
            const uint32_t to = path_nodes[i];
            if (from != to) edge_set.insert((uint64_t(from) << 32) | to);
        }
    }
    progress.finish();
    pool.stop();

    const std::vector<uint64_t> edges(edge_set.begin(), edge_set.end());
    const CycleIndex cycles = find_cycles_(static_cast<uint32_t>(node_ids.size()), edges);
    for (uint32_t node = 0; node < cycles.component.size(); ++node) {
        if (!cycles.cyclic[cycles.component[node]]) continue;
        const auto found = node_rules.find(node);
        if (found != node_rules.end()) bad.insert(found->second.begin(), found->second.end());
    }
    return bad;
}

GfaPathCycleGuard::PathScan GfaPathCycleGuard::scan_path_(
    const GfaPath& path,
    const SegReplace::Expander& expander
) const {
    struct LastInterval {
        SegReplace::Seg owner{0};
        uint32_t beg{0};
        uint32_t end{0};
        size_t step{0};
        bool reverse{false};
        bool replaced{false};
    };

    std::unordered_map<uint64_t, LastInterval> last;
    PathScan result;
    std::vector<std::pair<size_t, size_t>> bad_ranges;

    for (const PathSegment& segment : path.segments) {
        const uint32_t sid = static_cast<uint32_t>(segment.node_id);
        if (sid >= graph_.nodes_.size() || graph_.nodes_[sid].deleted) continue;

        const std::vector<uint32_t>* boundaries = nullptr;
        std::vector<uint32_t> fallback;
        if (sid < graph_.cuts_.size() && graph_.cuts_[sid].v.size() >= 2) {
            boundaries = &graph_.cuts_[sid].v;
        } else {
            fallback = {0, graph_.nodes_[sid].length};
            boundaries = &fallback;
        }

        for (size_t step = 0; step + 1 < boundaries->size(); ++step) {
            const size_t i = segment.is_reverse ? boundaries->size() - step - 2 : step;
            const SegReplace::Seg owner = SegReplace::Interval::pack(
                sid, (*boundaries)[i], (*boundaries)[i + 1], segment.is_reverse
            );
            const SegReplace::Expansion expansion = expander.query(owner);
            const bool replaced = expansion.size() != 1 || expansion.front() != owner;

            for (SegReplace::Seg interval : expansion) {
                const size_t current_step = result.steps.size();
                result.steps.push_back({interval, owner, replaced});
                const uint64_t target = SegReplace::Interval::seg_id(interval);
                const uint32_t beg = SegReplace::Interval::beg(interval);
                const uint32_t end = SegReplace::Interval::end(interval);
                const bool reverse = SegReplace::Interval::is_reverse(interval);
                auto found = last.find(target);
                if (found == last.end()) {
                    last.emplace(target, LastInterval{owner, beg, end, current_step, reverse, replaced});
                    continue;
                }

                const LastInterval& previous = found->second;
                const bool conflict = reverse != previous.reverse || (!reverse && beg < previous.end) || (reverse && end > previous.beg);
                if (!conflict) {
                    found->second = {owner, beg, end, current_step, reverse, replaced};
                    continue;
                }

                if (previous.replaced || replaced) {
                    bad_ranges.emplace_back(previous.step, current_step);
                }
                if (!previous.replaced || replaced) continue;
                found->second = {owner, beg, end, current_step, reverse, false};
            }
        }
    }

    std::sort(bad_ranges.begin(), bad_ranges.end());
    size_t range = 0;
    while (range < bad_ranges.size()) {
        size_t begin = bad_ranges[range].first;
        size_t end = bad_ranges[range].second;
        while (++range < bad_ranges.size() && bad_ranges[range].first <= end + 1) {
            end = std::max(end, bad_ranges[range].second);
        }
        for (size_t i = begin; i <= end; ++i) {
            if (result.steps[i].replaced) result.bad_rules.push_back(result.steps[i].owner);
        }
    }

    std::sort(result.bad_rules.begin(), result.bad_rules.end());
    result.bad_rules.erase(
        std::unique(result.bad_rules.begin(), result.bad_rules.end()),
        result.bad_rules.end()
    );
    return result;
}

GfaPathCycleGuard::CycleIndex GfaPathCycleGuard::find_cycles_(
    uint32_t node_count,
    const std::vector<uint64_t>& edges
) {
    CycleIndex result;
    result.component.assign(node_count, UINT32_MAX);
    if (node_count == 0) return result;

    std::vector<uint64_t> out_offset(node_count + 1, 0);
    std::vector<uint64_t> in_offset(node_count + 1, 0);
    for (uint64_t edge : edges) {
        ++out_offset[static_cast<uint32_t>(edge >> 32) + 1];
        ++in_offset[static_cast<uint32_t>(edge) + 1];
    }
    for (uint32_t i = 1; i <= node_count; ++i) {
        out_offset[i] += out_offset[i - 1];
        in_offset[i] += in_offset[i - 1];
    }

    std::vector<uint32_t> out(edges.size());
    std::vector<uint32_t> in(edges.size());
    std::vector<uint64_t> out_next = out_offset;
    std::vector<uint64_t> in_next = in_offset;
    for (uint64_t edge : edges) {
        const uint32_t from = static_cast<uint32_t>(edge >> 32);
        const uint32_t to = static_cast<uint32_t>(edge);
        out[out_next[from]++] = to;
        in[in_next[to]++] = from;
    }

    struct Frame {
        uint32_t node;
        uint64_t next;
    };
    std::vector<uint8_t> seen(node_count, 0);
    std::vector<uint32_t> order;
    std::vector<Frame> dfs;
    order.reserve(node_count);
    dfs.reserve(256);
    for (uint32_t start = 0; start < node_count; ++start) {
        if (seen[start]) continue;
        seen[start] = 1;
        dfs.push_back({start, out_offset[start]});
        while (!dfs.empty()) {
            Frame& frame = dfs.back();
            if (frame.next < out_offset[frame.node + 1]) {
                const uint32_t next = out[frame.next++];
                if (!seen[next]) {
                    seen[next] = 1;
                    dfs.push_back({next, out_offset[next]});
                }
                continue;
            }
            order.push_back(frame.node);
            dfs.pop_back();
        }
    }

    std::vector<uint32_t> stack;
    std::vector<uint32_t> component_size;
    stack.reserve(256);
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const uint32_t start = *it;
        if (result.component[start] != UINT32_MAX) continue;

        const uint32_t component = static_cast<uint32_t>(component_size.size());
        uint32_t size = 0;
        result.component[start] = component;
        stack.push_back(start);
        while (!stack.empty()) {
            const uint32_t node = stack.back();
            stack.pop_back();
            ++size;
            for (uint64_t i = in_offset[node]; i < in_offset[node + 1]; ++i) {
                const uint32_t next = in[i];
                if (result.component[next] == UINT32_MAX) {
                    result.component[next] = component;
                    stack.push_back(next);
                }
            }
        }
        component_size.push_back(size);
    }

    result.cyclic.resize(component_size.size(), 0);
    for (uint32_t i = 0; i < component_size.size(); ++i) {
        result.cyclic[i] = component_size[i] > 1;
    }
    for (uint64_t edge : edges) {
        const uint32_t from = static_cast<uint32_t>(edge >> 32);
        const uint32_t to = static_cast<uint32_t>(edge);
        if (from == to) result.cyclic[result.component[from]] = 1;
    }
    return result;
}
/* ================================================================================================================
 *                                          PATH CYCLE GUARD END
 * ================================================================================================================ */
