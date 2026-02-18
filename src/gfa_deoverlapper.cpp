#include "../include/gfa_deoverlapper.hpp"
#include "../include/progress_tracker.hpp"

#include <set>
#include <tuple>
#include <iostream>


void GfaDeoverlapper::prune_overlaps_() {
    log_stream() << "Detecting complementary edges with inconsistent overlap lengths ...\n";

    struct ComplementaryOverlapInfo {
        uint32_t forward_len{0};
        uint32_t reverse_len{0};
    };

    std::unordered_map<uint64_t, ComplementaryOverlapInfo> comp_overlap_map;  // key: (min << 32) | max

    for (const auto& e : arcs_) {
        if (e.get_del() || e.get_comp()) continue;

        const uint32_t v_seg_id = e.get_source_segment_id();
        const uint32_t w_seg_id = e.get_target_segment_id();
        const bool     v_rev    = e.get_source_is_reverse();
        const bool     w_rev    = e.get_target_is_reverse();

        const uint32_t v_len = nodes_[v_seg_id].length;
        const uint32_t w_len = nodes_[w_seg_id].length;

        const uint32_t OV = (e.ov == INT32_MAX ? 0u : clamp32_(e.ov, v_len));
        const uint32_t OW = (e.ow == INT32_MAX ? 0u : clamp32_(e.ow, w_len));
        const uint32_t L  = std::min(OV, OW);
        if (L == 0) continue;

        const uint64_t key    = (static_cast<uint64_t>(std::min(v_seg_id, w_seg_id)) << 32) | std::max(v_seg_id, w_seg_id);
        const bool small_rev  = (v_seg_id < w_seg_id) ? v_rev : w_rev;

        auto& ov = comp_overlap_map.try_emplace(key, ComplementaryOverlapInfo{}).first->second;
        if (small_rev) {
            if (L > ov.reverse_len) ov.reverse_len = L;
        } else {
            if (L > ov.forward_len) ov.forward_len = L;
        }
    }

    // filter: remove equal overlaps
    {
        std::unordered_map<uint64_t, ComplementaryOverlapInfo> tmp;

        for (auto& kv : comp_overlap_map) {
            const auto& ov = kv.second;
            if (ov.forward_len != 0 && ov.reverse_len != 0 && ov.forward_len != ov.reverse_len) {
                tmp.emplace(kv.first, kv.second);
            }
        }
        comp_overlap_map.swap(tmp);

        // remove from the arcs
        for (auto& e : arcs_) {
            if (e.get_del() || e.get_comp()) continue;

            const uint32_t v_seg_id = e.get_source_segment_id();
            const uint32_t w_seg_id = e.get_target_segment_id();

            const uint64_t key = (static_cast<uint64_t>(std::min(v_seg_id, w_seg_id)) << 32) | std::max(v_seg_id, w_seg_id);
            auto find_it = comp_overlap_map.find(key);
            if (find_it != comp_overlap_map.end()) {
                if (std::min(e.ov, e.ow) != std::max(find_it->second.forward_len, find_it->second.reverse_len)) {
                    e.set_del(true);
                }
            }
        }
    }

    // print mismatches
    {
        std::vector<std::pair<uint64_t, ComplementaryOverlapInfo> > mismatches(comp_overlap_map.begin(), comp_overlap_map.end());
        std::sort(mismatches.begin(), mismatches.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

        warning_stream() << "   - Total edges with unequal overlap: " << mismatches.size() << "\n";

        for (const auto& [key, ov] : mismatches) {
            const uint32_t small_seg_id = static_cast<uint32_t>(key >> 32);
            const uint32_t large_seg_id = static_cast<uint32_t>(key & 0xFFFFFFFFu);

            std::string small_name, large_name;
            if (small_seg_id < nodes_.size()) small_name = nodes_[small_seg_id].name;
            if (large_seg_id < nodes_.size()) large_name = nodes_[large_seg_id].name;

            warning_stream() << "      - " << small_name << " and " << large_name << ": " << ov.forward_len << " != " << ov.reverse_len << "\n";
        }
    }
}

std::vector<GfaDeoverlapper::MmWfaHit> GfaDeoverlapper::filter_aligns_(std::vector<MmWfaHit> a)
{
    if (a.size() <= 1) return a;

    auto rb = [](const MmWfaHit& x){ return std::min(x.r_beg, x.r_end); };
    auto re = [](const MmWfaHit& x){ return std::max(x.r_beg, x.r_end); };
    auto qb = [](const MmWfaHit& x){ return std::min(x.q_beg, x.q_end); };
    auto qe = [](const MmWfaHit& x){ return std::max(x.q_beg, x.q_end); };

    auto len = [&](const MmWfaHit& x)->uint32_t {
        uint32_t b = rb(x), e = re(x);
        return (e > b) ? (e - b) : 0u;
    };

    auto ovlp = [](uint32_t l1, uint32_t r1, uint32_t l2, uint32_t r2)->uint32_t {
        uint32_t l = (l1 > l2) ? l1 : l2;
        uint32_t r = (r1 < r2) ? r1 : r2;
        return (r > l) ? (r - l) : 0u;
    };

    std::sort(a.begin(), a.end(), [&](const MmWfaHit& A, const MmWfaHit& B){
        uint32_t la = len(A), lb = len(B);
        if (la != lb) return la > lb;
        return rb(A) < rb(B);
    });

    std::vector<MmWfaHit> kept;
    kept.reserve(a.size());

    for (auto& cand : a) {
        uint32_t crb = rb(cand), cre = re(cand);
        uint32_t cqb = qb(cand), cqe = qe(cand);
        if (cre <= crb || cqe <= cqb) continue;

        int64_t crk = int64_t(crb);
        int64_t cqk = int64_t(cqb);

        bool ok = true;
        for (const auto& k : kept) {
            if (ovlp(crb, cre, rb(k), re(k)) > 0) { ok = false; break; }
            if (ovlp(cqb, cqe, qb(k), qe(k)) > 0) { ok = false; break; }

            int64_t krk = int64_t(rb(k));
            int64_t kqk = int64_t(qb(k));
            if ((crk > krk && cqk <= kqk) || (crk < krk && cqk >= kqk)) { ok = false; break; }
        }

        if (ok) kept.push_back(std::move(cand));
    }

    return kept;
}

std::vector<GfaDeoverlapper::MmWfaHit> GfaDeoverlapper::align_wfa_(
    const std::string& v_name,
    const std::string& w_name,
    const std::string& v_seq_slice,
    const std::string& w_seq_slice
) {
    std::vector<MmWfaHit> hits;

    // ---- Build a tiny minimizer index for v (reference) ----
    std::vector<std::string> v_names{v_name};
    std::vector<std::string_view> v_seqs{v_seq_slice};
    std::vector<std::vector<std::string>> v_right_seqs{{""}};

    mmidx::MinimizerIndex idx(v_names, v_seqs, v_right_seqs, chainOpts_, anchorOpts_);
    idx.build_mm();

    // ---- Single-read aligner over this local index ----
    aligner::Alignmenter aln(idx, v_names, v_seqs, extendOpts_, alignOpts_);

    // ---- Run the per-read pipeline on v (ref) vs w (read) ----
    auto aligns = aln.produce_read(w_name, w_seq_slice, /*keep same strand only=*/true);
    if (aligns.empty()) {
        log_stream() << "   - No alignment produced for " << v_name << " vs " << w_name << "\n";
        return hits;
    }

    hits.reserve(aligns.size());
    for (auto& a : aligns) {
        MmWfaHit h;
        h.r_beg = (uint32_t)a.r_beg;
        h.r_end = (uint32_t)a.r_end;
        h.q_beg = (uint32_t)a.q_beg;
        h.q_end = (uint32_t)a.q_end;

        h.cigar = std::move(a.cigar);
        hits.emplace_back(std::move(h));
    }

    return hits;
}

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "ketopt.h"

static inline std::string mm2_cigar_to_string_(const uint32_t* cigar, uint32_t n_cigar) {
    // mm2 op: M I D N S H P = X B
    static const char op_table[] = "MIDNSHP=XB";
    std::string s;
    s.reserve(n_cigar * 3);
    for (uint32_t i = 0; i < n_cigar; ++i) {
        uint32_t len = cigar[i] >> 4;
        uint32_t op  = cigar[i] & 0x0f;
        if (op >= 10) continue;
        s += std::to_string(len);
        s.push_back(op_table[op]);
    }
    return s;
}

std::vector<GfaDeoverlapper::MmWfaHit> GfaDeoverlapper::align_mm2_(
    const std::string& v_name,
    const std::string& w_name,
    const std::string& v_seq_slice,
    const std::string& w_seq_slice
) {
    std::vector<MmWfaHit> hits;

    mm_idxopt_t ipt;
    mm_mapopt_t opt;
    mm_set_opt(0, &ipt, &opt);
    mm_set_opt("asm5", &ipt, &opt);

    ipt.k = (short)chainOpts_.k;
    ipt.w = (short)chainOpts_.w;

    opt.flag |= MM_F_CIGAR;
    opt.flag |= MM_F_EQX;  // output =/X
    opt.best_n = (short)anchorOpts_.max_kept;

    const char* ref_seqs[1]  = { v_seq_slice.c_str() };
    const char* ref_names[1] = { v_name.c_str() };

    const int is_hpc = (ipt.flag & MM_I_HPC) ? 1 : 0;
    mm_idx_t* mi = mm_idx_str(ipt.w, ipt.k, is_hpc, ipt.bucket_bits, 1, ref_seqs, ref_names);
    if (!mi) return hits;

    mm_mapopt_update(&opt, mi);

    static thread_local mm_tbuf_t* tbuf = nullptr;
    if (!tbuf) tbuf = mm_tbuf_init();

    int n_regs = 0;
    mm_reg1_t* regs = mm_map(mi, (int)w_seq_slice.size(), w_seq_slice.c_str(), &n_regs, tbuf, &opt, w_name.c_str());

    for (int i = 0; i < n_regs; ++i) {
        const mm_reg1_t& r = regs[i];
        if (!r.p) continue;
        if (r.rev) continue;  // keep same strand only
        const mm_extra_t* ex = (const mm_extra_t*)r.p;
        if (!ex || ex->n_cigar == 0) continue;

        std::string cigar = mm2_cigar_to_string_(ex->cigar, ex->n_cigar);
        if (cigar.empty()) continue;

        MmWfaHit h;
        h.r_beg = (uint32_t)r.rs;
        h.r_end = (uint32_t)r.re;
        h.q_beg = (uint32_t)r.qs;
        h.q_end = (uint32_t)r.qe;

        h.cigar = std::move(cigar);
        hits.emplace_back(std::move(h));
    }

    if (regs) {
        for (int i = 0; i < n_regs; ++i) free(regs[i].p);
        free(regs);
    }
    mm_idx_destroy(mi);

    return hits;
}

// Align v-slice (reference) against w-slice (query) and pack as BubbleAlignment.
std::vector<GfaDeoverlapper::BubbleAlignment> GfaDeoverlapper::align_and_pack_(
    uint32_t            v_vertex,           // packed vertex id for A
    uint32_t            w_vertex,           // packed vertex id for B
    const std::string&  v_name,             // segment name of A
    const std::string&  w_name,             // segment name of B
    const std::string&  v_seq_slice,        // oriented slice sequence of A (already rev-comp if needed)
    const std::string&  w_seq_slice,        // oriented slice sequence of B (already rev-comp if needed)
    uint32_t            vb, uint32_t ve,    // slice coordinates on A in plus strand: [vb, ve)
    uint32_t            wb, uint32_t we     // slice coordinates on B in plus strand: [wb, we)
) {
    std::vector<BubbleAlignment> outs;

    if (v_seq_slice.empty() || w_seq_slice.empty() || v_seq_slice == "*" || w_seq_slice == "*") {
        warning_stream() << "   - Empty/unknown sequence for " << v_name << " vs " << w_name << "\n";
        return outs;
    }

    const bool v_rev = NodeHandle::get_is_reverse(v_vertex);
    const bool w_rev = NodeHandle::get_is_reverse(w_vertex);

    auto add_alignment = [&](
        uint32_t r_beg, uint32_t r_end,
        uint32_t q_beg, uint32_t q_end,
        std::string cigar,
        std::vector<CIGAR::COp>&& ops
    ) {
        // Convert to original segment coordinates
        //   forward:  [vb + r_beg, vb + r_end)
        //   reverse:  [ve - r_end, ve - r_beg)
        const uint32_t beg_a = v_rev ? (ve - r_end) : (vb + r_beg);
        const uint32_t end_a = v_rev ? (ve - r_beg) : (vb + r_end);
        const uint32_t beg_b = w_rev ? (we - q_end) : (wb + q_beg);
        const uint32_t end_b = w_rev ? (we - q_beg) : (wb + q_end);

        debug_stream() << "   - ref: " << v_name << "(" << (v_rev ? "-" : "+") << ")\n";
        debug_stream() << "   - qry: " << w_name << "(" << (w_rev ? "-" : "+") << ")\n";
        debug_stream() << "     - [r_beg, r_end, r_len) = [" << beg_a << ", " << end_a << ", " << v_seq_slice.length() << ")\n";
        debug_stream() << "     - [q_beg, q_end, q_len) = [" << beg_b << ", " << end_b << ", " << w_seq_slice.length() << ")\n";
        debug_stream() << "     - CIGAR = " << cigar << "\n";

        outs.emplace_back(
            v_name, w_name,
            beg_a, end_a,
            beg_b, end_b,
            v_vertex, w_vertex,
            std::move(ops)
        );
    };

    if (v_seq_slice == w_seq_slice) {
        const uint32_t L = (uint32_t)v_seq_slice.size();
        if (L == 0) return outs;
        std::string cigar = std::to_string(L) + "=";
        auto ops = CIGAR::parse(cigar);
        add_alignment(0, L, 0, L, std::move(cigar), std::move(ops));
        return outs;
    }

    std::vector<MmWfaHit> hits = use_wfa_
        ? align_wfa_(v_name, w_name, v_seq_slice, w_seq_slice)
        : align_mm2_(v_name, w_name, v_seq_slice, w_seq_slice);

    hits = filter_aligns_(std::move(hits));
    
    for (auto& h : hits) {
        if (h.cigar.empty()) continue;
        auto ops = CIGAR::parse(h.cigar);
        add_alignment(h.r_beg, h.r_end, h.q_beg, h.q_end, std::move(h.cigar), std::move(ops));
    }

    return outs;
}

void GfaDeoverlapper::initialize_cuts_() {
    log_stream() << "Initializing cut points ...\n";

    const uint32_t orig_seg_n = static_cast<uint32_t>(nodes_.size());

    cuts_.assign(orig_seg_n, Cuts{});
    for (uint32_t s = 0; s < orig_seg_n; ++s) {
        cuts_[s].name = nodes_[s].name;
        cuts_[s].v    = {0u, nodes_[s].length};
    }
}


void GfaDeoverlapper::overlaps_align_() {
    log_stream() << "Aligning edge overlaps ...\n";

    using namespace wfa;  // WFAlignerGapAffine

    // Progress tracker for merging unitigs
    ProgressTracker prog(arcs_.size());

    // to avoid redundant alignment of (v,w) and (w,v)
    std::unordered_set<uint64_t> seen_pairs;
    seen_pairs.reserve(arcs_.size() * 2);
    auto pair_key = [](uint32_t a, uint32_t b) -> uint64_t {
        uint32_t x = std::min(a, b);
        uint32_t y = std::max(a, b);
        return (uint64_t(x) << 32) | uint64_t(y);
    };

    // edges contribute window boundaries
    for (const auto &e : arcs_) {
        prog.hit();  // Update progress

        if (e.get_del() || e.get_comp()) continue;
        uint32_t v_vertex = e.get_source_vertex_id();
        uint32_t w_vertex = e.get_target_vertex_id();
        uint32_t v_seg_id = e.get_source_segment_id();
        uint32_t w_seg_id = e.get_target_segment_id();
        bool     v_rev    = e.get_source_is_reverse();
        bool     w_rev    = e.get_target_is_reverse();
        std::string v_name = nodes_[v_seg_id].name, w_name = nodes_[w_seg_id].name;

        uint32_t v_len = nodes_[v_seg_id].length, w_len = nodes_[w_seg_id].length;
        uint32_t OV = (e.ov == INT32_MAX ? 0u : clamp32_(e.ov, v_len));
        uint32_t OW = (e.ow == INT32_MAX ? 0u : clamp32_(e.ow, w_len));
        const uint32_t L = std::min(OV, OW);
        if (L == 0) continue;

        auto [vb, ve] = v_overlap_pos_(v_len, OV, v_rev);
        auto [wb, we] = w_overlap_pos_(w_len, OW, w_rev);

        // record cut points (overlap endpoints) 2025-10-14
        cuts_[v_seg_id].v.push_back(vb); cuts_[v_seg_id].v.push_back(ve);
        cuts_[w_seg_id].v.push_back(wb); cuts_[w_seg_id].v.push_back(we);

        // avoid redundant alignment of (v,w) and (w,v).
        // Don't change the order of code, because the cut points need to be recorded first.
        uint64_t key = pair_key(v_seg_id, w_seg_id);
        if (!seen_pairs.insert(key).second) continue;

        std::string v_seq_slice = slice_seq_or_star_(nodes_, v_seg_id, vb, ve, v_rev);
        std::string w_seq_slice = slice_seq_or_star_(nodes_, w_seg_id, wb, we, w_rev);

        if (v_seq_slice == "*" || w_seq_slice == "*" || v_seq_slice.empty() || w_seq_slice.empty()) {
            warning_stream()
                << "   - Overlap sequences unknown (\"*\"): "
                << v_name << "(" << vb << "-" << ve << (v_rev?":-":":+") << ") vs "
                << w_name << "(" << wb << "-" << we << (w_rev?":-":":+") << "), len=" << L
                << "\n";
            continue;
        }

        std::vector<BubbleAlignment> alns = align_and_pack_(
            v_vertex, w_vertex,
            v_name, w_name,
            v_seq_slice, w_seq_slice,
            vb, ve, wb, we
        );

        if (alns.empty()) continue;
        bubble_aligns_.insert(bubble_aligns_.end(), std::make_move_iterator(alns.begin()), std::make_move_iterator(alns.end()));
    }

    log_stream() << "   - Total alignments produced: " << bubble_aligns_.size() << "\n";
}

void GfaDeoverlapper::dedup_aligns_()
{
    gfaName namer;

    auto seg_name = [this](uint32_t v) -> const std::string& {
        return nodes_[NodeHandle::get_segment_id(v)].name;
    };

    auto norm_interval = [] (uint32_t a, uint32_t b) {
        uint32_t l = std::min(a, b);
        uint32_t r = std::max(a, b);
        return std::pair<uint32_t,uint32_t>(l, r);
    };

    auto overlap = [] (uint32_t l1, uint32_t r1, uint32_t l2, uint32_t r2) -> uint32_t {
        uint32_t l = std::max(l1, l2);
        uint32_t r = std::min(r1, r2);
        return (r > l) ? (r - l) : 0u;
    };

    struct SpanItem {
        size_t      aln_idx;   // Index in bubble_aligns_
        std::string base;      // The base segment name
        uint32_t    beg;       // The real start position on the base
        uint32_t    end;       // The real end position on the base
        uint32_t    len;       // end - beg
    };

    const size_t n_aln = bubble_aligns_.size();
    std::unordered_map<std::string, std::vector<SpanItem>> groups;
    groups.reserve(n_aln * 4);

    for (size_t i = 0; i < n_aln; ++i) {
        const auto& aln = bubble_aligns_[i];

        auto [a_l, a_r] = norm_interval(aln.beg_a, aln.end_a);
        auto [b_l, b_r] = norm_interval(aln.beg_b, aln.end_b);
        if (a_r <= a_l || b_r <= b_l) continue;

        const std::string& sA = seg_name(aln.v_a);
        const std::string& sB = seg_name(aln.v_b);

        // key=ref, value=qry
        {
            std::string fullB = namer.format_interval_name(sB, b_l, b_r, false);
            auto piecesB = namer.parse_composite_with_dir(fullB);

            for (const auto& pc : piecesB) {
                if (pc.len == 0) continue;
                uint32_t lo = static_cast<uint32_t>(pc.lo);
                uint32_t hi = static_cast<uint32_t>(pc.hi);
                if (hi <= lo) continue;
                uint32_t len = hi - lo;

                groups[sA].push_back(SpanItem{i, pc.root, lo, hi, len});
            }
        }

        // key=qry, value=ref
        {
            std::string fullA = namer.format_interval_name(sA, a_l, a_r, false);
            auto piecesA = namer.parse_composite_with_dir(fullA);

            for (const auto& pc : piecesA) {
                if (pc.len == 0) continue;
                uint32_t lo = static_cast<uint32_t>(pc.lo);
                uint32_t hi = static_cast<uint32_t>(pc.hi);
                if (hi <= lo) continue;
                uint32_t len = hi - lo;

                groups[sB].push_back(SpanItem{i, pc.root, lo, hi, len});
            }
        }
    }

    // Sort and dedup within each ref group
    std::vector<uint8_t> keep(n_aln, 1);
    for (auto &kv : groups) {
        const std::string& ref = kv.first;
        auto &spans            = kv.second;

        if (spans.size() <= 1) continue;

        std::sort(spans.begin(), spans.end(),
            [] (const SpanItem &x, const SpanItem &y) {
                if (x.base != y.base) return x.base < y.base;
                if (x.len  != y.len ) return x.len  > y.len;
                if (x.beg  != y.beg ) return x.beg  < y.beg;
                return x.end < y.end;
            }
        );

        std::unordered_map<std::string, std::vector<SpanItem>> chosen_by_base;
        chosen_by_base.reserve(8);

        for (const auto& s : spans) {
            if (!keep[s.aln_idx]) continue;

            auto &chosen = chosen_by_base[s.base];
            bool conflict = false;

            for (const auto& c : chosen) {
                uint32_t ov = overlap(s.beg, s.end, c.beg, c.end);
                if (ov > 0) {
                    keep[s.aln_idx] = 0;
                    conflict = true;
                    break;
                }
            }

            if (!conflict) {
                chosen.push_back(s);
            }
        }
    }


    // Rebuild bubble_aligns_
    std::vector<BubbleAlignment> new_aligns;
    new_aligns.reserve(n_aln);

    for (size_t i = 0; i < n_aln; ++i) {
        if (keep[i]) {
            new_aligns.push_back(std::move(bubble_aligns_[i]));
        }
    }

    bubble_aligns_ = std::move(new_aligns);

    // Print results
    for (const auto& aln : bubble_aligns_) {
        auto [b_l_raw, b_r_raw] = norm_interval(aln.beg_b, aln.end_b);
        if (b_r_raw <= b_l_raw) continue;
        const std::string& snameB = seg_name(aln.v_b);
        std::string fullB = namer.format_interval_name(snameB, b_l_raw, b_r_raw, false);
    }

    log_stream() << "   - After base-overlap dedup, alignments kept: " << bubble_aligns_.size() << " / " << n_aln << "\n";
}

void GfaDeoverlapper::build_cuts_from_cigar_() {
    log_stream() << "Building cut points through alignments ...\n";

    // Progress tracker for merging unitigs
    ProgressTracker prog(bubble_aligns_.size());

    for (const auto& align : bubble_aligns_) {
        prog.hit();  // Update progress

        uint32_t seg_a = NodeHandle::get_segment_id(align.v_a);
        uint32_t seg_b = NodeHandle::get_segment_id(align.v_b);
        bool     rev_a = NodeHandle::get_is_reverse(align.v_a);
        bool     rev_b = NodeHandle::get_is_reverse(align.v_b);
        if (seg_a >= cuts_.size() || seg_b >= cuts_.size()) continue;

        uint32_t len_a = nodes_[seg_a].length;
        uint32_t len_b = nodes_[seg_b].length;

        uint32_t pos_a = rev_a ? align.end_a : align.beg_a;
        uint32_t pos_b = rev_b ? align.end_b : align.beg_b;

        // guard coordinates
        if (!(align.beg_a < align.end_a && align.end_a <= len_a)) continue;
        if (!(align.beg_b < align.end_b && align.end_b <= len_b)) continue;

        // step function to move along the sequences
        auto step = [](uint32_t& p, bool rev, uint32_t n) {
            if (n == 0) return;
            p = rev ? (p - n) : (p + n);
        };

        for (const auto& op : align.ops) {
            switch (op.op) {
                case 'M':
                case '=': {  // match or sequence match
                    const uint32_t prev_a = pos_a;
                    const uint32_t prev_b = pos_b;

                    step(pos_a, rev_a, op.len);
                    step(pos_b, rev_b, op.len);

                    if (op.len >= MIN_EQ_FOR_CUT_) {
                        if (prev_a > 0 && prev_a < len_a) cuts_[seg_a].v.push_back(prev_a);
                        if (prev_b > 0 && prev_b < len_b) cuts_[seg_b].v.push_back(prev_b);
                        if (pos_a  > 0 && pos_a  < len_a) cuts_[seg_a].v.push_back(pos_a);
                        if (pos_b  > 0 && pos_b  < len_b) cuts_[seg_b].v.push_back(pos_b);
                    }
                } break;

                case 'X': { // mismatch
                    step(pos_a, rev_a, op.len);
                    step(pos_b, rev_b, op.len);
                } break;

                case 'I': { // insertion
                    step(pos_b, rev_b, op.len);
                } break;

                case 'D': { // deletion
                    step(pos_a, rev_a, op.len);
                } break;

                case 'S': // soft clip
                case 'H': // hard clip
                    break;

                default:
                    warning_stream() << "   - Unknown CIGAR op: " << op.op << "\n";
                    break;
            }

            if (pos_a > len_a || pos_b > len_b) {
                warning_stream() 
                    << "   - CIGAR stepping overflow: pos_a=" << pos_a
                    << " len_a=" << len_a << " pos_b=" << pos_b
                    << " len_b=" << len_b << "\n";
                break;
            }
        }
    }

    // sort unique + ensure boundary
    for (uint32_t s = 0; s < cuts_.size(); ++s) {
        auto &C = cuts_[s].v;
        std::sort(C.begin(), C.end());
        C.erase(std::unique(C.begin(), C.end()), C.end());
        if (C.empty() || C.front() != 0u) C.insert(C.begin(), 0u);
        if (C.back() != nodes_[s].length) C.push_back(nodes_[s].length);
    }
}


void GfaDeoverlapper::propagate_cuts_() {
    log_stream() << "Propagating cut points through alignments ...\n";

    int iter = 0;
    while (iter < MAX_PROPAGATION_ITERS_) {
        uint64_t aligns_changed = 0;
        ++iter;

        for (const auto& align : bubble_aligns_) {
            uint32_t seg_a = NodeHandle::get_segment_id(align.v_a);
            uint32_t seg_b = NodeHandle::get_segment_id(align.v_b);
            bool     rev_a = NodeHandle::get_is_reverse(align.v_a);
            bool     rev_b = NodeHandle::get_is_reverse(align.v_b);
            std::string name_a = align.name_a;
            std::string name_b = align.name_b;
            if (seg_a >= cuts_.size() || seg_b >= cuts_.size()) continue;

            uint32_t len_a = nodes_[seg_a].length;
            uint32_t len_b = nodes_[seg_b].length;

            if (!(align.beg_a < align.end_a && align.end_a <= len_a)) continue;
            if (!(align.beg_b < align.end_b && align.end_b <= len_b)) continue;

            uint32_t pos_a = rev_a ? align.end_a : align.beg_a;
            uint32_t pos_b = rev_b ? align.end_b : align.beg_b;

            auto& a_cuts = cuts_[seg_a].v;
            auto& b_cuts = cuts_[seg_b].v;

            std::vector<uint32_t> new_a_cuts, new_b_cuts;

            bool leaf_is_a;
            if (!prefer_v_as_leaf_(seg_a, seg_b, leaf_is_a)) continue;

            auto step = [](uint32_t& p, bool rev, uint32_t n) {
                if (n == 0) return;
                p = rev ? (p - n) : (p + n);
            };

            for (const auto& op : align.ops) {
                switch (op.op) {
                    case 'M':
                    case '=': {
                        const uint32_t prev_a = pos_a;
                        const uint32_t prev_b = pos_b;

                        step(pos_a, rev_a, op.len);
                        step(pos_b, rev_b, op.len);

                        const auto [beg_a, end_a] = rev_a
                            ? std::make_pair(pos_a, prev_a)
                            : std::make_pair(prev_a, pos_a);
                        const auto [beg_b, end_b] = rev_b
                            ? std::make_pair(pos_b, prev_b)
                            : std::make_pair(prev_b, pos_b);

                        auto it_a_1 = std::lower_bound(a_cuts.begin(), a_cuts.end(), beg_a);
                        auto it_a_2 = std::upper_bound(a_cuts.begin(), a_cuts.end(), end_a);
                        auto it_b_1 = std::lower_bound(b_cuts.begin(), b_cuts.end(), beg_b);
                        auto it_b_2 = std::upper_bound(b_cuts.begin(), b_cuts.end(), end_b);
                        std::vector<uint32_t> cut_in_align_a(it_a_1, it_a_2);
                        std::vector<uint32_t> cut_in_align_b(it_b_1, it_b_2);

                        if (leaf_is_a) {
                            for (const uint32_t cut_b_pos : cut_in_align_b) {
                                const uint32_t offset       = rev_b ? (end_b - cut_b_pos) : (cut_b_pos - beg_b);
                                const uint32_t mapped_a_pos = rev_a ? (end_a - offset)    : (beg_a + offset);
                                if (mapped_a_pos > 0 && mapped_a_pos < len_a) new_a_cuts.push_back(mapped_a_pos);
                            }
                        } else {
                            for (const uint32_t cut_a_pos : cut_in_align_a) {
                                const uint32_t offset       = rev_a ? (end_a - cut_a_pos) : (cut_a_pos - beg_a);
                                const uint32_t mapped_b_pos = rev_b ? (end_b - offset)    : (beg_b + offset);
                                if (mapped_b_pos > 0 && mapped_b_pos < len_b) new_b_cuts.push_back(mapped_b_pos);
                            }
                        }
                    } break;

                    case 'X': { // mismatch
                        step(pos_a, rev_a, op.len);
                        step(pos_b, rev_b, op.len);
                    } break;

                    case 'I': { // insertion
                        step(pos_b, rev_b, op.len);
                    } break;

                    case 'D': { // deletion
                        step(pos_a, rev_a, op.len);
                    } break;

                    case 'S': // soft clip
                    case 'H': // hard clip
                        break;

                    default:
                        warning_stream() << "   - Unknown CIGAR op: " << op.op << "\n";
                        break;
                }

                if (pos_a > len_a || pos_b > len_b) {
                    warning_stream() 
                        << "   - CIGAR stepping overflow: pos_a=" << pos_a
                        << " len_a=" << len_a << " pos_b=" << pos_b
                        << " len_b=" << len_b << "\n";
                    break;
                }
            }

            bool increase_a = merge_and_dedup_cuts_(a_cuts, new_a_cuts);
            bool increase_b = merge_and_dedup_cuts_(b_cuts, new_b_cuts);
            if (increase_a || increase_b) {
                aligns_changed++;
            }
        }
        log_stream() << "   - [Iter " << iter << "] " << aligns_changed << " aligns had cut granularity changes\n";
        if (aligns_changed == 0) break;
    }
}

void GfaDeoverlapper::print_cuts(const std::vector<Cuts>& cuts) {
    log_stream() << "Segment cuts:\n";
    for (size_t i = 0; i < cuts.size(); ++i) {
        const auto& cut = cuts[i];
        std::string log = "   - " + cut.name + ": ";
        for (const auto& pos : cut.v) log += std::to_string(pos) + " ";
        log += "\n";
        log_stream() << log;
    }
}


void GfaDeoverlapper::build_rules_from_pair_windows_(
    SegReplace::RuleMap& rulemap,
    // trunk
    uint32_t trunk_seg_id, bool trunk_is_rev,
    const std::vector<uint32_t>& cut_in_overlap_trunk,
    // leaf
    uint32_t leaf_seg_id, bool leaf_is_rev,
    const std::vector<uint32_t>& new_cuts_leaf,
    const std::vector<int32_t>&  new_offsets_leaf,
    // leaf older cuts (in overlap region only)
    const std::vector<uint32_t>& cut_in_overlap_leaf
) const {
    // ---------- A) trunk -> leaf ----------
    if (cut_in_overlap_trunk.size() >= 2 && new_cuts_leaf.size() >= 2 && new_offsets_leaf.size() >= 2) {
        for (size_t j = 1; j < cut_in_overlap_trunk.size(); ++j) {
            const uint32_t trunk_beg = trunk_is_rev ? cut_in_overlap_trunk[j]     : cut_in_overlap_trunk[j - 1];
            const uint32_t trunk_end = trunk_is_rev ? cut_in_overlap_trunk[j - 1] : cut_in_overlap_trunk[j];

            const int32_t trunk_beg_off = trunk_is_rev
                ? (static_cast<int32_t>(cut_in_overlap_trunk[0]) - static_cast<int32_t>(trunk_end))
                : (static_cast<int32_t>(trunk_beg) - static_cast<int32_t>(cut_in_overlap_trunk[0]));
            const int32_t trunk_end_off = trunk_is_rev
                ? (static_cast<int32_t>(cut_in_overlap_trunk[0]) - static_cast<int32_t>(trunk_beg))
                : (static_cast<int32_t>(trunk_end) - static_cast<int32_t>(cut_in_overlap_trunk[0]));

            const uint32_t leaf_beg_idx = find_index_by_offset(trunk_beg_off, new_offsets_leaf);
            const uint32_t leaf_end_idx = find_index_by_offset(trunk_end_off, new_offsets_leaf);
            if (leaf_beg_idx == UINT32_MAX || leaf_end_idx == UINT32_MAX || leaf_beg_idx == leaf_end_idx) continue;

            std::vector<SegReplace::Seg> leaf_parts;
            leaf_parts.reserve(leaf_end_idx - leaf_beg_idx);
            for (uint32_t k = leaf_beg_idx; k <= leaf_end_idx - 1; ++k) {
                uint32_t lb = new_cuts_leaf[k];
                uint32_t le = new_cuts_leaf[k + 1];
                if (lb > le) std::swap(lb, le);
                leaf_parts.push_back(SegReplace::Interval::pack(leaf_seg_id, lb, le, leaf_is_rev));
            }

            SegReplace::Seg trunk_u128 = SegReplace::Interval::pack(trunk_seg_id, trunk_beg, trunk_end, trunk_is_rev);
            if (trunk_is_rev) {
                trunk_u128 = SegReplace::Interval::toggle_strand(trunk_u128);
                SegReplace::Expander::reverse_and_toggle(leaf_parts);
            }
            rulemap[trunk_u128] = std::move(leaf_parts);
        }
    }

    // ---------- B) Update leaf cuts ----------
    if (cut_in_overlap_leaf.size() >= 2 && new_cuts_leaf.size() >= 2) {
        std::vector<uint32_t> leaf_win_sorted = cut_in_overlap_leaf;
        std::sort(leaf_win_sorted.begin(), leaf_win_sorted.end());

        std::vector<uint32_t> cuts_leaf_sorted = new_cuts_leaf;
        std::sort(cuts_leaf_sorted.begin(), cuts_leaf_sorted.end());

        for (size_t j = 1; j < leaf_win_sorted.size(); ++j) {
            const uint32_t leaf_win_beg = leaf_win_sorted[j - 1];
            const uint32_t leaf_win_end = leaf_win_sorted[j];

            const auto it1 = std::lower_bound(cuts_leaf_sorted.begin(), cuts_leaf_sorted.end(), leaf_win_beg);
            const auto it2 = std::upper_bound(cuts_leaf_sorted.begin(), cuts_leaf_sorted.end(), leaf_win_end);
            std::vector<uint32_t> sub(it1, it2);

            if (sub.size() < 3) continue;

            std::vector<SegReplace::Seg> leaf_parts;
            leaf_parts.reserve(sub.size() - 1);
            for (size_t k = 1; k < sub.size(); ++k) {
                uint32_t lb = sub[k - 1];
                uint32_t le = sub[k];
                if (lb > le) std::swap(lb, le);
                leaf_parts.push_back(SegReplace::Interval::pack(leaf_seg_id, lb, le, false));
            }

            SegReplace::Seg leaf_win_u128 = SegReplace::Interval::pack(leaf_seg_id, leaf_win_beg, leaf_win_end, false);
            rulemap[leaf_win_u128] = std::move(leaf_parts);
        }
    }
}


void GfaDeoverlapper::record_rulemap(
    uint32_t v_seg_id, uint32_t w_seg_id,
    bool v_is_rev, bool w_is_rev,
    const std::string& v_name, const std::string& w_name,
    uint32_t vb, uint32_t ve, uint32_t wb, uint32_t we,
    SegReplace::RuleMap& rulemap
) {
    // Helper to calculate offset from cut point
    auto offset_from = [](uint32_t cut, uint32_t beg, uint32_t end, bool rev) -> uint64_t {
        if (beg >= end) {
            error_stream() << "Invalid bounds: [" << beg << ", " << end << "]\n";
            std::exit(1);
        }
        if (cut < beg || cut > end) {
            error_stream() << "Cut point " << cut << " out of bounds [" << beg << ", " << end << "]\n";
            std::exit(1);
        }
        return rev ? (uint64_t)end - cut : (uint64_t)cut - beg;
    };
    auto offset_to = [](uint64_t off, uint32_t beg, uint32_t end, bool rev) -> uint32_t {
        if (beg >= end) {
            error_stream() << "Invalid bounds: [" << beg << ", " << end << "]\n";
            std::exit(1);
        }
        if (off > (uint64_t)(end - beg)) {
            error_stream() << "Offset " << off << " out of bounds [0, " << (end - beg) << "]\n";
            std::exit(1);
        }
        uint64_t pos = rev ? (uint64_t)end - off : (uint64_t)beg + off;
        return (uint32_t)pos;
    };
    auto build_offsets = [&](const std::vector<uint32_t>& cuts, uint32_t beg, uint32_t end, bool rev) -> std::vector<int32_t> {
        std::vector<int32_t> offsets;
        offsets.reserve(cuts.size());
        for (uint32_t c : cuts) {
            offsets.push_back((int32_t)offset_from(c, beg, end, rev));
        }
        return offsets;
    };

    auto& v_cuts = cuts_[v_seg_id].v;
    auto& w_cuts = cuts_[w_seg_id].v;

    auto it_v_1 = std::lower_bound(v_cuts.begin(), v_cuts.end(), vb);
    auto it_v_2 = std::upper_bound(v_cuts.begin(), v_cuts.end(), ve);
    auto it_w_1 = std::lower_bound(w_cuts.begin(), w_cuts.end(), wb);
    auto it_w_2 = std::upper_bound(w_cuts.begin(), w_cuts.end(), we);
    std::vector<uint32_t> cut_in_overlap_v(it_v_1, it_v_2);
    std::vector<uint32_t> cut_in_overlap_w(it_w_1, it_w_2);

    bool v_is_leaf;
    if (!prefer_v_as_leaf_(v_seg_id, w_seg_id, v_is_leaf)) return;

    auto&        leaf_cuts     = v_is_leaf ? v_cuts : w_cuts;
    auto&        trunk_cuts    = v_is_leaf ? w_cuts : v_cuts;
    auto&        leaf_vec      = v_is_leaf ? cut_in_overlap_v : cut_in_overlap_w;
    auto&        trunk_vec     = v_is_leaf ? cut_in_overlap_w : cut_in_overlap_v;
    bool         leaf_is_rev   = v_is_leaf ? v_is_rev : w_is_rev;
    bool         trunk_is_rev  = v_is_leaf ? w_is_rev : v_is_rev;
    uint32_t     leaf_seg_id   = v_is_leaf ? v_seg_id : w_seg_id;
    uint32_t     trunk_seg_id  = v_is_leaf ? w_seg_id : v_seg_id;
    const auto&  leaf_name     = v_is_leaf ? v_name   : w_name;
    const auto&  trunk_name    = v_is_leaf ? w_name   : v_name;
    uint32_t     lb            = v_is_leaf ? vb       : wb;
    uint32_t     le            = v_is_leaf ? ve       : we;
    uint32_t     tb            = v_is_leaf ? wb       : vb;
    uint32_t     te            = v_is_leaf ? we       : ve;

    // reverse according to direction
    if (leaf_is_rev)  std::reverse(leaf_vec.begin(),  leaf_vec.end());
    if (trunk_is_rev) std::reverse(trunk_vec.begin(), trunk_vec.end());

    // Compute leaf offsets
    std::vector<int32_t> trunk_offsets = build_offsets(trunk_vec, tb, te, trunk_is_rev);
    std::vector<int32_t> leaf_offsets  = build_offsets(leaf_vec, lb, le, leaf_is_rev);

    for (size_t l_idx = 0; l_idx < trunk_offsets.size() - 1; ++l_idx) {
        size_t r_idx = l_idx + 1;
        uint32_t trunk_beg = trunk_is_rev ? trunk_vec[r_idx] : trunk_vec[l_idx];
        uint32_t trunk_end = trunk_is_rev ? trunk_vec[l_idx] : trunk_vec[r_idx];

        int32_t trunk_beg_off = trunk_offsets[l_idx];
        int32_t trunk_end_off = trunk_offsets[r_idx];
        uint32_t leaf_beg_idx = find_index_by_offset(trunk_beg_off, leaf_offsets);
        uint32_t leaf_end_idx = find_index_by_offset(trunk_end_off, leaf_offsets);

        if (leaf_beg_idx == UINT32_MAX) { continue; }

        // 2026-01-08
        //  Example:
        //     Trunk offsets: 0 30 50 70 120
        //     Leaf  offsets: 0 30    70 120
        // merge_and_dedup_cuts_ will add new cuts, so there will be some cuts in trunk_vec that are not in leaf_vec.
        // In that case, we will ignore the missing cuts by extending the trunk window until we find a matching cut in trunk_vec.
        while (leaf_end_idx == UINT32_MAX && (r_idx + 1) < trunk_offsets.size()) {
            ++r_idx;
            if (trunk_is_rev) { trunk_beg = trunk_vec[r_idx]; }
            else              { trunk_end = trunk_vec[r_idx]; }
            trunk_end_off = trunk_offsets[r_idx];
            leaf_end_idx = find_index_by_offset(trunk_end_off, leaf_offsets);
        }

        if (leaf_end_idx == UINT32_MAX) {
            warning_stream() << "Cut offset not found!\n";
            warning_stream() << "   - " << leaf_name << " (leaf) and " << trunk_name << " (trunk)\n";
            warning_stream() << "   - Requested offset=" << trunk_end_off << "\n";
            std::string log_tmp = "   - Available leaf offsets=(";
            for (size_t i = 0; i < leaf_offsets.size(); i++) {
                log_tmp += std::to_string(leaf_offsets[i]);
                if (i + 1 < leaf_offsets.size()) log_tmp += ", ";
            }
            log_tmp += ")\n";
            warning_stream() << log_tmp;
            continue;
        }
        if (leaf_end_idx <= leaf_beg_idx) {
            warning_stream() << "Unexpected index detected!\n";
            warning_stream() << "   - trunk_offset=" << trunk_beg_off << " - " << trunk_end_off << "\n";
            warning_stream() << "   - leaf_idx=" << leaf_beg_idx << " - " << leaf_end_idx << "\n";
            std::string log_tmp = "   - leaf_offsets=(";
            for (size_t i = 0; i < leaf_offsets.size(); i++) {
                log_tmp += std::to_string(leaf_offsets[i]);
                if (i + 1 < leaf_offsets.size()) log_tmp += ", ";
            }
            log_tmp += ")\n";
            warning_stream() << log_tmp;
            continue;
        }


        if (DEBUG_ENABLED) {
            std::cerr << "Processing edge " << v_name << "(" << v_seg_id << ",rev=" << v_is_rev << ") <-> " << w_name << "(" << w_seg_id << ",rev=" << w_is_rev << ")\n"
                    << "   - Trunk segment overlap window: [" << trunk_vec.front() << " .. " << trunk_vec.back() << "]\n"
                    << "   - Leaf  segment overlap window: [" << leaf_vec.front()  << " .. " << leaf_vec.back()  << "]\n"
                    << "   - Current trunk position: beg=" << trunk_beg << " end=" << trunk_end << "\n"
                    << "   - Current trunk offsets: beg=" << trunk_beg_off << " end=" << trunk_end_off << "\n"
                    << "   - Mapped leaf idx range: " << leaf_beg_idx << " .. " << leaf_end_idx << "\n";

            std::cerr << "   - trunk_vec (" << trunk_name << "): ";
            for (auto p : trunk_vec) std::cerr << p << " ";
            std::cerr << "\n";

            std::cerr << "   - leaf_offsets: ";
            for (auto off : leaf_offsets) std::cerr << off << " ";
            std::cerr << "\n";
        }

        std::vector<SegReplace::Seg> leaf_segs;
        for (uint32_t i = leaf_beg_idx; i <= leaf_end_idx - 1; ++i) {
            uint32_t leaf_beg = leaf_vec[i];
            uint32_t leaf_end = leaf_vec[i + 1];
            if (leaf_beg > leaf_end) std::swap(leaf_beg, leaf_end);
            SegReplace::Seg seg_u128 = SegReplace::Interval::pack(leaf_seg_id, leaf_beg, leaf_end, leaf_is_rev);
            leaf_segs.push_back(seg_u128);
        }
        SegReplace::Seg trunk_u128 = SegReplace::Interval::pack(trunk_seg_id, trunk_beg, trunk_end, trunk_is_rev);
        if (trunk_is_rev) {
            trunk_u128 = SegReplace::Interval::toggle_strand(trunk_u128);
            SegReplace::Expander::reverse_and_toggle(leaf_segs);
        }

        // Add to rulemap
        auto it_rule = rulemap.find(trunk_u128);
        if (it_rule == rulemap.end()) {
            rulemap[trunk_u128] = leaf_segs;
        } else {
            auto exist_segs = it_rule->second;

            if (exist_segs.size() == 0 || leaf_segs.size() == 0) {
                error_stream() << "Empty segs detected!\n";
                error_stream() << "   - Edge: " << v_name << "(" << v_seg_id << ",rev=" << v_is_rev << ") <-> " << w_name << "(" << w_seg_id << ",rev=" << w_is_rev << ")\n";

                std::string log_tmp;

                log_tmp = "   - exist segs: ";
                if (exist_segs.empty()) log_tmp += "(EMPTY)";
                else for (auto& seg : exist_segs) log_tmp += SegReplace::Interval::format(seg) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                log_tmp = "   - leaf_segs: ";
                if (leaf_segs.empty()) log_tmp += "(EMPTY)";
                else for (auto& seg : leaf_segs) log_tmp += SegReplace::Interval::format(seg) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                log_tmp = "   - trunk_vec (" + trunk_name + "): ";
                for (auto p : trunk_vec) log_tmp += std::to_string(p) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                error_stream() << "   - Trunk range: beg=" << trunk_beg << " end=" << trunk_end << " offsets=(" << trunk_beg_off << "," << trunk_end_off << ")\n";

                log_tmp = "   - leaf_vec (" + leaf_name + "): ";
                for (auto p : leaf_vec) log_tmp += std::to_string(p) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                log_tmp = "   - leaf_offsets: ";
                for (auto off : leaf_offsets) log_tmp += std::to_string(off) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                error_stream() << "   - Leaf index range: " << leaf_beg_idx << " .. " << leaf_end_idx
                                << " (index_len=" << (leaf_end_idx >= leaf_beg_idx ? (leaf_end_idx - leaf_beg_idx) : -1) << ")\n";

                if (!leaf_vec.empty() && leaf_beg_idx < leaf_vec.size() && leaf_end_idx < leaf_vec.size()) {
                    error_stream() << "   - Leaf range: beg=" << leaf_vec[leaf_beg_idx] << " end=" << leaf_vec[leaf_end_idx] << "\n";
                }
                exit(1);
            }

            const auto& exist_first = exist_segs.front();
            const auto& exist_last  = exist_segs.back();
            const auto& leaf_first = leaf_segs.front();
            const auto& leaf_last  = leaf_segs.back();

            uint32_t exist_first_beg = SegReplace::Interval::beg(exist_first);
            uint32_t exist_first_end = SegReplace::Interval::end(exist_first);
            uint32_t exist_last_beg  = SegReplace::Interval::beg(exist_last);
            uint32_t exist_last_end  = SegReplace::Interval::end(exist_last);
            uint32_t leaf_first_beg  = SegReplace::Interval::beg(leaf_first);
            uint32_t leaf_first_end  = SegReplace::Interval::end(leaf_first);
            uint32_t leaf_last_beg   = SegReplace::Interval::beg(leaf_last);
            uint32_t leaf_last_end   = SegReplace::Interval::end(leaf_last);
            uint32_t exist_beg       = std::min({exist_first_beg, exist_first_end, exist_last_beg, exist_last_end});
            uint32_t exist_end       = std::max({exist_first_beg, exist_first_end, exist_last_beg, exist_last_end});
            uint32_t leaf_beg        = std::min({leaf_first_beg, leaf_first_end, leaf_last_beg, leaf_last_end});
            uint32_t leaf_end        = std::max({leaf_first_beg, leaf_first_end, leaf_last_beg, leaf_last_end});
            bool     exist_is_rev    = SegReplace::Interval::is_reverse(exist_first);
            bool     leaf_is_rev1    = SegReplace::Interval::is_reverse(leaf_first);
            uint64_t exist_seg_id    = SegReplace::Interval::seg_id(exist_first);

            // Avoid self-mapping, But when the size of exist_segs > 1, this check is not sufficient. (2025-12-05)
            // Can be packed into a function later if needed.
            std::vector<std::string> leaf_base_names, exist_base_names;
            gfaName::collect_roots_from_name(leaf_name, leaf_base_names);
            for (const auto& seg : exist_segs) {
                std::string exist_name = nodes_[SegReplace::Interval::seg_id(seg)].name;
                std::vector<std::string> exist_names_tmp;
                gfaName::collect_roots_from_name(exist_name, exist_names_tmp);
                exist_base_names.insert(exist_base_names.end(), exist_names_tmp.begin(), exist_names_tmp.end());
            }
            if (std::any_of(leaf_base_names.begin(), leaf_base_names.end(),
                            [&](const std::string& s) {
                                return std::find(exist_base_names.begin(),
                                                exist_base_names.end(), s)
                                    != exist_base_names.end();
                            })) {
                continue;
            }

            auto& exist_cuts = cuts_[exist_seg_id].v;

            auto it_exist_1 = std::lower_bound(exist_cuts.begin(), exist_cuts.end(), exist_beg);
            auto it_exist_2 = std::upper_bound(exist_cuts.begin(), exist_cuts.end(), exist_end);
            auto it_leaf_1  = std::lower_bound(leaf_cuts.begin(),  leaf_cuts.end(),  leaf_beg);
            auto it_leaf_2  = std::upper_bound(leaf_cuts.begin(),  leaf_cuts.end(),  leaf_end);
            auto exist_diff = it_exist_2 - it_exist_1;
            auto leaf_diff  = it_leaf_2  - it_leaf_1;
            std::vector<uint32_t> cut_in_overlap_exist(it_exist_1, it_exist_2);
            std::vector<uint32_t> cut_in_overlap_leaf(it_leaf_1, it_leaf_2);

            std::vector<uint32_t> new_cuts_exist = cut_in_overlap_exist;
            std::vector<uint32_t> new_cuts_leaf  = cut_in_overlap_leaf;

            // exist -> leaf
            for (uint32_t cut : cut_in_overlap_exist) {
                uint64_t off = offset_from(cut, exist_beg, exist_end, exist_is_rev);
                uint32_t mapped = offset_to(off, leaf_beg, leaf_end, leaf_is_rev1);
                if (mapped < leaf_beg || mapped > leaf_end) continue;
                new_cuts_leaf.push_back(mapped);
            }
            std::sort(new_cuts_leaf.begin(), new_cuts_leaf.end());
            new_cuts_leaf.erase(std::unique(new_cuts_leaf.begin(), new_cuts_leaf.end()), new_cuts_leaf.end());

            // leaf -> exist
            for (uint32_t cut : cut_in_overlap_leaf) {
                uint64_t off = offset_from(cut, leaf_beg, leaf_end, leaf_is_rev1);
                uint32_t mapped = offset_to(off, exist_beg, exist_end, exist_is_rev);
                if (mapped < exist_beg || mapped > exist_end) continue;
                new_cuts_exist.push_back(mapped);
            }
            std::sort(new_cuts_exist.begin(), new_cuts_exist.end());
            new_cuts_exist.erase(std::unique(new_cuts_exist.begin(), new_cuts_exist.end()), new_cuts_exist.end());

            // The original cuts are all forward, now after adding the cutting points
            // they need to be reversed according to their respective directions
            if (exist_is_rev) {
                std::reverse(cut_in_overlap_exist.begin(), cut_in_overlap_exist.end());
                std::reverse(new_cuts_exist.begin(),  new_cuts_exist.end());
            }
            if (leaf_is_rev1) {
                std::reverse(cut_in_overlap_leaf.begin(), cut_in_overlap_leaf.end());
                std::reverse(new_cuts_leaf.begin(), new_cuts_leaf.end());
            }

            // calculate offsets
            std::vector<int32_t> new_exist_offsets = build_offsets(new_cuts_exist, exist_beg, exist_end, exist_is_rev);
            std::vector<int32_t> new_leaf_offsets  = build_offsets(new_cuts_leaf,  leaf_beg,  leaf_end,  leaf_is_rev1);

            if (DEBUG_ENABLED && (exist_diff == 1 || leaf_diff == 1)) {  // debug 2025-10-10
                std::cerr << "Skip: insufficient cut points in overlap window.\n";
                std::cerr << "   - Trunk: " << trunk_name << " (seg=" << exist_seg_id
                          << ", rev=" << exist_is_rev << ") window=[" << exist_beg << "," << exist_end
                          << "] #cuts=" << exist_diff << "\n   - cuts: ";
                for (auto c : cut_in_overlap_exist) std::cerr << c << " ";
                std::cerr << "\n";

                std::cerr << "   - Leaf : " << leaf_name  << " (seg=" << leaf_seg_id
                          << ", rev=" << leaf_is_rev1 << ") window=[" << leaf_beg << "," << leaf_end
                          << "] #cuts=" << leaf_diff << "\n   - cuts: ";
                for (auto c : cut_in_overlap_leaf)  std::cerr << c << " ";
                std::cerr << "\n";
                continue;
            }

            if (DEBUG_ENABLED && (exist_diff > 2 || leaf_diff > 2)) {
                std::string log_tmp = "   - exist segs: ";
                for (auto& seg : exist_segs) log_tmp += SegReplace::Interval::format(seg) + " ";
                log_tmp += "\n";
                log_tmp += "   - leaf_segs: ";
                for (auto& seg : leaf_segs) log_tmp += SegReplace::Interval::format(seg) + " ";
                log_tmp += "\n";
                std::cerr << log_tmp;

                std::cerr << "   - Existing segment cuts(" << exist_seg_id << "-" << exist_diff << "): ";
                for (auto cut : cut_in_overlap_exist) std::cerr << cut << " ";
                std::cerr << "\n";
                std::cerr << "   - Leaf segment cuts(" << leaf_seg_id << "-" << leaf_diff << "): ";
                for (auto cut : cut_in_overlap_leaf) std::cerr << cut << " ";
                std::cerr << "\n";

                std::cerr << "     - exist original\n";
                for (const auto& cut : cut_in_overlap_exist) std::cerr << cut << " ";
                std::cerr << "\n";
                std::cerr << "     - leaf original\n";
                for (const auto& cut : cut_in_overlap_leaf) std::cerr << cut << " ";
                std::cerr << "\n";
                std::cerr << "     - exist new\n";
                for (const auto& cut : new_cuts_exist) std::cerr << cut << " ";
                std::cerr << "\n";
                std::cerr << "     - leaf new\n";
                for (const auto& cut : new_cuts_leaf) std::cerr << cut << " ";
                std::cerr << "\n";
            }

            SegReplace::Seg exist_u128 = SegReplace::Interval::pack(exist_seg_id, exist_beg, exist_end, exist_is_rev);
            SegReplace::Seg leaf_u128  = SegReplace::Interval::pack(leaf_seg_id,  leaf_beg,  leaf_end,  leaf_is_rev1);

            bool leaf_is_exist;
            if (!prefer_v_as_leaf_(exist_seg_id, leaf_seg_id, leaf_is_exist)) continue;

            if (leaf_is_exist) {
                merge_and_dedup_cuts_(exist_cuts, new_cuts_exist);  // debug 2025-10-10, Update cut list

                build_rules_from_pair_windows_(
                    rulemap,
                    /*trunk*/ leaf_seg_id, leaf_is_rev1, cut_in_overlap_leaf,
                    /*leaf*/ exist_seg_id, exist_is_rev, new_cuts_exist, new_exist_offsets,
                    /*leaf older*/ cut_in_overlap_exist
                );
            } else {
                merge_and_dedup_cuts_(leaf_cuts, new_cuts_leaf);  // debug 2025-10-10, Update cut list

                build_rules_from_pair_windows_(
                    rulemap,
                    /*trunk*/ exist_seg_id, exist_is_rev, cut_in_overlap_exist,
                    /*leaf*/ leaf_seg_id, leaf_is_rev1, new_cuts_leaf, new_leaf_offsets,
                    /*leaf older*/ cut_in_overlap_leaf
                );
            }
        }

        if (DEBUG_ENABLED) {
            std::cerr << "   - Rule added: " << SegReplace::Interval::format(trunk_u128) << " -> ";
            for (const auto& leaf_seg_u : leaf_segs) {
                std::cerr << SegReplace::Interval::format(leaf_seg_u) << " ";
            }
            std::cerr << "\n";
        }
    }
}

void GfaDeoverlapper::build_rulemap_() {
    log_stream() << "Building replacement rules from alignments ...\n";

    rulemap_.clear();

    // Progress tracker
    ProgressTracker prog(bubble_aligns_.size());

    for (const auto& align : bubble_aligns_) {
        prog.hit();  // Update progress

        const uint32_t seg_a = NodeHandle::get_segment_id(align.v_a);
        const uint32_t seg_b = NodeHandle::get_segment_id(align.v_b);
        const bool     rev_a = NodeHandle::get_is_reverse(align.v_a);
        const bool     rev_b = NodeHandle::get_is_reverse(align.v_b);
        if (seg_a >= cuts_.size() || seg_b >= cuts_.size()) continue;
        if (seg_a == seg_b) continue;

        const std::string& name_a = nodes_[seg_a].name;
        const std::string& name_b = nodes_[seg_b].name;

        const uint32_t len_a = nodes_[seg_a].length;
        const uint32_t len_b = nodes_[seg_b].length;

        if (!(align.beg_a < align.end_a && align.end_a <= len_a)) continue;
        if (!(align.beg_b < align.end_b && align.end_b <= len_b)) continue;

        uint32_t pos_a = rev_a ? align.end_a : align.beg_a;
        uint32_t pos_b = rev_b ? align.end_b : align.beg_b;

        auto step = [](uint32_t& p, bool rev, uint32_t n) {
            if (n == 0) return;
            p = rev ? (p - n) : (p + n);
        };

        for (const auto& op : align.ops) {
            switch (op.op) {
                case 'M':
                case '=': {
                    const uint32_t prev_a = pos_a;
                    const uint32_t prev_b = pos_b;

                    step(pos_a, rev_a, op.len);
                    step(pos_b, rev_b, op.len);

                    const auto [beg_a, end_a] = rev_a
                        ? std::make_pair(pos_a, prev_a)
                        : std::make_pair(prev_a, pos_a);
                    const auto [beg_b, end_b] = rev_b
                        ? std::make_pair(pos_b, prev_b)
                        : std::make_pair(prev_b, pos_b);

                    if (op.len >= MIN_EQ_FOR_CUT_) {
                        GfaDeoverlapper::record_rulemap(
                            seg_a, seg_b,
                            rev_a, rev_b,
                            name_a, name_b,
                            beg_a, end_a, beg_b, end_b,
                            rulemap_
                        );
                    }
                } break;

                case 'X': // mismatch
                    step(pos_a, rev_a, op.len);
                    step(pos_b, rev_b, op.len);
                    break;

                case 'I': // insertion
                    step(pos_b, rev_b, op.len);
                    break;

                case 'D': // deletion
                    step(pos_a, rev_a, op.len);
                    break;

                case 'S': // soft clip
                case 'H': // hard clip
                    break;

                default:
                    warning_stream() << "   - Unknown CIGAR op: " << op.op << "\n";
                    break;
            }

            if (pos_a > len_a || pos_b > len_b) {
                warning_stream() 
                    << "   - CIGAR stepping overflow: pos_a=" << pos_a
                    << " len_a=" << len_a << " pos_b=" << pos_b
                    << " len_b=" << len_b << "\n";
                break;
            }
        }
    }
}

SegReplace::Expander GfaDeoverlapper::build_SegReplace_(bool filter_abnormal)
{
    SegReplace::Expander ex(rulemap_, getAllSegmentNames());
    ex.build_index();
    const SegReplace::RuleMap& idx = ex.index_view();
    
    if (filter_abnormal) {
        auto seg_less = [](const SegReplace::Seg& a, const SegReplace::Seg& b) {
            return std::tuple<uint64_t, uint32_t, uint32_t, bool>(SegReplace::Interval::seg_id(a), SegReplace::Interval::beg(a), SegReplace::Interval::end(a), SegReplace::Interval::is_reverse(a)) <
                std::tuple<uint64_t, uint32_t, uint32_t, bool>(SegReplace::Interval::seg_id(b), SegReplace::Interval::beg(b), SegReplace::Interval::end(b), SegReplace::Interval::is_reverse(b));
        };

        bool changed = true;
        while (changed) {
            changed = false;

            // Sort keys
            std::vector< SegReplace::Seg> keys;
            keys.reserve(idx.size());
            for (const auto& kv : idx) keys.push_back(kv.first);
            std::sort(keys.begin(), keys.end(), seg_less);

            std::vector<SegReplace::Seg> to_del;
            to_del.reserve(keys.size() / 8 + 8);

            for (size_t i = 1; i < keys.size(); ++i) {
                const SegReplace::Seg& prev = keys[i - 1];
                const SegReplace::Seg& cur  = keys[i];
                const bool prev_is_rev = SegReplace::Interval::is_reverse(prev);
                const bool cur_is_rev  = SegReplace::Interval::is_reverse(cur);

                // seg_id and adjacency check
                if (SegReplace::Interval::seg_id(prev) != SegReplace::Interval::seg_id(cur)) continue;
                if (SegReplace::Interval::end(prev) != SegReplace::Interval::beg(cur)) continue;

                auto it0 = idx.find(prev);
                auto it1 = idx.find(cur);
                if (it0 == idx.end() || it1 == idx.end()) continue;

                auto touches = [&](const SegReplace::Expansion& prev_v, const SegReplace::Expansion& cur_v) -> bool {
                    if (prev_v.empty() || cur_v.empty()) return false;
                    const SegReplace::Seg& a = prev_is_rev ? SegReplace::Interval::toggle_strand(prev_v.front()) : prev_v.back();
                    const SegReplace::Seg& b = cur_is_rev  ? SegReplace::Interval::toggle_strand(cur_v.back())   : cur_v.front();
                    return a == b;
                };

                const SegReplace::Expansion& v0 = it0->second;
                const SegReplace::Expansion& v1 = it1->second;

                if (!touches(v0, v1)) continue;

                const uint64_t L0 = v0.size();
                const uint64_t L1 = v1.size();

                bool del_prev = false;
                if (L0 != L1) del_prev = (L0 < L1);
                else {
                    const uint32_t b0 = SegReplace::Interval::beg(prev), e0 = SegReplace::Interval::end(prev);
                    const uint32_t b1 = SegReplace::Interval::beg(cur),  e1 = SegReplace::Interval::end(cur);
                    const uint32_t K0 = (e0 > b0) ? (e0 - b0) : 0;
                    const uint32_t K1 = (e1 > b1) ? (e1 - b1) : 0;
                    del_prev = (K0 < K1);
                }

                to_del.push_back(del_prev ? prev : cur);
            }

            // Delete marked keys
            if (!to_del.empty()) {
                std::sort(to_del.begin(), to_del.end(), seg_less);
                to_del.erase(std::unique(to_del.begin(), to_del.end()), to_del.end());

                for (const auto& k : to_del) ex.remove_from_index_by_key(k);
                changed = true;
            }
        }
    }

    if (DEBUG_ENABLED) ex.print_index();
    if (DEBUG_ENABLED) rulemap_verify(nodes_, ex.index_view());

    return ex;
}


void GfaDeoverlapper::rulemap_verify(const std::vector<GfaNode>& nodes, const SegReplace::RuleMap& idx) {
    wfa::WFAlignerGapAffine aligner(
        /*match*/ 0,
        /*gap_open*/ 4,
        /*gap_extend*/ 2,
        /*mismatch*/ 4,
        wfa::WFAligner::Alignment,
        wfa::WFAligner::MemoryHigh
    );

    std::vector<std::pair<SegReplace::Seg, SegReplace::Expansion>> items(idx.begin(), idx.end());
    auto seg_less = [](const std::pair<SegReplace::Seg, SegReplace::Expansion>& a,
                       const std::pair<SegReplace::Seg, SegReplace::Expansion>& b) {
        return std::tuple<uint64_t, uint32_t, uint32_t, bool>(
                SegReplace::Interval::seg_id(a.first),
                SegReplace::Interval::beg(a.first),
                SegReplace::Interval::end(a.first),
                SegReplace::Interval::is_reverse(a.first)
            ) < std::tuple<uint64_t, uint32_t, uint32_t, bool>(
                SegReplace::Interval::seg_id(b.first),
                SegReplace::Interval::beg(b.first),
                SegReplace::Interval::end(b.first),
                SegReplace::Interval::is_reverse(b.first)
            );
    };
    std::sort(items.begin(), items.end(), seg_less);

    auto slice = [&](SegReplace::Seg s)->std::string{
        return slice_seq_or_star_(
            nodes,
            (uint32_t)SegReplace::Interval::seg_id(s),
            SegReplace::Interval::beg(s),
            SegReplace::Interval::end(s),
            SegReplace::Interval::is_reverse(s)
        );
    };

    auto getNodeName = [&](SegReplace::Interval::u128 interId)->std::string {
        uint64_t seg_id = SegReplace::Interval::seg_id(interId);
        if (seg_id < nodes.size()) return nodes[seg_id].name + ":" + std::to_string(SegReplace::Interval::beg(interId)) + "-" + std::to_string(SegReplace::Interval::end(interId)) + (SegReplace::Interval::is_reverse(interId) ? "-" : "+");
        else return std::to_string(seg_id);
    };

    for (const auto& kv : items) {
        std::string first_seq, second_seq;
        for (const auto& x : kv.second) {
            if (first_seq.empty()) {
                first_seq = slice(kv.first);
            }
            second_seq += slice(x);
        }
        if (first_seq != second_seq) {
            std::string log;
            log += getNodeName(kv.first) + " =>";
            for (const auto& x : kv.second) log += " " + getNodeName(x);
            log += "\n";
            debug_stream() << log;
            debug_stream() << "F: " << first_seq << "\n";
            debug_stream() << "S: " << second_seq << "\n";
            aligner.alignEnd2End(first_seq.data(), (int)first_seq.size(), second_seq.data(), (int)second_seq.size());
            const std::string cigar = aligner.getCIGAR(true);
            debug_stream() << "CIGAR: " << cigar << "\n";
        }
    }
}

void GfaDeoverlapper::refine_chain_(SegReplace::Expansion& chain, const SegReplace::Expander& ex) const {
    if (chain.size() < 2) return;

    auto query = [&](SegReplace::Seg big)->SegReplace::Expansion {
        return ex.query(big);
    };

    SegReplace::Expansion refined;
    refined.reserve(chain.size());

    size_t i = 0;
    while (i < chain.size()) {
        size_t j = i;
        size_t k = i;
        uint64_t chr = SegReplace::Interval::seg_id(chain[i]);
        bool     st  = SegReplace::Interval::is_reverse(chain[i]);
        uint32_t low = SegReplace::Interval::beg(chain[i]);
        uint32_t high= SegReplace::Interval::end(chain[i]);

        SegReplace::Expansion refined_tmp;
        while (j + 1 < chain.size() && touchable_(chain[j], chain[j+1])) {
            ++j;
            uint32_t nb = SegReplace::Interval::beg(chain[j]);
            uint32_t ne = SegReplace::Interval::end(chain[j]);
            if (low == ne) low = nb;
            else if (high == nb) high = ne;
            else break;

            SegReplace::Seg big = SegReplace::Interval::pack(chr, low, high, st);
            const auto e = query(big);
            if (!is_identity_(e, big)) {
                refined_tmp = e;
                k = j;
                break;
            }
        }
        if (!refined_tmp.empty()) {
            refined.insert(refined.end(), refined_tmp.begin(), refined_tmp.end());
            i = k + 1;
        } else {
            refined.push_back(chain[i]);
            i++;
        }
    }

    chain.swap(refined);
}

void GfaDeoverlapper::emit_node_expansion_edges_(
    const SegReplace::Expansion& node_expansion,
    std::unordered_set<uint64_t>& seen_edges,
    gfaName& namer,
    const std::string& v_name,
    const std::string& w_name,
    bool v_is_rev,
    bool w_is_rev,
    std::optional<std::pair<uint32_t,uint32_t>> v_span,
    std::optional<std::pair<uint32_t,uint32_t>> w_span
) {
    if (v_span && w_span) {
        debug_stream() << "Add arc: " << v_name << "("
                       << v_span->first << "-" << v_span->second << ":"
                       << (v_is_rev ? "-" : "+") << ") and "
                       << w_name << "("
                       << w_span->first << "-" << w_span->second << ":"
                       << (w_is_rev ? "-" : "+") << ")\n";
    } else {
        debug_stream() << "Add arc: " << v_name << "(" << (v_is_rev ? "-" : "+") << ") and " << w_name << "(" << (w_is_rev ? "-" : "+") << ")\n";
    }

    auto pretty_name = [&](SegReplace::Seg s) -> std::string {
        const std::string& base = getNodeName(SegReplace::Interval::seg_id(s));
        uint32_t beg = SegReplace::Interval::beg(s);
        uint32_t end = SegReplace::Interval::end(s);
        return namer.format_interval_name(base, beg, end, false);
    };

    for (size_t i = 0; i + 1 < node_expansion.size(); ++i) {
        std::string v_name_tmp = pretty_name(node_expansion[i]);
        std::string w_name_tmp = pretty_name(node_expansion[i + 1]);

        bool v_rev = SegReplace::Interval::is_reverse(node_expansion[i]);
        uint32_t v_sid = get_or_add_segment(v_name_tmp);
        std::string v_seg_seq = slice_seq_or_star_(
            nodes_,
            SegReplace::Interval::seg_id(node_expansion[i]),
            SegReplace::Interval::beg(node_expansion[i]),
            SegReplace::Interval::end(node_expansion[i]),
            false
        );
        GfaNode& v_n = nodes_[v_sid];
        v_n.name = v_name_tmp;
        v_n.sequence = v_seg_seq;
        v_n.length = v_seg_seq != "*" ? (uint32_t)v_seg_seq.length() : 0;
        total_segment_length_ += v_seg_seq.length();
        uint32_t v_vid = (v_sid << 1) | (v_rev ? 1 : 0);

        uint32_t w_sid = get_or_add_segment(w_name_tmp);
        bool w_rev = SegReplace::Interval::is_reverse(node_expansion[i + 1]);
        std::string w_seg_seq = slice_seq_or_star_(
            nodes_,
            SegReplace::Interval::seg_id(node_expansion[i + 1]),
            SegReplace::Interval::beg(node_expansion[i + 1]),
            SegReplace::Interval::end(node_expansion[i + 1]),
            false
        );
        GfaNode& w_n = nodes_[w_sid];
        w_n.name = w_name_tmp;
        w_n.sequence = w_seg_seq;
        w_n.length = w_seg_seq != "*" ? (uint32_t)w_seg_seq.length() : 0;
        total_segment_length_ += w_seg_seq.length();
        uint32_t w_vid = (w_sid << 1) | (w_rev ? 1 : 0);

        seen_edges.insert(encode_edge_u64_(v_vid, w_vid));

        debug_stream() << "   - " << v_name_tmp << (v_rev ? "-" : "+") << " -> " << w_name_tmp << (w_rev ? "-" : "+") << "\n";
    }
    debug_stream() << "\n";
}

void GfaDeoverlapper::expand_and_rewire_edges_(
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
        uint32_t v_length = nodes_[v_seg_id].length, w_length = nodes_[w_seg_id].length;

        if (arc.ov != arc.ow) {
            error_stream() << "Overlap values are not equal for " << v_name << " and " << w_name;
            std::exit(1);
        }
        auto [vb, ve] = v_overlap_pos_(v_length, arc.ov, v_is_rev);
        auto [wb, we] = w_overlap_pos_(w_length, arc.ow, w_is_rev);

        auto& v_cuts = const_cast<std::vector<uint32_t>&>(cuts_[v_seg_id].v);
        auto& w_cuts = const_cast<std::vector<uint32_t>&>(cuts_[w_seg_id].v);

        auto it_v_1 = std::lower_bound(v_cuts.begin(), v_cuts.end(), vb);
        auto it_v_2 = std::upper_bound(v_cuts.begin(), v_cuts.end(), ve);
        auto it_w_1 = std::lower_bound(w_cuts.begin(), w_cuts.end(), wb);
        auto it_w_2 = std::upper_bound(w_cuts.begin(), w_cuts.end(), we);
        std::vector<uint32_t> cut_in_overlap_v(it_v_1, it_v_2);
        std::vector<uint32_t> cut_in_overlap_w(it_w_1, it_w_2);

        // outside overlap
        std::vector<uint32_t> cut_outside_v, cut_outside_w;
        cut_outside_v.insert(cut_outside_v.end(), v_cuts.begin(), it_v_1);
        cut_outside_v.insert(cut_outside_v.end(), it_v_2, v_cuts.end());
        cut_outside_w.insert(cut_outside_w.end(), w_cuts.begin(), it_w_1);
        cut_outside_w.insert(cut_outside_w.end(), it_w_2, w_cuts.end());

        bool v_is_leaf;
        if (!prefer_v_as_leaf_(v_seg_id, w_seg_id, v_is_leaf)) continue;

        auto&        trunk_cuts          = v_is_leaf ? w_cuts : v_cuts;
        auto&        leaf_cuts           = v_is_leaf ? v_cuts : w_cuts;
        auto&        trunk_in_overlap_v  = v_is_leaf ? cut_in_overlap_w : cut_in_overlap_v;
        auto&        leaf_in_overlap_v   = v_is_leaf ? cut_in_overlap_v : cut_in_overlap_w;
        bool         trunk_is_rev        = v_is_leaf ? w_is_rev : v_is_rev;
        bool         leaf_is_rev         = v_is_leaf ? v_is_rev : w_is_rev;
        uint32_t     trunk_seg_id        = v_is_leaf ? w_seg_id : v_seg_id;
        uint32_t     leaf_seg_id         = v_is_leaf ? v_seg_id : w_seg_id;
        const auto&  trunk_name          = v_is_leaf ? w_name : v_name;
        const auto&  leaf_name           = v_is_leaf ? v_name : w_name;
        uint32_t     tb                  = v_is_leaf ? wb : vb;
        uint32_t     te                  = v_is_leaf ? we : ve;
        uint32_t     lb                  = v_is_leaf ? vb : wb;
        uint32_t     le                  = v_is_leaf ? ve : we;

        // ============= Collect the expansion of source/target/overlap =============
        auto collect_non_overlap = [&](
            uint32_t seg_id, bool seg_rev,
            const std::vector<uint32_t>& cuts_vec,
            uint32_t win_beg, uint32_t win_end,
            bool left_side /* true: Only the left side of the window, false: Only take the right side */
        ) -> std::vector<SegReplace::Expansion> {
            std::vector<SegReplace::Expansion> res;
            for (size_t i = 0; i + 1 < cuts_vec.size();) {
                uint32_t beg = cuts_vec[i], end = cuts_vec[i + 1];
                bool overlap_with_window =
                    (!seg_rev && ((left_side && end > win_beg) || (!left_side && beg < win_end))) ||
                    ( seg_rev && ((left_side && beg < win_end) || (!left_side && end > win_beg)));

                if (overlap_with_window) {
                    ++i;
                    continue;
                }

                SegReplace::Seg window_u128 = SegReplace::Interval::pack(seg_id, beg, end, seg_rev);
                const auto leafs_u128 = query(window_u128);

                size_t j = i + 1;
                size_t k = i + 1;
                if (!is_identity_(leafs_u128, window_u128)) {
                    res.push_back(leafs_u128);
                } else {
                    uint32_t acc_beg = beg;
                    uint32_t acc_end = end;
                    SegReplace::Expansion chosen = leafs_u128;
                    if (!seg_rev) {
                        while ((j + 1) < cuts_vec.size()) {
                            ++j;
                            acc_end = cuts_vec[j];
                            if (left_side && acc_end > win_beg) break;
                            SegReplace::Seg big_window_u128 = SegReplace::Interval::pack(seg_id, acc_beg, acc_end, seg_rev);
                            const auto big_leafs_u128 = query(big_window_u128);
                            if (!is_identity_(big_leafs_u128, big_window_u128)) {
                                chosen = big_leafs_u128;
                                k = j;
                            }
                        }
                    } else {
                        while ((j + 1) < cuts_vec.size()) {
                            ++j;
                            acc_end = cuts_vec[j];
                            if (!left_side && acc_end > win_beg) break;
                            SegReplace::Seg big_window_u128 = SegReplace::Interval::pack(seg_id, acc_beg, acc_end, seg_rev);
                            const auto big_leafs_u128 = query(big_window_u128);
                            if (!is_identity_(big_leafs_u128, big_window_u128)) {
                                chosen = big_leafs_u128;
                                k = j;
                            }
                        }
                    }
                    res.push_back(chosen);
                }
                i = k;
            }
            return res;
        };

        auto collect_overlap = [&](
            uint32_t seg_id, bool seg_rev, const std::vector<uint32_t>& vec
        ) -> std::vector<SegReplace::Expansion> {
            std::vector<SegReplace::Expansion> res;
            for (size_t i = 0; i + 1 < vec.size();) {
                uint32_t beg = vec[i], end = vec[i + 1];
                SegReplace::Seg window_u128 = SegReplace::Interval::pack(seg_id, beg, end, seg_rev);
                const auto leafs_u128 = query(window_u128);

                size_t j = i + 1;
                size_t k = i + 1;
                if (!is_identity_(leafs_u128, window_u128)) {
                    res.push_back(leafs_u128);
                } else {
                    uint32_t acc_beg = beg;
                    uint32_t acc_end = end;
                    SegReplace::Expansion chosen = leafs_u128;
                    while ((j + 1) < vec.size()) {
                        ++j;
                        acc_end = vec[j];
                        SegReplace::Seg big_window_u128 = SegReplace::Interval::pack(seg_id, acc_beg, acc_end, seg_rev);
                        const auto big_leafs_u128 = query(big_window_u128);
                        if (!is_identity_(big_leafs_u128, big_window_u128)) {
                            chosen = big_leafs_u128;
                            k = j;
                        }
                    }
                    res.push_back(chosen);
                }
                i = k;
            }
            return res;
        };

        // no-overlap
        std::vector<SegReplace::Expansion> source_expansions = collect_non_overlap(v_seg_id, v_is_rev, v_cuts, vb, ve, /*left_side=*/true);
        std::vector<SegReplace::Expansion> target_expansions = collect_non_overlap(w_seg_id, w_is_rev, w_cuts, wb, we, /*left_side=*/false);

        // overlap
        std::vector<SegReplace::Expansion> source_overlap_expansions = collect_overlap(v_seg_id, v_is_rev, cut_in_overlap_v);
        std::vector<SegReplace::Expansion> target_overlap_expansions = collect_overlap(w_seg_id, w_is_rev, cut_in_overlap_w);

        // ========================= chain base on source =========================
        {
            auto src = source_expansions;
            auto tgt = target_expansions;

            // insert the overlap expansions back to the corresponding side (direction decides front/back)
            if (v_is_rev) src.insert(src.begin(), source_overlap_expansions.begin(), source_overlap_expansions.end());
            else src.insert(src.end(), source_overlap_expansions.begin(), source_overlap_expansions.end());

            // inversion according to the edge direction
            if (v_is_rev) std::reverse(src.begin(), src.end());
            if (w_is_rev) std::reverse(tgt.begin(), tgt.end());

            // chaining the whole chain
            SegReplace::Expansion node_expansion_v;
            for (const auto& seq : src) node_expansion_v.insert(node_expansion_v.end(), seq.begin(), seq.end());
            for (const auto& seq : tgt) node_expansion_v.insert(node_expansion_v.end(), seq.begin(), seq.end());

            // refinement (cross-boundary merging and expansion)
            refine_chain_(node_expansion_v, ex);

            // emit edges (seen_edges)
            emit_node_expansion_edges_(
                node_expansion_v,
                seen_edges,
                namer,
                v_name, w_name,
                v_is_rev, w_is_rev,
                std::make_optional(std::pair<uint32_t,uint32_t>{vb, ve}),
                std::make_optional(std::pair<uint32_t,uint32_t>{wb, we})
            );
        }

        // ========================= chain base on target =========================
        {
            auto src = source_expansions;
            auto tgt = target_expansions;

            // insert the overlap expansions back to the corresponding side (direction decides front/back)
            if (w_is_rev) tgt.insert(tgt.end(), target_overlap_expansions.begin(), target_overlap_expansions.end());
            else tgt.insert(tgt.begin(), target_overlap_expansions.begin(), target_overlap_expansions.end());

            // inversion according to the edge direction
            if (v_is_rev) std::reverse(src.begin(), src.end());
            if (w_is_rev) std::reverse(tgt.begin(), tgt.end());

            // chaining the whole chain
            SegReplace::Expansion node_expansion_w;
            for (const auto& seq : src) node_expansion_w.insert(node_expansion_w.end(), seq.begin(), seq.end());
            for (const auto& seq : tgt) node_expansion_w.insert(node_expansion_w.end(), seq.begin(), seq.end());

            // refinement (cross-boundary merging and expansion)
            refine_chain_(node_expansion_w, ex);

            // emit edges (seen_edges)
            emit_node_expansion_edges_(
                node_expansion_w,
                seen_edges,
                namer,
                v_name, w_name,
                v_is_rev, w_is_rev,
                std::make_optional(std::pair<uint32_t,uint32_t>{vb, ve}),
                std::make_optional(std::pair<uint32_t,uint32_t>{wb, we})
            );
        }
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

void GfaDeoverlapper::remove_unused_nodes_() {
    log_stream() << "Removing unused nodes ...\n";

    std::vector<bool> used(nodes_.size(), false);

    for (const auto& e : arcs_) {
        if (e.get_del() || e.get_comp()) continue;
        used[e.get_source_segment_id()] = true;
        used[e.get_target_segment_id()] = true;
    }

    for (size_t i = 0; i < nodes_.size(); i++) {
        if (!used[i]) {
            delete_segment(i);
        }
    }

    rebuild_after_edits();
}

// Main process
void GfaDeoverlapper::deoverlap(const std::string& prefix) {
    log_stream() << "Normalizing overlapping GFA to non-overlapping GFA ...\n";
    // 0. Basic finishing
    finalize_();

    // 1. Handle inconsistencies in forward and reverse overlap and filter corresponding edges
    prune_overlaps_();

    // 2. Initialize cuts_: boundaries + window endpoints
    initialize_cuts_();

    // 3.1. Align overlaps
    overlaps_align_();

    // 3.2. Deduplicate alignments
    dedup_aligns_();

    // 4. Build cuts_ from alignments
    build_cuts_from_cigar_();

    // 5. Propagate cuts_ within windows
    propagate_cuts_();
    if (DEBUG_ENABLED) print_cuts(cuts_);

    // 6. Build rulemap_
    build_rulemap_();
    SegReplace::Expander ex = build_SegReplace_(true);
    ex.save_map(prefix + ".deoverlap.map");  // save rulemap

    // 7. Apply rules: expand and rewire 0M edges
    expand_and_rewire_edges_(ex);

    // 8. Remove unused nodes
    remove_unused_nodes_();

    log_stream() << "Finished normalizing to non-overlapping GFA.\n";
}