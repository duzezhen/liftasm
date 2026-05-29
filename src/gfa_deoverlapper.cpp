#include "../include/gfa_deoverlapper.hpp"
#include "../include/progress_tracker.hpp"

#include <set>
#include <tuple>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "ketopt.h"


// ------------------------------------------------ Helper functions ------------------------------------------------
static inline uint32_t find_index_by_offset(int32_t offset, const std::vector<int32_t>& offsets) {
    auto it = std::lower_bound(offsets.begin(), offsets.end(), offset);
    if (it != offsets.end() && *it == offset) {
        return static_cast<uint32_t>(it - offsets.begin());
    }
    return UINT32_MAX;
}

static inline uint64_t encode_edge_u64_(uint32_t v, uint32_t w) {
    return (uint64_t(v) << 32) | w;
};
static inline std::pair<uint32_t, uint32_t> decode_edge_u64_(uint64_t e) {
    uint32_t v = (uint32_t)(e >> 32);
    uint32_t w = (uint32_t)(e & 0xFFFFFFFF);
    return {v, w};
};

static inline uint32_t clamp32_(int64_t x, uint32_t L) {
    if (x < 0) return 0u;
    if ((uint64_t)x > (uint64_t)L) return L;
    return (uint32_t)x;
}

// Check if an expansion is an identity mapping
static inline bool is_identity_(const SegReplace::Expansion& e, SegReplace::Seg k) {
    return e.size() == 1 && e[0] == k;
}

// touchable: same chr & same direction, and end-to-end
static inline bool touchable_(SegReplace::Seg a, SegReplace::Seg b) {
    if (SegReplace::Interval::seg_id(a) != SegReplace::Interval::seg_id(b)) return false;
    if (SegReplace::Interval::is_reverse(a) != SegReplace::Interval::is_reverse(b)) return false;
    uint32_t ab = SegReplace::Interval::beg(a), ae = SegReplace::Interval::end(a);
    uint32_t bb = SegReplace::Interval::beg(b), be = SegReplace::Interval::end(b);
    return (ae == bb) || (be == ab);
}

// On segment v, take the overlap window [vb, ve) based on the direction
static inline std::pair<uint32_t,uint32_t> v_overlap_pos_(uint32_t Lseg, uint32_t OV, bool rev) {
    if (!rev) return {Lseg - OV, Lseg};
    else      return {0u, OV};
}

// On segment w, take the overlap window [wb, we) based on the direction
static inline std::pair<uint32_t,uint32_t> w_overlap_pos_(uint32_t Lseg, uint32_t OW, bool rev) {
    if (!rev) return {0u, OW};
    else      return {Lseg - OW, Lseg};
}


// ------------------------------------------------ GFA deoverlapper ------------------------------------------------
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
        warning_stream() << "  ! Total edges with unequal overlap: " << mismatches.size() << "\n\n";
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
        log_stream() << "  - No alignment produced for " << v_name << " vs " << w_name << "\n";
        return hits;
    }

    hits.reserve(aligns.size());
    for (auto& a : aligns) {
        MmWfaHit h;
        h.r_beg = (uint32_t)a.r_beg;
        h.r_end = (uint32_t)a.r_end;
        h.q_beg = (uint32_t)a.q_beg;
        h.q_end = (uint32_t)a.q_end;

        if (CIGAR::match_ratio(a.cigar) < MIN_MATCH_RATIO_) continue;

        h.cigar = std::move(a.cigar);
        hits.emplace_back(std::move(h));
    }

    return hits;
}


std::vector<GfaDeoverlapper::MmWfaHit> GfaDeoverlapper::align_short_mm2_(
    const std::string& v_name,
    const std::string& w_name,
    const std::string& v_seq_slice,
    const std::string& w_seq_slice
) {
    std::vector<MmWfaHit> hits;

    if (v_seq_slice.empty() || w_seq_slice.empty()) return hits;

    if (v_seq_slice == w_seq_slice) {
        MmWfaHit h;
        h.r_beg = 0;
        h.r_end = (uint32_t)v_seq_slice.size();
        h.q_beg = 0;
        h.q_end = (uint32_t)w_seq_slice.size();
        h.cigar = std::to_string(v_seq_slice.size()) + "=";
        hits.emplace_back(std::move(h));
        return hits;
    }

    const size_t min_len = std::min(v_seq_slice.size(), w_seq_slice.size());
    if (min_len < 4) return hits;

    mm_idxopt_t ipt;
    mm_mapopt_t opt;

    mm_set_opt(nullptr, &ipt, &opt);
    mm_set_opt("asm5", &ipt, &opt);

    ipt.k = (short)std::min<size_t>(std::max<size_t>(4, min_len), 11);
    ipt.w = 1;

    opt.flag |= MM_F_CIGAR;
    opt.flag |= MM_F_EQX;
    opt.best_n = (short)anchorOpts_.max_kept;
    opt.zdrop = extendOpts_.dyn_zdrop;

    const char* ref_seqs[1] = { v_seq_slice.c_str() };
    const char* ref_names[1] = { v_name.c_str() };

    mm_idx_t* mi = mm_idx_str(ipt.w, ipt.k, 0, ipt.bucket_bits, 1, ref_seqs, ref_names);
    if (!mi) return hits;

    mm_mapopt_update(&opt, mi);

    mm_tbuf_t* tbuf = mm_tbuf_init();

    int n_regs = 0;
    mm_reg1_t* regs = mm_map(
        mi,
        (int)w_seq_slice.size(),
        w_seq_slice.c_str(),
        &n_regs,
        tbuf,
        &opt,
        w_name.c_str()
    );

    for (int i = 0; i < n_regs; ++i) {
        const mm_reg1_t& r = regs[i];
        if (!r.p || r.rev) continue;

        const mm_extra_t* ex = (const mm_extra_t*)r.p;
        if (!ex || ex->n_cigar == 0) continue;

        std::string cigar = CIGAR::to_string(ex->cigar, ex->n_cigar);
        if (cigar.empty() || cigar == "*") continue;
        if (CIGAR::match_ratio(cigar) < MIN_MATCH_RATIO_) continue;

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

    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);

    return hits;
}

struct TBufHolder {
    mm_tbuf_t* p = nullptr;

    ~TBufHolder() {
        if (p) mm_tbuf_destroy(p);
    }
};


std::vector<GfaDeoverlapper::MmWfaHit> GfaDeoverlapper::align_mm2_(
    const std::string& v_name,
    const std::string& w_name,
    const std::string& v_seq_slice,
    const std::string& w_seq_slice
) {
    std::vector<MmWfaHit> hits;

    if (v_seq_slice.empty() || w_seq_slice.empty()) {
        return hits;
    }

    mm_idxopt_t ipt;
    mm_mapopt_t opt;

    mm_set_opt(nullptr, &ipt, &opt);
    mm_set_opt("asm5", &ipt, &opt);

    ipt.k = (short)chainOpts_.k;
    ipt.w = (short)chainOpts_.w;

    opt.flag |= MM_F_CIGAR;
    opt.flag |= MM_F_EQX;  // output =/X

    opt.best_n = (short)anchorOpts_.max_kept;
    opt.zdrop = extendOpts_.dyn_zdrop;

    const char* ref_seqs[1]  = { v_seq_slice.c_str() };
    const char* ref_names[1] = { v_name.c_str() };

    const int is_hpc = (ipt.flag & MM_I_HPC) ? 1 : 0;
    mm_idx_t* mi = mm_idx_str(ipt.w, ipt.k, is_hpc, ipt.bucket_bits, 1, ref_seqs, ref_names);
    if (!mi) return hits;

    mm_mapopt_update(&opt, mi);

    static thread_local TBufHolder holder;
    if (!holder.p) holder.p = mm_tbuf_init();
    mm_tbuf_t* tbuf = holder.p;

    int n_regs = 0;
    mm_reg1_t* regs = mm_map(mi, (int)w_seq_slice.size(), w_seq_slice.c_str(), &n_regs, tbuf, &opt, w_name.c_str());

    for (int i = 0; i < n_regs; ++i) {
        const mm_reg1_t& r = regs[i];
        if (!r.p) continue;
        if (r.rev) continue;  // keep same strand only
        const mm_extra_t* ex = (const mm_extra_t*)r.p;
        if (!ex || ex->n_cigar == 0) continue;

        std::string cigar = CIGAR::to_string(ex->cigar, ex->n_cigar);
        if (cigar.empty() || cigar == "*") continue;

        if (CIGAR::match_ratio(cigar) < MIN_MATCH_RATIO_) continue;

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

    if (v_name == w_name) { return outs; }

    if (v_seq_slice.empty() || v_seq_slice == "*") {
        error_stream() << "Empty/unknown sequence for " << v_name << "\n";
        std::exit(1);
    }

    if (w_seq_slice.empty() || w_seq_slice == "*") {
        error_stream() << "Empty/unknown sequence for " << w_name << "\n";
        std::exit(1);
    }

    const bool v_rev = NodeHandle::get_is_reverse(v_vertex);
    const bool w_rev = NodeHandle::get_is_reverse(w_vertex);

    auto print_alignment = [&](
        uint32_t r_beg, uint32_t r_end,
        uint32_t q_beg, uint32_t q_end,
        std::string cigar,
        std::vector<CIGAR::COp>& ops
    ) {
        // Convert to original segment coordinates
        //   forward:  [vb + r_beg, vb + r_end)
        //   reverse:  [ve - r_end, ve - r_beg)
        const uint32_t beg_a = v_rev ? (ve - r_end) : (vb + r_beg);
        const uint32_t end_a = v_rev ? (ve - r_beg) : (vb + r_end);
        const uint32_t beg_b = w_rev ? (we - q_end) : (wb + q_beg);
        const uint32_t end_b = w_rev ? (we - q_beg) : (wb + q_end);

        debug_stream() << "    - ref: " << v_name << "(" << (v_rev ? "-" : "+") << ")\n";
        debug_stream() << "    - qry: " << w_name << "(" << (w_rev ? "-" : "+") << ")\n";
        debug_stream() << "      - [r_beg, r_end, r_len) = [" << beg_a << ", " << end_a << ", " << v_seq_slice.length() << ")\n";
        debug_stream() << "      - [q_beg, q_end, q_len) = [" << beg_b << ", " << end_b << ", " << w_seq_slice.length() << ")\n";
        debug_stream() << "      - CIGAR = " << cigar << "\n";
    };

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

        outs.emplace_back(
            v_name, w_name,
            beg_a, end_a,
            beg_b, end_b,
            v_vertex, w_vertex,
            std::move(ops)
        );
    };

    if (DEBUG_ENABLED) {
        debug_stream() << "----------------------------------------\n";
        debug_stream() << "Aligning " << v_name << " vs " << w_name << " ...\n";
    }

    if (v_seq_slice == w_seq_slice) {
        const uint32_t L = (uint32_t)v_seq_slice.size();
        if (L == 0) return outs;
        std::string cigar = std::to_string(L) + "=";
        auto ops = CIGAR::parse(cigar);
        if (DEBUG_ENABLED) print_alignment(0, L, 0, L, cigar, ops);
        add_alignment(0, L, 0, L, std::move(cigar), std::move(ops));
        return outs;
    }

    std::vector<MmWfaHit> hits;
    if (std::min(v_seq_slice.size(), w_seq_slice.size()) < std::max(chainOpts_.k, chainOpts_.w)) {
        hits = align_short_mm2_(v_name, w_name, v_seq_slice, w_seq_slice);
    } else {
        hits = use_wfa_ ? align_wfa_(v_name, w_name, v_seq_slice, w_seq_slice) : align_mm2_(v_name, w_name, v_seq_slice, w_seq_slice);
    }

    if (DEBUG_ENABLED && !hits.empty()) {
        debug_stream() << "  - Total raw alignments: " << hits.size() << "\n";
        for (auto& h : hits) {
            if (h.cigar.empty()) continue;
            auto ops = CIGAR::parse(h.cigar);
            if (DEBUG_ENABLED) print_alignment(h.r_beg, h.r_end, h.q_beg, h.q_end, h.cigar, ops);
        }
    }

    hits = filter_aligns_(std::move(hits));

    if (DEBUG_ENABLED && !hits.empty()) debug_stream() << "  - Total filtered alignments: " << hits.size() << "\n";
    for (auto& h : hits) {
        if (h.cigar.empty()) continue;
        auto ops = CIGAR::parse(h.cigar);
        if (DEBUG_ENABLED) print_alignment(h.r_beg, h.r_end, h.q_beg, h.q_end, h.cigar, ops);
        add_alignment(h.r_beg, h.r_end, h.q_beg, h.q_end, std::move(h.cigar), std::move(ops));
    }

    return outs;
}


void GfaDeoverlapper::initialize_cuts_() {
    log_stream() << "Initializing cut points ...\n" << "\n";

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

    // ---------------- Thread pool ----------------
    ThreadPool pool(alignOpts_.threads);
    std::vector<std::future<std::vector<BubbleAlignment>>> futs;
    futs.reserve(arcs_.size());

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

        std::string& v_name = nodes_[v_seg_id].name;
        std::string& w_name = nodes_[w_seg_id].name;

        uint32_t v_len = nodes_[v_seg_id].length;
        uint32_t w_len = nodes_[w_seg_id].length;

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
        uint64_t key = pair_key(v_seg_id, w_seg_id);
        if (!seen_pairs.insert(key).second) continue;

        std::string v_seq_slice = slice_seq_or_star_(nodes_, v_seg_id, vb, ve, v_rev);
        std::string w_seq_slice = slice_seq_or_star_(nodes_, w_seg_id, wb, we, w_rev);

        if (v_seq_slice == "*" || w_seq_slice == "*" || v_seq_slice.empty() || w_seq_slice.empty()) {
            warning_stream() << "  ! Overlap sequences unknown (\"*\"): "
                        << v_name << "(" << vb << "-" << ve << (v_rev?":-":":+") << ") vs "
                        << w_name << "(" << wb << "-" << we << (w_rev?":-":":+") << "), len=" << L
                        << "\n";
            continue;
        }

        futs.emplace_back(
            pool.submit([this, v_vertex, w_vertex, v_name = std::string(v_name), w_name = std::string(w_name), v_seq_slice = std::move(v_seq_slice), w_seq_slice = std::move(w_seq_slice), vb, ve, wb, we]() {
                return align_and_pack_(
                    v_vertex, w_vertex,
                    v_name, w_name,
                    v_seq_slice, w_seq_slice,
                    vb, ve, wb, we
                );
            })
        );
    }

    // Collect results
    for (auto &f : futs) {
        auto vec = f.get();
        if (vec.empty()) continue;

        bubble_aligns_.insert(
            bubble_aligns_.end(),
            std::make_move_iterator(vec.begin()),
            std::make_move_iterator(vec.end())
        );
    }

    pool.stop();

    log_stream() << "  - Total alignments produced: " << bubble_aligns_.size() << "\n" << "\n";
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
        int32_t     idx;       // Index in bubble_aligns_
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

                groups[sA].push_back(SpanItem{aln.idx, pc.root, lo, hi, len});
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

                groups[sB].push_back(SpanItem{aln.idx, pc.root, lo, hi, len});
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
                if (x.idx  != y.idx ) return x.idx  < y.idx;
                if (x.len  != y.len ) return x.len  > y.len;
                if (x.beg  != y.beg ) return x.beg  < y.beg;
                return x.end < y.end;
            }
        );

        std::unordered_map<std::string, std::vector<SpanItem>> chosen_by_base;
        chosen_by_base.reserve(8);

        for (const auto& s : spans) {
            if (!keep[s.idx]) continue;

            auto &chosen = chosen_by_base[s.base];
            bool conflict = false;

            for (const auto& c : chosen) {
                uint32_t ov = overlap(s.beg, s.end, c.beg, c.end);
                if (ov > 0) {
                    keep[s.idx] = 0;
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

    log_stream() << "After base-overlap dedup, alignments kept: " << bubble_aligns_.size() << " / " << n_aln << "\n" << "\n";
}


std::vector<std::vector<size_t>> GfaDeoverlapper::build_align_groups_() const {
    log_stream() << "Grouping alignments by segment overlaps ...\n";

    const uint32_t nseg = static_cast<uint32_t>(cuts_.size());

    std::vector<uint32_t> parent(nseg);
    std::vector<uint8_t> rank(nseg, 0);

    for (uint32_t i = 0; i < nseg; ++i) parent[i] = i;

    auto find_root = [&](uint32_t x) -> uint32_t {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    for (const auto& aln : bubble_aligns_) {
        const uint32_t seg_a = NodeHandle::get_segment_id(aln.v_a);
        const uint32_t seg_b = NodeHandle::get_segment_id(aln.v_b);
        if (seg_a >= nseg || seg_b >= nseg) continue;

        uint32_t ra = find_root(seg_a);
        uint32_t rb = find_root(seg_b);
        if (ra == rb) continue;

        if (rank[ra] < rank[rb]) std::swap(ra, rb);
        parent[rb] = ra;
        if (rank[ra] == rank[rb]) ++rank[ra];
    }

    std::vector<std::vector<size_t>> groups;
    std::unordered_map<uint32_t, size_t> group_of_root;
    group_of_root.reserve(nseg);

    for (size_t i = 0; i < bubble_aligns_.size(); ++i) {
        const auto& aln = bubble_aligns_[i];

        const uint32_t seg_a = NodeHandle::get_segment_id(aln.v_a);
        const uint32_t seg_b = NodeHandle::get_segment_id(aln.v_b);
        if (seg_a >= nseg || seg_b >= nseg) continue;

        const uint32_t root = find_root(seg_a);

        auto it = group_of_root.find(root);
        if (it == group_of_root.end()) {
            const size_t gid = groups.size();
            group_of_root.emplace(root, gid);
            groups.emplace_back();
            groups.back().push_back(i);
        } else {
            groups[it->second].push_back(i);
        }
    }

    log_stream() << "  - Total groups generated: " << groups.size() << "\n" << "\n";

    return groups;
}


void GfaDeoverlapper::normalize_all_cuts_() {
    for (uint32_t sid = 0; sid < cuts_.size(); ++sid) {
        auto& C = cuts_[sid].v;

        std::sort(C.begin(), C.end());
        C.erase(std::unique(C.begin(), C.end()), C.end());

        if (C.empty() || C.front() != 0u) {
            C.insert(C.begin(), 0u);
        }

        const uint32_t len = nodes_[sid].length;
        if (C.back() != len) {
            C.push_back(len);
        }
    }
}


bool GfaDeoverlapper::add_cuts_from_one_alignment_(const BubbleAlignment& align) {
    uint32_t seg_a = NodeHandle::get_segment_id(align.v_a);
    uint32_t seg_b = NodeHandle::get_segment_id(align.v_b);
    bool     rev_a = NodeHandle::get_is_reverse(align.v_a);
    bool     rev_b = NodeHandle::get_is_reverse(align.v_b);

    if (seg_a >= cuts_.size() || seg_b >= cuts_.size()) return false;

    uint32_t len_a = nodes_[seg_a].length;
    uint32_t len_b = nodes_[seg_b].length;

    if (!(align.beg_a < align.end_a && align.end_a <= len_a)) return false;
    if (!(align.beg_b < align.end_b && align.end_b <= len_b)) return false;

    uint32_t pos_a = rev_a ? align.end_a : align.beg_a;
    uint32_t pos_b = rev_b ? align.end_b : align.beg_b;

    const bool is_single_exact_match = (align.ops.size() == 1 && align.ops[0].op == '=');
    bool changed = false;

    auto add_cut = [&](uint32_t seg, uint32_t p, uint32_t len) {
        if (p > 0 && p < len) {
            cuts_[seg].v.push_back(p);
            changed = true;
        }
    };

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

                if (op.len >= MIN_EQ_FOR_CUT_ || is_single_exact_match) {
                    add_cut(seg_a, prev_a, len_a);
                    add_cut(seg_b, prev_b, len_b);
                    add_cut(seg_a, pos_a,  len_a);
                    add_cut(seg_b, pos_b,  len_b);
                }
            } break;

            case 'X':
                step(pos_a, rev_a, op.len);
                step(pos_b, rev_b, op.len);
                break;

            case 'I':
                step(pos_b, rev_b, op.len);
                break;

            case 'D':
                step(pos_a, rev_a, op.len);
                break;

            case 'S':
            case 'H':
                break;

            default:
                warning_stream() << "  ! Unknown CIGAR op: " << op.op << "\n";
                break;
        }

        if (pos_a > len_a || pos_b > len_b) {
            warning_stream() << "  ! CIGAR stepping overflow: pos_a=" << pos_a << " len_a=" << len_a << " pos_b=" << pos_b << " len_b=" << len_b << "\n";
            break;
        }
    }

    return changed;
}

void GfaDeoverlapper::normalize_group_cuts_(const std::vector<size_t>& group) {
    std::vector<uint32_t> segs;
    segs.reserve(group.size() * 2);

    for (const size_t aln_id : group) {
        const auto& aln = bubble_aligns_[aln_id];

        const uint32_t seg_a = NodeHandle::get_segment_id(aln.v_a);
        const uint32_t seg_b = NodeHandle::get_segment_id(aln.v_b);

        if (seg_a < cuts_.size()) segs.push_back(seg_a);
        if (seg_b < cuts_.size()) segs.push_back(seg_b);
    }

    std::sort(segs.begin(), segs.end());
    segs.erase(std::unique(segs.begin(), segs.end()), segs.end());

    for (uint32_t sid : segs) {
        auto& C = cuts_[sid].v;

        std::sort(C.begin(), C.end());
        C.erase(std::unique(C.begin(), C.end()), C.end());

        if (C.empty() || C.front() != 0u) C.insert(C.begin(), 0u);
        if (C.back() != nodes_[sid].length) C.push_back(nodes_[sid].length);
    }
}


uint64_t GfaDeoverlapper::propagate_cuts_run_(const std::vector<size_t>& group) {
    std::unordered_map<uint32_t, std::vector<size_t>> incident;
    incident.reserve(group.size() * 2 + 8);

    std::vector<uint32_t> dirty;
    dirty.reserve(group.size() * 2);

    for (const size_t aln_id : group) {
        const auto& aln = bubble_aligns_[aln_id];

        const uint32_t seg_a = NodeHandle::get_segment_id(aln.v_a);
        const uint32_t seg_b = NodeHandle::get_segment_id(aln.v_b);

        if (seg_a >= cuts_.size() || seg_b >= cuts_.size()) continue;

        incident[seg_a].push_back(aln_id);
        if (seg_b != seg_a) incident[seg_b].push_back(aln_id);

        dirty.push_back(seg_a);
        dirty.push_back(seg_b);
    }

    std::sort(dirty.begin(), dirty.end());
    dirty.erase(std::unique(dirty.begin(), dirty.end()), dirty.end());

    std::vector<uint32_t> next_dirty;
    std::vector<size_t> active;

    std::unordered_set<size_t> active_seen;
    std::unordered_set<uint32_t> dirty_seen;
    active_seen.reserve(group.size() + 8);
    dirty_seen.reserve(dirty.size() + 8);

    uint64_t total_changed = 0;

    auto step = [](uint32_t& p, bool rev, uint32_t n) {
        if (n == 0) return;
        p = rev ? (p - n) : (p + n);
    };

    for (int iter = 0; iter < MAX_PROPAGATION_ITERS_ && !dirty.empty(); ++iter) {
        active.clear();
        active_seen.clear();

        for (uint32_t sid : dirty) {
            auto it = incident.find(sid);
            if (it == incident.end()) continue;

            for (size_t aln_id : it->second) {
                if (active_seen.insert(aln_id).second) {
                    active.push_back(aln_id);
                }
            }
        }

        next_dirty.clear();
        dirty_seen.clear();

        for (size_t aln_id : active) {
            const auto& align = bubble_aligns_[aln_id];

            uint32_t seg_a = NodeHandle::get_segment_id(align.v_a);
            uint32_t seg_b = NodeHandle::get_segment_id(align.v_b);
            bool     rev_a = NodeHandle::get_is_reverse(align.v_a);
            bool     rev_b = NodeHandle::get_is_reverse(align.v_b);

            if (seg_a >= cuts_.size() || seg_b >= cuts_.size()) continue;

            uint32_t len_a = nodes_[seg_a].length;
            uint32_t len_b = nodes_[seg_b].length;

            if (!(align.beg_a < align.end_a && align.end_a <= len_a)) continue;
            if (!(align.beg_b < align.end_b && align.end_b <= len_b)) continue;

            uint32_t pos_a = rev_a ? align.end_a : align.beg_a;
            uint32_t pos_b = rev_b ? align.end_b : align.beg_b;

            auto& a_cuts = cuts_[seg_a].v;
            auto& b_cuts = cuts_[seg_b].v;

            std::vector<uint32_t> new_a_cuts;
            std::vector<uint32_t> new_b_cuts;

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

                        auto it_b_1 = std::lower_bound(b_cuts.begin(), b_cuts.end(), beg_b);
                        auto it_b_2 = std::upper_bound(b_cuts.begin(), b_cuts.end(), end_b);

                        for (auto it = it_b_1; it != it_b_2; ++it) {
                            const uint32_t cut_b_pos = *it;
                            const uint32_t offset = rev_b ? (end_b - cut_b_pos) : (cut_b_pos - beg_b);
                            const uint32_t mapped_a_pos = rev_a ? (end_a - offset) : (beg_a + offset);

                            if (mapped_a_pos > 0 && mapped_a_pos < len_a) {
                                new_a_cuts.push_back(mapped_a_pos);
                            }
                        }

                        auto it_a_1 = std::lower_bound(a_cuts.begin(), a_cuts.end(), beg_a);
                        auto it_a_2 = std::upper_bound(a_cuts.begin(), a_cuts.end(), end_a);

                        for (auto it = it_a_1; it != it_a_2; ++it) {
                            const uint32_t cut_a_pos = *it;
                            const uint32_t offset = rev_a ? (end_a - cut_a_pos) : (cut_a_pos - beg_a);
                            const uint32_t mapped_b_pos = rev_b ? (end_b - offset) : (beg_b + offset);

                            if (mapped_b_pos > 0 && mapped_b_pos < len_b) {
                                new_b_cuts.push_back(mapped_b_pos);
                            }
                        }
                    } break;

                    case 'X':
                        step(pos_a, rev_a, op.len);
                        step(pos_b, rev_b, op.len);
                        break;

                    case 'I':
                        step(pos_b, rev_b, op.len);
                        break;

                    case 'D':
                        step(pos_a, rev_a, op.len);
                        break;

                    case 'S':
                    case 'H':
                        break;

                    default:
                        warning_stream() << "  ! Unknown CIGAR op: " << op.op << "\n";
                        break;
                }

                if (pos_a > len_a || pos_b > len_b) {
                    warning_stream()
                        << "  ! CIGAR stepping overflow: pos_a=" << pos_a
                        << " len_a=" << len_a << " pos_b=" << pos_b
                        << " len_b=" << len_b << "\n";
                    break;
                }
            }

            const bool increase_a = merge_and_dedup_cuts_(a_cuts, new_a_cuts);
            const bool increase_b = merge_and_dedup_cuts_(b_cuts, new_b_cuts);

            if (increase_a && dirty_seen.insert(seg_a).second) {
                next_dirty.push_back(seg_a);
            }
            if (increase_b && dirty_seen.insert(seg_b).second) {
                next_dirty.push_back(seg_b);
            }

            if (increase_a || increase_b) {
                ++total_changed;
            }
        }

        dirty.swap(next_dirty);
    }

    return total_changed;
}

uint64_t GfaDeoverlapper::prune_cuts_run_(const std::vector<size_t>& group) {
    if (MAX_ABNORMAL_CUT_LEN_ <= 0 || MIN_ABNORMAL_CUT_COUNT_ <= 0) {
        return 0;
    }

    std::vector<uint32_t> segs;
    segs.reserve(group.size() * 2);

    for (const size_t aln_id : group) {
        const auto& aln = bubble_aligns_[aln_id];

        const uint32_t seg_a = NodeHandle::get_segment_id(aln.v_a);
        const uint32_t seg_b = NodeHandle::get_segment_id(aln.v_b);

        if (seg_a < cuts_.size()) segs.push_back(seg_a);
        if (seg_b < cuts_.size()) segs.push_back(seg_b);
    }

    std::sort(segs.begin(), segs.end());
    segs.erase(std::unique(segs.begin(), segs.end()), segs.end());

    uint64_t n_removed = 0;

    for (uint32_t sid : segs) {
        auto& cuts = cuts_[sid].v;
        if (cuts.size() <= 2) continue;

        std::vector<uint32_t> filtered;
        filtered.reserve(cuts.size());
        filtered.push_back(cuts.front());

        size_t i = 0;
        while (i + 1 < cuts.size()) {
            const uint32_t seg_len = cuts[i + 1] - cuts[i];

            if (seg_len <= MAX_ABNORMAL_CUT_LEN_) {
                const size_t abnormal_beg = i;

                while (i + 1 < cuts.size() && cuts[i + 1] - cuts[i] <= MAX_ABNORMAL_CUT_LEN_) {
                    ++i;
                }

                const size_t abnormal_seg_count = i - abnormal_beg;

                if (abnormal_seg_count >= MIN_ABNORMAL_CUT_COUNT_) {
                    n_removed += abnormal_seg_count - 1;

                    if (filtered.back() != cuts[i]) {
                        filtered.push_back(cuts[i]);
                    }
                } else {
                    for (size_t j = abnormal_beg + 1; j <= i; ++j) {
                        filtered.push_back(cuts[j]);
                    }
                }
            } else {
                filtered.push_back(cuts[i + 1]);
                ++i;
            }
        }

        if (filtered.size() != cuts.size()) {
            cuts.swap(filtered);
        }
    }

    return n_removed;
}


std::pair<uint64_t, uint64_t> GfaDeoverlapper::build_propagate_prune_cuts_run_(const std::vector<size_t>& group) {
    uint64_t init_changed = 0;

    for (size_t aln_id : group) {
        if (add_cuts_from_one_alignment_(bubble_aligns_[aln_id])) {
            ++init_changed;
        }
    }

    normalize_group_cuts_(group);
    const uint64_t propagated_changed = propagate_cuts_run_(group);
    normalize_group_cuts_(group);
    const uint64_t n_removed = prune_cuts_run_(group);
    normalize_group_cuts_(group);

    return {
        init_changed + propagated_changed,
        n_removed
    };
}

void GfaDeoverlapper::build_propagate_prune_cuts_(const std::vector<std::vector<size_t>>& groups) {
    log_stream() << "Building, propagating, and pruning cut points ...\n";

    if (bubble_aligns_.empty()) {
        for (uint32_t s = 0; s < cuts_.size(); ++s) {
            auto& C = cuts_[s].v;
            if (C.empty() || C.front() != 0u) C.insert(C.begin(), 0u);
            if (C.back() != nodes_[s].length) C.push_back(nodes_[s].length);
        }

        log_stream() << "  - Cut point processing completed.\n\n";
        return;
    }

    ThreadPool pool(alignOpts_.threads);
    std::vector<std::future<std::pair<uint64_t, uint64_t>>> futs;
    futs.reserve(groups.size());

    for (const auto& group : groups) {
        const auto* group_ptr = &group;
        futs.emplace_back(pool.submit([this, group_ptr]() {
            return build_propagate_prune_cuts_run_(*group_ptr);
        }));
    }

    ProgressTracker prog(groups.size());

    uint64_t total_changed = 0;
    uint64_t total_removed = 0;

    for (auto& f : futs) {
        auto [changed, removed] = f.get();
        prog.hit();

        total_changed += changed;
        total_removed += removed;
    }

    prog.finish();
    pool.stop();

    normalize_all_cuts_();

    log_stream() << "  - Alignments with cut changes: " << total_changed << "\n";
    log_stream() << "  - Removed abnormal cut points: " << total_removed << "\n\n";
}


void GfaDeoverlapper::print_cuts(const std::vector<Cuts>& cuts) {
    log_stream() << "Segment cuts:\n";
    for (size_t i = 0; i < cuts.size(); ++i) {
        const auto& cut = cuts[i];
        std::string log = "  - " + cut.name + ": ";
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
    const std::vector<uint32_t>& cut_in_overlap_leaf, 
    PairRuleIndex_& index
) const {
    /**
     * @brief Merge new_value into old_value if they have the same granularity (i.e., same number of segments and same segment lengths). Return true if old_value is updated, false otherwise.
     * @date 2026-05-18
     * @version 0.1.3-r2
     */
    auto merge_same_granularity_value = [&](std::vector<SegReplace::Seg>& old_value, const std::vector<SegReplace::Seg>& new_value) -> bool
    {
        if (old_value.size() != new_value.size()) {
            return new_value.size() > old_value.size();
        }

        for (size_t i = 0; i < old_value.size(); ++i) {
            const uint32_t old_beg = SegReplace::Interval::beg(old_value[i]);
            const uint32_t old_end = SegReplace::Interval::end(old_value[i]);
            const uint32_t new_beg = SegReplace::Interval::beg(new_value[i]);
            const uint32_t new_end = SegReplace::Interval::end(new_value[i]);

            const uint32_t old_len = std::max(old_beg, old_end) - std::min(old_beg, old_end);
            const uint32_t new_len = std::max(new_beg, new_end) - std::min(new_beg, new_end);

            if (old_len != new_len) {
                return false;
            }
        }

        for (size_t i = 0; i < old_value.size(); ++i) {
            const uint64_t old_sid = SegReplace::Interval::seg_id(old_value[i]);
            const uint64_t new_sid = SegReplace::Interval::seg_id(new_value[i]);

            if (old_sid == new_sid) continue;

            if (old_sid < new_sid) {
                rulemap.emplace(new_value[i], std::vector<SegReplace::Seg>{old_value[i]});
            } else {
                rulemap.emplace(old_value[i], std::vector<SegReplace::Seg>{new_value[i]});
                old_value[i] = new_value[i];
            }
        }

        return false;
    };


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

            auto it = rulemap.find(trunk_u128);
            if (it == rulemap.end()) {
                add_pair_rule_index_(index, trunk_u128, leaf_parts);
                rulemap.emplace(trunk_u128, std::move(leaf_parts));
            } else if (merge_same_granularity_value(it->second, leaf_parts)) {
                add_pair_rule_index_(index, trunk_u128, leaf_parts);
                it->second = std::move(leaf_parts);
            }
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
            
            auto it = rulemap.find(leaf_win_u128);
            if (it == rulemap.end()) {
                add_pair_rule_index_(index, leaf_win_u128, leaf_parts);
                rulemap.emplace(leaf_win_u128, std::move(leaf_parts));
            } else if (merge_same_granularity_value(it->second, leaf_parts)) {
                add_pair_rule_index_(index, leaf_win_u128, leaf_parts);
                it->second = std::move(leaf_parts);
            }
        }
    }
}


GfaDeoverlapper::PairRuleIndex_ GfaDeoverlapper::build_pair_rule_index_(
    const SegReplace::RuleMap& rulemap
) const {
    PairRuleIndex_ index;
    index.reserve(rulemap.size() * 2 + 1);

    for (const auto& kv : rulemap) {
        add_pair_rule_index_(index, kv.first, kv.second);
    }

    return index;
}


void GfaDeoverlapper::add_pair_rule_index_(
    PairRuleIndex_& index,
    SegReplace::Seg key,
    const std::vector<SegReplace::Seg>& value
) const {
    if (value.empty()) return;

    auto make_side = [&](const std::vector<SegReplace::Seg>& segs, uint32_t& sid, uint32_t& beg, uint32_t& end, bool& rev) -> bool
    {
        if (segs.empty()) return false;

        const uint64_t raw_sid = SegReplace::Interval::seg_id(segs.front());
        if (raw_sid > UINT32_MAX || raw_sid >= nodes_.size()) return false;

        sid = static_cast<uint32_t>(raw_sid);
        rev = SegReplace::Interval::is_reverse(segs.front());

        uint32_t lo = UINT32_MAX;
        uint32_t hi = 0;

        for (const auto& seg : segs) {
            if (SegReplace::Interval::seg_id(seg) != raw_sid) return false;
            if (SegReplace::Interval::is_reverse(seg) != rev) return false;

            const uint32_t b = SegReplace::Interval::beg(seg);
            const uint32_t e = SegReplace::Interval::end(seg);

            lo = std::min(lo, std::min(b, e));
            hi = std::max(hi, std::max(b, e));
        }

        if (lo >= hi || hi > nodes_[sid].length) return false;

        beg = lo;
        end = hi;
        return true;
    };

    PairRuleWin_ w;

    {
        const uint64_t raw_sid = SegReplace::Interval::seg_id(key);
        if (raw_sid > UINT32_MAX || raw_sid >= nodes_.size()) return;

        const uint32_t b = SegReplace::Interval::beg(key);
        const uint32_t e = SegReplace::Interval::end(key);

        w.a_sid = static_cast<uint32_t>(raw_sid);
        w.a_beg = std::min(b, e);
        w.a_end = std::max(b, e);
        w.a_rev = SegReplace::Interval::is_reverse(key);

        if (w.a_beg >= w.a_end || w.a_end > nodes_[w.a_sid].length) return;
    }

    if (!make_side(value, w.b_sid, w.b_beg, w.b_end, w.b_rev)) return;
    if (w.a_sid == w.b_sid) return;

    if (w.a_sid > w.b_sid) {
        std::swap(w.a_sid, w.b_sid);
        std::swap(w.a_beg, w.b_beg);
        std::swap(w.a_end, w.b_end);
        std::swap(w.a_rev, w.b_rev);
    }

    const uint32_t len_a = nodes_[w.a_sid].length;
    const uint32_t len_b = nodes_[w.b_sid].length;

    const uint32_t a_obeg = w.a_rev ? (len_a - w.a_end) : w.a_beg;
    const uint32_t b_obeg = w.b_rev ? (len_b - w.b_end) : w.b_beg;
    w.delta = int64_t(b_obeg) - int64_t(a_obeg);

    const uint64_t pair_key = (uint64_t(w.a_sid) << 32) | uint64_t(w.b_sid);
    auto& vec = index[pair_key];

    for (auto& old : vec) {
        if (old.a_rev != w.a_rev || old.b_rev != w.b_rev || old.delta != w.delta) continue;

        const bool touch_a = w.a_beg <= old.a_end && old.a_beg <= w.a_end;
        const bool touch_b = w.b_beg <= old.b_end && old.b_beg <= w.b_end;
        if (!touch_a || !touch_b) continue;

        old.a_beg = std::min(old.a_beg, w.a_beg);
        old.a_end = std::max(old.a_end, w.a_end);
        old.b_beg = std::min(old.b_beg, w.b_beg);
        old.b_end = std::max(old.b_end, w.b_end);
        return;
    }

    vec.push_back(w);
}


bool GfaDeoverlapper::has_conflicting_pair_rule_(
    PairRuleIndex_& index,
    const std::vector<SegReplace::Seg>& exist_segs,
    const std::vector<SegReplace::Seg>& leaf_segs
) const {
    if (exist_segs.empty() || leaf_segs.empty()) return false;

    PairRuleIndex_ tmp;
    add_pair_rule_index_(tmp, exist_segs.front(), leaf_segs);
    if (tmp.empty()) return false;

    const PairRuleWin_& cur = tmp.begin()->second.front();
    const uint64_t pair_key = (uint64_t(cur.a_sid) << 32) | uint64_t(cur.b_sid);

    auto it = index.find(pair_key);
    if (it == index.end()) return false;

    auto overlap = [](uint32_t a0, uint32_t a1, uint32_t b0, uint32_t b1) -> bool {
        return a0 < a1 && b0 < b1 && std::max(a0, b0) < std::min(a1, b1);
    };

    for (const auto& old : it->second) {
        const bool ov_a = overlap(cur.a_beg, cur.a_end, old.a_beg, old.a_end);
        const bool ov_b = overlap(cur.b_beg, cur.b_end, old.b_beg, old.b_end);

        if (!ov_a && !ov_b) continue;

        const bool same =
            cur.a_rev == old.a_rev &&
            cur.b_rev == old.b_rev &&
            cur.delta == old.delta;

        if (!same) return true;
    }

    return false;
}


void GfaDeoverlapper::record_rulemap(
    uint32_t v_seg_id, uint32_t w_seg_id,
    bool v_is_rev, bool w_is_rev,
    const std::string& v_name, const std::string& w_name,
    uint32_t vb, uint32_t ve, uint32_t wb, uint32_t we,
    SegReplace::RuleMap& rulemap, 
    std::vector<std::unordered_set<uint32_t>>& cut_sets, 
    PairRuleIndex_& index
) {
    // ------------------------------------------------ Helper functions ------------------------------------------------
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

    auto& v_cut_set = cut_sets[v_seg_id];
    auto& w_cut_set = cut_sets[w_seg_id];

    auto it_v_1 = std::lower_bound(v_cuts.begin(), v_cuts.end(), vb);
    auto it_v_2 = std::upper_bound(v_cuts.begin(), v_cuts.end(), ve);
    auto it_w_1 = std::lower_bound(w_cuts.begin(), w_cuts.end(), wb);
    auto it_w_2 = std::upper_bound(w_cuts.begin(), w_cuts.end(), we);
    std::vector<uint32_t> cut_in_overlap_v(it_v_1, it_v_2);
    std::vector<uint32_t> cut_in_overlap_w(it_w_1, it_w_2);

    bool v_is_leaf;
    if (!prefer_v_as_leaf_(v_seg_id, w_seg_id, v_is_leaf)) return;

    auto&        leaf_cuts     = v_is_leaf ? v_cuts : w_cuts;
    auto&        leaf_cut_set  = v_is_leaf ? v_cut_set : w_cut_set;
    // auto&        trunk_cuts    = v_is_leaf ? w_cuts : v_cuts;
    // auto&        trunk_cut_set = v_is_leaf ? w_cut_set : v_cut_set;
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

    // leaf name roots
    bool leaf_base_names_ready = false;
    std::vector<std::string> leaf_base_names_cache;

    // reverse according to direction
    if (leaf_is_rev)  std::reverse(leaf_vec.begin(),  leaf_vec.end());
    if (trunk_is_rev) std::reverse(trunk_vec.begin(), trunk_vec.end());

    // Compute leaf offsets
    std::vector<int32_t> trunk_offsets = build_offsets(trunk_vec, tb, te, trunk_is_rev);
    std::vector<int32_t> leaf_offsets  = build_offsets(leaf_vec, lb, le, leaf_is_rev);

    if (trunk_offsets.size() < 2 || leaf_offsets.empty()) return;

    for (size_t l_idx = 0; l_idx + 1 < trunk_offsets.size(); ++l_idx) {
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

        if (leaf_end_idx == UINT32_MAX) {  // No matching cut found in leaf_vec, skip this trunk window
            continue;
        }

        if (leaf_end_idx <= leaf_beg_idx) {
            warning_stream() << "Unexpected index detected!\n";
            warning_stream() << "  - trunk_offset=" << trunk_beg_off << " - " << trunk_end_off << "\n";
            warning_stream() << "  - leaf_idx=" << leaf_beg_idx << " - " << leaf_end_idx << "\n";
            std::string log_tmp = "  - leaf_offsets=(";
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
                    << "  - Trunk segment overlap window: [" << trunk_vec.front() << " .. " << trunk_vec.back() << "]\n"
                    << "  - Leaf  segment overlap window: [" << leaf_vec.front()  << " .. " << leaf_vec.back()  << "]\n"
                    << "  - Current trunk position: beg=" << trunk_beg << " end=" << trunk_end << "\n"
                    << "  - Current trunk offsets: beg=" << trunk_beg_off << " end=" << trunk_end_off << "\n"
                    << "  - Mapped leaf idx range: " << leaf_beg_idx << " .. " << leaf_end_idx << "\n";

            std::cerr << "  - trunk_vec (" << trunk_name << "): ";
            for (auto p : trunk_vec) std::cerr << p << " ";
            std::cerr << "\n";

            std::cerr << "  - leaf_offsets: ";
            for (auto off : leaf_offsets) std::cerr << off << " ";
            std::cerr << "\n";
        }

        std::vector<SegReplace::Seg> leaf_segs;
        leaf_segs.reserve(leaf_end_idx - leaf_beg_idx);
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
            add_pair_rule_index_(index, trunk_u128, leaf_segs);
        } else {
            const auto& exist_segs = it_rule->second;

            if (exist_segs.size() == 0 || leaf_segs.size() == 0) {
                error_stream() << "Empty segs detected!\n";
                error_stream() << "  - Edge: " << v_name << "(" << v_seg_id << ",rev=" << v_is_rev << ") <-> " << w_name << "(" << w_seg_id << ",rev=" << w_is_rev << ")\n";

                std::string log_tmp;

                log_tmp = "  - exist segs: ";
                if (exist_segs.empty()) log_tmp += "(EMPTY)";
                else for (auto& seg : exist_segs) log_tmp += SegReplace::Interval::format(seg) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                log_tmp = "  - leaf_segs: ";
                if (leaf_segs.empty()) log_tmp += "(EMPTY)";
                else for (auto& seg : leaf_segs) log_tmp += SegReplace::Interval::format(seg) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                log_tmp = "  - trunk_vec (" + trunk_name + "): ";
                for (auto p : trunk_vec) log_tmp += std::to_string(p) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                error_stream() << "  - Trunk range: beg=" << trunk_beg << " end=" << trunk_end << " offsets=(" << trunk_beg_off << "," << trunk_end_off << ")\n";

                log_tmp = "  - leaf_vec (" + leaf_name + "): ";
                for (auto p : leaf_vec) log_tmp += std::to_string(p) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                log_tmp = "  - leaf_offsets: ";
                for (auto off : leaf_offsets) log_tmp += std::to_string(off) + " ";
                log_tmp += "\n";
                error_stream() << log_tmp;

                error_stream() << "  - Leaf index range: " << leaf_beg_idx << " .. " << leaf_end_idx
                                << " (index_len=" << (leaf_end_idx >= leaf_beg_idx ? (leaf_end_idx - leaf_beg_idx) : -1) << ")\n";

                if (!leaf_vec.empty() && leaf_beg_idx < leaf_vec.size() && leaf_end_idx < leaf_vec.size()) {
                    error_stream() << "  - Leaf range: beg=" << leaf_vec[leaf_beg_idx] << " end=" << leaf_vec[leaf_end_idx] << "\n";
                }
                exit(1);
            }

            /**
             * @brief Check if there is a conflicting pair rule in the existing rulemap for the current trunk and leaf segments. A conflicting pair rule is defined as follows:
             *   - If exist_segs and leaf_segs already have an old relation in rulemap, but that relation is different from the current candidate mapping, skip the current mapping to avoid introducing conflicting replacement cycles.
             * @date 2026-05-18
             * @version 0.1.3-r2
             * @note This check is necessary to prevent cycles from being introduced into the graph.
             */
            if (has_conflicting_pair_rule_(index, exist_segs, leaf_segs)) {
                if (DEBUG_ENABLED) {
                    std::cerr << "  - Skip conflicting propagated rule:\n";
                    std::cerr << "    - trunk: " << SegReplace::Interval::format(trunk_u128) << "\n";

                    std::cerr << "    - exist segs: ";
                    for (const auto& seg : exist_segs) {
                        std::cerr << SegReplace::Interval::format(seg) << " ";
                    }
                    std::cerr << "\n";

                    std::cerr << "    - leaf segs : ";
                    for (const auto& seg : leaf_segs) {
                        std::cerr << SegReplace::Interval::format(seg) << " ";
                    }
                    std::cerr << "\n";
                }
                continue;
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
            if (!leaf_base_names_ready) {
                gfaName::collect_roots_from_name(leaf_name, leaf_base_names_cache);
                leaf_base_names_ready = true;
            }
            std::vector<std::string> exist_base_names;
            exist_base_names.reserve(exist_segs.size() * 2);

            for (const auto& seg : exist_segs) {
                const std::string& exist_name = nodes_[SegReplace::Interval::seg_id(seg)].name;
                std::vector<std::string> exist_names_tmp;
                gfaName::collect_roots_from_name(exist_name, exist_names_tmp);
                exist_base_names.insert(exist_base_names.end(), exist_names_tmp.begin(), exist_names_tmp.end());
            }
            if (std::any_of(leaf_base_names_cache.begin(), leaf_base_names_cache.end(), [&](const std::string& s) {return std::find(exist_base_names.begin(), exist_base_names.end(), s) != exist_base_names.end();})) {
                continue;
            }

            auto& exist_cuts = cuts_[exist_seg_id].v;
            auto& exist_cut_set = cut_sets[exist_seg_id];

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

            std::vector<uint32_t> extra_leaf_cuts;
            std::vector<uint32_t> extra_exist_cuts;
            extra_leaf_cuts.reserve(cut_in_overlap_exist.size());
            extra_exist_cuts.reserve(cut_in_overlap_leaf.size());

            // exist -> leaf
            for (uint32_t cut : cut_in_overlap_exist) {
                uint64_t off = offset_from(cut, exist_beg, exist_end, exist_is_rev);
                uint32_t mapped = offset_to(off, leaf_beg, leaf_end, leaf_is_rev1);
                if (mapped < leaf_beg || mapped > leaf_end) continue;

                if (leaf_cut_set.find(mapped) == leaf_cut_set.end()) {
                    extra_leaf_cuts.push_back(mapped);
                }
            }

            if (!extra_leaf_cuts.empty()) {
                std::sort(extra_leaf_cuts.begin(), extra_leaf_cuts.end());
                extra_leaf_cuts.erase(std::unique(extra_leaf_cuts.begin(), extra_leaf_cuts.end()), extra_leaf_cuts.end());

                new_cuts_leaf.reserve(cut_in_overlap_leaf.size() + extra_leaf_cuts.size());
                new_cuts_leaf.insert(new_cuts_leaf.end(), extra_leaf_cuts.begin(), extra_leaf_cuts.end());

                std::inplace_merge(
                    new_cuts_leaf.begin(),
                    new_cuts_leaf.begin() + cut_in_overlap_leaf.size(),
                    new_cuts_leaf.end()
                );
            }

            // leaf -> exist
            for (uint32_t cut : cut_in_overlap_leaf) {
                uint64_t off = offset_from(cut, leaf_beg, leaf_end, leaf_is_rev1);
                uint32_t mapped = offset_to(off, exist_beg, exist_end, exist_is_rev);
                if (mapped < exist_beg || mapped > exist_end) continue;

                if (exist_cut_set.find(mapped) == exist_cut_set.end()) {
                    extra_exist_cuts.push_back(mapped);
                }
            }

            if (!extra_exist_cuts.empty()) {
                std::sort(extra_exist_cuts.begin(), extra_exist_cuts.end());
                extra_exist_cuts.erase(std::unique(extra_exist_cuts.begin(), extra_exist_cuts.end()), extra_exist_cuts.end());

                new_cuts_exist.reserve(cut_in_overlap_exist.size() + extra_exist_cuts.size());
                new_cuts_exist.insert(new_cuts_exist.end(), extra_exist_cuts.begin(), extra_exist_cuts.end());

                std::inplace_merge(
                    new_cuts_exist.begin(),
                    new_cuts_exist.begin() + cut_in_overlap_exist.size(),
                    new_cuts_exist.end()
                );
            }

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
                std::cerr << "  - Trunk: " << trunk_name << " (seg=" << exist_seg_id
                          << ", rev=" << exist_is_rev << ") window=[" << exist_beg << "," << exist_end
                          << "] #cuts=" << exist_diff << "\n   - cuts: ";
                for (auto c : cut_in_overlap_exist) std::cerr << c << " ";
                std::cerr << "\n";

                std::cerr << "  - Leaf : " << leaf_name  << " (seg=" << leaf_seg_id
                          << ", rev=" << leaf_is_rev1 << ") window=[" << leaf_beg << "," << leaf_end
                          << "] #cuts=" << leaf_diff << "\n   - cuts: ";
                for (auto c : cut_in_overlap_leaf)  std::cerr << c << " ";
                std::cerr << "\n";
                continue;
            }

            if (DEBUG_ENABLED && (exist_diff > 2 || leaf_diff > 2)) {
                std::string log_tmp = "  - exist segs: ";
                for (auto& seg : exist_segs) log_tmp += SegReplace::Interval::format(seg) + " ";
                log_tmp += "\n";
                log_tmp += "  - leaf_segs: ";
                for (auto& seg : leaf_segs) log_tmp += SegReplace::Interval::format(seg) + " ";
                log_tmp += "\n";
                std::cerr << log_tmp;

                std::cerr << "  - Existing segment cuts(" << exist_seg_id << "-" << exist_diff << "): ";
                for (auto cut : cut_in_overlap_exist) std::cerr << cut << " ";
                std::cerr << "\n";
                std::cerr << "  - Leaf segment cuts(" << leaf_seg_id << "-" << leaf_diff << "): ";
                for (auto cut : cut_in_overlap_leaf) std::cerr << cut << " ";
                std::cerr << "\n";

                std::cerr << "     - exist original: ";
                for (const auto& cut : cut_in_overlap_exist) std::cerr << cut << " ";
                std::cerr << "\n";
                std::cerr << "     -      exist new: ";
                for (const auto& cut : new_cuts_exist) std::cerr << cut << " ";
                std::cerr << "\n";
                std::cerr << "     -  leaf original: ";
                for (const auto& cut : cut_in_overlap_leaf) std::cerr << cut << " ";
                std::cerr << "\n";
                std::cerr << "     -       leaf new: ";
                for (const auto& cut : new_cuts_leaf) std::cerr << cut << " ";
                std::cerr << "\n";
            }

            // SegReplace::Seg exist_u128 = SegReplace::Interval::pack(exist_seg_id, exist_beg, exist_end, exist_is_rev);
            // SegReplace::Seg leaf_u128  = SegReplace::Interval::pack(leaf_seg_id,  leaf_beg,  leaf_end,  leaf_is_rev1);

            bool leaf_is_exist;
            if (!prefer_v_as_leaf_(exist_seg_id, leaf_seg_id, leaf_is_exist)) continue;

            if (leaf_is_exist) {
                // merge_and_dedup_cuts_(exist_cuts, new_cuts_exist);  // debug 2025-10-10, Update cut list
                merge_and_dedup_cuts_with_set_(exist_cuts, exist_cut_set, extra_exist_cuts);  // 2025.04.10 update both cut list and cut set

                build_rules_from_pair_windows_(
                    rulemap,
                    /*trunk*/ leaf_seg_id, leaf_is_rev1, cut_in_overlap_leaf,
                    /*leaf*/ exist_seg_id, exist_is_rev, new_cuts_exist, new_exist_offsets,
                    /*leaf older*/ cut_in_overlap_exist, index
                );
            } else {
                // merge_and_dedup_cuts_(leaf_cuts, new_cuts_leaf);  // debug 2025-10-10, Update cut list
                merge_and_dedup_cuts_with_set_(leaf_cuts, leaf_cut_set, extra_leaf_cuts);  // 2025.04.10 update both cut list and cut set

                build_rules_from_pair_windows_(
                    rulemap,
                    /*trunk*/ exist_seg_id, exist_is_rev, cut_in_overlap_exist,
                    /*leaf*/ leaf_seg_id, leaf_is_rev1, new_cuts_leaf, new_leaf_offsets,
                    /*leaf older*/ cut_in_overlap_leaf, index
                );
            }
        }

        if (DEBUG_ENABLED) {
            std::cerr << "  - Rule added: " << SegReplace::Interval::format(trunk_u128) << " -> ";
            for (const auto& leaf_seg_u : leaf_segs) {
                std::cerr << SegReplace::Interval::format(leaf_seg_u) << " ";
            }
            std::cerr << "\n";
        }
    }
}


SegReplace::RuleMap GfaDeoverlapper::build_rulemap_run_(
    const std::vector<size_t>& align_ids,
    std::vector<std::unordered_set<uint32_t>>& cut_sets
) {
    SegReplace::RuleMap rulemap;
    PairRuleIndex_ pair_rule_index = build_pair_rule_index_(rulemap);  // For quick lookup of existing cuts when building rules, reduce the time complexity in sort and deduo setp.

    for (const size_t align_id : align_ids) {
        if (align_id >= bubble_aligns_.size()) continue;
        const auto& align = bubble_aligns_[align_id];

        // If the alignment result length is 1 and all are "=", directly add the cut point without filtering with MIN_EQ_FOR_CUT_
        const bool is_single_exact_match = (align.ops.size() == 1 && align.ops[0].op == '=');

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

                    if (op.len >= MIN_EQ_FOR_CUT_ || is_single_exact_match) {
                        GfaDeoverlapper::record_rulemap(
                            seg_a, seg_b,
                            rev_a, rev_b,
                            name_a, name_b,
                            beg_a, end_a, beg_b, end_b,
                            rulemap, cut_sets, pair_rule_index
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
                    warning_stream() << "  - Unknown CIGAR op: " << op.op << "\n";
                    break;
            }

            if (pos_a > len_a || pos_b > len_b) {
                warning_stream() 
                    << "  - CIGAR stepping overflow: pos_a=" << pos_a
                    << " len_a=" << len_a << " pos_b=" << pos_b
                    << " len_b=" << len_b << "\n";
                break;
            }
        }
    }

    return rulemap;
}


void GfaDeoverlapper::build_rulemap_(const std::vector<std::vector<size_t>>& groups) {
    log_stream() << "Building replacement rules from alignments ...\n";

    rulemap_.clear();
    if (bubble_aligns_.empty()) {
        log_stream() << "  - Total rules generated: 0\n\n";
        return;
    }

    std::vector<std::unordered_set<uint32_t>> cut_sets = build_cut_sets_();

    ThreadPool pool(alignOpts_.threads);
    std::vector<std::future<SegReplace::RuleMap>> futs;
    futs.reserve(groups.size());

    for (const auto& group : groups) {
        const auto* group_ptr = &group;
        futs.emplace_back(pool.submit([this, group_ptr, &cut_sets]() {
            return build_rulemap_run_(*group_ptr, cut_sets);
        }));
    }

    ProgressTracker prog(groups.size());

    for (auto& f : futs) {
        SegReplace::RuleMap local = f.get();
        prog.hit();

        for (auto& kv : local) {
            rulemap_.emplace(kv.first, std::move(kv.second));
        }
    }

    prog.finish();
    pool.stop();

    normalize_all_cuts_();

    log_stream() << "  - Total rules generated: " << rulemap_.size() << "\n\n";
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
            const auto e = ex.query(big);
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


bool GfaDeoverlapper::prefer_v_as_leaf_(uint32_t v_seg_id, uint32_t w_seg_id, bool& v_is_leaf, bool by_name) const {
    if (v_seg_id == w_seg_id) return false;

    if (!by_name) {
        v_is_leaf = v_seg_id < w_seg_id;
        return true;
    }

    // by name
    if (v_seg_id >= nodes_.size() || w_seg_id >= nodes_.size()) {
        error_stream() << "Invalid segment IDs provided: " 
                    << v_seg_id << ", " << w_seg_id 
                    << ". Segment ID exceeds current node list size (" 
                    << nodes_.size() << ").\n";
        std::exit(1);
    }
    const std::string& vn = nodes_[v_seg_id].name;
    const std::string& wn = nodes_[w_seg_id].name;

    v_is_leaf = vn < wn;

    return true;
}


bool GfaDeoverlapper::merge_and_dedup_cuts_(
    std::vector<uint32_t>& cuts,
    const std::vector<uint32_t>& new_cuts
) {
    if (new_cuts.empty()) {
        return false;
    }

    std::vector<uint32_t> normalized = new_cuts;

    if (normalized.size() >= 2) {
        bool asc = true;
        bool desc = true;

        for (size_t i = 1; i < normalized.size(); ++i) {
            if (normalized[i] < normalized[i - 1]) asc = false;
            if (normalized[i] > normalized[i - 1]) desc = false;
            if (!asc && !desc) break;
        }

        if (desc && !asc) {
            std::reverse(normalized.begin(), normalized.end());
        } else if (!asc && !desc) {
            std::sort(normalized.begin(), normalized.end());
        }
    }

    normalized.erase(
        std::unique(normalized.begin(), normalized.end()),
        normalized.end()
    );

    if (normalized.empty()) {
        return false;
    }

    {
        size_t i = 0;
        size_t j = 0;
        bool has_new_element = false;

        while (i < cuts.size() && j < normalized.size()) {
            if (cuts[i] < normalized[j]) {
                ++i;
            } else if (cuts[i] > normalized[j]) {
                has_new_element = true;
                break;
            } else {
                ++i;
                ++j;
            }
        }

        if (!has_new_element && j < normalized.size()) {
            has_new_element = true;
        }

        if (!has_new_element) {
            return false;
        }
    }

    std::vector<uint32_t> merged;
    merged.reserve(cuts.size() + normalized.size());

    size_t i = 0;
    size_t j = 0;

    while (i < cuts.size() && j < normalized.size()) {
        if (cuts[i] < normalized[j]) {
            merged.push_back(cuts[i]);
            ++i;
        } else if (cuts[i] > normalized[j]) {
            merged.push_back(normalized[j]);
            ++j;
        } else {
            merged.push_back(cuts[i]);
            ++i;
            ++j;
        }
    }

    while (i < cuts.size()) {
        merged.push_back(cuts[i]);
        ++i;
    }

    while (j < normalized.size()) {
        merged.push_back(normalized[j]);
        ++j;
    }

    cuts.swap(merged);
    return true;
}


std::vector<std::unordered_set<uint32_t>> GfaDeoverlapper::build_cut_sets_() const {
    std::vector<std::unordered_set<uint32_t>> cut_sets(cuts_.size());

    for (size_t seg_id = 0; seg_id < cuts_.size(); ++seg_id) {
        const auto& cuts = cuts_[seg_id].v;
        auto& cut_set = cut_sets[seg_id];
        cut_set.reserve(cuts.size() * 2 + 1);
        for (uint32_t c : cuts) {
            cut_set.insert(c);
        }
    }
    return cut_sets;
}


bool GfaDeoverlapper::merge_and_dedup_cuts_with_set_(
    std::vector<uint32_t>& cuts,
    std::unordered_set<uint32_t>& cut_set,
    const std::vector<uint32_t>& candidate_cuts
) {
    if (candidate_cuts.empty()) return false;

    std::vector<uint32_t> truly_new;
    truly_new.reserve(candidate_cuts.size());

    for (uint32_t c : candidate_cuts) {
        if (cut_set.insert(c).second) {
            truly_new.push_back(c);
        }
    }

    if (truly_new.empty()) return false;

    return merge_and_dedup_cuts_(cuts, truly_new);
}


void GfaDeoverlapper::emit_node_expansion_edges_(
    const SegReplace::Expansion& node_expansion,
    std::unordered_set<uint64_t>& seen_edges,
    uint64_t& new_node_count,
    uint64_t& new_edge_count,
    gfaName& namer,
    const std::string& v_name,
    const std::string& w_name,
    bool v_is_rev,
    bool w_is_rev,
    std::optional<std::pair<uint32_t,uint32_t>> v_span,
    std::optional<std::pair<uint32_t,uint32_t>> w_span
) {
    if (DEBUG_ENABLED) {
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
    }

    auto pretty_name = [&](SegReplace::Seg s) -> std::string {
        const std::string& base = getNodeName(SegReplace::Interval::seg_id(s));
        uint32_t beg = SegReplace::Interval::beg(s);
        uint32_t end = SegReplace::Interval::end(s);
        return namer.format_interval_name(base, beg, end, false);
    };

    auto ensure_segment = [&](SegReplace::Seg s, const std::string& seg_name) -> uint32_t {
        bool is_new = false;
        uint32_t sid = get_or_add_segment(seg_name, &is_new);

        if (is_new) {
            std::string seg_seq = slice_seq_or_star_(
                nodes_,
                SegReplace::Interval::seg_id(s),
                SegReplace::Interval::beg(s),
                SegReplace::Interval::end(s),
                false
            );

            GfaNode& n = nodes_[sid];
            n.name = seg_name;
            n.sequence = seg_seq;
            n.length = seg_seq != "*" ? (uint32_t)seg_seq.length() : 0;
            total_segment_length_ += seg_seq.length();
            ++new_node_count;
        }

        return sid;
    };

    if (node_expansion.size() < 2) {
        if (DEBUG_ENABLED) { debug_stream() << "\n"; }
        return;
    }

    for (size_t i = 0; i + 1 < node_expansion.size(); ++i) {
        const SegReplace::Seg vs = node_expansion[i];
        const SegReplace::Seg ws = node_expansion[i + 1];

        std::string v_name_tmp = pretty_name(vs);
        std::string w_name_tmp = pretty_name(ws);

        bool v_rev = SegReplace::Interval::is_reverse(vs);
        bool w_rev = SegReplace::Interval::is_reverse(ws);

        uint32_t v_sid = ensure_segment(vs, v_name_tmp);
        uint32_t w_sid = ensure_segment(ws, w_name_tmp);

        uint32_t v_vid = (v_sid << 1) | (v_rev ? 1 : 0);
        uint32_t w_vid = (w_sid << 1) | (w_rev ? 1 : 0);

        if (seen_edges.insert(encode_edge_u64_(v_vid, w_vid)).second) {
            ++new_edge_count;
        }

        if (DEBUG_ENABLED) {
            debug_stream() << "  - " << v_name_tmp << (v_rev ? "-" : "+") << " -> " << w_name_tmp << (w_rev ? "-" : "+") << "\n";
        }
    }

    if (DEBUG_ENABLED) { debug_stream() << "\n"; }
}


void GfaDeoverlapper::expand_and_rewire_edges_(
    const SegReplace::Expander& ex
) {
    log_stream() << "Expanding segments and rewiring edges ...\n";

    // ------------------------------------------------ Initialization ------------------------------------------------
    // make sure each edge is only added once
    std::unordered_set<uint64_t> seen_edges;
    seen_edges.reserve(arcs_.size() * 2);
    seen_edges.max_load_factor(0.7f);

    uint64_t new_node_count = 0;
    uint64_t new_edge_count = 0;

    gfaName namer;  // for generating new segment names

    /**
     * @brief Record the node that don't have any changes.
     * @date 2026-05-19
     * @version 0.1.3-r3
     */
    std::vector<uint32_t> active_degree(nodes_.size(), 0);
    for (const auto& e : arcs_) {
        if (e.get_del() || e.get_comp()) continue;
        ++active_degree[e.get_source_segment_id()];
        ++active_degree[e.get_target_segment_id()];
    }

    // ------------------------------------------------ Helper functions ------------------------------------------------
    auto trivial_segment = [&](uint32_t sid) -> bool {
        if (sid >= nodes_.size() || nodes_[sid].deleted) return false;
        const uint32_t len = nodes_[sid].length;
        if (sid >= cuts_.size()) return true;
        if (!(cuts_[sid].v.size() == 2 && cuts_[sid].v[0] == 0 && cuts_[sid].v[1] == len)) return false;
        if (len == 0) return true;

        auto f = SegReplace::Interval::pack(sid, 0, len, false);
        auto r = SegReplace::Interval::pack(sid, 0, len, true);
        return is_identity_(ex.query(f), f) && is_identity_(ex.query(r), r);
    };

    // original = original window_u128 (2026-05-22 v0.1.3-r4)
    // query    = replacement expansion from ex.query(window_u128)
    // If query produces a repeated interval later, fallback to the affected piece range to original windows.
    struct ExpansionPiece_ {
        SegReplace::Seg original;
        SegReplace::Expansion query;
    };

    // Rule (2026-05-22 v0.1.3-r4):
    //   If the same SegReplace::Seg appears multiple times in xs[*].query, then fallback all pieces from the first appearance to the current piece.
    //
    // Example:
    //   piece[2].query contains A
    //   piece[5].query contains A again
    //
    // Then:
    //   piece[2], piece[3], piece[4], piece[5]
    // are all replaced by their original windows.
    auto sanitize_expansions = [&](std::vector<ExpansionPiece_>& xs) -> std::vector<SegReplace::Expansion> {
        // ************************************************** Have bug here, need to rewrite later (Start) (2026-05-22 v0.1.3-r4) **************************************************
        // std::unordered_map<SegReplace::Seg, size_t, SegReplace::U128Hash, SegReplace::U128Eq> first_seen;
        // first_seen.reserve(xs.size() * 2 + 1);

        // std::vector<uint8_t> use_original(xs.size(), 0);

        // for (size_t i = 0; i < xs.size(); ++i) {
        //     for (const auto& seg : xs[i].query) {
        //         auto it = first_seen.find(seg);

        //         if (it == first_seen.end()) {
        //             first_seen.emplace(seg, i);
        //             continue;
        //         }

        //         for (size_t j = it->second; j <= i; ++j) {
        //             use_original[j] = 1;
        //         }
        //     }
        // }

        // std::vector<SegReplace::Expansion> out;
        // out.reserve(xs.size());

        // for (size_t i = 0; i < xs.size(); ++i) {
        //     if (use_original[i]) {
        //         out.push_back(SegReplace::Expansion{xs[i].original});
        //     } else {
        //         out.push_back(std::move(xs[i].query));
        //     }
        // }
        // ************************************************** Have bug here, need to rewrite later (End) (2026-05-22 v0.1.3-r4) **************************************************

        // May introduce loop in the graph (2026-05-22 v0.1.3-r4)
        std::vector<SegReplace::Expansion> out;
        out.reserve(xs.size());
        for (size_t i = 0; i < xs.size(); ++i) {
            out.push_back(std::move(xs[i].query));
        }

        return out;
    };

    // Collect non-overlap pieces on one side of the overlap window.
    // left_side = true: collect source-side sequence before the overlap window
    // left_side = false: collect target-side sequence after the overlap window
    auto collect_non_overlap = [&ex](
        uint32_t seg_id, bool seg_rev,
        const std::vector<uint32_t>& cuts_vec,
        uint32_t win_beg, uint32_t win_end,
        bool left_side  // true: Only the left side of the window, false: Only take the right side
    ) -> std::vector<ExpansionPiece_> {
        std::vector<ExpansionPiece_> res;

        if (cuts_vec.size() > 1) res.reserve(cuts_vec.size() - 1);

        for (size_t i = 0; i + 1 < cuts_vec.size(); ++i) {
            const uint32_t beg = cuts_vec[i];
            const uint32_t end = cuts_vec[i + 1];

            const bool overlap_with_window =
                (!seg_rev && ((left_side && end > win_beg) || (!left_side && beg < win_end))) ||
                ( seg_rev && ((left_side && beg < win_end) || (!left_side && end > win_beg)));

            if (overlap_with_window) continue;

            SegReplace::Seg window_u128 = SegReplace::Interval::pack(seg_id, beg, end, seg_rev);
            res.push_back(ExpansionPiece_{window_u128, ex.query(window_u128)});
        }

        return res;
    };

    // Collect pieces inside an overlap window.
    auto collect_overlap = [&ex](
        uint32_t seg_id, bool seg_rev, const std::vector<uint32_t>& vec
    ) -> std::vector<ExpansionPiece_> {
        std::vector<ExpansionPiece_> res;
        if (vec.size() > 1) res.reserve(vec.size() - 1);

        for (size_t i = 0; i + 1 < vec.size(); ++i) {
            const uint32_t beg = vec[i];
            const uint32_t end = vec[i + 1];

            SegReplace::Seg window_u128 = SegReplace::Interval::pack(seg_id, beg, end, seg_rev);
            res.push_back(ExpansionPiece_{window_u128, ex.query(window_u128)});
        }

        return res;
    };

    // ------------------------------------------------ Preserve isolated unchanged nodes ------------------------------------------------
    keep_unused_nodes_.assign(nodes_.size(), 0);
    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        if (active_degree[sid] == 0 && trivial_segment(sid)) {
            keep_unused_nodes_[sid] = 1;
        }
    }

    // ------------------------------------------------ Expand and collect new edge chains ------------------------------------------------
    const uint64_t arc_num = arcs_.size();

    ProgressTracker prog(arc_num);  // Progress tracker for merging unitigs

    for (size_t ai = 0; ai < arc_num; ++ai) {
        prog.hit();

        const GfaArc& arc = arcs_[ai];
        if (arc.get_del() || arc.get_comp()) continue;

        const uint32_t v_seg_id = arc.get_source_segment_id();
        const uint32_t w_seg_id = arc.get_target_segment_id();
        const bool     v_is_rev = arc.get_source_is_reverse();
        const bool     w_is_rev = arc.get_target_is_reverse();

        const std::string& v_name = nodes_[v_seg_id].name;
        const std::string& w_name = nodes_[w_seg_id].name;

        const uint32_t v_length = nodes_[v_seg_id].length;
        const uint32_t w_length = nodes_[w_seg_id].length;

        if (arc.ov != arc.ow) {
            error_stream() << "Overlap values are not equal for " << v_name << " and " << w_name << "\n";
            std::exit(1);
        }

        /**
         * @brief Record the node that don't have any changes.
         * @date 2026-05-19
         * @version 0.1.3-r3
         */
        if ((arc.ov == 0 || arc.ov == INT32_MAX) && (arc.ow == 0 || arc.ow == INT32_MAX) && trivial_segment(v_seg_id) && trivial_segment(w_seg_id)) {
            seen_edges.insert(encode_edge_u64_(arc.get_source_vertex_id(), arc.get_target_vertex_id()));
            continue;
        }

        // ------------------------------------------------ Compute overlap windows ------------------------------------------------
        // Compute overlap windows on source and target
        auto [vb, ve] = v_overlap_pos_(v_length, arc.ov, v_is_rev);
        auto [wb, we] = w_overlap_pos_(w_length, arc.ow, w_is_rev);

        // Cuts positions
        const auto& v_cuts = cuts_[v_seg_id].v;
        const auto& w_cuts = cuts_[w_seg_id].v;

        // Extract cuts within overlap windows
        auto it_v_1 = std::lower_bound(v_cuts.begin(), v_cuts.end(), vb);
        auto it_v_2 = std::upper_bound(v_cuts.begin(), v_cuts.end(), ve);
        auto it_w_1 = std::lower_bound(w_cuts.begin(), w_cuts.end(), wb);
        auto it_w_2 = std::upper_bound(w_cuts.begin(), w_cuts.end(), we);

        std::vector<uint32_t> cut_in_overlap_v(it_v_1, it_v_2);
        std::vector<uint32_t> cut_in_overlap_w(it_w_1, it_w_2);

        // ------------------------------------------------ Collect source/target pieces ------------------------------------------------
        // no-overlap
        std::vector<ExpansionPiece_> source_pieces = collect_non_overlap(v_seg_id, v_is_rev, v_cuts, vb, ve, /*left_side=*/true);
        std::vector<ExpansionPiece_> target_pieces = collect_non_overlap(w_seg_id, w_is_rev, w_cuts, wb, we, /*left_side=*/false);

        // overlap
        std::vector<ExpansionPiece_> source_overlap_pieces = collect_overlap(v_seg_id, v_is_rev, cut_in_overlap_v);
        // std::vector<ExpansionPiece_> target_overlap_pieces = collect_overlap(w_seg_id, w_is_rev, cut_in_overlap_w);

        // Merge source expansions with overlap
        if (v_is_rev) {
            source_pieces.insert(
                source_pieces.begin(),
                source_overlap_pieces.begin(),
                source_overlap_pieces.end()
            );
            std::reverse(source_pieces.begin(), source_pieces.end());
        } else {
            source_pieces.insert(
                source_pieces.end(),
                source_overlap_pieces.begin(),
                source_overlap_pieces.end()
            );
        }

        if (w_is_rev) {
            std::reverse(target_pieces.begin(), target_pieces.end());
        }

        // ------------------------------------------------ Sanitize repeated replacement pieces (2026-05-22 v0.1.3-r4) ------------------------------------------------
        // If one replacement interval appears more than once inside the same source/target side,
        // fallback the whole affected piece range to original windows instead of query results.
        std::vector<SegReplace::Expansion> source_expansions = sanitize_expansions(source_pieces);
        std::vector<SegReplace::Expansion> target_expansions = sanitize_expansions(target_pieces);

        // ------------------------------------------------ Build expansion chain ------------------------------------------------
        // Build expansion chain
        SegReplace::Expansion node_expansion;
        {
            size_t total_sz = 0;
            for (const auto& seq : source_expansions) total_sz += seq.size();
            for (const auto& seq : target_expansions) total_sz += seq.size();
            node_expansion.reserve(total_sz);
        }

        for (const auto& seq : source_expansions) {
            node_expansion.insert(node_expansion.end(), seq.begin(), seq.end());
        }
        for (const auto& seq : target_expansions) {
            node_expansion.insert(node_expansion.end(), seq.begin(), seq.end());
        }

        // Disabled due to high time complexity (O(N^2)). But it may improve results slightly, it will be rewrited later. (2026-04-07)
        // refine_chain_(node_expansion, ex);

        // ------------------------------------------------ Emit edges from sanitized chain ------------------------------------------------
        // Add edges
        emit_node_expansion_edges_(
            node_expansion,
            seen_edges,
            new_node_count,
            new_edge_count,
            namer,
            v_name, w_name,
            v_is_rev, w_is_rev,
            std::make_optional(std::pair<uint32_t,uint32_t>{vb, ve}),
            std::make_optional(std::pair<uint32_t,uint32_t>{wb, we})
        );
    }

    // ------------------------------------------------ Replace graph edges ------------------------------------------------
    // Mark all non-complementary arcs for deletion
    for (auto& e : arcs_) {
        if (!e.get_comp()) e.set_del(true);
    }
    rebuild_after_edits();

    // Add edges
    for (const auto& k : seen_edges) {
        auto [v, w] = decode_edge_u64_(k);
        add_arc(v, w, 0, 0, -1, false);
    }
    rebuild_after_edits();

    // ------------------------------------------------ Summary ------------------------------------------------
    log_stream() << "  - New nodes added: " << new_node_count << "\n";
    log_stream() << "  - New edges added: " << new_edge_count << "\n" << "\n";
}


void GfaDeoverlapper::remove_unused_nodes_() {
    log_stream() << "Removing unused nodes ...\n";

    std::vector<bool> used(nodes_.size(), false);

    for (const auto& e : arcs_) {
        if (e.get_del() || e.get_comp()) continue;
        used[e.get_source_segment_id()] = true;
        used[e.get_target_segment_id()] = true;
    }

    uint64_t removed_num = 0;
    for (size_t i = 0; i < nodes_.size(); i++) {
        if (!used[i]) {
            /**
             * @brief Retain the node that don't have any changes.
             * @date 2026-05-19
             * @version 0.1.3-r3
             */
            if (i < keep_unused_nodes_.size() && keep_unused_nodes_[i]) continue;
            delete_segment(i);
            removed_num++;
        }
    }

    rebuild_after_edits();

    log_stream() << "  - Removed " << removed_num << " unused nodes.\n" << "\n";
}


// Main process
void GfaDeoverlapper::deoverlap(const std::string& prefix) {
    log_stream() << "Normalizing overlapping GFA to non-overlapping GFA ...\n" << "\n";
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

    // 4. group alignments by connected segments
    const auto align_groups = build_align_groups_();

    // 5. build cuts, propagate cuts, and prune abnormal cuts
    build_propagate_prune_cuts_(align_groups);
    if (DEBUG_ENABLED) print_cuts(cuts_);

    // 6. build rulemap from final cuts
    build_rulemap_(align_groups);
    SegReplace::Expander ex = build_SegReplace_(true);
    ex.save_map(prefix + ".deoverlap.map");  // save rulemap

    // 7. Apply rules: expand and rewire 0M edges
    expand_and_rewire_edges_(ex);

    // 8. Remove unused nodes
    remove_unused_nodes_();

    log_stream() << "Finished normalizing to non-overlapping GFA.\n" << "\n";
}
