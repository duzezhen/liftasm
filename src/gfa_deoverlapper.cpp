#include "../include/gfa_deoverlapper.hpp"
#include "../include/progress_tracker.hpp"

#include <set>
#include <deque>
#include <tuple>
#include <unordered_map>
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


bool GfaDeoverlapper::trim_small_overlap_(MmWfaHit& hit, const MmWfaHit& kept) const {
    const uint32_t hrb = std::min(hit.r_beg, hit.r_end);
    const uint32_t hre = std::max(hit.r_beg, hit.r_end);
    const uint32_t hqb = std::min(hit.q_beg, hit.q_end);
    const uint32_t hqe = std::max(hit.q_beg, hit.q_end);
    const uint32_t krb = std::min(kept.r_beg, kept.r_end);
    const uint32_t kre = std::max(kept.r_beg, kept.r_end);
    const uint32_t kqb = std::min(kept.q_beg, kept.q_end);
    const uint32_t kqe = std::max(kept.q_beg, kept.q_end);

    const uint32_t hit_len = std::min(hre - hrb, hqe - hqb);
    const uint32_t kept_len = std::min(kre - krb, kqe - kqb);
    const uint32_t r_overlap = std::min(hre, kre) > std::max(hrb, krb) ? std::min(hre, kre) - std::max(hrb, krb) : 0;
    const uint32_t q_overlap = std::min(hqe, kqe) > std::max(hqb, kqb) ? std::min(hqe, kqe) - std::max(hqb, kqb) : 0;
    const uint32_t overlap = std::max(r_overlap, q_overlap);

    if (overlap == 0 || hit_len < params_.trim_min_len || kept_len < params_.trim_min_len) return false;
    if (static_cast<double>(overlap) / hit_len > params_.trim_max_overlap) return false;

    const bool trim_front = hrb >= krb && hqb >= kqb;
    const bool trim_back = hre <= kre && hqe <= kqe;
    if (trim_front == trim_back) return false;

    std::vector<CIGAR::COp> ops = CIGAR::parse(hit.cigar);
    if (ops.empty()) return false;

    if (trim_back) std::reverse(ops.begin(), ops.end());

    const uint32_t need_r = r_overlap ? (trim_front ? kre - hrb : hre - krb) : 0;
    const uint32_t need_q = q_overlap ? (trim_front ? kqe - hqb : hqe - kqb) : 0;
    uint32_t used_r = 0, used_q = 0;
    size_t i = 0;

    while (i < ops.size() && (used_r < need_r || used_q < need_q)) {
        CIGAR::COp& op = ops[i];
        const bool rc = op.op == 'M' || op.op == 'D' || op.op == 'N' || op.op == '=' || op.op == 'X';
        const bool qc = op.op == 'M' || op.op == 'I' || op.op == '=' || op.op == 'X';
        const uint32_t left_r = used_r < need_r ? need_r - used_r : 0;
        const uint32_t left_q = used_q < need_q ? need_q - used_q : 0;
        const bool can_finish = (left_r == 0 || (rc && left_r <= op.len)) && (left_q == 0 || (qc && left_q <= op.len));
        const uint32_t take = can_finish ? std::max(left_r, left_q) : op.len;
        used_r += rc ? take : 0;
        used_q += qc ? take : 0;
        op.len -= take;
        if (op.len == 0) ++i;
    }
    if (used_r < need_r || used_q < need_q) return false;

    ops.erase(ops.begin(), ops.begin() + std::min(i, ops.size()));
    if (trim_back) std::reverse(ops.begin(), ops.end());

    if (trim_front) {
        hit.r_beg += used_r;
        hit.q_beg += used_q;
    } else {
        hit.r_end -= used_r;
        hit.q_end -= used_q;
    }

    if (ops.empty() || hit.r_end <= hit.r_beg || hit.q_end <= hit.q_beg) return false;
    if (std::min(hit.r_end, kre) > std::max(hit.r_beg, krb)) return false;
    if (std::min(hit.q_end, kqe) > std::max(hit.q_beg, kqb)) return false;
    hit.cigar = CIGAR::pack(ops);
    return hit.cigar != "*";
}

std::vector<GfaDeoverlapper::MmWfaHit> GfaDeoverlapper::filter_aligns_(
    std::vector<MmWfaHit> a,
    uint32_t ref_len,
    uint32_t qry_len
) const {
    auto rb = [](const MmWfaHit& x) { return std::min(x.r_beg, x.r_end); };
    auto re = [](const MmWfaHit& x) { return std::max(x.r_beg, x.r_end); };
    auto qb = [](const MmWfaHit& x) { return std::min(x.q_beg, x.q_end); };
    auto qe = [](const MmWfaHit& x) { return std::max(x.q_beg, x.q_end); };

    auto rspan = [&](const MmWfaHit& x) -> uint32_t {
        const uint32_t b = rb(x), e = re(x);
        return e > b ? e - b : 0u;
    };

    auto qspan = [&](const MmWfaHit& x) -> uint32_t {
        const uint32_t b = qb(x), e = qe(x);
        return e > b ? e - b : 0u;
    };

    auto len = [&](const MmWfaHit& x) -> uint32_t {
        return std::min(rspan(x), qspan(x));
    };

    auto ovlp = [](uint32_t l1, uint32_t r1, uint32_t l2, uint32_t r2) -> uint32_t {
        const uint32_t l = std::max(l1, l2);
        const uint32_t r = std::min(r1, r2);
        return r > l ? r - l : 0u;
    };

    const uint32_t denom = std::min(ref_len, qry_len);

    std::vector<MmWfaHit> passed;
    passed.reserve(a.size());

    for (auto& x : a) {
        if (x.cigar.empty() || x.cigar == "*") continue;
        if (x.mapq < params_.min_mapq) continue;
        if (CIGAR::match_ratio(x.cigar) < params_.min_match_ratio) continue;

        const uint32_t aln_len = len(x);
        if (aln_len == 0) continue;

        if (denom > 0) {
            const double ali_ratio = static_cast<double>(aln_len) / static_cast<double>(denom);
            if (ali_ratio < params_.min_ali_ratio) continue;
        }

        passed.emplace_back(std::move(x));
    }

    if (passed.size() <= 1) return passed;

    std::sort(passed.begin(), passed.end(), [&](const MmWfaHit& A, const MmWfaHit& B) {
        const uint32_t la = len(A);
        const uint32_t lb = len(B);
        if (la != lb) return la > lb;
        if (A.mapq != B.mapq) return A.mapq > B.mapq;
        return rb(A) < rb(B);
    });

    std::vector<MmWfaHit> kept;
    kept.reserve(passed.size());

    for (auto& cand : passed) {
        if (re(cand) <= rb(cand) || qe(cand) <= qb(cand)) continue;

        bool ok = true;
        for (const auto& k : kept) {
            const uint32_t r_overlap = ovlp(rb(cand), re(cand), rb(k), re(k));
            const uint32_t q_overlap = ovlp(qb(cand), qe(cand), qb(k), qe(k));
            if (r_overlap > 0 || q_overlap > 0) {
                if (!trim_small_overlap_(cand, k)) { ok = false; break; }
                continue;
            }

            const int64_t crk = int64_t(rb(cand));
            const int64_t cqk = int64_t(qb(cand));
            const int64_t krk = int64_t(rb(k));
            const int64_t kqk = int64_t(qb(k));

            if ((crk > krk && cqk <= kqk) || (crk < krk && cqk >= kqk)) {
                ok = false;
                break;
            }
        }

        if (!ok || cand.cigar.empty() || cand.cigar == "*") continue;
        if (CIGAR::match_ratio(cand.cigar) < params_.min_match_ratio) continue;
        const uint32_t aln_len = len(cand);
        if (aln_len == 0) continue;
        if (denom > 0 && static_cast<double>(aln_len) / static_cast<double>(denom) < params_.min_ali_ratio) continue;
        kept.emplace_back(std::move(cand));
    }

    return kept;
}

std::vector<GfaDeoverlapper::MmWfaHit> GfaDeoverlapper::align_short_wfa_(
    const std::string& v_name,
    const std::string& w_name,
    const std::string& v_seq_slice,
    const std::string& w_seq_slice
) {
    std::vector<MmWfaHit> hits;

    if (v_seq_slice.empty() || w_seq_slice.empty()) {
        return hits;
    }

    MmWfaHit h;
    h.r_beg = 0;
    h.r_end = static_cast<uint32_t>(v_seq_slice.size());
    h.q_beg = 0;
    h.q_end = static_cast<uint32_t>(w_seq_slice.size());

    if (v_seq_slice == w_seq_slice) {
        h.mapq = 60;
        h.cigar = std::to_string(v_seq_slice.size()) + "=";
        hits.emplace_back(std::move(h));
        return hits;
    }

    std::string cigar;
    int score = 0;

    seedExtend::run_wfa_fragment(
        std::string_view(v_seq_slice),
        std::string_view(w_seq_slice),
        cigar,
        score,
        alignment_options_.extend
    );

    if (cigar.empty() || cigar == "*") {
        return hits;
    }

    h.cigar = std::move(cigar);

    h.mapq = 60;

    hits.emplace_back(std::move(h));
    return hits;
}

std::vector<GfaDeoverlapper::MmWfaHit> GfaDeoverlapper::align_wfa_(
    const std::string& v_name,
    const std::string& w_name,
    const std::string& v_seq_slice,
    const std::string& w_seq_slice
) {
    std::vector<MmWfaHit> hits;

    // ---- Build a tiny minimizer index for v (reference) ----
    ExpandedSeqs sequences;
    sequences.names = {v_name};
    sequences.seqs = {v_seq_slice};
    sequences.right_seqs = {{""}};

    mmidx::MinimizerIndex idx(sequences, alignment_options_);
    idx.build_mm();

    aligner::Alignmenter aln(idx, sequences, alignment_options_);

    auto aligns = aln.produce_read(w_name, w_seq_slice, /*keep same strand only=*/true);
    if (aligns.empty()) {
        return hits;
    }

    hits.reserve(aligns.size());
    for (auto& a : aligns) {
        MmWfaHit h;
        h.r_beg = (uint32_t)a.r_beg;
        h.r_end = (uint32_t)a.r_end;
        h.q_beg = (uint32_t)a.q_beg;
        h.q_end = (uint32_t)a.q_end;

        h.mapq = a.MAPQ;

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
        h.mapq = 60;
        h.cigar = std::to_string(v_seq_slice.size()) + "=";
        hits.emplace_back(std::move(h));
        return hits;
    }

    const size_t min_len = std::min(v_seq_slice.size(), w_seq_slice.size());
    if (min_len < 4) return hits;

    mm_idxopt_t ipt;
    mm_mapopt_t opt;

    mm_set_opt(nullptr, &ipt, &opt);
    if (mm_set_opt("sr", &ipt, &opt) < 0) {
        error_stream() << "Invalid mm2 preset: sr" << "\n";
        std::exit(1);
    }

    ipt.k = (short)std::min<size_t>(std::max<size_t>(4, min_len), 11);
    ipt.w = 1;

    opt.flag |= MM_F_CIGAR;
    opt.flag |= MM_F_EQX;
    opt.best_n = (short)alignment_options_.anchor.max_kept;
    opt.zdrop = alignment_options_.extend.dyn_zdrop;

    const char* ref_seqs[1] = { v_seq_slice.c_str() };
    const char* ref_names[1] = { v_name.c_str() };

    mm_idx_t* mi = mm_idx_str(ipt.w, ipt.k, 0, ipt.bucket_bits, 1, ref_seqs, ref_names);
    if (!mi) return hits;

    // mm_mapopt_update(&opt, mi);

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

        MmWfaHit h;
        h.r_beg = (uint32_t)r.rs;
        h.r_end = (uint32_t)r.re;
        h.q_beg = (uint32_t)r.qs;
        h.q_end = (uint32_t)r.qe;

        h.mapq = r.mapq;

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

    const bool retry_short = (v_seq_slice.size() < 500 && w_seq_slice.size() < 500);

    mm_idxopt_t ipt;
    mm_mapopt_t opt;

    mm_set_opt(nullptr, &ipt, &opt);
    if (mm_set_opt(params_.mm2_preset.c_str(), &ipt, &opt) < 0) {
        error_stream() << "Invalid mm2 preset: " << params_.mm2_preset << "\n";
        std::exit(1);
    }

    ipt.k = (short)alignment_options_.chain.k;
    ipt.w = (short)alignment_options_.chain.w;

    opt.flag |= MM_F_CIGAR;
    opt.flag |= MM_F_EQX;  // output =/X

    opt.best_n = (short)alignment_options_.anchor.max_kept;
    opt.zdrop = alignment_options_.extend.dyn_zdrop;

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

        MmWfaHit h;
        h.r_beg = (uint32_t)r.rs;
        h.r_end = (uint32_t)r.re;
        h.q_beg = (uint32_t)r.qs;
        h.q_end = (uint32_t)r.qe;

        h.mapq = r.mapq;

        h.cigar = std::move(cigar);
        hits.emplace_back(std::move(h));
    }

    if (regs) {
        for (int i = 0; i < n_regs; ++i) free(regs[i].p);
        free(regs);
    }

    mm_idx_destroy(mi);

    // Retry short sequences with a smaller k-mer and window.
    if (hits.empty() && retry_short) {
        return align_short_mm2_(v_name, w_name, v_seq_slice, w_seq_slice);
    }

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
        uint8_t mapq,
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
        debug_stream() << "      - mapq = " << (int)mapq << "\n";
        debug_stream() << "      - CIGAR = " << cigar << "\n";
    };

    auto add_alignment = [&](
        uint32_t r_beg, uint32_t r_end,
        uint32_t q_beg, uint32_t q_end,
        std::string cigar,
        uint8_t mapq,
        std::vector<CIGAR::COp>&& ops,
        bool min_eq_checked,
        bool direct_equal
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
            mapq,
            std::move(ops)
        );
        outs.back().min_eq_checked = min_eq_checked;
        outs.back().direct_equal = direct_equal;
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
        if (DEBUG_ENABLED) print_alignment(0, L, 0, L, cigar, 60, ops);
        add_alignment(0, L, 0, L, std::move(cigar), 60, std::move(ops), false, true);
        return outs;
    }

    std::vector<MmWfaHit> hits;
    hits = align_mm2_(v_name, w_name, v_seq_slice, w_seq_slice);

    if (DEBUG_ENABLED && !hits.empty()) {
        debug_stream() << "  - Total raw alignments: " << hits.size() << "\n";
        for (auto& h : hits) {
            if (h.cigar.empty()) continue;
            auto ops = CIGAR::parse(h.cigar);
            if (DEBUG_ENABLED) print_alignment(h.r_beg, h.r_end, h.q_beg, h.q_end, h.cigar, h.mapq, ops);
        }
    }

    hits = filter_aligns_(std::move(hits), static_cast<uint32_t>(v_seq_slice.size()), static_cast<uint32_t>(w_seq_slice.size()));

    if (DEBUG_ENABLED && !hits.empty()) debug_stream() << "  - Total filtered alignments: " << hits.size() << "\n";
    for (auto& h : hits) {
        if (h.cigar.empty()) continue;
        auto ops = CIGAR::parse(h.cigar);
        if (!CIGAR::mask_short_matches(ops, static_cast<uint32_t>(params_.min_eq))) continue;
        if (DEBUG_ENABLED) print_alignment(h.r_beg, h.r_end, h.q_beg, h.q_end, h.cigar, h.mapq, ops);
        add_alignment(h.r_beg, h.r_end, h.q_beg, h.q_end, std::move(h.cigar), h.mapq, std::move(ops), true, false);
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
    ThreadPool pool(alignment_options_.align.threads);
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

void GfaDeoverlapper::dedup_aligns_(size_t begin)
{
    const size_t n_aln = bubble_aligns_.size();
    if (begin >= n_aln) return;

    struct SeenSpan {
        uint32_t a_beg, a_end;
        uint32_t b_beg, b_end;
    };

    struct CandidateBatch {
        int32_t idx;
        std::vector<size_t> alignments;
    };

    std::vector<BubbleAlignment> kept;
    kept.reserve(n_aln);
    for (size_t i = 0; i < begin; ++i) {
        kept.emplace_back(std::move(bubble_aligns_[i]));
    }

    log_stream() << "Deduplicating alignments ..." << "\n";
    log_stream() << "  - Number of alignments: " << n_aln << "\n";
    log_stream() << "  - Deduplication start: " << begin << "\n";

    std::vector<CandidateBatch> batches;
    std::unordered_map<int32_t, size_t> batch_of;
    batch_of.reserve(n_aln - begin);
    for (size_t i = begin; i < n_aln; ++i) {
        const int32_t idx = bubble_aligns_[i].idx;
        auto inserted = batch_of.emplace(idx, batches.size());
        if (inserted.second) batches.push_back({idx, {}});
        batches[inserted.first->second].alignments.push_back(i);
    }

    std::unordered_map<uint64_t, std::vector<SeenSpan>> accepted;
    accepted.reserve((n_aln - begin) * 2 + 8);
    size_t rejected_alignments = 0;
    size_t rejected_batches = 0;

    for (const CandidateBatch& batch : batches) {
        size_t kept_in_batch = 0;
        size_t rejected_in_batch = 0;

        for (size_t i : batch.alignments) {
            BubbleAlignment& aln = bubble_aligns_[i];
            const uint32_t seg_a = NodeHandle::get_segment_id(aln.v_a);
            const uint32_t seg_b = NodeHandle::get_segment_id(aln.v_b);

            if (seg_a >= nodes_.size() || seg_b >= nodes_.size()) {
                if (DEBUG_ENABLED) {
                    debug_stream() << "  - INVALID_SEGMENT_ID: " << "idx=" << aln.idx << " " << aln.name_a << ":" << aln.beg_a << "-" << aln.end_a << " " << aln.name_b << ":" << aln.beg_b << "-" << aln.end_b << "\n";
                }
                continue;
            }

            std::vector<std::string> roots_a;
            std::vector<std::string> roots_b;

            gfaName::collect_roots_from_name(nodes_[seg_a].name, roots_a);
            gfaName::collect_roots_from_name(nodes_[seg_b].name, roots_b);

            bool same_root = false;

            for (const auto& a : roots_a) {
                for (const auto& b : roots_b) {
                    if (a == b) {
                        same_root = true;
                        break;
                    }
                }
                if (same_root) break;
            }

            if (same_root) {
                if (DEBUG_ENABLED) {
                    debug_stream() << "  - SHARED_ROOT: " << "idx=" << aln.idx << " " << aln.name_a << ":" << aln.beg_a << "-" << aln.end_a << " " << aln.name_b << ":" << aln.beg_b << "-" << aln.end_b << "\n";
                }

                continue;
            }

            uint32_t a = seg_a;
            uint32_t b = seg_b;
            SeenSpan span{aln.beg_a, aln.end_a, aln.beg_b, aln.end_b};
            if (a > b) {
                std::swap(a, b);
                std::swap(span.a_beg, span.b_beg);
                std::swap(span.a_end, span.b_end);
            }
            const uint64_t key = (uint64_t(a) << 32) | uint64_t(b);

            bool conflict = false;
            const auto found = accepted.find(key);
            if (found != accepted.end()) {
                for (const SeenSpan& previous : found->second) {
                    const bool overlap_a = std::max(span.a_beg, previous.a_beg) < std::min(span.a_end, previous.a_end);
                    const bool overlap_b = std::max(span.b_beg, previous.b_beg) < std::min(span.b_end, previous.b_end);
                    if (overlap_a && overlap_b) {
                        conflict = true;
                        break;
                    }
                }
            }
            if (conflict) {
                ++rejected_alignments;
                ++rejected_in_batch;
                if (DEBUG_ENABLED) {
                    debug_stream() << "  - DROP_OVERLAPPING_ALIGNMENT: idx=" << aln.idx << " " << aln.name_a << ":" << aln.beg_a << "-" << aln.end_a << " " << aln.name_b << ":" << aln.beg_b << "-" << aln.end_b << "\n";
                }
                continue;
            }

            accepted[key].push_back(span);
            if (DEBUG_ENABLED) {
                debug_stream() << "  - KEEP: idx=" << aln.idx << " " << aln.name_a << ":" << aln.beg_a << "-" << aln.end_a << " " << aln.name_b << ":" << aln.beg_b << "-" << aln.end_b << "\n";
            }
            kept.emplace_back(std::move(aln));
            ++kept_in_batch;
        }

        if (kept_in_batch == 0 && rejected_in_batch > 0) ++rejected_batches;
    }

    bubble_aligns_ = std::move(kept);

    log_stream() << "  - Overlapping alignments removed: " << rejected_alignments << "\n";
    log_stream() << "  - Candidates fully removed: " << rejected_batches << "\n";
    log_stream() << "  - Number of alignments after deduplication: " << bubble_aligns_.size() << "\n\n";
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


bool GfaDeoverlapper::match_can_collapse_(const BubbleAlignment& align, const CIGAR::COp& op) const {
    return op.len >= params_.min_eq ||
           align.min_eq_checked ||
           align.direct_equal ||
           align.force_cuts;
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

                if (match_can_collapse_(align, op)) {
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

    for (int iter = 0; iter < params_.max_iters && !dirty.empty(); ++iter) {
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

                        // Only trusted matches may carry cuts between segments.
                        if (!match_can_collapse_(align, op)) break;

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
    if (params_.max_abnormal_cut_len <= 0 || params_.min_abnormal_cut_count <= 0) {
        return 0;
    }
    for (const size_t aln_id : group) {
        if (aln_id < bubble_aligns_.size() && bubble_aligns_[aln_id].force_cuts) return 0;
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

            if (seg_len <= params_.max_abnormal_cut_len) {
                const size_t abnormal_beg = i;

                while (i + 1 < cuts.size() && cuts[i + 1] - cuts[i] <= params_.max_abnormal_cut_len) {
                    ++i;
                }

                const size_t abnormal_seg_count = i - abnormal_beg;

                if (abnormal_seg_count >= params_.min_abnormal_cut_count) {
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

    ThreadPool pool(alignment_options_.align.threads);
    std::vector<std::future<std::pair<uint64_t, uint64_t>>> futs;
    futs.reserve(groups.size());
    ProgressTracker prog(groups.size());

    for (const auto& group : groups) {
        const auto* group_ptr = &group;
        futs.emplace_back(pool.submit([this, group_ptr, &prog]() {
            auto result = build_propagate_prune_cuts_run_(*group_ptr);
            prog.hit();
            return result;
        }));
    }

    uint64_t total_changed = 0;
    uint64_t total_removed = 0;

    for (auto& f : futs) {
        auto [changed, removed] = f.get();

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
                merge_and_dedup_cuts_with_set_(exist_cuts, exist_cut_set, extra_exist_cuts);  // 2026.04.10 update both cut list and cut set

                build_rules_from_pair_windows_(
                    rulemap,
                    /*trunk*/ leaf_seg_id, leaf_is_rev1, cut_in_overlap_leaf,
                    /*leaf*/ exist_seg_id, exist_is_rev, new_cuts_exist, new_exist_offsets,
                    /*leaf older*/ cut_in_overlap_exist, index
                );
            } else {
                // merge_and_dedup_cuts_(leaf_cuts, new_cuts_leaf);  // debug 2025-10-10, Update cut list
                merge_and_dedup_cuts_with_set_(leaf_cuts, leaf_cut_set, extra_leaf_cuts);  // 2026.04.10 update both cut list and cut set

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

        uint64_t consumed_a = 0;
        uint64_t consumed_b = 0;
        for (const CIGAR::COp& op : align.ops) {
            if (op.op == 'M' || op.op == '=' || op.op == 'X' || op.op == 'D') consumed_a += op.len;
            if (op.op == 'M' || op.op == '=' || op.op == 'X' || op.op == 'I') consumed_b += op.len;
        }
        if (consumed_a != align.end_a - align.beg_a || consumed_b != align.end_b - align.beg_b) continue;

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

                    if (match_can_collapse_(align, op)) {
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

    ThreadPool pool(alignment_options_.align.threads);
    std::deque<std::future<SegReplace::RuleMap>> futs;
    const size_t max_pending = std::max<size_t>(1, alignment_options_.align.threads * 2);
    ProgressTracker prog(groups.size());

    for (const auto& group : groups) {
        const auto* group_ptr = &group;
        futs.emplace_back(pool.submit([this, group_ptr, &cut_sets, &prog]() {
            SegReplace::RuleMap result = build_rulemap_run_(*group_ptr, cut_sets);
            prog.hit();
            return result;
        }));

        if (futs.size() >= max_pending) {
            SegReplace::RuleMap local = futs.front().get();
            futs.pop_front();
            for (auto& kv : local) {
                rulemap_.emplace(kv.first, std::move(kv.second));
            }
        }
    }
    while (!futs.empty()) {
        SegReplace::RuleMap local = futs.front().get();
        futs.pop_front();
        for (auto& kv : local) {
            rulemap_.emplace(kv.first, std::move(kv.second));
        }
    }

    prog.finish();
    pool.stop();

    normalize_all_cuts_();

    log_stream() << "  - Total rules generated: " << rulemap_.size() << "\n\n";
}

SegReplace::Expander GfaDeoverlapper::build_SegReplace_()
{
    log_stream() << "Building segment replacement rules ...\n";

    SegReplace::Expander ex(rulemap_, getAllSegmentNames(), params_.min_trans_len);
    ex.build_index();

    const size_t raw_rules = ex.index_view().size();
    log_stream() << "  - Raw expanded rules: " << raw_rules << "\n";

    const size_t removed_nonmonotonic = ex.filter_nonmonotonic_index();

    log_stream() << "  - Removed non-monotonic rules: " << removed_nonmonotonic << "\n";
    log_stream() << "  - Final rules: " << ex.index_view().size() << "\n\n";

    if (DEBUG_ENABLED) ex.print_index();
    if (DEBUG_ENABLED) rulemap_verify(nodes_, ex.index_view());

    return ex;
}

uint32_t GfaDeoverlapper::materialized_segment_(SegReplace::Seg interval) const {
    const uint32_t sid = static_cast<uint32_t>(SegReplace::Interval::seg_id(interval));
    if (sid >= nodes_.size()) return UINT32_MAX;

    const uint32_t beg = SegReplace::Interval::beg(interval);
    const uint32_t end = SegReplace::Interval::end(interval);
    if (beg == 0 && end == nodes_[sid].length) return sid;

    const gfaName namer;
    const std::string name = namer.format_interval_name(nodes_[sid].name, beg, end, false);
    const auto found = name_to_id_map_.find(name);
    return found == name_to_id_map_.end() ? UINT32_MAX : static_cast<uint32_t>(found->second);
}

void GfaDeoverlapper::rewrite_paths_(const SegReplace::Expander& expander) {
    std::vector<GfaPath> rewritten;
    rewritten.reserve(paths_.size());

    for (const GfaPath& source : paths_) {
        GfaPath path;
        path.name = source.name;
        for (const PathSegment& segment : source.segments) {
            const uint32_t sid = static_cast<uint32_t>(segment.node_id);
            if (sid >= nodes_.size()) continue;

            std::vector<uint32_t> boundaries;
            if (sid < cuts_.size() && cuts_[sid].v.size() >= 2) boundaries = cuts_[sid].v;
            else boundaries = {0, nodes_[sid].length};

            SegReplace::Expansion pieces;
            pieces.reserve(boundaries.size() - 1);
            for (size_t i = 1; i < boundaries.size(); ++i) {
                pieces.push_back(SegReplace::Interval::pack(
                    sid, boundaries[i - 1], boundaries[i], segment.is_reverse
                ));
            }
            if (segment.is_reverse) std::reverse(pieces.begin(), pieces.end());

            for (SegReplace::Seg piece : pieces) {
                for (SegReplace::Seg interval : expander.query(piece)) {
                    const uint32_t mapped = materialized_segment_(interval);
                    if (mapped == UINT32_MAX || mapped >= nodes_.size() || nodes_[mapped].deleted) {
                        error_stream() << "Path replacement was not materialized for segment '" << nodes_[sid].name << "'\n";
                        std::exit(1);
                    }
                    const bool reverse = SegReplace::Interval::is_reverse(interval);
                    if (!path.segments.empty() && path.segments.back().node_id == mapped && path.segments.back().is_reverse == reverse) continue;
                    path.segments.push_back({mapped, reverse});
                }
            }
        }
        if (!path.segments.empty()) rewritten.push_back(std::move(path));
    }
    paths_ = std::move(rewritten);
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


/* ================================================================================================================
 *                                         EXPAND REWIRER START
 * ================================================================================================================ */
GfaDeoverlapExpandRewirer::GfaDeoverlapExpandRewirer(
    GfaDeoverlapper& graph,
    const SegReplace::Expander& expander
) : graph_(graph), expander_(expander) {}

void GfaDeoverlapExpandRewirer::run() {
    log_stream() << "Expanding segments and rewiring edges ...\n";

    seen_edges_.reserve(graph_.arcs_.size() * 2 + 1);
    seen_edges_.max_load_factor(0.7f);

    build_active_degree_();
    build_piece_queries_();
    sanitize_local_replacements_();
    preserve_isolated_identity_nodes_();
    emit_isolated_expanded_nodes_();
    rewire_existing_edges_();
    replace_graph_edges_();

    log_stream() << "  - Local replacement pieces reverted: " << reverted_piece_count_ << "\n";
    log_stream() << "  - New nodes added: " << new_node_count_ << "\n";
    log_stream() << "  - New edges added: " << new_edge_count_ << "\n" << "\n";
}

void GfaDeoverlapExpandRewirer::build_active_degree_() {
    active_degree_.assign(graph_.nodes_.size(), 0);

    for (const auto& e : graph_.arcs_) {
        if (e.get_del() || e.get_comp()) continue;
        const uint32_t v = e.get_source_segment_id();
        const uint32_t w = e.get_target_segment_id();

        ++active_degree_[v];
        ++active_degree_[w];
    }
}

void GfaDeoverlapExpandRewirer::build_piece_queries_() {
    pieces_.clear();
    pieces_.reserve(graph_.nodes_.size() * 4 + 1);
    std::vector<SegReplace::Seg> repeated_pieces;

    for (uint32_t sid = 0; sid < graph_.nodes_.size(); ++sid) {
        if (graph_.nodes_[sid].deleted) continue;

        for (SegReplace::Seg piece : collect_all_pieces_(sid, false)) {
            SegReplace::Expansion query = expander_.query(piece);
            if (piece_query_has_repeat_(query)) repeated_pieces.push_back(piece);
            pieces_.try_emplace(piece, PieceQuery_{piece, std::move(query), false});
        }
        for (SegReplace::Seg piece : collect_all_pieces_(sid, true)) {
            SegReplace::Expansion query = expander_.query(piece);
            if (piece_query_has_repeat_(query)) repeated_pieces.push_back(piece);
            pieces_.try_emplace(piece, PieceQuery_{piece, std::move(query), false});
        }
    }

    for (SegReplace::Seg piece : repeated_pieces) {
        revert_piece_(piece);
    }
}

void GfaDeoverlapExpandRewirer::sanitize_local_replacements_() {
    for (const GfaArc& arc : graph_.arcs_) {
        if (arc.get_del() || arc.get_comp()) continue;
        sanitize_arc_repeats_(arc);
    }
}

void GfaDeoverlapExpandRewirer::sanitize_arc_repeats_(const GfaArc& arc) {
    if (arc.ov != arc.ow) return;

    const uint32_t v_seg_id = arc.get_source_segment_id();
    const uint32_t w_seg_id = arc.get_target_segment_id();
    const bool v_is_rev = arc.get_source_is_reverse();
    const bool w_is_rev = arc.get_target_is_reverse();

    const uint32_t v_length = graph_.nodes_[v_seg_id].length;
    const uint32_t w_length = graph_.nodes_[w_seg_id].length;

    auto [vb, ve] = v_overlap_pos_(v_length, arc.ov, v_is_rev);
    auto [wb, we] = w_overlap_pos_(w_length, arc.ow, w_is_rev);

    std::vector<SegReplace::Seg> source = collect_non_overlap_pieces_(v_seg_id, v_is_rev, vb, ve, true);
    std::vector<SegReplace::Seg> target = collect_non_overlap_pieces_(w_seg_id, w_is_rev, wb, we, false);
    std::vector<SegReplace::Seg> source_overlap = collect_overlap_pieces_(v_seg_id, v_is_rev, collect_overlap_cuts_(v_seg_id, vb, ve));
    std::vector<SegReplace::Seg> target_overlap = collect_overlap_pieces_(w_seg_id, w_is_rev, collect_overlap_cuts_(w_seg_id, wb, we));

    std::vector<SegReplace::Seg> local;
    local.reserve(source.size() + std::max(source_overlap.size(), target_overlap.size()) + target.size());

    local.insert(local.end(), source.begin(), source.end());
    local.insert(local.end(), source_overlap.begin(), source_overlap.end());
    local.insert(local.end(), target.begin(), target.end());
    sanitize_piece_repeats_(local);

    local.clear();
    local.insert(local.end(), source.begin(), source.end());
    local.insert(local.end(), target_overlap.begin(), target_overlap.end());
    local.insert(local.end(), target.begin(), target.end());
    sanitize_piece_repeats_(local);
}

void GfaDeoverlapExpandRewirer::sanitize_piece_repeats_(
    const std::vector<SegReplace::Seg>& local_pieces
) {
    if (local_pieces.empty()) return;

    std::unordered_map<SegReplace::Seg, SegReplace::Seg, SegReplace::U128Hash, SegReplace::U128Eq> first_owner;
    first_owner.reserve(local_pieces.size() * 4 + 1);

    for (SegReplace::Seg owner : local_pieces) {
        auto pit = pieces_.find(owner);
        if (pit == pieces_.end() || pit->second.use_original) continue;

        const SegReplace::Expansion& query = pit->second.query;
        if (is_identity_(query, owner)) continue;

        for (SegReplace::Seg s : query) {
            if (s == owner) {
                revert_piece_(owner);
                break;
            }

            auto it = first_owner.find(s);
            if (it == first_owner.end()) {
                first_owner.emplace(s, owner);
                continue;
            }

            revert_piece_(it->second);
            revert_piece_(owner);
            break;
        }
    }
}

void GfaDeoverlapExpandRewirer::preserve_isolated_identity_nodes_() {
    graph_.keep_unused_nodes_.assign(graph_.nodes_.size(), 0);

    for (uint32_t sid = 0; sid < graph_.nodes_.size(); ++sid) {
        if (sid < active_degree_.size() && active_degree_[sid] == 0 && trivial_segment_(sid)) {
            graph_.keep_unused_nodes_[sid] = 1;
        }
    }
}

void GfaDeoverlapExpandRewirer::emit_isolated_expanded_nodes_() {
    for (uint32_t sid = 0; sid < graph_.nodes_.size(); ++sid) {
        if (graph_.nodes_[sid].deleted) continue;
        if (sid >= active_degree_.size() || active_degree_[sid] != 0) continue;
        if (trivial_segment_(sid)) continue;

        std::vector<uint32_t> sample_ids;
        graph_.merge_sample_ids_into(sample_ids, graph_.nodes_[sid].sample_ids);

        Chain_ chain = build_chain_(
            collect_all_pieces_(sid, false),
            sample_ids,
            graph_.nodes_[sid].aux,
            {},
            {},
            {},
            {},
            {},
            {}
        );

        std::vector<uint32_t> vids = materialize_chain_(chain);
        for (uint32_t vid : vids) {
            const uint32_t keep_sid = vid >> 1;
            if (keep_sid >= graph_.keep_unused_nodes_.size()) {
                graph_.keep_unused_nodes_.resize(keep_sid + 1, 0);
            }
            graph_.keep_unused_nodes_[keep_sid] = 1;
        }

        for (size_t i = 0; i + 1 < vids.size(); ++i) {
            emit_bridge_edge_(vids[i], vids[i + 1]);
        }
    }
}

void GfaDeoverlapExpandRewirer::rewire_existing_edges_() {
    ProgressTracker prog(graph_.arcs_.size());

    for (const GfaArc& arc : graph_.arcs_) {
        prog.hit();

        if (arc.get_del() || arc.get_comp()) continue;

        const uint32_t v_seg_id = arc.get_source_segment_id();
        const uint32_t w_seg_id = arc.get_target_segment_id();
        const bool v_is_rev = arc.get_source_is_reverse();
        const bool w_is_rev = arc.get_target_is_reverse();

        const std::string& v_name = graph_.nodes_[v_seg_id].name;
        const std::string& w_name = graph_.nodes_[w_seg_id].name;

        if (arc.ov != arc.ow) {
            error_stream() << "Overlap values are not equal for " << v_name << " and " << w_name << "\n";
            std::exit(1);
        }

        if ((arc.ov == 0 || arc.ov == INT32_MAX) && (arc.ow == 0 || arc.ow == INT32_MAX) && trivial_segment_(v_seg_id) && trivial_segment_(w_seg_id)) {
            emit_bridge_edge_(arc.get_source_vertex_id(), arc.get_target_vertex_id());
            continue;
        }

        const uint32_t v_length = graph_.nodes_[v_seg_id].length;
        const uint32_t w_length = graph_.nodes_[w_seg_id].length;

        auto [vb, ve] = v_overlap_pos_(v_length, arc.ov, v_is_rev);
        auto [wb, we] = w_overlap_pos_(w_length, arc.ow, w_is_rev);

        std::vector<uint32_t> cut_in_overlap_v = collect_overlap_cuts_(v_seg_id, vb, ve);
        std::vector<uint32_t> cut_in_overlap_w = collect_overlap_cuts_(w_seg_id, wb, we);

        std::vector<uint32_t> source_sample_ids;
        std::vector<uint32_t> target_sample_ids;
        std::vector<uint32_t> overlap_sample_ids;

        graph_.merge_sample_ids_into(source_sample_ids, graph_.nodes_[v_seg_id].sample_ids);
        graph_.merge_sample_ids_into(overlap_sample_ids, graph_.nodes_[v_seg_id].sample_ids);
        graph_.merge_sample_ids_into(target_sample_ids, graph_.nodes_[w_seg_id].sample_ids);
        graph_.merge_sample_ids_into(overlap_sample_ids, graph_.nodes_[w_seg_id].sample_ids);
        GfaAux overlap_tags = graph_.nodes_[v_seg_id].aux;
        graph_.merge_variant_tags_into(overlap_tags, graph_.nodes_[w_seg_id].aux);

        std::vector<SegReplace::Seg> source = collect_non_overlap_pieces_(v_seg_id, v_is_rev, vb, ve, true);
        std::vector<SegReplace::Seg> target = collect_non_overlap_pieces_(w_seg_id, w_is_rev, wb, we, false);

        Chain_ chain = build_chain_(
            source,
            source_sample_ids,
            graph_.nodes_[v_seg_id].aux,
            collect_overlap_pieces_(v_seg_id, v_is_rev, cut_in_overlap_v),
            overlap_sample_ids,
            overlap_tags,
            target,
            target_sample_ids,
            graph_.nodes_[w_seg_id].aux
        );

        if (DEBUG_ENABLED) {
            debug_stream() << "Add arc: " << v_name << "(" << vb << "-" << ve << ":" << (v_is_rev ? "-" : "+") << ") and " << w_name << "(" << wb << "-" << we << ":" << (w_is_rev ? "-" : "+") << ")\n";
        }

        emit_chain_edges_(chain);

        Chain_ target_overlap_chain = build_chain_(
            source,
            source_sample_ids,
            graph_.nodes_[v_seg_id].aux,
            collect_overlap_pieces_(w_seg_id, w_is_rev, cut_in_overlap_w),
            overlap_sample_ids,
            overlap_tags,
            target,
            target_sample_ids,
            graph_.nodes_[w_seg_id].aux
        );

        emit_chain_edges_(target_overlap_chain);

        if (DEBUG_ENABLED) debug_stream() << "\n";
    }
}

void GfaDeoverlapExpandRewirer::replace_graph_edges_() {
    for (auto& e : graph_.arcs_) {
        if (!e.get_comp()) e.set_del(true);
    }
    graph_.rebuild_after_edits();

    for (const auto& k : seen_edges_) {
        auto [v, w] = decode_edge_u64_(k);
        graph_.add_arc(v, w, 0, 0, -1, false);
    }
    graph_.rebuild_after_edits();
}

bool GfaDeoverlapExpandRewirer::trivial_segment_(uint32_t sid) const {
    if (sid >= graph_.nodes_.size() || graph_.nodes_[sid].deleted) return false;

    const uint32_t len = graph_.nodes_[sid].length;
    if (sid >= graph_.cuts_.size()) return true;
    if (!(graph_.cuts_[sid].v.size() == 2 && graph_.cuts_[sid].v[0] == 0 && graph_.cuts_[sid].v[1] == len)) return false;

    SegReplace::Seg f = SegReplace::Interval::pack(sid, 0, len, false);
    SegReplace::Seg r = SegReplace::Interval::pack(sid, 0, len, true);

    return is_identity_(safe_expansion_(f), f) && is_identity_(safe_expansion_(r), r);
}

bool GfaDeoverlapExpandRewirer::piece_query_has_repeat_(const SegReplace::Expansion& query) const {
    std::unordered_set<SegReplace::Seg, SegReplace::U128Hash, SegReplace::U128Eq> seen;
    seen.reserve(query.size() * 2 + 1);

    for (SegReplace::Seg s : query) {
        if (!seen.insert(s).second) return true;
    }
    return false;
}

SegReplace::Expansion GfaDeoverlapExpandRewirer::safe_expansion_(SegReplace::Seg piece) const {
    auto it = pieces_.find(piece);
    if (it == pieces_.end() || it->second.use_original) {
        return SegReplace::Expansion{piece};
    }
    return it->second.query;
}

std::vector<SegReplace::Seg> GfaDeoverlapExpandRewirer::collect_all_pieces_(uint32_t sid, bool rev) const {
    std::vector<SegReplace::Seg> out;

    if (sid >= graph_.nodes_.size() || graph_.nodes_[sid].deleted) return out;

    const std::vector<uint32_t>* cuts = nullptr;
    std::vector<uint32_t> fallback;

    if (sid < graph_.cuts_.size() && graph_.cuts_[sid].v.size() >= 2) {
        cuts = &graph_.cuts_[sid].v;
    } else {
        fallback = {0, graph_.nodes_[sid].length};
        cuts = &fallback;
    }

    out.reserve(cuts->size() - 1);

    for (size_t i = 0; i + 1 < cuts->size(); ++i) {
        out.push_back(SegReplace::Interval::pack(sid, (*cuts)[i], (*cuts)[i + 1], rev));
    }

    if (rev) std::reverse(out.begin(), out.end());
    return out;
}

std::vector<SegReplace::Seg> GfaDeoverlapExpandRewirer::collect_non_overlap_pieces_(
    uint32_t sid,
    bool rev,
    uint32_t win_beg,
    uint32_t win_end,
    bool left_side
) const {
    std::vector<SegReplace::Seg> out;
    if (sid >= graph_.cuts_.size()) return collect_all_pieces_(sid, rev);

    const auto& cuts = graph_.cuts_[sid].v;
    if (cuts.size() > 1) out.reserve(cuts.size() - 1);

    for (size_t i = 0; i + 1 < cuts.size(); ++i) {
        const uint32_t beg = cuts[i];
        const uint32_t end = cuts[i + 1];

        bool keep = false;
        if (left_side) {
            keep = !rev ? (end <= win_beg) : (beg >= win_end);
        } else {
            keep = !rev ? (beg >= win_end) : (end <= win_beg);
        }

        if (keep) {
            out.push_back(SegReplace::Interval::pack(sid, beg, end, rev));
        }
    }

    if (rev) std::reverse(out.begin(), out.end());
    return out;
}

std::vector<SegReplace::Seg> GfaDeoverlapExpandRewirer::collect_overlap_pieces_(
    uint32_t sid,
    bool rev,
    const std::vector<uint32_t>& cuts
) const {
    std::vector<SegReplace::Seg> out;
    if (cuts.size() > 1) out.reserve(cuts.size() - 1);

    for (size_t i = 0; i + 1 < cuts.size(); ++i) {
        out.push_back(SegReplace::Interval::pack(sid, cuts[i], cuts[i + 1], rev));
    }

    if (rev) std::reverse(out.begin(), out.end());
    return out;
}

std::vector<uint32_t> GfaDeoverlapExpandRewirer::collect_overlap_cuts_(
    uint32_t sid,
    uint32_t beg,
    uint32_t end
) const {
    std::vector<uint32_t> cuts;
    cuts.reserve(8);
    cuts.push_back(beg);

    if (sid < graph_.cuts_.size()) {
        const auto& src = graph_.cuts_[sid].v;
        auto first = std::lower_bound(src.begin(), src.end(), beg);
        auto last = std::upper_bound(src.begin(), src.end(), end);
        cuts.insert(cuts.end(), first, last);
    }

    cuts.push_back(end);
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    return cuts;
}

GfaDeoverlapExpandRewirer::Chain_ GfaDeoverlapExpandRewirer::build_chain_(
    const std::vector<SegReplace::Seg>& source,
    const std::vector<uint32_t>& source_sample_ids,
    const GfaAux& source_tags,
    const std::vector<SegReplace::Seg>& overlap,
    const std::vector<uint32_t>& overlap_sample_ids,
    const GfaAux& overlap_tags,
    const std::vector<SegReplace::Seg>& target,
    const std::vector<uint32_t>& target_sample_ids,
    const GfaAux& target_tags
) const {
    Chain_ chain;

    append_piece_group_(source, source_sample_ids, source_tags, chain);
    append_piece_group_(overlap, overlap_sample_ids, overlap_tags, chain);
    append_piece_group_(target, target_sample_ids, target_tags, chain);

    return chain;
}

void GfaDeoverlapExpandRewirer::append_piece_group_(
    const std::vector<SegReplace::Seg>& pieces,
    const std::vector<uint32_t>& sample_ids,
    const GfaAux& variant_tags,
    Chain_& chain
) const {
    for (SegReplace::Seg piece : pieces) {
        SegReplace::Expansion safe = safe_expansion_(piece);
        for (SegReplace::Seg s : safe) {
            chain.expansion.push_back(s);
            chain.sample_ids.push_back(sample_ids);
            chain.variant_tags.push_back(variant_tags);
        }
    }
}

bool GfaDeoverlapExpandRewirer::revert_piece_(SegReplace::Seg piece) {
    bool changed = mark_piece_original_(piece);
    changed = mark_piece_original_(SegReplace::Interval::toggle_strand(piece)) || changed;

    if (changed) {
        expander_.remove_from_index_by_key(piece);
    }

    return changed;
}

bool GfaDeoverlapExpandRewirer::mark_piece_original_(SegReplace::Seg piece) {
    auto it = pieces_.find(piece);
    if (it == pieces_.end()) return false;
    if (it->second.use_original) return false;
    if (is_identity_(it->second.query, it->second.original)) return false;

    it->second.use_original = true;
    ++reverted_piece_count_;
    return true;
}

std::vector<uint32_t> GfaDeoverlapExpandRewirer::materialize_chain_(const Chain_& chain) {
    if (chain.expansion.size() != chain.sample_ids.size() || chain.expansion.size() != chain.variant_tags.size()) {
        error_stream() << "metadata size mismatch with node expansion\n";
        std::exit(1);
    }

    std::vector<uint32_t> vids;
    vids.reserve(chain.expansion.size());

    for (size_t i = 0; i < chain.expansion.size(); ++i) {
        SegReplace::Seg s = chain.expansion[i];
        uint32_t sid = materialize_segment_(s, chain.sample_ids[i], chain.variant_tags[i]);
        bool rev = SegReplace::Interval::is_reverse(s);
        vids.push_back((sid << 1) | (rev ? 1u : 0u));
    }

    return vids;
}

uint32_t GfaDeoverlapExpandRewirer::materialize_segment_(
    SegReplace::Seg s,
    const std::vector<uint32_t>& sample_ids,
    const GfaAux& variant_tags
) {
    const uint32_t original_sid = static_cast<uint32_t>(SegReplace::Interval::seg_id(s));

    if (original_sid >= graph_.nodes_.size()) {
        error_stream() << "Invalid replacement segment id: " << original_sid << "\n";
        std::exit(1);
    }

    const bool full_interval = SegReplace::Interval::beg(s) == 0 && SegReplace::Interval::end(s) == graph_.nodes_[original_sid].length;

    if (full_interval) {
        graph_.merge_sample_ids_into(graph_.nodes_[original_sid].sample_ids, sample_ids);
        graph_.merge_variant_tags_into(graph_.nodes_[original_sid].aux, variant_tags);
        return original_sid;
    }

    const std::string& base = graph_.getNodeName(original_sid);
    const uint32_t beg = SegReplace::Interval::beg(s);
    const uint32_t end = SegReplace::Interval::end(s);
    const std::string seg_name = namer_.format_interval_name(base, beg, end, false);
    const bool inherited_complex = graph_.nodes_[original_sid].is_complex;

    bool is_new = false;
    uint32_t sid = graph_.get_or_add_segment(seg_name, &is_new);

    graph_.merge_sample_ids_into(graph_.nodes_[sid].sample_ids, sample_ids);
    graph_.nodes_[sid].is_complex = graph_.nodes_[sid].is_complex || inherited_complex;
    graph_.merge_variant_tags_into(graph_.nodes_[sid].aux, variant_tags);

    if (is_new) {
        std::string seg_seq = GfaGraph::slice_seq_or_star_(graph_.nodes_, original_sid, beg, end, false);

        GfaNode& n = graph_.nodes_[sid];
        n.name = seg_name;
        n.sequence = seg_seq;
        n.length = seg_seq != "*" ? static_cast<uint32_t>(seg_seq.length()) : 0;
        n.deleted = false;
        graph_.total_segment_length_ += n.length;
        ++new_node_count_;
    }

    return sid;
}

void GfaDeoverlapExpandRewirer::emit_chain_edges_(const Chain_& chain) {
    std::vector<uint32_t> vids = materialize_chain_(chain);

    for (size_t i = 0; i + 1 < vids.size(); ++i) {
        emit_bridge_edge_(vids[i], vids[i + 1]);

        if (DEBUG_ENABLED) {
            const uint32_t v_sid = vids[i] >> 1;
            const uint32_t w_sid = vids[i + 1] >> 1;
            debug_stream() << "  - " << graph_.getNodeName(v_sid) << ((vids[i] & 1u) ? "-" : "+") << " -> " << graph_.getNodeName(w_sid) << ((vids[i + 1] & 1u) ? "-" : "+") << "\n";
        }
    }
}

void GfaDeoverlapExpandRewirer::emit_bridge_edge_(uint32_t from, uint32_t to) {
    if (seen_edges_.insert(encode_edge_u64_(from, to)).second) {
        ++new_edge_count_;
    }
}
/* ================================================================================================================
 *                                         EXPAND REWIRER END
 * ================================================================================================================ */


void GfaDeoverlapper::expand_and_rewire_edges_(
    const SegReplace::Expander& ex
) {
    GfaDeoverlapExpandRewirer rewirer(*this, ex);
    rewirer.run();
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
        if (nodes_[i].deleted) continue;
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
    log_stream() << "Normalizing overlapping GFA to non-overlapping GFA ...\n\n";
    
    finalize_();
    prune_overlaps_();
    initialize_cuts_();
    overlaps_align_();
    dedup_aligns_(bubble_aligns_.size());

    const auto align_groups = build_align_groups_();
    build_propagate_prune_cuts_(align_groups);
    if (DEBUG_ENABLED) print_cuts(cuts_);
    build_rulemap_(align_groups);

    SegReplace::Expander ex = build_SegReplace_();
    ex.save_map(prefix + ".deoverlap.map");
    expand_and_rewire_edges_(ex);
    remove_unused_nodes_();

    log_stream() << "Finished normalizing to non-overlapping GFA.\n" << "\n";
}
