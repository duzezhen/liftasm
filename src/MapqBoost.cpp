#include "../include/MapqBoost.hpp"
#include "../include/ThreadPool.hpp"
#include "../include/coordmap.hpp"
#include "../include/restriction_sites.hpp"
#include "../include/CIGAR.hpp"
#include "../include/progress_tracker.hpp"
#include "../include/logger.hpp"

#include <htslib/sam.h>
#include <htslib/hts.h>

#include <cctype>
#include <cstring>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cctype>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

namespace mapqboost {

static inline const uint32_t* bam_cigar_data(const bam1_t* b)
{
    return b ? bam_get_cigar(b) : nullptr;
}

static inline uint32_t bam_cigar_size(const bam1_t* b)
{
    return b ? b->core.n_cigar : 0u;
}

// 0 = single-end，1 = first-in-pair，2 = second-in-pair
static inline int read_end_id(const bam1_t* b)
{
    const uint16_t f = b->core.flag;
    if (!(f & BAM_FPAIRED))          return 0;
    if (f & BAM_FREAD1)              return 1;
    if (f & BAM_FREAD2)              return 2;
    return 0;
}

static inline int qname_cmp(const char* a, const char* b)
{
    const unsigned char* s1 = (const unsigned char*)a;
    const unsigned char* s2 = (const unsigned char*)b;

    while (*s1 && *s2) {
        const bool d1 = std::isdigit(*s1);
        const bool d2 = std::isdigit(*s2);

        if (d1 && d2) {
            // skip leading zeros
            const unsigned char* z1 = s1; while (*z1 == '0') ++z1;
            const unsigned char* z2 = s2; while (*z2 == '0') ++z2;

            // find end of digit runs
            const unsigned char* e1 = z1; while (std::isdigit(*e1)) ++e1;
            const unsigned char* e2 = z2; while (std::isdigit(*e2)) ++e2;

            const int n1 = (int)(e1 - z1);
            const int n2 = (int)(e2 - z2);

            if (n1 != n2) return (n1 < n2) ? -1 : 1;
            if (n1) {
                const int r = std::memcmp(z1, z2, (size_t)n1);
                if (r) return (r < 0) ? -1 : 1;
            }

            const int lz1 = (int)(z1 - s1);
            const int lz2 = (int)(z2 - s2);
            if (lz1 != lz2) return (lz1 < lz2) ? -1 : 1;

            s1 = e1;
            s2 = e2;
            continue;
        }

        if (*s1 != *s2) return (*s1 < *s2) ? -1 : 1;
        ++s1; ++s2;
    }
    if (*s1) return 1;
    if (*s2) return -1;
    return 0;
}

static inline std::string detect_write_mode(const std::string& out)
{
    if (out.empty() || out == "-") return "w";
    auto lower_suffix = [&](size_t k) {
        if (out.size() < k) return std::string();
        std::string suf = out.substr(out.size() - k);
        for (auto& c : suf) c = (char)std::tolower((unsigned char)c);
        return suf;
    };
    if (lower_suffix(5) == ".cram") return "wc";
    if (lower_suffix(4) == ".sam")  return "w";
    if (lower_suffix(4) == ".bam")  return "wb";
    return "wb";
}

static inline bool is_primary_like(const bam1_t* b) {
    const uint16_t f = b->core.flag;
    return ( (f & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0 );
}

static inline bool is_mapped(const bam1_t* b) {
    return b->core.tid >= 0 && (b->core.flag & BAM_FUNMAP) == 0;
}

static inline int32_t get_nm(const bam1_t* b) {
    if (const uint8_t* nm = bam_aux_get(b, "NM")) {
        int32_t v = (int32_t)bam_aux2i(nm);
        return (v >= 0) ? v : -1;
    } else if (const uint8_t* nm = bam_aux_get(b, "nM")) {
        int32_t v = (int32_t)bam_aux2i(nm);
        return (v >= 0) ? v : -1;
    }
    return -1;
}

// Lift exon-like blocks independently and accumulate all mapped intervals into eq
static inline void lift_blocks_to_eq(
    const coordmap::CoordMap& coormap_idx,
    std::string_view rname,
    const std::vector<std::pair<uint32_t,uint32_t>>& blocks,
    int cm_max_hops,
    uint32_t cm_max_fanout,
    uint32_t cm_min_len,
    double cm_min_frac,
    uint32_t cm_max_total_hits,
    std::unordered_map<std::string, std::vector<std::pair<uint32_t,uint32_t>>>& eq
) {
    // add source blocks
    {
        auto& src = eq[std::string(rname)];
        src.reserve(src.size() + blocks.size());
        for (const auto& iv : blocks) {
            if (iv.first < iv.second) src.push_back(iv);
        }
    }

    // lift each block separately
    for (const auto& iv : blocks) {
        const uint32_t beg = iv.first;
        const uint32_t end = iv.second;
        if (beg >= end) continue;

        auto hits = coormap_idx.map_range(rname, beg, end, cm_max_hops, cm_max_fanout, cm_min_len, cm_min_frac, cm_max_total_hits, true);

        for (auto& h : hits) {
            eq[coormap_idx.contig_name(h.ctg)].push_back({h.beg, h.end});
        }
    }
}

// Merge all interval lists in-place
static inline void merge_eq_inplace(
    std::unordered_map<std::string, std::vector<std::pair<uint32_t,uint32_t>>>& eq
) {
    for (auto& kv : eq) {
        auto& ivs = kv.second;
        if (ivs.empty()) continue;  // 2026-03-22
        std::sort(ivs.begin(), ivs.end());
        uint32_t cb = ivs[0].first, ce = ivs[0].second;
        size_t out_i = 0;
        for (size_t i = 1; i < ivs.size(); ++i) {
            auto [b, e] = ivs[i];
            if (b <= ce) ce = std::max(ce, e);
            else { ivs[out_i++] = {cb, ce}; cb = b; ce = e; }
        }
        ivs[out_i++] = {cb, ce};
        ivs.resize(out_i);
    }
}

static inline uint32_t overlapped_len(const std::vector<std::pair<uint32_t,uint32_t>>& ivs, uint32_t b0, uint32_t e0)
{
    if (b0 >= e0 || ivs.empty()) return 0;

    // lower_bound by iv.second > b0
    size_t lo = 0, hi = ivs.size();
    while (lo < hi) {
        size_t mid = (lo + hi) >> 1;
        if (ivs[mid].second <= b0) lo = mid + 1;
        else hi = mid;
    }

    uint32_t tot = 0;
    for (size_t i = lo; i < ivs.size(); ++i) {
        auto [ib, ie] = ivs[i];
        if (ib >= e0) break;
        uint32_t ob = (b0 > ib ? b0 : ib);
        uint32_t oe = (e0 < ie ? e0 : ie);
        if (ob < oe) tot += (oe - ob);
    }
    return tot;
}

static inline uint32_t overlapped_len_blocks(
    const std::vector<std::pair<uint32_t,uint32_t>>& eq_ivs,
    const std::vector<std::pair<uint32_t,uint32_t>>& query_blocks
) {
    uint32_t tot = 0;
    for (const auto& iv : query_blocks) {
        if (iv.first < iv.second) {
            tot += overlapped_len(eq_ivs, iv.first, iv.second);
        }
    }
    return tot;
}

static inline bool query_overlap(uint32_t b1, uint32_t e1, uint32_t b2, uint32_t e2)
{
    return (b1 < e2) && (b2 < e1);
}

static inline double query_overlap_frac(uint32_t b1, uint32_t e1, uint32_t b2, uint32_t e2)
{
    if (!(b1 < e2 && b2 < e1)) return 0.0;

    const uint32_t ob = std::max(b1, b2);
    const uint32_t oe = std::min(e1, e2);
    const uint32_t ov = oe - ob;

    const uint32_t l1 = e1 - b1;
    const uint32_t l2 = e2 - b2;
    const uint32_t den = std::min(l1, l2);

    return den ? (double)ov / (double)den : 0.0;
}

MapqBooster::MapqBooster(
    const coordmap::CoordMap& coormap_idx, const rsite::Index& rsite_idx,
    std::size_t batch_size, uint8_t mapq_low, uint8_t mapq_cap, bool name_check, 
    int cm_max_hops, uint32_t cm_max_fanout, uint32_t cm_min_len, double cm_min_frac, uint32_t cm_max_total_hits, 
    double sub_ovlp_frac,
    double K_mapq, double K_rs, double K_as, double K_ml, double K_nm, double W_mapq, double W_rs, double W_as, double W_ml, double W_nm, double close_as_eps,
    int threads, int io_threads
)
    : coormap_idx_(coormap_idx), rsite_idx_(rsite_idx)
    , batch_size_(batch_size >= 1 ? batch_size : 20000), mapq_low_(mapq_low), mapq_cap_(mapq_cap), name_check_(name_check)
    , cm_max_hops_(cm_max_hops), cm_max_fanout_(cm_max_fanout), cm_min_len_(cm_min_len), cm_min_frac_(cm_min_frac), cm_max_total_hits_(cm_max_total_hits)
    , sub_ovlp_frac_(sub_ovlp_frac)
    , K_mapq_(K_mapq), K_rs_(K_rs), K_as_(K_as), K_ml_(K_ml), K_nm_(K_nm), W_mapq_(W_mapq), W_rs_(W_rs), W_as_(W_as), W_ml_(W_ml), W_nm_(W_nm), close_as_eps_(close_as_eps)
    , threads_((threads <= 0) ? (int)std::max(1u, std::thread::hardware_concurrency()) : threads) , io_threads_(std::max(1, io_threads))
{}

// Score formula:
//   score = W_mapq * mapq_norm
//         + W_rs   * rs_good    (only if both flanking restriction sites are found)
//         + W_as   * as_norm
//         + W_ml   * ml_norm
//         - W_nm   * nm_pen
//
// Michaelis-Menten normalization:
//   mapq_norm = MAPQ  / (MAPQ  + K_mapq)
//   as_norm   = AS    / (AS    + K_as)    (0 if AS unavailable)
//   ml_norm   = ML    / (ML    + K_ml)
//   nm_pen    = NM    / (NM    + K_nm)
//
// restriction-site term:
//   closer  = min(dist_L, dist_R)
//   dvd     = 0.5*(dist_L + dist_R) - closer
//   rs_good = 0.7 / (1 + closer/K_rs) + 0.3 / (1 + dvd/K_rs)
//
static inline double score_formula(
    const SubRec& r, const sam_hdr_t* hdr, 
    const double K_mapq, const double K_rs, const double K_as, const double K_ml, const double K_nm,
    const double W_mapq, const double W_rs, const double W_as, const double W_ml, const double W_nm
) {
    if (K_mapq == 0.0 || K_rs == 0.0 || K_as == 0.0 || K_ml == 0.0 || K_nm == 0.0) {
        return (double)r.mapq;
    }

    const double mapq_norm = (double)r.mapq / ((double)r.mapq + K_mapq);
    const double mapq_term = W_mapq * mapq_norm;

    double rs_term = 0.0;
    if (r.rsite_found_L && r.rsite_found_R) {
        const double L = (double)r.rsite_dist_L;
        const double R = (double)r.rsite_dist_R;

        const double both   = 0.5 * (L + R);
        const double closer = std::min(L, R);
        const double dvd    = both - closer;

        const double closer_good = 1.0 / (1.0 + closer / K_rs);
        const double dvd_good    = 1.0 / (1.0 + dvd    / K_rs);

        const double rs_good = 0.7 * closer_good + 0.3 * dvd_good;
        rs_term = W_rs * rs_good;
    }

    const double as_norm = (r.as == INT32_MIN) ? 0.0 : ((double)std::max(0, r.as) / ((double)std::max(0, r.as) + K_as));

    const double ml_norm = (double)r.ml / ((double)r.ml + K_ml);

    const double nm_pen  = (double)std::max<int32_t>(0, r.nm) / ((double)std::max<int32_t>(0, r.nm) + K_nm);

    return  mapq_term + rs_term + W_as * as_norm + W_ml * ml_norm - W_nm * nm_pen;;
}

void MapqBooster::build_subgroups_by_query_overlap(
    const std::vector<bam1_t*>& group,
    std::vector<Subgroup>& out_subgroups, 
    const sam_hdr_t* hdr
) const
{
    struct Item { bam1_t* b; int end_id; uint32_t qbeg, qend; };

    std::vector<Item> items;
    items.reserve(group.size());

    for (bam1_t* b : group) {
        if (!b) continue;
        if (!is_mapped(b)) continue;

        uint32_t qbeg = 0, qend = 0;
        if (!CIGAR::query_interval_fwd(bam_cigar_data(b), bam_cigar_size(b), bam_is_rev(b), qbeg, qend)) continue;

        items.push_back(Item{b, read_end_id(b), qbeg, qend});
    }

    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b){
        if (a.end_id != b.end_id) return a.end_id < b.end_id;
        if (a.qbeg != b.qbeg) return a.qbeg < b.qbeg;
        return a.qend < b.qend;
    });

    out_subgroups.clear();
    if (items.empty()) return;

    std::vector<char> used(items.size(), 0);

    for (size_t i = 0; i < items.size(); ++i) {
        if (used[i]) continue;

        const int cur_end_id = items[i].end_id;

        used[i] = 1;
        Subgroup sub;
        sub.end_id = cur_end_id;
        sub.recs.reserve(8);

        uint32_t sub_beg = items[i].qbeg;
        uint32_t sub_end = items[i].qend;

        auto make_rec = [&](bam1_t* b, uint32_t qbeg, uint32_t qend) -> SubRec {
            SubRec r;
            r.b = b;
            r.qbeg = qbeg; r.qend = qend;

            r.tid = b->core.tid;
            if (r.tid >= 0) {
                r.rbeg = (uint32_t)b->core.pos;
                uint32_t span = CIGAR::ref_span(bam_cigar_data(b), bam_cigar_size(b));
                r.rend = r.rbeg + span;
            }

            r.mapq = b->core.qual;
            r.nm = get_nm(b);
            if (r.nm < 0) r.nm = INT32_MAX/2;
            if (const uint8_t* as_tmp = bam_aux_get(b, "AS")) r.as = (int32_t)bam_aux2i(as_tmp);
            r.ml = CIGAR::match_len(bam_cigar_data(b), bam_cigar_size(b));

            r.primary_like = is_primary_like(b);
            
            // restriction sites
            if (r.tid >= 0 && r.rend > r.rbeg) {
                const char* chr = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), r.tid);
                if (chr) {
                    rsite::Nearest L = rsite_idx_.nearest(std::string_view(chr), r.rbeg);
                    rsite::Nearest R = rsite_idx_.nearest(std::string_view(chr), r.rend - 1);

                    if (L.found) {
                        r.rsite_found_L = true;
                        r.rsite_pos_L  = (uint32_t)L.pos;
                        r.rsite_dist_L = (uint32_t)std::llabs((long long)L.dist);
                    }
                    if (R.found) {
                        r.rsite_found_R = true;
                        r.rsite_pos_R  = (uint32_t)R.pos;
                        r.rsite_dist_R = (uint32_t)std::llabs((long long)R.dist);
                    }
                }
            }

            r.score = score_formula(r, hdr, K_mapq_, K_rs_, K_as_, K_ml_, K_nm_, W_mapq_, W_rs_, W_as_, W_ml_, W_nm_);

            return r;
        };

        // First
        sub.recs.push_back(make_rec(items[i].b, items[i].qbeg, items[i].qend));

        // Collect
        for (size_t j = i + 1; j < items.size(); ++j) {
            if (used[j]) continue;

            if (items[j].end_id != cur_end_id) {
                if (items[j].end_id > cur_end_id) break;
                continue;
            }

            if (items[j].qbeg >= sub_end) break;

            if (sub_ovlp_frac_ > 0.0) {
                const double f = query_overlap_frac(sub_beg, sub_end, items[j].qbeg, items[j].qend);
                if (f < sub_ovlp_frac_) continue;
            }

            used[j] = 1;
            sub.recs.push_back(make_rec(items[j].b, items[j].qbeg, items[j].qend));
        }

        out_subgroups.push_back(std::move(sub));
    }

    // Re-rank each subgroup and put the best alignment at sub[0] (2026-03-03)
    MapqBooster::select_best_per_subgroup(out_subgroups, hdr);

    if (DEBUG_ENABLED) {print_group(group, out_subgroups, hdr);}
}

void MapqBooster::select_best_per_subgroup(std::vector<Subgroup>& out_subgroups, const sam_hdr_t* hdr) const
{
    if (!hdr) return;

    auto find_primary_for_end = [&](int end_id) -> bam1_t* {
        for (auto& sub : out_subgroups) {
            for (auto& r : sub.recs) {
                bam1_t* b = r.b;
                if (!b || !is_mapped(b)) continue;
                if (b->core.flag & BAM_FSUPPLEMENTARY) continue;
                if (read_end_id(b) != end_id) continue;
                if (is_primary_like(b)) return b;
            }
        }
        return nullptr;
    };

    auto tid_hap = [&](int32_t tid) -> int {
        if (tid < 0) return 0;
        const char* rn = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), tid);
        if (!rn) return 0;
        if (rn[0]=='h' && (rn[1]=='1' || rn[1]=='2')) return rn[1]-'0';
        return 0;
    };

    for (auto& sub : out_subgroups) {
        if (sub.recs.empty()) continue;

        const int end_id = sub.end_id;
        bam1_t* pri = find_primary_for_end(end_id);

        const int32_t pri_tid = pri ? pri->core.tid : -1;
        const int32_t pri_pos = pri ? pri->core.pos : 0;
        const int     pri_hap = pri ? tid_hap(pri_tid) : 0;

        auto better_than = [&](const SubRec& x, const SubRec& y) -> bool {
            // 0. Primary
            if (pri) {
                if (x.b == pri) return true;
                if (y.b == pri) return false;
            }

            // 1. Score
            if (x.score != y.score) return x.score > y.score;

            // 2. Same hap as primary
            if (pri_hap != 0) {
                const bool x_h = (x.tid >= 0 && tid_hap(x.tid) == pri_hap);
                const bool y_h = (y.tid >= 0 && tid_hap(y.tid) == pri_hap);
                if (x_h != y_h) return x_h;
            }

            // 3. Closer to primary if on same contig
            if (pri) {
                const bool x_on = (x.tid == pri_tid);
                const bool y_on = (y.tid == pri_tid);
                if (x_on && y_on) {
                    const long long dx = llabs((long long)x.rbeg - (long long)pri_pos);
                    const long long dy = llabs((long long)y.rbeg - (long long)pri_pos);
                    if (dx != dy) return dx < dy;
                }
            }

            // 4. Smaller tid/pos
            if (x.tid != y.tid) return x.tid < y.tid;
            return x.rbeg < y.rbeg;
        };

        std::stable_sort(sub.recs.begin(), sub.recs.end(), [&](const SubRec& a, const SubRec& b){ return better_than(a, b); });
    }
}

bool MapqBooster::should_boost_with_XA_(SubRec& cand, const sam_hdr_t* hdr) const
{
    const bam1_t* b = cand.b;
    if (!b || !hdr) return false;

    if (cand.mapq > mapq_low_) return false;

    const uint8_t* xa = bam_aux_get(b, "XA");
    if (!xa || bam_aux_type(xa) != 'Z') return false;

    const int32_t nm_pri  = cand.nm;
    const uint32_t ml_pri = cand.ml;

    const char* r_pri_c = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), b->core.tid);
    if (!r_pri_c) return false;
    const std::string_view r_pri(r_pri_c);

    const auto pri_blocks = CIGAR::ref_blocks(cand.rbeg, bam_cigar_data(b), bam_cigar_size(b));
    const uint32_t span_pri = CIGAR::total_block_bases(pri_blocks);
    if (!span_pri) return false;

    std::unordered_map<std::string, std::vector<std::pair<uint32_t,uint32_t>>> eq;
    eq.reserve(16);
    lift_blocks_to_eq(coormap_idx_, r_pri, pri_blocks, cm_max_hops_, cm_max_fanout_, cm_min_len_, cm_min_frac_, cm_max_total_hits_, eq);
    merge_eq_inplace(eq);

    const char* p = bam_aux2Z(xa);

    while (*p) {
        const char* rbeg = p;
        while (*p && *p != ',') ++p;
        if (!*p) break;
        std::string_view rname(rbeg, (size_t)(p - rbeg));
        ++p;

        if (*p == '+' || *p == '-') ++p;

        long pos1 = 0;
        while (std::isdigit((unsigned char)*p)) {
            pos1 = pos1 * 10 + (*p - '0');
            ++p;
        }
        if (*p != ',') break;
        ++p;

        const char* cbeg = p;
        while (*p && *p != ',') ++p;
        std::string_view cig(cbeg, (size_t)(p - cbeg));
        if (*p != ',') break;
        ++p;

        long nm_xa = 0;
        while (std::isdigit((unsigned char)*p)) {
            nm_xa = nm_xa * 10 + (*p - '0');
            ++p;
        }
        while (*p == ';') ++p;

        if (rname == r_pri) continue;  // 2026-03-22
        if (nm_pri < INT32_MAX/2 && nm_xa != nm_pri) continue;
        if (pos1 <= 0) continue;
        if (CIGAR::match_len(cig) != ml_pri) continue;

        const uint32_t beg = (uint32_t)(pos1 - 1);
        const auto xa_blocks = CIGAR::ref_blocks(beg, cig);
        if (xa_blocks.empty()) continue;

        auto it = eq.find(std::string(rname));
        cand.close_score_num++;  // 2026-03-22, This function (XA) don't have the score filed, so we use the NM implacitly as the socre.
        if (it == eq.end()) continue;  // 2026-03-22

        const uint32_t ov = overlapped_len_blocks(it->second, xa_blocks);
        cand.homol_frac += (double)ov / (double)span_pri;  // 2026-03-22, This function (XA) don't have the score filed, so we use the NM implacitly as the socre.
    }

    return true;
}

static inline bool has_close_as(const SubRec& best, const SubRec& r, const double eps = 0.01)
{
    const double as_best = (best.as == INT32_MIN) ? 0.0 : best.as;

    const double as_r = (r.as == INT32_MIN) ? 0.0 : r.as;

    return std::fabs(as_best - as_r) < (eps * as_best);
}

bool MapqBooster::should_boost_from_group_(Subgroup& sub, const sam_hdr_t* hdr) const
{
    if (!hdr) return false;
    if (sub.recs.empty()) return false;

    SubRec& cand = sub.recs[0];
    const bam1_t* b = cand.b;
    if (!b || !is_mapped(b)) return false;
    if (cand.mapq > mapq_low_) return false;

    const uint8_t* xa = bam_aux_get(b, "XA");
    if (xa && bam_aux_type(xa) == 'Z') {
        return should_boost_with_XA_(cand, hdr);
    }

    const char* ref_c_c = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), cand.tid);
    if (!ref_c_c) return false;
    const std::string_view ref_c(ref_c_c);

    const auto cand_blocks = CIGAR::ref_blocks(cand.rbeg, bam_cigar_data(b), bam_cigar_size(b));
    const uint32_t span = CIGAR::total_block_bases(cand_blocks);
    if (!span) return false;

    std::unordered_map<std::string,std::vector<std::pair<uint32_t,uint32_t>>> eq;
    eq.reserve(16);
    lift_blocks_to_eq(coormap_idx_, ref_c, cand_blocks, cm_max_hops_, cm_max_fanout_, cm_min_len_, cm_min_frac_, cm_max_total_hits_, eq);
    merge_eq_inplace(eq);

    const int end_id_c = sub.end_id;

    for (const auto& r : sub.recs)
    {
        const bam1_t* o = r.b;
        if (!o || o == b) continue;
        if (!is_mapped(o)) continue;
        if (read_end_id(o) != end_id_c) continue;

        const char* r2_c = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), r.tid);
        if (!r2_c) continue;
        const std::string_view r2(r2_c);

        const auto other_blocks = CIGAR::ref_blocks(r.rbeg, bam_cigar_data(o), bam_cigar_size(o));
        if (other_blocks.empty()) continue;

        bool close_as = has_close_as(cand, r, close_as_eps_);  // 2026-03-22
        if (close_as) cand.close_score_num++;   // 2026-03-22

        auto it = eq.find(std::string(r2));
        if (it == eq.end()) continue;

        const uint32_t ov = overlapped_len_blocks(it->second, other_blocks);
        if (close_as) cand.homol_frac += (double)ov / (double)span;  // 2026-03-22
    }

    return true;
}

// support = homol_frac / close_score_num
//   - p_wrong = 1 - support
//   - base_q  = -10 * log10(p_wrong)
static inline uint8_t recal_MAPQ(
    const SubRec& best,
    uint8_t mapq_cap = 60
) {
    const uint8_t oldq = best.mapq;

    if (best.close_score_num <= 0) return oldq;

    double support = best.homol_frac / (double)best.close_score_num;

    if (support < 0.0) support = 0.0;
    if (support > 1.0) support = 1.0;

    double p_wrong = 1.0 - support;

    const double eps = 1e-12;
    if (p_wrong < eps) p_wrong = eps;

    int newq = (int)std::llround(-10.0 * std::log10(p_wrong));

    if (newq < 0) newq = 0;
    if (newq > (int)mapq_cap) newq = (int)mapq_cap;

    return (uint8_t)newq;
}

void MapqBooster::process_group_(std::vector<bam1_t*>& group, const sam_hdr_t* hdr, Stats& st_local) const {
    if (!hdr) return;

    // 1. split into subgroups by (same end_id) and (query-overlap)
    std::vector<Subgroup> subgroups;  // First record in each subgroup is the "best" record by match length, MAPQ, and NM, and closest to primary if possible
    subgroups.reserve(8);
    build_subgroups_by_query_overlap(group, subgroups, hdr);

    // 2. Boost for each subgroup
    for (auto& sub : subgroups) {
        if (sub.recs.empty()) continue;
        SubRec& best = sub.recs[0];
        bam1_t* b = best.b;

        if (should_boost_from_group_(sub, hdr)) {
            const uint8_t newq = recal_MAPQ(best, mapq_cap_);
            if (newq != b->core.qual) {
                b->core.qual = newq;
                ++st_local.changed;
            }
        }
    }
}

// ---------------- run ----------------
int MapqBooster::run(const std::string& in_path_, const std::string& out_path_) const
{
    const std::string in_path  = (in_path_.empty()  ? "-" : in_path_);
    const std::string out_path = (out_path_.empty() ? "-" : out_path_);

    log_stream() << "Boosting MAPQ ...\n";

    // Overall stats
    Stats st;

    // Input
    htsFile* in_raw = sam_open(in_path.c_str(), "r");
    if (!in_raw) { error_stream() << in_path << ": No such file or directory\n"; std::exit(3); }
    std::unique_ptr<htsFile, void(*)(htsFile*)> in(in_raw, [](htsFile* f){ if (f) sam_close(f); });
    hts_set_threads(in.get(), io_threads_);

    sam_hdr_t* hdr_raw = sam_hdr_read(in.get());
    if (!hdr_raw) { error_stream() << "Header read failed"; std::exit(3); }
    std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t*)> hdr(hdr_raw, [](sam_hdr_t* p){ if (p) bam_hdr_destroy(p); });

    // Output
    std::string mode = detect_write_mode(out_path);
    htsFile* out_raw = sam_open(out_path.c_str(), mode.c_str());
    if (!out_raw) { error_stream() << out_path << ": No such file or directory\n"; std::exit(3); }
    std::unique_ptr<htsFile, void(*)(htsFile*)> out(out_raw, [](htsFile* f){ if (f) sam_close(f); });
    hts_set_threads(out.get(), io_threads_);

    if (sam_hdr_write(out.get(), hdr.get()) < 0) { error_stream() << "Header write failed"; std::exit(3); }

    struct Group {
        std::string qname;
        std::vector<bam1_t*> recs;
    };

    auto free_group = [](Group& g){
        for (auto* b : g.recs) bam_destroy1(b);
        g.recs.clear();
    };

    const std::size_t B = std::max<std::size_t>(1, batch_size_);
    std::vector<Group> batch;
    batch.reserve(1024);

    std::unique_ptr<bam1_t, void(*)(bam1_t*)> rec(bam_init1(), [](bam1_t* p){ if (p) bam_destroy1(p); });

    bool has_last_name = false;
    std::string last_name;

    Group cur;
    cur.recs.reserve(8);

    std::size_t batch_records = 0;

    const int T = std::max(1, threads_);
    ThreadPool pool(T);

    auto flush_batch = [&]() -> bool {
        if (batch.empty()) return true;

        std::vector<std::future<Stats>> futs;
        futs.reserve(batch.size());

        for (auto& g : batch) {
            futs.emplace_back(
                pool.submit([this, &g, hdr_ptr = hdr.get()]() -> Stats {
                    Stats local;
                    process_group_(g.recs, hdr_ptr, local);
                    return local;
                })
            );
        }

        for (auto& f : futs) {
            Stats x = f.get();
            st.changed += x.changed;
        }

        for (auto& g : batch) {
            for (auto* b : g.recs) {
                if (sam_write1(out.get(), hdr.get(), b) < 0) {
                    error_stream() << "Write error: " << g.qname << "\n";
                    std::exit(3);
                }
                ++st.written;
                bam_destroy1(b);
            }
            g.recs.clear();
        }

        batch.clear();
        batch_records = 0;
        return true;
    };

    auto tracker = ProgressTracker::Every(3000000);

    while (true) {
        int ret = sam_read1(in.get(), hdr.get(), rec.get());
        if (ret < 0) break;

        tracker.hit();
        ++st.total_in;

        bam1_t* dup = bam_dup1(rec.get());
        const char* qn_c = bam_get_qname(dup);
        const std::string qn = qn_c ? std::string(qn_c) : std::string();

        // Auto-detect QNAME-sorted input
        if (name_check_) {
            if (has_last_name) {
                if (qname_cmp(qn.c_str(), last_name.c_str()) < 0) {
                    error_stream() << "Input is not read-name sorted (QNAME not non-decreasing). Offending order: prev='" + last_name + "', curr='" + qn + "'. Please sort by name first, e.g.: samtools sort -n ...\n";
                    bam_destroy1(dup);
                    free_group(cur);
                    for (auto& g : batch) free_group(g);
                    batch.clear();
                    std::exit(4);
                }
            } else {
                has_last_name = true;
            }
            last_name = qn;
        }

        if (cur.qname.empty()) {
            cur.qname = qn;
        } else if (qn != cur.qname) {
            // finalize current group into batch
            batch_records += cur.recs.size();
            batch.push_back(std::move(cur));
            cur = Group{};
            cur.qname = qn;
            cur.recs.reserve(8);

            // flush if batch is large enough
            if (batch_records >= B) { flush_batch(); }
        }

        cur.recs.push_back(dup);
    }

    if (!cur.qname.empty()) {
        batch_records += cur.recs.size();
        batch.push_back(std::move(cur));
        cur = Group{};
    }

    flush_batch();
    tracker.finish();

    log_stream() << "  - processed=" << st.written << " changed=" << st.changed << "\n" << "\n";

    return 0;
}

void MapqBooster::print_group(
    const std::vector<bam1_t*>& group,
    const std::vector<Subgroup>& subs,
    const sam_hdr_t* hdr
) const
{
    if (!hdr || group.empty()) return;

    constexpr int W_IDX   = 10;
    constexpr int W_QINT  = 16;
    constexpr int W_FLAG  = 8;
    constexpr int W_MAPQ  = 6;
    constexpr int W_PRI   = 10;
    constexpr int W_AS    = 8;
    constexpr int W_ML    = 10;
    constexpr int W_NM    = 8;
    constexpr int W_RNAME = 15;
    constexpr int W_RPOS  = 24;
    constexpr int W_RDIS  = 24;
    constexpr int W_SCORE = 9;

    auto print_header = [&]() {
        debug_stream() << "  "
           << std::left
           << std::setw(W_IDX)   << "group"
           << std::setw(W_QINT)  << "q_interval"
           << std::setw(W_FLAG)  << "flag"
           << std::setw(W_MAPQ)  << "MAPQ"
           << std::setw(W_PRI)   << "primary"
           << std::setw(W_AS)    << "AS"
           << std::setw(W_ML)    << "ML"
           << std::setw(W_NM)    << "NM"
           << std::setw(W_RNAME) << "r_name"
           << std::setw(W_RPOS)  << "r_interval"
           << std::setw(W_RDIS)  << "rsite_dist(L,R)"
           << std::setw(W_SCORE) << "score"
           << "cigar"
           << "\n";
    };

    auto print_one = [&](const SubRec& r, int idx, int subidx){
        const bam1_t* b = r.b;
        if (!b) return;

        const char* rname = (r.tid >= 0)
            ? sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), r.tid)
            : "*";

        const std::string cig = CIGAR::to_string(bam_cigar_data(b), bam_cigar_size(b));

        if (idx == 0 && subidx == 0) print_header();

        std::string idxss_str = "(" + std::to_string(idx)  + "," + std::to_string(subidx) + ")";
        std::string qiss_str  = "[" + std::to_string(r.qbeg) + "," + std::to_string(r.qend) + ")";
        std::string riss_str  = "[" + std::to_string((int32_t)r.rbeg) + "," + std::to_string((int32_t)r.rend) + ")";
        std::string rdis_str  = ((r.rsite_found_L ) ? std::to_string(r.rsite_dist_L) : "*") + "," + (r.rsite_found_R ? std::to_string(r.rsite_dist_R) : "*");

        debug_stream() << "  " << std::left
           << std::setw(W_IDX)   << idxss_str
           << std::setw(W_QINT)  << qiss_str
           << std::setw(W_FLAG)  << b->core.flag
           << std::setw(W_MAPQ)  << (int)r.mapq
           << std::setw(W_PRI)   << (r.primary_like ? "T" : "F")
           << std::setw(W_AS)    << r.as
           << std::setw(W_ML)    << r.ml
           << std::setw(W_NM)    << r.nm
           << std::setw(W_RNAME) << (rname ? rname : "*")
           << std::setw(W_RPOS)  << riss_str
           << std::setw(W_RDIS)  << rdis_str
           << std::setw(W_SCORE) << int(r.score)
           << cig
           << "\n";
    };

    const char* qn0 = bam_get_qname(group[0]);
    debug_stream() << "QNAME: " << (qn0 ? qn0 : "") << "  (records=" << group.size() << ")" << "  (subgroups=" << subs.size() << ")" << "\n";

    for (size_t si = 0; si < subs.size(); ++si) {
        const auto& sg = subs[si];
        for (size_t k = 0; k < sg.recs.size(); ++k) {
            print_one(sg.recs[k], (int)si, (int)k);
        }
    }

    debug_stream() << "\n";
}

} // namespace mapqboost