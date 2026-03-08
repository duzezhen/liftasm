#include "../include/MapqBoost.hpp"
#include "../include/coordmap.hpp"
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
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

namespace mapqboost {

// ---------------- helpers ----------------

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
    const coordmap::CoordMap& map,
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

        auto hits = map.map_range(rname, beg, end, cm_max_hops, cm_max_fanout, cm_min_len, cm_min_frac, cm_max_total_hits, true);

        for (auto& h : hits) {
            eq[map.contig_name(h.ctg)].push_back({h.beg, h.end});
        }
    }
}

// Merge all interval lists in-place
static inline void merge_eq_inplace(
    std::unordered_map<std::string, std::vector<std::pair<uint32_t,uint32_t>>>& eq
) {
    for (auto& kv : eq) {
        auto& ivs = kv.second;
        if (ivs.empty()) return;
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
    const coordmap::CoordMap& map, 
    std::size_t batch_size, uint8_t mapq_low, uint8_t mapq_new, double min_frac, int min_equiv_contigs, bool name_check, 
    int cm_max_hops, uint32_t cm_max_fanout, uint32_t cm_min_len, double cm_min_frac, uint32_t cm_max_total_hits, 
    double sub_ovlp_frac,
    int threads, int io_threads
)
    : map_(map)
    , batch_size_(batch_size >= 1 ? batch_size : 20000), mapq_low_(mapq_low), mapq_new_(mapq_new), min_frac_(min_frac), min_equiv_contigs_(min_equiv_contigs), name_check_(name_check)
    , cm_max_hops_(cm_max_hops), cm_max_fanout_(cm_max_fanout), cm_min_len_(cm_min_len), cm_min_frac_(cm_min_frac), cm_max_total_hits_(cm_max_total_hits)
    , sub_ovlp_frac_(sub_ovlp_frac)
    , threads_((threads <= 0) ? (int)std::max(1u, std::thread::hardware_concurrency()) : threads) , io_threads_(std::max(1, io_threads))
{}

void MapqBooster::build_subgroups_by_query_overlap(
    const std::vector<bam1_t*>& group,
    std::vector<std::vector<bam1_t*>>& out_subgroups, 
    const sam_hdr_t* hdr
) const
{
    struct Item { bam1_t* b; int end_id; uint32_t qb, qe; };

    std::vector<Item> items;
    items.reserve(group.size());

    for (bam1_t* b : group) {
        if (!b) continue;
        if (!is_mapped(b)) continue;

        uint32_t qb = 0, qe = 0;
        if (!CIGAR::query_interval_fwd(bam_cigar_data(b), bam_cigar_size(b), bam_is_rev(b), qb, qe)) continue;

        items.push_back(Item{b, read_end_id(b), qb, qe});
    }

    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b){
        if (a.end_id != b.end_id) return a.end_id < b.end_id;
        if (a.qb != b.qb) return a.qb < b.qb;
        return a.qe < b.qe;
    });

    auto find_primary_for_end = [&](int end_id) -> bam1_t* {
        for (bam1_t* b : group) {
            if (!b) continue;
            if (!is_mapped(b)) continue;
            if (b->core.flag & BAM_FSUPPLEMENTARY) continue;
            if (read_end_id(b) != end_id) continue;
            if (is_primary_like(b)) return b;
        }
        return nullptr;
    };

    auto better = [&](bam1_t* a, bam1_t* b) -> bool {
        // Match length
        uint32_t mla = CIGAR::match_len(bam_cigar_data(a), bam_cigar_size(a));
        uint32_t mlb = CIGAR::match_len(bam_cigar_data(b), bam_cigar_size(b));
        if (mla != mlb) return mla > mlb;

        // MAPQ
        if (a->core.qual != b->core.qual) return a->core.qual > b->core.qual;

        // NM
        int32_t na = get_nm(a); if (na < 0) na = INT32_MAX/2;
        int32_t nb = get_nm(b); if (nb < 0) nb = INT32_MAX/2;
        if (na != nb) return na < nb;

        // final stable tiebreak
        if (a->core.tid != b->core.tid) return a->core.tid < b->core.tid;
        return a->core.pos < b->core.pos;
    };

    out_subgroups.clear();
    if (items.empty()) return;

    std::vector<char> used(items.size(), 0);

    bam1_t* pri = nullptr;

    for (size_t i = 0; i < items.size(); ++i) {
        if (used[i]) continue;

        const int cur_end_id = items[i].end_id;

        if (!pri || read_end_id(pri) != cur_end_id) {
            pri = find_primary_for_end(cur_end_id);
        }

        used[i] = 1;
        std::vector<bam1_t*> sub;
        sub.reserve(8);
        sub.push_back(items[i].b);

        uint32_t sub_beg = items[i].qb;
        uint32_t sub_end = items[i].qe;

        for (size_t j = i + 1; j < items.size(); ++j) {
            if (used[j]) continue;

            if (items[j].end_id != cur_end_id) {
                if (items[j].end_id > cur_end_id) break;
                continue;
            }

            if (items[j].qb >= sub_end) break;

            if (sub_ovlp_frac_ > 0.0) {
                const double f = query_overlap_frac(sub_beg, sub_end, items[j].qb, items[j].qe);
                if (f < sub_ovlp_frac_) continue;
            }

            sub.push_back(items[j].b);
            used[j] = 1;
        }

        out_subgroups.push_back(std::move(sub));
    }

    // Re-rank each subgroup and put the best alignment at sub[0] (2026-03-03)
    MapqBooster::select_best_per_subgroup(out_subgroups, hdr);

    if (DEBUG_ENABLED) {print_group(group, out_subgroups, hdr);}
}

void MapqBooster::select_best_per_subgroup(std::vector<std::vector<bam1_t*>>& out_subgroups, const sam_hdr_t* hdr) const
{
    if (!hdr) return;

    auto find_primary_for_end = [&](int end_id) -> bam1_t* {
        for (auto& sub : out_subgroups) {
            for (bam1_t* b : sub) {
                if (!b || !is_mapped(b)) continue;
                if (b->core.flag & BAM_FSUPPLEMENTARY) continue;
                if (read_end_id(b) != end_id) continue;
                if (is_primary_like(b)) return b;
            }
        }
        return nullptr;
    };

    auto get_as = [&](bam1_t* b) -> int32_t {
        if (!b) return INT32_MIN;
        if (const uint8_t* as_tmp = bam_aux_get(b, "AS")) return (int32_t)bam_aux2i(as_tmp);
        return INT32_MIN;
    };

    auto get_nm_norm = [&](bam1_t* b) -> int32_t {
        int32_t nm = get_nm(b);
        return (nm >= 0) ? nm : INT32_MAX/2;
    };

    auto tid_hap = [&](int32_t tid) -> int {
        if (tid < 0) return 0;
        const char* rn = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), tid);
        if (!rn) return 0;
        if (rn[0]=='h' && (rn[1]=='1' || rn[1]=='2')) return rn[1]-'0';
        return 0;
    };

    for (auto& sub : out_subgroups) {
        if (sub.empty()) continue;

        const int end_id = read_end_id(sub[0]);
        bam1_t* pri = find_primary_for_end(end_id);

        const int32_t pri_tid = pri ? pri->core.tid : -1;
        const int32_t pri_pos = pri ? pri->core.pos : 0;
        const int     pri_hap = pri ? tid_hap(pri_tid) : 0;

        auto better_than = [&](bam1_t* x, bam1_t* y) -> bool {
            if (!x) return false;
            if (!y) return true;

            // 0. primary
            if (pri) {
                if (x == pri) return true;
                if (y == pri) return false;
            }

            // 1. AS
            const int32_t ax = get_as(x), ay = get_as(y);
            if (ax != ay) return ax > ay;

            // 2. MAPQ
            if (x->core.qual != y->core.qual) return x->core.qual > y->core.qual;

            // 3. ML
            const uint32_t mlx = CIGAR::match_len(bam_cigar_data(x), bam_cigar_size(x));
            const uint32_t mly = CIGAR::match_len(bam_cigar_data(y), bam_cigar_size(y));
            if (mlx != mly) return mlx > mly;

            // 4. NM
            const int32_t nx = get_nm_norm(x), ny = get_nm_norm(y);
            if (nx != ny) return nx < ny;

            // 5. same hap as primary
            if (pri_hap != 0) {
                const bool x_h = (is_mapped(x) && tid_hap(x->core.tid) == pri_hap);
                const bool y_h = (is_mapped(y) && tid_hap(y->core.tid) == pri_hap);
                if (x_h != y_h) return x_h;
            }

            // 6. closer to primary
            if (pri) {
                const bool x_on = (x->core.tid == pri_tid);
                const bool y_on = (y->core.tid == pri_tid);
                if (x_on && y_on) {
                    const long long dx = llabs((long long)x->core.pos - (long long)pri_pos);
                    const long long dy = llabs((long long)y->core.pos - (long long)pri_pos);
                    if (dx != dy) return dx < dy;
                } else if (x_on != y_on) {
                    return x_on;
                }
            }

            // 7. smaller tid/pos
            if (x->core.tid != y->core.tid) return x->core.tid < y->core.tid;
            return x->core.pos < y->core.pos;
        };

        std::stable_sort(sub.begin(), sub.end(), [&](bam1_t* a, bam1_t* b){ return better_than(a, b); });
    }
}

bool MapqBooster::should_boost_with_XA_(const bam1_t* b, const sam_hdr_t* hdr) const
{
    if (!is_primary_like(b)) return false;

    const uint8_t* xa = bam_aux_get(b, "XA");
    if (!xa || bam_aux_type(xa) != 'Z') return false;

    const int32_t nm_pri = get_nm(b);
    const uint32_t ml_pri = CIGAR::match_len(bam_cigar_data(b), bam_cigar_size(b));

    const char* r_pri_c = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), b->core.tid);
    if (!r_pri_c) return false;
    const std::string_view r_pri(r_pri_c);

    const uint32_t beg_pri = b->core.pos;

    const auto pri_blocks = CIGAR::ref_blocks(beg_pri, bam_cigar_data(b), bam_cigar_size(b));
    const uint32_t span_pri = CIGAR::total_block_bases(pri_blocks);
    if (!span_pri) return false;

    std::unordered_map<std::string, std::vector<std::pair<uint32_t,uint32_t>>> eq;
    eq.reserve(16);
    lift_blocks_to_eq(map_, r_pri, pri_blocks, cm_max_hops_, cm_max_fanout_, cm_min_len_, cm_min_frac_, cm_max_total_hits_, eq);
    merge_eq_inplace(eq);

    const uint32_t min_cov = (uint32_t)std::ceil(min_frac_ * span_pri);

    const char* p = bam_aux2Z(xa);
    int ok_ctg = 0;

    while (*p) {
        const char* rbeg = p;
        while (*p && *p != ',') ++p;
        if (!*p) break;
        std::string_view rname(rbeg, (size_t)(p - rbeg));
        ++p;

        if (rname == r_pri) return false;

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

        if (nm_pri >= 0 && nm_xa != nm_pri) continue;
        if (pos1 <= 0) continue;
        if (CIGAR::match_len(cig) != ml_pri) continue;

        const uint32_t beg = (uint32_t)(pos1 - 1);
        const auto xa_blocks = CIGAR::ref_blocks(beg, cig);
        if (xa_blocks.empty()) continue;

        auto it = eq.find(std::string(rname));
        if (it == eq.end()) return false;

        const uint32_t ov = overlapped_len_blocks(it->second, xa_blocks);
        if (ov < min_cov) return false;

        ++ok_ctg;
    }

    return ok_ctg >= min_equiv_contigs_;
}

bool MapqBooster::should_boost_from_group_(const bam1_t* cand, const std::vector<bam1_t*>& group, const sam_hdr_t* hdr) const
{
    if (!cand || !hdr)                return false;
    if (!is_mapped(cand))             return false;
    if (cand->core.qual > mapq_low_)  return false;

    const uint8_t* xa = bam_aux_get(cand, "XA");
    if (xa && bam_aux_type(xa) == 'Z') {
        return should_boost_with_XA_(cand, hdr);
    }

    const char* ref_c_c = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), cand->core.tid);
    if (!ref_c_c) return false;
    const std::string_view ref_c(ref_c_c);

    const uint32_t beg_c = cand->core.pos;
    const auto cand_blocks = CIGAR::ref_blocks(beg_c, bam_cigar_data(cand), bam_cigar_size(cand));
    const uint32_t span = CIGAR::total_block_bases(cand_blocks);
    if (!span) return false;

    std::unordered_map<std::string,std::vector<std::pair<uint32_t,uint32_t>>> eq;
    eq.reserve(16);
    lift_blocks_to_eq(map_, ref_c, cand_blocks, cm_max_hops_, cm_max_fanout_, cm_min_len_, cm_min_frac_, cm_max_total_hits_, eq);
    merge_eq_inplace(eq);

    const uint32_t min_cov = (uint32_t)std::ceil(min_frac_ * span);
    const int end_id_c = read_end_id(cand);

    std::unordered_set<std::string_view> supported;
    supported.reserve(16);

    for (const bam1_t* o : group)
    {
        if (!o || o == cand) continue;
        if (!is_mapped(o)) continue;
        if (read_end_id(o) != end_id_c) continue;

        const char* r2_c = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), o->core.tid);
        if (!r2_c) continue;
        const std::string_view r2(r2_c);

        const uint32_t obeg = o->core.pos;
        const auto other_blocks = CIGAR::ref_blocks(obeg, bam_cigar_data(o), bam_cigar_size(o));
        if (other_blocks.empty()) continue;

        auto it = eq.find(std::string(r2));
        if (it == eq.end()) continue;

        const uint32_t ov = overlapped_len_blocks(it->second, other_blocks);
        if (ov < min_cov) continue;

        supported.insert(r2);
    }

    return (int)supported.size() >= min_equiv_contigs_;
}

void MapqBooster::process_group_(std::vector<bam1_t*>& group, const sam_hdr_t* hdr, Stats& st_local) const {
    if (!hdr) return;

    // 1. split into subgroups by (same end_id) and (query-overlap)
    std::vector<std::vector<bam1_t*>> subgroups;  // First record in each subgroup is the "best" record by match length, MAPQ, and NM, and closest to primary if possible
    subgroups.reserve(8);
    build_subgroups_by_query_overlap(group, subgroups, hdr);

    // 2. Boost for each subgroup
    for (auto& sub : subgroups) {
        if (sub.empty()) continue;

        bam1_t* b = sub[0];  // best record in subgroup
        if (b->core.qual > mapq_low_) continue;  // Record which need to be boost.

        if (should_boost_from_group_(b, sub, hdr)) {
            if (b->core.qual != mapq_new_) {
                b->core.qual = mapq_new_;
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
    auto tracker = ProgressTracker::Every(3000000);

    auto flush_batch = [&]() -> bool {
        if (batch.empty()) return true;

        std::atomic<size_t> idx{0};
        const int T = std::max(1, threads_);

        // thread-local stats (avoid atomic hot path)
        std::vector<Stats> tls((size_t)T);

        auto worker = [&](int tid) {
            size_t i;
            while ((i = idx.fetch_add(1)) < batch.size()) {
                process_group_(batch[i].recs, hdr.get(), tls[(size_t)tid]);
            }
        };

        std::vector<std::thread> pool;
        pool.reserve((size_t)T);
        for (int t = 0; t < T; ++t) pool.emplace_back(worker, t);
        for (auto& th : pool) th.join();

        for (const auto& x : tls) st.changed += x.changed;

        // Write in original order
        for (auto& g : batch) {
            for (auto* b : g.recs) {
                if (sam_write1(out.get(), hdr.get(), b) < 0) { error_stream() << "Write error: " << g.qname << "\n"; std::exit(3); }
                ++st.written;
                bam_destroy1(b);
            }
            g.recs.clear();
        }
        batch.clear();
        batch_records = 0;
        return true;
    };

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

    log_stream() << "processed=" << st.written << " changed=" << st.changed << "\n";
    return 0;
}

void MapqBooster::print_group(
    const std::vector<bam1_t*>& group,
    const std::vector<std::vector<bam1_t*>>& subs,
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
    constexpr int W_RPOS  = 13;

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
           << std::setw(W_RPOS)  << "r_pos"
           << "cigar"
           << "\n";
    };

    auto print_one = [&](const bam1_t* b, int idx, int subidx){
        if (!b) return;

        const bool mapped = is_mapped(b);
        const bool pri_like = is_primary_like(b);

        const char* rname = mapped
            ? sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), b->core.tid)
            : "*";
        const int32_t pos0 = mapped ? b->core.pos : -1;

        uint32_t qb = 0, qe = 0;
        CIGAR::query_interval_fwd(bam_cigar_data(b), bam_cigar_size(b), bam_is_rev(b), qb, qe);

        const uint8_t mq = b->core.qual;

        const std::string cig = CIGAR::to_string(bam_cigar_data(b), bam_cigar_size(b));
        int32_t as = INT32_MIN;
        if (const uint8_t* as_tmp = bam_aux_get(b, "AS")) as = (int32_t)bam_aux2i(as_tmp);
        const int32_t nm = get_nm(b);
        const uint32_t ml = CIGAR::match_len(bam_cigar_data(b), bam_cigar_size(b));

        if (idx == 0 && subidx == 0) print_header();

        std::string idxss_str = "(" + std::to_string(idx) + "," + std::to_string(subidx) + ")";
        std::string qiss_str = "[" + std::to_string(qb) + "," + std::to_string(qe) + ")";

        debug_stream() << "  " << std::left
           << std::setw(W_IDX)   << idxss_str
           << std::setw(W_QINT)  << qiss_str
           << std::setw(W_FLAG)  << b->core.flag
           << std::setw(W_MAPQ)  << (int)mq
           << std::setw(W_PRI)   << (pri_like ? "T" : "F")
           << std::setw(W_AS)    << as
           << std::setw(W_ML)    << ml
           << std::setw(W_NM)    << nm
           << std::setw(W_RNAME) << (rname ? rname : "*")
           << std::setw(W_RPOS)  << pos0
           << cig
           << "\n";
    };

    const char* qn0 = bam_get_qname(group[0]);
    debug_stream() << "QNAME: " << (qn0 ? qn0 : "") << "  (records=" << group.size() << ")" << "  (subgroups=" << subs.size() << ")" << "\n";

    for (size_t si = 0; si < subs.size(); ++si) {
        const auto& sub = subs[si];
        for (size_t k = 0; k < sub.size(); ++k) {
            print_one(sub[k], (int)si, (int)k);
        }
    }

    debug_stream() << "\n";
}

} // namespace mapqboost