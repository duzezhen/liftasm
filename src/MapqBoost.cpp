#include "../include/MapqBoost.hpp"
#include "../include/coordmap.hpp"
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

namespace mapqboost {

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

// ---------------- Write mode ----------------
static inline std::string detect_write_mode(const std::string& out)
{
    if (out.empty() || out == "-") return "w"; // stdout => SAM
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

// ---------------- CIGAR: reference span (M/D/N/=/X) ----------------
static inline uint32_t ref_span_from_cigar(const bam1_t* b)
{
    const uint32_t* cig = bam_get_cigar(b);
    const uint32_t n = b->core.n_cigar;
    uint32_t span = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const int op = bam_cigar_op(cig[i]);
        const uint32_t len = bam_cigar_oplen(cig[i]);
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF)
            span += len;
    }
    return span;
}

static inline uint32_t ref_span_from_cigar(std::string_view s)
{
    uint32_t span = 0, num = 0;
    for (char c : s) {
        if (std::isdigit((unsigned char)c)) {
            num = num * 10 + (c - '0');
        } else {
            if (c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X') span += num;
            num = 0;
        }
    }
    return span;
}

static inline uint32_t query_len(const bam1_t* b)
{
    const uint32_t* c = bam_get_cigar(b);
    uint32_t n = b->core.n_cigar, len = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t op32 = c[i];
        const uint32_t op   = bam_cigar_op(op32);
        const uint32_t l    = bam_cigar_oplen(op32);
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF)
            len += l;
    }
    return len;
}

static inline uint32_t match_len_cigar(const bam1_t* r)
{
    const uint32_t* c = bam_get_cigar(r);
    uint32_t n = r->core.n_cigar, m = 0;
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t len = bam_cigar_oplen(c[i]);
        switch (bam_cigar_op(c[i])) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF: m += len; break;
            default: break;
        }
    }
    return m;
}

static inline uint32_t match_len_cigar(std::string_view s)
{
    uint32_t m = 0, num = 0;
    for (char c : s) {
        if (std::isdigit((unsigned char)c)) {
            num = num * 10 + (c - '0');
        } else {
            if (c == 'M' || c == '=' || c == 'X') m += num;
            num = 0;
        }
    }
    return m;
}

// ---------------- ctor ----------------
MapqBooster::MapqBooster(
    const coordmap::CoordMap& map, 
    std::size_t batch_size, uint8_t mapq_low, uint8_t mapq_new, double min_frac, int min_equiv_contigs, bool name_check, 
    int cm_max_hops, uint32_t cm_max_fanout, uint32_t cm_min_len, double cm_min_frac, uint32_t cm_max_total_hits, 
    int threads, int io_threads
)
    : map_(map)
    , batch_size_(batch_size >= 1 ? batch_size : 20000), mapq_low_(mapq_low), mapq_new_(mapq_new), min_frac_(min_frac), min_equiv_contigs_(min_equiv_contigs), name_check_(name_check)
    , cm_max_hops_(cm_max_hops), cm_max_fanout_(cm_max_fanout), cm_min_len_(cm_min_len), cm_min_frac_(cm_min_frac), cm_max_total_hits_(cm_max_total_hits)
    , threads_((threads <= 0) ? (int)std::max(1u, std::thread::hardware_concurrency()) : threads) , io_threads_(std::max(1, io_threads))
{}

// ---------------- helpers ----------------
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

static inline void merge_inplace(std::vector<std::pair<uint32_t,uint32_t>>& ivs) {
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

// XA-based boosting
bool MapqBooster::should_boost_with_XA_(const bam1_t* b, const sam_hdr_t* hdr) const
{
    const uint8_t* xa = bam_aux_get(b, "XA");
    if (!xa || bam_aux_type(xa) != 'Z') return false;

    const int32_t nm_pri = get_nm(b);  // NM
    const uint32_t ml_pri = match_len_cigar(b);  // Match length
    const char* r_pri = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), b->core.tid);
    if (!r_pri) return false;

    const uint32_t beg_pri  = b->core.pos;
    const uint32_t span_pri = ref_span_from_cigar(b);
    if (!span_pri) return false;
    const uint32_t end_pri  = beg_pri + span_pri;

    /* build equivalence map for the primary interval */
    std::unordered_map<std::string, std::vector<std::pair<uint32_t,uint32_t>>> eq;
    eq[r_pri].push_back({beg_pri, end_pri});
    auto hits = map_.map_range(r_pri, beg_pri, end_pri, cm_max_hops_, cm_max_fanout_, cm_min_len_, cm_min_frac_, cm_max_total_hits_, true);
    for (auto& h : hits) eq[map_.contig_name(h.ctg)].push_back({h.beg, h.end});
    for (auto& kv : eq) merge_inplace(kv.second);

    const uint32_t min_cov = (uint32_t)std::ceil(min_frac_ * span_pri);

    /* scan XA segments */
    const char* p = bam_aux2Z(xa);
    int ok_ctg = 0;

    while (*p) {
        /* <rname>,<+|-><pos>,<cigar>,<NM>; */

        // rname
        const char* rbeg = p;
        while (*p && *p != ',') ++p;
        if (!*p) break;
        std::string rname(rbeg, p - rbeg);
        ++p;

        // Filter same ref as primary
        if (std::string_view(rname) == std::string_view(r_pri)) return false;

        // strand (+/-)
        if (*p == '+' || *p == '-') ++p;

        // Convert 1-based pos to 0-based
        long pos1 = 0;
        while (std::isdigit((unsigned char)*p)) { pos1 = pos1 * 10 + (*p - '0'); ++p; }
        if (*p != ',') break; ++p;

        // CIGAR
        const char* cbeg = p;
        while (*p && *p != ',') ++p;
        std::string_view cig(cbeg, p - cbeg);
        if (*p != ',') break; ++p;

        // NM
        long nm_xa = 0;
        while (std::isdigit((unsigned char)*p)) { nm_xa = nm_xa * 10 + (*p - '0'); ++p; }
        while (*p == ';') ++p;

        // Filter by NM and match length
        if (nm_pri >= 0 && nm_xa != nm_pri) continue;  // NM
        if (pos1 <= 0) continue;
        if (match_len_cigar(cig) != ml_pri) continue;  // Match length

        const uint32_t span = ref_span_from_cigar(cig);
        if (!span) continue;
        const uint32_t beg = (uint32_t)(pos1 - 1), end = beg + span;

        auto it = eq.find(rname);
        if (it == eq.end()) return false;

        const uint32_t ov = overlapped_len(it->second, beg, end);
        if (ov < min_cov) return false;

        ++ok_ctg;
    }
    return ok_ctg >= min_equiv_contigs_;
}

// Sorted group-based boosting
bool MapqBooster::should_boost_from_group_(const bam1_t* cand, const std::vector<bam1_t*>& group, const sam_hdr_t* hdr) const
{
    // filter
    if (!cand || !hdr)                return false;
    if (!is_mapped(cand))             return false;
    if (!is_primary_like(cand))       return false;
    if (cand->core.qual > mapq_low_)  return false;

    // XA-based check
    const uint8_t* xa = bam_aux_get(cand, "XA");
    if (xa && bam_aux_type(xa) == 'Z') {
        return should_boost_with_XA_(cand, hdr);
    }

    /* candidate reference interval & NM */
    const char* ref_c = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), cand->core.tid);
    if (!ref_c) return false;

    const uint32_t beg_c = cand->core.pos;
    const uint32_t span  = ref_span_from_cigar(cand);
    if (!span) return false;
    const uint32_t end_c = beg_c + span;
    const int32_t  nm_c  = get_nm(cand);  // NM
    const uint32_t ml_c   = match_len_cigar(cand);  // Match length

    /* equivalence intervals on every contig */
    std::unordered_map<std::string,std::vector<std::pair<uint32_t,uint32_t>>> eq;
    eq[ref_c].push_back({beg_c, end_c});
    auto hits = map_.map_range(ref_c, beg_c, end_c, cm_max_hops_, cm_max_fanout_, cm_min_len_, cm_min_frac_, cm_max_total_hits_, true);
    for (auto& h : hits) eq[map_.contig_name(h.ctg)].push_back({h.beg, h.end});
    for (auto& kv : eq) merge_inplace(kv.second);

    const uint32_t min_cov = (uint32_t)std::ceil(min_frac_ * span);

    // query interval of candidate (0-based, closed-open)
    auto query_interval = [](const bam1_t* r, uint32_t& qb, uint32_t& qe){
        qb = 0; qe = 0;
        if (!r) return;
        const uint32_t* cg = bam_get_cigar(r);
        const uint32_t n   = r->core.n_cigar;
        auto aln_q  = [](int op) {return op==BAM_CMATCH || op==BAM_CINS || op==BAM_CEQUAL || op==BAM_CDIFF;};
        auto read_q = [](int op) {return op==BAM_CMATCH || op==BAM_CINS || op==BAM_CEQUAL || op==BAM_CDIFF || op==BAM_CSOFT_CLIP || op==BAM_CHARD_CLIP;};

        uint32_t qpos = 0;
        uint32_t L0   = 0;
        bool started  = false;
        uint32_t qb_st = 0, qe_st = 0;

        for (uint32_t i = 0; i < n; ++i){
            const uint32_t len = bam_cigar_oplen(cg[i]);
            const int op       = bam_cigar_op(cg[i]);
            if (read_q(op)) L0 += len;
            if (aln_q(op)) {
                if (!started) { qb_st = qpos; started = true; }
                qpos += len;
                qe_st = qpos;
            } else {
                if (read_q(op)) qpos += len;
            }
        }
        if (!started) { qb = 0; qe = 0; return; }
        if (bam_is_rev(r)) {
            qb = L0 - qe_st;
            qe = L0 - qb_st;
        } else {
            qb = qb_st;
            qe = qe_st;
        }
    };
    uint32_t qbeg_c=0,qend_c=0; query_interval(cand,qbeg_c,qend_c);

    const int end_id_c = read_end_id(cand);   /* 0/1/2 */

    /* contigs that support homology */
    std::unordered_set<std::string_view> supported; supported.reserve(16);

    for (const bam1_t* o : group)
    {
        if (!o || o==cand) continue;
        if (!is_mapped(o) || (o->core.flag & BAM_FSUPPLEMENTARY)) continue;
        if (read_end_id(o) != end_id_c) continue;  // opposite mate

        // Filtering by NM and match length
        if (nm_c >= 0 && get_nm(o) != nm_c) continue;
        if (match_len_cigar(o) != ml_c) continue;

        // query-overlap
        uint32_t qb=0,qe=0; query_interval(o, qb, qe);
        if (!(qbeg_c < qe && qb < qend_c)) continue;  // no overlap

        // reference homology check
        const char* r2 = sam_hdr_tid2name(const_cast<sam_hdr_t*>(hdr), o->core.tid);
        if (!r2) return false;  // mapped to unknown

        // Filter same ref as candidate
        if (std::string_view(r2) == std::string_view(ref_c)) return false;

        const uint32_t obeg = o->core.pos;
        const uint32_t ospan= ref_span_from_cigar(o);
        if (!ospan) return false;
        const uint32_t oend = obeg + ospan;

        auto it = eq.find(r2);
        if (it == eq.end()) return false;  // ref not equivalent

        if (overlapped_len(it->second, obeg, oend) < min_cov) return false;  // non-homologous

        supported.insert(std::string_view(r2));
    }

    return (int)supported.size() >= min_equiv_contigs_;
}

void MapqBooster::process_group_(std::vector<bam1_t*>& group, const sam_hdr_t* hdr, Stats& st_local) const {
    for (bam1_t* b : group) {
        if (!b) continue;
        if (!is_primary_like(b)) continue;
        if (!is_mapped(b)) continue;
        if (b->core.qual > mapq_low_) continue;

        if (should_boost_from_group_(b, group, hdr)) {
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

} // namespace mapqboost