#pragma once
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <ostream>
#include <fstream>
#include <regex>
#include <algorithm>
#include <cctype>
#include <limits>

#include <zlib.h>
#include "logger.hpp"
#include "kio.hpp"
#include "kseq.h"
#include "coordmap.hpp"
#include "kpaf.hpp"

KSEQ_INIT(gzFile, gzread)

namespace liftover {

// Helpers
namespace detail {
    bool parse_u32(std::string_view sv, uint32_t& out);
    bool split3_bed(std::string_view line, std::string_view& chrom, std::string_view& s0, std::string_view& s1, std::string_view& rest);
    std::string hit_label(const coordmap::CoordMap& M, const coordmap::Hit& h);
    void rc_inplace(std::string& s);
} // namespace detail

// ============================ High-level runner ============================
struct RunOpts {
    std::string map_file;  // required
    std::string bed_file;  // required if !do_check
    std::string ref_file;  // required if do_check
    std::string out_file;  // empty => stdout
    
    std::string regex;

    // normal liftover
    double min_frac = 0.9;

    // extend liftover
    uint32_t flank_win  = 5000;   // window size for each extension (bp)
    uint32_t max_flank  = 20000;  // max extension on each side (bp)
    uint32_t max_gap    = 20000;  // max allowed difference between target and query gap
    uint16_t max_hit    = 5;      // max extended liftover results to keep

    // Optional PAF liftover
    std::string paf_file;
    uint32_t min_mapq  = 5;
    uint32_t min_len   = 50000;

    // Check-mode
    bool do_check   = false;
    uint32_t win = 500, step = 500, max_examples = 20;

    // CoordMapOpts
    int      cm_max_hops;
    uint32_t cm_max_fanout;
    uint32_t cm_min_len;
    double   cm_min_frac;
    uint32_t cm_max_total_hits;

    int threads = 1;
};

RunOpts set_opts(
    std::string map_file, std::string bed_file, std::string ref_file, std::string out_file, 
    std::string regex, double min_frac, 
    uint32_t flank_win, uint32_t max_flank, uint32_t max_gap, uint16_t max_hit,
    std::string paf_file, uint32_t min_mapq, uint32_t min_len,
    bool do_check, uint32_t win, uint32_t step, uint32_t max_examples, 
    int cm_max_hops, uint32_t cm_max_fanout, uint32_t cm_min_len, double cm_min_frac, uint32_t cm_max_total_hits, int threads
);


struct BEDinfo { uint32_t beg, end; std::string rest; };
struct Pafinfo {
    std::string tname;
    uint32_t qlen=0, qstart=0, qend=0;
    char strand = '+';
    uint32_t tlen=0, tstart=0, tend=0;
    uint32_t nmatch=0, alen=0, mapq=0;

    // tags
    char tp = 0;          // 'P' or 'I'
    std::string cg;
};
int file_reader(
    const RunOpts& opt, 
    std::unordered_map<std::string, std::vector<BEDinfo>>& bed, 
    coordmap::CoordMap& M, 
    std::unordered_map<std::string, std::vector<Pafinfo>>& paf,
    std::unordered_map<std::string, std::vector<uint32_t>>& paf_end,
    std::unordered_map<std::string, std::vector<uint32_t>>& paf_idx_end
);


struct LIFTresult {
    std::string tname="";
    int64_t tbeg=0, tend=0, qbeg=0, qend=0;
    char strand='+';
    std::string rest="";  // qname_qbeg_qend_suffix

    LIFTresult() = default;
    LIFTresult(std::string tn, int64_t tb, int64_t te, int64_t qb, int64_t qe, char s, std::string r) : tname(std::move(tn)), tbeg(tb), tend(te), qbeg(qb), qend(qe), strand(s), rest(std::move(r)) {}
};
std::vector<std::string> liftover(const RunOpts& opt);

// PAF liftover
std::vector<LIFTresult> paf_liftover(
    const std::unordered_map<std::string, std::vector<Pafinfo>>& pafs,
    const std::unordered_map<std::string, std::vector<uint32_t>>& paf_end,
    const std::unordered_map<std::string, std::vector<uint32_t>>& paf_idx_end,
    std::string qname, uint32_t qbeg, uint32_t qend
);
// CoordMap liftover
std::vector<LIFTresult> map_liftover(
    const RunOpts& opt, 
    const coordmap::CoordMap& M, 
    std::string qname, uint32_t qbeg, uint32_t qend
);

/*
 * @brief Deduplicate liftover results, keeping those with the smallest gap coverage on query.
 * @param: fixed:    always kept results
 * @param: cand:     candidate results to deduplicate against fixed and themselves (must be sorted by score, the first being the best)
 * @param: qbeg:     query begin position, used to compute gap
 * @param: qend:     query end position, used to compute gap
 * 
 * @return de-duplicated list of LIFTresult
*/
std::vector<LIFTresult> dedup_keep_best_gap(
    std::vector<LIFTresult> fixed,
    std::vector<LIFTresult> cand,
    uint32_t qbeg, uint32_t qend
);

void save_liftover_results(const std::string out_file, const std::vector<std::string>& blocks);

// FASTA sequence database for liftover checking
struct SeqDB {
    struct SvHash {
        using is_transparent = void;
        size_t operator()(std::string_view s) const noexcept { return std::hash<std::string_view>{}(s); }
        size_t operator()(const std::string& s) const noexcept { return (*this)(std::string_view{s}); }
    };
    struct SvEq {
        using is_transparent = void;
        bool operator()(std::string_view a, std::string_view b) const noexcept { return a == b; }
        bool operator()(const std::string& a, const std::string& b) const noexcept { return a == b; }
        bool operator()(const std::string& a, std::string_view b) const noexcept { return std::string_view{a} == b; }
        bool operator()(std::string_view a, const std::string& b) const noexcept { return a == std::string_view{b}; }
    };

    std::unordered_map<std::string, std::string, SvHash, SvEq> seq;
    std::vector<std::string> names;

    bool load(const std::string& path) {
        gzFile fp = gzopen(path.c_str(), "rb");
        if (!fp) {
            error_stream() << path << ": No such file or directory" << std::endl;
            std::exit(1);
        }

        kseq_t* ks = kseq_init(fp);
        if (!ks) {
            gzclose(fp);
            error_stream() << "kseq_init failed" << std::endl;
            std::exit(1);
        }

        while (kseq_read(ks) >= 0) {
            if (!ks->name.s || ks->name.l == 0) continue;
            std::string name(ks->name.s, ks->name.l);

            std::string s;
            s.assign(ks->seq.s ? ks->seq.s : "", (size_t)ks->seq.l);
            for (char& c : s) c = (char)std::toupper((unsigned char)c);

            names.push_back(name);
            seq.emplace(std::move(name), std::move(s));
        }

        kseq_destroy(ks);
        gzclose(fp);
        return true;
    }

    const std::string* get(std::string_view name) const {
        auto it = seq.find(name);
        return (it == seq.end()) ? nullptr : &it->second;
    }
};

void check(const coordmap::CoordMap& M, const RunOpts& opt);

} // namespace liftover