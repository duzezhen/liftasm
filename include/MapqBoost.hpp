#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

struct bam1_t;
struct sam_hdr_t;
namespace coordmap { class CoordMap; }
namespace rsite { class Index; }

namespace mapqboost {

struct Stats {
    std::uint64_t total_in    = 0;
    std::uint64_t written     = 0;
    std::uint64_t changed     = 0;
    std::uint64_t unmapped    = 0;
    std::uint64_t skipped_ref = 0;
};

struct SubRec {
    bam1_t*  b{nullptr};

    uint32_t qbeg{0}, qend{0};

    int32_t  tid{-1};
    uint32_t rbeg{0}, rend{0};
    uint8_t  mapq{0};
    int32_t  nm{INT32_MAX/2};
    int32_t  as{INT32_MIN};
    uint32_t ml{0};

    bool primary_like{false};

    // restriction site info
    bool rsite_found_L{false}, rsite_found_R{false};
    uint32_t rsite_pos_L{0}, rsite_pos_R{0};
    uint32_t rsite_dist_L{UINT32_MAX}, rsite_dist_R{UINT32_MAX};

    // combined score, calculate from MAPQ, AS, NM, rsite info, etc. Higher is better.
    double score{0.0};
    // Close score to primary
    int    close_score_num{0};
    double homol_frac{0.0};
};

struct Subgroup {
    int end_id{0};
    std::vector<SubRec> recs;   // sorted, best-first
};

class MapqBooster {
public:
    MapqBooster(
        const coordmap::CoordMap& coormap_idx, const rsite::Index& rsite_idx,
        std::size_t batch_size, std::uint8_t mapq_low, std::uint8_t mapq_cap, bool name_check, 
        int cm_max_hops, uint32_t cm_max_fanout, uint32_t cm_min_len, double cm_min_frac, uint32_t cm_max_total_hits,
        double sub_ovlp_frac, 
        double K_mapq, double K_rs, double K_as, double K_ml, double K_nm, double W_mapq, double W_rs, double W_as, double W_ml, double W_nm, double close_as_eps,
        int threads, int io_threads
    );

    // in_path/out_path: "-" or empty => stdin/stdout
    int run(const std::string& in_path, const std::string& out_path) const;

private:
    const coordmap::CoordMap& coormap_idx_;
    const rsite::Index& rsite_idx_;
    const std::size_t   batch_size_;
    const std::uint8_t  mapq_low_;
    const std::uint8_t  mapq_cap_;
    const bool          name_check_;

    // CoordMapOpts
    const int      cm_max_hops_;
    const uint32_t cm_max_fanout_;
    const uint32_t cm_min_len_;
    const double   cm_min_frac_;
    const uint32_t cm_max_total_hits_;

    // Record cluster parameters
    const double sub_ovlp_frac_;
    const double K_mapq_, K_rs_, K_as_, K_ml_, K_nm_;
    const double W_mapq_, W_rs_, W_as_, W_ml_, W_nm_;
    const double close_as_eps_;

    const int threads_;
    const int io_threads_;

private:
    // Build subgroups for one read group:
    // 1. same read_end_id (0/1/2), 2. query intervals overlap (on forward read coordinates)
    void build_subgroups_by_query_overlap(
        const std::vector<bam1_t*>& group,
        std::vector<Subgroup>& out_subgroups, 
        const sam_hdr_t* hdr
    ) const;
    // Re-rank each subgroup and put the chosen alignment at sub[0].
    void select_best_per_subgroup(std::vector<Subgroup>& out_subgroups, const sam_hdr_t* hdr) const;

    // XA-based decision; returns true = boost, false = undecided/fail
    bool should_boost_with_XA_(SubRec& cand, const sam_hdr_t* hdr) const;

    // Decide boost for ONE candidate primary-like record using other records in the same QNAME group.
    bool should_boost_from_group_(Subgroup& sub, const sam_hdr_t* hdr) const;

    // Boost all candidate primary-like records in one group (may be multiple primaries).
    void process_group_(std::vector<bam1_t*>& group, const sam_hdr_t* hdr, Stats& st_local) const;

    MapqBooster(const MapqBooster&) = delete;
    MapqBooster& operator=(const MapqBooster&) = delete;
    MapqBooster(MapqBooster&&) = delete;
    MapqBooster& operator=(MapqBooster&&) = delete;

    // Debug
    void print_group(
        const std::vector<bam1_t*>& group,
        const std::vector<Subgroup>& subs,
        const sam_hdr_t* hdr
    ) const;
};

} // namespace mapqboost