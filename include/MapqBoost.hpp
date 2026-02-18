#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

struct bam1_t;
struct sam_hdr_t;
namespace coordmap { class CoordMap; }

namespace mapqboost {

struct Stats {
    std::uint64_t total_in    = 0;
    std::uint64_t written     = 0;
    std::uint64_t changed     = 0;
    std::uint64_t unmapped    = 0;
    std::uint64_t skipped_ref = 0;
};

class MapqBooster {
public:
    MapqBooster(
        const coordmap::CoordMap& map, 
        std::size_t batch_size, std::uint8_t mapq_low, std::uint8_t mapq_new, double min_frac, int min_equiv_contigs, bool name_check, 
        int cm_max_hops, uint32_t cm_max_fanout, uint32_t cm_min_len, double cm_min_frac, uint32_t cm_max_total_hits,
        int threads, int io_threads
    );

    // in_path/out_path: "-" or empty => stdin/stdout
    int run(const std::string& in_path, const std::string& out_path) const;

private:
    const coordmap::CoordMap& map_;
    const std::size_t  batch_size_;
    const std::uint8_t mapq_low_;
    const std::uint8_t mapq_new_;
    const double       min_frac_;
    const int          min_equiv_contigs_;
    const bool         name_check_;

    // CoordMapOpts
    int      cm_max_hops_;
    uint32_t cm_max_fanout_;
    uint32_t cm_min_len_;
    double   cm_min_frac_;
    uint32_t cm_max_total_hits_;

    const int          threads_;
    const int          io_threads_;

private:
    // XA-based decision; returns true = boost, false = undecided/fail
    bool should_boost_with_XA_(const bam1_t* b, const sam_hdr_t* hdr) const;

    // Decide boost for ONE candidate primary-like record using other records in the same QNAME group.
    bool should_boost_from_group_(const bam1_t* cand, const std::vector<bam1_t*>& group, const sam_hdr_t* hdr) const;

    // Boost all candidate primary-like records in one group (may be multiple primaries).
    void process_group_(std::vector<bam1_t*>& group, const sam_hdr_t* hdr, Stats& st_local) const;

    MapqBooster(const MapqBooster&) = delete;
    MapqBooster& operator=(const MapqBooster&) = delete;
    MapqBooster(MapqBooster&&) = delete;
    MapqBooster& operator=(MapqBooster&&) = delete;
};

} // namespace mapqboost