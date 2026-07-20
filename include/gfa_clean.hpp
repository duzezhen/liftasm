#pragma once

#include "gfa_parser.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

struct GfaCleanOptions {
    uint32_t min_reads{5};
    double min_purity{0.8};
    double min_comp_overlap{0.05};
};

class GfaCleaner : public GfaGraph {
public:
    struct ReadPlacement {
        uint32_t component{UINT32_MAX};
        uint32_t contig_id{UINT32_MAX};
        std::string contig_name;
        uint32_t contig_position{0};
        bool is_reverse{false};
        std::string read_name;
        uint32_t read_start{0};
        uint32_t read_end{0};
        std::string optional_tags;
    };

    explicit GfaCleaner(GfaCleanOptions options);

    size_t clean(const std::vector<std::string>& ctg_files);
    void save_read_placements(const std::string& output_file) const;
    const std::vector<ReadPlacement>& read_placements() const { return read_placements_; }

private:
    struct ReadInfo {
        uint32_t id{UINT32_MAX};
        std::vector<uint32_t> components;
    };
    struct LabelCount {
        uint32_t label{0};
        uint32_t count{0};
    };
    struct NodeEvidence {
        std::vector<LabelCount> labels;
        uint32_t total{0};
        uint32_t dominant_label{UINT32_MAX};
        uint32_t dominant_count{0};
    };
    struct ComponentOverlapIndex {
        std::unordered_map<uint32_t, uint32_t> node_counts;
        std::unordered_map<uint64_t, uint32_t> shared_nodes;
    };

    using ReadIndex = std::unordered_map<std::string, ReadInfo>;

    void build_ctg_read_index_(const std::vector<std::string>& ctg_files, ReadIndex& read_index);
    std::vector<NodeEvidence> build_utg_evidence_(const ReadIndex& read_index) const;

    bool resolved_different_(const NodeEvidence& a, const NodeEvidence& b) const;
    ComponentOverlapIndex build_component_overlap_index_(
        const std::vector<NodeEvidence>& evidence,
        const std::unordered_set<uint64_t>& candidate_pairs
    ) const;
    static double component_overlap_(const ComponentOverlapIndex& index, uint32_t a, uint32_t b);
    static uint64_t component_pair_key_(uint32_t a, uint32_t b);

    GfaCleanOptions options_;
    std::vector<ReadPlacement> read_placements_;
};
