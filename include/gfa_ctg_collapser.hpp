#pragma once

#include "gfa_collapser.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

class CtgCollapseDebugger;

/* ================================================================================================================
 *                                        CTG COLLAPSER START
 * ================================================================================================================ */
class GfaCtgCollapser : public GfaCollapser {
public:
    using GfaCollapser::GfaCollapser;

    void load_ctgs(
        const std::vector<std::string>& hap1_files,
        const std::vector<std::string>& hap2_files,
        const std::vector<std::string>& sample_names,
        uint32_t min_contig_len
    );

    void collapse_ctgs(
        const std::vector<std::string>& vcf_files,
        const std::string& prefix,
        double min_anchor_coverage,
        uint32_t short_contig_len,
        bool anchor_only,
        bool write_paf
    );

    std::vector<std::string> collapse_samples(
        const std::vector<std::string>& hap1_files,
        const std::vector<std::string>& hap2_files,
        const std::vector<std::string>& sample_names,
        const std::vector<std::string>& vcf_files,
        const std::string& prefix,
        double min_anchor_coverage,
        uint32_t short_contig_len,
        uint32_t min_contig_len,
        bool anchor_only,
        bool write_paf,
        const std::string& command_line
    ) const;

    void set_component_filter(double min_coverage, double end_fraction) noexcept;
    void rebuild_component_paths();

private:
    friend class CtgCollapseDebugger;

    double min_component_coverage_{0.5};
    double component_end_fraction_{0.01};

    struct SegmentOrigin {
        uint32_t sample{UINT32_MAX};
        uint8_t hap{0};
    };

    struct Anchor {
        uint32_t a_sid{0}, b_sid{0};
        uint32_t a_beg{0}, a_end{0};
        uint32_t b_beg{0}, b_end{0};
        bool reverse{false};
        const std::string* read_name{nullptr};
    };

    struct AnchorBlock {
        uint32_t a_sid{0}, b_sid{0};
        uint32_t a_beg{0}, a_end{0};
        uint32_t b_beg{0}, b_end{0};
        uint32_t anchors{0};
        bool reverse{false};
    };

    struct Backbone {
        std::string name;
        uint32_t sid{0};
    };

    struct ReadOccurrence {
        const GfaAlignment* hap1{nullptr};
        const GfaAlignment* hap2{nullptr};
        uint32_t hap1_count{0};
        uint32_t hap2_count{0};
    };

    struct PathPiece {
        uint32_t vertex{0};
        uint32_t beg{0}, end{0};
    };

    struct ComponentBackbone {
        std::string name;
        std::string sequence;
        std::vector<PathPiece> pieces;
        uint32_t id{0};
        uint32_t sample{0};
    };

    struct ComponentHit {
        uint32_t rid{0};
        uint32_t ref_beg{0}, ref_end{0};
        uint32_t query_beg{0}, query_end{0};
        uint32_t raw_query_beg{0}, raw_query_end{0};
        uint32_t matches{0}, block_len{0};
        uint8_t mapq{0};
        bool reverse{false};
        std::string cigar;
        std::vector<CIGAR::COp> ops;
    };

    std::vector<SegmentOrigin> origins_;
    std::vector<Backbone> backbones_;
    std::vector<std::string> sample_names_;

    static std::string file_prefix_(const std::string& path);
    static uint32_t alignment_span_(const GfaAlignment& alignment, uint32_t segment_len);
    static bool intervals_overlap_(uint32_t a_beg, uint32_t a_end, uint32_t b_beg, uint32_t b_end);

    void index_origins_(const std::vector<std::string>& sample_names);
    void filter_short_contigs_(uint32_t min_contig_len);
    void namespace_segment_names_(const std::string& suffix);
    void initialize_backbones_();
    std::vector<Anchor> collect_unique_shared_anchors_() const;
    std::vector<AnchorBlock> build_anchor_blocks_(std::vector<Anchor> anchors) const;
    std::vector<AnchorBlock> select_nonconflicting_blocks_(
        std::vector<AnchorBlock> blocks,
        uint32_t short_contig_len
    ) const;
    std::vector<AnchorBlock> filter_spanning_blocks_(
        std::vector<AnchorBlock> blocks,
        double min_coverage
    ) const;
    void align_unanchored_contigs_(
        const std::vector<AnchorBlock>& blocks,
        uint32_t short_contig_len,
        const std::string& prefix,
        bool write_paf
    );
    std::vector<AnchorBlock> build_alignment_regions_(
        const std::vector<AnchorBlock>& anchors,
        bool anchor_only
    ) const;
    std::vector<BubbleAlignment> align_block_(const AnchorBlock& block);
    void align_regions_(const std::vector<AnchorBlock>& regions);
    void rebuild_backbone_paths_(const SegReplace::Expander& expander);
    std::vector<GfaPath> build_component_paths_(
        const std::vector<size_t>& path_ids,
        const std::string& name_prefix,
        size_t* component_count = nullptr,
        uint64_t* removed_cycle_edges = nullptr
    ) const;
    void append_component_paths_();
    void collapse_backbone_samples_(
        const std::vector<std::string>& files,
        const std::vector<std::string>& sample_names,
        const std::string& prefix,
        uint32_t short_contig_len,
        bool write_paf,
        const std::string& command_line
    );
    std::vector<ComponentBackbone> component_backbones_(const std::vector<std::string>& sample_names) const;
    static uint32_t best_reference_sample_(
        const std::vector<ComponentBackbone>& backbones,
        uint64_t& reference_n50
    );
    std::vector<BubbleAlignment> align_component_backbones_(
        const std::vector<ComponentBackbone>& backbones,
        const std::vector<std::string>& sample_names,
        const std::string& prefix,
        uint32_t short_contig_len,
        bool write_paf
    );
    std::vector<std::vector<ComponentHit>> select_component_hits_(
        std::vector<std::vector<ComponentHit>> hits,
        const std::vector<ComponentBackbone>& backbones,
        uint32_t short_contig_len,
        std::vector<uint32_t>& component_groups
    ) const;
    std::vector<BubbleAlignment> split_component_alignment_(
        const ComponentBackbone& ref,
        const ComponentBackbone& query,
        bool query_reverse,
        uint32_t ref_beg,
        uint32_t query_beg,
        uint8_t mapq,
        const std::vector<CIGAR::COp>& ops
    ) const;
};
/* ================================================================================================================
 *                                         CTG COLLAPSER END
 * ================================================================================================================ */

/* ================================================================================================================
 *                                      CTG COLLAPSE DEBUG START
 * ================================================================================================================ */
class CtgCollapseDebugger {
public:
    static void inputs(
        const std::vector<std::string>& hap1_files,
        const std::vector<std::string>& hap2_files,
        const std::vector<std::string>& sample_names
    );
    static void anchor_decision(
        const GfaCtgCollapser& collapser,
        const GfaCtgCollapser::Anchor& anchor,
        const char* decision
    );
    static void blocks(
        const GfaCtgCollapser& collapser,
        const std::vector<GfaCtgCollapser::AnchorBlock>& blocks,
        const char* stage
    );
    static void pair_score(
        const GfaCtgCollapser& collapser,
        uint32_t a_sid,
        uint32_t b_sid,
        uint64_t anchors,
        uint64_t span,
        bool selected
    );
    static void block_decision(
        const GfaCtgCollapser& collapser,
        const GfaCtgCollapser::AnchorBlock& block,
        const char* decision
    );
    static void alignment(
        const GfaCtgCollapser& collapser,
        const GfaCtgCollapser::AnchorBlock& block,
        size_t alignment_count
    );
    static void paths(const GfaCtgCollapser& collapser);
};
/* ================================================================================================================
 *                                       CTG COLLAPSE DEBUG END
 * ================================================================================================================ */
