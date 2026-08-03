#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_set>

#include "gfa_parser.hpp"
#include "logger.hpp"
#include "CIGAR.hpp"
#include "gfa_augment.hpp"
#include "gfa_bubble.hpp"
#include "path_pair_split.hpp"
#include "seg_replace.hpp"
#include "progress_tracker.hpp"
#include "options.hpp"

#include "bindings/cpp/WFAligner.hpp"

class LocalGraphComplexity;
class GfaPathCycleGuard;

class GfaCollapser : public GfaAugmenter {

public:
    GfaCollapser(
        double min_jaccard, 
        int min_eq,
        double min_match_ratio,
        double min_ali_ratio,
        uint8_t min_mapq,
        uint32_t trim_min_len,
        double trim_max_overlap,
        int max_prop_iters, 
        uint32_t homo_k, 
        uint32_t homo_w,
        uint32_t all_pair_len,
        uint32_t repeat_mask_min_len,
        uint32_t repeat_mask_max_period,
        uint32_t repeat_mask_max_mismatch, 
        uint32_t max_abnormal_cut_len,
        uint32_t min_abnormal_cut_count, 
        uint32_t min_trans_len,
        std::string mm2_preset
    )
        : GfaAugmenter(
            min_eq,
            min_match_ratio,
            min_ali_ratio,
            min_mapq,
            trim_min_len,
            trim_max_overlap,
            max_prop_iters,
            max_abnormal_cut_len,
            min_abnormal_cut_count,
            min_trans_len,
            mm2_preset
        ),
        mm_opt_{/*k=*/homo_k, /*w=*/homo_w, /*seek_reverse=*/false, /*seed=*/0x8a5cd789635d2dffULL},
        ALL_PAIR_LEN_(all_pair_len)
    {
        MIN_JACCARD_FOR_ALIGN_ = min_jaccard;
        REPEAT_MASK_MIN_LEN_ = repeat_mask_min_len;
        REPEAT_MASK_MAX_PERIOD_ = repeat_mask_max_period;
        REPEAT_MASK_MAX_MISMATCH_ = repeat_mask_max_mismatch;
        // set_forbid_overlap(true);
    }

    /**
     * @brief Merge linear chains in a de-overlap graph.
     * 
     *  Example:
     *
     *   Before: X -> A -> B -> C -> D -> Y
     *    After: X -> (A-B-C-D) -> Y
     *
     * Overlap graphs are not supported. If any active arc has positive overlap, this step is skipped.
    */
    void merge_linear_chains();

    void normalize_homopolymer_bubbles(
        const std::vector<GfaBubble::Bubble>& bubbles,
        uint32_t max_path_len
    );

    void disable_repeat_mask() {
        REPEAT_MASK_MIN_LEN_ = 0;
        REPEAT_MASK_MAX_PERIOD_ = 0;
        REPEAT_MASK_MAX_MISMATCH_ = 0;
    }

    void set_complex_marking(bool enabled) noexcept {
        mark_rejected_complex_ = enabled;
    }

    /**
     * @brief Collapse homologous sequences within bubbles/paths into a single node.
     * @date 2026-04-03
     * @version 0.1.3
     * 
     * @param bubbles          list of bubbles detected in the graph
     * @param homologous_paths list of homologous paths detected in the graph
     * @param prefix           prefix for output files
     */
    void collapse_homologous_seq(
        const std::vector<GfaBubble::Bubble>& bubbles, 
        const std::vector<GfaBubble::HomologousPath>& homologous_paths, 
        const std::string& prefix
    );

protected:

    float MIN_JACCARD_FOR_ALIGN_ = 0.8;      // Threshold: only when Jaccard >= 0.8, perform alignment
    uint32_t REPEAT_MASK_MIN_LEN_ = 24;      // minimum tandem-repeat span to mask
    uint32_t REPEAT_MASK_MAX_PERIOD_ = 12;   // maximum tandem-repeat period to test
    uint32_t REPEAT_MASK_MAX_MISMATCH_ = 2;  // maximum mismatches allowed in a repeat window
    const minimizerdna::Options mm_opt_;     // Used to calculate Jaccard similarity between two homologous paths
    const uint32_t ALL_PAIR_LEN_;

private:
    friend class LocalGraphComplexity;
    friend class GfaPathCycleGuard;

    /* -------------------------------------------- merge_linear_chains -------------------------------------------- */
    /**
     * @brief Extend from a starting vertex to build a linear chain (unitig).
     *
     *   - each internal node has indegree = 1 and outdegree = 1
     *   - extension stops when branching or dead-end is encountered
     *
     * @param start          starting vertex ID
     * @param chain          output vector of vertex IDs forming the chain
     * @param edge_idx_seq   output indices of arcs used in the chain
     *
     * @return true if a valid chain (length >= 2) is found, false otherwise
     */
    bool try_build_chain_(
        uint32_t start,
        std::vector<uint32_t>& chain,
        std::vector<size_t>& edge_idx_seq
    );

    /* -------------------------------------------- collapse bubbles -------------------------------------------- */
    std::vector<BubbleAlignment> align_subpaths_(
        const std::vector<uint32_t>& subA,
        const std::vector<uint32_t>& subB,
        const path_pair_split::Span& spanA,
        const path_pair_split::Span& spanB
    );

    std::vector<BubbleAlignment> align_paths_(
        const std::vector<std::vector<uint32_t>>& paths, 
        const std::vector<uint32_t>& clusters,
        GfaBubble::Type path_type,
        const bool skip_same_start,
        path_pair_split::ComparedIndex& compared,
        uint32_t idx = UINT32_MAX
    );

    struct CandidateTopoSpan_ {
        uint32_t component{UINT32_MAX};
        uint32_t begin{UINT32_MAX};
        uint32_t end{0};
    };

    struct AlignmentCandidate_ {
        bool homologous{false};
        size_t source_index{0};
        GfaBubble::Type type{GfaBubble::Type::Normal};
        bool skip_same_start{false};
        uint32_t idx{UINT32_MAX};
        uint32_t input_paths{0};
        uint64_t priority{0};
        std::vector<CandidateTopoSpan_> spans;
        std::vector<uint32_t> region_vertices;
    };

    std::vector<AlignmentCandidate_> alignment_candidates_;
    bool mark_rejected_complex_{true};

    std::vector<AlignmentCandidate_> build_alignment_candidates_(
        const std::vector<GfaBubble::Bubble>& bubbles,
        const std::vector<GfaBubble::HomologousPath>& homologous_paths
    ) const;

    std::vector<std::vector<size_t>> group_alignment_candidates_(
        const std::vector<AlignmentCandidate_>& candidates
    ) const;

    std::vector<uint32_t> collect_candidate_region_(
        const std::vector<std::vector<uint32_t>>& paths
    ) const;

    std::vector<BubbleAlignment> align_candidate_group_(
        const std::vector<size_t>& group,
        const std::vector<AlignmentCandidate_>& candidates,
        const std::vector<GfaBubble::Bubble>& bubbles,
        const std::vector<GfaBubble::HomologousPath>& homologous_paths
    );

    void reset_collapse_cuts_();

    SegReplace::Expander build_safe_expander_(size_t alignment_begin);

    /**
     * @brief Align forks, bubble and homologous paths to find homologous regions for cutting.
     * @date 2026-04-03
     * @version 0.1.3
     * 
     * @param bubbles          list of bubbles detected in the graph
     * @param homologous_paths list of homologous paths detected in the graph
     * @param bubble_finder    reference to the bubble finder object
     */
    void homologous_align_(
        const std::vector<GfaBubble::Bubble>& bubbles, 
        const std::vector<GfaBubble::HomologousPath>& homologous_paths
    );
};

/* ================================================================================================================
 *                                          LOCAL GRAPH COMPLEXITY START
 * ================================================================================================================ */
class LocalGraphComplexity {
public:
    LocalGraphComplexity(GfaCollapser& collapser, size_t alignment_begin);

    std::unordered_set<int32_t> find_increasing_candidates() const;
    size_t reject_candidates(const std::unordered_set<int32_t>& candidates);

private:
    struct LocalGraph {
        uint32_t node_count{0};
        std::unordered_set<uint64_t> directed_edges;
    };

    struct Metrics {
        uint64_t max_block_cycle_excess{0};
        uint64_t max_block_branch_excess{0};
    };

    // Ignore small cut-induced bubbles; reject only order-of-magnitude growth in one inseparable block.
    static constexpr uint64_t MAX_COMPLEXITY_GROWTH_ = 8;
    GfaCollapser& collapser_;
    const size_t alignment_begin_;

    LocalGraph build_local_graph_(const std::vector<uint32_t>& vertices, const SegReplace::Expander* expander) const;
    bool candidate_increases_(size_t candidate_id, const SegReplace::Expander& expander) const;
    bool increases_(uint32_t input_paths, const LocalGraph& before, const LocalGraph& after) const;
    Metrics measure_(uint32_t node_count, const std::unordered_set<uint64_t>& directed_edges) const;
    static void measure_block_(const std::vector<std::pair<uint32_t, uint32_t>>& edges, const std::vector<uint32_t>& block_edges, Metrics& metrics);
};
/* ================================================================================================================
 *                                          LOCAL GRAPH COMPLEXITY END
 * ================================================================================================================ */

/* ================================================================================================================
 *                                         PATH CYCLE GUARD START
 * ================================================================================================================ */
class GfaPathCycleGuard {
public:
    explicit GfaPathCycleGuard(GfaCollapser& graph);

    size_t filter_path_rules(SegReplace::Expander& expander) const;
    size_t filter_rules(const std::vector<GfaPath>& paths, SegReplace::Expander& expander) const;

private:
    struct ExpandedStep {
        SegReplace::Seg interval{0};
        SegReplace::Seg owner{0};
        bool replaced{false};
    };

    struct PathScan {
        std::vector<ExpandedStep> steps;
        std::vector<SegReplace::Seg> bad_rules;
    };

    struct CycleIndex {
        std::vector<uint32_t> component;
        std::vector<uint8_t> cyclic;
    };

    GfaCollapser& graph_;

    size_t filter_(const std::vector<const GfaPath*>& paths, SegReplace::Expander& expander) const;
    std::unordered_set<SegReplace::Seg, SegReplace::U128Hash, SegReplace::U128Eq> find_bad_rules_(const std::vector<const GfaPath*>& paths, const SegReplace::Expander& expander) const;
    PathScan scan_path_(const GfaPath& path, const SegReplace::Expander& expander) const;
    static CycleIndex find_cycles_(uint32_t node_count, const std::vector<uint64_t>& edges);
};
/* ================================================================================================================
 *                                          PATH CYCLE GUARD END
 * ================================================================================================================ */
