#pragma once
#include <vector>
#include <string>
#include <cstdint>

#include "gfa_parser.hpp"
#include "logger.hpp"
#include "CIGAR.hpp"
#include "gfa_deoverlapper.hpp"
#include "gfa_bubble.hpp"
#include "path_pair_split.hpp"
#include "seg_replace.hpp"
#include "progress_tracker.hpp"
#include "options.hpp"

#include "bindings/cpp/WFAligner.hpp"

class GfaCollapser : public GfaDeoverlapper {

public:
    GfaCollapser(
        double min_jaccard, 
        int min_eq, double 
        min_match_ratio, 
        int max_prop_iters, 
        uint32_t homo_k, 
        uint32_t homo_w,
        uint32_t repeat_mask_min_len,
        uint32_t repeat_mask_max_period,
        uint32_t repeat_mask_max_mismatch, 
        uint32_t max_abnormal_cut_len,
        uint32_t min_abnormal_cut_count
    )
    : GfaDeoverlapper(min_eq, min_match_ratio, max_prop_iters, max_abnormal_cut_len, min_abnormal_cut_count), mm_opt_{/*k=*/homo_k, /*w=*/homo_w, /*seek_reverse=*/false, /*seed=*/0x8a5cd789635d2dffULL}
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

    /**
     * @brief Collapse homologous sequences within bubbles/paths into a single node.
     * @date 2026-04-03
     * @version 0.1.3
     * 
     * @param bubbles          list of bubbles detected in the graph
     * @param forks            list of fork groups detected in the graph
     * @param homologous_paths list of homologous paths detected in the graph
     * @param prefix           prefix for output files
     * @param bubble_finder    reference to the bubble finder object
     */
    void collapse_homologous_seq(
        const std::vector<GfaBubble::Bubble>& bubbles, 
        const std::vector<GfaBubble::ForkGroup>& forks, 
        const std::vector<GfaBubble::HomologousPath>& homologous_paths, 
        const std::string& prefix, 
        const GfaBubble::GfaBubbleFinder& bubble_finder
    );

protected:

    float MIN_JACCARD_FOR_ALIGN_ = 0.8;      // Threshold: only when Jaccard >= 0.8, perform alignment
    uint32_t REPEAT_MASK_MIN_LEN_ = 24;      // minimum tandem-repeat span to mask
    uint32_t REPEAT_MASK_MAX_PERIOD_ = 12;   // maximum tandem-repeat period to test
    uint32_t REPEAT_MASK_MAX_MISMATCH_ = 2;  // maximum mismatches allowed in a repeat window
    const minimizerdna::Options mm_opt_;     // Used to calculate Jaccard similarity between two homologous paths

private:
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
        const bool skip_same_start,
        std::shared_ptr<path_pair_split::ComparedIndex> shared_cmp,
        std::shared_ptr<std::shared_mutex> shared_cmp_mutex, 
        uint32_t idx = UINT32_MAX
    );

    /**
     * @brief Align forks, bubble and homologous paths to find homologous regions for cutting.
     * @date 2026-04-03
     * @version 0.1.3
     * 
     * @param bubbles          list of bubbles detected in the graph
     * @param forks            list of fork groups detected in the graph
     * @param homologous_paths list of homologous paths detected in the graph
     * @param bubble_finder    reference to the bubble finder object
     */
    void homologous_align_(
        const std::vector<GfaBubble::Bubble>& bubbles, 
        const std::vector<GfaBubble::ForkGroup>& forks, 
        const std::vector<GfaBubble::HomologousPath>& homologous_paths, 
        const GfaBubble::GfaBubbleFinder& bubble_finder
    );
};