#pragma once
#include <vector>
#include <string>
#include <cstdint>

#include "gfa_parser.hpp"
#include "logger.hpp"
#include "CIGAR.hpp"
#include "gfa_deoverlapper.hpp"
#include "gfa_bubble.hpp"
#include "seg_replace.hpp"
#include "progress_tracker.hpp"
#include "options.hpp"

#include "bindings/cpp/WFAligner.hpp"

class GfaCollapser : public GfaDeoverlapper {
public:
    GfaCollapser(double min_jaccard, double min_new_frac, int min_eq, int max_prop_iters)
    : GfaDeoverlapper(min_eq, max_prop_iters)
    {
        MIN_JACCARD_FOR_ALIGN_ = min_jaccard;
        MIN_NEW_LEN_FRAC_ = min_new_frac;
        set_forbid_overlap(true);
    }

    /**
     * @brief Collapse linear unitigs in the GFA graph (no overlaps GFA as input).
     * @date 2025-08-27
     * @author Zezhen Du
     *
     * This function searches for all linear chains (unitigs) in the directed
     * segment graph that satisfy:
     *   - The head (first vertex) has exactly one outgoing non-deletion edge.
     *   - The tail (last vertex) has exactly one incoming non-deletion edge.
     *   - Every internal vertex has indegree = 1 and outdegree = 1.
     *
     * For each such chain:
     *   1. Concatenate the oriented sequences of all nodes to form a new merged segment.
     *   2. Redirect all external edges:
     *        - In-edges of the chain head are rewired to the new segment.
     *        - Out-edges of the chain tail are rewired from the new segment.
     *   3. Delete all internal edges of the chain.
     *   4. Delete all original segments in the chain.
     */
    void collapse_unitigs();

    /**
     * @brief Collapse homologous sequences within bubbles into a single node (no overlaps GFA as input).
     * @date 2025-11-17
     * @author Zezhen Du
     */
    void collapse_bubbles(const std::vector<GfaBubble::Bubble>& bubbles, const std::vector<GfaBubble::ForkGroup>& forks, const std::string& prefix, const GfaBubble::GfaBubbleFinder& bubble_finder);

protected:

    float MIN_JACCARD_FOR_ALIGN_;    // Threshold: only when Jaccard >= 0.8, perform alignment
    double MIN_NEW_LEN_FRAC_;        // --min-new-frac: submit align only if the unaligned fraction ≥ threshold

private:
    /* -------------------------------------------- collapse unitigs -------------------------------------------- */
    // Extend from start forward to get a linear chain; if successful, chain.size() >= 2
    bool try_build_chain_(
        uint32_t start,
        std::vector<uint32_t>& chain,
        std::vector<size_t>& edge_idx_seq
    );

    /* -------------------------------------------- collapse bubbles -------------------------------------------- */
    std::vector<std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> split_paths_(const std::vector<uint32_t> &path1, const std::vector<uint32_t> &path2);
    std::string vertex_name_(uint32_t vtx) const;
    void build_concat_seq_(const std::vector<uint32_t>& path, std::string &seq, std::vector<uint32_t> &offsets);
    std::pair<size_t,uint32_t> find_node_idx_(uint32_t pos, const std::vector<uint32_t> &offs) const;
    std::vector<BubbleAlignment> split_cigar_(
        const std::vector<uint32_t> &pathA,
        const std::vector<uint32_t> &pathB,
        const std::vector<uint32_t> &offsA,
        const std::vector<uint32_t> &offsB,
        const BubbleAlignment &big
    );
    std::vector<BubbleAlignment> align_subpaths_(
        const std::vector<uint32_t> &subA,
        const std::vector<uint32_t> &subB
    );
    /**
     * @brief Align bubble branches pairwise to find homologous regions for cutting.
     * @date 2025-11-15
     * @author Zezhen Du
     */
    void bubbles_align_(const std::vector<GfaBubble::Bubble>& bubbles, const std::vector<GfaBubble::ForkGroup>& forks, const GfaBubble::GfaBubbleFinder& bubble_finder);
    

    void expand_and_rewire_edges_(const SegReplace::Expander& ex);  // Expand and replace edges based on rule index (retain original logs and structure updates), and finally rebuild
};