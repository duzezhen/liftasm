#pragma once
#include <unordered_set>
#include <iterator>
#include <utility>
#include <optional>

#include "gfa_parser.hpp"
#include "seg_replace.hpp"
#include "logger.hpp"
#include "seq_utils.hpp"
#include "options.hpp"
#include "gfa_name.hpp"
#include "aligner.hpp"

// GfaDeoverlapper
// ------------------------------------------------------------------
// Purpose:
//   Transform an "overlap-based" GFA into a "non-overlap" GFA by
//   cutting segments at consistent boundaries and rebuilding edges.
//
// Workflow:
//   1. Prune inconsistent overlaps (forward vs reverse mismatches).
//   2. Initialize cut points for all segments (0, L, and overlap windows).
//   3. Iteratively propagate cut points between connected segments
//      until convergence (or max iterations).
//   4. Build canonical replacement rules (SegReplace).
//   5. Expand and rewire edges based on rules to eliminate overlaps.
//   6. Reconstruct the GFA with non-overlapping segments.
//
// Usage:
//   GfaDeoverlapper deoverlapper;
//   deoverlapper.load_from_GFA(FILE);
//   deoverlapper.deoverlap();
//   deoverlapper.write_to_GFA(FILE);

class GfaDeoverlapExpandRewirer;

class GfaDeoverlapper : public GfaGraph {
public:
    friend class GfaDeoverlapExpandRewirer;

    struct Cuts {
        std::vector<uint32_t> v;
        std::string           name;
    };

public:

    GfaDeoverlapper(
        int min_eq, 
        double min_match_ratio,
        double min_ali_ratio,
        uint8_t min_mapq,
        uint32_t trim_min_len,
        double trim_max_overlap,
        uint16_t max_prop_iters, 
        uint32_t max_abnormal_cut_len, 
        uint32_t min_abnormal_cut_count, 
        uint32_t min_trans_len,
        std::string mm2_preset
    ) {
        set_forbid_overlap(false);
        MIN_EQ_FOR_CUT_ = min_eq;
        MIN_MATCH_RATIO_ = min_match_ratio;
        MIN_ALI_RATIO_ = min_ali_ratio;
        MIN_MAPQ_ = min_mapq;
        TRIM_MIN_ALIGN_LEN_ = trim_min_len;
        TRIM_MAX_OVERLAP_ = trim_max_overlap;
        MAX_PROPAGATION_ITERS_ = max_prop_iters;
        MAX_ABNORMAL_CUT_LEN_ = max_abnormal_cut_len;
        MIN_ABNORMAL_CUT_COUNT_ = min_abnormal_cut_count;
        MIN_TRANS_LEN_ = min_trans_len;
        if (!mm2_preset.empty()) MM2_PRESET_ = mm2_preset;
    }

    void set_opts(
        const opt::ChainOpts&  chain, const opt::AnchorOpts& anchor,
        const opt::ExtendOpts& extend, const opt::AlignOpts&  align, 
        bool use_wfa
    ) {
        chainOpts_ = chain; anchorOpts_ = anchor; extendOpts_ = extend; alignOpts_ = align; use_wfa_ = use_wfa;
    }

    void deoverlap(const std::string& prefix);

protected:  // data will be used in rulemap building
    /* -------------------------------------------- collapse bubbles -------------------------------------------- */
    struct MmWfaHit {
        uint32_t r_beg = 0, r_end = 0;     // alignment start/end on reference (0-based, end-exclusive)
        uint32_t q_beg = 0, q_end = 0;     // alignment start/end on query (0-based, end-exclusive, excluded soft-clip)
        uint8_t mapq = 0;
        std::string cigar;
    };
    struct BubbleAlignment {
        std::string             name_a;   // segment name of A
        std::string             name_b;   // segment name of B
        // Coordinates are on the forward strand of the segments
        std::uint32_t           beg_a;    // alignment start on A (0-based, inclusive)
        std::uint32_t           end_a;    // alignment end on A (0-based, exclusive)
        std::uint32_t           beg_b;    // alignment start on B (0-based, inclusive), excluded soft-clip
        std::uint32_t           end_b;    // alignment end on B (0-based, exclusive), excluded soft-clip
        uint32_t                v_a;      // [31:1]=segment id, [0]=rev, A vertex id
        uint32_t                v_b;      // [31:1]=segment id, [0]=rev, B vertex id
        std::vector<CIGAR::COp> ops;      // alignment cigar (sega vs segb)
        uint8_t                 mapq;     // mapping quality
        int32_t                 idx = INT32_MAX;    // Index of this alignment in bubbles, homologous and forks. Used for deduplication.
        bool                    force_cuts = false; // preserve exact boundaries for externally supplied alignments (augment)

        BubbleAlignment() = default;

        BubbleAlignment(std::string na, std::string nb, std::uint32_t ba, std::uint32_t ea, std::uint32_t bb, std::uint32_t eb, uint32_t va, uint32_t vb, uint8_t mq, std::vector<CIGAR::COp>&& o, int32_t i = UINT32_MAX)
            : name_a(std::move(na)), name_b(std::move(nb)), beg_a(ba), end_a(ea), beg_b(bb), end_b(eb), v_a(va), v_b(vb), mapq(mq), ops(std::move(o)), idx(i) {}
    };

    std::vector<BubbleAlignment> bubble_aligns_;  // Pairwise alignments between bubble branches

    std::vector<Cuts> cuts_;  // Cutting points for each segment
    SegReplace::RuleMap rulemap_;  // Replacement rules for segments
    std::vector<uint8_t> keep_unused_nodes_;  // Whether to keep segments that are not involved in any alignment-based replacement rule.

    int MIN_EQ_FOR_CUT_;              // Only when the length of '=' in CIGAR >= threshold, a cut point is made / rulemap is established
    double MIN_MATCH_RATIO_;          // Only when the ratio of '=/M' in CIGAR >= threshold, the alignment will be saved for rulemap building
    double MIN_ALI_RATIO_;            // minimum aligned length ratio against shorter sequence
    uint8_t MIN_MAPQ_;                // Minimum mapping quality to be considered
    uint32_t TRIM_MIN_ALIGN_LEN_ = 5'000'000;
    double TRIM_MAX_OVERLAP_ = 0.05;
    uint16_t MAX_PROPAGATION_ITERS_;  // Maximum number of propagation iterations

    uint32_t MAX_ABNORMAL_CUT_LEN_ = 4;     // Maximum short segment length (bp) considered abnormal
    uint32_t MIN_ABNORMAL_CUT_COUNT_ = 5;   // Minimum consecutive abnormal short segments required for pruning

    uint32_t MIN_TRANS_LEN_ = 3;  // Minimum interval length to allow transitive replacement expansion (i.e. if A->B and B->C, then A->C is allowed only if the replacement interval is >= this length)

    std::string MM2_PRESET_ = "asm5";

    opt::ChainOpts  chainOpts_;
    opt::AnchorOpts anchorOpts_;
    opt::ExtendOpts extendOpts_;
    opt::AlignOpts  alignOpts_;
    bool            use_wfa_ = false;  // whether to use WFA for alignment (default: use minimap2)


// Used to prevent generating cycle in graph.
protected:
    struct PairRuleWin_ {
        uint32_t a_sid = UINT32_MAX;
        uint32_t a_beg = 0, a_end = 0;
        bool     a_rev = false;

        uint32_t b_sid = UINT32_MAX;
        uint32_t b_beg = 0, b_end = 0;
        bool     b_rev = false;

        int64_t  delta = 0;
    };

    using PairRuleIndex_ = std::unordered_map<uint64_t, std::vector<PairRuleWin_>>;


protected:
    /**
     * @brief Decide which segment should be treated as the leaf between two segments.
     *
     * If @p by_name is false, the decision is based on segment IDs; otherwise, it is
     * based on lexicographical order of segment names. If both segment IDs are equal,
     * no decision is made.
     *
     * @param v_seg_id   Segment ID of v.
     * @param w_seg_id   Segment ID of w.
     * @param v_is_leaf  Output flag; set to true if v is chosen as the leaf.
     * @param by_name    Whether to compare by segment name instead of ID.
     *
     * @return True if a decision was made, false if v and w refer to the same segment.
     */
    bool prefer_v_as_leaf_(uint32_t v_seg_id, uint32_t w_seg_id, bool& v_is_leaf, bool by_name=false) const;

    /**
     * @brief Merge new cut positions into an existing sorted cut list with deduplication.
     *
     * The input @p new_cuts may be unsorted or reversed; it will be normalized (sorted ascending
     * and deduplicated) before merging. If all elements already exist in @p cuts, no modification
     * is performed.
     *
     * @param cuts      Existing sorted and deduplicated cut vector (will be updated in-place).
     * @param new_cuts  Candidate cut positions to merge.
     *
     * @return True if @p cuts was modified (i.e., at least one new cut was added), false otherwise.
     */
    bool merge_and_dedup_cuts_(std::vector<uint32_t>& cuts, const std::vector<uint32_t>& new_cuts);

    /**
     * @brief Build hash sets of cut positions for all segments.
     *
     * The returned sets are used for fast cut existence lookup during rule generation.
     *
     * @return One unordered_set per segment, containing all cut positions on segments.
     */
    std::vector<std::unordered_set<uint32_t>> build_cut_sets_() const;

    /**
     * @brief Merge candidate cuts into an existing cut list using a hash set for deduplication.
     *
     * Only cuts not already present in @p cut_set are appended to the merge input. If no truly new
     * cut exists, the function returns false.
     *
     * @param cuts            Sorted/deduplicated cut vector to update.
     * @param cut_set         Hash set tracking existing cuts.
     * @param candidate_cuts  Candidate cut positions to insert.
     *
     * @return True if new cuts were added to @p cuts, false otherwise.
     */
    bool merge_and_dedup_cuts_with_set_(
        std::vector<uint32_t>& cuts,
        std::unordered_set<uint32_t>& cut_set,
        const std::vector<uint32_t>& candidate_cuts
    );

protected:
    // Generate and filter the "forward and reverse inconsistent overlap" mapping, and mark the edges that need to be removed in arcs_
    void prune_overlaps_();

    // Initialize all segment cut points (0, L, and "window endpoints"), and sort and deduplicate
    void initialize_cuts_();

    // Align v-slice (reference) against w-slice (query) and pack as BubbleAlignment.
    bool trim_small_overlap_(MmWfaHit& hit, const MmWfaHit& kept) const;
    std::vector<MmWfaHit> filter_aligns_(std::vector<MmWfaHit> a, uint32_t ref_len, uint32_t qry_len) const;

    std::vector<MmWfaHit> align_short_wfa_(
        const std::string& v_name,
        const std::string& w_name,
        const std::string& v_seq_slice,
        const std::string& w_seq_slice
    );

    std::vector<MmWfaHit> align_wfa_(
        const std::string& v_name,
        const std::string& w_name,
        const std::string& v_seq_slice,
        const std::string& w_seq_slice
    );
    std::vector<MmWfaHit> align_short_mm2_(
        const std::string& v_name,
        const std::string& w_name,
        const std::string& v_seq_slice,
        const std::string& w_seq_slice
    );

    std::vector<MmWfaHit> align_mm2_(
        const std::string& v_name,
        const std::string& w_name,
        const std::string& v_seq_slice,
        const std::string& w_seq_slice
    );
    std::vector<BubbleAlignment> align_and_pack_(
        uint32_t            v_vertex,           // packed vertex id for A
        uint32_t            w_vertex,           // packed vertex id for B
        const std::string&  v_name,             // segment name of A
        const std::string&  w_name,             // segment name of B
        const std::string&  v_seq_slice,        // oriented slice sequence of A (already rev-comp if needed)
        const std::string&  w_seq_slice,        // oriented slice sequence of B (already rev-comp if needed)
        uint32_t            vb, uint32_t ve,    // slice coordinates on A in plus strand: [vb, ve)
        uint32_t            wb, uint32_t we     // slice coordinates on B in plus strand: [wb, we)
    );
    // Align overlaps between edges
    void overlaps_align_();

    // Deduplicate alignments (keep the longest among identical segment pairs)
    void dedup_aligns_(size_t begin = 0);

    // Build alignment groups by connected segments. (2026-05-25, v0.1.3-r7)
    std::vector<std::vector<size_t>> build_align_groups_() const;

    // Build, propagate and prune cuts based on alignment results. (2026-05-25, v0.1.3-r7)
    void normalize_all_cuts_();
    bool add_cuts_from_one_alignment_(const BubbleAlignment& align);
    void normalize_group_cuts_(const std::vector<size_t>& group);
    uint64_t propagate_cuts_run_(const std::vector<size_t>& group);
    uint64_t prune_cuts_run_(const std::vector<size_t>& group);
    std::pair<uint64_t, uint64_t> build_propagate_prune_cuts_run_(const std::vector<size_t>& group);
    void build_propagate_prune_cuts_(const std::vector<std::vector<size_t>>& groups);

    // Build replacement rules rulemap
    void build_rulemap_(const std::vector<std::vector<size_t>>& groups);
    SegReplace::RuleMap build_rulemap_run_(
        const std::vector<size_t>& align_ids,
        std::vector<std::unordered_set<uint32_t>>& cut_sets
    );
    SegReplace::Expander build_SegReplace_();

    // Print/verify expansion index
    static void rulemap_verify(const std::vector<GfaNode>& nodes, const SegReplace::RuleMap& idx);

    /**
     * @brief Expand segment intervals around each arc and rewire the graph with new edges.
     * 
     * Example:
     *
     *   Suppose overlap and de-overlap graph:
     *
     *       v ----overlap----> w
     *
     *   and the cuts split the segments into:
     *
     *       v: [left non-overlap]  [overlap pieces]
     *       w:                     [overlap pieces]  [right non-overlap]
     *
     *   Expands the corresponding interval pieces and produces a single rewritten chain like:
     *
     *       [v_left_1] | [v_left_2] | [v_ov_1] | [v_ov_2] | [w_right_1] | [w_right_2]
     *
     * @param ex  the expander used to query replacement expansions for segment intervals
     * 
     */
    void expand_and_rewire_edges_(const SegReplace::Expander& ex);

    // Remove unused nodes
    void remove_unused_nodes_();

    // Print cut points (use after initialize_cuts_ and propagate_cuts_)
    static void print_cuts(const std::vector<Cuts>& cuts);

    void build_rules_from_pair_windows_(
        SegReplace::RuleMap& rulemap,
        // trunk
        uint32_t trunk_seg_id, bool trunk_is_rev,
        const std::vector<uint32_t>& cut_in_overlap_trunk,
        // leaf
        uint32_t leaf_seg_id, bool leaf_is_rev,
        const std::vector<uint32_t>& new_cuts_leaf,
        const std::vector<int32_t>&  new_offsets_leaf,
        // leaf older cuts (in overlap region only)
        const std::vector<uint32_t>& cut_in_overlap_leaf, 
        PairRuleIndex_& index
    ) const;

    /**
     * @brief Convert an overlap window and the cuts at both ends into replacement rules and write them into the rulemap.
     * 
     * The overlap of v/w is [vb,ve), [wb,we) (half-open interval), and cuts are already sorted and deduplicated (the function will still sort/unique as needed internally).
     * 
     * @param v_seg_id    segment ID of v
     * @param w_seg_id    segment ID of w
     * @param v_is_rev    orientation of v
     * @param w_is_rev    orientation of w
     * @param v_name      segment name of v (for logging only)
     * @param w_name      segment name of w (for logging only)
     * @param vb, ve      overlap window on v in forward strand coordinates [vb, ve)
     * @param wb, we      overlap window on w in forward strand coordinates [wb, we)
     * @param rulemap     rule table to write the generated rules into (output)
     * @param cut_sets    pre-built cut sets for quick lookup during rule generation
     * 
     */
    void record_rulemap(
        uint32_t v_seg_id, uint32_t w_seg_id,
        bool v_is_rev, bool w_is_rev,
        const std::string& v_name, const std::string& w_name,
        uint32_t vb, uint32_t ve, uint32_t wb, uint32_t we,
        SegReplace::RuleMap& rulemap, 
        std::vector<std::unordered_set<uint32_t>>& cut_sets, 
        PairRuleIndex_& index
    );

    /**
     * @brief Refine an expansion chain by merging adjacent intervals when possible.
     *
     * Example:
     *
     *   Input chain:
     *
     *       [A:0-10] | [A:10-20] | [A:20-30] | [B:0-10]
     *
     *   Suppose:
     *       - The first three intervals are touchable and form [A:0-30]
     *       - ex.query(A:0-30) returns a non-identity expansion, e.g.: [X:0-25]
     *
     *       Before refinement:
     *           [A:0-10] | [A:10-20] | [A:20-30] | [B:0-10]
     *
     *       After refinement:
     *           [X:0-25] | [B:0-10]
     *
     * @param chain   the input expansion chain to be refined (modified in-place)
     * @param ex      the expander to query for expansions of merged intervals
     * 
     */
    void refine_chain_(SegReplace::Expansion& chain, const SegReplace::Expander& ex) const;


// Used to prevent generating cycle in graph.
protected:
    GfaDeoverlapper::PairRuleIndex_ build_pair_rule_index_(const SegReplace::RuleMap& rulemap) const;
    void add_pair_rule_index_(PairRuleIndex_& index, SegReplace::Seg key, const std::vector<SegReplace::Seg>& value) const;

    /**
     * @brief Check if there is a conflicting pair rule in the existing rulemap for the current trunk and leaf segments. A conflicting pair rule is defined as follows:
     *   - If exist_segs and leaf_segs already have an old relation in rulemap, but that relation is different from the current candidate mapping, skip the current mapping to avoid introducing conflicting replacement cycles.
     * @date 2026-05-18
     * @version 0.1.3-r2
     * @note This check is necessary to prevent cycles from being introduced into the graph.
     */
    bool has_conflicting_pair_rule_(
        PairRuleIndex_& index,
        const std::vector<SegReplace::Seg>& exist_segs,
        const std::vector<SegReplace::Seg>& leaf_segs
    ) const;
};


/* ================================================================================================================
 *                                         EXPAND REWIRER START
 * ================================================================================================================ */
class GfaDeoverlapExpandRewirer {
public:
    GfaDeoverlapExpandRewirer(GfaDeoverlapper& graph, const SegReplace::Expander& expander);

    void run();

private:
    struct PieceQuery_ {
        SegReplace::Seg original{};
        SegReplace::Expansion query;
        bool use_original{false};
    };

    struct Chain_ {
        SegReplace::Expansion expansion;
        std::vector<std::vector<uint32_t>> sample_ids;
        std::vector<GfaAux> variant_tags;
    };

    GfaDeoverlapper& graph_;
    const SegReplace::Expander& expander_;

    std::vector<uint32_t> active_degree_;
    std::unordered_map<SegReplace::Seg, PieceQuery_, SegReplace::U128Hash, SegReplace::U128Eq> pieces_;
    std::unordered_set<uint64_t> seen_edges_;

    uint64_t new_node_count_{0};
    uint64_t new_edge_count_{0};
    uint64_t reverted_piece_count_{0};
    gfaName namer_;

private:
    void build_active_degree_();
    void build_piece_queries_();
    void sanitize_local_replacements_();
    void preserve_isolated_identity_nodes_();
    void emit_isolated_expanded_nodes_();
    void rewire_existing_edges_();
    void replace_graph_edges_();

    bool trivial_segment_(uint32_t sid) const;
    bool piece_query_has_repeat_(const SegReplace::Expansion& query) const;
    SegReplace::Expansion safe_expansion_(SegReplace::Seg piece) const;
    void sanitize_arc_repeats_(const GfaArc& arc);
    void sanitize_piece_repeats_(const std::vector<SegReplace::Seg>& local_pieces);

    std::vector<SegReplace::Seg> collect_all_pieces_(uint32_t sid, bool rev) const;
    std::vector<SegReplace::Seg> collect_non_overlap_pieces_(
        uint32_t sid,
        bool rev,
        uint32_t win_beg,
        uint32_t win_end,
        bool left_side
    ) const;
    std::vector<SegReplace::Seg> collect_overlap_pieces_(
        uint32_t sid,
        bool rev,
        const std::vector<uint32_t>& cuts
    ) const;
    std::vector<uint32_t> collect_overlap_cuts_(uint32_t sid, uint32_t beg, uint32_t end) const;

    Chain_ build_chain_(
        const std::vector<SegReplace::Seg>& source,
        const std::vector<uint32_t>& source_sample_ids,
        const GfaAux& source_tags,
        const std::vector<SegReplace::Seg>& overlap,
        const std::vector<uint32_t>& overlap_sample_ids,
        const GfaAux& overlap_tags,
        const std::vector<SegReplace::Seg>& target,
        const std::vector<uint32_t>& target_sample_ids,
        const GfaAux& target_tags
    ) const;
    void append_piece_group_(
        const std::vector<SegReplace::Seg>& pieces,
        const std::vector<uint32_t>& sample_ids,
        const GfaAux& variant_tags,
        Chain_& chain
    ) const;
    bool revert_piece_(SegReplace::Seg piece);
    bool mark_piece_original_(SegReplace::Seg piece);

    std::vector<uint32_t> materialize_chain_(const Chain_& chain);
    uint32_t materialize_segment_(SegReplace::Seg s, const std::vector<uint32_t>& sample_ids, const GfaAux& variant_tags);
    void emit_chain_edges_(const Chain_& chain);
    void emit_bridge_edge_(uint32_t from, uint32_t to);
};
/* ================================================================================================================
 *                                         EXPAND REWIRER END
 * ================================================================================================================ */
