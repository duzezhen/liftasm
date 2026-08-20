#pragma once

#include "gfa_parser.hpp"
#include "gfa_bubble_types.hpp"
#include "options.hpp"

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class GfaGapfill;


/* =================================================================================================================
 *                                             GAP BOUNDARY START
 * ================================================================================================================= */
// Finds ordered shared-node anchors and caches local graph evidence.
class GfaGapfillBoundary {
public:
    struct PhaseSupport {
        uint64_t target_bp{0};   // Bubble bp examined on the target contig.
        uint64_t bridge_bp{0};   // Bubble bp examined on the bridge contig.
        uint64_t shared_bp{0};   // Bubble bp shared by both local paths.
        double similarity{-1.0}; // shared_bp / min(target_bp, bridge_bp); -1 if unavailable.
    };

    struct Boundary {
        uint32_t target_pos{0}; // Outermost ordered shared node in the target contig.
        uint32_t bridge_pos{0}; // Same oriented node in the spanning bridge.
        uint64_t target_cut_bp{0}; // Base coordinate where left/right contig is cut.
        uint64_t bridge_cut_bp{0}; // Base coordinate where bridge is cut on this side.
        PhaseSupport phase; // Phase support around this boundary.
    };

    explicit GfaGapfillBoundary(const GfaGapfill& gapfill) : gapfill_(gapfill) {}

    // Cache ordered anchors and all inexpensive evidence used by candidate selection.
    bool prepare(uint32_t left, uint32_t right, uint32_t bridge, Boundary& left_boundary, Boundary& right_boundary) const;
    // Test whether a bridge spans a target pair without running base alignment.
    bool spans(uint32_t left, uint32_t right, uint32_t bridge) const;
    // Inspect one known overlap boundary for debug output.
    Boundary inspect(uint32_t target, uint32_t target_pos, uint32_t bridge, uint32_t bridge_pos, bool before) const;

private:
    const GfaGapfill& gapfill_;

    // Anchor discovery and graph-order validation.
    bool find_anchors_(uint32_t left, uint32_t right, uint32_t bridge, Boundary& left_boundary, Boundary& right_boundary) const;
    bool is_cycle_free_(uint32_t left, uint32_t right, uint32_t bridge, const Boundary& left_boundary, const Boundary& right_boundary) const;

    // Local phase evidence extending inward from the final graph boundary.
    std::pair<uint32_t, uint32_t> phase_region_(uint32_t fragment, uint32_t pos, bool before) const;
    PhaseSupport phase_similarity_(uint32_t target, uint32_t target_pos, uint32_t bridge, uint32_t bridge_pos, bool before) const;

    uint32_t find_position_(uint32_t fragment, uint32_t vertex) const;
};
/* =================================================================================================================
 *                                              GAP BOUNDARY END
 * ================================================================================================================= */


class GfaGapfill : public GfaGraph {
public:
    using PrimaryPaths = std::array<std::vector<GfaBubble::ReferencePath>, 2>;

    struct MisassemblyContig {
        std::string name, sample, source, side;
        uint8_t hap{0};
        uint32_t source_component{UINT32_MAX};
    };

    struct Params {
        // Bubble phase evidence and graph-walk limits.
        uint32_t max_depth{100'000};
        uint64_t dfs_guard{1'000'000};
        double path_difference{0.05};
        std::vector<std::string> full_phase_samples;
        uint32_t phase_path_len{10};
        uint64_t phase_window_bp{500'000};

        // Gap anchors and spanning evidence.
        uint64_t min_contig_bp{1'000'000};
        uint64_t max_gap_bp{10'000'000};
        double max_target_overlap{0.1};
        uint64_t min_overlap_bp{1'000'000};
        double min_similarity{0.7};

        // Misassembly detection and output deduplication.
        uint64_t misassembly_check_bp{10'000'000};
        double misassembly_similarity{0.7};
        double dedup_similarity{0.95};
        uint64_t dedup_component_bp{10'000'000};
        uint32_t threads{1};
        std::string html_file;
    };

    explicit GfaGapfill(Params params);
    ~GfaGapfill() = default;

    // Copy minimap2/index settings used by relocation and deduplication checks.
    void set_alignment_options(
        const opt::ChainOpts& chain,
        const opt::AnchorOpts& anchor,
        const opt::ExtendOpts& extend,
        const opt::AlignOpts& align,
        const std::string& preset,
        double min_match,
        double min_ali_ratio,
        uint8_t min_mapq
    );

    size_t gapfill(const std::vector<GfaBubble::Bubble>& bubbles);
    std::vector<MisassemblyContig> prepare_misassemblies(const std::vector<GfaBubble::Bubble>& bubbles);
    void set_misassemblies(const std::vector<MisassemblyContig>& records, const std::unordered_set<std::string>& relocated);
    PrimaryPaths save_samples(const std::string& prefix, const std::string& command_line);

private:
    friend class GfaGapfillBoundary;

    /* Graph layout records used to order contigs without reference coordinates. */
    struct LayoutAnchor {
        uint32_t segment{0}, length{0}; // graph node
        uint64_t offset{0};             // position on contig
    };

    struct LayoutContig {
        uint32_t id{0}, component{0};  // fragment and component
        uint64_t bp{0};  // contig length
        int64_t start{0};  // layout position
        bool reverse{false};  // need reverse
        std::vector<LayoutAnchor> anchors;
    };

    // Normalized input P path with graph, phase and output state cached once.
    struct Fragment {
        size_t path_id{0};  // original P path id
        std::string sample;  // sample name
        uint32_t component{UINT32_MAX};  // graph component
        uint8_t hap{0};  // haplotype
        uint64_t length{0};  // path length
        int64_t layout_start{0};  // position after layout
        bool reverse{false};  // path was reversed
        bool eligible{true};  // long enough for gapfill
        bool redundant{false};  // covered by a selected gap
        int32_t relocation{-1};  // relocation record id
        std::vector<uint32_t> vertices;  // oriented path nodes
        std::vector<std::pair<uint32_t, uint32_t>> vertex_index;  // node -> path position
        std::vector<uint32_t> bubble_positions;  // bubble node positions
        std::vector<uint64_t> path_bp;  // cumulative path bp
        std::vector<uint64_t> bubble_bp;  // cumulative bubble bp
    };

    struct PathOverlap {
        uint64_t bp{0};  // Total overlapping base between two contig
        uint32_t nodes{0};  // Total overlapping nodes between two contig
        uint32_t a_begin{UINT32_MAX}, a_end{0};  // Overlapping node range on contig A
        uint32_t b_begin{UINT32_MAX}, b_end{0};  // Overlapping node range on contig B
    };
    using OverlapIndex = std::vector<std::unordered_map<uint32_t, PathOverlap>>;

    struct OverlapSupport {
        uint64_t shared_bp{0}; // shared node bp
        uint64_t overlap_bp{0}; // shorter covered region
        double similarity{0.0}; // shared_bp / overlap_bp
    };

    /* Candidate records keep all evidence from graph anchors through final selection. */
    using Boundary = GfaGapfillBoundary::Boundary;

    // Node-level sequence selected between the two Boundary anchors.
    struct FillPath {
        std::vector<uint32_t> vertices; // complete source-to-sink walk
        std::string sequence;           // walk sequence
        uint64_t length{0}, left_cut{0}, right_cut{0}; // full walk and copied interval
    };

    struct Candidate {
        // Why a candidate was kept or dropped during selection.
        enum Status : uint8_t {
            PENDING,             // not selected yet
            KEPT,                // selected
            NO_PATH,             // no bounded source-to-sink graph walk
            PAIR_CONFLICT,       // incomplete full-phase reciprocal pair
            USED_END,            // left/right end already used
            COORDINATE_CONFLICT, // overlaps an already selected connection
            CYCLE                // would create a loop
        };

        uint32_t left{UINT32_MAX}, right{UINT32_MAX}, bridge{UINT32_MAX}; // left -> bridge -> right
        Boundary left_boundary, right_boundary;
        uint64_t left_shared{0}, right_shared{0}; // shared bp on both sides
        uint64_t left_unplaced{0}, right_unplaced{0};
        uint64_t target_distance{0}; // gap or overlap size
        bool target_overlaps{false};  // whether the target contigs overlap in the graph
        double left_unplaced_similarity{-1.0}, right_unplaced_similarity{-1.0};
        bool left_misassembly{false}, right_misassembly{false};
        bool homolog_span{false};
        uint32_t phase_group{UINT32_MAX}; // reciprocal full-phase connections selected as one unit
        uint32_t phase_locus{UINT32_MAX}; // alternative straight/cross units for one overlapping diploid gap
        uint8_t walk_hap{0}; // output haplotype track after full-phase exchanges
        std::shared_ptr<FillPath> fill; // final node-level gap sequence
        double phase_score{0.5}; // haplotype support
        Status status{PENDING}; // selection result
    };

    // Candidate search result produced independently for one sample/thread.
    struct CandidateBatch {
        std::vector<Candidate> candidates; // within-haplotype gap candidates from one sample
        std::vector<Candidate> phase_edges; // cross-haplotype assignment evidence for full-phase samples
        uint64_t tested{0}; // pairs tested
        uint64_t low_support{0}; // not enough shared bp
        uint64_t wrong_group{0}; // wrong sample/component/layout
        uint64_t anchor_failed{0}; // boundary prepare failed
    };

    /* Records retained for final gap and relocation reports. */
    struct GapRecord {
        std::string sample, contig, left, right, bridge;
        uint8_t hap{0};
        uint32_t component{UINT32_MAX};
        uint64_t gap_beg{0}, gap_end{0};
        double phase_left{0.0}, phase_right{0.0};
    };

    struct RelocationRecord {
        std::string sample, source, side, status, relocated_contig, target_backbone, strand;
        uint8_t hap{0};
        uint32_t source_component{UINT32_MAX}, target_component{UINT32_MAX};
        uint64_t target_beg{0}, target_end{0};
        double probability{0.0};
        bool used_for_gap{false};
    };

    // Final sample contig assembled from target fragments and node walks.
    struct Chain {
        // One sequence piece copied from a target fragment or node walk.
        struct Part {
            uint32_t fragment{UINT32_MAX}; // source fragment
            uint32_t gap{UINT32_MAX}; // gap record for a node-walk piece
            uint64_t begin{0}, end{0}; // bp range on the fragment
            bool filled{false}; // copied from a node walk
        };
        std::string sample; // sample name
        uint8_t hap{0}; // assigned haplotype
        uint32_t component{UINT32_MAX}; // graph component
        std::vector<uint32_t> vertices; // assembled graph path
        std::vector<uint32_t> source_fragments; // fragments used by this chain
        std::vector<Part> parts; // bp pieces used for output
        std::vector<Candidate> gaps; // gap connections in this chain
    };
    using ChainRefs = std::array<std::vector<const Chain*>, 2>;
    struct PrimarySelection {
        ChainRefs primary;
        ChainRefs low_quality;
    };

    /* Configuration followed by pipeline state in production order. */
    Params params_;
    struct Mm2Options {
        uint16_t k{25}, w{30};
        uint16_t best_n{5};
        int zdrop{5000};
        double min_match{0.9};
        double min_ali_ratio{0.05};
        uint8_t min_mapq{30};
        std::string preset{"asm5"};
    } mm2_;

    std::vector<uint8_t> bubble_nodes_;
    std::unordered_set<std::string> full_phase_samples_;  // Samples are assembled with Hi-C, Pore-C or trio-binning, and the contigs are fully phased within a single ctg.
    std::unordered_map<std::string, uint32_t> source_ids_;  // source_ids_["HG002.hap1"] = 0; source_ids_["HG002.hap2"] = 1;
    std::unordered_set<uint32_t> source_scope_;
    std::unordered_map<uint64_t, std::vector<uint32_t>> source_transitions_; // oriented P-path edge -> source IDs
    std::vector<Fragment> fragments_;
    std::map<std::string, std::array<std::vector<Chain>, 2>> sample_chains_;  // sample -> haplotype -> chains
    std::vector<GapRecord> records_;
    std::vector<RelocationRecord> relocations_;
    std::vector<MisassemblyContig> misassemblies_;
    std::unordered_map<std::string, uint32_t> misassembly_index_;
    std::unordered_set<std::string> relocated_misassemblies_;
    std::vector<Candidate> prepared_candidates_;
    bool candidates_prepared_{false};

    // ================================================= Clear state =================================================
    void clear_state_();

    // ================================================= Contig preparation and graph layout =================================================
    void index_bubble_nodes_(const std::vector<GfaBubble::Bubble>& bubbles);
    void build_fragments_();
    void build_graph_order_();
    void index_source_transitions_();
    static bool parse_path_name_(const std::string& name, std::string& sample, uint8_t& hap);
    void rebuild_fragment_index_(Fragment& fragment) const;
    static void place_contigs_(std::vector<LayoutContig>& contigs);

    // ================================================= Shared-node candidate evidence =================================================
    std::vector<Candidate> build_candidates_(bool include_long_gaps = false) const;
    CandidateBatch build_sample_candidates_(const std::vector<uint32_t>& sample_fragments, const OverlapIndex& overlaps, bool include_long_gaps) const;
    bool ordered_pair_(const Fragment& left, const Fragment& right, bool include_long_gaps) const;
    bool skips_backbone_(uint32_t left, uint32_t right) const;
    OverlapSupport overlap_support_(const Fragment& a, const Fragment& b, const PathOverlap& overlap) const;
    bool homolog_spans_(uint32_t left, uint32_t right, const OverlapIndex& overlaps) const;
    static double raw_phase_score_(const Candidate& candidate);
    void group_full_phase_candidates_(std::vector<Candidate>& candidates, const std::vector<Candidate>& phase_edges) const;

    // ================================================= Candidate selection =================================================
    std::vector<Candidate> select_and_walk_(std::vector<Candidate>& candidates);
    std::vector<size_t> select_candidates_(std::vector<Candidate>& candidates);
    static bool candidate_better_(const Candidate& a, const Candidate& b);
    static bool connection_better_(const Candidate& a, const Candidate& b);

    // ================================================= Node-level gap filling =================================================
    size_t resolve_fill_paths_(std::vector<Candidate>& candidates, std::vector<size_t>& selected) const;
    size_t resolve_sample_fill_paths_(std::vector<Candidate>& candidates, const std::vector<size_t>& selected, size_t begin, size_t end, const std::vector<uint32_t>& sample_fragments) const;
    bool resolve_fill_path_(Candidate& candidate, const std::vector<uint8_t>& soft_blocked, uint8_t soft_block_mask) const;

    // ================================================= Misassembly detection and relocation =================================================
    void check_unplaced_sequence_(std::vector<Candidate>& candidates) const;
    std::vector<MisassemblyContig> split_misassemblies_(const std::vector<Candidate>& candidates);
    void mark_used_relocations_(const std::vector<Candidate>& selected);

    // ================================================= Assign haplotypes =================================================
    void build_sample_chains_(const std::vector<Candidate>& selected);
    Chain build_chain_(uint32_t start, const std::vector<int32_t>& incoming, const std::vector<int32_t>& outgoing, const std::vector<Candidate>& selected) const;

    // ================================================= Deduplication =================================================
    void mark_redundant_fragments_(const std::vector<Candidate>& selected);
    double mm2_similarity_(const std::string& reference, const std::string& query) const;

    // ================================================= Select the primary chains =================================================
    PrimarySelection select_primary_chains_() const;
    ChainRefs unmatched_small_primary_chains_(const std::map<uint32_t, ChainRefs>& components, const std::unordered_map<const Chain*, uint64_t>& chain_lengths) const;

    // ================================================= Sequence, chain, and output helpers =================================================
    void print_summary_(size_t selected) const;
    void write_html_(const std::vector<Candidate>& candidates) const;
    void save_haplotype_(const std::string& label, uint8_t hap, const std::vector<const Chain*>& chains, const std::string& prefix, const std::string& command_line, bool primary = false, const char* primary_kind = "pr", std::vector<GfaBubble::ReferencePath>* reference_paths = nullptr);
    
    // ================================================= Sequence helpers =================================================
    std::string fragment_sequence_(const Fragment& fragment, uint32_t begin, uint32_t end) const;
    std::string fragment_subsequence_(const Fragment& fragment, uint64_t begin, uint64_t end) const;
    uint64_t chain_length_(const Chain& chain) const;
    static uint64_t ng50_(std::vector<uint64_t> lengths, uint64_t component_bp);
};
