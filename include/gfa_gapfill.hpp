#pragma once

#include "gfa_parser.hpp"
#include "gfa_bubble_types.hpp"
#include "options.hpp"

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class GfaGapfill;


/* =================================================================================================================
 *                                             GAP BOUNDARY START
 * ================================================================================================================= */
// Builds graph anchors and local evidence, then refines selected cuts with mm2.
class GfaGapfillBoundary {
public:
    struct PhaseSupport {
        uint64_t target_bp{0};   // Bubble bp examined on the target contig.
        uint64_t bridge_bp{0};   // Bubble bp examined on the bridge contig.
        uint64_t shared_bp{0};   // Bubble bp shared by both local paths.
        uint64_t windows{1};     // Phase windows needed for sufficient evidence.
        double similarity{-1.0}; // shared_bp / min(target_bp, bridge_bp); -1 if unavailable.
    };

    struct Alignment {
        uint64_t reference_cut{0}; // Cut relative to the bridge window.
        uint64_t query_cut{0};     // Cut relative to the target window.
        uint32_t query_begin{0}, query_end{0};         // Aligned target interval.
        uint32_t reference_begin{0}, reference_end{0}; // Aligned bridge interval.
        uint32_t hits{0};          // Number of minimap2 hits.
        uint8_t mapq{0};           // Mapping quality of the reported hit.
        double identity{0.0};      // Matching bases / alignment block length.
        double coverage{0.0};      // Aligned target bases / target window length.
        bool accepted{false};      // Identity, coverage and cut mapping passed.
        bool reverse{false};       // Reported hit is reverse-complemented.
        bool cut_mapped{false};    // A match operation defines an exact cut.
    };

    struct Boundary {
        uint32_t target_pos{0};          // Trusted shared-node index on the target.
        uint32_t bridge_pos{0};          // Matching shared-node index on the bridge.
        uint32_t junction_target_pos{0}; // Target junction before skipping its unreliable end.
        uint32_t junction_bridge_pos{0}; // Bridge node paired with the original junction.
        uint64_t target_cut_bp{0};       // Final target cut, initially a node-boundary fallback.
        uint64_t bridge_cut_bp{0};       // Final bridge cut, initially a node-boundary fallback.
        uint64_t target_begin{0}, target_end{0}; // Target minimap2 interval.
        uint64_t bridge_begin{0}, bridge_end{0}; // Bridge minimap2 interval.
        PhaseSupport phase;              // Local bubble-path phase evidence.
        double minimizer{-1.0};          // Local target/bridge minimizer Jaccard.
        double identity{-1.0};           // Boundary alignment identity; -1 if unavailable.
        double coverage{-1.0};           // Boundary alignment coverage; -1 if unavailable.
        bool alignable{false};           // Both cached alignment intervals are non-empty.
    };

    explicit GfaGapfillBoundary(const GfaGapfill& gapfill) : gapfill_(gapfill) {}

    // Cache ordered anchors and all inexpensive evidence used by candidate selection.
    bool prepare(uint32_t left, uint32_t right, uint32_t bridge, Boundary& left_boundary, Boundary& right_boundary) const;
    // Test whether a bridge spans a target pair without running base alignment.
    bool spans(uint32_t left, uint32_t right, uint32_t bridge) const;
    // Inspect one known overlap boundary for debug output.
    Boundary inspect(uint32_t target, uint32_t target_pos, uint32_t bridge, uint32_t bridge_pos, bool before, bool expand_minimizer = false) const;
    // Refine cached node cuts to base-pair cuts when minimap2 passes both filters.
    bool refine(uint32_t left, uint32_t right, uint32_t bridge, Boundary& left_boundary, Boundary& right_boundary) const;

private:
    const GfaGapfill& gapfill_;

    // Anchor discovery and graph-order validation.
    bool find_anchors_(uint32_t left, uint32_t right, uint32_t bridge, Boundary& left_boundary, Boundary& right_boundary) const;
    bool expand_anchors_(uint32_t left, uint32_t right, uint32_t bridge, const std::vector<std::pair<uint32_t, uint32_t>>& left_matches, const std::vector<std::pair<uint32_t, uint32_t>>& right_matches, Boundary& left_boundary, Boundary& right_boundary) const;
    bool path_unique_(uint32_t left, uint32_t right, uint32_t bridge, const Boundary& left_boundary, const Boundary& right_boundary) const;

    // Local phase and minimizer evidence inside trusted boundary sequence.
    uint64_t phase_available_(uint32_t fragment, uint32_t pos, bool before) const;
    std::pair<uint32_t, uint32_t> phase_region_(uint32_t fragment, uint32_t pos, bool before, uint64_t windows) const;
    PhaseSupport phase_similarity_(uint32_t target, uint32_t target_pos, uint32_t bridge, uint32_t bridge_pos, bool before) const;
    double minimizer_similarity_(uint32_t target, uint32_t target_pos, uint32_t bridge, uint32_t bridge_pos, bool before, uint64_t windows) const;

    // Alignment-window construction and minimap2 cut refinement.
    bool shared_anchor_(uint32_t target, uint32_t bridge, const Boundary& left_boundary, const Boundary& right_boundary, uint32_t target_pos, bool before, uint32_t& bridge_pos) const;
    bool prepare_alignment_window_(uint32_t target, uint32_t bridge, const Boundary& left_boundary, const Boundary& right_boundary, bool before, Boundary& boundary) const;
    bool align_(const std::string& reference, const std::string& query, bool before, Alignment& alignment) const;
    uint32_t find_position_(uint32_t fragment, uint32_t vertex) const;
};
/* =================================================================================================================
 *                                              GAP BOUNDARY END
 * ================================================================================================================= */


class GfaGapfill : public GfaGraph {
public:
    struct MisassemblyContig {
        std::string name, sample, source, side;
        uint8_t hap{0};
        uint32_t source_component{UINT32_MAX};
    };

    struct Params {
        uint64_t min_overlap_bp{1'000'000};
        uint64_t min_contig_bp{1'000'000};
        uint64_t max_gap_bp{10'000'000};
        double min_similarity{0.7};
        double max_target_overlap{0.1};
        double min_confidence{0.5};
        uint64_t misassembly_check_bp{10'000'000};
        double misassembly_similarity{0.7};
        double path_difference{0.05};
        uint32_t phase_path_len{10};
        uint64_t phase_skip_bp{500'000};
        uint64_t phase_window_bp{500'000};
        uint64_t phase_min_bp{1'000};
        double dedup_similarity{0.95};
        uint64_t dedup_component_bp{20'000'000};
        uint32_t threads{1};
        std::string html_file;
    };

    explicit GfaGapfill(Params params);
    ~GfaGapfill() = default;

    // Copy minimap2/index settings shared by boundary refinement and relocation checks.
    void set_alignment_options(
        const opt::ChainOpts& chain,
        const opt::AnchorOpts& anchor,
        const opt::ExtendOpts& extend,
        const opt::AlignOpts& align,
        const std::string& preset,
        uint32_t max_occ,
        double min_match,
        double min_ali_ratio,
        uint8_t min_mapq,
        double boundary_identity,
        double boundary_coverage
    );

    // Run candidate discovery, refinement, selection, relocation and chain construction.
    size_t gapfill(const std::vector<GfaBubble::Bubble>& bubbles);
    // Detect suspect terminal sequence and split each interval into one standalone ms segment/path.
    std::vector<MisassemblyContig> prepare_misassemblies(const std::vector<GfaBubble::Bubble>& bubbles);
    // Restore first-pass relocation metadata after the recollapsed graph is loaded.
    void set_misassemblies(
        const std::vector<MisassemblyContig>& records,
        const std::unordered_set<std::string>& relocated
    );
    // Write per-sample, primary and low-quality GFA files plus TSV reports.
    void save_samples(const std::string& prefix, const std::string& command_line);

private:
    friend class GfaGapfillBoundary;

    /* Graph layout records used to order contigs without reference coordinates. */
    struct LayoutAnchor {
        uint32_t segment{0}, length{0}; // graph node
        uint64_t offset{0};             // position on contig
    };

    struct LayoutContig {
        uint32_t id{0}, component{0};  // fragment and component
        uint8_t hap{0};
        uint64_t bp{0};  // contig length
        int64_t start{0};  // layout position
        bool reverse{false};  // need reverse
        std::string name, sample;
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
        std::vector<uint32_t> rank_drops;  // topology rank drops
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
    using BoundaryPhaseSupport = GfaGapfillBoundary::PhaseSupport;
    using BoundaryAlignment = GfaGapfillBoundary::Alignment;

    struct Candidate {
        // Why a candidate was kept or dropped during selection.
        enum Status : uint8_t {
            PENDING,             // not selected yet
            KEPT,                // selected
            LOW_CONFIDENCE,      // score too low
            USED_END,            // left/right end already used
            PHASE_CONFLICT,      // conflicts with another haplotype path
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
        uint32_t sample_support{0}; // samples supporting this bridge path
        uint32_t spanning_samples{0}; // samples that provide a bridge
        double sample_score{0.0}; // sample_support / spanning_samples
        double phase_score{0.5}; // haplotype support
        double alignment_score{1.0}; // boundary support
        double confidence{0.0}; // combined score
        Status status{PENDING}; // selection result
    };

    // Candidate search result produced independently for one sample/thread.
    struct CandidateBatch {
        std::vector<Candidate> candidates; // candidates from one sample
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
        double identity_left{-1.0}, coverage_left{-1.0};
        double identity_right{-1.0}, coverage_right{-1.0};
        uint32_t sample_support{0}, spanning_samples{0};
        double confidence{0.0};
    };

    struct RelocationRecord {
        std::string sample, source, side, status, relocated_contig, target_backbone, strand;
        uint8_t hap{0};
        uint32_t source_component{UINT32_MAX}, target_component{UINT32_MAX};
        uint64_t target_beg{0}, target_end{0};
        double probability{0.0};
        bool used_for_gap{false};
    };

    // Final sample contig assembled from target fragments and bridge intervals.
    struct Chain {
        struct Part {
            uint32_t fragment{UINT32_MAX};
            uint64_t begin{0}, end{0};
            bool filled{false};
        };
        std::string sample;
        uint8_t hap{0};
        uint32_t component{UINT32_MAX};
        std::vector<uint32_t> vertices;
        std::vector<uint32_t> source_fragments;
        std::vector<Part> parts;
        std::vector<Candidate> gaps;
        std::array<uint64_t, 2> hap_bp{0, 0};
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
        uint32_t max_occ{10};
        double min_match{0.9};
        double min_ali_ratio{0.05};
        uint8_t min_mapq{30};
        double boundary_identity{0.9};
        double boundary_coverage{0.8};
        std::string preset{"asm5"};
    } mm2_;
    std::vector<uint8_t> bubble_nodes_;
    std::vector<Fragment> fragments_;
    std::map<std::string, std::array<std::vector<Chain>, 2>> sample_chains_;
    std::vector<GapRecord> records_;
    std::vector<RelocationRecord> relocations_;
    std::vector<MisassemblyContig> misassemblies_;
    std::unordered_map<std::string, uint32_t> misassembly_index_;
    std::unordered_set<std::string> relocated_misassemblies_;
    std::vector<Candidate> prepared_candidates_, prepared_selected_;
    bool candidates_prepared_{false};

    // Pipeline.
    void clear_state_();

    // Contig preparation and graph layout.
    void index_bubble_nodes_(const std::vector<GfaBubble::Bubble>& bubbles);
    void build_fragments_();
    void build_graph_order_();
    static bool parse_path_name_(const std::string& name, std::string& sample, uint8_t& hap);
    void rebuild_fragment_index_(Fragment& fragment) const;
    static void place_contigs_(std::vector<LayoutContig>& contigs);

    // Shared-node candidate evidence.
    std::vector<Candidate> build_candidates_(bool include_long_gaps = false) const;
    CandidateBatch build_sample_candidates_(const std::vector<uint32_t>& sample_fragments, const OverlapIndex& overlaps, const std::vector<std::string>& group_keys, bool include_long_gaps) const;
    bool ordered_pair_(const Fragment& left, const Fragment& right, bool include_long_gaps) const;
    OverlapSupport overlap_support_(const Fragment& a, const Fragment& b, const PathOverlap& overlap) const;
    bool homolog_spans_(uint32_t left, uint32_t right, const OverlapIndex& overlaps) const;
    void update_sample_support_(std::vector<Candidate>& evidence) const;
    static double raw_phase_score_(const Candidate& candidate);

    // Refine left/bridge and bridge/right cut positions with local alignment.
    void refine_boundaries_(std::vector<Candidate>& selected, std::vector<Candidate>& candidates) const;
    static void update_confidence_(Candidate& candidate);

    // Candidate selection.
    std::vector<Candidate> select_candidates_(std::vector<Candidate>& candidates);
    static bool candidate_better_(const Candidate& a, const Candidate& b);
    static bool connection_better_(const Candidate& a, const Candidate& b);

    // Misassembly detection and relocation.
    void check_unplaced_sequence_(std::vector<Candidate>& candidates) const;
    void check_unplaced_candidate_(Candidate& candidate, bool left, bool right) const;
    std::vector<MisassemblyContig> split_misassemblies_(const std::vector<Candidate>& candidates);
    void mark_used_relocations_(const std::vector<Candidate>& selected);

    // Deduplication
    void mark_redundant_fragments_(const std::vector<Candidate>& selected);
    double mm2_similarity_(const std::string& reference, const std::string& query) const;

    // Sequence, chain, and output helpers.
    Chain build_chain_(uint32_t start, const std::vector<int32_t>& incoming, const std::vector<int32_t>& outgoing, const std::vector<Candidate>& selected) const;
    void build_sample_chains_(const std::vector<Candidate>& selected);
    void print_summary_(size_t selected) const;
    void write_html_(const std::vector<Candidate>& candidates) const;
    PrimarySelection select_primary_chains_() const;
    ChainRefs unmatched_small_primary_chains_(const std::map<uint32_t, ChainRefs>& components, const std::unordered_map<const Chain*, uint64_t>& chain_lengths) const;
    void save_haplotype_(const std::string& label, uint8_t hap, const std::vector<const Chain*>& chains, const std::string& prefix, const std::string& command_line, bool primary = false, const char* primary_kind = "pr");
    std::string fragment_sequence_(const Fragment& fragment, uint32_t begin, uint32_t end) const;
    std::string fragment_subsequence_(const Fragment& fragment, uint64_t begin, uint64_t end) const;
    uint64_t chain_length_(const Chain& chain) const;
    static uint64_t ng50_(std::vector<uint64_t> lengths, uint64_t component_bp);
};
