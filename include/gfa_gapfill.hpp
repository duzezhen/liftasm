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

class GfaGapfill : public GfaGraph {
public:
    struct Params {
        uint64_t min_overlap_bp{1'000'000};
        uint64_t min_contig_bp{1'000'000};
        uint64_t max_gap_bp{10'000'000};
        double min_similarity{0.7};
        double max_target_overlap{0.1};
        double min_probability{0.5};
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

    size_t gapfill(const std::vector<GfaBubble::Bubble>& bubbles);
    void save_samples(const std::string& prefix, const std::string& command_line);

private:
    struct Fragment {
        size_t path_id{0};
        std::string sample;
        uint32_t component{UINT32_MAX};
        uint8_t hap{0};
        uint64_t length{0};
        int64_t layout_start{0};
        bool reverse{false};
        bool eligible{true};
        bool redundant{false};
        int32_t relocation{-1};
        std::vector<uint32_t> vertices;
        std::vector<std::pair<uint32_t, uint32_t>> vertex_index;
        std::vector<uint32_t> rank_drops;
        std::vector<uint32_t> bubble_positions;
        std::vector<uint64_t> path_bp;
        std::vector<uint64_t> bubble_bp;
    };

    struct PathOverlap {
        uint64_t bp{0};
        uint64_t bubble_bp{0};
        uint32_t nodes{0};
        uint32_t a_begin{UINT32_MAX}, a_end{0};
        uint32_t b_begin{UINT32_MAX}, b_end{0};
    };
    using OverlapIndex = std::vector<std::unordered_map<uint32_t, PathOverlap>>;

    struct BoundaryPhaseSupport {
        uint64_t target_bp{0}, bridge_bp{0}, shared_bp{0};
        uint64_t windows{1};
        double similarity{-1.0};
    };

    struct BoundaryAlignment {
        uint64_t reference_cut{0};
        uint64_t query_cut{0};
        uint32_t query_begin{0}, query_end{0};
        uint32_t reference_begin{0}, reference_end{0};
        uint32_t hits{0};
        uint8_t mapq{0};
        double identity{0.0};
        double coverage{0.0};
        bool accepted{false};
        bool reverse{false};
        bool cut_mapped{false};
    };

    struct Candidate {
        enum Status : uint8_t {
            PENDING, KEPT, LOW_PROBABILITY, USED_END,
            PHASE_CONFLICT, COORDINATE_CONFLICT, CYCLE
        };

        uint32_t left{UINT32_MAX}, right{UINT32_MAX}, bridge{UINT32_MAX};
        uint32_t left_pos{0}, bridge_left{0}, bridge_right{0}, right_pos{0};
        uint32_t junction_left_pos{0}, junction_bridge_left{0};
        uint32_t junction_bridge_right{0}, junction_right_pos{0};
        uint64_t left_cut_bp{0}, bridge_left_cut_bp{0};
        uint64_t bridge_right_cut_bp{0}, right_cut_bp{0};
        uint64_t left_shared{0}, right_shared{0};
        uint64_t left_unplaced{0}, right_unplaced{0};
        uint64_t target_distance{0};
        bool target_overlaps{false};
        double left_similarity{0.0}, right_similarity{0.0};
        BoundaryPhaseSupport left_phase, right_phase;
        double left_boundary_similarity{-1.0}, right_boundary_similarity{-1.0};
        double left_alignment_similarity{-1.0}, right_alignment_similarity{-1.0};
        double left_unplaced_similarity{-1.0}, right_unplaced_similarity{-1.0};
        bool left_misassembly{false}, right_misassembly{false};
        bool phase_consistent{false};
        bool homolog_span{false};
        uint32_t sample_support{0};
        uint32_t informative_samples{0};
        uint32_t spanning_samples{0};
        double score{0.0}, phase_score{0.0};
        double probability{0.0};
        Status status{PENDING};
    };

    struct GapRecord {
        std::string sample, contig, left, right, bridge;
        uint8_t hap{0};
        uint32_t component{UINT32_MAX};
        uint64_t gap_beg{0}, gap_end{0};
        double phase_left{0.0}, phase_right{0.0};
        double alignment_left{-1.0}, alignment_right{-1.0};
        double probability{0.0};
    };

    struct RelocationRecord {
        std::string sample, source, side, status, relocated_contig, target_backbone, strand;
        uint8_t hap{0};
        uint32_t source_component{UINT32_MAX}, target_component{UINT32_MAX};
        uint64_t target_beg{0}, target_end{0};
        double probability{0.0};
        bool used_for_gap{false};
        std::vector<uint32_t> unplaced_vertices;
    };

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
        uint32_t hap1_votes{0}, hap2_votes{0};
    };
    using ChainRefs = std::array<std::vector<const Chain*>, 2>;
    struct PrimarySelection {
        ChainRefs primary;
        ChainRefs low_quality;
    };

    struct OverlapSupport {
        uint64_t shared_bp{0};
        uint64_t overlap_bp{0};
        uint64_t bubble_shared_bp{0};
        uint64_t bubble_overlap_bp{0};
        double similarity{0.0};
        double phase_similarity{0.0};
    };

    struct BoundaryWindow {
        uint64_t target_begin{0}, target_end{0};
        uint64_t bridge_begin{0}, bridge_end{0};
    };

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
    std::vector<Fragment> fragments_;
    std::vector<uint8_t> bubble_nodes_;
    std::map<std::string, std::array<std::vector<Chain>, 2>> sample_chains_;
    std::vector<GapRecord> records_;
    std::vector<RelocationRecord> relocations_;

    static bool parse_path_name_(const std::string& name, std::string& sample, uint8_t& hap);
    static uint32_t find_position_(const Fragment& fragment, uint32_t vertex);
    static std::pair<uint32_t, uint64_t> path_region_(const Fragment& fragment, uint32_t beg, uint32_t end);
    static bool candidate_better_(const Candidate& a, const Candidate& b);
    static bool connection_better_(const Candidate& a, const Candidate& b);
    std::string vertex_name_(uint32_t vertex) const;

    bool build_fragment_(size_t path_id, Fragment& fragment);
    void index_bubble_nodes_(const std::vector<GfaBubble::Bubble>& bubbles);
    void build_fragments_();
    void build_graph_order_();
    void rebuild_fragment_index_(Fragment& fragment) const;
    bool ordered_pair_(const Fragment& left, const Fragment& right) const;
    OverlapSupport overlap_support_(const Fragment& a, const Fragment& b, const PathOverlap& overlap) const;
    bool gap_boundary_supported_(
        const Fragment& a, const Fragment& b,
        const PathOverlap& overlap, bool right_boundary
    ) const;
    std::string boundary_sequence_(const Fragment& fragment, uint32_t pos, bool before, uint64_t windows) const;
    double boundary_similarity_(
        const Fragment& target, uint32_t target_pos,
        const Fragment& bridge, uint32_t bridge_pos,
        bool before, uint64_t windows
    ) const;
    std::pair<uint32_t, uint32_t> boundary_phase_region_(
        const Fragment& fragment, uint32_t pos, bool before, uint64_t windows
    ) const;
    std::pair<uint64_t, uint64_t> boundary_alignment_interval_(
        const Fragment& fragment, bool before
    ) const;
    uint64_t boundary_phase_available_(
        const Fragment& fragment, uint32_t pos, bool before
    ) const;
    BoundaryPhaseSupport boundary_phase_similarity_(
        const Fragment& target, uint32_t target_pos,
        const Fragment& bridge, uint32_t bridge_pos,
        bool before
    ) const;
    bool homolog_spans_(uint32_t left, uint32_t right, const OverlapIndex& overlaps) const;
    bool find_anchors_(const Fragment& left, const Fragment& right, const Fragment& bridge, Candidate& candidate) const;
    bool candidate_path_unique_(
        const Fragment& left, const Fragment& right,
        const Fragment& bridge, const Candidate& candidate
    ) const;
    bool expand_boundary_(
        const Fragment& left, const Fragment& right, const Fragment& bridge,
        const std::vector<std::pair<uint32_t, uint32_t>>& left_matches,
        const std::vector<std::pair<uint32_t, uint32_t>>& right_matches,
        Candidate& candidate
    ) const;
    void check_unplaced_sequence_(std::vector<Candidate>& candidates, bool retry = false) const;
    void check_unplaced_candidate_(Candidate& candidate) const;
    std::string fragment_sequence_(const Fragment& fragment, uint32_t begin, uint32_t end) const;
    std::string fragment_subsequence_(
        const Fragment& fragment, uint64_t begin, uint64_t end
    ) const;
    double mm2_similarity_(const std::string& reference, const std::string& query) const;
    bool align_boundary_(
        const std::string& reference,
        const std::string& query,
        bool before,
        BoundaryAlignment& alignment
    ) const;
    bool boundary_alignment_window_(
        const Fragment& target,
        const Fragment& bridge,
        const Candidate& candidate,
        bool before,
        BoundaryWindow& window
    ) const;
    bool shared_boundary_anchor_(
        const Fragment& target,
        const Fragment& bridge,
        const Candidate& candidate,
        uint32_t target_pos,
        bool before,
        uint32_t& bridge_pos
    ) const;
    bool refine_boundary_(Candidate& candidate) const;
    void refine_boundaries_(
        std::vector<Candidate>& selected, std::vector<Candidate>& candidates
    ) const;
    void update_support_(Candidate& candidate, const std::vector<Candidate>& evidence) const;
    bool relocate_misassemblies_(const std::vector<Candidate>& selected, std::unordered_set<uint32_t>& affected_components);
    void mark_used_relocations_(const std::vector<Candidate>& selected);
    void mark_redundant_fragments_(const std::vector<Candidate>& selected);
    std::vector<Candidate> build_candidates_(const std::unordered_set<uint32_t>* components = nullptr) const;
    std::vector<Candidate> select_candidates_(std::vector<Candidate>& candidates);
    Chain build_chain_(uint32_t start, const std::vector<int32_t>& incoming, const std::vector<int32_t>& outgoing, const std::vector<Candidate>& selected) const;
    void build_sample_chains_(const std::vector<Candidate>& selected);
    void print_summary_(size_t selected) const;
    void write_html_(const std::vector<Candidate>& candidates) const;
    uint64_t chain_length_(const Chain& chain) const;
    static uint64_t ng50_(std::vector<uint64_t> lengths, uint64_t component_bp);
    ChainRefs unmatched_small_primary_chains_(
        const std::map<uint32_t, ChainRefs>& components,
        const std::unordered_map<const Chain*, uint64_t>& chain_lengths
    ) const;
    PrimarySelection select_primary_chains_() const;
    void save_haplotype_(const std::string& label, uint8_t hap, const std::vector<const Chain*>& chains, const std::string& prefix, const std::string& command_line, bool primary = false, const char* primary_kind = "pr");
    std::string vertices_sequence_(const std::vector<uint32_t>& vertices) const;
};

class GfaGapfillPlotter {
public:
    struct Anchor {
        uint32_t segment{0}, length{0};
        uint64_t offset{0};
    };

    struct Contig {
        uint32_t id{0}, component{0};
        uint8_t hap{0};
        uint64_t bp{0};
        int64_t start{0};
        bool reverse{false};
        std::string name, sample, first, last;
        std::vector<Anchor> anchors;
    };

    struct Connection {
        uint32_t component{0}, left{0}, right{0}, bridge{0};
        double phase_left{0.0}, phase_right{0.0};
        uint64_t phase_left_target_bp{0}, phase_left_bridge_bp{0};
        uint64_t phase_right_target_bp{0}, phase_right_bridge_bp{0};
        double boundary_left{0.0}, boundary_right{0.0};
        double alignment_left{-1.0}, alignment_right{-1.0};
        double probability{0.0};
        bool homolog_span{false};
        uint32_t sample_support{0};
        uint32_t informative_samples{0};
        uint32_t spanning_samples{0};
        uint64_t left_cut_bp{0}, right_cut_bp{0};
        uint64_t bridge_left_bp{0}, bridge_right_bp{0};
        std::string status;
    };

    static void save(
        const std::string& file,
        const std::vector<Contig>& contigs,
        const std::vector<Connection>& connections
    );
};

class GfaGapfillDebugger {
public:
    struct BoundarySupport {
        double minimizer{-1.0}, phase{-1.0};
        uint64_t first_bubble_bp{0}, second_bubble_bp{0};
    };

    struct Alignment {
        uint64_t target_begin{0}, target_end{0};
        uint64_t bridge_begin{0}, bridge_end{0};
        uint64_t query_cut{0}, reference_cut{0};
        uint32_t query_begin{0}, query_end{0};
        uint32_t reference_begin{0}, reference_end{0};
        uint32_t hits{0};
        uint8_t mapq{0};
        double identity{0.0}, coverage{0.0};
        bool accepted{false};
        bool reverse{false};
    };

    static void path(const std::string& name, const std::string& decision, const std::string& sample = "", uint8_t hap = 0, uint32_t component = UINT32_MAX, size_t nodes = 0, uint64_t bp = 0, const std::string& first = "", const std::string& last = "");
    static void component(const std::string& sample, uint32_t component, size_t paths);
    static void overlap(const std::string& a, const std::string& b, uint64_t shared_bp, uint64_t overlap_bp, double similarity, double a_fraction, double b_fraction, const BoundarySupport& left, const BoundarySupport& right);
    static void search(const std::string& sample, size_t paths, uint64_t tested, uint64_t low_support, uint64_t wrong_group, uint64_t anchor_failed, size_t candidates);
    static void candidate(const std::string& left, const std::string& right, const std::string& bridge, uint32_t component, double left_phase, double right_phase, uint64_t left_target_bp, uint64_t left_bridge_bp, uint64_t right_target_bp, uint64_t right_bridge_bp, double left_boundary, double right_boundary, bool phase_consistent, bool homolog_span, uint32_t sample_support, uint32_t informative_samples, uint32_t spanning_samples, double probability);
    static void boundary_alignment(const std::string& left, const std::string& right, const std::string& bridge, const Alignment& left_alignment, const Alignment& right_alignment);
    static void refined_boundary(
        const std::string& left, uint64_t left_cut,
        const std::string& right, uint64_t right_cut, uint64_t right_length,
        const std::string& bridge, uint64_t bridge_begin, uint64_t bridge_end,
        const BoundarySupport& left_support, const BoundarySupport& right_support,
        double left_alignment, double right_alignment
    );
    static void selection(const std::string& left, const std::string& right, const std::string& bridge, const std::string& decision, double probability);
    static void assignment(const std::string& sample, uint32_t component, uint8_t hap, uint32_t hap1_votes, uint32_t hap2_votes, bool conflict);
    static void primary_candidate(uint32_t component, const std::string& sample, uint64_t component_bp, uint64_t hap1_ng50, uint64_t hap2_ng50);
    static void chain(const std::string& name, uint8_t hap, uint32_t component, size_t contigs, size_t gaps, uint64_t bp);
};
