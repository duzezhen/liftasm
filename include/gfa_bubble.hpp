#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <queue>
#include <stack>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <mutex>

#include "logger.hpp"
#include "gfa_parser.hpp"
#include "gfa_bubble_types.hpp"

#include "bloom_filter.hpp"

#include "liftover.hpp"
#include "gfa_name.hpp"
#include "seq_utils.hpp"
#include "save.hpp"


namespace GfaBubble {

class HomologousPathEnumerator;

class GfaBubbleFinder {

    friend class HomologousPathEnumerator;

public:
    GfaBubbleFinder(
        const GfaGraph& graph, uint32_t max_depth, uint16_t max_paths, uint64_t dfs_guard, 
        double path_diff, uint32_t stall_round_limit_, bool skip_comp, bool keep_nested, 
        double same_sim, uint32_t same_min_len,
        uint32_t diff_min_src, double diff_sim, uint32_t diff_min_len, 
        uint32_t homo_num, uint32_t homo_k, uint32_t homo_w, 
        uint32_t homo_extend_bp, uint64_t homo_bloom_bits, uint32_t homo_bloom_hash, 
        uint32_t thread
    );
    
    void find_bubbles();
    void filter_nested_bubbles();
    const std::vector<Bubble>& get_bubbles() const noexcept { return bubbles_; }
    void save_bubble_as_gfa(const std::string& output_file, const uint32_t min_len, const uint32_t min_num, bool write_seq = true, const std::string& command_line = "") const;
    void save_bubble_as_vcf(const std::string& output_prefix, const std::string& paf_file, const std::string& ref_file, uint32_t min_mapq, uint32_t min_aln_len) const;
    void print_bubbles() const;

    void find_forks();
    const std::vector<ForkGroup>& get_forks() const { return forks_; }
    void print_forks() const;

    void find_homologous_paths();
    const std::vector<HomologousPath>& get_homologous_paths() const { return homologous_paths_; }
    void print_homologous_paths() const;
    void label_path_clusters();

    std::vector<uint32_t> collect_orientation_conflict_segments() const;

    uint64_t hash_vec(std::vector<uint32_t> vertexs, bool use_direction = true) { return hash_vec_(std::move(vertexs), use_direction); }
    static std::vector<uint32_t> trim_path_cluster_anchors(const std::vector<uint32_t>& path, Type type);

private:
    uint64_t hash_vec_(std::vector<uint32_t> vertexs, bool use_direction = true) const;
    bool is_strict_source_(Vertex vtx, std::vector<Vertex>& tgt) const;  // Collect first‑step targets of v (unique by segment id)
    void collect_strict_sources_tips_();
    bool is_strict_sink_(Vertex vtx, Vertex w) const;  // Check whether w is a strict sink

    bool find_common_nearest_sink_(Vertex source, const std::vector<Vertex>& seeds, uint64_t bfs_limit, uint32_t& best_sink, std::unordered_set<uint32_t>& local_nodes);

    Bubble detect_closed_bubble_from_source_(Vertex src, uint64_t bfs_limit);
    GfaBubble::BubbleBranchPairs_ build_bubble_branch_pairs_() const;

    std::vector<Vertex> collect_sources_for_homo_() const;
    std::vector<minimizerdna::Sketch> build_source_mm_library_(const std::vector<Vertex>& sources) const;
    std::vector<NodeGroupEntry> cluster_sources_by_mm_(const std::vector<Vertex>& sources, const std::vector<minimizerdna::Sketch>& node_sketches) const;
    std::vector<NodeGroupEntry> dedup_adjacent_source_groups_(std::vector<NodeGroupEntry> groups) const;
    std::vector<HomoTask_> build_homo_tasks_(std::vector<NodeGroupEntry> groups) const;
    void sort_and_dedup_homologous_paths_(const HomologousParam& same_params, const HomologousParam& diff_params);
    void filter_diff_source_context_();
    HomologousPath detect_homologous_paths_from_sources_(const std::vector<Vertex>& srcs, const HomologousParam& params, const GfaBubble::BubbleBranchPairs_& bubble_pairs) const;

private:
    const uint32_t    thread_;
    const GfaGraph&   graph_;
    const uint32_t    DEPTH_MARGIN_ = 50;  // Margin added to BFS max sink depth to limit DFS path enumeration

    const uint32_t    max_depth_;
    const uint16_t    max_paths_;
    const uint64_t    dfs_guard_;           // --DFS-guard: max DFS states (stop if exceeded)
    const double      path_diff_;           // --path-diff: diff ratio threshold (<= is same cluster)
    const uint32_t    stall_round_limit_;   // stop DFS path search after this many rounds without a new unique path

    const bool        skip_comp_;

    const bool        keep_nested_;         // whether to keep bubbles contained within larger bubbles (i.e. non-local bubbles)

    // Homologous path parameters
    double same_sim_;                       // --same_sim: minimum similarity for homologous paths from the same source
    uint32_t same_min_len_;                 // --same_len: minimum total length of each same-source homologous path
    uint32_t diff_min_src_;                 // --diff_src_len: minimum source length for searching homologous paths between different sources
    double diff_sim_;                       // --diff_sim: minimum similarity for homologous paths between different sources
    uint32_t diff_min_len_;                 // --diff_len: minimum total length of each different-source homologous path
    uint32_t homo_num_{1};                  // --homo_num: max paths per source branch
    uint32_t homo_extend_bp_{1000000};      // --hextend: bp extended per homologous-path extension round
    uint64_t homo_bloom_bits_{1ULL << 27};  // --hbits: Bloom filter size in bits
    uint32_t homo_bloom_hash_{4};           // --hhash: Bloom filter hash count

    // Used to calculate Jaccard similarity between two homologous paths
    const minimizerdna::Options mm_opt_;

private:
    std::vector<Vertex> sources_, tips_;
    std::vector<Bubble> bubbles_;
    GfaBubble::BubbleBranchPairs_ bubble_branch_pairs_;
    std::vector<ForkGroup> forks_;
    std::vector<HomologousPath> homologous_paths_;
};


/* ================================================================================================================
 *                                         HOMOLOGOUS PATH ENUMERATION START
 * ================================================================================================================ */
class HomologousPathEnumeratorLogger {
public:
    explicit HomologousPathEnumeratorLogger(const GfaGraph& graph);

    void header(const std::vector<Vertex>& srcs) const;
    void starts(const std::vector<Vertex>& starts) const;
    void compressed_starts(size_t before, const std::vector<Vertex>& starts) const;

    void similarity(
        double min_first,
        double max_first,
        double min_second,
        double max_second,
        double min_hi,
        double max_hi,
        bool any_pair_good
    ) const;

    void return_reason(const std::string& reason) const;
    void break_reason(const std::string& reason) const;
    void sequence(const std::string& seq) const;

    void kept(
        const std::vector<Vertex>& starts,
        const std::vector<Vertex>& ends,
        const std::vector<uint64_t>& lens,
        const std::vector<uint8_t>& keep
    ) const;

    void enumerate(
        Vertex start,
        Vertex end,
        const std::vector<std::vector<uint32_t>>& paths
    ) const;

    void stop_final_intersected_pair(
        Vertex a_start,
        Vertex a_end,
        Vertex b_start,
        Vertex b_end,
        uint32_t rank
    ) const;

    void stop_intersected_branch(
        Vertex start,
        Vertex end,
        Vertex by,
        uint32_t rank
    ) const;

private:
    const GfaGraph& graph_;
};

class HomologousPathEnumerator {
public:
    using PathSketch = fastbloom::BloomSketch<minimizerdna::Sketch>;
    using RegionSet = std::unordered_set<uint32_t>;

    HomologousPathEnumerator(
        const GfaBubbleFinder& finder,
        const std::vector<Vertex>& srcs,
        const HomologousParam& params,
        const BubbleBranchPairs_& bubble_pairs
    );

    HomologousPath run();

private:
    static constexpr uint8_t PATH_DEAD   = 0;  // stop moving and don't participate in similarity checks
    static constexpr uint8_t PATH_ACTIVE = 1;  // movable
    static constexpr uint8_t PATH_FROZEN = 2;  // stop moving, but still participate in similarity checks

    bool share_one_bubble_group_(const std::vector<Vertex>& vs) const;
    bool prepare_starts_();
    bool initialize_states_();
    bool bootstrap_();

    bool extend_path_(size_t i, bool include_current);

    std::pair<double, double> pair_containments_(size_t i, size_t j) const;
    bool pair_supported_(size_t i, size_t j, const std::pair<double, double>& c) const;

    bool pair_can_still_change_(
        size_t i,
        size_t j,
        const std::vector<uint8_t>& still_supported,
        const std::vector<uint8_t>& need_move
    ) const;

    size_t pair_cache_index_(size_t i, size_t j) const;

    bool compute_pair_stats_(double& min_hi, double& max_hi);

    void collect_alive_active_(std::vector<size_t>& alive_idx, std::vector<size_t>& active_idx) const;

    bool plan_moves_(
        const std::vector<size_t>& alive_idx,
        std::vector<uint8_t>& need_move,
        std::vector<uint8_t>& still_supported,
        std::vector<std::pair<double, double>>& pair_cache,
        std::vector<uint8_t>& pair_cache_valid
    );

    bool extend_needed_paths_(
        const std::vector<size_t>& active_idx,
        const std::vector<uint8_t>& need_move,
        std::vector<uint8_t>& can_move,
        std::vector<Vertex>& next_vtxs,
        std::vector<uint64_t>& next_bp,
        std::vector<std::unique_ptr<PathSketch>>& next_sketch,
        std::vector<std::unique_ptr<std::unordered_set<uint32_t>>>& next_visited,
        std::vector<std::unique_ptr<std::unordered_set<uint32_t>>>& next_visited_seg,
        std::vector<std::unique_ptr<std::vector<uint32_t>>>& next_walks
    );

    bool evaluate_next_state_(
        const std::vector<size_t>& alive_idx,
        const std::vector<uint8_t>& can_move,
        const std::vector<std::unique_ptr<PathSketch>>& next_sketch,
        const std::vector<std::pair<double, double>>& pair_cache,
        const std::vector<uint8_t>& pair_cache_valid,
        std::vector<uint8_t>& next_supported,
        double& next_min_hi,
        double& next_max_hi
    ) const;

    void commit_next_state_(
        const std::vector<size_t>& active_idx,
        const std::vector<uint8_t>& can_move,
        const std::vector<Vertex>& next_vtxs,
        const std::vector<uint64_t>& next_bp,
        std::vector<std::unique_ptr<PathSketch>>& next_sketch,
        std::vector<std::unique_ptr<std::unordered_set<uint32_t>>>& next_visited,
        std::vector<std::unique_ptr<std::unordered_set<uint32_t>>>& next_visited_seg,
        std::vector<std::unique_ptr<std::vector<uint32_t>>>& next_walks,
        double next_min_hi,
        double next_max_hi
    );

    void kill_unsupported_(
        const std::vector<size_t>& alive_idx,
        const std::vector<uint8_t>& supported
    );

    void stop_intersected_same_source_paths_();

    HomologousPath emit_result_();

private:
    uint32_t DFS_DEPTH_;
    uint64_t DFS_GUARD_;

    const GfaBubbleFinder& finder_;
    const GfaGraph& graph_;
    HomologousPathEnumeratorLogger log_;

    const std::vector<Vertex>& srcs_;
    const HomologousParam& params_;
    const BubbleBranchPairs_& bubble_pairs_;

    HomologousPath out_;

    std::vector<Vertex> starts_;
    std::vector<Vertex> owners_;

    std::vector<Vertex> cur_;
    std::vector<Vertex> ends_;
    std::vector<uint64_t> acc_bp_;
    std::vector<uint8_t> active_;
    std::vector<uint8_t> keep_;
    std::vector<std::unordered_set<uint32_t>> visited_;
    std::vector<std::unordered_set<uint32_t>> visited_seg_;
    std::vector<std::vector<uint32_t>> walks_;
    std::vector<PathSketch> sketches_;
};
/* ================================================================================================================
 *                                          HOMOLOGOUS PATH ENUMERATION END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                         HOMOLOGOUS PATH DEDUPLICATION START
 * ================================================================================================================ */
class HomologousPathDeduplicator {
public:
    HomologousPathDeduplicator(
        const GfaGraph& graph,
        const HomologousParam& same_params,
        const HomologousParam& diff_params
    );

    void run(std::vector<HomologousPath>& paths) const;

private:
    struct Candidate {
        uint64_t size{0};
        uint32_t tie_vertex{UINT32_MAX};
        HomologousPath hp;
    };

    struct RunState {
        uint32_t last_pos{UINT32_MAX};
        uint32_t last_old_pos{UINT32_MAX};
        uint32_t cur_beg{0};
        uint32_t cur_len{0};
        uint32_t best_beg{0};
        uint32_t best_len{0};
    };

    struct CoveredRun {
        uint32_t group_id{UINT32_MAX};
        uint32_t path_id{UINT32_MAX};
        uint32_t beg{0};
        uint32_t len{0};

        uint32_t end() const { return beg + len; }
    };

    struct TrimPlan {
        uint32_t group_id{UINT32_MAX};
        std::vector<CoveredRun> runs;

        bool valid() const {
            return group_id != UINT32_MAX && !runs.empty();
        }
    };

    struct NodeHit {
        uint32_t group_id{UINT32_MAX};
        uint32_t path_id{UINT32_MAX};
        uint32_t pos{UINT32_MAX};
    };

private:
    uint64_t group_size_(const HomologousPath& hp) const;
    uint32_t tie_vertex_(const HomologousPath& hp) const;
    uint64_t path_bp_(const std::vector<Vertex>& path) const;

    const HomologousParam& params_for_(Type type) const;

    static bool candidate_better_(const Candidate& a, const Candidate& b);

    std::vector<CoveredRun> covered_runs_for_path_(
        const std::vector<Vertex>& path,
        uint32_t path_id,
        const std::unordered_map<uint32_t, std::vector<NodeHit>>& node_index
    ) const;

    bool find_trim_plan_(
        const HomologousPath& hp,
        const std::unordered_map<uint32_t, std::vector<NodeHit>>& node_index,
        TrimPlan& plan
    ) const;

    bool build_trimmed_piece_(
        const HomologousPath& hp,
        const TrimPlan& plan,
        bool keep_prefix,
        HomologousPath& out
    ) const;

    bool append_trimmed_path_(
        const HomologousPath& old_hp,
        size_t i,
        const std::vector<Vertex>& path,
        HomologousPath& out
    ) const;

    bool passes_thresholds_(
        const HomologousPath& hp
    ) const;

    void index_group_(
        const HomologousPath& hp,
        uint32_t group_id,
        std::unordered_map<uint32_t, std::vector<NodeHit>>& node_index
    ) const;

private:
    const GfaGraph& graph_;
    const HomologousParam& same_params_;
    const HomologousParam& diff_params_;
};
/* ================================================================================================================
 *                                          HOMOLOGOUS PATH DEDUPLICATION END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                         PATH CLUSTERING START
 * ================================================================================================================ */
class PathClusterer {
public:
    PathClusterer(const GfaGraph& graph, double max_difference, uint32_t kmer, uint32_t window);

    std::vector<uint32_t> cluster(const std::vector<std::vector<uint32_t>>& paths, Type type) const;

private:
    const GfaGraph& graph_;
    const double max_difference_;
    const uint32_t kmer_;
    const uint32_t window_;
};
/* ================================================================================================================
 *                                         PATH CLUSTERING END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                         BUBBLE WRITER START
 * ================================================================================================================ */
class BubbleWriter {
public:
    BubbleWriter(
        const GfaGraph& graph,
        const std::vector<Bubble>& bubbles,
        const std::vector<HomologousPath>& homologous_paths
    );

    void save_gfa(
        const std::string& output_file,
        uint32_t min_len,
        uint32_t min_num,
        bool write_seq = true,
        const std::string& command_line = ""
    ) const;

    void save_vcf(
        const std::string& output_prefix,
        const std::string& paf_file,
        const std::string& ref_file,
        uint32_t min_mapq,
        uint32_t min_aln_len
    ) const;

private:
    struct BoundaryWindow {
        std::string chrom;
        uint32_t win_beg{0};
        uint32_t win_end{0};
        uint32_t boundary{0};
        bool boundary_is_left{true};
    };

    std::string node_variant_types_(uint32_t sid) const;
    std::string collect_variant_types_(const std::unordered_set<uint32_t>& segments) const;

    static std::string join_strings_(const std::vector<std::string>& values, std::string_view separator);
    std::string path_name_u32_(const std::vector<uint32_t>& path, uint32_t skip_a = UINT32_MAX, uint32_t skip_b = UINT32_MAX) const;
    std::string path_name_vtx_(const std::vector<Vertex>& path) const;
    static std::string overlap_string_(const std::vector<uint32_t>& overlaps);
    uint64_t path_span_u32_(const std::vector<uint32_t>& path, const std::vector<uint32_t>& overlaps, uint32_t skip_a, uint32_t skip_b) const;
    uint64_t path_span_vtx_(const std::vector<Vertex>& path, const std::vector<uint32_t>& overlaps) const;

    bool boundary_window_(const std::string& name_with_sign, bool is_source, uint64_t plain_len, BoundaryWindow& out) const;
    static bool map_boundary_(const liftover::LIFTresult& hit, const BoundaryWindow& window, uint32_t& ref_pos);
    static bool first_base_all_same_(const std::string& ref, const std::vector<std::string>& alleles);
    static bool can_trim_prefix_(const std::string& ref, const std::vector<std::string>& alleles);
    static bool can_trim_suffix_(const std::string& ref, const std::vector<std::string>& alleles);
    std::string path_root_(const std::string& path_name) const;
    std::vector<uint32_t> path_sample_ids_(const std::vector<uint32_t>& path, uint32_t src_seg, uint32_t sink_seg) const;
    std::vector<std::string> sample_gt_tokens_(const std::vector<std::string>& path_gt_tokens, const std::vector<std::vector<uint32_t>>& path_sample_ids) const;
    static std::string vcf_record_string_(const VcfRecord& record);
    void normalize_vcf_record_(std::vector<VcfRecord>& records, size_t index) const;

    const GfaGraph& graph_;
    const std::vector<Bubble>& bubbles_;
    const std::vector<HomologousPath>& homologous_paths_;
};
/* ================================================================================================================
 *                                         BUBBLE WRITER END
 * ================================================================================================================ */

} // namespace GfaBubble
