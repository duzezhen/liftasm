#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <queue>
#include <stack>
#include <string>
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

class GfaBubbleFinder {
public:
    GfaBubbleFinder(
        const GfaGraph& graph, std::size_t max_depth, std::size_t max_paths, uint64_t dfs_guard, 
        uint16_t cx_depth, uint32_t cx_nodes, uint32_t cx_branches, int cx_deg_branch, int cx_deg_hub,
        double path_diff, uint32_t stall_round_limit_, bool skip_comp, bool keep_nested, 
        double same_sim, uint32_t same_min_num, uint32_t same_min_len,
        uint32_t diff_min_src, double diff_sim, uint32_t diff_min_num, uint32_t diff_min_len, 
        uint32_t homo_num, uint32_t homo_k, uint32_t homo_w, 
        uint32_t homo_boot_bp, uint32_t homo_cnt_len, uint64_t homo_bloom_bits, uint32_t homo_bloom_hash,
        uint32_t thread
    );
    
    void find_bubbles();
    const std::vector<Bubble>& get_bubbles() const noexcept { return bubbles_; }
    void save_bubble_as_gfa(const std::string& output_file, const uint32_t min_len, const uint32_t min_num, bool write_seq = true) const;
    void save_bubble_as_vcf(const std::string& output_prefix, const std::string& paf_file, const std::string& ref_file, uint32_t min_mapq, uint32_t min_aln_len) const;
    void print_bubbles() const;

    void find_forks();
    const std::vector<ForkGroup>& get_forks() const { return forks_; }
    void print_forks() const;

    void find_homologous_paths();
    const std::vector<HomologousPath>& get_homologous_paths() const { return homologous_paths_; }
    void print_homologous_paths() const;

    // cluster paths by node-length symmetric-difference ratio, pick longest rep per cluster
    std::vector<std::vector<uint32_t>> pick_representative_paths(const std::vector<std::vector<uint32_t>>& paths) const ;
    std::vector<std::vector<uint32_t>> pick_representative_paths(const std::vector<std::vector<Vertex>>& paths) const ;

    // Pack a vector of vertex IDs into a hash value.
    uint64_t hash_vec(std::vector<uint32_t> vertexs, bool use_direction = true) { return hash_vec_(std::move(vertexs), use_direction); }

private:
    uint64_t hash_vec_(std::vector<uint32_t> vertexs, bool use_direction = true) const;
    bool is_single_node_(Vertex v) const;
    bool is_strict_source_(Vertex vtx, std::vector<Vertex>& tgt) const;  // Collect first‑step targets of v (unique by segment id)
    bool is_complex_source_(Vertex vtx, const std::vector<Vertex>& seed_buf);
    void collect_strict_sources_tips_();
    bool is_strict_sink_(Vertex vtx, Vertex w) const;  // Check whether w is a strict sink

    // Mark nodes in 'nodes' as complex
    void mark_nodes_complex_(const std::vector<Vertex>& nodes);
    void mark_nodes_simple_(const std::vector<Vertex>& nodes);

    /**
     * @brief Find the nearest common sink reachable by all branches from a source.
     * @date 2026-04-16
     * @version 0.1.3
     *
     * @param source          the source vertex of the bubble
     * @param seeds           starting vertices (branches) from the source
     * @param bfs_limit       maximum BFS depth limit
     * @param best_sink       (output) selected common sink vertex id
     * @param best_sink_depth (output) depth of the selected sink
     * @param local_nodes     (output) nodes visited during BFS (for DFS use)
     *
     * @return true if a valid common sink is found; false otherwise
     */
    bool find_common_nearest_sink_(
        Vertex source,
        const std::vector<Vertex>& seeds,
        uint64_t bfs_limit,
        uint32_t& best_sink,
        uint32_t& best_sink_depth,
        std::unordered_set<uint32_t>& local_nodes
    );

    /**
     * @brief Detect closed bubbles starting from a source vertex.
     * @date 2026-04-16
     * @version 0.1.3
     *
     * @param src        the source vertex of the bubble
     * @param bfs_limit  maximum BFS depth limit
     *
     * @return a Bubble object containing the detected bubble information; empty if no valid bubble is found
     */
    Bubble detect_closed_bubble_from_source_(Vertex src, uint64_t bfs_limit);

    // Filter out non-local bubbles
    void filter_nonlocal_bubbles_();

    /**
     * @brief Build a set of all bubble branch pairs (fitst and second vertex), 
     */
    std::unordered_set<uint64_t> build_bubble_branch_pairs_() const;


    std::vector<Vertex> collect_sources_for_homo_();


    std::vector<minimizerdna::Sketch> build_all_node_mm_library_() const;

    /**
     * @param source         the source vertex for clustering
     * @param node_sketches  vector of all node sketches
     * @param NodeMMLibs     library of node sketches
     * @param sim_thr        similarity threshold for clustering (0-1)
     * 
     * @return a list of NodeGroupEntry, each containing a group of source nodes clustered by minimizer similarity
     */
    std::vector<NodeGroupEntry> cluster_sources_by_mm_(
        const std::vector<Vertex>& sources,
        const std::vector<minimizerdna::Sketch>& node_sketches,
        uint32_t min_len,
        double sim_thr
    ) const;


    bool homo_vertices_already_used_(
        const std::vector<Vertex>& vertices,
        const HomoUsedIntervals_& used_intervals
    ) const;

    bool homo_path_already_used_(
        const HomologousPath& hp,
        const HomoUsedIntervals_& used_intervals
    ) const;

    void insert_homo_used_interval_(
        HomoUsedIntervals_& used_intervals,
        const HomoUsedKey_& key,
        HomoUsedInterval_ interval
    ) const;

    void mark_homo_used_intervals_(
        HomoUsedIntervals_& used_intervals,
        const HomologousPath& hp
    ) const;

    std::vector<HomoTaskGroup_> build_homo_tasks_(
        const std::vector<NodeGroupEntry>& groups, 
        bool no_grouping = false
    ) const;

    HomoRunResult_ find_homologous_paths_run_(
        const HomoTaskGroup_& group,
        const HomologousParam& same_params,
        const HomologousParam& diff_params,
        const std::unordered_set<uint64_t>& bubble_pairs,
        const std::vector<minimizerdna::Sketch>& node_sketches
    ) const;

    /**
     * @brief Detect homologous paths grown from one or more source vertices.
     *
     * @param srcs            Source vertices used to seed homologous path detection.
     * @param params          Homologous detection parameters.
     * @param bubble_pairs    Hash set of already known bubble/source pairs used to skip or compress redundant starts.
     * @param node_sketches   Precomputed minimizer sketches indexed by oriented vertex ID.
     * @param used_intervals  Previously used homologous topological intervals.
     *
     * @return HomologousPath containing detected homologous path metadata and enumerated vertex paths.
     */
    HomologousPath detect_homologous_paths_from_sources_(
        const std::vector<Vertex>& srcs,
        const HomologousParam& params,
        const std::unordered_set<uint64_t>& bubble_pairs,
        const std::vector<minimizerdna::Sketch>& node_sketches,
        const HomoUsedIntervals_& used_intervals
    ) const;

private:
    const uint32_t    thread_;
    const GfaGraph&   graph_;
    const uint32_t    DEPTH_MARGIN_ = 50;  // Margin added to BFS max sink depth to limit DFS path enumeration

    const std::size_t max_depth_;
    const std::size_t max_paths_;
    const uint64_t    dfs_guard_;          // --DFS-guard: max DFS states (stop if exceeded)
    const uint16_t    cx_depth_;           // --cx-depth: BFS depth
    const uint32_t    cx_nodes_;           // --cx-nodes: visited-node cap (complex if exceeded)
    const uint32_t    cx_branches_;        // --cx-branches: branching-node cap (complex if exceeded)
    const int         cx_deg_branch_;      // --cx-deg-branch: degree>=this counts as branching
    const int         cx_deg_hub_;         // --cx-deg-hub: degree>=this is hub => complex
    const double      path_diff_;          // --path-diff: diff ratio threshold (<= is same cluster)
    const uint32_t    stall_round_limit_;  // stop DFS path search after this many rounds without a new unique path

    const bool        skip_comp_;

    const bool        keep_nested_;        // whether to keep bubbles contained within larger bubbles (i.e. non-local bubbles)

    // Homologous path parameters
    double same_sim_;              // --same_sim: minimum similarity for homologous paths from the same source
    uint32_t same_min_num_;        // --same_num: minimum number of nodes in each same-source homologous path
    uint32_t same_min_len_;        // --same_len: minimum total length of each same-source homologous path
    uint32_t diff_min_src_;        // --diff_src_len: minimum source length for searching homologous paths between different sources
    double diff_sim_;              // --diff_sim: minimum similarity for homologous paths between different sources
    uint32_t diff_min_num_;        // --diff_num: minimum number of nodes in each different-source homologous path
    uint32_t diff_min_len_;        // --diff_len: minimum total length of each different-source homologous path
    uint32_t homo_num_{1};         // --homo_num: max paths per source branch
    uint32_t homo_boot_bp_{1000};           // --hboot: bootstrap bp before first similarity check
    uint32_t homo_cnt_len_{1000};           // --hcnt: minimum node length counted as one extension step
    uint64_t homo_bloom_bits_{1ULL << 27};  // --hbits: Bloom filter size in bits
    uint32_t homo_bloom_hash_{4};           // --hhash: Bloom filter hash count


    // Used to calculate Jaccard similarity between two homologous paths
    const minimizerdna::Options mm_opt_;
    const uint32_t overlap_len_;

    // 0 = unknown, 1 = complex, 2 = simple
    static constexpr uint8_t UNKNOWN_FLAG_  = 0;
    static constexpr uint8_t COMPLEX_FLAG_  = 1;
    static constexpr uint8_t SIMPLE_FLAG_   = 2;
    mutable std::vector<uint8_t> complex_flag_;
    std::mutex complex_mtx_;

    // scratch buffers
    static constexpr uint8_t UNVIS_FLAG_  = 0xFFu;  // unvisited flag
    static constexpr uint8_t MERGED_FLAG_ = 0xFEu;  // merged flag, used to find larger bubbles
    mutable std::vector<uint8_t> reach_buf_;        // BFS branch tag (UNVIS_FLAG_ == unvisited, MERGED_FLAG_ == merged)

private:
    std::vector<Vertex> sources_, tips_;
    std::vector<Bubble> bubbles_;
    std::vector<ForkGroup> forks_;
    std::vector<HomologousPath> homologous_paths_;
};

} // namespace GfaBubble
