#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <queue>
#include <stack>
#include <string>
#include <unordered_set>
#include <vector>

#include "logger.hpp"
#include "gfa_parser.hpp"
#include "gfa_bubble_types.hpp"
#include "save.hpp"

namespace GfaBubble {

class GfaBubbleFinder {
public:
    /**
     * @param graph      Reference to an immutable GfaGraph.
     * @param max_depth  DFS/BFS depth cap when enumerating each path.
     * @param max_paths  Per‑bubble path enumeration cap.
     * @param skip_comp  If true, ignore arcs marked as "comp" (compressed).
     * @param thread     Number of threads to use for parallel processing.
     */
    GfaBubbleFinder(
        const GfaGraph& graph, std::size_t max_depth, std::size_t max_paths, uint64_t dfs_guard, 
        uint16_t cx_depth, uint32_t cx_nodes, uint32_t cx_branches, int cx_deg_branch, int cx_deg_hub, int cx_deg_cap, 
        double path_diff, bool skip_comp, uint32_t thread
    );
    
    void find_bubbles(const bool filter_nonlocal = true);  // Run bubble detection over the whole graph.
    const std::vector<Bubble>& get_bubbles() const noexcept { return bubbles_; }  // getters
    void save_bubble_as_gfa(const std::string& filename, const uint32_t min_len, const uint32_t min_num) const;  // Save all detected bubbles as an independent GFA sub‑graph.
    void print_bubbles() const;  // Print bubbles

    void find_forks();  // Find all strict forks in the graph
    const std::vector<ForkGroup>& get_forks() const { return forks_; }
    void print_forks() const;  // Print strict forks

    // cluster paths by node-length symmetric-difference ratio, pick longest rep per cluster
    std::vector<std::vector<uint32_t>> pick_representative_paths(
        const std::vector<std::vector<uint32_t>>& paths
    ) const ;

    static uint64_t hash_vec(std::vector<uint32_t> vec, bool use_direction = true) { return hash_vec_(std::move(vec), use_direction); }

private:
    static uint64_t hash_vec_(std::vector<uint32_t> vec, bool use_direction = true) {
        if (!use_direction) {
            for (auto &x : vec) x >>= 1;
        }
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

        uint64_t h = 14695981039346656037ULL;  // FNV offset basis
        for (uint32_t v : vec) {
            h ^= static_cast<uint64_t>(v);
            h *= 1099511628211ULL;  // FNV prime
        }
        return h;
    }

    // Collect first‑step targets of v (unique by segment id).
    bool is_strict_source_(uint32_t v, std::vector<uint32_t>& tgt) const;
    // Check whether w is a strict sink
    bool is_strict_sink_(uint32_t v, uint32_t w) const;

    // Mark nodes in 'nodes' as complex
    void mark_nodes_complex_(const std::vector<uint32_t>& nodes);
    void mark_nodes_simple_(const std::vector<uint32_t>& nodes);
    bool is_complex_source_(uint32_t v, const std::vector<uint32_t>& seed_buf);

    // Closed-bubble
    Bubble detect_closed_bubble_from_source_(uint32_t src, uint64_t bfs_limit);
    // Open-bubble (Tip)
    void find_open_bubbles_from_tips_(const std::vector<uint32_t>& tip_nodes);
    
    // Filter out non-local bubbles
    void filter_nonlocal_bubbles_();

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
    const int         cx_deg_cap_;         // --cx-deg-cap: stop counting degree after this
    const double      path_diff_;          // --path-diff: diff ratio threshold (<= is same cluster)

    const bool        skip_comp_;

    // 0 = unknown, 1 = complex, 2 = simple
    static constexpr uint8_t UNKNOWN_FLAG_  = 0;
    static constexpr uint8_t COMPLEX_FLAG_  = 1;
    static constexpr uint8_t SIMPLE_FLAG_   = 2;
    mutable std::vector<uint8_t> complex_flag_;
    std::mutex complex_mtx_;

    // scratch buffers
    static constexpr uint8_t UNVIS_FLAG_  = 0xFFu;  // unvisited flag
    static constexpr uint8_t MERGED_FLAG_ = 0xFEu;  // merged flag, used to find larger bubbles
    mutable std::vector<uint8_t> reach_buf_;   // BFS branch tag (UNVIS_FLAG_ == unvisited, MERGED_FLAG_ == merged)

    std::vector<Bubble> bubbles_;
    std::vector<ForkGroup> forks_;
};

} // namespace GfaBubble
