#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "gfa_parser.hpp"
#include "gfa_parser_types.hpp"

class GfaWalkerLogger;

struct GreedyFrame {
    uint32_t v{UINT32_MAX};
    std::vector<uint32_t> cand;
    size_t next{0};
};


/* ================================================================================================================
 *                                      GREEDY BRANCH TRACKER START
 * ================================================================================================================ */
class GreedyBranchTracker {
public:
    static constexpr uint64_t NO_BRANCH = UINT64_MAX;

    explicit GreedyBranchTracker(uint32_t stall_limit);

    bool disabled() const;
    bool has_branch(uint64_t branch) const;
    uint32_t last_stall() const;

    uint64_t make_branch(uint32_t src, uint32_t next) const;
    void remove_closed(uint32_t src, std::vector<uint32_t>& cand) const;
    uint32_t stall(uint64_t branch);
    void reset(uint64_t branch);

private:
    uint32_t stall_limit_{0};
    uint32_t last_stall_{0};
    std::unordered_map<uint64_t, uint32_t> stall_rounds_;
    std::unordered_set<uint64_t> closed_branches_;
};
/* ================================================================================================================
 *                                       GREEDY BRANCH TRACKER END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                        GREEDY DFS CONTEXT START
 * ================================================================================================================ */
class GreedyDfsContext {
public:
    GreedyDfsContext(
        const GfaGraph& graph,
        uint32_t src,
        uint32_t sink,
        const std::unordered_set<uint32_t>& region_set,
        bool skip_comp,
        const GfaWalkerLogger& log
    );

    bool use_topo_filter() const;
    uint32_t src_rank() const;
    uint32_t sink_rank() const;

    bool in_topo_window(uint32_t v) const;
    bool valid_next(uint32_t from, uint32_t w, const std::unordered_set<uint32_t>& on_path_seg) const;
    bool path_has_new_node(const std::vector<uint32_t>& path, const std::unordered_set<uint32_t>& covered) const;

    std::vector<uint32_t> collect_candidates(
        uint32_t from,
        const std::unordered_set<uint32_t>& on_path_seg,
        const std::unordered_set<uint32_t>& covered,
        const std::unordered_map<uint64_t, uint32_t>& edge_try_count
    ) const;

private:
    std::string reject_reason(uint32_t from, uint32_t w, const std::unordered_set<uint32_t>& on_path_seg) const;
    void build_sink_reachable_();

    void sort_candidates(
        uint32_t from,
        std::vector<uint32_t>& cand,
        const std::unordered_set<uint32_t>& covered,
        const std::unordered_map<uint64_t, uint32_t>& edge_try_count
    ) const;

private:
    const GfaGraph& graph_;
    const GfaWalkerLogger& log_;
    uint32_t src_{UINT32_MAX};
    uint32_t sink_{UINT32_MAX};
    const std::unordered_set<uint32_t>& region_set_;
    bool skip_comp_{false};

    bool use_topo_filter_{false};
    uint32_t src_rank_{UINT32_MAX};
    uint32_t sink_rank_{UINT32_MAX};

    std::unordered_set<uint32_t> sink_reachable_;
};
/* ================================================================================================================
 *                                         GREEDY DFS CONTEXT END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                            GFA WALKER START
 * ================================================================================================================ */
class GfaWalker {
public:
    explicit GfaWalker(const GfaGraph& graph);

    std::vector<std::vector<uint32_t>> enumerate_paths_greedy_DFS(
        uint32_t src,
        uint32_t sink,
        const std::unordered_set<uint32_t>& region_set,
        uint32_t max_depth,
        uint32_t max_paths,
        bool skip_comp,
        bool& hit_limits,
        uint64_t DFS_guard,
        uint32_t stall_round_limit = 2
    ) const;

    std::vector<std::vector<uint32_t>> open_walk(
        uint32_t src,
        const std::unordered_set<uint32_t>& region_set,
        uint32_t max_depth,
        bool skip_comp,
        uint64_t DFS_guard,
        uint64_t walk_bp,
        const std::unordered_set<uint32_t>* blocked_seg = nullptr
    ) const;

private:
    static uint64_t hash_vec_(std::vector<uint32_t> vec, bool use_direction = true);

private:
    const GfaGraph& graph_;
};
/* ================================================================================================================
 *                                            GFA WALKER END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                       SOURCE-AWARE WALKER START
 * ================================================================================================================ */
class GfaSourceWalker {
public:
    struct Options {
        // Sample IDs come from GfaGraph::getSampleNames(); source_scope excludes non-haplotype SN labels.
        uint32_t preferred_source{UINT32_MAX};
        const std::unordered_set<uint32_t>* source_scope{nullptr};
        const std::unordered_map<uint64_t, std::vector<uint32_t>>* transition_sources{nullptr}; // oriented P-path edge -> source IDs
        const std::vector<uint8_t>* soft_blocked{nullptr};  // dense bit flags indexed by unoriented segment ID
        uint8_t soft_block_mask{1};
        uint32_t max_depth{4096};
        uint64_t max_bp{UINT64_MAX}, max_states{1000000};
        bool skip_comp{false};
    };

    struct Result {
        std::vector<uint32_t> vertices;  // source and sink included
        uint64_t states{0};
        bool truncated{false};

        explicit operator bool() const { return !vertices.empty(); }
    };

    explicit GfaSourceWalker(const GfaGraph& graph);
    Result walk(uint32_t source, uint32_t sink, const Options& options) const;
    Result walk(uint32_t source, uint32_t sink) const { return walk(source, sink, Options{}); }

private:
    const GfaGraph& graph_;
};
/* ================================================================================================================
 *                                        SOURCE-AWARE WALKER END
 * ================================================================================================================ */
