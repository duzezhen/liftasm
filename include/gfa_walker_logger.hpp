#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../include/gfa_parser.hpp"
#include "../include/gfa_parser_types.hpp"

class GfaWalkerLogger {
public:
    explicit GfaWalkerLogger(const GfaGraph& graph, bool enabled = DEBUG_ENABLED);

    bool enabled() const;
    std::string indent(int level) const;
    std::string node_name(uint32_t v) const;

    void path(const std::vector<uint32_t>& path, int level, const std::string& title) const;

    void greedy_header(
        uint32_t src,
        uint32_t sink,
        uint32_t max_depth,
        uint32_t max_paths,
        uint64_t dfs_guard,
        uint32_t stall_limit,
        bool use_topo,
        uint32_t src_rank,
        uint32_t sink_rank,
        size_t region_size
    ) const;

    void greedy_round(
        uint32_t round,
        size_t path_count,
        size_t covered_count,
        uint32_t stall_rounds,
        uint64_t dfs_states
    ) const;

    void greedy_start(uint32_t src) const;
    void greedy_expand(uint32_t v, size_t depth, size_t out_arcs) const;
    void greedy_candidate(uint32_t w, uint32_t len, bool used, uint32_t edge_try) const;
    void greedy_skip(uint32_t from, uint32_t to, const std::string& reason) const;
    void greedy_sorted_candidates(
        const std::vector<uint32_t>& cand,
        uint32_t from,
        const std::unordered_set<uint32_t>& covered,
        const std::unordered_map<uint64_t, uint32_t>& edge_try_count
    ) const;
    void greedy_choose(uint32_t w, size_t rank, size_t total) const;
    void greedy_backtrack(uint32_t v, const std::string& reason) const;
    void greedy_record(const std::vector<uint32_t>& path, size_t path_id, bool has_new_node) const;
    void greedy_discard(const std::vector<uint32_t>& path, const std::string& reason, uint32_t next_stall) const;
    void greedy_summary(size_t path_count, size_t covered_count, uint64_t dfs_states, uint32_t rounds) const;

    void open_header(uint32_t src, uint64_t walk_bp, uint32_t max_depth, uint64_t dfs_guard) const;
    void open_step(uint32_t from, uint32_t to) const;
    void open_summary(const std::vector<uint32_t>& path, uint64_t bp) const;

private:
    uint32_t edge_try(uint32_t from, uint32_t to, const std::unordered_map<uint64_t, uint32_t>& edge_try_count) const;

private:
    const GfaGraph& graph_;
    bool enabled_{false};
};
