#include "../include/gfa_walker_logger.hpp"

#include "../include/logger.hpp"


static inline uint64_t edge_key(uint32_t a, uint32_t b) {
    return (uint64_t(a) << 32) | uint64_t(b);
}

GfaWalkerLogger::GfaWalkerLogger(const GfaGraph& graph, bool enabled)
    : graph_(graph),
      enabled_(enabled)
{}

bool GfaWalkerLogger::enabled() const {
    return enabled_;
}

std::string GfaWalkerLogger::indent(int level) const {
    return std::string(level * 2, ' ');
}

std::string GfaWalkerLogger::node_name(uint32_t v) const {
    const uint32_t sid = NodeHandle::get_segment_id(v);
    const GfaNode* node = graph_.getNode(sid);
    return node ? node->name : "*";
}

void GfaWalkerLogger::path(const std::vector<uint32_t>& p, int level, const std::string& title) const {
    if (!enabled_) return;

    debug_stream() << indent(level) << title << "\n";
    for (uint32_t v : p) {
        debug_stream() << indent(level + 1) << node_name(v) << "\n";
    }
}

void GfaWalkerLogger::greedy_header(
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
) const {
    if (!enabled_) return;

    debug_stream() << "=========================================\n";
    debug_stream() << "Enumerate greedy DFS paths\n";
    debug_stream() << indent(1) << "src       : " << node_name(src) << "\n";
    debug_stream() << indent(1) << "sink      : " << node_name(sink) << "\n";
    debug_stream() << indent(1) << "max_depth : " << max_depth << "\n";
    debug_stream() << indent(1) << "max_paths : " << max_paths << "\n";
    debug_stream() << indent(1) << "DFS_guard : " << dfs_guard << "\n";
    debug_stream() << indent(1) << "stall_lim : " << stall_limit << "\n";
    debug_stream() << indent(1) << "region_sz : " << region_size << "\n";

    if (use_topo) {
        debug_stream() << indent(1) << "topo_win  : " << src_rank << " .. " << sink_rank << "\n";
    }
}

void GfaWalkerLogger::greedy_round(
    uint32_t round,
    size_t path_count,
    size_t covered_count,
    uint32_t stall_rounds,
    uint64_t dfs_states
) const {
    if (!enabled_) return;

    debug_stream() << "\n";
    debug_stream() << indent(1) << "round " << round << "\n";
    debug_stream() << indent(2) << "path_count   : " << path_count << "\n";
    debug_stream() << indent(2) << "covered      : " << covered_count << "\n";
    debug_stream() << indent(2) << "stall_rounds : " << stall_rounds << "\n";
    debug_stream() << indent(2) << "dfs_states   : " << dfs_states << "\n";
}

void GfaWalkerLogger::greedy_start(uint32_t src) const {
    if (!enabled_) return;
    debug_stream() << indent(2) << "start greedy DFS from " << node_name(src) << "\n";
}

void GfaWalkerLogger::greedy_expand(uint32_t v, size_t depth, size_t out_arcs) const {
    if (!enabled_) return;
    debug_stream() << indent(2) << "expand " << node_name(v) << "  depth=" << depth << "  out_arcs=" << out_arcs << "\n";
}

void GfaWalkerLogger::greedy_candidate(uint32_t w, uint32_t len, bool used, uint32_t edge_try) const {
    if (!enabled_) return;
    debug_stream() << indent(3) << "candidate " << node_name(w) << "  len=" << len << "  used=" << (used ? "yes" : "no") << "  edge_try=" << edge_try << "\n";
}

void GfaWalkerLogger::greedy_skip(uint32_t from, uint32_t to, const std::string& reason) const {
    if (!enabled_) return;
    debug_stream() << indent(3) << "skip " << node_name(from) << " -> " << node_name(to) << "  reason=" << reason << "\n";
}

uint32_t GfaWalkerLogger::edge_try(
    uint32_t from,
    uint32_t to,
    const std::unordered_map<uint64_t, uint32_t>& edge_try_count
) const {
    const auto it = edge_try_count.find(edge_key(from, to));
    return it == edge_try_count.end() ? 0u : it->second;
}

void GfaWalkerLogger::greedy_sorted_candidates(
    const std::vector<uint32_t>& cand,
    uint32_t from,
    const std::unordered_set<uint32_t>& covered,
    const std::unordered_map<uint64_t, uint32_t>& edge_try_count
) const {
    if (!enabled_ || cand.empty()) return;

    debug_stream() << indent(2) << "sorted candidates:\n";

    for (uint32_t v : cand) {
        const uint32_t sid = NodeHandle::get_segment_id(v);
        debug_stream() << indent(3) << node_name(v) << "  len=" << graph_.getNodeLength(sid) << "  used=" << (covered.find(v) != covered.end() ? "yes" : "no") << "  edge_try=" << edge_try(from, v, edge_try_count) << "\n";
    }
}

void GfaWalkerLogger::greedy_choose(uint32_t w, size_t rank, size_t total) const {
    if (!enabled_) return;
    debug_stream() << indent(2) << "choose candidate: " << node_name(w) << "  rank=" << rank << "/" << total << "\n";
}

void GfaWalkerLogger::greedy_backtrack(uint32_t v, const std::string& reason) const {
    if (!enabled_) return;
    debug_stream() << indent(2) << reason << " at " << node_name(v) << ", backtrack\n";
}

void GfaWalkerLogger::greedy_record(const std::vector<uint32_t>& p, size_t path_id, bool has_new_node) const {
    if (!enabled_) return;

    debug_stream() << indent(2) << "record unique path #" << path_id << "  len=" << p.size() << "  new_node=" << (has_new_node ? "yes" : "no") << "\n";
    path(p, 3, "full path:");
}

void GfaWalkerLogger::greedy_discard(const std::vector<uint32_t>& p, const std::string& reason, uint32_t next_stall) const {
    if (!enabled_) return;

    debug_stream() << indent(2) << "path discarded, " << reason << "\n";
    path(p, 3, "discarded path:");
    debug_stream() << indent(2) << "stall_rounds -> " << next_stall << "\n";
}

void GfaWalkerLogger::greedy_summary(size_t path_count, size_t covered_count, uint64_t dfs_states, uint32_t rounds) const {
    if (!enabled_) return;

    debug_stream() << "\n";
    debug_stream() << "  Greedy enumeration finished\n";
    debug_stream() << indent(2) << "total paths : " << path_count << "\n";
    debug_stream() << indent(2) << "covered     : " << covered_count << "\n";
    debug_stream() << indent(2) << "dfs_states  : " << dfs_states << "\n";
    debug_stream() << indent(2) << "rounds      : " << rounds << "\n";
    debug_stream() << "=========================================\n";
}

void GfaWalkerLogger::open_header(uint32_t src, uint64_t walk_bp, uint32_t max_depth, uint64_t dfs_guard) const {
    if (!enabled_) return;

    debug_stream() << "=========================================\n";
    debug_stream() << "Open walk by topo sinks\n";
    debug_stream() << indent(1) << "src       : " << node_name(src) << "\n";
    debug_stream() << indent(1) << "walk_bp   : " << walk_bp << "\n";
    debug_stream() << indent(1) << "max_depth : " << max_depth << "\n";
    debug_stream() << indent(1) << "DFS_guard : " << dfs_guard << "\n";
}

void GfaWalkerLogger::open_sink(uint32_t from, uint32_t sink) const {
    if (!enabled_) return;
    debug_stream() << indent(1) << "sink: " << node_name(from) << " -> " << node_name(sink) << "\n";
}

void GfaWalkerLogger::open_reject(uint32_t sink, const std::string& reason) const {
    if (!enabled_) return;
    debug_stream() << indent(1) << "reject " << node_name(sink) << "  reason=" << reason << "\n";
}

void GfaWalkerLogger::open_summary(const std::vector<uint32_t>& p, uint64_t bp) const {
    if (!enabled_) return;

    debug_stream() << "  Open walk finished\n";
    debug_stream() << indent(2) << "bp   : " << bp << "\n";
    debug_stream() << indent(2) << "nodes: " << p.size() << "\n";
    path(p, 2, "path:");
    debug_stream() << "=========================================\n";
}
