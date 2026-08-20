#include "../include/gfa_walker.hpp"
#include "../include/gfa_walker_logger.hpp"

#include <algorithm>
#include <cmath>
#include <queue>


static inline uint64_t edge_key(uint32_t a, uint32_t b) {
    return (uint64_t(a) << 32) | uint64_t(b);
}


/* ================================================================================================================
 *                                      GREEDY BRANCH TRACKER START
 * ================================================================================================================ */
GreedyBranchTracker::GreedyBranchTracker(uint32_t stall_limit)
    : stall_limit_(stall_limit)
{
    stall_rounds_.reserve(64);
    closed_branches_.reserve(64);
}

bool GreedyBranchTracker::disabled() const {
    return stall_limit_ == 0;
}

bool GreedyBranchTracker::has_branch(uint64_t branch) const {
    return branch != NO_BRANCH;
}

uint32_t GreedyBranchTracker::last_stall() const {
    return last_stall_;
}

uint64_t GreedyBranchTracker::make_branch(uint32_t src, uint32_t next) const {
    return edge_key(src, next);
}

void GreedyBranchTracker::remove_closed(uint32_t src, std::vector<uint32_t>& cand) const {
    if (disabled()) return;

    size_t keep = 0;
    for (uint32_t w : cand) {
        if (closed_branches_.find(make_branch(src, w)) == closed_branches_.end()) {
            cand[keep++] = w;
        }
    }
    cand.resize(keep);
}

uint32_t GreedyBranchTracker::stall(uint64_t branch) {
    if (disabled()) return last_stall_;
    if (!has_branch(branch)) return last_stall_;

    uint32_t& count = stall_rounds_[branch];
    ++count;
    last_stall_ = count;

    if (count >= stall_limit_) {
        closed_branches_.insert(branch);
    }

    return last_stall_;
}

void GreedyBranchTracker::reset(uint64_t branch) {
    if (disabled()) return;
    if (!has_branch(branch)) return;
    stall_rounds_[branch] = 0;
    last_stall_ = 0;
}
/* ================================================================================================================
 *                                       GREEDY BRANCH TRACKER END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                        GREEDY DFS CONTEXT START
 * ================================================================================================================ */
GreedyDfsContext::GreedyDfsContext(
    const GfaGraph& graph,
    uint32_t src,
    uint32_t sink,
    const std::unordered_set<uint32_t>& region_set,
    bool skip_comp,
    const GfaWalkerLogger& log
)
    : graph_(graph),
      log_(log),
      src_(src),
      sink_(sink),
      region_set_(region_set),
      skip_comp_(skip_comp)
{
    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);

    src_rank_ = graph_.vertex_topo_rank(Vertex(src_));
    sink_rank_ = graph_.vertex_topo_rank(Vertex(sink_));

    use_topo_filter_ = graph_.get_vertex_topological_index().size() >= V && src_rank_ != UINT32_MAX && sink_rank_ != UINT32_MAX;

    if (src_rank_ > sink_rank_) {
        std::swap(src_rank_, sink_rank_);
    }

    build_sink_reachable_();
}

bool GreedyDfsContext::use_topo_filter() const {
    return use_topo_filter_;
}

uint32_t GreedyDfsContext::src_rank() const {
    return src_rank_;
}

uint32_t GreedyDfsContext::sink_rank() const {
    return sink_rank_;
}

bool GreedyDfsContext::in_topo_window(uint32_t v) const {
    if (!use_topo_filter_) return true;

    const uint32_t r = graph_.vertex_topo_rank(Vertex(v));
    return r != UINT32_MAX && r >= src_rank_ && r <= sink_rank_;
}

void GreedyDfsContext::build_sink_reachable_() {
    if (region_set_.empty()) return;

    std::queue<uint32_t> q;
    sink_reachable_.reserve(region_set_.size() + 2);

    sink_reachable_.insert(sink_);
    q.push(sink_);

    while (!q.empty()) {
        const uint32_t v = q.front();
        q.pop();

        const auto& arcs = graph_.getArcsToVertex(v);
        for (const GfaArc* a : arcs) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;

            const uint32_t u = a->get_source_vertex_id();
            if (u != src_ && region_set_.find(u) == region_set_.end()) continue;
            if (!in_topo_window(u)) continue;

            const uint32_t sid = NodeHandle::get_segment_id(u);
            if (graph_.getNodeDeleted(sid)) continue;

            if (sink_reachable_.insert(u).second) {
                q.push(u);
            }
        }
    }
}

std::string GreedyDfsContext::reject_reason(uint32_t from, uint32_t w, const std::unordered_set<uint32_t>& on_path_seg) const {
    if (w == from) return "self-loop";
    if (!in_topo_window(w)) return "outside-topo-window";

    if (!region_set_.empty() && w != sink_ && region_set_.find(w) == region_set_.end()) {
        return "outside-region";
    }

    if (!sink_reachable_.empty() && sink_reachable_.find(w) == sink_reachable_.end()) {
        return "cannot-reach-sink";
    }

    const uint32_t sid = NodeHandle::get_segment_id(w);
    if (graph_.getNodeDeleted(sid)) return "deleted-node";
    if (on_path_seg.find(sid) != on_path_seg.end()) return "cycle-on-path";

    return "";
}

bool GreedyDfsContext::valid_next(uint32_t from, uint32_t w, const std::unordered_set<uint32_t>& on_path_seg) const {
    return reject_reason(from, w, on_path_seg).empty();
}

bool GreedyDfsContext::path_has_new_node(const std::vector<uint32_t>& path, const std::unordered_set<uint32_t>& covered) const {
    for (uint32_t v : path) {
        if (covered.find(v) == covered.end()) return true;
    }

    return false;
}

std::vector<uint32_t> GreedyDfsContext::collect_candidates(
    uint32_t from,
    const std::unordered_set<uint32_t>& on_path_seg,
    const std::unordered_set<uint32_t>& covered,
    const std::unordered_map<uint64_t, uint32_t>& edge_try_count
) const {
    std::vector<uint32_t> cand;
    const auto& arcs = graph_.getArcsFromVertex(from);
    cand.reserve(arcs.size());

    log_.greedy_expand(from, on_path_seg.size() > 0 ? on_path_seg.size() - 1 : 0, arcs.size());

    std::unordered_set<uint32_t> seen_target;
    seen_target.reserve(arcs.size() * 2 + 1);

    for (const GfaArc* a : arcs) {
        if (!a) {
            continue;
        }

        if (a->get_del()) {
            continue;
        }

        if (skip_comp_ && a->get_comp()) {
            continue;
        }

        const uint32_t w = a->get_target_vertex_id();

        if (!seen_target.insert(w).second) {
            log_.greedy_skip(from, w, "duplicate-target");
            continue;
        }

        const std::string reason = reject_reason(from, w, on_path_seg);
        if (!reason.empty()) {
            log_.greedy_skip(from, w, reason);
            continue;
        }

        const uint32_t sid = NodeHandle::get_segment_id(w);
        const auto it_try = edge_try_count.find(edge_key(from, w));
        const uint32_t tries = it_try == edge_try_count.end() ? 0u : it_try->second;

        cand.push_back(w);
        log_.greedy_candidate(
            w,
            graph_.getNodeLength(sid),
            covered.find(w) != covered.end(),
            tries
        );
    }

    sort_candidates(from, cand, covered, edge_try_count);
    log_.greedy_sorted_candidates(cand, from, covered, edge_try_count);

    return cand;
}

void GreedyDfsContext::sort_candidates(
    uint32_t from,
    std::vector<uint32_t>& cand,
    const std::unordered_set<uint32_t>& covered,
    const std::unordered_map<uint64_t, uint32_t>& edge_try_count
) const {
    std::sort(cand.begin(), cand.end(),
        [&](uint32_t a, uint32_t b) {
            const bool sa = graph_.intersect_vertex_sample(from, a);
            const bool sb = graph_.intersect_vertex_sample(from, b);
            if (sa != sb) return sa > sb;

            const auto ita = edge_try_count.find(edge_key(from, a));
            const auto itb = edge_try_count.find(edge_key(from, b));

            const uint32_t ta = ita == edge_try_count.end() ? 0u : ita->second;
            const uint32_t tb = itb == edge_try_count.end() ? 0u : itb->second;
            if (ta != tb) return ta < tb;

            const bool ua = covered.find(a) == covered.end();
            const bool ub = covered.find(b) == covered.end();
            if (ua != ub) return ua > ub;

            const uint32_t oa = graph_.get_edge_ow(from, a);
            const uint32_t ob = graph_.get_edge_ow(from, b);
            if (oa != ob) return oa > ob;

            const uint32_t la = graph_.getNodeLength(NodeHandle::get_segment_id(a));
            const uint32_t lb = graph_.getNodeLength(NodeHandle::get_segment_id(b));
            if (la != lb) return la > lb;

            return a < b;
        }
    );
}
/* ================================================================================================================
 *                                        GREEDY DFS CONTEXT END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                            GFA WALKER START
 * ================================================================================================================ */
GfaWalker::GfaWalker(const GfaGraph& graph)
    : graph_(graph)
{}

uint64_t GfaWalker::hash_vec_(std::vector<uint32_t> vec, bool use_direction) {
    if (!use_direction) {
        for (uint32_t& v : vec) v >>= 1;
    }

    uint64_t h = 14695981039346656037ULL;
    for (uint32_t v : vec) {
        h ^= static_cast<uint64_t>(v);
        h *= 1099511628211ULL;
    }
    return h;
}

std::vector<std::vector<uint32_t>> GfaWalker::enumerate_paths_greedy_DFS(
    uint32_t src,
    uint32_t sink,
    const std::unordered_set<uint32_t>& region_set,
    uint32_t max_depth,
    uint32_t max_paths,
    bool skip_comp,
    bool& hit_limits,
    uint64_t DFS_guard,
    uint32_t stall_round_limit
) const {
    std::vector<std::vector<uint32_t>> out;
    hit_limits = false;

    if (src == sink) {
        out.push_back({src});
        return out;
    }

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (V == 0 || src >= V || sink >= V) return out;

    GfaWalkerLogger log(graph_);
    GreedyDfsContext ctx(graph_, src, sink, region_set, skip_comp, log);

    log.greedy_header(
        src,
        sink,
        max_depth,
        max_paths,
        DFS_guard,
        stall_round_limit,
        ctx.use_topo_filter(),
        ctx.src_rank(),
        ctx.sink_rank(),
        region_set.size()
    );

    std::unordered_set<uint32_t> covered;
    covered.reserve(region_set.size() + 8);
    covered.insert(src);

    std::unordered_set<uint64_t> seen_path;
    seen_path.reserve(64);

    std::unordered_map<uint64_t, uint32_t> edge_try_count;
    edge_try_count.reserve(256);

    GreedyBranchTracker branches(stall_round_limit);

    uint64_t dfs_states = 0;
    uint32_t round_id = 0;

    while (true) {
        if (max_paths > 0 && out.size() >= max_paths) break;

        log.greedy_round(round_id, out.size(), covered.size(), branches.last_stall(), dfs_states);
        ++round_id;

        std::vector<uint32_t> path;
        path.reserve(max_depth + 1);
        path.push_back(src);

        std::unordered_set<uint32_t> on_path_seg;
        on_path_seg.reserve(max_depth + 1);
        on_path_seg.insert(NodeHandle::get_segment_id(src));

        std::vector<GreedyFrame> stack;
        stack.reserve(max_depth + 1);
        stack.push_back(GreedyFrame{src, {}, 0});

        bool reached_sink = false;
        uint64_t branch = GreedyBranchTracker::NO_BRANCH;

        log.greedy_start(src);

        auto pop_node = [&]() {
            const uint32_t v = stack.back().v;
            stack.pop_back();

            if (path.size() > 1) {
                on_path_seg.erase(NodeHandle::get_segment_id(path.back()));
                path.pop_back();
            }

            (void)v;
        };

        while (!stack.empty()) {
            if (++dfs_states > DFS_guard) {
                hit_limits = true;
                log.greedy_discard(path, "exceed DFS_guard", branches.last_stall() + 1);
                return out;
            }

            GreedyFrame& fr = stack.back();
            const uint32_t v = fr.v;
            const uint32_t depth = static_cast<uint32_t>(path.size() - 1);

            if (v == sink) {
                reached_sink = true;
                break;
            }

            if (depth >= max_depth) {
                log.greedy_backtrack(v, "reach max_depth");
                pop_node();
                continue;
            }

            if (fr.cand.empty() && fr.next == 0) {
                fr.cand = ctx.collect_candidates(v, on_path_seg, covered, edge_try_count);

                if (v == src) {
                    branches.remove_closed(src, fr.cand);
                }

                if (fr.cand.empty()) {
                    log.greedy_backtrack(v, "dead end");
                    pop_node();
                    continue;
                }
            }

            if (v == src && branches.has_branch(branch)) {
                log.greedy_backtrack(v, "branch exhausted");
                stack.clear();
                continue;
            }

            if (fr.next >= fr.cand.size()) {
                log.greedy_backtrack(v, "all candidates tried");
                pop_node();
                continue;
            }

            const uint32_t w = fr.cand[fr.next++];

            if (!ctx.valid_next(v, w, on_path_seg)) {
                log.greedy_skip(v, w, "stale-candidate");
                continue;
            }

            log.greedy_choose(w, fr.next, fr.cand.size());

            if (v == src) {
                branch = branches.make_branch(src, w);
            }

            ++edge_try_count[edge_key(v, w)];

            path.push_back(w);
            on_path_seg.insert(NodeHandle::get_segment_id(w));
            stack.push_back(GreedyFrame{w, {}, 0});
        }

        const bool has_new_node = ctx.path_has_new_node(path, covered);

        for (uint32_t v : path) {
            covered.insert(v);
        }

        if (!branches.has_branch(branch)) {
            break;
        }

        if (!reached_sink) {
            log.greedy_discard(path, "did not reach sink", branches.stall(branch));
            continue;
        }

        const uint64_t hv = hash_vec_(path);
        if (!seen_path.insert(hv).second) {
            log.greedy_discard(path, "duplicate path", branches.stall(branch));
            continue;
        }

        out.push_back(path);
        log.greedy_record(path, out.size(), has_new_node);

        if (has_new_node) {
            branches.reset(branch);
        } else {
            branches.stall(branch);
        }
    }

    log.greedy_summary(out.size(), covered.size(), dfs_states, round_id);
    return out;
}

std::vector<std::vector<uint32_t>> GfaWalker::open_walk(
    uint32_t src,
    const std::unordered_set<uint32_t>& region_set,
    uint32_t max_depth,
    bool skip_comp,
    uint64_t DFS_guard,
    uint64_t walk_bp,
    const std::unordered_set<uint32_t>* blocked_seg
) const {
    std::vector<std::vector<uint32_t>> out;

    const uint64_t V = graph_.getNumNodes() * 2;
    if (V == 0 || src >= V || walk_bp == 0) return out;

    const uint32_t depth_limit = max_depth == 0 ? static_cast<uint32_t>(std::min<uint64_t>(V, UINT32_MAX)) : max_depth;

    GfaWalkerLogger log(graph_);
    log.open_header(src, walk_bp, depth_limit, DFS_guard);

    std::vector<uint32_t> path;
    path.reserve(std::min<size_t>(depth_limit + uint64_t(1), 4096));
    path.push_back(src);

    std::unordered_set<uint32_t> used_seg;
    used_seg.reserve(std::min<size_t>(depth_limit + uint64_t(1), 4096) + (blocked_seg ? blocked_seg->size() : 0));

    if (blocked_seg) {
        used_seg.insert(blocked_seg->begin(), blocked_seg->end());
    }
    used_seg.insert(NodeHandle::get_segment_id(src));

    uint32_t cur = src;
    uint32_t steps = 0;
    uint64_t total_bp = 0;
    uint64_t examined_arcs = 0;

    auto has_continuation = [&](uint32_t from, bool& guard_hit) {
        const uint32_t from_sid = NodeHandle::get_segment_id(from);

        for (const GfaArc* arc : graph_.getArcsFromVertex(from)) {
            if (DFS_guard > 0 && examined_arcs >= DFS_guard) {
                guard_hit = true;
                return false;
            }
            ++examined_arcs;

            if (!arc || arc->get_del() || (skip_comp && arc->get_comp())) continue;

            const uint32_t next = arc->get_target_vertex_id();
            if (next == from || next >= V) continue;
            if (!region_set.empty() && !region_set.count(next)) continue;

            const uint32_t sid = NodeHandle::get_segment_id(next);
            if (sid == from_sid || graph_.getNodeDeleted(sid) || used_seg.count(sid)) continue;

            return true;
        }

        return false;
    };

    while (total_bp < walk_bp && steps < depth_limit) {
        uint32_t best = UINT32_MAX;
        uint32_t best_len = 0;
        uint32_t best_ow = 0;
        bool best_same_sample = false;
        bool best_dead_end = true;
        bool guard_hit = false;

        for (const GfaArc* arc : graph_.getArcsFromVertex(cur)) {
            if (DFS_guard > 0 && examined_arcs >= DFS_guard) {
                guard_hit = true;
                break;
            }
            ++examined_arcs;

            if (!arc || arc->get_del()) continue;
            if (skip_comp && arc->get_comp()) continue;

            const uint32_t next = arc->get_target_vertex_id();
            if (next == cur || next >= V) continue;

            if (!region_set.empty() && !region_set.count(next)) continue;

            const uint32_t sid = NodeHandle::get_segment_id(next);
            if (graph_.getNodeDeleted(sid) || used_seg.count(sid)) continue;

            const uint32_t len = graph_.getNodeLength(sid);
            const bool dead_end = !has_continuation(next, guard_hit);
            if (guard_hit) break;

            if (best != UINT32_MAX) {
                if (dead_end && !best_dead_end) continue;
                if (dead_end == best_dead_end && len < best_len) continue;
            }

            const uint32_t ow = arc->ow > 0 && arc->ow != INT32_MAX ? static_cast<uint32_t>(arc->ow) : 0u;
            const bool same_sample = graph_.intersect_vertex_sample(cur, next);

            const bool better =
                best == UINT32_MAX ||
                (best_dead_end && !dead_end) ||
                (dead_end == best_dead_end && len > best_len) ||
                (dead_end == best_dead_end && len == best_len && same_sample > best_same_sample) ||
                (dead_end == best_dead_end && len == best_len && same_sample == best_same_sample && ow > best_ow) ||
                (dead_end == best_dead_end && len == best_len && same_sample == best_same_sample && ow == best_ow && next < best);

            if (better) {
                best = next;
                best_len = len;
                best_ow = ow;
                best_same_sample = same_sample;
                best_dead_end = dead_end;
            }
        }

        if (guard_hit || best == UINT32_MAX) break;

        log.open_step(cur, best);

        total_bp += best_len > best_ow ? uint64_t(best_len - best_ow) : 0;

        cur = best;
        path.push_back(best);
        used_seg.insert(NodeHandle::get_segment_id(best));
        ++steps;
    }

    log.open_summary(path, total_bp);
    out.push_back(std::move(path));

    return out;
}
/* ================================================================================================================
 *                                            GFA WALKER END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                       SOURCE-AWARE WALKER START
 * ================================================================================================================ */
namespace {

struct SourceWalkCandidate {
    uint32_t vertex{UINT32_MAX}, distance{UINT32_MAX};
    uint64_t added_bp{0};
    double source_weight{0.0};
    bool blocked{false}, preferred{false}, continuous{false}, supported{false};
    std::vector<uint32_t> sources;
};

struct SourceWalkStep {
    uint64_t added_bp{0};
    uint32_t previous_source{UINT32_MAX};
    std::vector<uint32_t> sources;
    double source_share{0.0};
};

struct SourceWalkFrame {
    uint32_t vertex{UINT32_MAX};
    std::vector<SourceWalkCandidate> candidates;
    size_t next{0};
    SourceWalkStep step;
};

bool source_in_scope(uint32_t source, const GfaSourceWalker::Options& options) {
    return !options.source_scope || options.source_scope->count(source);
}

std::vector<uint32_t> node_sources(
    const GfaGraph& graph,
    uint32_t vertex,
    const GfaSourceWalker::Options& options
) {
    std::vector<uint32_t> sources;
    for (uint32_t source : graph.getNodeSampleIds(Vertex::get_segment_id(vertex))) {
        if (source < graph.getSampleNames().size() && source_in_scope(source, options)) sources.push_back(source);
    }
    return sources;
}

std::vector<uint32_t> transition_sources(
    const GfaGraph& graph,
    uint32_t from,
    uint32_t to,
    const GfaSourceWalker::Options& options
) {
    if (!options.transition_sources) return node_sources(graph, to, options);

    std::vector<uint32_t> sources;
    const uint64_t key = edge_key(from, to);
    const auto found = options.transition_sources->find(key);
    if (found == options.transition_sources->end()) return sources;

    sources.reserve(found->second.size());
    for (uint32_t source : found->second) {
        if (source < graph.getSampleNames().size() && source_in_scope(source, options)) sources.push_back(source);
    }
    return sources;
}

uint32_t choose_source(
    const std::vector<uint32_t>& sources,
    const std::vector<double>& source_weight,
    uint32_t preferred
) {
    uint32_t best = UINT32_MAX;
    for (uint32_t source : sources) {
        if (source >= source_weight.size()) continue;
        if (best == UINT32_MAX || source_weight[source] > source_weight[best] ||
            (source_weight[source] == source_weight[best] && source == preferred && best != preferred) ||
            (source_weight[source] == source_weight[best] && best != preferred && source < best)) {
            best = source;
        }
    }
    return best;
}

}  // namespace

GfaSourceWalker::GfaSourceWalker(const GfaGraph& graph)
    : graph_(graph)
{}

GfaSourceWalker::Result GfaSourceWalker::walk(
    uint32_t source,
    uint32_t sink,
    const Options& options
) const {
    Result result;
    const uint64_t vertex_count = graph_.getNumNodes() * 2;
    if (source >= vertex_count || sink >= vertex_count ||
        graph_.getNodeDeleted(Vertex::get_segment_id(source)) ||
        graph_.getNodeDeleted(Vertex::get_segment_id(sink))) return result;

    if (source == sink) {
        result.vertices.push_back(source);
        return result;
    }

    const uint32_t depth_limit = options.max_depth == 0
        ? static_cast<uint32_t>(std::min<uint64_t>(vertex_count, UINT32_MAX))
        : options.max_depth;
    const uint64_t state_limit = options.max_states == 0 ? UINT64_MAX : options.max_states;

    const bool have_topology = graph_.get_vertex_topological_index().size() >= vertex_count;
    const uint32_t source_rank = have_topology ? graph_.vertex_topo_rank(Vertex(source)) : UINT32_MAX;
    const uint32_t sink_rank = have_topology ? graph_.vertex_topo_rank(Vertex(sink)) : UINT32_MAX;
    const bool use_topology = source_rank != UINT32_MAX && sink_rank != UINT32_MAX;
    if (use_topology && source_rank > sink_rank) return result;

    // Reverse reachability removes branches that cannot reach the requested sink.
    std::unordered_map<uint32_t, uint32_t> sink_distance;
    sink_distance.reserve(1024);
    std::queue<uint32_t> pending;
    sink_distance.emplace(sink, 0);
    pending.push(sink);

    while (!pending.empty()) {
        const uint32_t vertex = pending.front();
        pending.pop();
        const uint32_t distance = sink_distance[vertex];
        if (distance >= depth_limit) continue;

        std::vector<uint32_t> predecessors;
        // In a bidirected GFA, v^1 -> u^1 is the reverse-complement arc of
        // u -> v.  Reading that adjacency directly avoids rebuilding all
        // incoming arcs through a second neighbour scan.
        for (const GfaArc* arc : graph_.getArcsFromVertex(vertex ^ 1u)) {
            if (++result.states > state_limit) {
                result.truncated = true;
                return result;
            }
            if (!arc || arc->get_del() || (options.skip_comp && arc->get_comp())) continue;

            const uint32_t previous = arc->get_target_vertex_id() ^ 1u;
            if (previous >= vertex_count || graph_.getNodeDeleted(Vertex::get_segment_id(previous))) continue;
            if (use_topology) {
                const uint32_t rank = graph_.vertex_topo_rank(Vertex(previous));
                if (rank == UINT32_MAX || rank < source_rank || rank > sink_rank) continue;
            }
            predecessors.push_back(previous);
        }
        std::sort(predecessors.begin(), predecessors.end());
        predecessors.erase(std::unique(predecessors.begin(), predecessors.end()), predecessors.end());
        for (uint32_t previous : predecessors) {
            if (sink_distance.emplace(previous, distance + 1).second) pending.push(previous);
        }
    }
    if (!sink_distance.count(source)) return result;

    // The DFS is deterministic. Dynamic SN support and soft blocking only order branches;
    // every reachable branch remains available for backtracking.
    std::vector<double> source_weight(graph_.getSampleNames().size(), 0.0);
    std::vector<uint32_t> path{source};
    std::unordered_set<uint32_t> path_segments{Vertex::get_segment_id(source)};
    std::vector<SourceWalkFrame> stack{{source, {}, 0, {}}};
    uint64_t path_bp = 0;
    uint32_t active_source = UINT32_MAX;

    const auto pop_frame = [&]() {
        SourceWalkStep& step = stack.back().step;
        for (uint32_t sample : step.sources) source_weight[sample] -= step.source_share;
        path_bp -= step.added_bp;
        active_source = step.previous_source;
        path_segments.erase(Vertex::get_segment_id(path.back()));
        path.pop_back();
        stack.pop_back();
    };

    while (!stack.empty()) {
        SourceWalkFrame& frame = stack.back();
        const uint32_t depth = static_cast<uint32_t>(path.size() - 1);

        if (frame.vertex == sink) {
            result.vertices = path;
            return result;
        }

        if (depth >= depth_limit) {
            if (stack.size() == 1) break;
            pop_frame();
            continue;
        }

        if (frame.candidates.empty() && frame.next == 0) {
            std::vector<uint32_t> targets;
            for (const GfaArc* arc : graph_.getArcsFromVertex(frame.vertex)) {
                if (++result.states > state_limit) {
                    result.truncated = true;
                    result.vertices.clear();
                    return result;
                }
                if (!arc || arc->get_del() || (options.skip_comp && arc->get_comp())) continue;
                const uint32_t next = arc->get_target_vertex_id();
                if (next >= vertex_count || next == frame.vertex || !sink_distance.count(next)) continue;
                if (depth + 1 + sink_distance[next] > depth_limit) continue;

                const uint32_t segment = Vertex::get_segment_id(next);
                if (graph_.getNodeDeleted(segment) || path_segments.count(segment)) continue;
                targets.push_back(next);
            }
            std::sort(targets.begin(), targets.end());
            targets.erase(std::unique(targets.begin(), targets.end()), targets.end());

            frame.candidates.reserve(targets.size());
            for (uint32_t next : targets) {
                const uint32_t segment = Vertex::get_segment_id(next);
                const uint64_t added = next == sink ? 0 : graph_.getNodeLength(segment);
                if (path_bp > options.max_bp || added > options.max_bp - path_bp) continue;

                SourceWalkCandidate candidate;
                candidate.vertex = next;
                candidate.distance = sink_distance[next];
                candidate.added_bp = added;
                candidate.blocked = next != sink && options.soft_blocked &&
                    segment < options.soft_blocked->size() &&
                    ((*options.soft_blocked)[segment] & options.soft_block_mask);

                candidate.sources = transition_sources(graph_, frame.vertex, next, options);
                candidate.supported = !candidate.sources.empty();
                for (uint32_t sample : candidate.sources) {
                    if (sample >= source_weight.size()) continue;
                    candidate.source_weight = std::max(candidate.source_weight, source_weight[sample]);
                    candidate.preferred |= sample == options.preferred_source;
                    candidate.continuous |= sample == active_source;
                }
                frame.candidates.push_back(candidate);
            }

            std::sort(frame.candidates.begin(), frame.candidates.end(), [](const SourceWalkCandidate& a, const SourceWalkCandidate& b) {
                if (a.blocked != b.blocked) return a.blocked < b.blocked;
                if (a.source_weight != b.source_weight) return a.source_weight > b.source_weight;
                if (a.continuous != b.continuous) return a.continuous > b.continuous;
                if (a.preferred != b.preferred) return a.preferred > b.preferred;
                if (a.supported != b.supported) return a.supported > b.supported;
                if (a.distance != b.distance) return a.distance < b.distance;
                return a.vertex < b.vertex;
            });
        }

        if (frame.next >= frame.candidates.size()) {
            if (stack.size() == 1) break;
            pop_frame();
            continue;
        }
        SourceWalkCandidate candidate = std::move(frame.candidates[frame.next++]);
        SourceWalkStep step;
        step.added_bp = candidate.added_bp;
        step.previous_source = active_source;

        if (candidate.vertex != sink) {
            step.sources = std::move(candidate.sources);
            if (!step.sources.empty()) {
                step.source_share = static_cast<double>(candidate.added_bp) / step.sources.size();
                for (uint32_t sample : step.sources) source_weight[sample] += step.source_share;

                if (!std::binary_search(step.sources.begin(), step.sources.end(), active_source)) {
                    const uint32_t next_source = choose_source(step.sources, source_weight, options.preferred_source);
                    active_source = next_source;
                }
            }
        }

        path_bp += step.added_bp;
        path.push_back(candidate.vertex);
        path_segments.insert(Vertex::get_segment_id(candidate.vertex));
        stack.push_back({candidate.vertex, {}, 0, std::move(step)});
    }

    return result;
}
/* ================================================================================================================
 *                                        SOURCE-AWARE WALKER END
 * ================================================================================================================ */
