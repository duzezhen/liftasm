#include "../include/gfa_walker.hpp"
#include "../include/gfa_walker_logger.hpp"

#include <algorithm>
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
 *                                         TOPO SINK FINDER START
 * ================================================================================================================ */
TopoSinkFinder::TopoSinkFinder(
    const GfaGraph& graph,
    const std::unordered_set<uint32_t>& region_set,
    const std::unordered_set<uint32_t>& used_seg,
    const std::unordered_set<uint32_t>& bad_sinks,
    bool skip_comp,
    uint32_t max_depth,
    uint64_t dfs_guard
)
    : graph_(graph),
      region_set_(region_set),
      used_seg_(used_seg),
      bad_sinks_(bad_sinks),
      skip_comp_(skip_comp),
      max_depth_(max_depth),
      dfs_guard_(dfs_guard)
{}

bool TopoSinkFinder::valid_vertex(uint32_t v) const {
    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (v >= V) return false;

    const uint32_t sid = NodeHandle::get_segment_id(v);
    if (graph_.getNodeDeleted(sid)) return false;

    if (!region_set_.empty() && region_set_.find(v) == region_set_.end()) {
        return false;
    }

    return true;
}

std::vector<uint32_t> TopoSinkFinder::collect_candidates(uint32_t from, std::unordered_set<uint32_t>& seen) const {
    std::vector<uint32_t> cand;
    const auto& arcs = graph_.getArcsFromVertex(from);
    cand.reserve(arcs.size());

    for (const GfaArc* a : arcs) {
        if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;

        const uint32_t w = a->get_target_vertex_id();
        if (!valid_vertex(w)) continue;
        if (!seen.insert(w).second) continue;
        if (bad_sinks_.find(w) != bad_sinks_.end()) continue;

        const uint32_t sid = NodeHandle::get_segment_id(w);
        if (used_seg_.find(sid) != used_seg_.end()) continue;

        cand.push_back(w);
    }

    return cand;
}

void TopoSinkFinder::sort_candidates(uint32_t from, uint32_t from_rank, std::vector<uint32_t>& cand) const {
    std::sort(cand.begin(), cand.end(),
        [&](uint32_t a, uint32_t b) {
            const bool sa = graph_.intersect_vertex_sample(from, a);
            const bool sb = graph_.intersect_vertex_sample(from, b);
            if (sa != sb) return sa > sb;

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

std::vector<uint32_t> TopoSinkFinder::find_path(uint32_t from, uint64_t target_bp, bool& hit_limits) const {
    hit_limits = false;

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (from >= V || graph_.get_vertex_topological_index().size() < V) return {};
    if (!valid_vertex(from)) return {};

    const uint32_t from_rank = graph_.vertex_topo_rank(Vertex(from));
    if (from_rank == UINT32_MAX) return {};

    struct State {
        uint32_t v{UINT32_MAX};
        uint32_t depth{0};
        uint64_t bp{0};
    };

    std::queue<State> q;
    q.push(State{from, 0, 0});

    std::unordered_set<uint32_t> seen;
    seen.reserve(256);
    seen.insert(from);

    std::unordered_map<uint32_t, uint32_t> parent;
    parent.reserve(256);

    auto build_path = [&](uint32_t sink) {
        std::vector<uint32_t> path;
        for (uint32_t v = sink;;) {
            path.push_back(v);
            if (v == from) break;
            const auto it = parent.find(v);
            if (it == parent.end()) return std::vector<uint32_t>{};
            v = it->second;
        }
        std::reverse(path.begin(), path.end());
        return path;
    };

    uint64_t states = 0;
    uint32_t best_sink = UINT32_MAX;
    uint64_t best_bp = 0;

    while (!q.empty()) {
        if (++states > dfs_guard_) {
            if (best_sink != UINT32_MAX) return build_path(best_sink);
            hit_limits = true;
            return {};
        }

        const State cur = q.front();
        q.pop();

        if (cur.depth >= max_depth_) continue;

        std::vector<uint32_t> cand = collect_candidates(cur.v, seen);
        sort_candidates(cur.v, from_rank, cand);

        for (uint32_t w : cand) {
            parent.emplace(w, cur.v);
            const uint32_t sid = NodeHandle::get_segment_id(w);
            const uint32_t len = graph_.getNodeLength(sid);
            const uint32_t ow = graph_.get_edge_ow(cur.v, w);
            const uint64_t next_bp = cur.bp + (len > ow ? uint64_t(len - ow) : 0);
            const uint32_t wrank = graph_.vertex_topo_rank(Vertex(w));

            if (wrank != UINT32_MAX && wrank > from_rank && !graph_.vertex_in_cycle(Vertex(w))) {
                if (best_sink == UINT32_MAX || next_bp > best_bp) {
                    best_sink = w;
                    best_bp = next_bp;
                }
                if (next_bp >= target_bp) return build_path(w);
            }

            q.push(State{w, cur.depth + 1, next_bp});
        }
    }

    return best_sink == UINT32_MAX ? std::vector<uint32_t>{} : build_path(best_sink);
}
/* ================================================================================================================
 *                                         TOPO SINK FINDER END
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

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    const std::vector<GfaTopoIndex>& topo = graph_.get_vertex_topological_index();

    if (V == 0 || src >= V || walk_bp == 0) return out;
    if (topo.size() < V) return out;

    const uint32_t depth_limit = max_depth == 0 ? V : max_depth;

    GfaWalkerLogger log(graph_);
    log.open_header(src, walk_bp, depth_limit, DFS_guard);

    auto added_bp = [&](uint32_t from, uint32_t to) -> uint64_t {
        const uint32_t sid = NodeHandle::get_segment_id(to);
        const uint32_t len = graph_.getNodeLength(sid);
        const uint32_t ow = graph_.get_edge_ow(from, to);

        return len > ow ? uint64_t(len - ow) : 0;
    };

    std::unordered_set<uint32_t> used_seg;
    used_seg.reserve(depth_limit + 8);
    used_seg.insert(NodeHandle::get_segment_id(src));

    if (blocked_seg) {
        used_seg.insert(blocked_seg->begin(), blocked_seg->end());
    }

    std::unordered_set<uint32_t> bad_sinks;
    bad_sinks.reserve(32);

    std::vector<uint32_t> full_path;
    full_path.reserve(depth_limit + 1);
    full_path.push_back(src);

    uint32_t cur = src;
    uint64_t total_bp = 0;
    while (total_bp < walk_bp && full_path.size() - 1 < depth_limit) {
        const uint32_t remaining_depth = depth_limit - static_cast<uint32_t>(full_path.size() - 1);
        TopoSinkFinder finder(
            graph_,
            region_set,
            used_seg,
            bad_sinks,
            skip_comp,
            remaining_depth,
            DFS_guard
        );

        bool sink_hit_limits = false;
        std::vector<uint32_t> piece = finder.find_path(cur, walk_bp - total_bp, sink_hit_limits);

        if (sink_hit_limits || piece.size() < 2) {
            break;
        }

        const uint32_t sink = piece.back();
        log.open_sink(cur, sink);

        bool clean = true;
        std::unordered_set<uint32_t> piece_seg;
        piece_seg.reserve(piece.size() + 8);

        for (size_t k = 1; k < piece.size(); ++k) {
            const uint32_t sid = NodeHandle::get_segment_id(piece[k]);

            if (used_seg.find(sid) != used_seg.end()) {
                clean = false;
                break;
            }

            if (!piece_seg.insert(sid).second) {
                clean = false;
                break;
            }
        }

        if (!clean) {
            bad_sinks.insert(sink);
            log.open_reject(sink, "cycle or used segment");
            continue;
        }

        for (size_t k = 1; k < piece.size(); ++k) {
            total_bp += added_bp(full_path.back(), piece[k]);
            full_path.push_back(piece[k]);
            used_seg.insert(NodeHandle::get_segment_id(piece[k]));
        }

        cur = full_path.back();
    }

    log.open_summary(full_path, total_bp);

    if (!full_path.empty()) {
        out.push_back(std::move(full_path));
    }

    return out;
}
/* ================================================================================================================
 *                                            GFA WALKER END
 * ================================================================================================================ */
