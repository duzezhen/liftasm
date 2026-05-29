#pragma once
#include <vector>
#include <string>
#include <queue>
#include <stack>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "gfa_parser.hpp"
#include "gfa_parser_types.hpp"

enum class WalkDir { From, To };

struct DFSFrame {
    uint32_t v;          // current vertex id (orientation‑aware)
    uint64_t nextIdx;    // next outgoing arc index to visit
    uint64_t depthLeft;  // remaining depth budget
};

/* ============================================================
 *  GfaNodeSelector: helper for selecting the next node
 * ============================================================ */
class GfaNodeSelector {
public:
    explicit GfaNodeSelector(const GfaGraph& graph)
        : graph_(graph)
    {}

    struct SegNameInfo {
        // e.g. "utg000545l:177754-177755"
        std::string base;
        uint64_t    beg{0};
        uint64_t    end{0};
        bool        ok{false};
    };

    // Statistics for simple-base in the current path, maintained by DFS
    struct NodeSimpleStats {
        std::unordered_map<std::string, uint32_t> base_cnt;      // base -> count
        std::unordered_map<std::string, uint32_t> base_last_end; // base -> last end coordinate
    };

    struct ScoreContext {
        const std::unordered_set<uint32_t>*           region_nodes{nullptr};
        const std::unordered_set<uint32_t>*           covered_nodes{nullptr};
        const std::unordered_map<uint32_t, uint32_t>* visit_times{nullptr};
        const std::unordered_map<uint32_t, uint32_t>* fail_times{nullptr};
    };

    // Parse vertex name
    SegNameInfo get_name_info(uint32_t v) const {
        uint32_t sid = NodeHandle::get_segment_id(v);
        const GfaNode* nd = graph_.getNode(sid);
        if (!nd) return {};
        return parse_seg_name_(nd->name);
    }

    // Check if a vertex is simple (indegree=1 and outdegree=1)
    bool is_simple(uint32_t v) const {
        const auto& outs = graph_.getArcsFromVertex(v);
        const auto& ins = graph_.getArcsToVertex(v);
        return (outs.size() == 1 && ins.size() == 1);
    }

    // Update simple-base statistics for the current path when entering a new vertex
    void add_stats_for_vertex(uint32_t v, NodeSimpleStats& stats, bool skip_simple = true) const {
        if (skip_simple && !is_simple(v)) return;
        SegNameInfo info = get_name_info(v);
        if (!info.ok) return;
        uint32_t &cnt = stats.base_cnt[info.base];
        ++cnt;
        stats.base_last_end[info.base] = info.end;
    }
    NodeSimpleStats build_stats_for_path(const std::vector<uint32_t>& path) const {
        NodeSimpleStats st;
        for (uint32_t v : path) {
            add_stats_for_vertex(v, st, /*skip_simple=*/false);
        }
        return st;
    }

    int64_t compute_weight(
        uint32_t from_v,
        uint32_t w,
        const NodeSimpleStats& stats,
        const ScoreContext* sc = nullptr
    ) const {
        const int IN_REGION_NEW   =  +5;  // Node in region but not yet covered
        const int IN_REGION_USED  =  -2;  // Node in region and already covered (non-simple)
        const int SIMPLE_BONUS    =  +2;  // Simple node
        const int BASE_FREQ_UNIT  =  +1;  // Reward for base appearing once in current path
        const int DIST_STEP_PEN   =  -1;  // Penalty per distance bin (farther is worse)
        const int DEG_STEP_PEN    =  -1;  // Degree penalty
        const int LEN_BIG_BONUS   =  +2;  // Very long segment
        const int LEN_MED_BONUS   =  +1;  // Medium length

        const int VISIT_PEN_UNIT  =  -1;  // Penalty for one successful visit (light)
        const int FAIL_PEN_UNIT   =  -2;  // Penalty for one failed visit (heavy)

        // Distance bin size: 100bp
        const uint64_t DIST_BIN   = 100;

        int64_t score = 0;

        // --- region / covered / visit ---
        bool in_region = false;
        bool used      = false;
        uint32_t suc   = 0;
        uint32_t fail  = 0;

        if (sc) {
            if (sc->region_nodes) {
                in_region = (sc->region_nodes->find(w) != sc->region_nodes->end());
            }
            if (sc->covered_nodes) {
                used = (sc->covered_nodes->find(w) != sc->covered_nodes->end());
            }
            if (sc->visit_times) {
                auto it = sc->visit_times->find(w);
                if (it != sc->visit_times->end()) suc = it->second;
            }
            if (sc->fail_times) {
                auto it2 = sc->fail_times->find(w);
                if (it2 != sc->fail_times->end()) fail = it2->second;
            }
        }

        const bool simple = is_simple(w);

        if (in_region && !used) {
            score += IN_REGION_NEW;
        }
        if (in_region && used && !simple) {
            score += IN_REGION_USED;
        }

        if (simple) {
            score += SIMPLE_BONUS;
        }

        // --- base / coordinate / length ---
        SegNameInfo from_info = get_name_info(from_v);
        SegNameInfo w_info    = get_name_info(w);

        const GfaNode* nd = graph_.getNode(NodeHandle::get_segment_id(w));
        uint32_t len = nd ? nd->length : 0u;

        // Degree
        int indeg = 0, outdeg = 0;
        for (const GfaArc* a : graph_.getArcsToVertex(w)) {
            if (!a || a->get_del() || a->get_comp()) continue;
            ++indeg;
        }
        for (const GfaArc* a : graph_.getArcsFromVertex(w)) {
            if (!a || a->get_del() || a->get_comp()) continue;
            ++outdeg;
        }
        int degree = indeg + outdeg;

        int      base_freq = 0;
        uint64_t distance  = 0;

        if (w_info.ok) {
            auto itc = stats.base_cnt.find(w_info.base);
            if (itc != stats.base_cnt.end()) {
                base_freq = static_cast<int>(itc->second);
            }

            uint64_t last_end = 0;
            bool have_last = false;

            auto it_end = stats.base_last_end.find(w_info.base);
            if (it_end != stats.base_last_end.end()) {
                last_end  = it_end->second;
                have_last = true;
            }

            if (have_last) {
                distance = (w_info.beg > last_end) ? (w_info.beg - last_end) : (last_end - w_info.beg);
            }
        }

        // Frequency of base in current path
        if (base_freq > 0) {
            int capped = std::min(base_freq, 3);
            score += capped * BASE_FREQ_UNIT;
        }

        // Distance penalty: farther distance results in more penalty,
        // 1 bin (~100bp) per point, up to 3 points
        if (distance > 0 && DIST_BIN > 0) {
            uint64_t bin = distance / DIST_BIN;
            int dpen = static_cast<int>(std::min<uint64_t>(bin, 3));
            score += dpen * DIST_STEP_PEN;
        }

        // Degree penalty
        if (degree > 0) {
            int d = std::min(degree, 3);
            score += d * DEG_STEP_PEN;
        }

        // Length: roughly divided into 3 tiers
        if (len >= 10000) {
            score += LEN_BIG_BONUS;
        } else if (len >= 1000) {
            score += LEN_MED_BONUS;
        } // Less than 1kb gets no bonus

        // Visit times penalty: failures weigh more than successes, with the same cap
        if (suc > 0) {
            int s = std::min<uint32_t>(suc, 5);
            score += s * VISIT_PEN_UNIT;
        }
        if (fail > 0) {
            int f = std::min<uint32_t>(fail, 5);
            score += f * FAIL_PEN_UNIT;
        }

        return score;
    }

    // Sort cand_vs by weight (higher score first)
    void sort_candidates(
        uint32_t from_v,
        const std::vector<uint32_t>& path,
        std::vector<uint32_t>& cand_vs,
        const ScoreContext* sc = nullptr
    ) const {
        if (cand_vs.empty()) return;

        NodeSimpleStats stats = build_stats_for_path(path);

        std::sort(cand_vs.begin(), cand_vs.end(),
            [&](uint32_t a, uint32_t b) {
                int64_t wa = compute_weight(from_v, a, stats, sc);
                int64_t wb = compute_weight(from_v, b, stats, sc);
                if (wa != wb) return wa > wb;
                return a < b;
            }
        );
    }

private:
    static SegNameInfo parse_seg_name_(const std::string& name) {
        SegNameInfo info;
        auto pos_colon = name.find(':');
        if (pos_colon == std::string::npos) return info;

        auto pos_dash  = name.find('-', pos_colon + 1);
        if (pos_dash == std::string::npos) return info;

        info.base = name.substr(0, pos_colon);
        try {
            info.beg = std::stoull(name.substr(pos_colon + 1, pos_dash - pos_colon - 1));
            info.end = std::stoull(name.substr(pos_dash + 1));
            info.ok  = true;
        } catch (...) {
            info.ok = false;
        }
        return info;
    }

private:
    const GfaGraph& graph_;
};


class GfaWalker {
public:
    explicit GfaWalker(const GfaGraph& graph)
        : graph_(graph), selector_(graph), V_(static_cast<uint32_t>(graph.getNumNodes() * 2)) {}

    std::vector<PathSequence> dfs(uint32_t start_v, WalkDir dir, uint32_t max_paths = 0, uint32_t max_steps = 0) const {
        return dfs_(start_v, dir, max_paths, max_steps);
    }
    std::vector<PathSequence> bfs(uint32_t start_v, WalkDir dir, uint32_t max_paths = 0, uint32_t max_steps = 0) const {
        return bfs_(start_v, dir, max_paths, max_steps);
    }
    std::vector<std::vector<uint32_t>> enumerate_paths_DFS(
        const uint32_t src, const uint32_t sink, const std::unordered_set<uint32_t>& region_set,
        const uint32_t max_depth, const uint32_t max_paths,
        bool skip_comp, bool& hit_limits, const uint64_t DFS_guard, const uint32_t stall_round_limit
    ) const {
        return enumerate_paths_DFS_(src, sink, region_set, max_depth, max_paths, skip_comp, hit_limits, DFS_guard, stall_round_limit);
    }
    std::vector<std::vector<uint32_t>> enumerate_paths_greedy_DFS(
        const uint32_t src, const uint32_t sink, const std::unordered_set<uint32_t>& region_set,
        const uint32_t max_depth, const uint32_t max_paths,
        bool skip_comp, bool& hit_limits, const uint64_t DFS_guard, const uint32_t stall_round_limit
    ) const {
        return enumerate_paths_greedy_DFS_(src, sink, region_set, max_depth, max_paths, skip_comp, hit_limits, DFS_guard, stall_round_limit);
    }
    

private:
    static uint64_t hash_vec_(std::vector<uint32_t> vec, bool use_direction = true) {
        if (!use_direction) {
            for (auto &x : vec) x >>= 1;
        }

        uint64_t h = 14695981039346656037ULL;  // FNV offset basis
        for (uint32_t v : vec) {
            h ^= static_cast<uint64_t>(v);
            h *= 1099511628211ULL;  // FNV prime
        }
        return h;
    }

    struct SegNameInfo {
        std::string base;   // utg000545l
        uint64_t    beg{0};
        uint64_t    end{0};
        bool        ok{false};
    };
    // parse "utg000545l:177754-177755"
    SegNameInfo parse_seg_name_(const std::string& name) const {
        SegNameInfo info;
        auto pos_colon = name.find(':');
        if (pos_colon == std::string::npos) return info;

        auto pos_dash  = name.find('-', pos_colon + 1);
        if (pos_dash == std::string::npos) return info;

        info.base = name.substr(0, pos_colon);
        try {
            info.beg = std::stoull(name.substr(pos_colon + 1, pos_dash - pos_colon - 1));
            info.end = std::stoull(name.substr(pos_dash + 1));
            info.ok  = true;
        } catch (...) {
            info.ok = false;
        }
        return info;
    }

private:
    static uint32_t append_tail_from_node_(const GfaGraph& graph, std::string& seq, uint32_t v, int32_t ow) {
        const uint32_t seg = NodeHandle::get_segment_id(v);
        const GfaNode* nd  = graph.getNode(seg);
        if (!nd || nd->deleted || nd->length == 0) return 0u;

        const std::string oriented = graph.get_oriented_sequence(Vertex(v));
        if (oriented.empty() || oriented == "*") return 0u;

        uint32_t ow_u = (ow == INT32_MAX) ? 0u : (ow < 0 ? 0u : static_cast<uint32_t>(ow));
        if (ow_u > nd->length) ow_u = nd->length;

        const uint32_t take = (nd->length > ow_u) ? (nd->length - ow_u) : 0u;
        if (take) seq.append(oriented, ow_u, take);
        return take;
    }

    static uint32_t append_left_suffix_of_prefix_(const GfaGraph& graph, std::vector<std::string>& chunks, uint32_t v, uint32_t upto) {
        const uint32_t seg = NodeHandle::get_segment_id(v);
        const GfaNode* nd  = graph.getNode(seg);
        if (!nd || nd->deleted || nd->length == 0) return 0u;

        const std::string oriented = graph.get_oriented_sequence(Vertex(v));
        if (oriented.empty() || oriented == "*") return 0u;

        const uint32_t lim  = std::min<uint32_t>(upto, nd->length);
        if (lim == 0) return 0u;

        chunks.emplace_back(oriented.substr(0, lim));
        return lim;
    }

    static uint32_t instant_lv_to_(const GfaGraph& graph, const GfaArc* a) {
        if (!a) return 0u;
        const uint32_t u     = a->get_source_vertex_id();
        const uint32_t u_seg = NodeHandle::get_segment_id(u);
        const GfaNode* nd    = graph.getNode(u_seg);
        if (!nd || nd->deleted || nd->length == 0) return 0u;
        if (a->ov == INT32_MAX) return nd->length;
        const uint32_t ov_u = (a->ov < 0) ? 0u : static_cast<uint32_t>(a->ov);
        return (nd->length > ov_u) ? (nd->length - ov_u) : 0u;
    }

    static std::string join_chunks_rev_(const std::vector<std::string>& chunks) {
        std::string s;
        size_t sum = 0; for (auto& c : chunks) sum += c.size();
        s.reserve(sum);
        for (auto it = chunks.rbegin(); it != chunks.rend(); ++it) s += *it;
        return s;
    }

    // DFS
    std::vector<PathSequence> dfs_(uint32_t start_v, WalkDir dir, uint32_t max_paths, uint32_t max_steps) const {
        std::vector<PathSequence> results;
        if (start_v >= V_) return results;

        static constexpr int kHardStepCap = 100000;

        auto hit_cap  = [&](size_t n){ return (max_paths != 0) && (n >= max_paths); };
        auto stop_stp = [&](uint32_t stp){ return (max_steps != 0) && (stp >= max_steps); };

        std::string seq; seq.reserve(256);
        std::vector<std::string> chunks; chunks.reserve(8);
        std::vector<uint32_t> path; path.reserve(16);
        std::vector<uint8_t> onp(V_, 0);

        std::function<void(uint32_t,uint32_t,int)> dfs =
        [&](uint32_t v, uint32_t steps, int guard) {
            if (hit_cap(results.size()) || guard > kHardStepCap) return;

            std::vector<const GfaArc*> arcs =
                (dir == WalkDir::From) ? graph_.getArcsFromVertex(v) : graph_.getArcsToVertex(v);

            if (arcs.empty()) {
                if (dir == WalkDir::To) results.push_back(PathSequence{path, join_chunks_rev_(chunks)});
                else                    results.push_back(PathSequence{path, seq});
                return;
            }

            for (const GfaArc* a : arcs) {
                if (hit_cap(results.size())) break;
                if (!a || a->get_del()) continue;

                uint32_t nxt = (dir == WalkDir::From) ? a->get_target_vertex_id() : a->get_source_vertex_id();
                if (nxt >= V_ || onp[nxt]) continue;

                const size_t old_seq_len = seq.size();
                const size_t old_chunks  = chunks.size();

                if (dir == WalkDir::From) {
                    append_tail_from_node_(graph_, seq, nxt, a->ow);
                } else {
                    const uint32_t lv = instant_lv_to_(graph_, a);
                    append_left_suffix_of_prefix_(graph_, chunks, nxt, lv);
                }

                path.push_back(nxt);
                onp[nxt] = 1;

                const bool stop_here = stop_stp(steps + 1);
                std::vector<const GfaArc*> next_arcs;
                if (!stop_here) {
                    next_arcs = (dir == WalkDir::From) ? graph_.getArcsFromVertex(nxt) : graph_.getArcsToVertex(nxt);
                }

                if (stop_here || next_arcs.empty()) {
                    if (dir == WalkDir::To) results.push_back(PathSequence{path, join_chunks_rev_(chunks)});
                    else                    results.push_back(PathSequence{path, seq});
                } else {
                    dfs(nxt, steps + 1, guard + 1);
                }

                onp[nxt] = 0;
                path.pop_back();
                if (dir == WalkDir::From) seq.resize(old_seq_len);
                else                      chunks.resize(old_chunks);
            }
        };

        dfs(start_v, /*steps=*/0, /*guard=*/0);
        return results;
    }

    // BFS
    std::vector<PathSequence> bfs_(uint32_t start_v, WalkDir dir, uint32_t max_paths, uint32_t max_steps) const {
        std::vector<PathSequence> results;
        if (start_v >= V_) return results;

        auto hit_cap  = [&](size_t n){ return (max_paths != 0) && (n >= max_paths); };
        auto stop_stp = [&](uint32_t stp){ return (max_steps != 0) && (stp >= max_steps); };

        if (dir == WalkDir::From) {
            struct StateF {
                uint32_t v;
                uint32_t steps;
                std::string seq;
                std::vector<uint32_t> path;
                std::vector<uint8_t>  onp;
            };
            std::queue<StateF> q;

            for (const GfaArc* a : graph_.getArcsFromVertex(start_v)) {
                if (!a || a->get_del()) continue;
                const uint32_t w = a->get_target_vertex_id();

                StateF s; s.v = w; s.steps = 1;
                s.seq.reserve(256);
                s.path.reserve(8);
                s.onp.assign(V_, 0);

                append_tail_from_node_(graph_, s.seq, w, a->ow);

                s.path.push_back(w);
                s.onp[w] = 1;
                q.push(std::move(s));
            }
            if (q.empty()) return results;

            while (!q.empty() && !hit_cap(results.size())) {
                StateF cur = std::move(q.front()); q.pop();

                const bool stop_here = stop_stp(cur.steps);
                const auto outs = stop_here ? std::vector<const GfaArc*>() : graph_.getArcsFromVertex(cur.v);

                if (stop_here || outs.empty()) {
                    results.push_back(PathSequence{cur.path, cur.seq});
                    continue;
                }

                for (const GfaArc* a : outs) {
                    if (hit_cap(results.size())) break;
                    if (!a || a->get_del()) continue;
                    const uint32_t w = a->get_target_vertex_id();
                    if (w >= cur.onp.size() || cur.onp[w]) continue;

                    StateF nxt = cur;
                    nxt.v = w; nxt.steps = cur.steps + 1;

                    append_tail_from_node_(graph_, nxt.seq, w, a->ow);

                    nxt.path.push_back(w);
                    nxt.onp[w] = 1;
                    q.push(std::move(nxt));
                }
            }
        } else { // WalkDir::To
            struct StateT {
                uint32_t v;
                uint32_t steps;
                std::vector<std::string> chunks;
                std::vector<uint32_t> path;
                std::vector<uint8_t>  onp;
            };
            std::queue<StateT> q;

            for (const GfaArc* a : graph_.getArcsToVertex(start_v)) {
                if (!a || a->get_del()) continue;
                const uint32_t u = a->get_source_vertex_id();

                StateT s; s.v = u; s.steps = 1;
                s.chunks.reserve(8);
                s.path.reserve(8);
                s.onp.assign(V_, 0);

                const uint32_t lv = instant_lv_to_(graph_, a);
                append_left_suffix_of_prefix_(graph_, s.chunks, u, lv);

                s.path.push_back(u);
                s.onp[u] = 1;
                q.push(std::move(s));
            }
            if (q.empty()) return results;

            while (!q.empty() && !hit_cap(results.size())) {
                StateT cur = std::move(q.front()); q.pop();

                const bool stop_here = stop_stp(cur.steps);
                const auto ins = stop_here ? std::vector<const GfaArc*>() : graph_.getArcsToVertex(cur.v);

                if (stop_here || ins.empty()) {
                    results.push_back(PathSequence{cur.path, join_chunks_rev_(cur.chunks)});
                    continue;
                }

                for (const GfaArc* a : ins) {
                    if (hit_cap(results.size())) break;
                    if (!a || a->get_del()) continue;
                    const uint32_t u = a->get_source_vertex_id();
                    if (u >= cur.onp.size() || cur.onp[u]) continue;

                    StateT nxt = cur;
                    nxt.v = u; nxt.steps = cur.steps + 1;

                    const uint32_t lv = instant_lv_to_(graph_, a);
                    append_left_suffix_of_prefix_(graph_, nxt.chunks, u, lv);

                    nxt.path.push_back(u);
                    nxt.onp[u] = 1;
                    q.push(std::move(nxt));
                }
            }
        }

        return results;
    }

    const std::string& get_node_name(uint32_t v) const {
        uint32_t sid = NodeHandle::get_segment_id(v);
        const GfaNode* nd = graph_.getNode(sid);
        static const std::string empty = "";
        return nd ? nd->name : empty;
    }

    /**
     * @brief Enumerate representative DFS paths from src to sink within a local region.
     *
     * Two stages:
     *
     *   1. From 'src', it explores forward and stores the best prefix path for each reachable vertex under the current scoring rule.
     *   2. Using these prefix vertices as pivots, it performs DFS extensions toward 'sink' and records unique full paths.
     *
     * Stops when:
     *
     *   - 'max_paths' unique paths have been collected, or
     *   - no new unique path is found for 'stall_round_limit' consecutive pivot rounds, or
     *   - DFS exploration exceeds 'DFS_guard', or
     *   - no more pivot vertices remain to try.
     *
     * Traversal can be restricted to 'region_set' (except that 'sink' is always allowed), and cycles on the current path are avoided.
     *
     * @param src                source vertex id
     * @param sink               sink vertex id
     * @param region_set         allowed internal vertices for traversal; if empty, no region restriction is applied
     * @param max_depth          maximum DFS depth allowed for a full path
     * @param max_paths          maximum number of unique paths to return; 0 means no explicit path-count cap
     * @param skip_comp          whether to skip arcs marked as complementary
     * @param hit_limits         output flag; set to true if traversal is aborted due to 'DFS_guard'
     * @param DFS_guard          hard upper bound on explored DFS states
     * @param stall_round_limit  stop after this many consecutive pivot rounds without finding a new unique path
     *
     * @return a list of unique src-to-sink paths, each represented as a vertex-id vector
     */
    std::vector<std::vector<uint32_t>> enumerate_paths_DFS_(
        const uint32_t src,
        const uint32_t sink,
        const std::unordered_set<uint32_t>& region_set,
        const uint32_t max_depth,
        const uint32_t max_paths,
        const bool skip_comp,
        bool& hit_limits,
        const uint64_t DFS_guard, 
        const uint32_t stall_round_limit = 2
    ) const
    {
        std::vector<std::vector<uint32_t>> out;
        hit_limits = false;

        if (src == sink) {
            out.push_back({src});
            return out;
        }

        std::unordered_set<uint32_t> covered;
        covered.reserve(region_set.size() + 8);

        const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
        if (V == 0 || src >= V || sink >= V) return out;

        // Only vertices whose topo_rank is between src and sink are allowed.
        uint32_t src_rank = graph_.vertex_topo_rank(Vertex(src));
        uint32_t sink_rank = graph_.vertex_topo_rank(Vertex(sink));
        const bool use_topo_filter = graph_.get_vertex_topological_index().size() >= V && src_rank != UINT32_MAX && sink_rank != UINT32_MAX;

        if (src_rank > sink_rank) std::swap(src_rank, sink_rank);

        // Used to check if a vertex is within the topo window defined by src and sink. If topo ranks are not available, this always returns true.
        auto in_topo_window = [&](uint32_t v) -> bool {
            if (!use_topo_filter) return true;
            const uint32_t r = graph_.vertex_topo_rank(Vertex(v));
            return r != UINT32_MAX && r >= src_rank && r <= sink_rank;
        };

        // Used for debugging
        auto log_indent = [&](int level) -> std::string {
            return std::string(level * 2, ' ');
        };

        auto log_path = [&](const std::vector<uint32_t>& p, int level, const std::string& title) {
            if (!DEBUG_ENABLED) return;
            debug_stream() << log_indent(level) << title << "\n";
            for (uint32_t x : p) {
                debug_stream() << log_indent(level + 1) << get_node_name(x) << "\n";
            }
        };

        if (DEBUG_ENABLED) {
            debug_stream() << "=========================================\n";
            debug_stream() << "Enumerate DFS paths\n";
            debug_stream() << log_indent(1) << "src       : " << get_node_name(src) << "\n";
            debug_stream() << log_indent(1) << "sink      : " << get_node_name(sink) << "\n";
            debug_stream() << log_indent(1) << "max_depth : " << max_depth << "\n";
            debug_stream() << log_indent(1) << "max_paths : " << max_paths << "\n";
            if (use_topo_filter) {
                debug_stream() << log_indent(1) << "topo_win  : " << src_rank << " .. " << sink_rank << "\n";
            }
            debug_stream() << log_indent(1) << "region_set:\n";
            for (uint32_t v : region_set) {
                debug_stream() << log_indent(2) << get_node_name(v) << "\n";
            }
        }

        // ================== Step 1: From src, record best prefix paths (half DFS) ================== //

        struct PrefixInfo {
            std::vector<uint32_t> path;
            uint32_t depth{0};
            int64_t  score{0};
        };

        std::unordered_map<uint32_t, PrefixInfo> best_prefix;
        best_prefix.reserve(256);

        {
            PrefixInfo info;
            info.path  = { src };
            info.depth = 0;
            info.score = 0;
            best_prefix[src] = std::move(info);
        }

        std::queue<uint32_t> q;
        q.push(src);

        uint64_t prefix_states = 0;

        if (DEBUG_ENABLED) {
            debug_stream() << "\n";
            debug_stream() << "  Prefix exploration\n";
            debug_stream() << log_indent(2) << "start from " << get_node_name(src) << "\n";
        }

        while (!q.empty()) {
            uint32_t v = q.front();
            q.pop();

            auto it_v = best_prefix.find(v);
            if (it_v == best_prefix.end()) continue;
            const PrefixInfo &cur = it_v->second;

            if (++prefix_states > DFS_guard) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "stop: exceed DFS_guard = " << DFS_guard << "\n";
                }
                hit_limits = true;
                break;
            }

            if (cur.depth >= max_depth) continue;
            if (v == sink) continue;

            GfaNodeSelector::NodeSimpleStats stats = selector_.build_stats_for_path(cur.path);

            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(2) << "expand " << get_node_name(v)
                               << "  depth=" << cur.depth
                               << "  score=" << cur.score
                               << "  path_len=" << cur.path.size() << "\n";
            }

            const auto &arcs = graph_.getArcsFromVertex(v);
            for (const GfaArc* a : arcs) {
                if (!a || a->get_del() || (skip_comp && a->get_comp())) continue;

                uint32_t w = a->get_target_vertex_id();

                if (w == v) continue;  // self-loop
                if (!in_topo_window(w)) continue;  // skip vertices outside the topo window defined by src and sink
                // If region_set is provided, restrict traversal to region_set.
                if (!region_set.empty() && w != sink && region_set.find(w) == region_set.end()) continue;
                if (std::find(cur.path.begin(), cur.path.end(), w) != cur.path.end()) continue;  // avoid cycle in prefix path

                uint32_t new_depth = cur.depth + 1;
                if (new_depth > max_depth) continue;

                // Calculate score
                GfaNodeSelector::ScoreContext sc;
                sc.region_nodes  = &region_set;
                sc.covered_nodes = &covered;
                int64_t step_score = selector_.compute_weight(v, w, stats, &sc);
                int64_t new_score  = cur.score + step_score;

                // Compare with existing best
                auto it_w = best_prefix.find(w);
                bool better = false;
                if (it_w == best_prefix.end()) {
                    better = true;
                } else {
                    const PrefixInfo &old = it_w->second;
                    if (w == sink) {
                        if (new_depth < old.depth) better = true;
                        else if (new_depth == old.depth && new_score > old.score) better = true;
                    } else {
                        if (new_score > old.score) better = true;
                        else if (new_score == old.score && new_depth < old.depth) better = true;
                    }
                }

                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3)
                                   << get_node_name(v) << " -> " << get_node_name(w)
                                   << "  step=" << step_score
                                   << "  total=" << new_score
                                   << "  depth=" << new_depth
                                   << (better ? "  update" : "  keep_old")
                                   << "\n";
                }

                if (!better) continue;

                PrefixInfo np;
                np.depth = new_depth;
                np.score = new_score;
                np.path  = cur.path;
                np.path.push_back(w);

                best_prefix[w] = std::move(np);
                covered.insert(w);

                if (w != sink) {
                    q.push(w);
                }
            }
        }

        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "best prefixes:\n";
            for (const auto &kv : best_prefix) {
                debug_stream() << log_indent(3)
                               << get_node_name(kv.first)
                               << "  path_len=" << kv.second.path.size()
                               << "  depth=" << kv.second.depth
                               << "  score=" << kv.second.score << "\n";
            }
        }

        if (best_prefix.size() <= 1) {
            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(2) << "only src reachable, no paths\n";
                debug_stream() << "=========================================\n\n";
            }
            return out;
        }

        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "reachable nodes with prefix: " << best_prefix.size() << "\n";
        }

        // Construct pivot candidates from best_prefix
        struct NodeStart {
            uint32_t v;
            uint32_t depth;
            int64_t  score;
            uint32_t len;
        };
        std::vector<NodeStart> starts;
        starts.reserve(best_prefix.size());

        // All nodes need to be covered
        std::unordered_set<uint32_t> to_cover;
        to_cover.reserve(best_prefix.size() * 2 + 8);

        for (const auto &kv : best_prefix) {
            uint32_t v = kv.first;
            const PrefixInfo &info = kv.second;

            if (!in_topo_window(v)) continue;  // skip vertices outside the topo window defined by src and sink

            to_cover.insert(v);

            if (v == src) continue;
            starts.push_back(NodeStart{
                v, 
                info.depth, 
                info.score, 
                graph_.getNodeLength(NodeHandle::get_segment_id(v))
            });
        }

        // Sort starts by depth , score , v
        std::sort(starts.begin(), starts.end(),
            [](const NodeStart &a, const NodeStart &b) {
                if (a.len != b.len) return a.len > b.len;
                if (a.depth != b.depth) return a.depth < b.depth;
                if (a.score != b.score) return a.score > b.score;
                return a.v < b.v;
            }
        );

        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "pivot candidates:\n";
            for (const auto &ns : starts) {
                debug_stream() << log_indent(3)
                               << get_node_name(ns.v)
                               << "  depth=" << ns.depth
                               << "  score=" << ns.score << "\n";
            }
        }

        // ================== Step 2: From pivots, do DFS to sink ================== //

        covered.clear();
        covered.reserve(region_set.size() + 8);

        // Deduplication
        std::unordered_set<uint64_t> seen_path;
        seen_path.reserve(64 + starts.size());

        // record which pivots have been tried
        std::unordered_set<uint32_t> tried_pivots;
        tried_pivots.reserve(starts.size() * 2 + 8);

        uint32_t stall_rounds = 0;

        uint64_t dfs_states = 0;

        if (DEBUG_ENABLED) {
            debug_stream() << "\n";
            debug_stream() << "  Pivot-to-sink DFS\n";
            debug_stream() << log_indent(2) << "stall_round_limit: " << stall_round_limit << "\n";
        }

        while (true) {
            if (DEBUG_ENABLED) {
                debug_stream() << "\n";
                debug_stream() << log_indent(2) << "round\n";
                debug_stream() << log_indent(3) << "covered      : " << covered.size() << " / " << to_cover.size() << "\n";
                debug_stream() << log_indent(3) << "path_count   : " << out.size() << "\n";
                debug_stream() << log_indent(3) << "stall_rounds : " << stall_rounds << "\n";
                debug_stream() << log_indent(3) << "tried_pivots : " << tried_pivots.size() << " / " << starts.size() << "\n";
            }

            if (max_paths > 0 && out.size() >= max_paths) break;
            if (stall_rounds >= stall_round_limit) break;

            uint32_t pivot = UINT32_MAX;

            // First priority: not tried && not covered
            for (const auto &ns : starts) {
                if (ns.depth == 0) continue;
                if (tried_pivots.find(ns.v) != tried_pivots.end()) continue;
                if (covered.find(ns.v) != covered.end()) continue;
                pivot = ns.v;
                break;
            }

            // Second priority: not tried (even if already covered, continue trying)
            if (pivot == UINT32_MAX) {
                for (const auto &ns : starts) {
                    if (ns.depth == 0) continue;
                    if (tried_pivots.find(ns.v) != tried_pivots.end()) continue;
                    pivot = ns.v;
                    break;
                }
            }

            if (pivot == UINT32_MAX) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "no more pivots to try\n";
                }
                break;
            }

            tried_pivots.insert(pivot);

            auto it_pfx = best_prefix.find(pivot);
            if (it_pfx == best_prefix.end()) {
                warning_stream() << "  ! " << get_node_name(pivot) << " has no prefix info, skip.\n";
                covered.insert(pivot);
                ++stall_rounds;
                continue;
            }

            const PrefixInfo &pfx = it_pfx->second;

            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(2) << "pivot: " << get_node_name(pivot)
                            << "  depth=" << pfx.depth
                            << "  score=" << pfx.score << "\n";
                log_path(pfx.path, 3, "prefix path:");
            }

            std::vector<uint32_t> path = pfx.path;
            std::unordered_set<uint32_t> on_path(path.begin(), path.end());

            struct Frame2 {
                uint32_t v;
                uint32_t depth_left;
                std::vector<uint32_t> cand_vs;
                size_t idx{0};
                std::string state = "normal";
            };

            auto build_frame2 = [&](uint32_t v, uint32_t depth_left) -> Frame2 {
                Frame2 fr;
                fr.v = v;
                fr.depth_left = depth_left;
                fr.idx = 0;

                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3)
                                << "frame at " << get_node_name(v)
                                << "  depth_left=" << depth_left << "\n";
                }

                const auto &arcs = graph_.getArcsFromVertex(v);

                if (arcs.empty()) {
                    fr.state = "end";
                    return fr;
                }

                fr.cand_vs.reserve(arcs.size());
                for (const GfaArc* a : arcs) {
                    if (!a || a->get_del() || (skip_comp && a->get_comp())) continue;

                    uint32_t w = a->get_target_vertex_id();
                    if (w == v) continue;
                    if (depth_left == 0) continue;
                    if (!in_topo_window(w)) continue;  // skip vertices outside the topo window defined by src and sink
                    if (!region_set.empty() && w != sink && region_set.find(w) == region_set.end()) continue;
                    if (on_path.find(w) != on_path.end()) continue;  // avoid cycle on current path

                    fr.cand_vs.push_back(w);
                }

                if (!fr.cand_vs.empty()) {
                    GfaNodeSelector::ScoreContext sc;
                    sc.region_nodes  = &region_set;
                    sc.covered_nodes = &covered;
                    selector_.sort_candidates(v, path, fr.cand_vs, &sc);
                }

                return fr;
            };

            std::vector<Frame2> stk;
            stk.reserve(max_depth + 4);

            uint32_t depth_left0 = (max_depth > pfx.depth) ? (max_depth - pfx.depth) : 0;
            stk.push_back(build_frame2(pivot, depth_left0));

            bool reached_sink = false;
            bool legal_path = true;

            while (!stk.empty()) {
                Frame2 &fr = stk.back();
                uint32_t v = fr.v;

                if (++dfs_states > DFS_guard) {
                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(3)
                                    << "stop: region too complex from "
                                    << get_node_name(src) << " to " << get_node_name(sink)
                                    << " at " << get_node_name(v) << "\n";
                    }
                    hit_limits = true;
                    return out;
                }

                if (v == sink) {
                    reached_sink = true;
                    if (DEBUG_ENABLED) {
                        log_path(path, 3, "reached sink, full path:");
                    }
                    break;
                }

                if ((fr.state == "end" || fr.depth_left == 0) && v != sink) {
                    legal_path = false;
                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(3)
                                    << "dead end at " << get_node_name(v)
                                    << " before reaching sink\n";
                    }
                    break;
                }

                if (fr.idx >= fr.cand_vs.size()) {
                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(4)
                                    << "backtrack from " << get_node_name(v) << "\n";
                    }
                    on_path.erase(v);
                    path.pop_back();
                    stk.pop_back();
                    continue;
                }

                uint32_t w = fr.cand_vs[fr.idx++];
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(4)
                                << "try " << get_node_name(v)
                                << " -> " << get_node_name(w) << "\n";
                }

                on_path.insert(w);
                path.push_back(w);

                uint32_t next_depth = (fr.depth_left > 0) ? (fr.depth_left - 1) : 0;
                stk.push_back(build_frame2(w, next_depth));
            }

            if (!legal_path || !reached_sink) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3)
                                << "path discarded, mark pivot/path nodes as covered\n";
                }
                covered.insert(pivot);
                for (uint32_t v : path) {
                    covered.insert(v);
                }
                ++stall_rounds;
                continue;
            }

            // Record path and update covered
            bool is_new_path = false;
            {
                uint64_t hv = hash_vec_(path);
                is_new_path = seen_path.insert(hv).second;
            }

            if (is_new_path) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3)
                                << "record unique path #" << (out.size() + 1)
                                << "  len=" << path.size() << "\n";
                }
                out.push_back(path);
                stall_rounds = 0;  // New paths found, reset to zero
            } else {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3)
                                << "duplicate path, skip record\n";
                }
                ++stall_rounds;
            }

            for (uint32_t v : path) {
                if (to_cover.find(v) != to_cover.end()) {
                    covered.insert(v);
                }
            }
        }

        if (DEBUG_ENABLED) {
            debug_stream() << "\n";
            debug_stream() << "  Enumeration finished\n";
            debug_stream() << log_indent(2) << "total paths: " << out.size() << "\n";
            debug_stream() << "=========================================\n";
        }

        return out;
    }

    /**
     * @brief Enumerate representative paths between two vertices using greedy DFS.
     *
     * Strategy:
     *   1. Start from @p src.
     *   2. At each step, enumerate all valid outgoing neighbors.
     *   3. Rank candidates using the following priority:
     *        a. Prefer vertices not yet present in @p covered.
     *        b. Prefer larger edge overlap (get_edge_ow()).
     *        c. Prefer longer node sequence length.
     *        d. Prefer smaller vertex ID for deterministic tie-breaking.
     *   4. Choose the best candidate and continue until:
     *        - @p sink is reached,
     *        - maximum depth is reached,
     *        - no valid successor exists,
     *        - DFS_guard is exceeded.
     *   5. Record the path if it is unique.
     *   6. Mark all visited vertices as covered so subsequent rounds prefer unexplored branches.
     *   7. Stop when:
     *        - max_paths paths have been found,
     *        - too many consecutive rounds produce no new path
     *          (stall_round_limit),
     *        - DFS_guard is exceeded.
     *
     * Special case:
     *   - If @p src == @p sink, the function returns a single path containing only @p src.
     *
     * @param src                Source oriented vertex ID.
     * @param sink               Sink oriented vertex ID.
     * @param region_set         Optional set restricting allowed internal vertices. If empty, all vertices are allowed.
     * @param max_depth          Maximum number of edges traversed in a single path.
     * @param max_paths          Maximum number of unique paths to return. If 0, no explicit path-count limit is applied.
     * @param skip_comp          If true, edges marked as composite (GfaArc::get_comp()) are ignored.
     * @param hit_limits         Output flag set to true if DFS_guard is exceeded.
     * @param DFS_guard          Maximum total number of DFS state expansions across all rounds.
     * @param stall_round_limit  Maximum number of consecutive rounds that fail to produce a new unique path.
     *
     * @return A vector of unique paths from @p src to @p sink. Each path is represented as a vector of oriented vertex IDs.
     */
    std::vector<std::vector<uint32_t>> enumerate_paths_greedy_DFS_(
        const uint32_t src,
        const uint32_t sink,
        const std::unordered_set<uint32_t>& region_set,
        const uint32_t max_depth,
        const uint32_t max_paths,
        const bool skip_comp,
        bool& hit_limits,
        const uint64_t DFS_guard,
        const uint32_t stall_round_limit = 2
    ) const
    {
        // ------------------------------------------------ Initialization ------------------------------------------------
        std::vector<std::vector<uint32_t>> out;
        hit_limits = false;

        // source and sink are the same.
        if (src == sink) {
            out.push_back({src});
            return out;
        }

        const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
        if (V == 0 || src >= V || sink >= V) return out;

        // ------------------------------------------------ Topological window ------------------------------------------------
        // Restrict the DFS walk to the topological interval between src and sink.
        uint32_t src_rank = graph_.vertex_topo_rank(Vertex(src));
        uint32_t sink_rank = graph_.vertex_topo_rank(Vertex(sink));

        const bool use_topo_filter =
            graph_.get_vertex_topological_index().size() >= V &&
            src_rank != UINT32_MAX &&
            sink_rank != UINT32_MAX;

        if (src_rank > sink_rank) std::swap(src_rank, sink_rank);

        auto in_topo_window = [&](uint32_t v) -> bool {
            if (!use_topo_filter) return true;

            const uint32_t r = graph_.vertex_topo_rank(Vertex(v));
            return r != UINT32_MAX && r >= src_rank && r <= sink_rank;
        };

        // ------------------------------------------------ Helper functions ------------------------------------------------
        auto valid_next = [&](uint32_t from, uint32_t w, const std::unordered_set<uint32_t>& on_path_seg ) -> bool {
            if (w == from) return false;
            if (!in_topo_window(w)) return false;

            // If region_set is provided, internal vertices must stay inside it.
            // The sink is always allowed.
            if (!region_set.empty() && w != sink && region_set.find(w) == region_set.end()) {
                return false;
            }

            // Segment-level cycle prevention.
            const uint32_t w_seg = NodeHandle::get_segment_id(w);
            if (on_path_seg.find(w_seg) != on_path_seg.end()) return false;

            return true;
        };

        // Encode a directed edge as a 64-bit key.
        auto edge_key = [](uint32_t a, uint32_t b) -> uint64_t {
            return (uint64_t(a) << 32) | uint64_t(b);
        };

        auto log_indent = [&](int level) -> std::string {
            return std::string(level * 2, ' ');
        };

        auto log_path = [&](const std::vector<uint32_t>& p, int level, const std::string& title) {
            if (!DEBUG_ENABLED) return;

            debug_stream() << log_indent(level) << title << "\n";
            for (uint32_t x : p) {
                debug_stream() << log_indent(level + 1) << get_node_name(x) << "\n";
            }
        };

        // ------------------------------------------------ Debug header ------------------------------------------------
        if (DEBUG_ENABLED) {
            debug_stream() << "=========================================\n";
            debug_stream() << "Enumerate greedy DFS paths\n";
            debug_stream() << log_indent(1) << "src       : " << get_node_name(src) << "\n";
            debug_stream() << log_indent(1) << "sink      : " << get_node_name(sink) << "\n";
            debug_stream() << log_indent(1) << "max_depth : " << max_depth << "\n";
            debug_stream() << log_indent(1) << "max_paths : " << max_paths << "\n";
            debug_stream() << log_indent(1) << "DFS_guard : " << DFS_guard << "\n";
            debug_stream() << log_indent(1) << "stall_lim : " << stall_round_limit << "\n";

            if (use_topo_filter) {
                debug_stream() << log_indent(1) << "topo_win  : " << src_rank << " .. " << sink_rank << "\n";
            }

            debug_stream() << log_indent(1) << "region_set size: " << region_set.size() << "\n";

            if (!region_set.empty()) {
                debug_stream() << log_indent(1) << "region_set:\n";
                for (uint32_t v : region_set) {
                    debug_stream() << log_indent(2) << get_node_name(v) << "\n";
                }
            }
        }

        // ------------------------------------------------ Global state across greedy rounds ------------------------------------------------
        // Used to sort the candidates nodes at each step, to prefer those not yet covered by previous paths.
        std::unordered_set<uint32_t> covered;
        covered.reserve(region_set.size() + 8);
        covered.insert(src);

        // Deduplicate complete paths.
        std::unordered_set<uint64_t> seen_path;
        seen_path.reserve(64);

        /**
        * @brief Count how many times each directed edge was chosen.
        * @date 2026-05-18
        * @version 0.1.3-r2
        * @note Used to prevent the greedy DFS from repeatedly choosing the same edge (used in sorting) and getting stuck in a local region.
        */
        std::unordered_map<uint64_t, uint32_t> edge_try_count;
        edge_try_count.reserve(256);

        uint64_t dfs_states = 0;
        uint32_t stall_rounds = 0;
        uint32_t round_id = 0;

        // ------------------------------------------------ Greedy path enumeration ------------------------------------------------
        // Each round performs one greedy walk from src to sink:
        //   1. Start from src.
        //   2. At each vertex, collect valid outgoing candidates.
        //   3. Sort candidates by priority.
        //   4. Choose the best candidate only.
        //   5. If sink is reached, record the path.
        while (true) {
            if (DEBUG_ENABLED) {
                debug_stream() << "\n";
                debug_stream() << log_indent(1) << "round " << round_id << "\n";
                debug_stream() << log_indent(2) << "path_count   : " << out.size() << "\n";
                debug_stream() << log_indent(2) << "covered      : " << covered.size() << "\n";
                debug_stream() << log_indent(2) << "stall_rounds : " << stall_rounds << "\n";
                debug_stream() << log_indent(2) << "dfs_states   : " << dfs_states << "\n";
            }

            // Stop when enough paths have been collected.
            if (max_paths > 0 && out.size() >= max_paths) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "stop: reach max_paths\n";
                }
                break;
            }

            // Stop after repeated rounds fail to produce new valid paths.
            if (stall_rounds >= stall_round_limit) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "stop: reach stall_round_limit\n";
                }
                break;
            }

            ++round_id;

            // Current greedy path.
            std::vector<uint32_t> path;
            path.reserve(max_depth + 1);
            path.push_back(src);

            // Segment-level nodes already used in the current path.
            std::unordered_set<uint32_t> on_path_seg;
            on_path_seg.reserve(max_depth + 1);
            on_path_seg.insert(NodeHandle::get_segment_id(src));

            uint32_t v = src;
            uint32_t depth = 0;
            bool reached_sink = false;

            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(2) << "start greedy walk from " << get_node_name(src) << "\n";
            }

            // ------------------------------------------------ One greedy walk ------------------------------------------------
            while (depth < max_depth) {
                // Global guard to avoid complex local graph structures.
                if (++dfs_states > DFS_guard) {
                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(2) << "stop: exceed DFS_guard = " << DFS_guard << " at " << get_node_name(v) << "\n";
                        log_path(path, 3, "current path:");
                        debug_stream() << "=========================================\n";
                    }

                    hit_limits = true;
                    return out;
                }

                if (v == sink) {
                    reached_sink = true;

                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(2) << "reached sink\n";
                    }

                    break;
                }

                // ------------------------------------------------ Candidate collection ------------------------------------------------
                std::vector<uint32_t> cand;

                const auto& arcs = graph_.getArcsFromVertex(v);
                cand.reserve(arcs.size());

                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "expand " << get_node_name(v) << "  depth=" << depth << "  out_arcs=" << arcs.size() << "\n";
                }

                for (const GfaArc* a : arcs) {
                    if (!a) {
                        if (DEBUG_ENABLED) {
                            debug_stream() << log_indent(3) << "skip null arc\n";
                        }
                        continue;
                    }

                    if (a->get_del()) {
                        if (DEBUG_ENABLED) {
                            debug_stream() << log_indent(3) << "skip deleted arc\n";
                        }
                        continue;
                    }

                    if (skip_comp && a->get_comp()) {
                        if (DEBUG_ENABLED) {
                            debug_stream() << log_indent(3) << "skip comp arc\n";
                        }
                        continue;
                    }

                    const uint32_t w = a->get_target_vertex_id();

                    if (!valid_next(v, w, on_path_seg)) {
                        if (DEBUG_ENABLED) {
                            debug_stream() << log_indent(3) << "skip " << get_node_name(v) << " -> " << get_node_name(w);

                            if (w == v) {
                                debug_stream() << "  reason=self-loop";
                            } else if (!in_topo_window(w)) {
                                debug_stream() << "  reason=outside-topo-window";
                            } else if (!region_set.empty() && w != sink && region_set.find(w) == region_set.end()) {
                                debug_stream() << "  reason=outside-region";
                            } else if (on_path_seg.find(NodeHandle::get_segment_id(w)) != on_path_seg.end()) {
                                debug_stream() << "  reason=cycle-on-path";
                            }

                            debug_stream() << "\n";
                        }

                        continue;
                    }

                    cand.push_back(w);

                    if (DEBUG_ENABLED) {
                        const uint64_t ek = edge_key(v, w);

                        const auto it_try = edge_try_count.find(ek);
                        const uint32_t tried = (it_try == edge_try_count.end()) ? 0u : it_try->second;

                        debug_stream() << log_indent(3) << "candidate " << get_node_name(w) << "  len=" << graph_.getNodeLength(NodeHandle::get_segment_id(w)) << "  used=" << (covered.find(w) != covered.end() ? "yes" : "no") << "  edge_try=" << tried << "\n";
                    }
                }

                if (cand.empty()) {
                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(2) << "dead end at " << get_node_name(v) << ", no valid candidate\n";
                    }

                    break;
                }

                // ------------------------------------------------ Candidate selection ------------------------------------------------
                // Candidate priority:
                //   1. Edges used fewer times in previous greedy rounds.
                //   2. Vertices not covered by previous paths.
                //   3. Larger overlap weight.
                //   4. Longer node length.
                //   5. Smaller vertex ID as deterministic tie-breaker.
                uint32_t best = UINT32_MAX;

                if (cand.size() == 1) {
                    best = cand.front();

                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(2) << "choose only candidate: " << get_node_name(best) << "\n";
                    }
                } else {
                    std::sort(cand.begin(), cand.end(),
                        [&](uint32_t a, uint32_t b) {
                            const uint64_t eka = edge_key(v, a);
                            const uint64_t ekb = edge_key(v, b);

                            const auto ita = edge_try_count.find(eka);
                            const auto itb = edge_try_count.find(ekb);

                            const uint32_t ta = (ita == edge_try_count.end()) ? 0u : ita->second;
                            const uint32_t tb = (itb == edge_try_count.end()) ? 0u : itb->second;

                            if (ta != tb) return ta < tb;

                            const bool a_unused = covered.find(a) == covered.end();
                            const bool b_unused = covered.find(b) == covered.end();

                            if (a_unused != b_unused) return a_unused > b_unused;

                            const uint32_t oa = graph_.get_edge_ow(v, a);
                            const uint32_t ob = graph_.get_edge_ow(v, b);

                            if (oa != ob) return oa > ob;

                            const uint32_t la = graph_.getNodeLength(NodeHandle::get_segment_id(a));
                            const uint32_t lb = graph_.getNodeLength(NodeHandle::get_segment_id(b));

                            if (la != lb) return la > lb;

                            return a < b;
                        }
                    );

                    best = cand.front();

                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(2) << "sorted candidates:\n";

                        for (uint32_t x : cand) {
                            const uint64_t ek = edge_key(v, x);

                            const auto it_try = edge_try_count.find(ek);
                            const uint32_t tried = (it_try == edge_try_count.end()) ? 0u : it_try->second;

                            debug_stream() << log_indent(3) << get_node_name(x) << "  len=" << graph_.getNodeLength(NodeHandle::get_segment_id(x)) << "  used=" << (covered.find(x) != covered.end() ? "yes" : "no") << "  edge_try=" << tried << "\n";
                        }

                        debug_stream() << log_indent(2) << "choose best candidate: " << get_node_name(best) << "\n";
                    }
                }

                if (best == UINT32_MAX) {
                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(2) << "no best candidate, break\n";
                    }

                    break;
                }

                // Mark the selected directed edge as tried.
                ++edge_try_count[edge_key(v, best)];

                // Extend the current path.
                path.push_back(best);
                on_path_seg.insert(NodeHandle::get_segment_id(best));

                v = best;
                ++depth;
            }

            if (v == sink) reached_sink = true;

            // ------------------------------------------------ Handle failed path ------------------------------------------------
            // If this greedy walk did not reach sink, the partial path is discarded, but its vertices are still marked as covered.
            if (!reached_sink) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "path discarded, did not reach sink\n";
                    log_path(path, 3, "discarded path:");
                    debug_stream() << log_indent(2) << "mark path nodes as covered\n";
                }

                for (uint32_t x : path) {
                    covered.insert(x);
                }

                ++stall_rounds;

                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "stall_rounds -> " << stall_rounds << "\n";
                }

                continue;
            }

            // ------------------------------------------------ Deduplicate complete path ------------------------------------------------
            // Only exactly identical vertex sequences are considered duplicates.
            bool is_new_path = false;
            {
                const uint64_t hv = hash_vec_(path);
                is_new_path = seen_path.insert(hv).second;
            }

            for (uint32_t x : path) {
                covered.insert(x);
            }

            if (is_new_path) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "record unique path #" << (out.size() + 1) << "  len=" << path.size() << "\n";
                    log_path(path, 3, "full path:");
                }

                out.push_back(path);
                stall_rounds = 0;
            } else {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(2) << "duplicate path, skip record\n";
                    log_path(path, 3, "duplicate path:");
                }

                ++stall_rounds;
            }
        }

        // ------------------------------------------------ Summary ------------------------------------------------
        if (DEBUG_ENABLED) {
            debug_stream() << "\n";
            debug_stream() << "  Greedy enumeration finished\n";
            debug_stream() << log_indent(2) << "total paths : " << out.size() << "\n";
            debug_stream() << log_indent(2) << "covered     : " << covered.size() << "\n";
            debug_stream() << log_indent(2) << "dfs_states  : " << dfs_states << "\n";
            debug_stream() << log_indent(2) << "rounds      : " << round_id << "\n";
            debug_stream() << "=========================================\n";
        }

        return out;
    }


private:
    const GfaGraph& graph_;
    GfaNodeSelector selector_;
    const uint32_t  V_;  // 2 * #segments
};