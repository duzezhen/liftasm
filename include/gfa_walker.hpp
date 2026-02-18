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
        bool skip_comp, bool& hit_limits, const uint64_t DFS_guard
    ) const {
        return enumerate_paths_DFS_(src, sink, region_set, max_depth, max_paths, skip_comp, hit_limits, DFS_guard);
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

        const std::string oriented = graph.get_oriented_sequence(v);
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

        const std::string oriented = graph.get_oriented_sequence(v);
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

    // // Enumerate all paths from src to sink using DFS, with optional depth and path limits.
    // std::vector<std::vector<uint32_t>> enumerate_paths_DFS_(
    //     const uint32_t src, const uint32_t sink, 
    //     const uint32_t max_depth, const uint32_t max_paths,
    //     const bool skip_comp = false
    // ) const {
    //     std::vector<std::vector<uint32_t>> out;
        
    //     // Check
    //     const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    //     if (V == 0 || src >= V || sink >= V) return out;

    //     // used to detect cycles on the current path
    //     std::unordered_set<uint32_t> on_path;
    //     on_path.reserve(max_depth + 1);

    //     // used to avoid duplicate paths (by hash)
    //     std::unordered_set<uint64_t> seen_path; seen_path.reserve(32);

    //     // DFS stack
    //     std::vector<uint32_t> path; path.reserve(max_depth + 1);
    //     std::vector<DFSFrame> stk;  stk.reserve(max_depth + 1);

    //     on_path.insert(src);
    //     path.push_back(src);
    //     stk.push_back({src, 0, static_cast<uint64_t>(max_depth)});

    //     // Limit the number of DFS states to prevent excessive CPU time usage
    //     static constexpr std::uint64_t MAX_DFS_STATES = 500'000;
    //     std::uint64_t dfs_states = 0;

    //     while (!stk.empty() && out.size() < max_paths) {
    //         if (++dfs_states > MAX_DFS_STATES) {
    //             out.clear();
    //             return out;
    //         }

    //         DFSFrame& fr = stk.back();
    //         uint32_t  v  = fr.v;
    //         if (v == sink) {
    //             uint64_t hv = hash_vec_(path);  // Deduplicate by hash
    //             if (seen_path.insert(hv).second) {
    //                 out.push_back(path);
    //                 if (out.size() >= max_paths) break;
    //             }
    //             on_path.erase(v);
    //             path.pop_back();
    //             stk.pop_back();
    //             continue;
    //         }

    //         const auto& arcs = graph_.getArcsFromVertex(v);
    //         bool advanced = false;
    //         while (fr.nextIdx < arcs.size()) {
    //             const GfaArc* a = arcs[fr.nextIdx++];
    //             if (!a || a->get_del() || (skip_comp && a->get_comp())) continue;
    //             uint32_t w = a->get_target_vertex_id();
    //             if (w == v) continue;
    //             if (on_path.find(w) != on_path.end()) continue;
    //             if (fr.depthLeft == 0) continue;
    //             on_path.insert(w);
    //             path.push_back(w);
    //             stk.push_back({w, 0, static_cast<uint64_t>(fr.depthLeft - 1)});
    //             advanced = true;
    //             break;
    //         }
    //         if (!advanced) {
    //             on_path.erase(v);
    //             path.pop_back();
    //             stk.pop_back();
    //         }
    //     }

    //     return out;
    // }

    const std::string& get_node_name(uint32_t v) const {
        uint32_t sid = NodeHandle::get_segment_id(v);
        const GfaNode* nd = graph_.getNode(sid);
        static const std::string empty = "";
        return nd ? nd->name : empty;
    }

    std::vector<std::vector<uint32_t>> enumerate_paths_DFS_(
        const uint32_t src,
        const uint32_t sink,
        const std::unordered_set<uint32_t>& region_set,
        const uint32_t max_depth,
        const uint32_t max_paths,
        const bool skip_comp, 
        bool& hit_limits, 
        const uint64_t DFS_guard
    ) const
    {
        std::vector<std::vector<uint32_t>> out;

        std::unordered_set<uint32_t> covered;
        covered.reserve(region_set.size() + 8);

        const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
        if (V == 0 || src >= V || sink >= V) return out;

        debug_stream() << "=========================================\n";
        debug_stream() << "[ENUM-DFS] src=" << get_node_name(src)
                    << " sink=" << get_node_name(sink)
                    << " max_depth=" << max_depth
                    << " max_paths=" << max_paths
                    << "\n";

        // print region_set
        debug_stream() << "[ENUM-DFS] region_set nodes:\n";
        for (uint32_t v : region_set) {
            debug_stream() << "  " << get_node_name(v) << "\n";
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
            best_prefix[src] = info;
        }

        std::queue<uint32_t> q;
        q.push(src);

        uint64_t prefix_states = 0;

        debug_stream() << "[PREFIX] start forward exploration from src\n";

        while (!q.empty()) {
            uint32_t v = q.front();
            q.pop();

            auto it_v = best_prefix.find(v);
            if (it_v == best_prefix.end()) continue;
            const PrefixInfo &cur = it_v->second;

            if (++prefix_states > DFS_guard) {
                debug_stream() << "[PREFIX] STOP: exceed DFS_guard=" << DFS_guard << "\n";
                break;
            }

            if (cur.depth >= max_depth) continue;
            if (v == sink) continue;

            GfaNodeSelector::NodeSimpleStats stats = selector_.build_stats_for_path(cur.path);

            debug_stream() << "[PREFIX] expand from " << get_node_name(v)
                << " depth=" << cur.depth
                << " score=" << cur.score
                << " path_len=" << cur.path.size()
                << "\n";

            const auto &arcs = graph_.getArcsFromVertex(v);
            for (const GfaArc* a : arcs) {
                if (!a || a->get_del() || (skip_comp && a->get_comp())) continue;

                uint32_t w = a->get_target_vertex_id();

                if (w == v) continue;  // Self-loop
                if (w != sink && region_set.find(w) == region_set.end()) continue;  // Not in region
                if (std::find(cur.path.begin(), cur.path.end(), w) != cur.path.end()) continue;  // Avoid cycles in the prefix path

                uint32_t new_depth = cur.depth + 1;
                if (new_depth > max_depth) continue;

                // Canculate score
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

                debug_stream() << "  [PREFIX-CAND] from " << get_node_name(v)
                            << " to "   << get_node_name(w)
                            << " step_score=" << step_score
                            << " new_score="  << new_score
                            << " new_depth="  << new_depth
                            << (better ? " [UPDATE]\n" : " [KEEP-OLD]\n");

                if (!better) continue;

                PrefixInfo np;
                np.depth = new_depth;
                np.score = new_score;
                np.path  = cur.path;
                np.path.push_back(w);

                best_prefix[w] = std::move(np);

                covered.insert(w);

                // Push node into queue
                if (w != sink) {
                    q.push(w);
                }
            }
        }

        // print best_prefix
        for (const auto &kv : best_prefix) {
            debug_stream() << "[BEST_PREFIX]  node=" << get_node_name(kv.first) << " path_len=" << kv.second.path.size() << " depth=" << kv.second.depth << " score=" << kv.second.score << "\n";
        }

        if (best_prefix.size() <= 1) {
            debug_stream() << "[PREFIX] only src reachable, no paths.\n";
            debug_stream() << "=========================================\n\n";
            return out;
        }

        debug_stream() << "[PREFIX] total nodes with prefix=" << best_prefix.size() << "\n";

        // Construct pivot candidates from best_prefix
        struct NodeStart {
            uint32_t v;
            uint32_t depth;
            int64_t  score;
        };
        std::vector<NodeStart> starts;
        starts.reserve(best_prefix.size());

        // All nodes need to be covered
        std::unordered_set<uint32_t> to_cover;
        to_cover.reserve(best_prefix.size() * 2 + 8);

        for (const auto &kv : best_prefix) {
            uint32_t v = kv.first;
            const PrefixInfo &info = kv.second;
            to_cover.insert(v);

            if (v == src) continue;
            starts.push_back(NodeStart{v, info.depth, info.score});
        }

        // Sort starts by depth , score , v
        std::sort(starts.begin(), starts.end(),
            [](const NodeStart &a, const NodeStart &b) {
                if (a.depth != b.depth) return a.depth < b.depth;
                if (a.score != b.score) return a.score < b.score;
                return a.v < b.v;
            }
        );

        debug_stream() << "[PREFIX] pivot candidates (by depth desc):\n";
        for (const auto &ns : starts) {
            debug_stream() << "  pivot=" << get_node_name(ns.v) << " depth=" << ns.depth << " score=" << ns.score << "\n";
        }

        // ================== Step 2: From pivots, do DFS to sink, until all nodes covered ================== //

        covered.clear();
        covered.reserve(region_set.size() + 8);

        // Deduplication
        std::unordered_set<uint64_t> seen_path;
        seen_path.reserve(64);

        // Defensive upper bound
        uint64_t dfs_states = 0;

        auto all_covered = [&]() -> bool {
            if (to_cover.empty()) return true;
            for (uint32_t v : to_cover) {
                if (covered.find(v) == covered.end()) return false;
            }
            return true;
        };

        debug_stream() << "[DFS] start second-phase DFS from pivots.\n";

        while (true) {
            debug_stream() << "-------------------------------------------\n";
            debug_stream() << "[DFS-ROUND] covered=" << covered.size() << " / need_cover=" << to_cover.size() << "  paths=" << out.size() << "\n";

            if (all_covered()) break;
            if (max_paths > 0 && out.size() >= max_paths) break;

            // Select next pivot
            uint32_t pivot = UINT32_MAX;
            for (const auto &ns : starts) {
                if (to_cover.find(ns.v) == to_cover.end()) continue;
                if (covered.find(ns.v) != covered.end()) continue;
                if (ns.depth == 0) continue;
                pivot = ns.v;
                break;
            }

            if (pivot == UINT32_MAX) break;

            auto it_pfx = best_prefix.find(pivot);
            if (it_pfx == best_prefix.end()) {
                warning_stream() << get_node_name(pivot) << " has no prefix info, skip.\n";
                covered.insert(pivot);
                continue;
            }

            const PrefixInfo &pfx = it_pfx->second;

            debug_stream() << "[PIVOT] node=" << get_node_name(pivot) << " depth=" << pfx.depth << " score=" << pfx.score << "\n";

            std::vector<uint32_t> path = pfx.path;
            std::unordered_set<uint32_t> on_path(path.begin(), path.end());

            struct Frame2 {
                uint32_t v;
                uint32_t depth_left;
                std::vector<uint32_t> cand_vs;
                size_t idx;
                std::string state = "normal";
            };

            auto build_frame2 = [&](uint32_t v, uint32_t depth_left) -> Frame2 {
                Frame2 fr;
                fr.v = v;
                fr.depth_left = depth_left;
                fr.idx = 0;

                debug_stream() << "[FRAME] at " << get_node_name(v) << " depth_left=" << depth_left << "\n";

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
                    if (on_path.find(w) != on_path.end()) continue;  // Acoid cycles on current path

                    fr.cand_vs.push_back(w);
                }

                if (!fr.cand_vs.empty()) {
                    GfaNodeSelector::ScoreContext sc;
                    sc.region_nodes  = &region_set;
                    sc.covered_nodes = &covered;
                    selector_.sort_candidates(v, path, fr.cand_vs, &sc);
                }

                // Print log
                if (!fr.cand_vs.empty()) {
                    GfaNodeSelector::NodeSimpleStats stats = selector_.build_stats_for_path(path);
                    for (auto w : fr.cand_vs) {
                        GfaNodeSelector::ScoreContext sc;
                        sc.region_nodes  = &region_set;
                        sc.covered_nodes = &covered;
                        int64_t wt = selector_.compute_weight(v, w, stats, &sc);
                    }
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
                    debug_stream() << "   - Region too complex from " << get_node_name(src) << " to " << get_node_name(sink) << " at " << get_node_name(v) << ", skip enumeration.\n";
                    out.clear();
                    hit_limits = true;
                    return out;
                }

                if (v == sink) {
                    reached_sink = true;
                    debug_stream() << "[HIT-SINK] path:\n";
                    for (auto x : path) debug_stream() << "  +node " << get_node_name(x) << "\n";
                    break;
                }

                if ((fr.state == "end" || fr.depth_left == 0) && v != sink) {
                    legal_path = false;
                    break;
                }

                if (fr.idx >= fr.cand_vs.size()) {
                    debug_stream() << "[BACKTRACK] at " << get_node_name(v) << "\n";
                    on_path.erase(v);
                    path.pop_back();
                    stk.pop_back();
                    continue;
                }

                uint32_t w = fr.cand_vs[fr.idx++];
                debug_stream() << "[TRY] " << get_node_name(v) << " -> " << get_node_name(w) << "\n";

                on_path.insert(w);
                path.push_back(w);

                uint32_t next_depth = (fr.depth_left > 0) ? (fr.depth_left - 1) : 0;
                stk.push_back(build_frame2(w, next_depth));
            }

            if (!legal_path || !reached_sink) {
                covered.insert(pivot);
                for (uint32_t v : path) {
                    covered.insert(v);
                }
                continue;
            }

            // Record path
            uint64_t hv = hash_vec_(path);
            if (!seen_path.insert(hv).second) {
                debug_stream() << "[DUP-PATH] already seen, but still update coverage.\n";
            }

            debug_stream() << "[RECORD-PATH] idx=" << (out.size() + 1) << " len=" << path.size() << "\n";
            out.push_back(path);

            bool new_cover = false;
            for (uint32_t v : path) {
                if (to_cover.find(v) != to_cover.end()) {
                    if (covered.insert(v).second) {
                        new_cover = true;
                        debug_stream() << "  +cover " << get_node_name(v) << "\n";
                    }
                }
            }

            if (!new_cover) {
                debug_stream() << "[NO-NEW-COVER] this path does not cover new nodes.\n";
            }
        }

        debug_stream() << "[ENUM-DFS-END] total paths=" << out.size() << "\n";
        debug_stream() << "=========================================\n";

        return out;
    }


private:
    const GfaGraph& graph_;
    GfaNodeSelector selector_;
    const uint32_t  V_;  // 2 * #segments
};