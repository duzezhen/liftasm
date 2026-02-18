#include "../include/gfa_bubble.hpp"
#include "../include/progress_tracker.hpp"
#include "../include/ThreadPool.hpp"

#include <numeric>

namespace GfaBubble {

GfaBubbleFinder::GfaBubbleFinder(
    const GfaGraph& g, std::size_t md, std::size_t mp, uint64_t dfs_guard, 
        uint16_t cx_depth, uint32_t cx_nodes, uint32_t cx_branches, int cx_deg_branch, int cx_deg_hub, int cx_deg_cap, 
        double path_diff, bool sc, uint32_t thread)
    : graph_(g), max_depth_(md), max_paths_(mp), dfs_guard_(dfs_guard), 
    cx_depth_(cx_depth), cx_nodes_(cx_nodes), cx_branches_(cx_branches), cx_deg_branch_(cx_deg_branch), cx_deg_hub_(cx_deg_hub), cx_deg_cap_(cx_deg_cap), 
    path_diff_(path_diff), skip_comp_(sc), thread_(thread)
{
    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    reach_buf_.assign(V, UNVIS_FLAG_);
    complex_flag_.assign(V, UNKNOWN_FLAG_);
}

bool GfaBubbleFinder::is_strict_source_(uint32_t v, std::vector<uint32_t>& tgt) const
{
    tgt.clear();
    tgt.reserve(8);

    std::unordered_set<uint32_t> seen_seg; seen_seg.reserve(8);

    for (const GfaArc* a : graph_.getArcsFromVertex(v)) {
        if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
        uint32_t w = a->get_target_vertex_id();
        uint32_t seg = NodeHandle::get_segment_id(w);
        if (!seen_seg.insert(seg).second) continue;
        if (w == v) continue;
        tgt.push_back(w);
    }

    return tgt.size() >= 2;
}

bool GfaBubbleFinder::is_strict_sink_(uint32_t v, uint32_t w) const {
    auto leads_uniquely_to = [&](uint32_t u, uint32_t target, int MAX) -> bool {
        if (u == target) return true;
        for (int step = 0; step < MAX; ++step) {
            uint32_t next = UINT32_MAX;
            int fanout = 0;
            for (const GfaArc* a : graph_.getArcsFromVertex(u)) {
                if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
                uint32_t t = a->get_target_vertex_id();
                if (t == u) continue;
                if (++fanout > 1) return false;
                next = t;
            }
            if (fanout == 0) return false;
            u = next;
            if (u == target) return true;
        }
        return (u == target);
    };

    constexpr int MAX_TUNNEL_STEPS = 4;

    for (const GfaArc* in : graph_.getArcsToVertex(w)) {
        if (!in || in->get_del() || (skip_comp_ && in->get_comp())) continue;
        const uint32_t pv = in->get_source_vertex_id();
        if (pv == w || pv == v) continue;

        int outdeg = 0;
        bool all_funnel_to_w = true;

        for (const GfaArc* a : graph_.getArcsFromVertex(pv)) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
            uint32_t x = a->get_target_vertex_id();
            if (x == pv) continue;
            ++outdeg;

            if (!leads_uniquely_to(x, w, MAX_TUNNEL_STEPS)) {
                all_funnel_to_w = false;
                if (outdeg > 1) return false;
            }
        }

        if (outdeg > 1 && !all_funnel_to_w) return false;
    }

    return true;
}

void GfaBubbleFinder::mark_nodes_complex_(const std::vector<uint32_t>& nodes)
{
    const uint32_t V = static_cast<uint32_t>(complex_flag_.size() * 2);
    if (V == 0) return;

    std::lock_guard<std::mutex> lk(complex_mtx_);
    for (uint32_t x : nodes) {
        if (x >= V) continue;
        complex_flag_[x] = COMPLEX_FLAG_;
    }
}
void GfaBubbleFinder::mark_nodes_simple_(const std::vector<uint32_t>& nodes)
{
    const uint32_t V = static_cast<uint32_t>(complex_flag_.size() * 2);
    if (V == 0) return;

    std::lock_guard<std::mutex> lk(complex_mtx_);
    for (uint32_t x : nodes) {
        if (x >= V) continue;
        if (complex_flag_[x] == UNKNOWN_FLAG_) {
            complex_flag_[x] = SIMPLE_FLAG_;
        }
    }
}

bool GfaBubbleFinder::is_complex_source_(
    uint32_t v,
    const std::vector<uint32_t>& seed_buf
) {
    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (complex_flag_.size() != V) return false;

    auto get_flag = [&](uint32_t x) -> uint8_t {
        if (x >= V) return UNKNOWN_FLAG_;
        std::lock_guard<std::mutex> lk(complex_mtx_);
        return complex_flag_[x];
    };
    auto is_marked_complex = [&](uint32_t x) -> bool {
        return get_flag(x) == COMPLEX_FLAG_;
    };
    auto is_marked_simple = [&](uint32_t x) -> bool {
        return get_flag(x) == SIMPLE_FLAG_;
    };

    // Count the total degreeof u
    auto total_deg = [&](uint32_t u) -> int {
        int deg = 0;

        const auto& out_arcs = graph_.getArcsFromVertex(u);
        for (const GfaArc* a : out_arcs) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
            uint32_t w = a->get_target_vertex_id();
            if (w == u) continue;
            ++deg;
            if (deg >= cx_deg_cap_) return deg;
        }

        return deg;
    };

    // 0. Use cached flags
    if (is_marked_complex(v)) return true;

    bool any_seed_complex = false;
    bool all_seed_simple  = !seed_buf.empty();
    for (uint32_t s : seed_buf) {
        if (is_marked_complex(s)) {
            any_seed_complex = true;
            break;
        }
        if (!is_marked_simple(s)) {
            all_seed_simple = false;
        }
    }
    if (any_seed_complex) return true;
    if (all_seed_simple)  return false;

    // 1. Quick hub screen: src or any seed with total_deg>=DEG_HUB is complex
    {
        if (total_deg(v) >= cx_deg_hub_) {
            mark_nodes_complex_({v});
            return true;
        }
        for (uint32_t s : seed_buf) {
            if (total_deg(s) >= cx_deg_hub_) {
                mark_nodes_complex_({v, s});
                return true;
            }
        }
    }

    // 2. Local BFS: count nodes / branching nodes
    std::queue<std::pair<uint32_t,uint16_t>> q; // (vertex, depth)
    std::unordered_set<uint32_t> visited;
    visited.reserve(256);

    std::vector<uint32_t> local_nodes;
    local_nodes.reserve(256);

    auto push_node = [&](uint32_t x, uint16_t d) {
        if (!visited.insert(x).second) return;
        q.emplace(x, d);
        local_nodes.push_back(x);
    };

    push_node(v, 0);
    for (uint32_t s : seed_buf) {
        if (s == v) continue;
        push_node(s, 0);
    }

    std::size_t branching_nodes = 0;

    while (!q.empty()) {
        auto [u, d] = q.front();
        q.pop();

        // Too many nodes; local region is large, mark complex
        if (visited.size() > cx_nodes_) {
            mark_nodes_complex_(local_nodes);
            return true;
        }

        int deg = total_deg(u);
        if (deg >= cx_deg_branch_) {
            ++branching_nodes;
            if (branching_nodes > cx_branches_) {
                mark_nodes_complex_(local_nodes);
                return true;
            }
        }

        // Depth limit
        if (d >= cx_depth_) continue;

        // Expand one layer (outgoing only to avoid double-counting in)
        const auto& arcs = graph_.getArcsFromVertex(u);
        for (const GfaArc* a : arcs) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
            uint32_t w = a->get_target_vertex_id();
            if (w == u) continue;

            if (is_marked_complex(w)) {
                local_nodes.push_back(w);
                mark_nodes_complex_(local_nodes);
                return true;
            }

            push_node(w, static_cast<uint16_t>(d + 1));
        }
    }

    mark_nodes_simple_(local_nodes);
    return false;
}

void GfaBubbleFinder::find_bubbles(const bool filter_nonlocal)
{
    log_stream() << "Detecting bubbles ...\n";
    bubbles_.clear();

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (V == 0) return;

    // 1. Collect strict sources and tips
    std::vector<uint32_t> sources;
    sources.reserve(V / 16 + 8);

    std::vector<uint32_t> tips;
    tips.reserve(V / 32 + 4);

    {
        std::vector<uint32_t> tmp;
        tmp.reserve(16);

        for (uint32_t v = 0; v < V; ++v) {
            if (graph_.getNodeDeleted(NodeHandle::get_segment_id(v))) continue;
            // 1.1 TIP
            const auto& outs = graph_.getArcsFromVertex(v);
            bool has_out = false;
            for (const GfaArc* a : outs) {
                if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
                uint32_t w = a->get_target_vertex_id();
                if (w == v) continue;
                has_out = true;
                break;
            }
            if (!has_out) {
                tips.push_back(v);
            }

            // 1.2 strict source
            tmp.clear();
            if (!is_strict_source_(v, tmp)) continue;
            if (tmp.size() > 64) continue;
            sources.push_back(v);
        }
    }

    log_stream() << "   - strict sources: " << sources.size() << ", tips: " << tips.size() << "\n";

    const uint64_t bfs_limit = static_cast<uint64_t>(std::min<std::size_t>(max_depth_ ? max_depth_ : 64, UINT64_MAX));

    // 2. Closed bubbles from strict sources
    auto pair_key = [](uint32_t s0, uint32_t s1) -> uint64_t {
        if (s0 > s1) std::swap(s0, s1);
        return (static_cast<uint64_t>(s0) << 32) | static_cast<uint64_t>(s1);
    };

    const unsigned threads = std::max(1u, thread_);
    ThreadPool pool(threads);

    std::vector<std::future<Bubble>> futs;
    futs.reserve(sources.size());

    for (uint32_t src : sources) {
        futs.emplace_back(
            pool.submit([&, src]() -> Bubble {
                return detect_closed_bubble_from_source_(src, bfs_limit);
            })
        );
    }

    std::unordered_set<uint64_t> seen_pair;
    seen_pair.reserve(2048);

    ProgressTracker prog(futs.size());
    for (auto &f : futs) {
        prog.hit();
        Bubble bb = f.get();
        const auto &paths = bb.get_paths();
        if (paths.empty()) continue;

        uint32_t src  = bb.get_source();
        uint32_t sink = bb.get_sink();

        uint32_t s0 = NodeHandle::get_segment_id(src);
        uint32_t s1 = NodeHandle::get_segment_id(sink);
        uint64_t key = pair_key(s0, s1);

        if (!seen_pair.insert(key).second) continue;

        bubbles_.emplace_back(std::move(bb));
    }

    pool.stop();

    // 3. Open bubbles from tips
    if (tips.size() >= 2) {
        find_open_bubbles_from_tips_(tips);
    }

    // 4. Post-processing
    if (filter_nonlocal) filter_nonlocal_bubbles_();
    log_stream() << "   - Total bubbles detected: " << bubbles_.size() << "\n";
}

Bubble GfaBubbleFinder::detect_closed_bubble_from_source_(uint32_t v, uint64_t bfs_limit)
{
    Bubble bb;

    auto vtx_name = [&](uint32_t x) -> std::string {
        uint32_t sid = NodeHandle::get_segment_id(x);
        bool     rv  = NodeHandle::get_is_reverse(x);
        return graph_.getNodeName(sid) + (rv ? "-" : "+");
    };

    debug_stream() << "[bubble-src] start from " << vtx_name(v) << "\n";

    // 1) Collect seeds
    std::vector<uint32_t> seed_buf;
    seed_buf.reserve(16);
    if (!is_strict_source_(v, seed_buf)) return bb;
    if (seed_buf.size() > 64) {
        debug_stream() << "  - too many seeds (" << seed_buf.size() << ") at " << vtx_name(v) << ", skip.\n";
        return bb;
    }

    const int16_t seeds_cnt = static_cast<int16_t>(seed_buf.size());
    if (seeds_cnt < 2) {
        debug_stream() << "  - seeds_cnt<2 at " << vtx_name(v) << ", skip.\n";
        return bb;
    }

    // 2) Complexity filter
    if (is_complex_source_(v, seed_buf)) {
        debug_stream() << "  - local region too complex, skip source " << vtx_name(v) << "\n";
        return bb;
    }

    // 3) meet_mask: vertex -> bitmask of branches
    const bool use_mask = (seeds_cnt > 1 && seeds_cnt <= 64);
    const uint64_t FULL = use_mask
        ? (seeds_cnt == 64 ? ~0ULL : ((1ULL << seeds_cnt) - 1))
        : 0ULL;

    std::unordered_map<uint32_t, uint64_t> meet_mask;
    meet_mask.reserve(128);

    // 4) Local nodes for this source (limit DFS enumeration)
    std::unordered_set<uint32_t> local_nodes;
    local_nodes.reserve(256);

    // 5) Multi-source BFS queue: (vertex, branch_id, depth)
    std::queue<std::tuple<uint32_t,int16_t,uint16_t>> q;

    // 6) visited on (vertex, branch)
    std::unordered_set<uint64_t> visited;
    visited.reserve(seeds_cnt * 32);

    auto encode_state = [](uint32_t node, int16_t br) -> uint64_t {
        return (static_cast<uint64_t>(node) << 8) | static_cast<uint64_t>(static_cast<uint8_t>(br));
    };

    // 7) Initialize seeds
    for (int16_t b = 0; b < seeds_cnt; ++b) {
        uint32_t w = seed_buf[static_cast<std::size_t>(b)];

        if (use_mask && b < 64) {
            uint64_t &m = meet_mask[w];
            m |= (1ULL << b);
        }

        uint64_t st = encode_state(w, b);
        visited.insert(st);
        q.emplace(w, b, 0u);

        local_nodes.insert(w);

        debug_stream() << "  [seed] src=" << vtx_name(v) << " seed_v=" << vtx_name(w) << " branch=" << b << "\n";
    }

    bool     early_all_met      = false;
    uint32_t global_sink        = UINT32_MAX;
    uint32_t sink_depth_from_src = 0;

    // 8) BFS: only find sink + local_nodes
    while (!q.empty()) {
        auto [cur, b, d] = q.front();
        q.pop();

        debug_stream() << "    [BFS] src=" << vtx_name(v)
                       << " cur=" << vtx_name(cur)
                       << " depth=" << d
                       << " branch=" << b
                       << " qsize=" << q.size() << "\n";

        if (cur == v) continue;
        if (d >= bfs_limit) {
            debug_stream() << "      - reach bfs_limit at " << vtx_name(cur) << "\n";
            return bb;
        }

        const auto& arcs = graph_.getArcsFromVertex(cur);

        for (const GfaArc* a : arcs) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
            uint32_t nxt = a->get_target_vertex_id();
            if (nxt == cur) continue;

            if (use_mask && b >= 0 && b < seeds_cnt) {
                uint64_t &m = meet_mask[nxt];
                m |= (1ULL << b);

                if (m == FULL && nxt != v) {
                    debug_stream() << "      - src=" << vtx_name(v)
                                   << "  cand_sink=" << vtx_name(nxt)
                                   << "  depth=" << d
                                   << "  seeds=" << seeds_cnt
                                   << "  qsize=" << q.size() << "\n";

                    if (!is_strict_sink_(v, nxt)) {
                        debug_stream() << "        - not a strict sink\n";
                    } else {
                        debug_stream() << "    + accept as bubble: src="
                                       << vtx_name(v) << " sink="
                                       << vtx_name(nxt) << "\n";
                        global_sink        = nxt;
                        early_all_met      = true;
                        sink_depth_from_src = static_cast<uint32_t>(d + 1);
                        break;
                    }
                }
            }

            uint64_t st2 = encode_state(nxt, b);
            if (!visited.insert(st2).second) continue;

            local_nodes.insert(nxt);
            q.emplace(nxt, b, static_cast<uint16_t>(d + 1));
        }

        if (early_all_met) break;
    }

    if (global_sink == UINT32_MAX) {
        debug_stream() << "        - " << vtx_name(v) << " -> no closed sink found\n";
        return bb;
    }

    const uint32_t sink = global_sink;

    // 9) Use BFS sink depth to limit DFS depth
    uint32_t eff_max_depth = max_depth_;

    if (sink_depth_from_src > 0) {
        uint32_t limit = sink_depth_from_src + DEPTH_MARGIN_;
        if (eff_max_depth == 0 || eff_max_depth > limit) {
            eff_max_depth = limit;
        }
    }

    debug_stream() << "  [bubble-final-closed] src=" << vtx_name(v)
                   << " sink=" << vtx_name(sink)
                   << "  bfs_sink_depth=" << sink_depth_from_src
                   << "  eff_max_depth=" << eff_max_depth << "\n";

    bool hit_limits = false;
    auto paths = graph_.enumerate_paths_DFS(v, sink, local_nodes, eff_max_depth, max_paths_, skip_comp_, hit_limits, dfs_guard_);

    debug_stream() << "  [bubble-final-closed] src=" << vtx_name(v)
                   << " sink=" << vtx_name(sink)
                   << "  paths=" << paths.size() << "\n\n";

    if (paths.empty()) return bb;

    bb.set_source(v);
    bb.set_sink(sink);
    bb.set_type(Type::Normal);
    bb.set_paths(std::move(paths));
    return bb;
}

void GfaBubbleFinder::find_open_bubbles_from_tips_(const std::vector<uint32_t>& tip_nodes)
{
    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (V == 0 || tip_nodes.size() < 2) return;

    auto vtx_name = [&](uint32_t x) -> std::string {
        uint32_t sid = NodeHandle::get_segment_id(x);
        bool rv = NodeHandle::get_is_reverse(x);
        return graph_.getNodeName(sid) + (rv ? "-" : "+");
    };
    auto sid_of = [&](uint32_t v) -> uint32_t { return NodeHandle::get_segment_id(v); };
    auto pair_key = [](uint32_t s0, uint32_t s1) -> uint64_t {
        if (s0 > s1) std::swap(s0, s1);
        return (static_cast<uint64_t>(s0) << 32) | static_cast<uint64_t>(s1);
    };
    auto rev = [](uint32_t x) { return x ^ 1u; };

    debug_stream() << "[open-tip] start detecting open bubbles from TIPs\n";

    // ---------- 1) back-BFS to identify all nodes connected to the Tip and distances ---------- //
    struct TipBackInfo {
        uint32_t tip;
        std::unordered_map<uint32_t, uint16_t> dist;
    };

    const unsigned threads = std::max(1u, thread_);
    ThreadPool pool_back(threads);

    std::vector<std::future<TipBackInfo>> futs;
    futs.reserve(tip_nodes.size());

    const uint16_t BACK_LIMIT = 1000;

    for (uint32_t t : tip_nodes) {
        futs.emplace_back(pool_back.submit([&, t]() -> TipBackInfo {
            TipBackInfo info;
            info.tip = t;

            std::queue<std::pair<uint32_t, uint16_t>> q;
            q.emplace(t, 0);
            info.dist[t] = 0;

            uint64_t states = 0;
            while (!q.empty()) {
                auto [cur, d] = q.front();
                q.pop();

                if (++states > dfs_guard_) break;
                if (d >= BACK_LIMIT) continue;

                const auto& ins = graph_.getArcsToVertex(cur);
                for (const GfaArc* a : ins) {
                    if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
                    uint32_t u = a->get_source_vertex_id();
                    if (u == cur) continue;

                    uint16_t nd = static_cast<uint16_t>(d + 1);
                    auto it = info.dist.find(u);
                    if (it != info.dist.end() && it->second <= nd) continue;

                    info.dist[u] = nd;
                    q.emplace(u, nd);
                }
            }

            debug_stream() << "  [TIP-BACK] tip=" << vtx_name(t) << " back_nodes=" << info.dist.size() << "\n";
            return info;
        }));
    }

    std::vector<TipBackInfo> back_infos;
    back_infos.reserve(tip_nodes.size());
    for (auto& f : futs) back_infos.push_back(f.get());
    pool_back.stop();

    const std::size_t N = back_infos.size();
    if (N < 2) return;

    const uint32_t OPEN_MAX_PATHS = (max_paths_ > 0) ? max_paths_ : 32;

    // ---------- 2) Find nearest neighbor tip ---------- //
    std::vector<int>      best_j(N, -1);
    std::vector<uint32_t> best_src(N, UINT32_MAX);
    std::vector<uint32_t> best_max(N, UINT32_MAX);
    std::vector<uint32_t> best_sum(N, UINT32_MAX);

    auto better = [](uint32_t mx1, uint32_t sm1, uint32_t mx2, uint32_t sm2) {
        return (mx1 < mx2) || (mx1 == mx2 && sm1 < sm2);
    };

    for (std::size_t i = 0; i < N; ++i) {
        const auto& bi = back_infos[i];

        for (std::size_t j = 0; j < N; ++j) {
            if (i == j) continue;

            if (is_complex_source_(bi.tip, {back_infos[j].tip})) continue;

            const auto& bj = back_infos[j];

            // iterate smaller map
            const auto* small = &bi.dist;
            const auto* large = &bj.dist;
            if (small->size() > large->size()) std::swap(small, large);

            uint32_t cand_src = UINT32_MAX;
            uint32_t cand_mx  = UINT32_MAX;
            uint32_t cand_sm  = UINT32_MAX;

            for (const auto& kv : *small) {
                uint32_t node = kv.first;
                uint16_t d1   = kv.second;
                auto it2 = large->find(node);
                if (it2 == large->end()) continue;
                uint16_t d2 = it2->second;

                uint32_t mx = std::max<uint32_t>(d1, d2);
                uint32_t sm = uint32_t(d1) + uint32_t(d2);

                if (better(mx, sm, cand_mx, cand_sm)) {
                    cand_mx  = mx;
                    cand_sm  = sm;
                    cand_src = node;
                }
            }

            if (cand_src == UINT32_MAX) continue;

            if (best_j[i] < 0 || better(cand_mx, cand_sm, best_max[i], best_sum[i])) {
                best_j[i]   = int(j);
                best_src[i] = cand_src;
                best_max[i] = cand_mx;
                best_sum[i] = cand_sm;
            }
        }
    }

    // ---------- 3) perform DFS only on pairs that are nearest neighbor to each other ---------- //
    ThreadPool pool_pairs(threads);

    struct PairResult {
        uint64_t key{0};
        Bubble   bubble;
    };

    std::vector<std::future<PairResult>> pair_futs;
    pair_futs.reserve(N);

    for (std::size_t i = 0; i < N; ++i) {
        int j = best_j[i];
        if (j < 0) continue;
        if (std::size_t(j) <= i) continue;
        if (best_j[j] != int(i)) continue;

        pair_futs.emplace_back(pool_pairs.submit([&, i, j]() -> PairResult {
            PairResult pr;

            const auto& bi = back_infos[i];
            const auto& bj = back_infos[j];

            const uint32_t tip_i = bi.tip;
            const uint32_t tip_j = bj.tip;

            pr.key = pair_key(sid_of(tip_i), sid_of(tip_j));

            const uint32_t mx = std::max(best_max[i], best_max[j]);
            const uint32_t eff_max_depth = mx + DEPTH_MARGIN_;
            const uint32_t common_src = best_src[i];

            debug_stream() << "  [open-tip-pair] tips=(" << vtx_name(tip_i) << ", " << vtx_name(tip_j)
                           << ") common_src=" << vtx_name(common_src)
                           << " best_max=" << mx << " best_sum=" << std::min(best_sum[i], best_sum[j])
                           << " depth=" << eff_max_depth << "\n";

            // Enumerate paths using DFS
            std::unordered_set<uint32_t> region;
            region.reserve(bi.dist.size() + bj.dist.size() + 8);
            for (const auto& kv : bi.dist) region.insert(rev(kv.first));
            for (const auto& kv : bj.dist) region.insert(rev(kv.first));
            region.insert(rev(common_src));
            region.insert(rev(tip_i));
            region.insert(rev(tip_j));

            bool hit_limits = false;

            auto paths_i = graph_.enumerate_paths_DFS(
                rev(tip_i), rev(common_src), region, eff_max_depth, OPEN_MAX_PATHS, skip_comp_, hit_limits, dfs_guard_
            );
            auto paths_j = graph_.enumerate_paths_DFS(
                rev(tip_j), rev(common_src), region, eff_max_depth, OPEN_MAX_PATHS, skip_comp_, hit_limits, dfs_guard_
            );

            if (hit_limits) mark_nodes_complex_({common_src, tip_i, tip_j});

            if (paths_i.empty() && paths_j.empty()) return pr;

            std::vector<std::vector<uint32_t>> all_paths;
            all_paths.reserve(paths_i.size() + paths_j.size());
            for (auto& p : paths_i) all_paths.emplace_back(std::move(p));
            for (auto& p : paths_j) all_paths.emplace_back(std::move(p));

            pr.bubble.set_source(common_src);
            pr.bubble.set_sink(common_src);
            pr.bubble.set_type(Type::Tip);
            pr.bubble.set_paths(std::move(all_paths));
            return pr;
        }));
    }

    std::vector<PairResult> pair_results;
    pair_results.reserve(pair_futs.size());
    for (auto& f : pair_futs) pair_results.push_back(f.get());
    pool_pairs.stop();

    // ---------- 4) Record the result ---------- //
    std::unordered_set<uint64_t> seen_tip_pairs;
    seen_tip_pairs.reserve(pair_results.size() * 2 + 8);

    for (auto& pr : pair_results) {
        if (pr.key == 0) continue;
        if (pr.bubble.get_paths().empty()) continue;
        if (!seen_tip_pairs.insert(pr.key).second) continue;
        bubbles_.emplace_back(std::move(pr.bubble));
    }

    debug_stream() << "[open-tip] done. bubbles now=" << bubbles_.size() << "\n";
}

// cluster paths by node-length symmetric-difference ratio, pick longest rep per cluster
std::vector<std::vector<uint32_t>> GfaBubbleFinder::pick_representative_paths(
    const std::vector<std::vector<uint32_t>>& paths
) const {
    std::vector<std::vector<uint32_t>> reps;
    if (paths.empty()) return reps;
    if (paths.size() == 1) { reps.push_back(paths[0]); return reps; }

    // 1) build per-path profile: sid -> total_bp, and total length
    std::vector<std::unordered_map<uint32_t, uint32_t>> prof(paths.size());
    std::vector<uint64_t> plen(paths.size(), 0);

    for (size_t i = 0; i < paths.size(); ++i) {
        auto& mp = prof[i];
        mp.reserve(paths[i].size() * 2 + 8);

        uint64_t L = 0;
        for (uint32_t vtx : paths[i]) {
            uint32_t sid = NodeHandle::get_segment_id(vtx);
            uint32_t w   = graph_.getNodeLength(sid);
            mp[sid] += w;
            L += w;
        }
        plen[i] = L;
    }

    // 2) sort indices by length desc (so first in cluster is naturally the rep = longest)
    std::vector<size_t> order(paths.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b){ return plen[a] > plen[b]; });

    // 3) greedy clustering on reps
    std::vector<size_t> rep_ids; // store index of representative path per cluster
    rep_ids.reserve(paths.size());

    auto diff_ratio = [&](size_t ia, size_t ib) -> double {
        const uint64_t A = plen[ia], B = plen[ib];
        if (A == 0 || B == 0) return 1.0;

        const auto& ma = prof[ia];
        const auto& mb = prof[ib];

        // intersection = sum(min(bpA[sid], bpB[sid])) over common sids
        uint64_t inter = 0;
        const auto& sm = (ma.size() <= mb.size()) ? ma : mb;
        const auto& lg = (ma.size() <= mb.size()) ? mb : ma;

        for (const auto& kv : sm) {
            auto it = lg.find(kv.first);
            if (it != lg.end()) inter += (uint64_t)std::min(kv.second, it->second);
        }

        const uint64_t uni = A + B - inter;
        return uni ? (1.0 - (double)inter / (double)uni) : 0.0;
    };

    for (size_t idx : order) {
        if (plen[idx] == 0) continue;

        bool placed = false;
        for (size_t r : rep_ids) {
            if (diff_ratio(idx, r) <= path_diff_) {
                placed = true; // belong to this cluster; rep already longest, no need update
                break;
            }
        }
        if (!placed) rep_ids.push_back(idx);
    }

    // 4) output representative paths (longest of each cluster)
    reps.reserve(rep_ids.size());
    for (size_t r : rep_ids) reps.push_back(paths[r]);
    return reps;
}

void GfaBubbleFinder::filter_nonlocal_bubbles_()
{
    log_stream() << "Filtering non-local bubbles ...\n";

    if (bubbles_.empty()) return;

    std::vector<Bubble> kept; 
    kept.reserve(bubbles_.size());

    for (const Bubble& bb : bubbles_) {
        // 1) Collect the set of vertices covered by this bubble, and the set of all sinks (the last vertex of each path)
        std::unordered_set<uint32_t> segids; segids.reserve(128);
        std::unordered_set<uint32_t> nodes; nodes.reserve(128);
        std::unordered_set<uint32_t> all_sinks; all_sinks.reserve(8);
        for (const auto& p : bb.get_paths()) {
            for (uint32_t v : p) {
                nodes.insert(v);
                segids.insert(NodeHandle::get_segment_id(v));
            }
            if (!p.empty()) { all_sinks.insert(p.back()); }
        }
        const uint32_t src  = bb.get_source();

        // 2) Check if any internal node has outgoing/incoming edges to outside the bubble
        bool leaky = false;
        for (uint32_t u : nodes) {
            if (u == src) continue;
            if (all_sinks.count(u)) continue;
            // 2a) outgoing edges: internal -> external
            for (const GfaArc* a : graph_.getArcsFromVertex(u)) {
                if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
                uint32_t w = a->get_target_vertex_id();
                if (graph_.getNodeDeleted(NodeHandle::get_segment_id(w))) continue;
                if (w == u) continue;
                if (!segids.count(NodeHandle::get_segment_id(w))) { leaky = true; break; }
            }
            if (leaky) break;

            // 2b) incoming edges: external -> internal
            for (const GfaArc* a : graph_.getArcsToVertex(u)) {
                if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
                uint32_t w = a->get_source_vertex_id();
                if (graph_.getNodeDeleted(NodeHandle::get_segment_id(w))) continue;
                if (w == u) continue;
                if (!segids.count(NodeHandle::get_segment_id(w))) { leaky = true; break; }
            }

            if (leaky) break;
        }

        if (leaky) { continue; }
        kept.emplace_back(bb);  // 'Locally independent' bubble
    }

    // Replace with the kept bubbles
    bubbles_.swap(kept);
}

void GfaBubbleFinder::save_bubble_as_gfa(
    const std::string& filename,
    const uint32_t min_len,
    const uint32_t min_num
) const {
    SAVE saver(filename);
    std::ostringstream oss;

    oss << "H\tVN:Z:1.0\tCreator:liftasm bubble\n";
    saver.save(oss.str()); oss.str(""); oss.clear();

    const std::size_t nb = bubbles_.size();

    std::vector<uint8_t> keep(nb, 0u);

    std::unordered_set<uint32_t> segs;
    segs.reserve(nb * 8);

    // Filter bubble and export H lines
    std::size_t out_idx = 0;
    for (std::size_t i = 0; i < nb; ++i) {
        const Bubble& bb = bubbles_[i];

        const uint32_t src_v    = bb.get_source();
        const uint32_t sink_v   = bb.get_sink();
        const uint32_t src_seg  = NodeHandle::get_segment_id(src_v);
        const uint32_t sink_seg = NodeHandle::get_segment_id(sink_v);

        std::unordered_set<uint32_t> inner_segs;
        inner_segs.reserve(16);

        for (const auto& p : bb.get_paths()) {
            for (uint32_t v : p) {
                const uint32_t sid = NodeHandle::get_segment_id(v);
                if (sid == src_seg || sid == sink_seg) continue;
                inner_segs.insert(sid);
            }
        }

        const uint32_t inner_nodes = static_cast<uint32_t>(inner_segs.size());

        uint64_t inner_bp = 0;
        for (uint32_t sid : inner_segs) {
            const GfaNode* node = graph_.getNode(sid);
            if (!node || graph_.getNodeDeleted(sid)) continue;
            inner_bp += graph_.getNodeLength(sid);
        }

        const bool pass_nodes = (min_num == 0 || inner_nodes >= min_num);
        const bool pass_bp    = (min_len == 0 || inner_bp    >= min_len);
        if (!pass_nodes || !pass_bp) continue;

        keep[i] = 1u;
        const std::size_t bid = out_idx++;

        const bool src_rev  = NodeHandle::get_is_reverse(src_v);
        const bool sink_rev = NodeHandle::get_is_reverse(sink_v);

        oss << "H\tbubble:i:" << bid
            << "\tsource:" << graph_.getNodeName(src_seg) << (src_rev ? "-" : "+")
            << "\tsink:"   << graph_.getNodeName(sink_seg) << (sink_rev ? "-" : "+")
            << "\tBT:Z:"   << bb.get_type() << "\tBL:i:"   << inner_bp << "\tBN:i:" << inner_nodes << "\n";
        saver.save(oss.str()); oss.str(""); oss.clear();

        for (const auto& p : bb.get_paths()) {
            for (uint32_t v : p) {
                segs.insert(NodeHandle::get_segment_id(v));
            }
        }
    }

    if (segs.empty()) return;

    // S lines
    std::vector<uint32_t> seg_vec(segs.begin(), segs.end());
    std::sort(seg_vec.begin(), seg_vec.end());

    for (uint32_t sid : seg_vec) {
        const GfaNode* node = graph_.getNode(sid);
        if (!node || graph_.getNodeDeleted(sid)) continue;

        oss << "S\t" << node->name << "\t";
        if (node->sequence.empty()) oss << "*";
        else oss << node->sequence;
        oss << "\tLN:i:" << graph_.getNodeLength(sid) << "\n";

        saver.save(oss.str()); oss.str(""); oss.clear();
    }

    // L lines
    for (const auto& arc : graph_.getAllArcs()) {
        if (arc.get_del() || arc.get_comp()) continue;

        const uint32_t v = arc.get_source_vertex_id();
        const uint32_t w = arc.get_target_vertex_id();

        const uint32_t sv = NodeHandle::get_segment_id(v);
        const uint32_t sw = NodeHandle::get_segment_id(w);
        if (!segs.count(sv) || !segs.count(sw)) continue;

        oss << "L\t" << graph_.getNodeName(sv) << "\t" << ((v & 1) ? '-' : '+') << "\t"
            << graph_.getNodeName(sw) << "\t" << ((w & 1) ? '-' : '+') << "\t"
            << graph_.format_overlap_field(arc.ov, arc.ow) << "\n";

        saver.save(oss.str()); oss.str(""); oss.clear();
    }

    // P lines
    for (std::size_t i = 0; i < nb; ++i) {
        if (!keep[i]) continue;

        const Bubble& bb = bubbles_[i];
        const std::size_t bid = i;

        const auto& paths = bb.get_paths();
        for (std::size_t pid = 0; pid < paths.size(); ++pid) {
            const auto& p = paths[pid];

            std::string seg_list_str;
            std::vector<uint32_t> vertex_list;
            seg_list_str.reserve(p.size() * 16);
            for (std::size_t k = 0; k < p.size(); ++k) {
                const uint32_t v = p[k];
                const uint32_t sid = NodeHandle::get_segment_id(v);
                if (sid == NodeHandle::get_segment_id(bb.get_source()) || sid == NodeHandle::get_segment_id(bb.get_sink())) continue;
                const bool rev = NodeHandle::get_is_reverse(v);
                seg_list_str.push_back(rev ? '<' : '>');
                seg_list_str += graph_.getNodeName(sid);
                vertex_list.push_back(v);
            }

            if (vertex_list.empty()) continue;

            const std::string seq = graph_.get_path_sequence(vertex_list);

            oss << "P\t"
                << "b" << bid << ".p" << pid << "\t"
                << seg_list_str << "\t*\t"
                << "BI:i:" << bid << "\t"
                << "PI:i:" << pid << "\t"
                << "PL:i:" << (seq.empty() ? 0 : seq.size()) << "\t"
                << "SEQ:Z:" << (seq.empty() ? "*" : seq)
                << "\n";

            saver.save(oss.str()); oss.str(""); oss.clear();
        }
    }
}

void GfaBubbleFinder::print_bubbles() const
{
    for (const auto& bb : bubbles_) {
        uint32_t s_seg = NodeHandle::get_segment_id(bb.get_source());
        uint32_t t_seg = NodeHandle::get_segment_id(bb.get_sink());
        bool     s_rev = NodeHandle::get_is_reverse(bb.get_source());
        bool     t_rev = NodeHandle::get_is_reverse(bb.get_sink());

        const std::string src  = graph_.getNodeName(s_seg) + (s_rev? "-":"+");
        const std::string sink = graph_.getNodeName(t_seg) + (t_rev? "-":"+");

        std::unordered_set<uint32_t> inner_segs;
        inner_segs.reserve(16);

        std::size_t total_bp = 0;
        for (const auto& p : bb.get_paths()) {
            for (uint32_t vtx : p) {
                uint32_t seg = NodeHandle::get_segment_id(vtx);
                if (seg == s_seg || seg == t_seg) continue;
                if (!inner_segs.insert(seg).second) continue;

                const GfaNode* node = graph_.getNode(seg);
                if (!node || graph_.getNodeDeleted(seg)) continue;
                total_bp += graph_.getNodeLength(seg);
            }
        }

        log_stream() << "Bubble: source=" << src << " sink=" << sink << "\n";
        log_stream() << "   - Total Nodes Length: " << total_bp << "\n";
        log_stream() << "   - Total Nodes Count: " << inner_segs.size() << "\n";
        log_stream() << "   - type: " << bb.get_type() << "\n";

        const auto& ps = bb.get_paths();
        std::string path_str;
        for (size_t i = 0; i < ps.size(); ++i) {
            path_str = "   - path[" + std::to_string(i) + "]: ";
            for (uint32_t vtx : ps[i]) {
                uint32_t seg = NodeHandle::get_segment_id(vtx);
                bool     rev = NodeHandle::get_is_reverse(vtx);
                path_str += graph_.getNodeName(seg) + (rev? "-":"+") + ' ';
            }
            log_stream() << path_str << "\n";
        }
    }
}

void GfaBubbleFinder::find_forks()
{
    log_stream() << "Detecting forks ...\n";
    forks_.clear();

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (V == 0) return;

    // used to avoid duplicate forks (by hash)
    std::unordered_set<uint64_t> seen; 
    seen.reserve(1024);

    for (uint32_t v = 0; v < V; ++v) {
        std::vector<uint32_t> children;
        children.reserve(8);
        std::unordered_set<std::uint32_t> seen_sid; seen_sid.reserve(8);

        for (const GfaArc* a : graph_.getArcsFromVertex(v)) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
            uint32_t w = a->get_target_vertex_id();
            if (w == v) continue;

            // deduplicate by segment name (2025-10-12)
            const uint32_t sid_w = NodeHandle::get_segment_id(w);
            if (!seen_sid.insert(sid_w).second) continue;

            children.push_back(w);
        }
        if (children.size() < 2) continue;  // Less than 2 successors, cannot form a fork group

        // 2025-10-24
        std::sort(children.begin(), children.end());
        const uint64_t key = hash_vec_(children, /*use_direction=*/false);
        if (!seen.insert(key).second) { continue; }  // duplication fork
        ForkGroup g(v, std::move(children));
        forks_.push_back(std::move(g));

        // // For each successor w, check if its filtered incoming edges are unique and the only predecessor is v
        // std::vector<uint32_t> good;
        // good.reserve(children.size());

        // for (uint32_t w : children) {
        //     uint32_t indeg = 0;
        //     bool only_from_v = true;

        //     for (const GfaArc* in : graph_.getArcsToVertex(w)) {
        //         if (!in || in->get_del() || (skip_comp_ && in->get_comp())) continue;
        //         uint32_t u = in->get_source_vertex_id();
        //         if (u == w) continue;
        //         ++indeg;
        //         if (u != v) { only_from_v = false; break; }
        //         if (indeg > 1) break;
        //     }

        //     if (indeg == 1 && only_from_v) {
        //         good.push_back(w);
        //     }
        // }

        // if (good.size() >= 2) {
        //     std::sort(good.begin(), good.end());
        //     const uint64_t key = hash_vec_(good, /*use_direction=*/false);
        //     if (!seen.insert(key).second) { continue; }  // duplication fork

        //     ForkGroup g(v, std::move(good));
        //     forks_.push_back(std::move(g));
        // }
    }
    log_stream() << "   - Total forks detected: " << forks_.size() << "\n";
}

void GfaBubbleFinder::print_forks() const
{
    for (const auto& g : forks_) {
        uint32_t s_seg = NodeHandle::get_segment_id(g.get_source());
        bool     s_rev = NodeHandle::get_is_reverse(g.get_source());
        std::string s_name = graph_.getNodeName(s_seg) + (s_rev ? "-" : "+");

        log_stream() << "Fork: source=" << s_name << "  branches=" << g.size() << "\n";
        for (size_t i = 0; i < g.size(); ++i) {
            uint32_t w = g.get_branches()[i];
            uint32_t seg = NodeHandle::get_segment_id(w);
            bool     rev = NodeHandle::get_is_reverse(w);
            log_stream() << "   - leaf[" << i << "]: " << graph_.getNodeName(seg) << (rev ? "-" : "+") << "\n";
        }
    }
}

} // namespace GfaBubble