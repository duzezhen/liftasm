#include "../include/gfa_bubble.hpp"
#include "../include/progress_tracker.hpp"
#include "../include/ThreadPool.hpp"

#include <numeric>

namespace GfaBubble {

GfaBubbleFinder::GfaBubbleFinder(
    const GfaGraph& g, std::size_t md, std::size_t mp, uint64_t dfs_guard, 
        uint16_t cx_depth, uint32_t cx_nodes, uint32_t cx_branches, int cx_deg_branch, int cx_deg_hub, 
        double path_diff, uint32_t stall_round_limit, bool sc, bool keep_nested,
        double same_sim, uint32_t same_min_num, uint32_t same_min_len,
        uint32_t diff_min_src, double diff_sim, uint32_t diff_min_num, uint32_t diff_min_len, 
        uint32_t homo_num, uint32_t homo_k, uint32_t homo_w, 
        uint32_t homo_boot_bp, uint32_t homo_cnt_len, uint64_t homo_bloom_bits, uint32_t homo_bloom_hash,
        uint32_t thread
    )
    : graph_(g), max_depth_(md), max_paths_(mp), dfs_guard_(dfs_guard), 
    cx_depth_(cx_depth), cx_nodes_(cx_nodes), cx_branches_(cx_branches), cx_deg_branch_(cx_deg_branch), cx_deg_hub_(cx_deg_hub), 
    path_diff_(path_diff), stall_round_limit_(stall_round_limit), skip_comp_(sc), keep_nested_(keep_nested), 
    same_sim_(same_sim), same_min_num_(same_min_num), same_min_len_(same_min_len),
    diff_min_src_(diff_min_src), diff_sim_(diff_sim), diff_min_num_(diff_min_num), diff_min_len_(diff_min_len), homo_num_(std::max<uint32_t>(1, homo_num)),
    homo_boot_bp_(homo_boot_bp), homo_cnt_len_(homo_cnt_len), homo_bloom_bits_(homo_bloom_bits), homo_bloom_hash_(homo_bloom_hash),
    mm_opt_{/*k=*/homo_k, /*w=*/homo_w, /*seek_reverse=*/false, /*seed=*/0x8a5cd789635d2dffULL},
    overlap_len_{homo_w > 0 ? homo_w - 1 : 0}, 
    thread_(thread)
{
    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    reach_buf_.assign(V, UNVIS_FLAG_);
    complex_flag_.assign(V, UNKNOWN_FLAG_);
}

uint64_t GfaBubbleFinder::hash_vec_(std::vector<uint32_t> vertexs, bool use_direction) const {
    if (!use_direction) {
        for (auto &x : vertexs) x >>= 1;
    }
    std::sort(vertexs.begin(), vertexs.end());
    vertexs.erase(std::unique(vertexs.begin(), vertexs.end()), vertexs.end());

    uint64_t h = 14695981039346656037ULL;
    for (uint32_t v : vertexs) {
        h ^= static_cast<uint64_t>(v);
        h *= 1099511628211ULL;
    }
    return h;
}

bool GfaBubbleFinder::is_single_node_(Vertex v) const
{
    for (const GfaArc* a : graph_.getArcsFromVertex(v.vertex_id())) {
        if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;

        uint32_t t = a->get_target_vertex_id();
        if (t != v.vertex_id()) return false;
    }

    for (const GfaArc* a : graph_.getArcsToVertex(v.vertex_id())) {
        if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;

        uint32_t s = a->get_source_vertex_id();
        if (s != v.vertex_id()) return false;
    }

    return true;
}

bool GfaBubbleFinder::is_strict_source_(Vertex vtx, std::vector<Vertex>& tgt) const
{
    tgt.clear();
    tgt.reserve(8);

    std::unordered_set<uint32_t> seen_seg; seen_seg.reserve(8);

    for (const GfaArc* a : graph_.getArcsFromVertex(vtx.vertex_id())) {
        if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
        Vertex w = Vertex(a->get_target_vertex_id());
        uint32_t seg = w.segment_id();
        if (!seen_seg.insert(seg).second) continue;
        if (w == vtx) continue;
        tgt.push_back(w);
    }

    return tgt.size() >= 2;
}

bool GfaBubbleFinder::is_strict_sink_(Vertex v, Vertex w) const {
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

    for (const GfaArc* in : graph_.getArcsToVertex(w.vertex_id())) {
        if (!in || in->get_del() || (skip_comp_ && in->get_comp())) continue;
        const uint32_t pv = in->get_source_vertex_id();
        if (pv == w.vertex_id() || pv == v.vertex_id()) continue;

        int outdeg = 0;
        bool all_funnel_to_w = true;

        for (const GfaArc* a : graph_.getArcsFromVertex(pv)) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
            uint32_t x = a->get_target_vertex_id();
            if (x == pv) continue;
            ++outdeg;

            if (!leads_uniquely_to(x, w.vertex_id(), MAX_TUNNEL_STEPS)) {
                all_funnel_to_w = false;
                if (outdeg > 1) return false;
            }
        }

        if (outdeg > 1 && !all_funnel_to_w) return false;
    }

    return true;
}

void GfaBubbleFinder::mark_nodes_complex_(const std::vector<Vertex>& nodes)
{
    const uint32_t V = static_cast<uint32_t>(complex_flag_.size());
    if (V == 0) return;

    std::lock_guard<std::mutex> lk(complex_mtx_);
    for (const Vertex& x : nodes) {
        if (x.vertex_id() >= V) continue;
        complex_flag_[x.vertex_id()] = COMPLEX_FLAG_;
    }
}
void GfaBubbleFinder::mark_nodes_simple_(const std::vector<Vertex>& nodes)
{
    const uint32_t V = static_cast<uint32_t>(complex_flag_.size());
    if (V == 0) return;

    std::lock_guard<std::mutex> lk(complex_mtx_);
    for (const Vertex& x : nodes) {
        if (x.vertex_id() >= V) continue;
        if (complex_flag_[x.vertex_id()] == UNKNOWN_FLAG_) {
            complex_flag_[x.vertex_id()] = SIMPLE_FLAG_;
        }
    }
}

bool GfaBubbleFinder::is_complex_source_(
    Vertex v,
    const std::vector<Vertex>& seed_buf
) {
    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (complex_flag_.size() != V) return false;

    if (cx_depth_ == 0) return false;

    auto get_flag = [&](Vertex x) -> uint8_t {
        if (x.vertex_id() >= V) return UNKNOWN_FLAG_;
        std::lock_guard<std::mutex> lk(complex_mtx_);
        return complex_flag_[x.vertex_id()];
    };
    auto is_marked_complex = [&](Vertex x) -> bool {
        return get_flag(x) == COMPLEX_FLAG_;
    };
    auto is_marked_simple = [&](Vertex x) -> bool {
        return get_flag(x) == SIMPLE_FLAG_;
    };

    // Count the total degreeof u
    auto total_deg = [&](Vertex u) -> int {
        int deg = 0;

        const auto& out_arcs = graph_.getArcsFromVertex(u.vertex_id());
        for (const GfaArc* a : out_arcs) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
            uint32_t w = a->get_target_vertex_id();
            if (w == u.vertex_id()) continue;
            ++deg;
        }

        return deg;
    };

    // 0. Use cached flags
    if (is_marked_complex(v)) return true;

    bool any_seed_complex = false;
    bool all_seed_simple  = !seed_buf.empty();
    for (const Vertex& s : seed_buf) {
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

    // 1. src or any seed with total_deg>=DEG_HUB is complex
    {
        if (total_deg(v) >= cx_deg_hub_) {
            mark_nodes_complex_({v});
            return true;
        }
        for (const Vertex& s : seed_buf) {
            if (total_deg(s) >= cx_deg_hub_) {
                mark_nodes_complex_({v, s});
                return true;
            }
        }
    }

    // 2. Local BFS: count nodes / branching nodes
    std::queue<std::pair<Vertex,uint16_t>> q; // (vertex, depth)
    std::unordered_set<Vertex> visited;
    visited.reserve(256);

    std::vector<Vertex> local_nodes;
    local_nodes.reserve(256);

    auto push_node = [&](Vertex x, uint16_t d) {
        if (!visited.insert(x).second) return;
        q.emplace(x, d);
        local_nodes.push_back(x);
    };

    push_node(v, 0);
    for (const Vertex& s : seed_buf) {
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
        const auto& arcs = graph_.getArcsFromVertex(u.vertex_id());
        for (const GfaArc* a : arcs) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;
            uint32_t w = a->get_target_vertex_id();
            if (w == u.vertex_id()) continue;

            if (is_marked_complex(Vertex(w))) {
                local_nodes.push_back(Vertex(w));
                mark_nodes_complex_(local_nodes);
                return true;
            }

            push_node(Vertex(w), static_cast<uint16_t>(d + 1));
        }
    }

    mark_nodes_simple_(local_nodes);
    return false;
}

void GfaBubbleFinder::collect_strict_sources_tips_()
{
    sources_.clear();
    tips_.clear();

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (V == 0) return;

    sources_.reserve(V / 16 + 8);
    tips_.reserve(V / 32 + 4);

    std::vector<Vertex> target_vtxs;
    target_vtxs.reserve(16);

    for (uint32_t v = 0; v < V; ++v) {
        if (graph_.getNodeDeleted(NodeHandle::get_segment_id(v))) continue;

        // Single-node (Chromosome-like)
        if (is_single_node_(Vertex(v))) {
            sources_.emplace_back(Vertex(v));
            continue;
        }

        target_vtxs.clear();
        if (!is_strict_source_(Vertex(v), target_vtxs)) {
            if (target_vtxs.empty()) {
                tips_.emplace_back(Vertex(v));
            }
            continue;
        }
        if (target_vtxs.size() < 2 || target_vtxs.size() > 64) continue;
        if (is_complex_source_(Vertex(v), target_vtxs)) continue;

        sources_.emplace_back(Vertex(v));
    }

    log_stream() << "  - strict sources: " << sources_.size() << ", tips: " << tips_.size() << "\n";

    return;
}

void GfaBubbleFinder::find_bubbles()
{
    log_stream() << "Detecting bubbles ...\n";
    bubbles_.clear();

    // Identify sources and tips
    if (sources_.empty() || tips_.empty()) {
        collect_strict_sources_tips_();
    }
    if (sources_.empty() || tips_.empty()) {
        warning_stream() << "  ! No strict sources or tips found, skip bubble detection.\n";
        return;
    }

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (V == 0) return;

    const uint64_t bfs_limit = static_cast<uint64_t>(std::min<std::size_t>(max_depth_ ? max_depth_ : 64, UINT64_MAX));

    // Closed bubbles from strict sources
    auto pair_key = [](uint32_t s0, uint32_t s1) -> uint64_t {
        if (s0 > s1) std::swap(s0, s1);
        return (static_cast<uint64_t>(s0) << 32) | static_cast<uint64_t>(s1);
    };

    const unsigned threads = std::max(1u, thread_);
    ThreadPool pool(threads);

    std::vector<std::future<Bubble>> futs;
    futs.reserve(sources_.size());

    ProgressTracker prog(sources_.size());

    for (Vertex src : sources_) {
        prog.hit();
        futs.emplace_back(
            pool.submit([&, src]() -> Bubble {
                return detect_closed_bubble_from_source_(src, bfs_limit);
            })
        );
    }

    std::unordered_set<uint64_t> seen_pair;
    seen_pair.reserve(2048);

    for (auto &f : futs) {
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

    log_stream() << "  - Total bubbles detected: " << bubbles_.size() << "\n" << "\n";

    // Post-processing
    if (!keep_nested_) filter_nonlocal_bubbles_();

    std::sort(bubbles_.begin(), bubbles_.end(),
        [](const Bubble& a, const Bubble& b) {
            if (a.get_len() != b.get_len()) return a.get_len() > b.get_len();
            return a.get_paths().size() > b.get_paths().size();
        }
    );
}

bool GfaBubbleFinder::find_common_nearest_sink_(
    Vertex source,
    const std::vector<Vertex>& seeds,
    uint64_t bfs_limit,
    uint32_t& best_sink,
    uint32_t& best_sink_depth,
    std::unordered_set<uint32_t>& local_nodes
) {
    auto vtx_name = [&](Vertex x) -> std::string {
        uint32_t sid = x.segment_id();
        bool     rv  = x.is_reverse();
        return graph_.getNodeName(sid) + (rv ? "-" : "+");
    };

    auto log_indent = [&](int level) -> std::string {
        return std::string(level * 2, ' ');
    };

    best_sink = UINT32_MAX;
    best_sink_depth = 0;

    const int16_t seeds_cnt = static_cast<int16_t>(seeds.size());
    if (seeds_cnt < 2) {
        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "seed count < 2 in sink finder, skip\n";
        }
        return false;
    }

    struct SinkStat {
        uint32_t seen_count = 0;
        uint32_t sum_depth  = 0;
        uint32_t max_depth  = 0;
    };

    std::queue<std::tuple<uint32_t, int16_t, uint16_t>> q;
    std::unordered_set<uint64_t> visited;
    visited.reserve(static_cast<std::size_t>(seeds_cnt) * 32);

    std::vector<std::unordered_set<uint32_t>> branch_sink_seen(static_cast<std::size_t>(seeds_cnt));
    std::unordered_map<uint32_t, SinkStat> sink_stats;
    sink_stats.reserve(128);

    std::vector<uint32_t> branch_active(static_cast<std::size_t>(seeds_cnt), 0);  // 0-inactive, 1-active

    auto encode_state = [](uint32_t node_id, int16_t br) -> uint64_t {
        return (static_cast<uint64_t>(node_id) << 8) | static_cast<uint64_t>(static_cast<uint8_t>(br));
    };

    uint32_t best_max_depth = UINT32_MAX;
    uint32_t best_sum_depth = UINT32_MAX;

    for (int16_t b = 0; b < seeds_cnt; ++b) {
        const Vertex w = seeds[static_cast<std::size_t>(b)];
        visited.insert(encode_state(w.vertex_id(), b));
        q.emplace(w.vertex_id(), b, 0u);
        ++branch_active[static_cast<std::size_t>(b)];
        local_nodes.insert(w.vertex_id());

        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "seed " << b << ": " << vtx_name(w) << "\n";
        }

        if (w.vertex_id() != source.vertex_id() && is_strict_sink_(source, w)) {
            auto& seen_set = branch_sink_seen[static_cast<std::size_t>(b)];
            if (seen_set.insert(w.vertex_id()).second) {
                SinkStat& ss = sink_stats[w.vertex_id()];
                ++ss.seen_count;

                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3)
                                   << "candidate sink: " << vtx_name(w)
                                   << "  depth=0"
                                   << "  branch_seen=" << ss.seen_count << "/" << seeds_cnt
                                   << "  max_depth=" << ss.max_depth
                                   << "  sum_depth=" << ss.sum_depth << "\n";
                }

                if (ss.seen_count == static_cast<uint32_t>(seeds_cnt)) {
                    bool better = false;
                    if (best_sink == UINT32_MAX) {
                        better = true;
                    } else if (ss.max_depth < best_max_depth) {
                        better = true;
                    } else if (ss.max_depth == best_max_depth && ss.sum_depth < best_sum_depth) {
                        better = true;
                    }

                    if (better) {
                        best_sink = w.vertex_id();
                        best_max_depth = ss.max_depth;
                        best_sum_depth = ss.sum_depth;
                        best_sink_depth = ss.max_depth;

                        if (DEBUG_ENABLED) {
                            debug_stream() << log_indent(4)
                                           << "update best sink: " << vtx_name(w)
                                           << "  best_max_depth=" << best_max_depth
                                           << "  best_sum_depth=" << best_sum_depth << "\n";
                        }
                    }
                }
            }
        }
    }

    while (!q.empty()) {
        auto [cur, b, d] = q.front();
        q.pop();

        --branch_active[static_cast<std::size_t>(b)];
        bool branch_extended = false;

        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2)
                           << "visit "
                           << vtx_name(Vertex(cur))
                           << "  depth=" << d
                           << "  branch=" << b
                           << "  qsize=" << q.size() << "\n";
        }

        if (cur == source.vertex_id()) {
            if (best_sink == UINT32_MAX && branch_active[static_cast<std::size_t>(b)] == 0) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3) << "stop bfs: branch " << b << " returned to source and has no active nodes\n";
                }
                break;
            }
            continue;
        }

        if (d >= bfs_limit) {
            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(3) << "reach bfs_limit at " << vtx_name(Vertex(cur)) << "\n";
            }
            break;
        }

        if (best_sink != UINT32_MAX && static_cast<uint32_t>(d) > best_sink_depth) {
            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(3) << "stop bfs early at depth=" << d << " because current best sink depth=" << best_sink_depth << "\n";
            }
            break;
        }

        const auto& arcs = graph_.getArcsFromVertex(cur);
        for (const GfaArc* a : arcs) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;

            const uint32_t nxt = a->get_target_vertex_id();
            if (nxt == cur) continue;

            const uint64_t st = encode_state(nxt, b);
            if (!visited.insert(st).second) continue;

            local_nodes.insert(nxt);
            q.emplace(nxt, b, static_cast<uint16_t>(d + 1));
            ++branch_active[static_cast<std::size_t>(b)];
            branch_extended = true;

            if (nxt == source.vertex_id()) continue;
            if (!is_strict_sink_(source, Vertex(nxt))) continue;

            auto& seen_set = branch_sink_seen[static_cast<std::size_t>(b)];
            if (!seen_set.insert(nxt).second) continue;

            SinkStat& ss = sink_stats[nxt];
            ++ss.seen_count;
            ss.sum_depth += static_cast<uint32_t>(d + 1);
            ss.max_depth = std::max(ss.max_depth, static_cast<uint32_t>(d + 1));

            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(3)
                               << "candidate sink: " << vtx_name(Vertex(nxt))
                               << "  depth=" << (d + 1)
                               << "  branch_seen=" << ss.seen_count << "/" << seeds_cnt
                               << "  max_depth=" << ss.max_depth
                               << "  sum_depth=" << ss.sum_depth << "\n";
            }

            if (ss.seen_count == static_cast<uint32_t>(seeds_cnt)) {
                bool better = false;
                if (best_sink == UINT32_MAX) {
                    better = true;
                } else if (ss.max_depth < best_max_depth) {
                    better = true;
                } else if (ss.max_depth == best_max_depth && ss.sum_depth < best_sum_depth) {
                    better = true;
                }

                if (better) {
                    best_sink = nxt;
                    best_max_depth = ss.max_depth;
                    best_sum_depth = ss.sum_depth;
                    best_sink_depth = ss.max_depth;

                    if (DEBUG_ENABLED) {
                        debug_stream() << log_indent(4) << "update best sink: " << vtx_name(Vertex(nxt)) << "  best_max_depth=" << best_max_depth << "  best_sum_depth=" << best_sum_depth << "\n";
                    }
                }
            }
        }

        if (best_sink == UINT32_MAX && !branch_extended && branch_active[static_cast<std::size_t>(b)] == 0) {
            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(3) << "stop bfs: branch " << b << " has no active nodes and no common sink found\n";
            }
            break;
        }
    }

    if (DEBUG_ENABLED && best_sink != UINT32_MAX) {
        debug_stream() << log_indent(2) << "selected sink: " << vtx_name(Vertex(best_sink)) << "  sink_depth=" << best_sink_depth << "\n";
    }

    return best_sink != UINT32_MAX;
}

Bubble GfaBubbleFinder::detect_closed_bubble_from_source_(Vertex v, uint64_t bfs_limit)
{
    Bubble bb;

    auto vtx_name = [&](Vertex x) -> std::string {
        uint32_t sid = x.segment_id();
        bool     rv  = x.is_reverse();
        return graph_.getNodeName(sid) + (rv ? "-" : "+");
    };

    auto bubble_len = [&](const Bubble& bb) -> uint64_t {
        const uint32_t src_seg  = Vertex::get_segment_id(bb.get_source());
        const uint32_t sink_seg = Vertex::get_segment_id(bb.get_sink());

        std::unordered_set<uint32_t> inner_segs;
        inner_segs.reserve(128);

        for (const auto& p : bb.get_paths()) {
            for (uint32_t v : p) {
                const uint32_t sid = Vertex::get_segment_id(v);
                if (sid == src_seg || sid == sink_seg) continue;
                inner_segs.insert(sid);
            }
        }

        uint64_t len = 0;
        for (uint32_t sid : inner_segs) {
            const GfaNode* node = graph_.getNode(sid);
            if (!node || graph_.getNodeDeleted(sid)) continue;
            len += graph_.getNodeLength(sid);
        }
        return len;
    };

    auto log_indent = [&](int level) -> std::string {
        return std::string(level * 2, ' ');
    };

    if (DEBUG_ENABLED) {
        debug_stream() << "Detect closed bubble\n";
        debug_stream() << log_indent(1) << "source    : " << vtx_name(v) << "\n";
        debug_stream() << log_indent(1) << "bfs_limit : " << bfs_limit << "\n";
    }

    // 1. Collect seeds
    if (DEBUG_ENABLED) {
        debug_stream() << log_indent(1) << "seed collection\n";
    }

    std::vector<Vertex> seed_buf;
    seed_buf.reserve(16);

    if (!is_strict_source_(v, seed_buf)) return bb;

    const int16_t seeds_cnt = static_cast<int16_t>(seed_buf.size());
    if (seeds_cnt < 2) {
        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "seed count < 2, skip\n";
        }
        return bb;
    }

    if (DEBUG_ENABLED) {
        debug_stream() << log_indent(2) << "seed count: " << seeds_cnt << "\n";
    }

    // 2. Complexity filter
    if (is_complex_source_(v, seed_buf)) {
        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(1) << "complexity check\n";
            debug_stream() << log_indent(2) << "local region too complex, skip\n";
        }
        return bb;
    }

    if (DEBUG_ENABLED) {
        debug_stream() << log_indent(1) << "complexity check\n";
        debug_stream() << log_indent(2) << "passed\n";
    }

    // 3. Collect sink candidates (BFS)
    std::unordered_set<uint32_t> local_nodes;
    local_nodes.reserve(256);

    uint32_t sink = UINT32_MAX;
    uint32_t sink_depth_from_src = 0;

    if (DEBUG_ENABLED) {
        debug_stream() << log_indent(1) << "bfs to detect sink\n";
    }

    if (!find_common_nearest_sink_(v, seed_buf, bfs_limit, sink, sink_depth_from_src, local_nodes)) {
        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "no closed sink found\n";
        }
        return bb;
    }

    // 4. Use BFS sink depth to limit DFS depth
    uint32_t eff_max_depth = max_depth_;
    if (sink_depth_from_src > 0) {
        uint32_t limit = sink_depth_from_src + DEPTH_MARGIN_;
        if (eff_max_depth == 0 || eff_max_depth > limit) {
            eff_max_depth = limit;
        }
    }

    if (DEBUG_ENABLED) {
        debug_stream() << log_indent(1) << "finalize bubble region\n";
        debug_stream() << log_indent(2) << "sink           : " << vtx_name(Vertex(sink)) << "\n";
        debug_stream() << log_indent(2) << "bfs sink depth : " << sink_depth_from_src << "\n";
        debug_stream() << log_indent(2) << "eff max depth  : " << eff_max_depth << "\n";
        debug_stream() << log_indent(2) << "local nodes    : " << local_nodes.size() << "\n";
    }

    // 5. DFS to enumerate paths from source to sink within local nodes and depth limit
    bool hit_limits = false;
    auto paths = graph_.enumerate_paths_greedy_DFS(
        v.vertex_id(), sink, local_nodes, eff_max_depth,
        max_paths_, skip_comp_, hit_limits, dfs_guard_, stall_round_limit_
    );

    if (DEBUG_ENABLED) {
        debug_stream() << log_indent(1) << "enumerate paths\n";
        debug_stream() << log_indent(2) << "path count : " << paths.size() << "\n";
        debug_stream() << log_indent(2) << "hit limits : " << (hit_limits ? "true" : "false") << "\n";
        debug_stream() << "\n";
    }

    if (paths.empty()) return bb;

    bb.set_source(v.vertex_id());
    bb.set_sink(sink);
    bb.set_type(Type::Normal);
    bb.set_paths(std::move(paths));
    bb.set_len(bubble_len(bb));
    return bb;
}

// cluster paths by node-length symmetric-difference ratio, pick longest rep per cluster
std::vector<std::vector<uint32_t>> GfaBubbleFinder::pick_representative_paths(
    const std::vector<std::vector<uint32_t>>& paths
) const {
    std::vector<std::vector<uint32_t>> reps;
    if (paths.empty()) return reps;
    if (paths.size() == 1) { reps.push_back(paths[0]); return reps; }

    // 1. build per-path profile: sid -> total_bp, and total length
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

    // 2. sort indices by length desc (so first in cluster is naturally the rep = longest)
    std::vector<size_t> order(paths.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b){ return plen[a] > plen[b]; });

    // 3. greedy clustering on reps
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

    // 4. output representative paths (longest of each cluster)
    reps.reserve(rep_ids.size());
    for (size_t r : rep_ids) reps.push_back(paths[r]);
    return reps;
}

std::vector<std::vector<uint32_t>> GfaBubbleFinder::pick_representative_paths(
    const std::vector<std::vector<Vertex>>& paths
) const {
    std::vector<std::vector<uint32_t>> conv;
    conv.reserve(paths.size());
    for (const auto& p : paths) {
        std::vector<uint32_t> vp;
        vp.reserve(p.size());
        for (const Vertex& v : p) vp.push_back(v.vertex_id());
        conv.push_back(std::move(vp));
    }
    return pick_representative_paths(conv);
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

    log_stream() << "  - Kept " << bubbles_.size() << " local bubbles" << "\n" << "\n";
}

std::unordered_set<uint64_t> GfaBubbleFinder::build_bubble_branch_pairs_() const
{
    log_stream() << "Building bubble branch pairs ...\n";

    std::unordered_set<uint64_t> used;
    used.reserve(bubbles_.size() * 4 + 16);

    auto add_group = [&](std::vector<uint32_t> vs) {
        vs.erase(std::remove(vs.begin(), vs.end(), UINT32_MAX), vs.end());
        if (vs.size() < 2) return;
        used.insert(hash_vec_(vs, /*use_direction=*/false));
    };

    auto add_pair = [&](uint32_t a, uint32_t b) {
        if (a == UINT32_MAX || b == UINT32_MAX) return;
        if (Vertex::get_segment_id(a) == Vertex::get_segment_id(b)) return;

        used.insert(hash_vec_(
            std::vector<uint32_t>{a, b},
            /*use_direction=*/false
        ));
    };

    for (const auto& bb : bubbles_) {
        const uint32_t first_vtx = bb.get_source();
        const uint32_t last_vtx  = bb.get_sink();

        std::vector<uint32_t> src_group;
        std::vector<uint32_t> sink_group;

        std::vector<uint32_t> all_nodes;

        src_group.reserve(bb.get_paths().size() + 1);
        sink_group.reserve(bb.get_paths().size() + 1);
        all_nodes.reserve(bb.get_paths().size() * 4 + 2);

        src_group.push_back(first_vtx);
        sink_group.push_back(last_vtx);

        std::unordered_set<uint32_t> seen_src;
        std::unordered_set<uint32_t> seen_sink;
        std::unordered_set<uint32_t> seen_all;

        seen_src.reserve(bb.get_paths().size() * 2 + 4);
        seen_sink.reserve(bb.get_paths().size() * 2 + 4);
        seen_all.reserve(bb.get_paths().size() * 8 + 8);

        auto add_all_node = [&](uint32_t vtx) {
            if (vtx == UINT32_MAX) return;
            if (seen_all.insert(Vertex::get_segment_id(vtx)).second) {
                all_nodes.push_back(vtx);
            }
        };

        add_all_node(first_vtx);
        add_all_node(last_vtx);

        for (const auto& p : bb.get_paths()) {
            if (p.empty()) continue;

            for (uint32_t vtx : p) {
                add_all_node(vtx);
            }

            uint32_t second_vtx = UINT32_MAX;
            for (uint32_t vtx : p) {
                if (vtx == first_vtx) continue;
                second_vtx = vtx;
                break;
            }

            uint32_t prev_to_last_vtx = UINT32_MAX;
            for (auto it = p.rbegin(); it != p.rend(); ++it) {
                if (*it == last_vtx) continue;
                prev_to_last_vtx = *it;
                break;
            }

            if (second_vtx != UINT32_MAX && seen_src.insert(second_vtx).second) {
                src_group.push_back(second_vtx);
            }
            if (prev_to_last_vtx != UINT32_MAX && seen_sink.insert(prev_to_last_vtx).second) {
                sink_group.push_back(prev_to_last_vtx);
            }
        }

        add_group(std::move(src_group));
        add_group(std::move(sink_group));

        for (size_t i = 0; i < all_nodes.size(); ++i) {
            for (size_t j = i + 1; j < all_nodes.size(); ++j) {
                add_pair(all_nodes[i], all_nodes[j]);
            }
        }
    }

    log_stream() << "  - Built " << used.size() << " bubble branch pairs\n" << "\n";

    return used;
}


void GfaBubbleFinder::save_bubble_as_gfa(
    const std::string& output_file,
    const uint32_t min_len,
    const uint32_t min_num, 
    bool write_seq
) const {
    log_stream() << "Writing bubble to " << output_file << " ...\n";

    const std::size_t nb = bubbles_.size();
    const std::size_t nh = homologous_paths_.size();

    std::vector<uint8_t> keep_b(nb, 0u);
    std::vector<uint8_t> keep_h(nh, 0u);

    std::unordered_set<uint32_t> segs;
    segs.reserve((nb + nh) * 8);

    SAVE saver(output_file);
    std::ostringstream oss;

    oss << "H\tVN:Z:1.0\tCreator:liftasm bubble\n";
    oss << "H\tCO:Z:Tags:\n";
    oss << "H\tCO:Z:  ID:i   record id\n";
    oss << "H\tCO:Z:  TP:Z   record type (bubble | homo)\n";
    oss << "H\tCO:Z:  SR:Z   source node(s)\n";
    oss << "H\tCO:Z:  SK:Z   sink node(s)\n";
    oss << "H\tCO:Z:  TY:Z   subtype\n";
    oss << "H\tCO:Z:  LN:i   total length (bp) of inner nodes\n";
    oss << "H\tCO:Z:  NC:i   number of unique inner nodes\n";
    saver.save(oss.str()); oss.str(""); oss.clear();

    std::size_t rid = 0;

    // Bubble H lines
    for (std::size_t i = 0; i < nb; ++i) {
        const Bubble& bb = bubbles_[i];

        const uint32_t src_v    = bb.get_source();
        const uint32_t sink_v   = bb.get_sink();
        const uint32_t src_seg  = Vertex::get_segment_id(src_v);
        const uint32_t sink_seg = Vertex::get_segment_id(sink_v);

        std::unordered_set<uint32_t> inner_segs;
        inner_segs.reserve(16);

        for (const auto& p : bb.get_paths()) {
            for (uint32_t v : p) {
                const uint32_t sid = Vertex::get_segment_id(v);
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

        keep_b[i] = 1u;

        const bool src_rev  = Vertex::get_is_reverse(src_v);
        const bool sink_rev = Vertex::get_is_reverse(sink_v);

        const std::string src_string  = graph_.getNodeName(src_seg)  + (src_rev  ? "-" : "+");
        const std::string sink_string = graph_.getNodeName(sink_seg) + (sink_rev ? "-" : "+");

        oss << "H\tID:i:" << rid++
            << "\tTP:Z:bubble"
            << "\tSR:Z:" << src_string
            << "\tSK:Z:" << sink_string
            << "\tTY:Z:" << bb.get_type()
            << "\tLN:i:" << inner_bp
            << "\tNC:i:" << inner_nodes
            << "\n";
        saver.save(oss.str()); oss.str(""); oss.clear();

        for (const auto& p : bb.get_paths()) {
            for (uint32_t v : p) {
                segs.insert(Vertex::get_segment_id(v));
            }
        }
    }

    // Homologous paths H lines
    for (std::size_t i = 0; i < nh; ++i) {
        const HomologousPath& hp = homologous_paths_[i];

        if (hp.paths.empty() || hp.sources.empty()) continue;

        std::unordered_set<uint32_t> inner_segs;
        inner_segs.reserve(32);

        for (const auto& path : hp.paths) {
            for (const Vertex& v : path) {
                inner_segs.insert(v.segment_id());
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

        keep_h[i] = 1u;

        std::string src_string, sink_string;
        src_string.reserve(hp.sources.size() * 24);
        sink_string.reserve(hp.ends.size() * 24);

        uint32_t source_num;
        if (hp.type == Type::SameSource) {
            source_num = 1;
        } else {
            source_num = static_cast<uint32_t>(hp.sources.size());
        }
        for (std::size_t j = 0; j < source_num; ++j) {
            if (j) src_string.push_back(',');
            src_string += graph_.getNodeName(hp.sources[j].segment_id());
            src_string += (hp.sources[j].is_reverse() ? "-" : "+");
        }
        for (std::size_t j = 0; j < hp.ends.size(); ++j) {
            if (j) sink_string.push_back(',');
            sink_string += graph_.getNodeName(hp.ends[j].segment_id());
            sink_string += (hp.ends[j].is_reverse() ? "-" : "+");
        }

        oss << "H\tID:i:" << rid++
            << "\tTP:Z:homo"
            << "\tSR:Z:" << src_string
            << "\tSK:Z:" << sink_string
            << "\tTY:Z:" << hp.get_type()
            << "\tLN:i:" << inner_bp
            << "\tNC:i:" << inner_nodes
            << "\n";
        saver.save(oss.str()); oss.str(""); oss.clear();

        for (const auto& path : hp.paths) {
            for (const Vertex& v : path) {
                segs.insert(v.segment_id());
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
        if (node->sequence.empty() || !write_seq) oss << "*";
        else oss << node->sequence;
        oss << "\tLN:i:" << graph_.getNodeLength(sid) << "\n";

        saver.save(oss.str()); oss.str(""); oss.clear();
    }

    // L lines
    for (const auto& arc : graph_.getAllArcs()) {
        if (arc.get_del() || arc.get_comp()) continue;

        const uint32_t v = arc.get_source_vertex_id();
        const uint32_t w = arc.get_target_vertex_id();

        const uint32_t sv = Vertex::get_segment_id(v);
        const uint32_t sw = Vertex::get_segment_id(w);
        if (!segs.count(sv) || !segs.count(sw)) continue;

        oss << "L\t" << graph_.getNodeName(sv) << "\t" << ((v & 1) ? '-' : '+') << "\t"
            << graph_.getNodeName(sw) << "\t" << ((w & 1) ? '-' : '+') << "\t"
            << graph_.format_overlap_field(arc.ov, arc.ow) << "\n";

        saver.save(oss.str()); oss.str(""); oss.clear();
    }

    // P lines
    auto path_name_string_u32 = [&](const std::vector<uint32_t>& path) -> std::string {
        if (path.empty()) return "*";
        std::string s;
        s.reserve(path.size() * 16);
        for (size_t i = 0; i < path.size(); ++i) {
            const uint32_t sid = Vertex::get_segment_id(path[i]);
            const bool rev = Vertex::get_is_reverse(path[i]);
            s.push_back(rev ? '<' : '>');
            s += graph_.getNodeName(sid);
        }
        return s.empty() ? "*" : s;
    };

    auto path_name_string_vtx = [&](const std::vector<Vertex>& path) -> std::string {
        if (path.empty()) return "*";
        std::string s;
        s.reserve(path.size() * 16);
        for (size_t i = 0; i < path.size(); ++i) {
            s.push_back(path[i].is_reverse() ? '<' : '>');
            s += graph_.getNodeName(path[i].segment_id());
        }
        return s.empty() ? "*" : s;
    };

    auto overlap_vec_to_string = [](const std::vector<uint32_t>& ov) -> std::string {
        if (ov.empty()) return "*";

        std::string s;
        s.reserve(ov.size() * 6);

        for (size_t i = 0; i < ov.size(); ++i) {
            if (i) s.push_back(',');
            s += std::to_string(ov[i]);
            s.push_back('M');
        }

        return s;
    };

    auto path_span_u32 = [&](const std::vector<uint32_t>& path, const std::vector<uint32_t>& ov, uint32_t src_vtx, uint32_t sink_vtx) -> uint64_t {
        uint64_t total = 0;
        for (uint32_t v : path) {
            if (v == src_vtx || v == sink_vtx) continue;
            const uint32_t sid = Vertex::get_segment_id(v);
            total += graph_.getNodeLength(sid);
        }
        uint64_t ov_sum = 0;
        for (uint32_t x : ov) ov_sum += x;
        return (total >= ov_sum) ? (total - ov_sum) : 0;
    };

    auto path_span_vtx = [&](const std::vector<Vertex>& path, const std::vector<uint32_t>& ov, const std::unordered_set<uint32_t>& exclude) -> uint64_t {
        uint64_t total = 0;
        for (const auto& v : path) {
            const uint32_t sid = v.segment_id();
            if (exclude.count(sid)) continue;
            total += graph_.getNodeLength(sid);
        }
        uint64_t ov_sum = 0;
        for (uint32_t x : ov) ov_sum += x;
        return (total >= ov_sum) ? (total - ov_sum) : 0;
    };

    std::size_t rid_p = 0;

    // Bubble P lines
    for (std::size_t i = 0; i < nb; ++i) {
        if (!keep_b[i]) continue;

        const Bubble& bb = bubbles_[i];
        const std::size_t cur_id = rid_p++;

        const auto& paths = bb.get_paths();
        for (std::size_t pid = 0; pid < paths.size(); ++pid) {
            const auto& p = paths[pid];

            std::vector<uint32_t> vertex_list;
            vertex_list.reserve(p.size());

            for (uint32_t v : p) {
                const uint32_t sid = NodeHandle::get_segment_id(v);
                if (sid == NodeHandle::get_segment_id(bb.get_source()) || sid == NodeHandle::get_segment_id(bb.get_sink())) {
                    continue;
                }
                vertex_list.push_back(v);
            }

            if (vertex_list.empty()) continue;

            const std::string seg_list_str = path_name_string_u32(vertex_list);
            const auto overlaps = graph_.get_path_overlaps(vertex_list);
            const std::string ovlp_str = overlap_vec_to_string(overlaps);
            const std::string seq = graph_.get_path_sequence(vertex_list);

            uint64_t path_len = 0;
            if (!seq.empty() && seq!="*") {
                path_len = seq.size();
            } else {
                path_len = path_span_u32(vertex_list, overlaps, bb.get_source(), bb.get_sink());
            }

            oss << "P\t"
                << "b" << cur_id << ".p" << pid << "\t"
                << seg_list_str << "\t"
                << ovlp_str << "\t"
                << "ID:i:" << cur_id << "\t"
                << "TP:Z:bubble\t"
                << "PI:i:" << pid << "\t"
                << "LN:i:" << path_len << "\t"
                << "SEQ:Z:" << (seq.empty() || !write_seq ? "*" : seq)
                << "\n";

            saver.save(oss.str()); oss.str(""); oss.clear();
        }
    }

    // Homologous P lines
    for (std::size_t i = 0; i < nh; ++i) {
        if (!keep_h[i]) continue;

        const HomologousPath& hp = homologous_paths_[i];
        const std::size_t cur_id = rid_p++;

        for (std::size_t pid = 0; pid < hp.paths.size(); ++pid) {
            const auto& p = hp.paths[pid];
            if (p.empty()) continue;

            const std::string seg_list_str = path_name_string_vtx(p);
            const auto overlaps = graph_.get_path_overlaps(p);
            const std::string ovlp_str = overlap_vec_to_string(overlaps);
            const std::string seq = graph_.get_path_sequence(p);

            std::unordered_set<uint32_t> exclude;
            uint64_t path_len = 0;
            if (!seq.empty() && seq!="*") {
                path_len = seq.size();
            } else {
                path_len = path_span_vtx(p, overlaps, exclude);
            }

            oss << "P\t"
                << "h" << cur_id << ".p" << pid << "\t"
                << seg_list_str << "\t"
                << ovlp_str << "\t"
                << "ID:i:" << cur_id << "\t"
                << "TP:Z:homo\t"
                << "PI:i:" << pid << "\t"
                << "LN:i:" << path_len << "\t"
                << "SEQ:Z:" << (seq.empty() || !write_seq ? "*" : seq)
                << "\n";

            saver.save(oss.str()); oss.str(""); oss.clear();
        }
    }

    log_stream() << "  - Total records written: " << nb + nh << "\n\n";
}

void GfaBubbleFinder::save_bubble_as_vcf(
    const std::string& output_prefix,
    const std::string& paf_file,
    const std::string& ref_file,
    uint32_t min_mapq,
    uint32_t min_aln_len
) const
{
    const bool use_liftover = !paf_file.empty();
    const bool use_ref = !ref_file.empty();
    const bool use_ref_sequence = use_liftover && use_ref;
    static constexpr uint32_t LIFT_FLANK_BP = 100;

    log_stream() << "Writing bubble VCF to " << output_prefix << ".bubbles.vcf ...\n";
    if (!use_liftover) {
        log_stream() << "  - No PAF provided; using bubble positions from GFA\n";
    }
    if (!use_ref) {
        log_stream() << "  - No reference provided; REF will not be extracted\n";
    }
    if (!use_liftover && use_ref) {
        log_stream() << "  - No PAF provided; REF will not be extracted\n";
    }

    std::unordered_map<std::string, std::vector<liftover::Pafinfo>> pafs;
    std::unordered_map<std::string, std::vector<uint32_t>> paf_ends;
    std::unordered_map<std::string, std::vector<uint32_t>> paf_idx_ends;
    if (use_liftover) {
        paf_reader(paf_file, min_mapq, min_aln_len, pafs, paf_ends, paf_idx_ends);
    }

    liftover::SeqDB db;
    if (use_ref_sequence && !db.load(ref_file)) {
        error_stream() << "Load ref failed: " << ref_file << "\n";
        std::exit(1);
    }

    SAVE saver(output_prefix + ".bubbles.vcf");
    std::ostringstream oss;

    auto join_strings = [](const std::vector<std::string>& xs, const std::string& sep) -> std::string {
        if (xs.empty()) return "";
        std::string out = xs[0];
        for (size_t i = 1; i < xs.size(); ++i) {
            out += sep;
            out += xs[i];
        }
        return out;
    };

    auto path_name_string_u32 = [&](const std::vector<uint32_t>& path, uint32_t src_seg, uint32_t sink_seg) -> std::string {
        std::string s;
        s.reserve(path.size() * 24);
        for (uint32_t v : path) {
            const uint32_t sid = Vertex::get_segment_id(v);
            if (sid == src_seg || sid == sink_seg) continue;
            s.push_back(Vertex::get_is_reverse(v) ? '<' : '>');
            s += graph_.getNodeName(sid);
        }
        return s.empty() ? "*" : s;
    };

    struct BoundaryWindow {
        std::string chrom;
        uint32_t win_beg{0};      // [win_beg, win_end)
        uint32_t win_end{0};
        uint32_t boundary{0};     // boundary base on query/root
        bool boundary_is_left{true};
    };

    auto get_boundary_window_from_name = [&](const std::string& name_with_sign, bool is_source, uint64_t plain_len, BoundaryWindow& bw) -> bool
    {
        bw = BoundaryWindow{};
        if (name_with_sign.empty()) return false;

        bool is_rev = false;
        std::string base = name_with_sign;
        if (!base.empty() && (base.back() == '+' || base.back() == '-')) {
            is_rev = (base.back() == '-');
            base.pop_back();
        }

        if (base.empty()) return false;

        // utg000175l
        if (base.find(':') == std::string::npos) {
            if (plain_len == 0) return false;
            bw.chrom = base;
            bw.win_beg = 0;
            bw.win_end = static_cast<uint32_t>(std::min<uint64_t>(plain_len, uint64_t(UINT32_MAX)));
            if (bw.win_end == 0) return false;

            if (is_source) {
                bw.boundary = is_rev ? 0u : (bw.win_end - 1);
            } else {
                bw.boundary = is_rev ? (bw.win_end - 1) : 0u;
            }
            bw.boundary_is_left = (bw.boundary == bw.win_beg);
            return true;
        }

        // utg000001l:3377-3469
        gfaName gname(';');
        std::string norm = gname.force_name_dir(base, is_rev);
        auto pieces = gname.parse_composite_with_dir(norm);
        if (pieces.empty()) return false;

        const auto& p = is_source ? pieces.back() : pieces.front();
        if (p.hi <= p.lo) return false;

        bw.chrom = p.root;

        // In path direction:
        // - source boundary is at the END of the source piece
        // - sink boundary   is at the START of the sink piece
        const bool take_head = is_source ? p.rev : !p.rev;

        if (take_head) {
            bw.win_beg = static_cast<uint32_t>(p.lo);
            bw.win_end = static_cast<uint32_t>(std::min<uint64_t>(p.hi, p.lo + LIFT_FLANK_BP));
            bw.boundary = bw.win_beg;
            bw.boundary_is_left = true;
        } else {
            const uint64_t wb = (p.hi > LIFT_FLANK_BP) ? std::max<uint64_t>(p.lo, p.hi - LIFT_FLANK_BP) : p.lo;
            bw.win_beg = static_cast<uint32_t>(wb);
            bw.win_end = static_cast<uint32_t>(p.hi);
            if (bw.win_end <= bw.win_beg) return false;
            bw.boundary = bw.win_end - 1;
            bw.boundary_is_left = false;
        }

        return bw.win_end > bw.win_beg;
    };

    auto map_boundary_from_hit = [](const liftover::LIFTresult& hit, const BoundaryWindow& bw, uint32_t& ref_pos) -> bool
    {
        if (hit.tend <= hit.tbeg) return false;
        if (bw.boundary_is_left) {
            ref_pos = (hit.strand == '+') ? static_cast<uint32_t>(hit.tbeg) : static_cast<uint32_t>(hit.tend - 1);
        } else {
            ref_pos = (hit.strand == '+') ? static_cast<uint32_t>(hit.tend - 1) : static_cast<uint32_t>(hit.tbeg);
        }
        return true;
    };

    auto first_base_all_same = [](const std::string& ref, const std::vector<std::string>& alleles) -> bool {
        if (ref.empty()) return false;
        const char c = ref.front();
        for (const auto& a : alleles) {
            if (a.empty()) return false;
            if (a.front() != c) return false;
        }
        return true;
    };

    auto can_trim_common_prefix = [](const std::string& ref, const std::vector<std::string>& alleles) -> bool {
        if (ref.size() <= 1 || alleles.empty()) return false;
        const char c = ref.front();
        for (const auto& a : alleles) {
            if (a.size() <= 1 || a.front() != c) return false;
        }
        return true;
    };

    auto can_trim_common_suffix = [](const std::string& ref, const std::vector<std::string>& alleles) -> bool {
        if (ref.size() <= 1 || alleles.empty()) return false;
        const char c = ref.back();
        for (const auto& a : alleles) {
            if (a.size() <= 1 || a.back() != c) return false;
        }
        return true;
    };

    auto pn_root_or_empty = [&](const std::string& pn_token) -> std::string {
        if (pn_token.empty() || pn_token == "*") return "";
        std::string name = pn_token;

        if (!name.empty() && (name.front() == '<' || name.front() == '>')) {
            name.erase(name.begin());
        }
        size_t cut = name.find_first_of("<>,");
        if (cut != std::string::npos) {
            name = name.substr(0, cut);
        }

        BoundaryWindow bw;
        if (get_boundary_window_from_name(name, true, 1, bw)) {
            return bw.chrom;
        }
        return name;
    };

    auto record_to_string = [&](const VcfRecord& r) -> std::string {
        std::ostringstream out;
        const std::string pn_field = join_strings(r.path_names, ",");
        const std::string gt_field = join_strings(r.gt_tokens, std::string(1, r.gt_sep));

        out << r.chrom << "\t"
            << r.pos << "\t"
            << r.id << "\t"
            << r.ref << "\t"
            << r.alt << "\t"
            << ".\tPASS\t"
            << "END=" << r.end
            << ";SRC=" << r.src
            << ";SNK=" << r.snk
            << ";NS=" << r.ns
            << ";LN=" << r.ln
            << ";NC=" << r.nc
            << ";PN=" << pn_field
            << "\tGT\t"
            << gt_field
            << "\n";

        return out.str();
    };

    auto normalize_pn_gt_order = [&](std::vector<VcfRecord>& records, size_t idx) {
        VcfRecord& r = records[idx];
        r.id = std::to_string(idx);
        const size_t n = r.path_names.size();

        r.gt_sep = '/';
        if (n <= 1 || r.gt_tokens.size() != n) return;

        std::vector<std::string> roots(n);
        std::vector<size_t> star_idx;
        star_idx.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            roots[i] = pn_root_or_empty(r.path_names[i]);
            if (roots[i].empty()) star_idx.push_back(i);
        }

        std::vector<std::string> prev_roots;
        if (idx > 0) {
            const VcfRecord& prev = records[idx - 1];
            for (const auto& pn : prev.path_names) {
                std::string root = pn_root_or_empty(pn);
                if (root.empty()) continue;
                if (std::find(prev_roots.begin(), prev_roots.end(), root) == prev_roots.end()) {
                    prev_roots.push_back(std::move(root));
                }
            }
        }

        if (!prev_roots.empty() && !star_idx.empty()) {
            std::unordered_set<std::string> used;
            for (const auto& root : roots) {
                if (!root.empty()) used.insert(root);
            }

            std::vector<std::string> missing;
            for (const auto& root : prev_roots) {
                if (!used.count(root)) missing.push_back(root);
            }

            const size_t m = std::min(star_idx.size(), missing.size());
            for (size_t k = 0; k < m; ++k) {
                roots[star_idx[k]] = missing[k];
            }
        }

        for (size_t i : star_idx) {
            if (roots[i].empty()) roots[i] = "~";
        }

        std::vector<size_t> order(n);
        std::iota(order.begin(), order.end(), 0);

        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            if (roots[a] != roots[b]) return roots[a] < roots[b];
            return a < b;
        });

        bool same_order = true;
        for (size_t i = 0; i < n; ++i) {
            if (order[i] != i) {
                same_order = false;
                break;
            }
        }
        if (same_order) return;

        std::vector<std::string> new_names;
        std::vector<std::string> new_gts;
        new_names.reserve(n);
        new_gts.reserve(n);

        for (size_t i : order) {
            new_names.push_back(std::move(r.path_names[i]));
            new_gts.push_back(std::move(r.gt_tokens[i]));
        }

        r.path_names = std::move(new_names);
        r.gt_tokens  = std::move(new_gts);
    };

    oss << "##fileformat=VCFv4.2\n";
    oss << "##source=liftasm bubble\n";
    for (const auto& name : db.names) {
        const std::string* seq = db.get(name);
        if (!seq) continue;
        oss << "##contig=<ID=" << name << ",length=" << seq->size() << ">\n";
    }
    oss << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n";
    oss << "##INFO=<ID=SRC,Number=1,Type=String,Description=\"Source node\">\n";
    oss << "##INFO=<ID=SNK,Number=1,Type=String,Description=\"Sink node\">\n";
    oss << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of paths\">\n";
    oss << "##INFO=<ID=LN,Number=1,Type=Integer,Description=\"Total inner bp\">\n";
    oss << "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of inner nodes\">\n";
    oss << "##INFO=<ID=PN,Number=.,Type=String,Description=\"Path names in genotype order\">\n";
    oss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    oss << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << output_prefix << "\n";
    saver.save(oss.str());
    oss.str(""); oss.clear();

    std::vector<VcfRecord> vcf_records;
    vcf_records.reserve(bubbles_.size());

    for (const Bubble& bb : bubbles_) {
        if (bb.get_type() != "Normal") continue;

        const uint32_t src_v    = bb.get_source();
        const uint32_t sink_v   = bb.get_sink();
        const uint32_t src_seg  = Vertex::get_segment_id(src_v);
        const uint32_t sink_seg = Vertex::get_segment_id(sink_v);

        const std::string src_name = graph_.getNodeName(src_seg) + (Vertex::get_is_reverse(src_v) ? "-" : "+");
        const std::string sink_name = graph_.getNodeName(sink_seg) + (Vertex::get_is_reverse(sink_v) ? "-" : "+");

        const uint64_t src_len = graph_.getNodeLength(src_seg);
        const uint64_t sink_len = graph_.getNodeLength(sink_seg);
        if (src_len == 0 || sink_len == 0) continue;

        std::unordered_set<uint32_t> inner_segs;
        inner_segs.reserve(16);
        for (const auto& p : bb.get_paths()) {
            for (uint32_t v : p) {
                const uint32_t sid = Vertex::get_segment_id(v);
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

        BoundaryWindow src_bw, sink_bw;
        if (!get_boundary_window_from_name(src_name, true, src_len, src_bw)) continue;
        if (!get_boundary_window_from_name(sink_name, false, sink_len, sink_bw)) continue;

        std::string chrom;
        uint32_t ref_src_pos = 0;
        uint32_t ref_sink_pos = 0;
        bool bubble_need_rc_to_ref = false;

        if (!use_liftover) {
            chrom = src_bw.chrom;
            ref_src_pos = src_bw.boundary;
            ref_sink_pos = sink_bw.boundary;
            bubble_need_rc_to_ref = (ref_src_pos > ref_sink_pos);
        } else {
            std::vector<liftover::LIFTresult> src_hits = paf_liftover(
                pafs, paf_ends, paf_idx_ends, src_bw.chrom, src_bw.win_beg, src_bw.win_end
            );
            std::vector<liftover::LIFTresult> sink_hits = paf_liftover(
                pafs, paf_ends, paf_idx_ends, sink_bw.chrom, sink_bw.win_beg, sink_bw.win_end
            );

            if (src_hits.empty() || sink_hits.empty()) continue;

            bool found_pair = false;
            uint64_t best_span = std::numeric_limits<uint64_t>::max();

            for (const auto& sh : src_hits) {
                uint32_t sref = 0;
                if (!map_boundary_from_hit(sh, src_bw, sref)) continue;

                for (const auto& th : sink_hits) {
                    if (sh.tname != th.tname) continue;
                    if (sh.strand != th.strand) continue;

                    uint32_t tref = 0;
                    if (!map_boundary_from_hit(th, sink_bw, tref)) continue;

                    const uint64_t span = (sref > tref) ? (uint64_t(sref) - uint64_t(tref)) : (uint64_t(tref) - uint64_t(sref));
                    if (!found_pair || span < best_span) {
                        found_pair = true;
                        best_span = span;
                        chrom = sh.tname;
                        ref_src_pos = sref;
                        ref_sink_pos = tref;
                    }
                }
            }

            if (!found_pair) continue;
            bubble_need_rc_to_ref = (ref_src_pos > ref_sink_pos);
        }

        uint32_t beg = std::min(ref_src_pos, ref_sink_pos);
        uint32_t end = std::max(ref_src_pos, ref_sink_pos) + 1;
        if (end <= beg) end = beg + 1;

        if (DEBUG_ENABLED) {
            debug_stream() << "\n";
            debug_stream() << "src_name:" << src_name << "; sink_name:" << sink_name
                           << "; " << chrom << ":" << beg << "-" << end
                           << "; bubble_need_rc_to_ref:" << (bubble_need_rc_to_ref ? "true" : "false") << "\n";
        }

        std::string src_seq = graph_.get_path_sequence(std::vector<uint32_t>{src_v});
        std::string sink_seq = graph_.get_path_sequence(std::vector<uint32_t>{sink_v});
        const bool have_src_seq = !(src_seq.empty() || src_seq == "*");
        const bool have_sink_seq = !(sink_seq.empty() || sink_seq == "*");

        std::vector<std::string> path_names;
        std::vector<std::string> allele_seqs;   // sequence or ""
        std::vector<uint8_t> has_allele_seq;    // 1 = real sequence, 0 = missing

        path_names.reserve(bb.get_paths().size());
        allele_seqs.reserve(bb.get_paths().size());
        has_allele_seq.reserve(bb.get_paths().size());

        for (const auto& p : bb.get_paths()) {
            if (p.empty()) continue;

            path_names.push_back(path_name_string_u32(p, src_seg, sink_seg));

            if (!use_ref_sequence) {
                std::vector<uint32_t> inner_path;
                inner_path.reserve(p.size());
                for (uint32_t v : p) {
                    const uint32_t sid = Vertex::get_segment_id(v);
                    if (sid == src_seg || sid == sink_seg) continue;
                    inner_path.push_back(v);
                }

                if (inner_path.empty()) {
                    allele_seqs.emplace_back("");
                    has_allele_seq.push_back(1);
                } else {
                    std::string alt_seq = graph_.get_path_sequence(inner_path);
                    if (alt_seq.empty() || alt_seq == "*") {
                        allele_seqs.emplace_back("");
                        has_allele_seq.push_back(0);
                    } else {
                        if (bubble_need_rc_to_ref) alt_seq = seqUtils::revcomp(alt_seq);
                        allele_seqs.push_back(std::move(alt_seq));
                        has_allele_seq.push_back(1);
                    }
                }
            } else {
                if (!have_src_seq || !have_sink_seq) {
                    allele_seqs.emplace_back("");
                    has_allele_seq.push_back(0);
                    continue;
                }

                std::string full = graph_.get_path_sequence(p);
                if (full.empty() || full == "*") {
                    allele_seqs.emplace_back("");
                    has_allele_seq.push_back(0);
                    continue;
                }

                if (bubble_need_rc_to_ref) {
                    full = seqUtils::revcomp(full);
                }

                const size_t left_anchor = bubble_need_rc_to_ref ? (sink_seq.size() - 1) : (src_seq.size() - 1);
                const size_t right_anchor = bubble_need_rc_to_ref ? (full.size() - src_seq.size()) : (full.size() - sink_seq.size());

                if (left_anchor >= full.size() || right_anchor >= full.size() || right_anchor < left_anchor) {
                    allele_seqs.emplace_back("");
                    has_allele_seq.push_back(0);
                    continue;
                }

                allele_seqs.push_back(full.substr(left_anchor, right_anchor - left_anchor + 1));
                has_allele_seq.push_back(1);
            }
        }

        if (path_names.empty()) continue;

        const bool have_any_allele_seq = std::find(has_allele_seq.begin(), has_allele_seq.end(), uint8_t(1)) != has_allele_seq.end();

        if (!use_ref_sequence) {
            VcfRecord rec;
            rec.chrom = chrom;
            rec.pos = beg + 1;
            rec.end = end;
            rec.src = src_name;
            rec.snk = sink_name;
            rec.ns = static_cast<uint32_t>(bb.get_paths().size());
            rec.ln = inner_bp;
            rec.nc = inner_nodes;
            rec.path_names = std::move(path_names);

            if (!have_any_allele_seq) {
                rec.ref = ".";
                rec.alt = ".";
                rec.gt_tokens.assign(rec.path_names.size(), ".");
                vcf_records.push_back(std::move(rec));
                continue;
            }

            rec.ref = "N";

            std::vector<std::string> alt_alleles;
            std::unordered_map<std::string, int> alt_index;
            rec.gt_tokens.reserve(allele_seqs.size());

            for (size_t i = 0; i < allele_seqs.size(); ++i) {
                if (!has_allele_seq[i]) {
                    rec.gt_tokens.push_back(".");
                    continue;
                }

                const std::string& seq = allele_seqs[i];
                auto it = alt_index.find(seq);
                if (it == alt_index.end()) {
                    const int idx = static_cast<int>(alt_alleles.size()) + 1;
                    alt_alleles.push_back(seq.empty() ? "." : seq);
                    alt_index.emplace(seq, idx);
                    rec.gt_tokens.push_back(std::to_string(idx));
                } else {
                    rec.gt_tokens.push_back(std::to_string(it->second));
                }
            }

            if (alt_alleles.empty()) {
                rec.alt = ".";
                if (rec.gt_tokens.empty()) rec.gt_tokens.assign(rec.path_names.size(), ".");
            } else {
                rec.alt = join_strings(alt_alleles, ",");
            }

            vcf_records.push_back(std::move(rec));
            continue;
        }

        const std::string* ref_ptr = db.get(chrom);
        if (!ref_ptr) continue;
        if (end > ref_ptr->size()) continue;

        VcfRecord rec;
        rec.chrom = chrom;
        rec.pos = beg + 1;
        rec.end = end;
        rec.src = src_name;
        rec.snk = sink_name;
        rec.ns = static_cast<uint32_t>(bb.get_paths().size());
        rec.ln = inner_bp;
        rec.nc = inner_nodes;
        rec.path_names = std::move(path_names);

        if (!have_any_allele_seq) {
            rec.ref = ".";
            rec.alt = ".";
            rec.gt_tokens.assign(rec.path_names.size(), ".");
            vcf_records.push_back(std::move(rec));
            continue;
        }

        std::string ref_seq = ref_ptr->substr(beg, end - beg);
        if (ref_seq.empty()) continue;

        uint32_t out_pos1 = beg + 1;

        std::vector<std::string> called_alleles;
        std::vector<size_t> called_idx;
        called_alleles.reserve(allele_seqs.size());
        called_idx.reserve(allele_seqs.size());

        for (size_t i = 0; i < allele_seqs.size(); ++i) {
            if (!has_allele_seq[i]) continue;
            called_idx.push_back(i);
            called_alleles.push_back(allele_seqs[i]);
        }

        while (!called_alleles.empty() && !first_base_all_same(ref_seq, called_alleles)) {
            if (beg == 0) break;
            --beg;
            out_pos1 = beg + 1;
            ref_seq.insert(ref_seq.begin(), (*ref_ptr)[beg]);
            for (auto& a : called_alleles) {
                a.insert(a.begin(), (*ref_ptr)[beg]);
            }
        }

        while (can_trim_common_prefix(ref_seq, called_alleles)) {
            ref_seq.erase(ref_seq.begin());
            for (auto& a : called_alleles) a.erase(a.begin());
            ++out_pos1;
        }

        while (can_trim_common_suffix(ref_seq, called_alleles)) {
            ref_seq.pop_back();
            for (auto& a : called_alleles) a.pop_back();
        }

        std::vector<std::string> alt_alleles;
        std::unordered_map<std::string, int> alt_index;
        rec.gt_tokens.reserve(allele_seqs.size());

        size_t ci = 0;
        for (size_t i = 0; i < allele_seqs.size(); ++i) {
            if (!has_allele_seq[i]) {
                rec.gt_tokens.push_back(".");
                continue;
            }

            const std::string& seq = called_alleles[ci++];
            if (seq == ref_seq) {
                rec.gt_tokens.push_back("0");
            } else {
                auto it = alt_index.find(seq);
                if (it == alt_index.end()) {
                    const int idx = static_cast<int>(alt_alleles.size()) + 1;
                    alt_alleles.push_back(seq);
                    alt_index.emplace(seq, idx);
                    rec.gt_tokens.push_back(std::to_string(idx));
                } else {
                    rec.gt_tokens.push_back(std::to_string(it->second));
                }
            }
        }

        if (alt_alleles.empty()) continue;

        rec.pos = out_pos1;
        rec.end = out_pos1 + static_cast<uint32_t>(ref_seq.size()) - 1;
        rec.ref = std::move(ref_seq);
        rec.alt = join_strings(alt_alleles, ",");
        vcf_records.push_back(std::move(rec));
    }

    std::sort(vcf_records.begin(), vcf_records.end(),
        [](const VcfRecord& a, const VcfRecord& b) {
            if (a.chrom != b.chrom) return a.chrom < b.chrom;
            return a.pos < b.pos;
        }
    );

    for (size_t i = 0; i < vcf_records.size(); ++i) {
        normalize_pn_gt_order(vcf_records, i);
        saver.save(record_to_string(vcf_records[i]));
    }

    log_stream() << "  - Total VCF records written: " << vcf_records.size() << "\n\n";
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
        log_stream() << "  - Total Nodes Length: " << total_bp << "\n";
        log_stream() << "  - Total Nodes Count: " << inner_segs.size() << "\n";
        log_stream() << "  - type: " << bb.get_type() << "\n";

        const auto& ps = bb.get_paths();
        std::string path_str;
        for (size_t i = 0; i < ps.size(); ++i) {
            path_str = "  - path[" + std::to_string(i) + "]: ";
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

    // used to avoid duplicate forks
    std::unordered_set<uint64_t> seen; 
    seen.reserve(1024);

    ProgressTracker prog(V);

    for (uint32_t v = 0; v < V; ++v) {
        prog.hit();
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
    }
    log_stream() << "  - Total forks detected: " << forks_.size() << "\n" << "\n";
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
            log_stream() << "  - leaf[" << i << "]: " << graph_.getNodeName(seg) << (rev ? "-" : "+") << "\n";
        }
    }
}


std::vector<Vertex> GfaBubbleFinder::collect_sources_for_homo_()
{
    log_stream() << "Collecting sources for homologous path enumeration ...\n";

    std::vector<Vertex> sources;

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (V == 0) return sources;

    sources.reserve(V / 16 + 8);

    for (uint32_t v = 0; v < V; ++v) {
        if (graph_.getNodeDeleted(NodeHandle::get_segment_id(v))) continue;

        const uint32_t len = graph_.getNodeLength(NodeHandle::get_segment_id(v));

        // Single-node (Chromosome-like)
        if (is_single_node_(Vertex(v))) {
            sources.emplace_back(Vertex(v));
            continue;
        }

        if (len < diff_min_src_) continue;

        sources.emplace_back(Vertex(v));
    }

    // Print source names
    if (DEBUG_ENABLED) {
        for (const Vertex& src : sources) {
            uint32_t seg = src.segment_id();
            bool     rev = src.is_reverse();
            debug_stream() << "  - source: " << graph_.getNodeName(seg) << (rev ? "-" : "+") << "\n";
        }
    }

    log_stream() << "  - sources: " << sources.size() << "\n" << "\n";

    return sources;
}

std::vector<minimizerdna::Sketch> GfaBubbleFinder::build_all_node_mm_library_() const
{
    log_stream() << "Building minimizer library for all nodes ...\n";

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    std::vector<minimizerdna::Sketch> node_sketches(V);

    const unsigned threads = std::max(1u, thread_);
    ThreadPool pool(threads);

    std::vector<std::future<std::pair<uint32_t, minimizerdna::Sketch>>> futs;
    futs.reserve(V);

    ProgressTracker prog(V);

    for (uint32_t vid = 0; vid < V; ++vid) {
        prog.hit();

        futs.emplace_back(
            pool.submit([&, vid]() -> std::pair<uint32_t, minimizerdna::Sketch> {
                const Vertex v(vid);

                if (graph_.getNodeDeleted(v.segment_id())) {
                    return {vid, minimizerdna::Sketch{}};
                }

                const std::string seq = graph_.get_oriented_sequence(v);
                if (seq.empty() || seq == "*") {
                    return {vid, minimizerdna::Sketch{}};
                }

                const minimizerdna::MinimizerBuilder mzb(mm_opt_);
                return {vid, mzb.build(seq)};
            })
        );
    }

    size_t valid = 0;

    for (auto& f : futs) {
        auto [vid, sk] = f.get();

        if (vid < node_sketches.size()) {
            if (!sk.empty()) ++valid;
            node_sketches[vid] = std::move(sk);
        }
    }

    pool.stop();

    log_stream() << "  - Total valid node sketches: " << valid << "\n\n";

    return node_sketches;
}


std::vector<NodeGroupEntry> GfaBubbleFinder::cluster_sources_by_mm_(
    const std::vector<Vertex>& sources,
    const std::vector<minimizerdna::Sketch>& node_sketches,
    uint32_t min_len,
    double sim_thr
) const
{
    log_stream() << "Clustering sources by minimizer similarity ...\n";

    const minimizerdna::MinimizerBuilder mzb(mm_opt_);

    struct SourceRef {
        Vertex vtx;
        uint32_t len;
    };

    std::vector<SourceRef> refs;
    refs.reserve(sources.size());

    for (Vertex src : sources) {
        const uint32_t vid = src.vertex_id();

        if (vid >= node_sketches.size()) continue;
        if (node_sketches[vid].empty()) continue;
        if (graph_.getNodeDeleted(src.segment_id())) continue;

        const uint32_t len = graph_.getNodeLength(src.segment_id());
        if (len < min_len) continue;

        refs.push_back(SourceRef{src, len});
    }

    std::vector<NodeGroupEntry> groups;
    const size_t n = refs.size();

    const unsigned threads = std::max(1u, thread_);
    ThreadPool pool(threads);

    std::vector<uint8_t> used(n, 0);

    ProgressTracker prog(n);
    std::size_t skipped_connected = 0;

    if (n < 2) {
        pool.stop();
        log_stream() << "  - Total source groups: 0\n";
        log_stream() << "  - Node pairs skipped due to connectivity in the graph: 0\n\n";
        return groups;
    }

    for (size_t i = 0; i < n; ++i) {
        prog.hit();

        if (used[i]) continue;

        const SourceRef& a = refs[i];
        const uint32_t comp_a = graph_.node_connectivity_id(a.vtx);
        const minimizerdna::Sketch& sketch_a = node_sketches[a.vtx.vertex_id()];

        NodeGroupEntry g;
        g.vtxs.reserve(2);
        g.vtxs.push_back(a.vtx);
        used[i] = 1;

        std::vector<std::future<std::pair<size_t, double>>> futs;
        futs.reserve((n > i + 1) ? (n - i - 1) : 0);

        for (size_t j = i + 1; j < n; ++j) {
            if (used[j]) continue;

            const SourceRef& b = refs[j];
            const uint32_t comp_b = graph_.node_connectivity_id(b.vtx);

            if (comp_a != UINT32_MAX && comp_b != UINT32_MAX && comp_a == comp_b) {
                ++skipped_connected;
                continue;
            }

            futs.emplace_back(
                pool.submit([&, j]() -> std::pair<size_t, double> {
                    const uint32_t vid_b = refs[j].vtx.vertex_id();
                    return {j, mzb.max_containment(sketch_a, node_sketches[vid_b])};
                })
            );
        }

        struct BestHit {
            size_t idx{SIZE_MAX};
            double sim{0.0};
            uint64_t len{0};
        };

        std::unordered_map<uint32_t, BestHit> best_by_comp;
        best_by_comp.reserve(futs.size());

        for (auto& f : futs) {
            auto [j, sim] = f.get();

            if (j == SIZE_MAX) continue;
            if (used[j]) continue;
            if (sim < sim_thr) continue;

            const SourceRef& b = refs[j];
            const uint32_t comp_b = graph_.node_connectivity_id(b.vtx);

            if (comp_b == UINT32_MAX) continue;

            BestHit hit;
            hit.idx = j;
            hit.sim = sim;
            hit.len = b.len;

            auto it = best_by_comp.find(comp_b);

            if (it == best_by_comp.end()) {
                best_by_comp.emplace(comp_b, hit);
            } else {
                const BestHit& old = it->second;

                if (hit.sim > old.sim ||
                    (hit.sim == old.sim && hit.len > old.len) ||
                    (hit.sim == old.sim && hit.len == old.len && hit.idx < old.idx)) {
                    it->second = hit;
                }
            }
        }

        for (const auto& kv : best_by_comp) {
            const BestHit& hit = kv.second;
            const size_t j = hit.idx;

            if (j == SIZE_MAX) continue;
            if (used[j]) continue;

            if (g.min_similarity <= 0.0) {
                g.min_similarity = hit.sim;
            } else {
                g.min_similarity = std::min(g.min_similarity, hit.sim);
            }

            g.vtxs.push_back(refs[j].vtx);
            used[j] = 1;
        }

        if (g.vtxs.size() >= 2) {
            groups.emplace_back(std::move(g));
        }
    }

    pool.stop();

    if (DEBUG_ENABLED) {
        for (const auto& g : groups) {
            std::string group_str;

            for (const Vertex& v : g.vtxs) {
                group_str += graph_.getNodeName(v.segment_id());
                group_str += (v.is_reverse() ? "-" : "+");
                group_str += ' ';
            }

            debug_stream() << "  - group: " << group_str << "  min_sim: " << g.min_similarity << "\n";
        }

        debug_stream() << "\n\n";
    }

    log_stream() << "  - Total source groups: " << groups.size() << "\n";
    log_stream() << "  - Node pairs skipped due to connectivity in the graph: " << skipped_connected << "\n\n";

    return groups;
}



bool GfaBubbleFinder::homo_vertices_already_used_(
    const std::vector<Vertex>& vertices,
    const HomoUsedIntervals_& used_intervals
) const
{
    struct Item {
        uint32_t vid;
        uint32_t comp;
        uint32_t rank;
    };

    std::vector<Item> xs;
    xs.reserve(vertices.size());

    for (Vertex v : vertices) {
        const uint32_t vid = v.vertex_id();
        const uint32_t comp = graph_.node_connectivity_id(v);
        const uint32_t rank = graph_.vertex_topo_rank(v);

        if (vid == UINT32_MAX || comp == UINT32_MAX || rank == UINT32_MAX) continue;
        xs.push_back({vid, comp, rank});
    }

    std::sort(xs.begin(), xs.end(),
        [](const Item& a, const Item& b) {
            return a.vid < b.vid;
        }
    );

    xs.erase(
        std::unique(xs.begin(), xs.end(),
            [](const Item& a, const Item& b) {
                return a.vid == b.vid;
            }
        ),
        xs.end()
    );

    if (xs.size() < 2) return false;

    auto covered = [&](const HomoUsedKey_& key, uint32_t rank) -> bool {
        const auto it = used_intervals.find(key);
        if (it == used_intervals.end()) return false;

        for (const HomoUsedInterval_& x : it->second) {
            if (rank <= x.beg) return false;
            if (rank < x.end) return true;
        }

        return false;
    };

    for (size_t i = 0; i + 1 < xs.size(); ++i) {
        for (size_t j = i + 1; j < xs.size(); ++j) {
            uint32_t a = xs[i].comp;
            uint32_t b = xs[j].comp;
            if (a > b) std::swap(a, b);

            if (!covered(HomoUsedKey_{a, b, xs[i].comp}, xs[i].rank)) return false;
            if (!covered(HomoUsedKey_{a, b, xs[j].comp}, xs[j].rank)) return false;
        }
    }

    return true;
}

bool GfaBubbleFinder::homo_path_already_used_(
    const HomologousPath& hp,
    const HomoUsedIntervals_& used_intervals
) const
{
    std::vector<Vertex> anchors;

    if (hp.type == Type::SameSource) {
        anchors.insert(anchors.end(), hp.sources.begin(), hp.sources.end());
    } else {
        anchors.insert(anchors.end(), hp.starts.begin(), hp.starts.end());
    }

    return homo_vertices_already_used_(anchors, used_intervals);
}

void GfaBubbleFinder::insert_homo_used_interval_(
    HomoUsedIntervals_& used_intervals,
    const HomoUsedKey_& key,
    HomoUsedInterval_ x
) const
{
    if (key.pair_a == UINT32_MAX || key.pair_b == UINT32_MAX || key.comp == UINT32_MAX) return;
    if (x.beg == UINT32_MAX || x.end == UINT32_MAX) return;
    if (x.beg > x.end) std::swap(x.beg, x.end);

    auto& xs = used_intervals[key];
    xs.push_back(x);

    std::sort(xs.begin(), xs.end(),
        [](const HomoUsedInterval_& a, const HomoUsedInterval_& b) {
            if (a.beg != b.beg) return a.beg < b.beg;
            return a.end < b.end;
        }
    );

    size_t w = 0;

    for (const HomoUsedInterval_& y : xs) {
        if (w == 0 || y.beg > xs[w - 1].end + 1) {
            xs[w++] = y;
        } else {
            xs[w - 1].end = std::max(xs[w - 1].end, y.end);
        }
    }

    xs.resize(w);
}

void GfaBubbleFinder::mark_homo_used_intervals_(
    HomoUsedIntervals_& used_intervals,
    const HomologousPath& hp
) const
{
    struct Arm {
        uint32_t comp;
        uint32_t beg;
        uint32_t end;
    };

    const size_t n = std::min(hp.sources.size(), std::min(hp.starts.size(), hp.ends.size()));
    if (n < 2 || hp.paths.empty()) return;

    std::vector<Arm> arms;
    arms.reserve(n * 2);

    for (size_t i = 0; i < n; ++i) {
        const Vertex src = (hp.type == Type::SameSource) ? hp.sources[i] : hp.starts[i];
        const Vertex snk = hp.ends[i];

        const uint32_t comp = graph_.node_connectivity_id(src);
        uint32_t beg = graph_.vertex_topo_rank(src);
        uint32_t end = graph_.vertex_topo_rank(snk);

        if (comp == UINT32_MAX || beg == UINT32_MAX || end == UINT32_MAX) continue;
        if (beg > end) std::swap(beg, end);

        arms.push_back({comp, beg, end});

        const Vertex rev_src(Vertex::get_reverse(src));
        const Vertex rev_snk(Vertex::get_reverse(snk));

        const uint32_t rev_comp = graph_.node_connectivity_id(rev_snk);
        uint32_t rev_beg = graph_.vertex_topo_rank(rev_snk);
        uint32_t rev_end = graph_.vertex_topo_rank(rev_src);

        if (rev_comp == UINT32_MAX || rev_beg == UINT32_MAX || rev_end == UINT32_MAX) continue;
        if (rev_beg > rev_end) std::swap(rev_beg, rev_end);

        arms.push_back({rev_comp, rev_beg, rev_end});
    }

    if (arms.size() < 2) return;

    for (size_t i = 0; i + 1 < arms.size(); ++i) {
        for (size_t j = i + 1; j < arms.size(); ++j) {
            if (hp.type == Type::DiffSource && arms[i].comp == arms[j].comp) continue;

            uint32_t a = arms[i].comp;
            uint32_t b = arms[j].comp;
            if (a > b) std::swap(a, b);

            insert_homo_used_interval_(
                used_intervals,
                HomoUsedKey_{a, b, arms[i].comp},
                HomoUsedInterval_{arms[i].beg, arms[i].end}
            );

            insert_homo_used_interval_(
                used_intervals,
                HomoUsedKey_{a, b, arms[j].comp},
                HomoUsedInterval_{arms[j].beg, arms[j].end}
            );
        }
    }
}

std::vector<HomoTaskGroup_> GfaBubbleFinder::build_homo_tasks_(
    const std::vector<NodeGroupEntry>& groups, 
    bool no_grouping
) const
{
    log_stream() << "Collecting homo task groups ...\n";

    std::vector<HomoTask_> same_tasks;
    same_tasks.reserve(sources_.size());

    std::vector<HomoTask_> diff_tasks;
    diff_tasks.reserve(groups.size());

    for (Vertex src : sources_) {
        same_tasks.push_back(HomoTask_{
            Type::SameSource,
            graph_.vertex_topo_rank(src),
            src.vertex_id(),
            {src}
        });
    }

    for (const NodeGroupEntry& g : groups) {
        if (g.vtxs.size() < 2) continue;

        uint32_t best_rank = UINT32_MAX;
        uint32_t best_vid = UINT32_MAX;

        for (Vertex v : g.vtxs) {
            const uint32_t r = graph_.vertex_topo_rank(v);
            const uint32_t vid = v.vertex_id();

            if (r < best_rank || (r == best_rank && vid < best_vid)) {
                best_rank = r;
                best_vid = vid;
            }
        }

        diff_tasks.push_back(HomoTask_{
            Type::DiffSource,
            best_rank,
            best_vid,
            g.vtxs
        });
    }

    auto task_comps = [&](const HomoTask_& t) {
        std::vector<uint32_t> comps;
        comps.reserve(t.srcs.size());

        for (Vertex v : t.srcs) {
            const uint32_t c = graph_.node_connectivity_id(v);
            if (c != UINT32_MAX) comps.push_back(c);
        }

        std::sort(comps.begin(), comps.end());
        comps.erase(std::unique(comps.begin(), comps.end()), comps.end());

        return comps;
    };

    auto refresh_group_meta = [](HomoTaskGroup_& g) {
        g.topo_rank = UINT32_MAX;
        g.tie_id = UINT32_MAX;

        for (const HomoTask_& t : g.tasks) {
            g.topo_rank = std::min(g.topo_rank, t.topo_rank);
            g.tie_id = std::min(g.tie_id, t.tie_id);
        }

        std::sort(g.comps.begin(), g.comps.end());
        g.comps.erase(std::unique(g.comps.begin(), g.comps.end()), g.comps.end());

        std::sort(g.tasks.begin(), g.tasks.end(),
            [](const HomoTask_& a, const HomoTask_& b) {
                if (a.topo_rank != b.topo_rank) return a.topo_rank < b.topo_rank;
                if (a.type != b.type) return a.type == Type::SameSource;
                return a.tie_id < b.tie_id;
            }
        );
    };

    std::vector<HomoTaskGroup_> out;
    out.reserve(same_tasks.size() + diff_tasks.size());

    auto push_single_task_group = [&](const HomoTask_& t) {
        HomoTaskGroup_ g;
        g.tasks.push_back(t);
        g.comps = task_comps(t);

        if (g.comps.empty()) {
            g.comps.push_back(UINT32_MAX);
        }

        refresh_group_meta(g);
        out.emplace_back(std::move(g));
    };

    // Same-source tasks
    for (const HomoTask_& t : same_tasks) {
        push_single_task_group(t);
    }

    if (no_grouping) {
        for (const HomoTask_& t : diff_tasks) {
            push_single_task_group(t);
        }
    } else {
        // Diff-source tasks: grouped by overlapping connectivity components.
        const size_t n = diff_tasks.size();

        std::vector<size_t> parent(n);
        std::iota(parent.begin(), parent.end(), 0);

        auto find_parent = [&](size_t x) {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };

        auto unite = [&](size_t a, size_t b) {
            a = find_parent(a);
            b = find_parent(b);
            if (a != b) parent[b] = a;
        };

        std::unordered_map<uint32_t, size_t> comp_owner;
        comp_owner.reserve(n * 4 + 8);

        for (size_t i = 0; i < n; ++i) {
            const std::vector<uint32_t> comps = task_comps(diff_tasks[i]);

            for (uint32_t c : comps) {
                auto it = comp_owner.find(c);

                if (it == comp_owner.end()) {
                    comp_owner.emplace(c, i);
                } else {
                    unite(i, it->second);
                }
            }
        }

        std::unordered_map<size_t, size_t> diff_group_by_root;
        diff_group_by_root.reserve(n * 2 + 8);

        for (size_t i = 0; i < n; ++i) {
            const size_t root = find_parent(i);

            auto it = diff_group_by_root.find(root);
            if (it == diff_group_by_root.end()) {
                const size_t idx = out.size();
                diff_group_by_root.emplace(root, idx);
                out.push_back(HomoTaskGroup_{});
                it = diff_group_by_root.find(root);
            }

            HomoTaskGroup_& g = out[it->second];
            g.tasks.push_back(diff_tasks[i]);

            const std::vector<uint32_t> comps = task_comps(diff_tasks[i]);
            g.comps.insert(g.comps.end(), comps.begin(), comps.end());
        }

        for (const auto& kv : diff_group_by_root) {
            refresh_group_meta(out[kv.second]);
        }
    }

    std::sort(out.begin(), out.end(),
        [](const HomoTaskGroup_& a, const HomoTaskGroup_& b) {
            if (a.topo_rank != b.topo_rank) return a.topo_rank < b.topo_rank;
            return a.tie_id < b.tie_id;
        }
    );

    if (DEBUG_ENABLED) {
        debug_stream() << "  - Homo task groups:\n";

        for (size_t gi = 0; gi < out.size(); ++gi) {
            const HomoTaskGroup_& g = out[gi];
            std::ostringstream oss;

            oss << "    group[" << gi << "] topo_rank=" << g.topo_rank
                << " tie_id=" << g.tie_id
                << " comps:";

            for (uint32_t c : g.comps) {
                oss << " " << c;
            }

            debug_stream() << oss.str() << "\n";
            oss.str("");
            oss.clear();

            for (const HomoTask_& t : g.tasks) {
                oss << "      - type=" << (t.type == Type::SameSource ? "Same" : "Diff") << " topo_rank=" << t.topo_rank << " tie_id=" << t.tie_id << " srcs:";

                for (const Vertex& v : t.srcs) {
                    oss << " " << graph_.getNodeName(v.segment_id()) << (v.is_reverse() ? "-" : "+");
                }

                debug_stream() << oss.str() << "\n";
                oss.str("");
                oss.clear();
            }
        }

        debug_stream() << "\n";
    }

    size_t total_tasks = 0;
    for (const auto& g : out) {
        total_tasks += g.tasks.size();
    }

    log_stream() << "  - Groups: " << out.size() << "\n";
    log_stream() << "  - Tasks: " << total_tasks << "\n\n";

    return out;
}

HomoRunResult_ GfaBubbleFinder::find_homologous_paths_run_(
    const HomoTaskGroup_& group,
    const HomologousParam& same_params,
    const HomologousParam& diff_params,
    const std::unordered_set<uint64_t>& bubble_pairs,
    const std::vector<minimizerdna::Sketch>& node_sketches
) const
{
    HomoRunResult_ out;
    out.paths.reserve(group.tasks.size());

    HomoUsedIntervals_ used_intervals;
    used_intervals.reserve(group.tasks.size() * 8 + 64);

    auto keep_result = [&](HomologousPath&& hr, Type task_type) {
        if (hr.sources.empty()) return;
        if (hr.starts.size() < 2) return;
        if (hr.ends.size() != hr.starts.size()) return;
        if (hr.paths.empty()) return;
        if (homo_path_already_used_(hr, used_intervals)) return;

        if (task_type == Type::SameSource) ++out.same_source_num;
        else ++out.diff_source_num;

        mark_homo_used_intervals_(used_intervals, hr);
        out.paths.emplace_back(std::move(hr));
    };

    for (const HomoTask_& task : group.tasks) {
        ++out.done_tasks;

        if (homo_vertices_already_used_(task.srcs, used_intervals)) {
            continue;
        }

        const HomologousParam& params = (task.type == Type::SameSource) ? same_params : diff_params;

        HomologousPath hr = detect_homologous_paths_from_sources_(
            task.srcs,
            params,
            bubble_pairs,
            node_sketches,
            used_intervals
        );

        keep_result(std::move(hr), task.type);
    }

    return out;
}


void GfaBubbleFinder::find_homologous_paths()
{
    /**
     * @brief Change homo_sameSource_Param.min_source_len to diff_min_src_
     * @date 2026-04-15
     * @version 0.1.3
     * @details Issue identified by WYF
     */
    HomologousParam homo_sameSource_Param(diff_min_src_, same_sim_, same_min_num_, same_min_len_, Type::SameSource);
    HomologousParam homo_diffSource_Param(diff_min_src_, diff_sim_, diff_min_num_, diff_min_len_, Type::DiffSource);

    std::vector<minimizerdna::Sketch> node_sketches = build_all_node_mm_library_();

    std::vector<Vertex> sources_for_homo = collect_sources_for_homo_();

    std::vector<NodeGroupEntry> groups = cluster_sources_by_mm_(
        sources_for_homo,
        node_sketches,
        homo_diffSource_Param.min_source_len,
        homo_diffSource_Param.min_similarity
    );
    const std::vector<HomoTaskGroup_> task_groups = build_homo_tasks_(groups, true);

    const std::unordered_set<uint64_t> bubble_pairs = build_bubble_branch_pairs_();

    log_stream() << "Detecting homologous paths ...\n";

    homologous_paths_.clear();

    size_t total_tasks = 0;
    for (const HomoTaskGroup_& g : task_groups) {
        total_tasks += g.tasks.size();
    }

    ProgressTracker prog(total_tasks);

    const unsigned threads = std::max(1u, thread_);
    ThreadPool pool(threads);

    std::vector<std::future<HomoRunResult_>> futs;
    futs.reserve(task_groups.size());

    for (const HomoTaskGroup_& group : task_groups) {
        futs.emplace_back(
            pool.submit([&, group]() -> HomoRunResult_ {
                return find_homologous_paths_run_(
                    group,
                    homo_sameSource_Param,
                    homo_diffSource_Param,
                    bubble_pairs,
                    node_sketches
                );
            })
        );
    }

    uint64_t same_source_num = 0;
    uint64_t diff_source_num = 0;

    for (auto& f : futs) {
        HomoRunResult_ r = f.get();

        for (size_t i = 0; i < r.done_tasks; ++i) {
            prog.hit();
        }

        same_source_num += r.same_source_num;
        diff_source_num += r.diff_source_num;

        for (HomologousPath& hp : r.paths) {
            homologous_paths_.emplace_back(std::move(hp));
        }
    }

    pool.stop();

    auto hp_size = [&](const HomologousPath& hp) -> uint64_t {
        std::unordered_set<uint32_t> inner_segs;
        inner_segs.reserve(32);

        for (const std::vector<Vertex>& path : hp.paths) {
            for (const Vertex& v : path) {
                inner_segs.insert(v.segment_id());
            }
        }

        uint64_t total = 0;
        for (uint32_t sid : inner_segs) {
            const GfaNode* node = graph_.getNode(sid);
            if (!node || graph_.getNodeDeleted(sid)) continue;
            total += graph_.getNodeLength(sid);
        }

        return total;
    };

    std::sort(homologous_paths_.begin(), homologous_paths_.end(),
        [&](const HomologousPath& a, const HomologousPath& b) {
            const uint64_t la = hp_size(a);
            const uint64_t lb = hp_size(b);
            if (la != lb) return la > lb;
            return a.paths.size() > b.paths.size();
        }
    );

    log_stream() << "  - Same source number: " << same_source_num << "\n";
    log_stream() << "  - Different source number: " << diff_source_num << "\n\n";
}


// Bloom filter
HomologousPath GfaBubbleFinder::detect_homologous_paths_from_sources_(
    const std::vector<Vertex>& srcs,
    const HomologousParam& params,
    const std::unordered_set<uint64_t>& bubble_pairs,
    const std::vector<minimizerdna::Sketch>& node_sketches,
    const HomoUsedIntervals_& used_intervals
) const
{
    HomologousPath out;
    out.type = params.type;


    // ------------------------------------------------ Local variables ------------------------------------------------
    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);

    static constexpr uint8_t PATH_DEAD   = 0;  // stop moving and don't participate in similarity checks
    static constexpr uint8_t PATH_ACTIVE = 1;  // movable
    static constexpr uint8_t PATH_FROZEN = 2;  // stop moving, but still participate in similarity checks

    using PathSketch = fastbloom::BloomSketch<minimizerdna::Sketch>;


    // ------------------------------------------------ Helper functions ------------------------------------------------
    auto fmt = [&](const Vertex& v) -> std::string {
        return graph_.getNodeName(v.segment_id()) + (v.is_reverse() ? "-" : "+");
    };

    auto pair_containments = [](const PathSketch& si, const PathSketch& sj) -> std::pair<double, double> {
        return si.bit_containments(sj);
    };

    auto pair_supported_from_containment = [&](size_t i, size_t j, const std::pair<double, double>& c, const std::vector<uint8_t>& active) -> bool {
        const bool first_good  = (c.first  >= params.min_similarity);
        const bool second_good = (c.second >= params.min_similarity);

        if (active[i] == PATH_FROZEN || active[j] == PATH_FROZEN) {
            return first_good && second_good;
        }

        return first_good || second_good;
    };

    auto pair_can_still_change_flags = [&](
        size_t i, size_t j, const std::vector<uint8_t>& active, const std::vector<uint8_t>& still_supported, const std::vector<uint8_t>& need_move
    ) -> bool {
        if (!still_supported[i] || !still_supported[j]) return true;
        if (active[i] == PATH_FROZEN && active[j] == PATH_FROZEN) return false;
        if (active[i] == PATH_FROZEN && active[j] == PATH_ACTIVE) return !need_move[j];
        if (active[j] == PATH_FROZEN && active[i] == PATH_ACTIVE) return !need_move[i];
        if (active[i] == PATH_ACTIVE && active[j] == PATH_ACTIVE) return !need_move[i] || !need_move[j];
        return false;
    };

    auto pair_cache_index = [](size_t i, size_t j, size_t n) -> size_t {
        if (i > j) std::swap(i, j);
        return i * n + j;
    };

    auto extend_one_step = [&](
        size_t i, PathSketch& target_sketch, uint64_t& target_bp, Vertex& target_cur, Vertex& target_end,
        std::unordered_set<uint32_t>& target_visited, std::unordered_set<uint32_t>& target_visited_seg, 
        std::vector<uint8_t>& active, uint32_t max_extend_steps = 1
    ) -> bool {
        bool moved = false;

        uint32_t step = 0;

        while (step < max_extend_steps) {
            Vertex best(UINT32_MAX);
            const GfaArc* best_arc = nullptr;
            uint32_t best_len = 0;

            for (const GfaArc* a : graph_.getArcsFromVertex(target_cur.vertex_id())) {
                if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;

                const uint32_t w = a->get_target_vertex_id();
                if (w == target_cur.vertex_id()) continue;

                const uint32_t sid = Vertex::get_segment_id(w);

                /**
                 * @brief Fixed cycle issues
                 * @date 2026-04-15
                 * @version 0.1.3
                 * @details Issue identified by WYF
                 */
                if (target_visited_seg.count(sid)) continue;
                if (graph_.getNodeDeleted(sid)) continue;

                const uint32_t L = graph_.getNodeLength(sid);

                const uint32_t ow = (a->ow == INT32_MAX) ? 0u : static_cast<uint32_t>(std::max<int32_t>(0, a->ow));
                const uint32_t best_ow = (!best_arc || best_arc->ow == INT32_MAX) ? 0u : static_cast<uint32_t>(std::max<int32_t>(0, best_arc->ow));

                /**
                 * @brief Candidate priority: overlap length > node length > vertex ID.
                 * @date 2026-05-11
                 * @version 0.1.3
                 */
                if (ow > best_ow || (ow == best_ow && L > best_len) || (ow == best_ow && L == best_len && w < best.vertex_id())) {
                    best = Vertex(w);
                    best_arc = a;
                    best_len = L;
                }
            }

            if (best.vertex_id() == UINT32_MAX) {
                active[i] = PATH_FROZEN;
                target_end = target_cur;
                return moved;
            }

            const uint32_t best_vid = best.vertex_id();

            if (best_vid >= node_sketches.size()) {
                active[i] = PATH_FROZEN;
                target_end = target_cur;
                return moved;
            }

            target_sketch.add_sketch(node_sketches[best_vid]);

            const uint32_t ow = best_arc ? static_cast<uint32_t>(best_arc->ow) : 0;
            target_bp += (best_len > ow) ? (best_len - ow) : 0;

            target_cur = best;
            target_end = best;
            target_visited.insert(best.vertex_id());
            target_visited_seg.insert(best.segment_id());

            moved = true;

            if (best_len >= homo_cnt_len_) {
                ++step;
            }
        }

        return moved;
    };

    auto compute_pair_stats = [&](
        const std::vector<PathSketch>& ss, const std::vector<uint8_t>& active, std::vector<uint8_t>& keep, 
        size_t N, double& min_hi, double& max_hi, bool& any_pair_good
    ) {
        min_hi = 1.0;
        max_hi = 0.0;
        any_pair_good = false;

        std::vector<size_t> alive;
        alive.reserve(N);

        for (size_t i = 0; i < N; ++i) {
            if (active[i] != PATH_DEAD) {
                alive.push_back(i);
            }
        }

        if (alive.size() < 2) {
            return;
        }

        for (size_t ai = 0; ai + 1 < alive.size(); ++ai) {
            const size_t i = alive[ai];

            for (size_t aj = ai + 1; aj < alive.size(); ++aj) {
                const size_t j = alive[aj];

                const auto c = pair_containments(ss[i], ss[j]);
                const double hi = std::max(c.first, c.second);

                if (hi < min_hi) min_hi = hi;
                if (hi > max_hi) max_hi = hi;

                if (pair_supported_from_containment(i, j, c, active)) {
                    keep[i] = 1;
                    keep[j] = 1;
                    any_pair_good = true;
                }
            }
        }

        if (DEBUG_ENABLED) {
            debug_stream() << "  - similarity:\n";
            debug_stream() << "    - min: " << min_hi << "\n";
            debug_stream() << "    - max: " << max_hi << "\n";
            debug_stream() << "    - any_pair_good: " << any_pair_good << "\n";
        }
    };


    // ------------------------------------------------ Initial checks ------------------------------------------------
    if (srcs.empty()) return out;
    if (V == 0 || srcs.front().vertex_id() >= V) return out;

    if (DEBUG_ENABLED) {
        debug_stream() << "Detect homologous paths from sources ...\n";
        debug_stream() << "  - srcs:\n";
        for (const auto& s : srcs) {
            debug_stream() << "    - " << fmt(s) << "\n";
        }
    }


    // ------------------------------------------------ Prepare starts ------------------------------------------------
    std::vector<Vertex> starts;
    std::vector<Vertex> owners;
    starts.reserve(std::max<size_t>(8, srcs.size()));
    owners.reserve(std::max<size_t>(8, srcs.size()));

    if (srcs.size() == 1) {  // Same source
        const Vertex src = srcs.front();

        if (graph_.getNodeDeleted(src.segment_id())) return out;

        std::unordered_set<uint32_t> seen_seg;
        seen_seg.reserve(8);

        for (const GfaArc* a : graph_.getArcsFromVertex(src.vertex_id())) {
            if (!a || a->get_del() || (skip_comp_ && a->get_comp())) continue;

            Vertex w(a->get_target_vertex_id());
            if (w == src) continue;
            if (graph_.getNodeDeleted(w.segment_id())) continue;

            uint32_t sid = w.segment_id();
            if (!seen_seg.insert(sid).second) continue;

            starts.push_back(w);
            owners.push_back(src);
        }

        std::vector<uint32_t> src_group;
        src_group.reserve(starts.size() + 1);
        src_group.push_back(src.vertex_id());

        std::unordered_set<uint32_t> seen_start_vtx;
        seen_start_vtx.reserve(starts.size() * 2 + 4);

        for (const auto& w : starts) {
            if (seen_start_vtx.insert(w.vertex_id()).second) {
                src_group.push_back(w.vertex_id());
            }
        }

        bool have_new_pairs = true;
        if (src_group.size() >= 2) {
            uint64_t h = hash_vec_(src_group, false);
            have_new_pairs = (bubble_pairs.find(h) == bubble_pairs.end());
        }

        if (!have_new_pairs) {
            if (DEBUG_ENABLED) debug_stream() << "  - return: bubble_pairs hit\n\n";
            return out;
        }
    } else {  // Different sources
        std::unordered_set<uint32_t> seen_seg;
        seen_seg.reserve(srcs.size() * 2 + 8);

        for (const Vertex& s : srcs) {
            if (s.vertex_id() >= V) continue;
            if (graph_.getNodeDeleted(s.segment_id())) continue;
            if (!seen_seg.insert(s.segment_id()).second) continue;

            starts.push_back(s);
            owners.push_back(s);
        }
    }

    if (DEBUG_ENABLED) {
        debug_stream() << "  - starts:\n";
        for (const auto& s : starts) {
            debug_stream() << "    - " << fmt(s) << "\n";
        }
    }


    // ------------------------------------------------ Compressed starts according to bubble_pairs ------------------------------------------------
    if (!bubble_pairs.empty() && starts.size() >= 2) {
        const size_t M = starts.size();

        std::vector<size_t> parent(M);
        std::iota(parent.begin(), parent.end(), 0);

        auto find_parent = [&](size_t x) {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };

        auto unite = [&](size_t a, size_t b) {
            a = find_parent(a);
            b = find_parent(b);
            if (a != b) parent[b] = a;
        };

        for (size_t i = 0; i < M; ++i) {
            for (size_t j = i + 1; j < M; ++j) {
                std::vector<uint32_t> pair_vs{
                    starts[i].vertex_id(),
                    starts[j].vertex_id()
                };

                if (bubble_pairs.count(hash_vec_(pair_vs, false))) {
                    unite(i, j);
                }
            }
        }

        std::unordered_map<size_t, size_t> rep;
        rep.reserve(M);

        for (size_t i = 0; i < M; ++i) {
            const size_t r = find_parent(i);
            auto it = rep.find(r);

            if (it == rep.end()) {
                rep.emplace(r, i);
                continue;
            }

            const uint32_t old_len = graph_.getNodeLength(starts[it->second].segment_id());
            const uint32_t new_len = graph_.getNodeLength(starts[i].segment_id());

            if (new_len > old_len) {
                it->second = i;
            }
        }

        if (rep.size() < M) {
            std::vector<size_t> keep_idx;
            keep_idx.reserve(rep.size());

            for (const auto& kv : rep) {
                keep_idx.push_back(kv.second);
            }

            std::sort(keep_idx.begin(), keep_idx.end());

            std::vector<Vertex> new_starts;
            std::vector<Vertex> new_owners;
            new_starts.reserve(keep_idx.size());
            new_owners.reserve(keep_idx.size());

            for (size_t i : keep_idx) {
                new_starts.push_back(starts[i]);
                new_owners.push_back(owners[i]);
            }

            if (DEBUG_ENABLED) {
                debug_stream() << "  - compressed starts by bubble_pairs: " << M << " -> " << new_starts.size() << "\n";
                for (const auto& s : new_starts) {
                    debug_stream() << "    - " << fmt(s) << "\n";
                }
            }

            starts.swap(new_starts);
            owners.swap(new_owners);
        }
    }

    if (starts.size() < 2) {
        if (DEBUG_ENABLED) debug_stream() << "  - return: invalid starts size\n\n";
        return out;
    }

    if (homo_vertices_already_used_(starts, used_intervals)) {
        if (DEBUG_ENABLED) {
            debug_stream() << "  - return: homologous topo interval already used\n\n";
        }
        return out;
    }


    // ------------------------------------------------ Initialize path states ------------------------------------------------
    const size_t N = starts.size();

    std::vector<Vertex> cur(N, Vertex(UINT32_MAX));
    std::vector<Vertex> ends(N, Vertex(UINT32_MAX));
    std::vector<uint64_t> acc_bp(N, 0);
    std::vector<uint8_t> active(N, PATH_ACTIVE);
    std::vector<uint8_t> keep(N, 0);
    std::vector<std::unordered_set<uint32_t>> visited(N);
    std::vector<std::unordered_set<uint32_t>> visited_seg(N);

    // Determine if one source can contain another; if so, the other source doesn't need separate exploration.
    std::unordered_map<Vertex, size_t> source_idx;
    if (srcs.size() > 1) {
        source_idx.reserve(N);
        for (size_t i = 0; i < N; ++i) {
            source_idx.emplace(starts[i], i);
        }
    }


    // ------------------------------------------------ Bloom filter initialization ------------------------------------------------
    fastbloom::BloomOptions bloom_opt;
    bloom_opt.bits = homo_bloom_bits_;
    bloom_opt.hashes = homo_bloom_hash_;
    bloom_opt.use_simd = true;

    std::vector<PathSketch> sketches;
    sketches.reserve(N);

    for (size_t i = 0; i < N; ++i) {
        const Vertex v = starts[i];

        cur[i] = v;
        ends[i] = v;
        visited[i].insert(v.vertex_id());
        visited_seg[i].insert(v.segment_id());

        const uint32_t vid = v.vertex_id();
        if (vid >= node_sketches.size()) {
            if (DEBUG_ENABLED) debug_stream() << "  - return: invalid seed sketch index " << fmt(v) << "\n\n";
            return HomologousPath();
        }

        acc_bp[i] = graph_.getNodeLength(v.segment_id());
        sketches.emplace_back(node_sketches[vid], bloom_opt);
    }


    // ------------------------------------------------ Initial bootstrap steps ------------------------------------------------
    for (size_t i = 0; i < N; ++i) {
        if (active[i] != PATH_ACTIVE) continue;

        while (active[i] == PATH_ACTIVE && acc_bp[i] < homo_boot_bp_) {
            if (!extend_one_step(i, sketches[i], acc_bp[i], cur[i], ends[i], visited[i], visited_seg[i], active, 1)) {
                break;
            }
        }
    }


    // ------------------------------------------------ Check pair stats after bootstrap ------------------------------------------------
    bool any_pair_good = false;
    compute_pair_stats(sketches, active, keep, N, out.min_hi_similarity, out.max_hi_similarity, any_pair_good);
    if (!any_pair_good) {
        if (DEBUG_ENABLED) debug_stream() << "  - return: no supported pairs\n\n";
        return HomologousPath();
    }


    // ------------------------------------------------ DFS exploration ------------------------------------------------
    uint32_t dfs_step = 0;
    uint64_t dfs_states = 0;

    std::vector<size_t> alive_idx;  // indices of paths that are not dead (include both active and frozen)
    std::vector<size_t> active_idx; // indices of paths that are active (movable)
    alive_idx.reserve(N);
    active_idx.reserve(N);

    std::vector<uint8_t> need_move(N, 0);       // need to move in the current step
    std::vector<uint8_t> still_supported(N, 0); // still supported in the current step
    std::vector<uint8_t> can_move(N, 0);        // can move in the current step
    std::vector<uint8_t> next_supported(N, 0);  // supported in the next step (after moving)

    std::vector<Vertex> next_vtxs(N, Vertex(UINT32_MAX));
    std::vector<uint64_t> next_bp(N, 0);

    std::vector<std::unique_ptr<PathSketch>> next_sketch(N);
    std::vector<std::unique_ptr<std::unordered_set<uint32_t>>> next_visited(N);
    std::vector<std::unique_ptr<std::unordered_set<uint32_t>>> next_visited_seg(N);

    std::vector<std::pair<double, double>> cur_pair_cache(N * N, {0.0, 0.0});
    std::vector<uint8_t> cur_pair_cache_valid(N * N, 0);

    while (dfs_step < max_depth_) {
        alive_idx.clear();
        active_idx.clear();

        // Collect indices of alive and active paths
        for (size_t i = 0; i < N; ++i) {
            if (active[i] != PATH_DEAD) {
                alive_idx.push_back(i);
            }
            if (active[i] == PATH_ACTIVE) {
                active_idx.push_back(i);
            }
        }

        // Only one path left or no movable paths
        if (alive_idx.size() < 2) {
            if (DEBUG_ENABLED) debug_stream() << "    - break: alive < 2\n";
            break;
        }
        // No paths can move
        if (active_idx.empty()) {
            if (DEBUG_ENABLED) debug_stream() << "    - break: no movable\n";
            break;
        }

        if (++dfs_states > dfs_guard_) break;

        // Reset
        std::fill(need_move.begin(), need_move.end(), 0);
        std::fill(still_supported.begin(), still_supported.end(), 0);
        std::fill(can_move.begin(), can_move.end(), 0);
        std::fill(next_supported.begin(), next_supported.end(), 0);
        std::fill(cur_pair_cache_valid.begin(), cur_pair_cache_valid.end(), 0);

        for (size_t i = 0; i < N; ++i) {
            next_vtxs[i] = Vertex(UINT32_MAX);
            next_bp[i] = 0;
            next_sketch[i].reset();
            next_visited[i].reset();
            next_visited_seg[i].reset();
        }

        bool any_need_move = false;

        // Similarity check for current state and determine which paths need to move
        for (size_t ai = 0; ai + 1 < alive_idx.size(); ++ai) {
            const size_t i = alive_idx[ai];

            for (size_t aj = ai + 1; aj < alive_idx.size(); ++aj) {
                const size_t j = alive_idx[aj];

                if (!pair_can_still_change_flags(i, j, active, still_supported, need_move)) continue;

                const auto c = pair_containments(sketches[i], sketches[j]);

                const size_t cache_idx = pair_cache_index(i, j, N);
                cur_pair_cache[cache_idx] = c;
                cur_pair_cache_valid[cache_idx] = 1;

                const bool first_good  = (c.first  >= params.min_similarity);
                const bool second_good = (c.second >= params.min_similarity);
                const bool supported = pair_supported_from_containment(i, j, c, active);

                if (supported) {
                    still_supported[i] = 1;
                    still_supported[j] = 1;
                }

                if (active[i] == PATH_FROZEN && active[j] == PATH_ACTIVE) {
                    if (supported) {
                        need_move[j] = 1;
                        any_need_move = true;
                    }
                    continue;
                }

                if (active[j] == PATH_FROZEN && active[i] == PATH_ACTIVE) {
                    if (supported) {
                        need_move[i] = 1;
                        any_need_move = true;
                    }
                    continue;
                }

                if (active[i] == PATH_FROZEN && active[j] == PATH_FROZEN) {
                    continue;
                }

                if (first_good && !second_good) {
                    need_move[i] = 1;
                    any_need_move = true;
                } else if (!first_good && second_good) {
                    need_move[j] = 1;
                    any_need_move = true;
                } else if (first_good && second_good) {
                    need_move[i] = 1;
                    need_move[j] = 1;
                    any_need_move = true;
                }
            }
        }

        // Paths that are not supported by any pair
        for (size_t i : alive_idx) {
            if (!still_supported[i]) {
                active[i] = PATH_DEAD;
            }
        }

        if (!any_need_move) break;

        bool any_move = false;

        // Extend paths that need to move
        for (size_t i : active_idx) {
            if (active[i] != PATH_ACTIVE || !need_move[i]) continue;

            Vertex tmp_cur = cur[i];
            Vertex tmp_end = ends[i];
            uint64_t tmp_bp = acc_bp[i];

            auto tmp_sketch = std::make_unique<PathSketch>(sketches[i]);
            auto tmp_visited = std::make_unique<std::unordered_set<uint32_t>>(visited[i]);
            auto tmp_visited_seg = std::make_unique<std::unordered_set<uint32_t>>(visited_seg[i]);

            if (!extend_one_step(i, *tmp_sketch, tmp_bp, tmp_cur, tmp_end, *tmp_visited, *tmp_visited_seg, active, 1)) {
                ends[i] = cur[i];
                continue;
            }

            next_vtxs[i] = tmp_cur;
            next_bp[i] = tmp_bp;

            next_sketch[i] = std::move(tmp_sketch);
            next_visited[i] = std::move(tmp_visited);
            next_visited_seg[i] = std::move(tmp_visited_seg);

            can_move[i] = 1;
            any_move = true;
        }

        if (!any_move) break;

        // Collect indices of alive and active paths after moving step
        alive_idx.clear();
        active_idx.clear();
        for (size_t i = 0; i < N; ++i) {
            if (active[i] != PATH_DEAD) {
                alive_idx.push_back(i);
            }
            if (active[i] == PATH_ACTIVE) {
                active_idx.push_back(i);
            }
        }

        // Similarity check for next state and determine which paths are supported in the next state
        double next_min_hi = 1.0, next_max_hi = 0.0;
        for (size_t ai = 0; ai + 1 < alive_idx.size(); ++ai) {
            const size_t i = alive_idx[ai];

            for (size_t aj = ai + 1; aj < alive_idx.size(); ++aj) {
                const size_t j = alive_idx[aj];

                std::pair<double, double> c;

                if (!can_move[i] && !can_move[j]) {
                    const size_t cache_idx = pair_cache_index(i, j, N);

                    if (cur_pair_cache_valid[cache_idx]) {
                        c = cur_pair_cache[cache_idx];
                    } else {
                        c = pair_containments(sketches[i], sketches[j]);
                        cur_pair_cache[cache_idx] = c;
                        cur_pair_cache_valid[cache_idx] = 1;
                    }
                } else {
                    const auto& si = can_move[i] ? *next_sketch[i] : sketches[i];
                    const auto& sj = can_move[j] ? *next_sketch[j] : sketches[j];

                    c = pair_containments(si, sj);
                }

                const double hi = std::max(c.first, c.second);

                if (pair_supported_from_containment(i, j, c, active)) {
                    next_supported[i] = 1;
                    next_supported[j] = 1;
                }

                if (hi < next_min_hi) next_min_hi = hi;
                if (hi > next_max_hi) next_max_hi = hi;
            }
        }

        // Check if enough supported paths after moving step
        size_t next_alive_supported_cnt = 0;
        size_t next_movable_supported_cnt = 0;
        for (size_t i : alive_idx) {
            if (next_supported[i]) {
                ++next_alive_supported_cnt;

                if (active[i] == PATH_ACTIVE) {
                    ++next_movable_supported_cnt;
                }
            }
        }
        if (next_alive_supported_cnt < 2) break;
        if (next_movable_supported_cnt == 0) break;

        // Update states for next step
        out.min_hi_similarity = next_min_hi;
        out.max_hi_similarity = next_max_hi;
        for (size_t i : active_idx) {
            if (!can_move[i]) continue;

            cur[i] = next_vtxs[i];
            ends[i] = next_vtxs[i];

            visited[i] = std::move(*next_visited[i]);
            visited_seg[i] = std::move(*next_visited_seg[i]);

            sketches[i] = std::move(*next_sketch[i]);
            acc_bp[i] = next_bp[i];
        }

        // If one source can contain another, the contained one doesn't need to continue; if they are still similar after moving, stop the contained one.
        if (srcs.size() > 1) {
            for (size_t i : alive_idx) {
                if (active[i] == PATH_DEAD) continue;

                auto it = source_idx.find(cur[i]);
                if (it == source_idx.end()) continue;

                const size_t j = it->second;
                if (j == i || active[j] == PATH_DEAD) continue;

                if (acc_bp[i] >= acc_bp[j]) {
                    active[j] = PATH_DEAD;
                    keep[j] = 0;
                }
            }
        }

        // Paths that are not supported by any pair in the next state
        for (size_t i : alive_idx) {
            if (active[i] != PATH_DEAD && !next_supported[i]) {
                active[i] = PATH_DEAD;
            }
        }

        // Check whether still-alive paths intersect a node, if so, stop one of the pairs.
        if (params.type == Type::SameSource) {
            active_idx.clear();

            for (size_t i = 0; i < N; ++i) {
                if (active[i] == PATH_ACTIVE) {
                    active_idx.push_back(i);
                }
            }

            if (active_idx.size() >= 2) {
                int16_t active_count = static_cast<int16_t>(active_idx.size());

                for (size_t x = 0; x + 1 < active_idx.size(); ++x) {
                    const size_t i = active_idx[x];
                    if (active[i] != PATH_ACTIVE) continue;

                    for (size_t y = x + 1; y < active_idx.size(); ++y) {
                        const size_t j = active_idx[y];
                        if (active[j] != PATH_ACTIVE) continue;

                        const auto* small = &visited[i];
                        const auto* large = &visited[j];

                        if (small->size() > large->size()) {
                            std::swap(small, large);
                        }

                        bool found = false;
                        uint32_t best_vid = UINT32_MAX;
                        uint32_t best_rank = UINT32_MAX;

                        for (uint32_t vid : *small) {
                            if (!large->count(vid)) continue;

                            const uint32_t rank = graph_.vertex_topo_rank(Vertex(vid));
                            if (rank == UINT32_MAX) continue;

                            if (!found || rank < best_rank ||
                                (rank == best_rank && vid < best_vid)) {
                                found = true;
                                best_vid = vid;
                                best_rank = rank;
                            }
                        }

                        if (!found) continue;

                        if (active_count == 2) {
                            ends[i] = Vertex(best_vid);
                            ends[j] = Vertex(best_vid);
                            active[i] = PATH_DEAD;
                            active[j] = PATH_DEAD;

                            if (DEBUG_ENABLED) {
                                debug_stream() << "    - stop final intersected pair: " << fmt(starts[i]) << " -> " << fmt(ends[i]) << "  and " << fmt(starts[j]) << " -> " << fmt(ends[j]) << "  rank=" << best_rank << "\n";
                            }

                            break;
                        }

                        ends[j] = Vertex(best_vid);
                        active[j] = PATH_DEAD;
                        active_count--;

                        if (DEBUG_ENABLED) {
                            debug_stream() << "    - stop intersected branch: " << fmt(starts[j]) << " -> " << fmt(ends[j]) << "  by " << fmt(starts[i]) << "  rank=" << best_rank << "\n";
                        }
                    }
                }
            }
        }

        ++dfs_step;
    }

done:
    out.sources.clear();
    out.starts.clear();
    out.ends.clear();
    out.lens.clear();
    out.paths.clear();

    std::vector<std::unordered_set<uint32_t>> visited_sets;
    visited_sets.reserve(N);

    for (size_t i = 0; i < N; ++i) {
        if (!keep[i]) continue;
        if (acc_bp[i] < params.min_path_bp) continue;

        out.sources.push_back(owners[i]);
        out.starts.push_back(starts[i]);
        out.ends.push_back(ends[i]);
        out.lens.push_back(acc_bp[i]);
        /**
         * @brief Fixed index mismatch issues
         * @date 2026-04-27
         * @version 0.1.3
         * @details Issue identified by WYF
         */
        visited_sets.push_back(std::move(visited[i]));
    }

    if (DEBUG_ENABLED) {
        debug_stream() << "  - kept:\n";
        for (size_t i = 0; i < N; ++i) {
            if (!keep[i]) continue;
            debug_stream() << "    - " << fmt(starts[i]) << " -> " << fmt(ends[i]) << " len: " << acc_bp[i] << "\n";
        }
    }

    out.paths.clear();
    out.paths.reserve(out.starts.size());

    for (size_t i = 0; i < out.starts.size(); ++i) {
        std::unordered_set<uint32_t> dummy_visited;
        std::unordered_set<uint32_t>& visited_set = (homo_num_ == 1) ? visited_sets[i] : dummy_visited;

        uint32_t eff_max_depth = max_depth_;
        if (!visited_set.empty()) {
            uint32_t limit = visited_set.size() + DEPTH_MARGIN_;
            if (eff_max_depth == 0 || eff_max_depth > limit) {
                eff_max_depth = limit;
            }
        }

        bool hit_limits = false;

        std::vector<std::vector<uint32_t>> paths = graph_.enumerate_paths_greedy_DFS(
            out.starts[i].vertex_id(), out.ends[i].vertex_id(), visited_set,
            eff_max_depth, homo_num_, skip_comp_, hit_limits, dfs_guard_, stall_round_limit_
        );

        for (const std::vector<uint32_t>& path : paths) {
            std::vector<Vertex> vpath;
            vpath.reserve(path.size());

            if (out.type == Type::SameSource) {
                vpath.push_back(out.sources[i]);
            }

            for (uint32_t vid : path) {
                vpath.emplace_back(Vertex(vid));
            }

            if (vpath.size() < params.min_path_nodes) continue;

            out.paths.push_back(std::move(vpath));
        }

        if (DEBUG_ENABLED) {
            debug_stream() << "  - enumerate:\n";
            debug_stream() << "    - " << fmt(out.starts[i]) << " -> " << fmt(out.ends[i]) << " paths: " << paths.size() << "\n";

            for (size_t p = 0; p < paths.size(); ++p) {
                std::string s;
                for (uint32_t vid : paths[p]) {
                    s += fmt(Vertex(vid)) + " ";
                }
                debug_stream() << "      - " << s << "\n";
            }
        }
    }

    if (DEBUG_ENABLED) debug_stream() << "\n";

    return out;
}

void GfaBubbleFinder::print_homologous_paths() const
{
    auto vname = [&](Vertex v) -> std::string {
        return graph_.getNodeName(v.segment_id()) + (v.is_reverse() ? "-" : "+");
    };

    for (const auto& hp : homologous_paths_) {
        if (hp.sources.empty()) continue;

        std::string src_string, sink_string;
        src_string.reserve(hp.sources.size() * 24);
        sink_string.reserve(hp.ends.size() * 24);

        uint32_t source_num;
        if (hp.type == Type::SameSource) {
            source_num = 1;
        } else {
            source_num = static_cast<uint32_t>(hp.sources.size());
        }
        for (std::uint32_t j = 0; j < source_num; ++j) {
            if (j) src_string += ",";
            src_string += graph_.getNodeName(hp.sources[j].segment_id());
            src_string += (hp.sources[j].is_reverse() ? "-" : "+");
        }
        for (std::uint32_t j = 0; j < hp.ends.size(); ++j) {
            if (j) sink_string += ",";
            sink_string += graph_.getNodeName(hp.ends[j].segment_id());
            sink_string += (hp.ends[j].is_reverse() ? "-" : "+");
        }

        std::unordered_set<uint32_t> inner_segs;
        inner_segs.reserve(16);

        std::size_t total_bp = 0;

        for (const std::vector<Vertex>& path : hp.paths) {
            for (const Vertex& v : path) {
                uint32_t seg = v.segment_id();
                if (!inner_segs.insert(seg).second) continue;

                const GfaNode* node = graph_.getNode(seg);
                if (!node || graph_.getNodeDeleted(seg)) continue;
                total_bp += graph_.getNodeLength(seg);
            }
        }

        log_stream() << "Homologous: source=" << src_string << "  sink=" << sink_string << "\n";
        log_stream() << "  - Total Nodes Length: " << total_bp << "\n";
        log_stream() << "  - Total Nodes Count: " << inner_segs.size() << "\n";
        log_stream() << "  - Similarity (Min/Max): " << hp.min_hi_similarity << "/" << hp.max_hi_similarity << "\n";

        std::string path_str;
        for (size_t i = 0; i < hp.paths.size(); ++i) {
            path_str = "  - path[" + std::to_string(i) + "]: ";
            for (const Vertex& v : hp.paths[i]) {
                path_str += vname(v) + ' ';
            }
            log_stream() << path_str << "\n";
        }
    }
}

} // namespace GfaBubble
