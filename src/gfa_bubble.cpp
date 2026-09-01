#include "../include/gfa_bubble.hpp"
#include "../include/gfa_bubble_annotation.hpp"
#include "../include/progress_tracker.hpp"
#include "../include/ThreadPool.hpp"

#include <iomanip>
#include <iterator>
#include <numeric>
#include <optional>
#include <set>
#include <type_traits>

namespace GfaBubble {

static std::vector<uint32_t> common_samples(
    const GfaGraph& graph,
    uint32_t first,
    uint32_t second
) {
    const auto& a = graph.getNodeSampleIds(Vertex::get_segment_id(first));
    const auto& b = graph.getNodeSampleIds(Vertex::get_segment_id(second));

    std::vector<uint32_t> out;
    out.reserve(std::min(a.size(), b.size()));
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));
    return out;
}

static bool node_is_vcf_variant(const GfaGraph& graph, uint32_t sid)
{
    const GfaNode* node = graph.getNode(sid);
    if (!node || graph.getNodeDeleted(sid)) return false;

    const uint8_t* tag = GfaAuxParser::get_aux_data(node->aux.aux_data, "VT");
    return tag && *tag == 'Z' && *(tag + 1) != '\0';
}

static bool is_haplotype_sample(const GfaGraph& graph, uint32_t sample)
{
    if (sample >= graph.getSampleNames().size()) return false;
    const std::string& name = graph.getSampleNames()[sample];
    const size_t marker = name.rfind(".hap");
    return marker != std::string::npos && marker + 4 < name.size() &&
        std::all_of(name.begin() + marker + 4, name.end(), [](unsigned char c) {
            return std::isdigit(c);
        });
}

GfaBubbleFinder::GfaBubbleFinder(const GfaGraph& graph, BubbleOpts params)
    : params_(std::move(params)), graph_(graph), path_walker_(graph),
      mm_opt_{/*k=*/params_.homo_k, /*w=*/params_.homo_w, /*seek_reverse=*/false, /*seed=*/0x8a5cd789635d2dffULL}
{}

std::vector<uint32_t> GfaBubbleFinder::trim_path_cluster_anchors(const std::vector<uint32_t>& path, Type type) {
    if (type == Type::SameSource) {
        if (path.size() <= 1) return {};
        return std::vector<uint32_t>(path.begin() + 1, path.end());
    }

    if (type == Type::Normal) {
        if (path.size() <= 2) return {};
        return std::vector<uint32_t>(path.begin() + 1, path.end() - 1);
    }

    return path;
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

bool GfaBubbleFinder::is_strict_source_(Vertex vtx, std::vector<Vertex>& tgt) const
{
    tgt.clear();
    tgt.reserve(8);

    std::unordered_set<uint32_t> seen_seg; seen_seg.reserve(8);

    for (const GfaArc* a : graph_.getArcsFromVertex(vtx.vertex_id())) {
        if (!a || a->get_del() || (params_.skip_comp && a->get_comp())) continue;
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
                if (!a || a->get_del() || (params_.skip_comp && a->get_comp())) continue;
                uint32_t t = a->get_target_vertex_id();
                if (t == u) continue;
                if (++fanout > 1) return false;
                next = t;
            }
            if (fanout == 0) return true;  // Dead end (0.1.3-r12 2026/06/29)
            u = next;
            if (u == target) return true;
        }
        return (u == target);
    };

    constexpr int MAX_TUNNEL_STEPS = 4;

    for (const GfaArc* in : graph_.getArcsToVertex(w.vertex_id())) {
        if (!in || in->get_del() || (params_.skip_comp && in->get_comp())) continue;
        const uint32_t pv = in->get_source_vertex_id();
        if (pv == w.vertex_id() || pv == v.vertex_id()) continue;

        int outdeg = 0;
        bool all_funnel_to_w = true;

        for (const GfaArc* a : graph_.getArcsFromVertex(pv)) {
            if (!a || a->get_del() || (params_.skip_comp && a->get_comp())) continue;
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
        const uint32_t sid = NodeHandle::get_segment_id(v);
        if (graph_.getNodeDeleted(sid)) continue;

        // Single-node (Chromosome-like)
        if (graph_.node_connectivity_size(Vertex(v)) == 1) {
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

        sources_.emplace_back(Vertex(v));
    }

    log_stream() << "  - strict sources: " << sources_.size() << ", tips: " << tips_.size() << "\n";

    return;
}

void GfaBubbleFinder::find_bubbles()
{
    log_stream() << "Detecting bubbles ...\n";
    bubbles_.clear();
    bubble_branch_pairs_ = {};

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

    const uint64_t bfs_limit = static_cast<uint64_t>(std::min<std::size_t>(params_.max_depth ? params_.max_depth : 64, UINT64_MAX));

    // Closed bubbles from strict sources
    auto pair_key = [](uint32_t s0, uint32_t s1) -> uint64_t {
        if (s0 > s1) std::swap(s0, s1);
        return (static_cast<uint64_t>(s0) << 32) | static_cast<uint64_t>(s1);
    };

    const unsigned threads = std::max(1u, params_.threads);
    ThreadPool pool(threads);

    std::vector<std::future<Bubble>> futs;
    futs.reserve(sources_.size());

    ProgressTracker prog(sources_.size());

    for (Vertex src : sources_) {
        futs.emplace_back(
            pool.submit([&, src]() -> Bubble {
                Bubble bubble = detect_closed_bubble_from_source_(src, bfs_limit);
                prog.hit();
                return bubble;
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
    if (!params_.keep_nested) filter_nested_bubbles();
    else bubble_branch_pairs_ = build_bubble_branch_pairs_();

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

    auto is_source_segment = [&](Vertex v) -> bool {
        return v.segment_id() == source.segment_id();
    };

    auto is_valid_sink_candidate = [&](Vertex v) -> bool {
        if (v.vertex_id() == source.vertex_id()) return false;
        if (v.segment_id() == source.segment_id()) return false;
        return is_strict_sink_(source, v);
    };

    best_sink = UINT32_MAX;
    uint32_t best_sink_depth = 0;

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

    std::queue<std::tuple<uint32_t, int16_t, uint32_t>> q;
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

        if (is_source_segment(w)) {
            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(2) << "seed " << b << ": " << vtx_name(w) << " returns to source segment, skip sink finder\n";
            }
            return false;
        }

        visited.insert(encode_state(w.vertex_id(), b));

        q.emplace(w.vertex_id(), b, 0u);
        ++branch_active[static_cast<std::size_t>(b)];
        local_nodes.insert(w.vertex_id());

        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "seed " << b << ": " << vtx_name(w) << "\n";
        }

        if (is_valid_sink_candidate(w)) {
            auto& seen_set = branch_sink_seen[static_cast<std::size_t>(b)];
            if (seen_set.insert(w.vertex_id()).second) {
                SinkStat& ss = sink_stats[w.vertex_id()];
                ++ss.seen_count;

                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3) << "candidate sink: " << vtx_name(w) << "  depth=0" << "  branch_seen=" << ss.seen_count << "/" << seeds_cnt << "  max_depth=" << ss.max_depth << "  sum_depth=" << ss.sum_depth << "\n";
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
                            debug_stream() << log_indent(4) << "update best sink: " << vtx_name(w) << "  best_max_depth=" << best_max_depth << "  best_sum_depth=" << best_sum_depth << "\n";
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
            debug_stream() << log_indent(2) << "visit " << vtx_name(Vertex(cur)) << "  depth=" << d << "  branch=" << b << "  qsize=" << q.size() << "\n";
        }

        if (is_source_segment(Vertex(cur))) {
            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(3) << "drop cyclic state: branch " << b << " reaches source segment at " << vtx_name(Vertex(cur)) << "\n";
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
            if (!a || a->get_del() || (params_.skip_comp && a->get_comp())) continue;

            const uint32_t nxt = a->get_target_vertex_id();
            if (nxt == cur) continue;

            const Vertex nxt_v(nxt);
            const uint32_t nxt_sid = nxt_v.segment_id();

            // Do not pass through source segment
            if (nxt_sid == source.segment_id()) {
                if (DEBUG_ENABLED) {
                    debug_stream() << log_indent(3) << "skip cyclic edge: " << vtx_name(Vertex(cur)) << " -> " << vtx_name(nxt_v) << " returns to source segment\n";
                }
                continue;
            }

            const uint64_t st = encode_state(nxt, b);
            if (!visited.insert(st).second) continue;

            local_nodes.insert(nxt);
            q.emplace(nxt, b, d + 1);
            ++branch_active[static_cast<std::size_t>(b)];
            branch_extended = true;

            if (!is_valid_sink_candidate(nxt_v)) continue;

            auto& seen_set = branch_sink_seen[static_cast<std::size_t>(b)];
            if (!seen_set.insert(nxt).second) continue;

            SinkStat& ss = sink_stats[nxt];
            ++ss.seen_count;
            ss.sum_depth += static_cast<uint32_t>(d + 1);
            ss.max_depth = std::max(ss.max_depth, static_cast<uint32_t>(d + 1));

            if (DEBUG_ENABLED) {
                debug_stream() << log_indent(3) << "candidate sink: " << vtx_name(nxt_v) << "  depth=" << (d + 1) << "  branch_seen=" << ss.seen_count << "/" << seeds_cnt << "  max_depth=" << ss.max_depth << "  sum_depth=" << ss.sum_depth << "\n";
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
                        debug_stream() << log_indent(4) << "update best sink: " << vtx_name(nxt_v) << "  best_max_depth=" << best_max_depth << "  best_sum_depth=" << best_sum_depth << "\n";
                    }
                }
            }
        }

        if (best_sink == UINT32_MAX && !branch_extended && branch_active[static_cast<std::size_t>(b)] == 0 && DEBUG_ENABLED) {
            debug_stream() << log_indent(3) << "branch " << b << " has no active nodes; keep scanning other branches for shared sinks\n";
        }
    }

    if (DEBUG_ENABLED && best_sink != UINT32_MAX) {
        debug_stream() << log_indent(2) << "selected sink: " << vtx_name(Vertex(best_sink)) << "  sink_depth=" << best_sink_depth << "\n";
    }

    return best_sink != UINT32_MAX;
}

std::vector<std::vector<uint32_t>> GfaBubbleFinder::enumerate_bubble_paths_(
    uint32_t source,
    uint32_t sink,
    const std::unordered_set<uint32_t>& region_set,
    uint32_t max_depth,
    bool& hit_limits
) const {
    hit_limits = false;

    // ------------------------------------------------ 1. Enumerate paths from P-lines in the GFA graph ------------------------------------------------
    std::vector<std::vector<uint32_t>> guided_paths;
    if (!path_walker_.empty()) {
        guided_paths = path_walker_.paths_between(source, sink, params_.max_paths);
    }

    // ------------------------------------------------ 2. If there are any nodes is from a VCF variant (don't have P-line), so use enumerate paths with sample guidance ------------------------------------------------
    bool vcf_hit_limits = false;
    for (auto& path : enumerate_vcf_paths_(source, sink, region_set, max_depth, vcf_hit_limits)) {
        if (std::find(guided_paths.begin(), guided_paths.end(), path) == guided_paths.end()) {
            guided_paths.push_back(std::move(path));
        }
    }
    hit_limits = hit_limits || vcf_hit_limits;
    if (guided_paths.size() >= 2) return guided_paths;

    // ------------------------------------------------ 3. If don't have any guided paths, use the original graph walk ------------------------------------------------
    const std::vector<uint32_t> samples = common_samples(graph_, source, sink);
    if (samples.empty()) {
        return graph_.enumerate_paths_greedy_DFS(
            source, sink, region_set, max_depth,
            params_.max_paths, params_.skip_comp, hit_limits, params_.DFS_guard, params_.stall_round_limit
        );
    }

    // // ------------------------------------------------ 4. Enumerate paths based on samples guided ------------------------------------------------
    std::vector<std::vector<std::vector<uint32_t>>> candidates;
    candidates.reserve(samples.size());
    size_t rounds = 0;

    for (uint32_t sample : samples) {
        bool sample_hit_limits = false;
        candidates.push_back(graph_.enumerate_paths_greedy_DFS(
            source, sink, region_set, max_depth,
            params_.max_paths, params_.skip_comp, sample_hit_limits, params_.DFS_guard, params_.stall_round_limit, sample
        ));
        hit_limits = hit_limits || sample_hit_limits;
        rounds = std::max(rounds, candidates.back().size());
    }

    std::vector<std::vector<uint32_t>> paths;
    paths.reserve(params_.max_paths);

    for (size_t rank = 0; rank < rounds; ++rank) {
        for (const auto& sample_paths : candidates) {
            if (rank >= sample_paths.size()) continue;
            if (std::find(paths.begin(), paths.end(), sample_paths[rank]) != paths.end()) continue;

            paths.push_back(sample_paths[rank]);
            if (params_.max_paths > 0 && paths.size() >= params_.max_paths) {
                return paths.size() >= 2 ? std::move(paths) : std::vector<std::vector<uint32_t>>{};
            }
        }
    }

    if (paths.size() >= 2) return paths;

    // ------------------------------------------------ 5. No guide produced a complete bubble; keep the original graph walk as ------------------------------------------------
    bool legacy_hit_limits = false;
    std::vector<std::vector<uint32_t>> legacy = graph_.enumerate_paths_greedy_DFS(
        source, sink, region_set, max_depth,
        params_.max_paths, params_.skip_comp, legacy_hit_limits, params_.DFS_guard, params_.stall_round_limit
    );
    hit_limits = hit_limits || legacy_hit_limits;
    return legacy;
}

std::vector<std::vector<uint32_t>> GfaBubbleFinder::enumerate_vcf_paths_(
    uint32_t source,
    uint32_t sink,
    const std::unordered_set<uint32_t>& region_set,
    uint32_t max_depth,
    bool& hit_limits
) const {
    hit_limits = false;

    std::unordered_set<uint32_t> variants;
    std::vector<uint32_t> samples;
    variants.reserve(4);
    samples.reserve(4);

    for (uint32_t vertex : region_set) {
        const uint32_t sid = Vertex::get_segment_id(vertex);
        if (vertex == source || vertex == sink || !node_is_vcf_variant(graph_, sid)) continue;
        if (!variants.insert(sid).second) continue;

        for (uint32_t sample : graph_.getNodeSampleIds(sid)) samples.push_back(sample);
    }

    std::sort(samples.begin(), samples.end());
    samples.erase(std::unique(samples.begin(), samples.end()), samples.end());
    if (variants.empty() || samples.empty()) return {};

    const uint32_t variant_limit = static_cast<uint32_t>(variants.size()) + 1;
    const uint32_t max_paths = params_.max_paths == 0 ? 0 : std::max<uint32_t>(params_.max_paths, variant_limit);

    std::vector<std::vector<uint32_t>> paths;
    for (uint32_t sample : samples) {
        bool sample_hit_limits = false;
        std::vector<std::vector<uint32_t>> sample_paths = graph_.enumerate_paths_greedy_DFS(
            source, sink, region_set, max_depth, max_paths, params_.skip_comp,
            sample_hit_limits, params_.DFS_guard, params_.stall_round_limit, sample
        );
        hit_limits = hit_limits || sample_hit_limits;

        for (auto& path : sample_paths) {
            const bool has_variant = std::any_of(path.begin(), path.end(), [&](uint32_t vertex) {
                return variants.contains(Vertex::get_segment_id(vertex));
            });
            if (!has_variant || std::find(paths.begin(), paths.end(), path) != paths.end()) continue;
            paths.push_back(std::move(path));
        }
    }
    return paths;
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

    // 2. Collect sink candidates (BFS)
    std::unordered_set<uint32_t> local_nodes;
    local_nodes.reserve(256);

    uint32_t sink = UINT32_MAX;

    if (DEBUG_ENABLED) {
        debug_stream() << log_indent(1) << "bfs to detect sink\n";
    }

    if (!find_common_nearest_sink_(v, seed_buf, bfs_limit, sink, local_nodes)) {
        if (DEBUG_ENABLED) {
            debug_stream() << log_indent(2) << "no closed sink found\n";
        }
        return bb;
    }

    // 3. Use BFS sink depth to limit DFS depth
    uint32_t eff_max_depth = params_.max_depth;
    if (local_nodes.size() > 0) {
        uint32_t limit = local_nodes.size() + DEPTH_MARGIN_;
        if (eff_max_depth == 0 || eff_max_depth > limit) {
            eff_max_depth = limit;
        }
    }

    if (DEBUG_ENABLED) {
        debug_stream() << log_indent(1) << "finalize bubble region\n";
        debug_stream() << log_indent(2) << "sink           : " << vtx_name(Vertex(sink)) << "\n";
        debug_stream() << log_indent(2) << "eff max depth  : " << eff_max_depth << "\n";
        debug_stream() << log_indent(2) << "local nodes    : " << local_nodes.size() << "\n";
    }

    // 4. Resolve paths by P-line, fixed sample, or legacy DFS
    bool hit_limits = false;
    auto paths = enumerate_bubble_paths_(
        v.vertex_id(), sink, local_nodes, eff_max_depth, hit_limits
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

void GfaBubbleFinder::filter_nested_bubbles()
{
    log_stream() << "Filtering fully nested bubbles ...\n";
    if (bubbles_.size() < 2) {
        bubble_branch_pairs_ = build_bubble_branch_pairs_();
        return;
    }

    struct BubbleNodes {
        size_t bubble{0};
        std::vector<uint32_t> nodes;
    };

    std::vector<BubbleNodes> ordered;
    ordered.reserve(bubbles_.size());
    for (size_t i = 0; i < bubbles_.size(); ++i) {
        BubbleNodes item;
        item.bubble = i;
        if (bubbles_[i].get_source() != UINT32_MAX) {
            item.nodes.push_back(Vertex::get_segment_id(bubbles_[i].get_source()));
        }
        if (bubbles_[i].get_sink() != UINT32_MAX) {
            item.nodes.push_back(Vertex::get_segment_id(bubbles_[i].get_sink()));
        }
        for (const auto& path : bubbles_[i].get_paths()) {
            for (uint32_t vertex : path) item.nodes.push_back(Vertex::get_segment_id(vertex));
        }
        std::sort(item.nodes.begin(), item.nodes.end());
        item.nodes.erase(std::unique(item.nodes.begin(), item.nodes.end()), item.nodes.end());
        ordered.push_back(std::move(item));
    }

    std::sort(ordered.begin(), ordered.end(), [this](const BubbleNodes& a, const BubbleNodes& b) {
        if (a.nodes.size() != b.nodes.size()) return a.nodes.size() > b.nodes.size();
        if (bubbles_[a.bubble].get_len() != bubbles_[b.bubble].get_len()) {
            return bubbles_[a.bubble].get_len() > bubbles_[b.bubble].get_len();
        }
        return a.bubble < b.bubble;
    });

    std::unordered_map<uint32_t, std::vector<uint32_t>> node_to_kept;
    node_to_kept.reserve(ordered.size() * 4 + 16);
    std::vector<BubbleNodes> kept_index;
    kept_index.reserve(ordered.size());
    std::vector<Bubble> kept;
    kept.reserve(ordered.size());

    for (BubbleNodes& item : ordered) {
        const std::vector<uint32_t>* candidates = nullptr;
        for (uint32_t node : item.nodes) {
            const auto found = node_to_kept.find(node);
            if (found == node_to_kept.end()) {
                candidates = nullptr;
                break;
            }
            if (!candidates || found->second.size() < candidates->size()) candidates = &found->second;
        }

        bool nested = false;
        if (candidates) {
            for (uint32_t id : *candidates) {
                const auto& container = kept_index[id].nodes;
                if (container.size() >= item.nodes.size() &&
                    std::includes(container.begin(), container.end(), item.nodes.begin(), item.nodes.end())) {
                    nested = true;
                    break;
                }
            }
        }
        if (nested) continue;

        const uint32_t id = static_cast<uint32_t>(kept_index.size());
        for (uint32_t node : item.nodes) node_to_kept[node].push_back(id);
        kept_index.push_back(std::move(item));
        kept.push_back(std::move(bubbles_[kept_index.back().bubble]));
    }

    const size_t removed = bubbles_.size() - kept.size();
    bubbles_ = std::move(kept);
    log_stream() << "  - Fully nested bubbles removed: " << removed << "\n";
    log_stream() << "  - Bubbles retained: " << bubbles_.size() << "\n\n";

    bubble_branch_pairs_ = build_bubble_branch_pairs_();
}


GfaBubble::BubbleBranchPairs_ GfaBubbleFinder::build_bubble_branch_pairs_() const
{
    log_stream() << "Building bubble branch groups ...\n";

    GfaBubble::BubbleBranchPairs_ used;
    used.seg_to_groups.reserve(bubbles_.size() * 8 + 16);

    uint32_t group_id = 0;
    uint64_t total_nodes = 0;

    auto add_bubble_group = [&](const std::vector<uint32_t>& segs) {
        if (segs.size() < 2) return;

        const uint32_t gid = group_id++;

        for (uint32_t sid : segs) {
            used.seg_to_groups[sid].push_back(gid);
            ++total_nodes;
        }
    };

    for (const auto& bb : bubbles_) {
        const uint32_t first_vtx = bb.get_source();
        const uint32_t last_vtx  = bb.get_sink();
        const auto& paths = bb.get_paths();

        std::vector<uint32_t> all_segs;
        all_segs.reserve(paths.size() * 4 + 2);

        std::unordered_set<uint32_t> seen_all;
        seen_all.reserve(paths.size() * 8 + 8);

        auto add_all_seg = [&](uint32_t vtx) {
            if (vtx == UINT32_MAX) return;

            const uint32_t sid = Vertex::get_segment_id(vtx);
            if (seen_all.insert(sid).second) {
                all_segs.push_back(sid);
            }
        };

        add_all_seg(first_vtx);
        add_all_seg(last_vtx);

        for (const auto& p : paths) {
            for (uint32_t vtx : p) {
                add_all_seg(vtx);
            }
        }

        add_bubble_group(all_segs);
    }

    for (auto& kv : used.seg_to_groups) {
        auto& xs = kv.second;
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    }

    log_stream() << "  - Built " << group_id << " bubble branch groups with " << total_nodes << " nodes\n\n";

    return used;
}


void GfaBubbleFinder::save_bubble_as_gfa(
    const std::string& output_file,
    const uint32_t min_len,
    const uint32_t min_num,
    bool write_seq,
    const std::string& command_line
) const {
    BubbleWriter(graph_, bubbles_, homologous_paths_).save_gfa(
        output_file, min_len, min_num, write_seq, command_line
    );
}

std::vector<MosaicSite> GfaBubbleFinder::save_bubble_as_vcf(
    const std::string& output_prefix,
    const std::string& paf_file,
    const std::string& ref_file,
    uint32_t min_mapq,
    uint32_t min_aln_len,
    const std::vector<ReferencePath>* reference_paths
) const {
    return BubbleWriter(graph_, bubbles_, homologous_paths_).save_vcf(
        output_prefix, paf_file, ref_file, min_mapq, min_aln_len,
        reference_paths, params_.max_mosaic_freq
    );
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
            if (!a || a->get_del() || (params_.skip_comp && a->get_comp())) continue;
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


std::vector<Vertex> GfaBubbleFinder::collect_sources_for_homo_() const
{
    log_stream() << "Collecting sources for homologous path enumeration ...\n";

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);

    std::vector<Vertex> sources;
    sources.reserve(V);

    for (uint32_t vid = 0; vid < V; ++vid) {
        Vertex v(vid);
        const uint32_t sid = v.segment_id();

        if (graph_.getNodeDeleted(sid)) continue;

        const uint32_t len = graph_.getNodeLength(sid);
        if (len >= params_.diff_min_src) {
            sources.push_back(v);
        }
    }

    if (DEBUG_ENABLED) {
        for (Vertex source : sources) {
            debug_stream() << "  - source: " << graph_.getNodeName(source.segment_id()) << (source.is_reverse() ? "-" : "+") << " " << graph_.getNodeLength(source.segment_id()) << "\n";
        }
    }

    log_stream() << "  - sources: " << sources.size() << "\n\n";

    return sources;
}

std::vector<minimizerdna::Sketch> GfaBubbleFinder::build_source_mm_library_(
    const std::vector<Vertex>& sources
) const
{
    log_stream() << "Building minimizer library for source nodes ...\n";

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    std::vector<minimizerdna::Sketch> node_sketches(V);

    const unsigned threads = std::max(1u, params_.threads);
    ThreadPool pool(threads);

    std::vector<std::future<std::pair<uint32_t, minimizerdna::Sketch>>> futs;
    futs.reserve(sources.size());

    ProgressTracker prog(sources.size());

    for (Vertex v : sources) {
        futs.emplace_back(
            pool.submit([&, v]() -> std::pair<uint32_t, minimizerdna::Sketch> {
                const uint32_t vid = v.vertex_id();
                const uint32_t sid = v.segment_id();

                if (vid >= node_sketches.size()) {
                    prog.hit();
                    return {vid, minimizerdna::Sketch{}};
                }

                if (graph_.getNodeDeleted(sid)) {
                    prog.hit();
                    return {vid, minimizerdna::Sketch{}};
                }

                const std::string seq = graph_.get_oriented_sequence(v);
                if (seq.empty() || seq == "*") {
                    prog.hit();
                    return {vid, minimizerdna::Sketch{}};
                }

                const minimizerdna::MinimizerBuilder mzb(mm_opt_);
                auto result = std::make_pair(vid, mzb.build(seq));
                prog.hit();
                return result;
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

    log_stream() << "  - Source nodes: " << sources.size() << "\n";
    log_stream() << "  - Total valid node sketches: " << valid << "\n\n";

    return node_sketches;
}

std::vector<NodeGroupEntry> GfaBubbleFinder::cluster_sources_by_mm_(
    const std::vector<Vertex>& sources,
    const std::vector<minimizerdna::Sketch>& node_sketches
) const
{
    log_stream() << "Clustering sources by minimizer index ...\n";

    struct SourceRef {
        Vertex vtx{Vertex(UINT32_MAX)};
        uint32_t len{0};
        uint32_t comp{UINT32_MAX};
    };


    // Collect valid source nodes
    std::vector<SourceRef> refs;
    refs.reserve(sources.size());

    std::unordered_set<uint32_t> seen_vertices;
    seen_vertices.reserve(sources.size() * 2 + 1);

    for (Vertex src : sources) {
        const uint32_t vid = src.vertex_id();
        const uint32_t sid = src.segment_id();

        if (vid >= node_sketches.size()) continue;
        if (node_sketches[vid].empty()) continue;
        if (graph_.getNodeDeleted(sid)) continue;
        if (!seen_vertices.insert(vid).second) continue;

        const uint32_t len = graph_.getNodeLength(sid);

        const uint32_t comp = graph_.node_connectivity_id(src);
        if (comp == UINT32_MAX) continue;

        refs.push_back(SourceRef{src, len, comp});
    }


    // Sort sources by length
    std::sort(refs.begin(), refs.end(), [](const SourceRef& a, const SourceRef& b) {
            if (a.len != b.len) return a.len > b.len;
            return a.vtx.vertex_id() < b.vtx.vertex_id();
        }
    );

    const size_t n = refs.size();
    std::vector<NodeGroupEntry> groups;

    if (n < 2) {
        log_stream() << "  - Total source groups: 0\n\n";
        return groups;
    }


    // Build the minimizer inverted index
    std::unordered_map<uint64_t, std::vector<uint32_t>> postings;

    size_t total_minimizers = 0;
    for (const SourceRef& ref : refs) {
        total_minimizers += node_sketches[ref.vtx.vertex_id()].size();
    }

    postings.reserve(std::min<size_t>(total_minimizers, n * 128 + 1024));

    for (uint32_t i = 0; i < n; ++i) {
        const auto& sketch = node_sketches[refs[i].vtx.vertex_id()];

        for (uint64_t hash : sketch.vals) {
            postings[hash].push_back(i);
        }
    }


    // Buffers for clustering
    std::vector<uint32_t> shared_count(n, 0);
    std::vector<uint32_t> touched;
    touched.reserve(1024);

    std::vector<uint8_t> used(n, 0);

    ProgressTracker prog(n);

    groups.reserve(n / 2 + 1);


    // Cluster sources using greedy anchors
    for (uint32_t i = 0; i < n; ++i) {
        prog.hit();

        if (used[i]) continue;

        used[i] = 1;
        touched.clear();

        const SourceRef& a = refs[i];
        const auto& sketch_a = node_sketches[a.vtx.vertex_id()];

        if (params_.diff_sim <= 0.0) {
            for (uint32_t j = i + 1; j < n; ++j) {
                if (used[j]) continue;
                if (refs[j].comp == a.comp) continue;

                shared_count[j] = 0;
                touched.push_back(j);
            }
        } else {
            for (uint64_t hash : sketch_a.vals) {
                const auto it = postings.find(hash);
                if (it == postings.end()) continue;

                for (uint32_t j : it->second) {
                    if (j <= i || used[j]) continue;

                    if (refs[j].comp == a.comp) {
                        continue;
                    }

                    if (shared_count[j]++ == 0) {
                        touched.push_back(j);
                    }
                }
            }
        }

        // Keep the best hit from each connected component
        struct BestHit {
            uint32_t idx{UINT32_MAX};
            uint32_t shared{0};
            uint32_t len{0};
            double similarity{0.0};
        };

        std::unordered_map<uint32_t, BestHit> best_by_component;
        best_by_component.reserve(touched.size() + 1);

        for (uint32_t j : touched) {
            const uint32_t shared = shared_count[j];
            shared_count[j] = 0;

            if (used[j]) continue;

            const SourceRef& b = refs[j];
            const auto& sketch_b = node_sketches[b.vtx.vertex_id()];

            if (params_.diff_sim > 0.0) {
                const size_t min_size = std::min(sketch_a.size(), sketch_b.size());
                const size_t required_shared = static_cast<size_t>(std::ceil(params_.diff_sim * double(min_size)));
                if (shared < required_shared) continue;
            }

            const size_t min_size = std::min(sketch_a.size(), sketch_b.size());
            const double similarity = min_size == 0 ? 0.0 : static_cast<double>(shared) / static_cast<double>(min_size);

            if (similarity < params_.diff_sim) continue;

            BestHit hit;
            hit.idx = j;
            hit.shared = shared;
            hit.len = b.len;
            hit.similarity = similarity;

            auto it = best_by_component.find(b.comp);

            if (it == best_by_component.end()) {
                best_by_component.emplace(b.comp, hit);
                continue;
            }

            const BestHit& old = it->second;

            const bool better =
                hit.similarity > old.similarity ||
                (hit.similarity == old.similarity && hit.shared > old.shared) ||
                (hit.similarity == old.similarity && hit.shared == old.shared && hit.len > old.len) ||
                (hit.similarity == old.similarity && hit.shared == old.shared && hit.len == old.len && hit.idx < old.idx);

            if (better) it->second = hit;
        }

        // Rank
        if (best_by_component.empty()) continue;

        std::vector<BestHit> hits;
        hits.reserve(best_by_component.size());

        for (const auto& [comp, hit] : best_by_component) {
            (void)comp;
            hits.push_back(hit);
        }

        std::sort(
            hits.begin(), hits.end(), [](const BestHit& x, const BestHit& y) {
                if (x.similarity != y.similarity) {
                    return x.similarity > y.similarity;
                }
                if (x.shared != y.shared) {
                    return x.shared > y.shared;
                }
                if (x.len != y.len) {
                    return x.len > y.len;
                }
                return x.idx < y.idx;
            }
        );

        // Build source group
        NodeGroupEntry group;
        group.vtxs.reserve(hits.size() + 1);
        group.vtxs.push_back(a.vtx);

        double min_similarity = 1.0;

        for (const BestHit& hit : hits) {
            if (hit.idx == UINT32_MAX || used[hit.idx]) continue;

            used[hit.idx] = 1;
            group.vtxs.push_back(refs[hit.idx].vtx);

            min_similarity = std::min(min_similarity, hit.similarity);
        }

        if (group.vtxs.size() >= 2) {
            group.min_similarity = min_similarity;
            groups.emplace_back(std::move(group));
        }
    }

    log_stream() << "  - Total source groups: " << groups.size() << "\n\n";

    if (DEBUG_ENABLED) {
        for (const NodeGroupEntry& group : groups) {
            std::string group_str;
            group_str += "  - group sim=" + std::to_string(group.min_similarity) + ":";

            for (Vertex v : group.vtxs) {
                group_str += " " + graph_.getNodeName(v.segment_id()) + (v.is_reverse() ? "-" : "+");
            }

            debug_stream() << group_str << "\n";
        }

        debug_stream() << "\n";
    }

    return groups;
}

std::vector<NodeGroupEntry> GfaBubbleFinder::dedup_adjacent_source_groups_(
    std::vector<NodeGroupEntry> groups
) const
{
    log_stream() << "Deduplicating adjacent source groups ...\n";

    struct SourcePosition {
        uint32_t comp{UINT32_MAX};
        uint32_t rank{UINT32_MAX};
        uint32_t vertex{UINT32_MAX};
    };

    struct RankedGroup {
        size_t index{0};
        std::vector<SourcePosition> positions;
        bool valid{true};
    };

    std::vector<RankedGroup> ranked;
    ranked.reserve(groups.size());

    for (size_t i = 0; i < groups.size(); ++i) {
        RankedGroup item;
        item.index = i;
        item.positions.reserve(groups[i].vtxs.size());

        for (Vertex v : groups[i].vtxs) {
            SourcePosition pos{
                graph_.node_connectivity_id(v),
                graph_.vertex_topo_rank(v),
                v.vertex_id()
            };
            if (pos.comp == UINT32_MAX || pos.rank == UINT32_MAX) item.valid = false;
            item.positions.push_back(pos);
        }

        std::sort(item.positions.begin(), item.positions.end(),
            [](const SourcePosition& a, const SourcePosition& b) {
                if (a.comp != b.comp) return a.comp < b.comp;
                if (a.rank != b.rank) return a.rank < b.rank;
                return a.vertex < b.vertex;
            }
        );

        ranked.emplace_back(std::move(item));
    }

    std::sort(ranked.begin(), ranked.end(), [](const RankedGroup& a, const RankedGroup& b) {
        return std::lexicographical_compare(
            a.positions.begin(), a.positions.end(), b.positions.begin(), b.positions.end(),
            [](const SourcePosition& x, const SourcePosition& y) {
                if (x.comp != y.comp) return x.comp < y.comp;
                if (x.rank != y.rank) return x.rank < y.rank;
                return x.vertex < y.vertex;
            }
        );
    });

    std::vector<NodeGroupEntry> deduped;
    deduped.reserve(groups.size());
    const RankedGroup* previous = nullptr;

    for (const RankedGroup& current : ranked) {
        bool adjacent = previous && previous->valid && current.valid && previous->positions.size() == current.positions.size();

        if (adjacent) {
            for (size_t i = 0; i < current.positions.size(); ++i) {
                const SourcePosition& prev = previous->positions[i];
                const SourcePosition& curr = current.positions[i];

                if (curr.comp != prev.comp || curr.rank < prev.rank || curr.rank - prev.rank > 1) {
                    adjacent = false;
                    break;
                }
            }
        }

        if (!adjacent) deduped.emplace_back(std::move(groups[current.index]));
        previous = &current;
    }

    log_stream() << "  - Source groups after adjacent-topology deduplication: " << deduped.size() << " (removed " << groups.size() - deduped.size() << ")\n\n";

    return deduped;
}

std::vector<HomoTask_> GfaBubbleFinder::build_homo_tasks_(
    std::vector<NodeGroupEntry> groups
) const
{
    log_stream() << "Collecting homo tasks ...\n";

    std::vector<HomoTask_> tasks;
    tasks.reserve(sources_.size() + groups.size());

    auto has_sequence = [&](Vertex v) {
        const GfaNode* node = graph_.getNode(v.segment_id());
        return node && !node->deleted && !node->sequence.empty() && node->sequence != "*";
    };

    for (Vertex src : sources_) {
        std::unordered_set<uint32_t> usable_starts;
        for (const GfaArc* arc : graph_.getArcsFromVertex(src.vertex_id())) {
            if (!arc || arc->get_del() || (params_.skip_comp && arc->get_comp())) continue;
            const Vertex start(arc->get_target_vertex_id());
            if (start != src && has_sequence(start)) usable_starts.insert(start.segment_id());
            if (usable_starts.size() >= 2) break;
        }
        if (usable_starts.size() < 2) continue;

        tasks.push_back(HomoTask_{
            Type::SameSource,
            graph_.vertex_topo_rank(src),
            src.vertex_id(),
            {src}
        });
    }

    for (NodeGroupEntry& group : groups) {
        group.vtxs.erase(
            std::remove_if(group.vtxs.begin(), group.vtxs.end(), [&](Vertex v) { return !has_sequence(v); }),
            group.vtxs.end()
        );
        if (group.vtxs.size() < 2) continue;

        uint32_t best_rank = UINT32_MAX;
        uint32_t best_vid = UINT32_MAX;

        for (Vertex v : group.vtxs) {
            const uint32_t rank = graph_.vertex_topo_rank(v);
            const uint32_t vid = v.vertex_id();

            if (rank < best_rank || (rank == best_rank && vid < best_vid)) {
                best_rank = rank;
                best_vid = vid;
            }
        }

        tasks.push_back(HomoTask_{
            Type::DiffSource,
            best_rank,
            best_vid,
            std::move(group.vtxs)
        });
    }

    std::sort(tasks.begin(), tasks.end(), [](const HomoTask_& a, const HomoTask_& b) {
        if (a.topo_rank != b.topo_rank) return a.topo_rank < b.topo_rank;
        if (a.type != b.type) return a.type == Type::SameSource;
        return a.tie_id < b.tie_id;
    });

    if (DEBUG_ENABLED) {
        debug_stream() << "  - Homo tasks:\n";

        for (size_t i = 0; i < tasks.size(); ++i) {
            const HomoTask_& task = tasks[i];
            std::ostringstream oss;
            oss << "    task[" << i << "] type="
                << (task.type == Type::SameSource ? "Same" : "Diff")
                << " topo_rank=" << task.topo_rank
                << " tie_id=" << task.tie_id
                << " srcs:";

            for (Vertex v : task.srcs) {
                oss << " " << graph_.getNodeName(v.segment_id()) << (v.is_reverse() ? "-" : "+");
            }

            debug_stream() << oss.str() << "\n";
        }

        debug_stream() << "\n";
    }

    log_stream() << "  - Tasks: " << tasks.size() << "\n\n";
    return tasks;
}

void GfaBubbleFinder::sort_and_dedup_homologous_paths_(
    const HomologousParam& same_params,
    const HomologousParam& diff_params
)
{
    HomologousPathDeduplicator dedup(
        graph_,
        same_params,
        diff_params
    );

    dedup.run(homologous_paths_);
}

void GfaBubbleFinder::filter_diff_source_context_()
{
    struct Anchor {
        size_t path_id{0};
        uint64_t source_pair{0};
        std::array<uint32_t, 2> components{};
        std::array<Vertex, 2> sources{};
        std::array<Vertex, 2> starts{};
        std::array<Vertex, 2> ends{};
    };

    const uint32_t context_bp = params_.diff_min_len;
    if (context_bp == 0 || homologous_paths_.empty()) return;

    bool has_diff_source = false;
    for (const HomologousPath& hp : homologous_paths_) {
        if (hp.type == Type::DiffSource) {
            has_diff_source = true;
            break;
        }
    }
    if (!has_diff_source) return;

    const uint32_t vertex_count = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    std::vector<uint32_t> order;
    order.reserve(vertex_count);

    for (uint32_t v = 0; v < vertex_count; ++v) {
        if (!graph_.getNodeDeleted(Vertex(v).segment_id())) order.push_back(v);
    }

    std::sort(order.begin(), order.end(), [this](uint32_t a, uint32_t b) {
        const uint32_t rank_a = graph_.vertex_topo_rank(Vertex(a));
        const uint32_t rank_b = graph_.vertex_topo_rank(Vertex(b));
        return rank_a != rank_b ? rank_a < rank_b : a < b;
    });

    std::vector<uint32_t> suffix_bp(vertex_count, 0);
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const uint32_t v = *it;
        if (graph_.vertex_in_cycle(Vertex(v))) {
            suffix_bp[v] = context_bp;
            continue;
        }

        uint64_t best = 0;
        for (const GfaArc* arc : graph_.getArcsFromVertex(v)) {
            if (!arc || arc->get_del() || (params_.skip_comp && arc->get_comp())) continue;

            const uint32_t w = arc->get_target_vertex_id();
            if (w >= vertex_count) continue;

            const uint32_t length = graph_.getNodeLength(Vertex(w).segment_id());
            const uint32_t overlap = graph_.get_edge_ow(v, w);
            const uint64_t step = length - std::min(length, overlap);
            best = std::max(best, step + suffix_bp[w]);
            if (best >= context_bp) break;
        }
        suffix_bp[v] = static_cast<uint32_t>(std::min<uint64_t>(best, context_bp));
    }

    std::vector<Anchor> anchors;
    std::unordered_map<uint64_t, std::vector<size_t>> component_pairs;

    for (size_t path_id = 0; path_id < homologous_paths_.size(); ++path_id) {
        const HomologousPath& hp = homologous_paths_[path_id];
        if (hp.type != Type::DiffSource || hp.sources.size() < 2) continue;
        if (hp.starts.size() != hp.sources.size() || hp.ends.size() != hp.sources.size()) continue;

        for (size_t i = 0; i + 1 < hp.sources.size(); ++i) {
            for (size_t j = i + 1; j < hp.sources.size(); ++j) {
                const uint32_t component_i = graph_.node_connectivity_id(hp.sources[i]);
                const uint32_t component_j = graph_.node_connectivity_id(hp.sources[j]);
                if (component_i == UINT32_MAX || component_j == UINT32_MAX || component_i == component_j) continue;

                Anchor anchor;
                anchor.path_id = path_id;
                size_t first = i;
                size_t second = j;
                if (component_i > component_j) std::swap(first, second);

                anchor.components = {
                    graph_.node_connectivity_id(hp.sources[first]),
                    graph_.node_connectivity_id(hp.sources[second])
                };
                anchor.sources = {hp.sources[first], hp.sources[second]};
                anchor.starts = {hp.starts[first], hp.starts[second]};
                anchor.ends = {hp.ends[first], hp.ends[second]};

                const uint32_t source_a = std::min(
                    hp.sources[first].segment_id(), hp.sources[second].segment_id()
                );
                const uint32_t source_b = std::max(
                    hp.sources[first].segment_id(), hp.sources[second].segment_id()
                );
                anchor.source_pair = (uint64_t(source_a) << 32) | source_b;

                const uint64_t key = (uint64_t(anchor.components[0]) << 32) | anchor.components[1];
                component_pairs[key].push_back(anchors.size());
                anchors.push_back(anchor);
            }
        }
    }

    std::vector<uint8_t> rejected(homologous_paths_.size(), 0);
    size_t current_path = SIZE_MAX;
    std::unordered_map<uint64_t, bool> source_pair_supported;

    auto finish_path = [&]() {
        if (current_path == SIZE_MAX) return;
        for (const auto& pair : source_pair_supported) {
            if (!pair.second) {
                rejected[current_path] = 1;
                break;
            }
        }
        source_pair_supported.clear();
    };

    for (size_t anchor_id = 0; anchor_id < anchors.size(); ++anchor_id) {
        const Anchor& anchor = anchors[anchor_id];
        if (anchor.path_id != current_path) {
            finish_path();
            current_path = anchor.path_id;
        }

        const uint64_t key = (uint64_t(anchor.components[0]) << 32) | anchor.components[1];
        const std::vector<size_t>& peers = component_pairs.at(key);

        const uint32_t left_a = suffix_bp[Vertex::get_reverse(anchor.starts[0]).vertex_id()];
        const uint32_t left_b = suffix_bp[Vertex::get_reverse(anchor.starts[1]).vertex_id()];
        const uint32_t right_a = suffix_bp[anchor.ends[0].vertex_id()];
        const uint32_t right_b = suffix_bp[anchor.ends[1].vertex_id()];

        bool left_supported = left_a < context_bp || left_b < context_bp;
        bool right_supported = right_a < context_bp || right_b < context_bp;
        const bool relative_reverse = anchor.sources[0].is_reverse() != anchor.sources[1].is_reverse();

        auto check_peer = [&](size_t peer_id) {
            const Anchor& peer = anchors[peer_id];
            if ((peer.sources[0].is_reverse() != peer.sources[1].is_reverse()) != relative_reverse) return;

            std::array<Vertex, 2> peer_starts = peer.starts;
            std::array<Vertex, 2> peer_ends = peer.ends;
            if (peer.sources[0].is_reverse() != anchor.sources[0].is_reverse()) {
                for (size_t i = 0; i < 2; ++i) {
                    const Vertex start = Vertex::get_reverse(peer.ends[i]);
                    peer_ends[i] = Vertex::get_reverse(peer.starts[i]);
                    peer_starts[i] = start;
                }
            }

            const bool before =
                graph_.vertex_topo_rank(peer_ends[0]) < graph_.vertex_topo_rank(anchor.starts[0]) &&
                graph_.vertex_topo_rank(peer_ends[1]) < graph_.vertex_topo_rank(anchor.starts[1]);
            const bool after =
                graph_.vertex_topo_rank(peer_starts[0]) > graph_.vertex_topo_rank(anchor.ends[0]) &&
                graph_.vertex_topo_rank(peer_starts[1]) > graph_.vertex_topo_rank(anchor.ends[1]);

            left_supported = left_supported || before;
            right_supported = right_supported || after;
        };

        // Anchors are appended in path order. Skip this path as one range
        // instead of comparing every pair of its P-line occurrences.
        auto own_begin = std::lower_bound(
            peers.begin(), peers.end(), anchor.path_id,
            [&](size_t peer_id, size_t path_id) {
                return anchors[peer_id].path_id < path_id;
            }
        );
        auto own_end = std::upper_bound(
            own_begin, peers.end(), anchor.path_id,
            [&](size_t path_id, size_t peer_id) {
                return path_id < anchors[peer_id].path_id;
            }
        );

        for (auto it = peers.begin(); it != own_begin && !(left_supported && right_supported); ++it) {
            check_peer(*it);
        }
        for (auto it = own_end; it != peers.end() && !(left_supported && right_supported); ++it) {
            check_peer(*it);
        }

        // One logical source may occur on several P-lines. Keep the pair when
        // any occurrence combination has collinear context on both sides.
        source_pair_supported[anchor.source_pair] |= left_supported && right_supported;
    }
    finish_path();

    size_t next = 0;
    size_t removed = 0;
    for (size_t i = 0; i < homologous_paths_.size(); ++i) {
        if (rejected[i]) {
            const HomologousPath& path = homologous_paths_[i];
            warning_stream() << "  ! Possible misassembly or repetitive region " << ++removed << '\n';

            for (size_t j = 0; j < path.sources.size(); ++j) {
                std::ostringstream line;
                line << std::fixed << std::setprecision(2)
                     << "    ! path" << j << ": "
                     << graph_.getNodeName(path.sources[j].segment_id())
                     << (path.sources[j].is_reverse() ? '-' : '+')
                     << " -> "
                     << graph_.getNodeName(path.ends[j].segment_id())
                     << (path.ends[j].is_reverse() ? '-' : '+')
                     << " (" << static_cast<double>(path.lens[j]) / 1'000'000.0 << " Mb)";
                warning_stream() << line.str() << '\n';
            }
            continue;
        }

        if (next != i) homologous_paths_[next] = std::move(homologous_paths_[i]);
        ++next;
    }
    homologous_paths_.resize(next);

    log_stream() << "  - Internal different-source paths without collinear context removed: " << removed << '\n';
}

void GfaBubbleFinder::find_homologous_paths()
{
    bool has_sequence = false;
    for (uint32_t sid = 0; sid < graph_.getNumNodes(); ++sid) {
        const GfaNode* node = graph_.getNode(sid);
        if (node && !node->deleted && !node->sequence.empty() && node->sequence != "*") {
            has_sequence = true;
            break;
        }
    }

    if (!has_sequence) {
        homologous_paths_.clear();
        log_stream() << "Detecting homologous paths ...\n";
        log_stream() << "  - Skipped: input graph has no segment sequences\n";
        log_stream() << "  - Homologous paths before deduplication: 0\n";
        log_stream() << "  - Homologous paths after deduplication: 0\n";
        log_stream() << "  - Same source number: 0\n";
        log_stream() << "  - Different source number: 0\n\n";
        return;
    }

    HomologousParam homo_sameSource_Param(
        params_.diff_min_src,
        params_.same_sim,
        params_.same_min_len,
        Type::SameSource
    );
    HomologousParam homo_diffSource_Param(
        params_.diff_min_src,
        params_.diff_sim,
        params_.diff_min_len,
        Type::DiffSource
    );

    // ------------------------------------------------ Collect and cluster sources ------------------------------------------------
    std::vector<HomoTask_> tasks;
    {
        std::vector<Vertex> sources_for_homo = collect_sources_for_homo_();
        std::vector<minimizerdna::Sketch> node_sketches = build_source_mm_library_(sources_for_homo);
        std::vector<NodeGroupEntry> groups = cluster_sources_by_mm_(sources_for_homo, node_sketches);

        groups = dedup_adjacent_source_groups_(std::move(groups));
        tasks = build_homo_tasks_(std::move(groups));
    }

    // ------------------------------------------------ Detect homologous paths ------------------------------------------------
    log_stream() << "Detecting homologous paths ...\n";

    homologous_paths_.clear();
    homologous_paths_.reserve(tasks.size());

    ProgressTracker prog(tasks.size());

    ThreadPool pool(std::max(1u, params_.threads));

    std::vector<std::future<HomologousPath>> futs;
    futs.reserve(tasks.size());

    for (size_t task_id = 0; task_id < tasks.size(); ++task_id) {
        futs.emplace_back(
            pool.submit(
                [&, task_id]() {
                    const HomoTask_& task = tasks[task_id];
                    const HomologousParam& params = task.type == Type::SameSource ? homo_sameSource_Param : homo_diffSource_Param;

                    HomologousPath path = detect_homologous_paths_from_sources_(
                        task.srcs,
                        params,
                        bubble_branch_pairs_
                    );
                    prog.hit();
                    return path;
                }
            )
        );
    }

    for (auto& f : futs) {
        HomologousPath hp = f.get();

        const size_t path_count = hp.paths.size();
        if (path_count < 2) continue;
        if (hp.sources.size() != path_count || hp.starts.size() != path_count ||
            hp.ends.size() != path_count || hp.lens.size() != path_count) {
            continue;
        }

        homologous_paths_.emplace_back(std::move(hp));
    }

    pool.stop();

    // ------------------------------------------------ Sort and deduplicate ------------------------------------------------
    const size_t before_dedup = homologous_paths_.size();

    sort_and_dedup_homologous_paths_(
        homo_sameSource_Param,
        homo_diffSource_Param
    );

    filter_diff_source_context_();

    uint64_t same_source_num = 0;
    uint64_t diff_source_num = 0;

    for (const HomologousPath& hp : homologous_paths_) {
        if (hp.type == Type::SameSource) {
            ++same_source_num;
        } else if (hp.type == Type::DiffSource) {
            ++diff_source_num;
        }
    }

    log_stream() << "  - Homologous paths before deduplication: " << before_dedup << "\n";
    log_stream() << "  - Homologous paths after deduplication: " << homologous_paths_.size() << "\n";
    log_stream() << "  - Same source number: " << same_source_num << "\n";
    log_stream() << "  - Different source number: " << diff_source_num << "\n\n";

    label_path_clusters();
}


/* ================================================================================================================
 *                                         HOMOLOGOUS PATH ENUMERATION START
 * ================================================================================================================ */
HomologousPathEnumerator::HomologousPathEnumerator(
    const GfaBubbleFinder& finder,
    const std::vector<Vertex>& srcs,
    const HomologousParam& params,
    const BubbleBranchPairs_& bubble_pairs
)
    : finder_(finder),
      graph_(finder.graph_),
      log_(finder.graph_),
      srcs_(srcs),
      params_(params),
      bubble_pairs_(bubble_pairs)
{
    DFS_DEPTH_ = finder_.params_.max_depth * 10;
    DFS_GUARD_ = finder_.params_.DFS_guard * 10;
}

HomologousPath HomologousPathEnumerator::run()
{
    out_.type = params_.type;

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);
    if (srcs_.empty()) return out_;
    if (V == 0 || srcs_.front().vertex_id() >= V) return out_;

    log_.header(srcs_);

    if (!prepare_starts_()) {
        log_.return_reason("invalid starts size");
        return out_;
    }
    out_.path_guided = use_p_paths_;

    log_.starts(starts_);

    if (!initialize_states_()) {
        log_.return_reason("state initialization failed");
        return out_;
    }

    bootstrap_();

    bool any_pair_good = compute_pair_stats_(out_.min_hi_similarity, out_.max_hi_similarity);

    if (!any_pair_good) {
        log_.return_reason("no supported pairs");
        return HomologousPath();
    }

    const size_t N = starts_.size();

    std::vector<size_t> alive_idx;
    std::vector<size_t> active_idx;
    alive_idx.reserve(N);
    active_idx.reserve(N);

    std::vector<uint8_t> need_move(N, 0);
    std::vector<uint8_t> still_supported(N, 0);
    std::vector<uint8_t> can_move(N, 0);
    std::vector<uint8_t> next_supported(N, 0);

    std::vector<Vertex> next_vtxs(N, Vertex(UINT32_MAX));
    std::vector<uint64_t> next_bp(N, 0);

    std::vector<std::unique_ptr<PathSketch>> next_sketch(N);
    std::vector<std::unique_ptr<std::unordered_set<uint32_t>>> next_visited(N);
    std::vector<std::unique_ptr<std::unordered_set<uint32_t>>> next_visited_seg(N);
    std::vector<std::unique_ptr<std::vector<uint32_t>>> next_walks(N);
    std::vector<GfaPathWalker::Cursor> next_cursors(N);

    uint32_t dfs_step = 0;
    uint64_t dfs_states = 0;

    // ------------------------------------------------ DFS exploration ------------------------------------------------
    while (dfs_step < DFS_DEPTH_) {
        collect_alive_active_(alive_idx, active_idx);

        if (alive_idx.size() < 2) {
            log_.break_reason("alive < 2");
            break;
        }

        if (active_idx.empty()) {
            log_.break_reason("no movable");
            break;
        }

        if (++dfs_states > DFS_GUARD_) {
            break;
        }

        std::fill(need_move.begin(), need_move.end(), 0);
        std::fill(still_supported.begin(), still_supported.end(), 0);
        std::fill(can_move.begin(), can_move.end(), 0);
        std::fill(next_supported.begin(), next_supported.end(), 0);

        for (size_t i = 0; i < N; ++i) {
            next_vtxs[i] = Vertex(UINT32_MAX);
            next_bp[i] = 0;
            next_sketch[i].reset();
            next_visited[i].reset();
            next_visited_seg[i].reset();
            next_walks[i].reset();
            next_cursors[i] = GfaPathWalker::Cursor{};
        }

        if (!plan_moves_(alive_idx, need_move, still_supported)) {
            break;
        }

        if (all_active_paths_intersected_()) {
            log_.break_reason("all active paths intersected");
            break;
        }

        if (!extend_needed_paths_(
                active_idx,
                need_move,
                can_move,
                next_vtxs,
                next_bp,
                next_sketch,
                next_visited,
                next_visited_seg,
                next_walks,
                next_cursors
            )
        ) {
            break;
        }

        collect_alive_active_(alive_idx, active_idx);

        double next_min_hi = 1.0;
        double next_max_hi = 0.0;

        if (!evaluate_next_state_(
                alive_idx,
                can_move,
                next_sketch,
                next_supported,
                next_min_hi,
                next_max_hi
            )
        ) {
            break;
        }

        commit_next_state_(
            active_idx,
            can_move,
            next_vtxs,
            next_bp,
            next_sketch,
            next_visited,
            next_visited_seg,
            next_walks,
            next_cursors,
            next_min_hi,
            next_max_hi
        );

        // Paths that are not supported by any pair in the next state
        kill_unsupported_(alive_idx, next_supported);

        if (all_active_paths_intersected_()) {
            log_.break_reason("all active paths intersected");
            break;
        }

        ++dfs_step;
    }

    return emit_result_();
}

bool HomologousPathEnumerator::share_one_bubble_group_(const std::vector<Vertex>& vs) const
{
    std::vector<uint32_t> segs;
    segs.reserve(vs.size());

    for (const Vertex& v : vs) {
        segs.push_back(v.segment_id());
    }

    std::sort(segs.begin(), segs.end());
    segs.erase(std::unique(segs.begin(), segs.end()), segs.end());

    if (segs.size() < 2) return false;

    size_t base_i = SIZE_MAX;
    size_t base_n = SIZE_MAX;

    for (size_t i = 0; i < segs.size(); ++i) {
        std::unordered_map<uint32_t, std::vector<uint32_t>>::const_iterator it = bubble_pairs_.seg_to_groups.find(segs[i]);

        if (it == bubble_pairs_.seg_to_groups.end() || it->second.empty()) {
            return false;
        }

        if (it->second.size() < base_n) {
            base_i = i;
            base_n = it->second.size();
        }
    }

    const std::vector<uint32_t>& base_groups = bubble_pairs_.seg_to_groups.at(segs[base_i]);

    for (uint32_t gid : base_groups) {
        bool ok = true;

        for (size_t i = 0; i < segs.size(); ++i) {
            if (i == base_i) continue;

            const std::vector<uint32_t>& xs = bubble_pairs_.seg_to_groups.at(segs[i]);

            if (!std::binary_search(xs.begin(), xs.end(), gid)) {
                ok = false;
                break;
            }
        }

        if (ok) return true;
    }

    return false;
}

bool HomologousPathEnumerator::prepare_starts_()
{
    starts_.clear();
    owners_.clear();
    path_cursors_.clear();
    required_samples_.clear();
    use_p_paths_ = false;
    use_fixed_samples_ = false;

    std::vector<Vertex> raw_starts;
    std::vector<Vertex> raw_owners;
    raw_starts.reserve(std::max<size_t>(8, srcs_.size()));
    raw_owners.reserve(std::max<size_t>(8, srcs_.size()));

    const uint32_t V = static_cast<uint32_t>(graph_.getNumNodes() * 2);

    if (srcs_.size() == 1) {  // Same source
        const Vertex src = srcs_.front();

        if (graph_.getNodeDeleted(src.segment_id())) return false;

        std::unordered_set<uint32_t> seen_seg;
        seen_seg.reserve(8);

        for (const GfaArc* a : graph_.getArcsFromVertex(src.vertex_id())) {
            if (!a || a->get_del() || (finder_.params_.skip_comp && a->get_comp())) continue;

            Vertex w(a->get_target_vertex_id());
            if (w == src) continue;
            if (graph_.getNodeDeleted(w.segment_id())) continue;

            const uint32_t sid = w.segment_id();
            if (!seen_seg.insert(sid).second) continue;

            raw_starts.push_back(w);
            raw_owners.push_back(src);
        }

        if (raw_starts.size() >= 2 && share_one_bubble_group_(raw_starts)) {
            log_.return_reason("bubble branch group hit");
            return false;
        }
    } else {  // Different sources
        std::unordered_set<uint32_t> seen_seg;
        seen_seg.reserve(srcs_.size() * 2 + 8);

        for (const Vertex& s : srcs_) {
            if (s.vertex_id() >= V) continue;
            if (graph_.getNodeDeleted(s.segment_id())) continue;
            if (!seen_seg.insert(s.segment_id()).second) continue;

            raw_starts.push_back(s);
            raw_owners.push_back(s);
        }
    }

    if (raw_starts.size() < 2) return false;

    auto has_two_groups = [&]() {
        std::unordered_set<uint32_t> groups;
        groups.reserve(starts_.size());

        for (size_t i = 0; i < starts_.size(); ++i) {
            const Vertex anchor = params_.type == Type::SameSource ? starts_[i] : owners_[i];
            groups.insert(anchor.segment_id());
        }
        return groups.size() >= 2;
    };

    // With real P-lines, every state follows one exact oriented contig occurrence.
    if (!finder_.path_walker_.empty()) {
        for (size_t i = 0; i < raw_starts.size(); ++i) {
            for (const GfaPathWalker::Cursor& cursor : finder_.path_walker_.locate(raw_starts[i].vertex_id())) {
                if (params_.type == Type::SameSource &&
                    !finder_.path_walker_.preceded_by(cursor, raw_owners[i].vertex_id())) {
                    continue;
                }

                starts_.push_back(raw_starts[i]);
                owners_.push_back(raw_owners[i]);
                path_cursors_.push_back(cursor);
                required_samples_.push_back(UINT32_MAX);
            }
        }

        use_p_paths_ = true;
        if (has_two_groups()) return true;

        starts_.clear();
        owners_.clear();
        path_cursors_.clear();
        required_samples_.clear();
        use_p_paths_ = false;
    }

    // Without P-lines, duplicate each graph start by its fixed SN label.
    for (size_t i = 0; i < raw_starts.size(); ++i) {
        std::vector<uint32_t> samples;
        if (params_.type == Type::SameSource) {
            samples = common_samples(
                graph_, raw_owners[i].vertex_id(), raw_starts[i].vertex_id()
            );
        } else {
            samples = graph_.getNodeSampleIds(raw_starts[i].segment_id());
        }

        for (uint32_t sample : samples) {
            starts_.push_back(raw_starts[i]);
            owners_.push_back(raw_owners[i]);
            path_cursors_.emplace_back();
            required_samples_.push_back(sample);
        }
    }

    if (!starts_.empty()) {
        use_fixed_samples_ = true;
        if (has_two_groups()) return true;

        starts_.clear();
        owners_.clear();
        path_cursors_.clear();
        required_samples_.clear();
        use_fixed_samples_ = false;
    }

    // Legacy greedy walk is used only when the graph has no usable P or SN evidence.
    starts_ = std::move(raw_starts);
    owners_ = std::move(raw_owners);
    path_cursors_.resize(starts_.size());
    required_samples_.assign(starts_.size(), UINT32_MAX);
    return true;
}

bool HomologousPathEnumerator::initialize_states_()
{
    // ------------------------------------------------ Initialize path states ------------------------------------------------
    const size_t N = starts_.size();
    if (N < 2) return false;
    if (owners_.size() != N || path_cursors_.size() != N || required_samples_.size() != N) return false;

    cur_.assign(N, Vertex(UINT32_MAX));
    ends_.assign(N, Vertex(UINT32_MAX));
    acc_bp_.assign(N, 0);
    active_.assign(N, PATH_ACTIVE);
    keep_.assign(N, 0);

    visited_.clear();
    visited_seg_.clear();
    walks_.clear();

    visited_.resize(N);
    visited_seg_.resize(N);
    walks_.resize(N);

    // ------------------------------------------------ Bloom filter initialization ------------------------------------------------
    fastbloom::BloomOptions bloom_opt;
    bloom_opt.bits = finder_.params_.homo_bloom_bits;
    bloom_opt.hashes = finder_.params_.homo_bloom_hash;
    bloom_opt.use_simd = true;

    sketches_.clear();
    sketches_.reserve(N);

    intersection_scan_pos_.clear();
    paths_intersected_.clear();
    if (params_.type == Type::SameSource) {
        intersection_scan_pos_.assign(N, 0);
        paths_intersected_.assign(N * N, 0);
    }

    for (size_t i = 0; i < N; ++i) {
        const Vertex v = starts_[i];

        cur_[i] = v;
        ends_[i] = v;

        visited_[i].insert(v.vertex_id());
        visited_seg_[i].insert(owners_[i].segment_id());
        visited_seg_[i].insert(v.segment_id());
        walks_[i].push_back(v.vertex_id());

        acc_bp_[i] = graph_.getNodeLength(v.segment_id());

        sketches_.emplace_back(minimizerdna::Sketch{}, bloom_opt);
    }

    return true;
}

bool HomologousPathEnumerator::bootstrap_()
{
    for (size_t i = 0; i < starts_.size(); ++i) {
        if (active_[i] != PATH_ACTIVE) continue;

        extend_path_(i, true);
    }

    update_path_intersections_();

    return true;
}

bool HomologousPathEnumerator::extend_path_(size_t i, bool include_current)
{
    auto freeze = [&]() {
        active_[i] = PATH_FROZEN;
        ends_[i] = cur_[i];
    };

    auto add_to_sketch = [&](std::string_view seq) {
        if (seq.empty() || seq == "*") return;

        const minimizerdna::MinimizerBuilder mzb(finder_.mm_opt_);
        const minimizerdna::Sketch sketch = mzb.build(seq);
        if (!sketch.empty()) sketches_[i].add_sketch(sketch);
    };

    const uint32_t current_len = graph_.getNodeLength(cur_[i].segment_id());
    const std::string current_seq = graph_.get_oriented_sequence(cur_[i]);

    if (current_seq.empty() || current_seq == "*") {
        freeze();
        return include_current;
    }

    if (include_current && finder_.params_.homo_extend_bp > 0 && current_len >= finder_.params_.homo_extend_bp) {
        add_to_sketch(current_seq);
        ends_[i] = cur_[i];
        log_.sequence(std::to_string(current_seq.size()));
        return true;
    }

    std::unordered_set<uint32_t> region_set;

    const bool old_debug = DEBUG_ENABLED;

    if (old_debug) {
        DEBUG_ENABLED = false;
    }

    GfaPathWalker::Cursor next_cursor = path_cursors_[i];
    std::vector<std::vector<uint32_t>> paths;

    if (use_p_paths_) {
        std::vector<uint32_t> path = finder_.path_walker_.extend(
            next_cursor, finder_.params_.homo_extend_bp
        );
        if (!path.empty()) paths.push_back(std::move(path));
    } else {
        paths = graph_.open_walk(
            cur_[i].vertex_id(),
            region_set,
            DFS_DEPTH_,
            finder_.params_.skip_comp,
            DFS_GUARD_,
            finder_.params_.homo_extend_bp,
            &visited_seg_[i],
            required_samples_[i]
        );
    }

    if (old_debug) DEBUG_ENABLED = true;

    if (paths.empty() || paths.front().size() < 2) {
        if (include_current) {
            add_to_sketch(current_seq);
            log_.sequence(std::to_string(current_seq.size()));
        }

        freeze();
        return include_current;
    }

    const std::vector<uint32_t>& path = paths.front();
    const std::string full_seq = graph_.get_path_sequence(path);

    if (full_seq.empty() || full_seq == "*") {
        freeze();
        return include_current;
    }

    std::string_view added_seq(full_seq);
    if (!include_current) {
        if (added_seq.size() < current_len) {
            freeze();
            return false;
        }
        added_seq.remove_prefix(current_len);
    }

    add_to_sketch(added_seq);
    acc_bp_[i] += include_current ? full_seq.size() - current_len : added_seq.size();

    for (size_t p = 1; p < path.size(); ++p) {
        const Vertex v(path[p]);
        cur_[i] = v;
        ends_[i] = v;
        visited_[i].insert(v.vertex_id());
        visited_seg_[i].insert(v.segment_id());
        walks_[i].push_back(v.vertex_id());
    }
    if (use_p_paths_) path_cursors_[i] = next_cursor;

    log_.sequence(std::to_string(added_seq.size()));

    return true;
}

std::pair<double, double> HomologousPathEnumerator::pair_containments_(size_t i, size_t j) const
{
    return sketches_[i].bit_containments(sketches_[j]);
}

bool HomologousPathEnumerator::pair_eligible_(size_t i, size_t j) const
{
    if (params_.type == Type::SameSource) {
        return starts_[i].segment_id() != starts_[j].segment_id();
    }
    return owners_[i].segment_id() != owners_[j].segment_id();
}

bool HomologousPathEnumerator::pair_supported_(size_t i, size_t j, const std::pair<double, double>& c) const
{
    const bool first_good  = c.first  >= params_.min_similarity;
    const bool second_good = c.second >= params_.min_similarity;

    if (active_[i] == PATH_FROZEN && active_[j] == PATH_FROZEN) return first_good && second_good;
    if (active_[i] == PATH_FROZEN) return second_good;
    if (active_[j] == PATH_FROZEN) return first_good;

    return first_good || second_good;
}

bool HomologousPathEnumerator::pair_can_still_change_(
    size_t i,
    size_t j,
    const std::vector<uint8_t>& still_supported,
    const std::vector<uint8_t>& need_move
) const {
    if (!still_supported[i] || !still_supported[j]) return true;
    if (active_[i] == PATH_FROZEN && active_[j] == PATH_FROZEN) return false;
    if (active_[i] == PATH_FROZEN && active_[j] == PATH_ACTIVE) return !need_move[j];
    if (active_[j] == PATH_FROZEN && active_[i] == PATH_ACTIVE) return !need_move[i];
    if (active_[i] == PATH_ACTIVE && active_[j] == PATH_ACTIVE) return !need_move[i] || !need_move[j];

    return false;
}

bool HomologousPathEnumerator::compute_pair_stats_(
    double& min_hi,
    double& max_hi
) {
    bool any_pair_good = false;

    min_hi = 1.0;
    max_hi = 0.0;

    double min_first = 1.0;
    double max_first = 0.0;
    double min_second = 1.0;
    double max_second = 0.0;

    std::vector<size_t> alive;
    alive.reserve(starts_.size());

    for (size_t i = 0; i < starts_.size(); ++i) {
        if (active_[i] != PATH_DEAD) {
            alive.push_back(i);
        }
    }

    if (alive.size() < 2) {
        return any_pair_good;
    }

    for (size_t ai = 0; ai + 1 < alive.size(); ++ai) {
        const size_t i = alive[ai];

        for (size_t aj = ai + 1; aj < alive.size(); ++aj) {
            const size_t j = alive[aj];
            if (!pair_eligible_(i, j)) continue;

            const std::pair<double, double> c = pair_containments_(i, j);
            const double hi = std::max(c.first, c.second);

            if (c.first < min_first) min_first = c.first;
            if (c.first > max_first) max_first = c.first;

            if (c.second < min_second) min_second = c.second;
            if (c.second > max_second) max_second = c.second;

            if (hi < min_hi) min_hi = hi;
            if (hi > max_hi) max_hi = hi;

            if (hi >= params_.min_similarity) {
                keep_[i] = 1;
                keep_[j] = 1;
                any_pair_good = true;
            }
        }
    }

    log_.similarity(
        min_first,
        max_first,
        min_second,
        max_second,
        min_hi,
        max_hi,
        any_pair_good
    );

    return any_pair_good;
}

void HomologousPathEnumerator::collect_alive_active_(
    std::vector<size_t>& alive_idx,
    std::vector<size_t>& active_idx
) const {
    alive_idx.clear();
    active_idx.clear();

    for (size_t i = 0; i < starts_.size(); ++i) {
        if (active_[i] != PATH_DEAD) {
            alive_idx.push_back(i);
        }

        if (active_[i] == PATH_ACTIVE) {
            active_idx.push_back(i);
        }
    }
}

bool HomologousPathEnumerator::plan_moves_(
    const std::vector<size_t>& alive_idx,
    std::vector<uint8_t>& need_move,
    std::vector<uint8_t>& still_supported
) {
    bool any_need_move = false;

    // Similarity check for current state and determine which paths need to move
    for (size_t ai = 0; ai + 1 < alive_idx.size(); ++ai) {
        const size_t i = alive_idx[ai];

        for (size_t aj = ai + 1; aj < alive_idx.size(); ++aj) {
            const size_t j = alive_idx[aj];
            if (!pair_eligible_(i, j)) continue;

            if (!pair_can_still_change_(i, j, still_supported, need_move)) continue;

            const std::pair<double, double> c = pair_containments_(i, j);

            const bool first_good  = c.first  >= params_.min_similarity;
            const bool second_good = c.second >= params_.min_similarity;
            const bool supported = pair_supported_(i, j, c);

            if (supported) {
                still_supported[i] = 1;
                still_supported[j] = 1;
            }

            if (active_[i] == PATH_FROZEN && active_[j] == PATH_ACTIVE) {
                if (supported) {
                    need_move[j] = 1;
                    any_need_move = true;
                }

                continue;
            }

            if (active_[j] == PATH_FROZEN && active_[i] == PATH_ACTIVE) {
                if (supported) {
                    need_move[i] = 1;
                    any_need_move = true;
                }

                continue;
            }

            if (active_[i] == PATH_FROZEN && active_[j] == PATH_FROZEN) {
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
            active_[i] = PATH_DEAD;
        }
    }

    return any_need_move;
}

bool HomologousPathEnumerator::extend_needed_paths_(
    const std::vector<size_t>& active_idx,
    const std::vector<uint8_t>& need_move,
    std::vector<uint8_t>& can_move,
    std::vector<Vertex>& next_vtxs,
    std::vector<uint64_t>& next_bp,
    std::vector<std::unique_ptr<PathSketch>>& next_sketch,
    std::vector<std::unique_ptr<std::unordered_set<uint32_t>>>& next_visited,
    std::vector<std::unique_ptr<std::unordered_set<uint32_t>>>& next_visited_seg,
    std::vector<std::unique_ptr<std::vector<uint32_t>>>& next_walks,
    std::vector<GfaPathWalker::Cursor>& next_cursors
) {
    bool any_move = false;

    // Extend paths that need to move
    for (size_t i : active_idx) {
        if (active_[i] != PATH_ACTIVE || !need_move[i]) continue;

        Vertex old_cur = cur_[i];
        Vertex old_end = ends_[i];
        uint64_t old_bp = acc_bp_[i];
        const GfaPathWalker::Cursor old_cursor = path_cursors_[i];

        std::unordered_set<uint32_t> old_visited = std::move(visited_[i]);
        std::unordered_set<uint32_t> old_visited_seg = std::move(visited_seg_[i]);
        std::vector<uint32_t> old_walk = std::move(walks_[i]);
        PathSketch old_sketch = std::move(sketches_[i]);

        visited_[i] = old_visited;
        visited_seg_[i] = old_visited_seg;
        walks_[i] = old_walk;
        sketches_[i] = old_sketch;

        if (!extend_path_(i, false)) {
            cur_[i] = old_cur;
            ends_[i] = old_cur;
            acc_bp_[i] = old_bp;
            visited_[i] = std::move(old_visited);
            visited_seg_[i] = std::move(old_visited_seg);
            walks_[i] = std::move(old_walk);
            sketches_[i] = std::move(old_sketch);
            path_cursors_[i] = old_cursor;
            continue;
        }

        next_vtxs[i] = cur_[i];
        next_bp[i] = acc_bp_[i];

        next_sketch[i] = std::make_unique<PathSketch>(std::move(sketches_[i]));
        next_visited[i] = std::make_unique<std::unordered_set<uint32_t>>(std::move(visited_[i]));
        next_visited_seg[i] = std::make_unique<std::unordered_set<uint32_t>>(std::move(visited_seg_[i]));
        next_walks[i] = std::make_unique<std::vector<uint32_t>>(std::move(walks_[i]));
        next_cursors[i] = path_cursors_[i];

        cur_[i] = old_cur;
        ends_[i] = old_end;
        acc_bp_[i] = old_bp;
        visited_[i] = std::move(old_visited);
        visited_seg_[i] = std::move(old_visited_seg);
        walks_[i] = std::move(old_walk);
        sketches_[i] = std::move(old_sketch);
        path_cursors_[i] = old_cursor;

        can_move[i] = 1;
        any_move = true;
    }

    return any_move;
}

bool HomologousPathEnumerator::evaluate_next_state_(
    const std::vector<size_t>& alive_idx,
    const std::vector<uint8_t>& can_move,
    const std::vector<std::unique_ptr<PathSketch>>& next_sketch,
    std::vector<uint8_t>& next_supported,
    double& next_min_hi,
    double& next_max_hi
) const {
    next_min_hi = 1.0;
    next_max_hi = 0.0;

    // Similarity check for next state and determine which paths are supported in the next state
    for (size_t ai = 0; ai + 1 < alive_idx.size(); ++ai) {
        const size_t i = alive_idx[ai];

        for (size_t aj = ai + 1; aj < alive_idx.size(); ++aj) {
            const size_t j = alive_idx[aj];
            if (!pair_eligible_(i, j)) continue;

            const PathSketch& si = can_move[i] ? *next_sketch[i] : sketches_[i];
            const PathSketch& sj = can_move[j] ? *next_sketch[j] : sketches_[j];
            const std::pair<double, double> c = si.bit_containments(sj);

            const double hi = std::max(c.first, c.second);

            if (pair_supported_(i, j, c)) {
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
        if (!next_supported[i]) continue;

        ++next_alive_supported_cnt;

        if (active_[i] == PATH_ACTIVE) {
            ++next_movable_supported_cnt;
        }
    }

    if (next_alive_supported_cnt < 2) return false;
    if (next_movable_supported_cnt == 0) return false;

    return true;
}

void HomologousPathEnumerator::commit_next_state_(
    const std::vector<size_t>& active_idx,
    const std::vector<uint8_t>& can_move,
    const std::vector<Vertex>& next_vtxs,
    const std::vector<uint64_t>& next_bp,
    std::vector<std::unique_ptr<PathSketch>>& next_sketch,
    std::vector<std::unique_ptr<std::unordered_set<uint32_t>>>& next_visited,
    std::vector<std::unique_ptr<std::unordered_set<uint32_t>>>& next_visited_seg,
    std::vector<std::unique_ptr<std::vector<uint32_t>>>& next_walks,
    const std::vector<GfaPathWalker::Cursor>& next_cursors,
    double next_min_hi,
    double next_max_hi
) {
    // Update states for next step
    out_.min_hi_similarity = next_min_hi;
    out_.max_hi_similarity = next_max_hi;

    for (size_t i : active_idx) {
        if (!can_move[i]) continue;

        cur_[i] = next_vtxs[i];
        ends_[i] = next_vtxs[i];

        visited_[i] = std::move(*next_visited[i]);
        visited_seg_[i] = std::move(*next_visited_seg[i]);
        walks_[i] = std::move(*next_walks[i]);

        sketches_[i] = std::move(*next_sketch[i]);
        acc_bp_[i] = next_bp[i];
        path_cursors_[i] = next_cursors[i];
    }

    update_path_intersections_();
}

void HomologousPathEnumerator::kill_unsupported_(
    const std::vector<size_t>& alive_idx,
    const std::vector<uint8_t>& supported
) {
    for (size_t i : alive_idx) {
        if (active_[i] != PATH_DEAD && !supported[i]) {
            active_[i] = PATH_DEAD;
        }
    }
}

void HomologousPathEnumerator::update_path_intersections_()
{
    const size_t N = walks_.size();
    if (intersection_scan_pos_.size() != N || paths_intersected_.size() != N * N) return;

    // Index only newly committed nodes. A path pair stays intersected once it has shared any segment.
    for (size_t i = 0; i < N; ++i) {
        size_t& pos = intersection_scan_pos_[i];

        for (; pos < walks_[i].size(); ++pos) {
            const uint32_t segment = Vertex(walks_[i][pos]).segment_id();

            for (size_t j = 0; j < N; ++j) {
                if (i == j || visited_seg_[j].find(segment) == visited_seg_[j].end()) continue;
                paths_intersected_[i * N + j] = 1;
                paths_intersected_[j * N + i] = 1;
            }
        }
    }
}

bool HomologousPathEnumerator::all_active_paths_intersected_() const
{
    if (params_.type != Type::SameSource) return false;

    const size_t N = starts_.size();
    if (paths_intersected_.size() != N * N) return false;

    size_t pair_count = 0;

    // Frozen paths cannot reach a future intersection. Keep them for similarity,
    // but let only movable paths decide when extension can stop.
    for (size_t i = 0; i + 1 < N; ++i) {
        if (active_[i] != PATH_ACTIVE) continue;

        for (size_t j = i + 1; j < N; ++j) {
            if (active_[j] != PATH_ACTIVE || !pair_eligible_(i, j)) continue;
            ++pair_count;

            if (!paths_intersected_[i * N + j]) return false;
        }
    }

    return pair_count > 0;
}

HomologousPath HomologousPathEnumerator::emit_result_()
{
    out_.sources.clear();
    out_.starts.clear();
    out_.ends.clear();
    out_.lens.clear();
    out_.paths.clear();

    log_.kept(starts_, ends_, acc_bp_, keep_);
    std::set<std::vector<Vertex>> seen_paths;

    auto path_bp = [&](const std::vector<uint32_t>& path) {
        uint64_t bp = 0;

        for (size_t p = 0; p < path.size(); ++p) {
            const uint32_t segment = Vertex::get_segment_id(path[p]);
            if (graph_.getNodeDeleted(segment)) return uint64_t(0);

            const uint32_t length = graph_.getNodeLength(segment);
            if (p == 0) {
                bp += length;
            } else {
                const uint32_t overlap = graph_.get_edge_ow(path[p - 1], path[p]);
                bp += length - std::min(length, overlap);
            }
        }
        return bp;
    };

    auto append_path = [&](size_t state, const std::vector<uint32_t>& path) {
        if (path.empty()) return;

        const uint64_t bp = path_bp(path);
        if (bp < params_.min_path_bp) return;

        std::vector<Vertex> vertices;
        vertices.reserve(path.size() + (out_.type == Type::SameSource));
        if (out_.type == Type::SameSource) vertices.push_back(owners_[state]);
        for (uint32_t vertex : path) vertices.emplace_back(vertex);

        if (!seen_paths.insert(vertices).second) return;

        out_.sources.push_back(owners_[state]);
        out_.starts.emplace_back(path.front());
        out_.ends.emplace_back(path.back());
        out_.lens.push_back(bp);
        out_.paths.push_back(std::move(vertices));
    };

    for (size_t i = 0; i < starts_.size(); ++i) {
        if (!keep_[i]) continue;

        if (use_p_paths_) {
            append_path(i, walks_[i]);
            log_.enumerate(starts_[i], ends_[i], {walks_[i]});
            continue;
        }

        static const std::unordered_set<uint32_t> empty_region;

        const std::unordered_set<uint32_t>& visited_set = params_.type == Type::DiffSource ?
            visited_[i] : (finder_.params_.homo_num == 1 ? visited_[i] : empty_region);

        uint32_t eff_max_depth = DFS_DEPTH_;

        if (!visited_set.empty()) {
            const uint32_t limit = static_cast<uint32_t>(visited_set.size()) + finder_.DEPTH_MARGIN_;

            if (eff_max_depth == 0 || eff_max_depth > limit) {
                eff_max_depth = limit;
            }
        }

        bool hit_limits = false;

        std::vector<std::vector<uint32_t>> paths = graph_.enumerate_paths_greedy_DFS(
            starts_[i].vertex_id(),
            ends_[i].vertex_id(),
            visited_set,
            eff_max_depth,
            finder_.params_.homo_num,
            finder_.params_.skip_comp,
            hit_limits,
            DFS_GUARD_,
            finder_.params_.stall_round_limit,
            required_samples_[i]
        );

        for (const std::vector<uint32_t>& path : paths) append_path(i, path);
        log_.enumerate(starts_[i], ends_[i], paths);
    }

    if (DEBUG_ENABLED) {
        debug_stream() << "\n";
    }

    return out_;
}

HomologousPath GfaBubbleFinder::detect_homologous_paths_from_sources_(
    const std::vector<Vertex>& srcs,
    const HomologousParam& params,
    const GfaBubble::BubbleBranchPairs_& bubble_pairs
) const {
    HomologousPathEnumerator enumerator(
        *this,
        srcs,
        params,
        bubble_pairs
    );

    return enumerator.run();
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
/* ================================================================================================================
 *                                          HOMOLOGOUS PATH ENUMERATION END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                         HOMOLOGOUS PATH DEDUPLICATION START
 * ================================================================================================================ */
HomologousPathDeduplicator::HomologousPathDeduplicator(
    const GfaGraph& graph,
    const HomologousParam& same_params,
    const HomologousParam& diff_params
)
    : graph_(graph),
      same_params_(same_params),
      diff_params_(diff_params)
{}

void HomologousPathDeduplicator::run(std::vector<HomologousPath>& paths) const
{
    if (paths.empty()) return;

    std::vector<Candidate> candidates;
    candidates.reserve(paths.size());

    for (HomologousPath& hp : paths) {
        candidates.push_back(Candidate{
            group_size_(hp),
            tie_vertex_(hp),
            std::move(hp)
        });
    }

    std::sort(candidates.begin(), candidates.end(), candidate_better_);

    paths.clear();
    paths.reserve(candidates.size());

    std::unordered_map<uint32_t, std::vector<NodeHit>> node_index;
    node_index.reserve(candidates.size() * 8 + 1);

    for (Candidate& candidate : candidates) {
        std::vector<HomologousPath> pending;
        pending.push_back(std::move(candidate.hp));

        while (!pending.empty()) {
            HomologousPath hp = std::move(pending.back());
            pending.pop_back();

            TrimPlan plan;
            if (!find_trim_plan_(hp, node_index, plan)) {
                const uint32_t group_id = static_cast<uint32_t>(paths.size());
                paths.push_back(std::move(hp));
                index_group_(paths.back(), group_id, node_index);
                continue;
            }

            HomologousPath prefix;
            if (build_trimmed_piece_(hp, plan, true, prefix)) {
                pending.push_back(std::move(prefix));
            }

            HomologousPath suffix;
            if (build_trimmed_piece_(hp, plan, false, suffix)) {
                pending.push_back(std::move(suffix));
            }
        }
    }
}

uint64_t HomologousPathDeduplicator::group_size_(const HomologousPath& hp) const
{
    std::vector<uint32_t> segs;

    for (const auto& path : hp.paths) {
        segs.reserve(segs.size() + path.size());
        for (Vertex v : path) {
            segs.push_back(v.segment_id());
        }
    }

    std::sort(segs.begin(), segs.end());
    segs.erase(std::unique(segs.begin(), segs.end()), segs.end());

    uint64_t total = 0;
    for (uint32_t sid : segs) {
        if (!graph_.getNodeDeleted(sid)) {
            total += graph_.getNodeLength(sid);
        }
    }

    return total;
}

uint32_t HomologousPathDeduplicator::tie_vertex_(const HomologousPath& hp) const
{
    uint32_t tie_vertex = UINT32_MAX;

    for (const auto& path : hp.paths) {
        for (Vertex v : path) {
            tie_vertex = std::min(tie_vertex, v.vertex_id());
        }
    }

    return tie_vertex;
}

uint64_t HomologousPathDeduplicator::path_bp_(const std::vector<Vertex>& path) const
{
    if (path.empty()) return 0;

    uint64_t total = 0;

    for (size_t i = 0; i < path.size(); ++i) {
        const uint32_t sid = path[i].segment_id();
        if (graph_.getNodeDeleted(sid)) return 0;

        const uint32_t len = graph_.getNodeLength(sid);

        if (i == 0) {
            total += len;
            continue;
        }

        const uint32_t ow = graph_.get_edge_ow(path[i - 1], path[i]);
        total += len > ow ? uint64_t(len - ow) : 0;
    }

    return total;
}

const HomologousParam& HomologousPathDeduplicator::params_for_(Type type) const
{
    return type == Type::SameSource ? same_params_ : diff_params_;
}

bool HomologousPathDeduplicator::candidate_better_(
    const Candidate& a,
    const Candidate& b
) {
    if (a.size != b.size) return a.size > b.size;

    if (a.hp.paths.size() != b.hp.paths.size()) {
        return a.hp.paths.size() > b.hp.paths.size();
    }

    if (a.hp.max_hi_similarity != b.hp.max_hi_similarity) {
        return a.hp.max_hi_similarity > b.hp.max_hi_similarity;
    }

    return a.tie_vertex < b.tie_vertex;
}

std::vector<HomologousPathDeduplicator::CoveredRun> HomologousPathDeduplicator::covered_runs_for_path_(
    const std::vector<Vertex>& path,
    uint32_t path_id,
    const std::unordered_map<uint32_t, std::vector<NodeHit>>& node_index
) const {
    std::vector<CoveredRun> out;
    if (path.size() < 2) return out;

    std::unordered_map<uint64_t, std::array<RunState, 2>> states;
    states.reserve(path.size() * 2 + 1);

    for (uint32_t pos = 0; pos < path.size(); ++pos) {
        const uint32_t sid = path[pos].segment_id();

        auto it = node_index.find(sid);
        if (it == node_index.end()) continue;

        for (const NodeHit& hit : it->second) {
            const uint64_t key = (uint64_t(hit.group_id) << 32) | uint64_t(hit.path_id);
            std::array<RunState, 2>& runs = states[key];

            for (uint32_t direction = 0; direction < runs.size(); ++direction) {
                RunState& st = runs[direction];
                const bool forward = direction == 0;
                const bool contiguous = st.last_pos != UINT32_MAX && st.last_pos + 1 == pos && (forward ? st.last_old_pos + 1 == hit.pos : st.last_old_pos > 0 && st.last_old_pos - 1 == hit.pos);

                if (contiguous) {
                    ++st.cur_len;
                } else {
                    st.cur_beg = pos;
                    st.cur_len = 1;
                }

                st.last_pos = pos;
                st.last_old_pos = hit.pos;

                if (st.cur_len > st.best_len) {
                    st.best_len = st.cur_len;
                    st.best_beg = st.cur_beg;
                }
            }
        }
    }

    out.reserve(states.size());

    for (const auto& kv : states) {
        const uint32_t group_id = static_cast<uint32_t>(kv.first >> 32);
        for (const RunState& st : kv.second) {
            if (st.best_len >= 2) {
                out.push_back(CoveredRun{
                    group_id,
                    path_id,
                    st.best_beg,
                    st.best_len
                });
            }
        }
    }

    std::sort(out.begin(), out.end(),
        [](const CoveredRun& a, const CoveredRun& b) {
            if (a.group_id != b.group_id) return a.group_id < b.group_id;
            if (a.len != b.len) return a.len > b.len;
            return a.beg < b.beg;
        }
    );

    out.erase(
        std::unique(out.begin(), out.end(),
            [](const CoveredRun& a, const CoveredRun& b) {
                return a.group_id == b.group_id;
            }
        ),
        out.end()
    );

    return out;
}

bool HomologousPathDeduplicator::find_trim_plan_(
    const HomologousPath& hp,
    const std::unordered_map<uint32_t, std::vector<NodeHit>>& node_index,
    TrimPlan& plan
) const {
    if (hp.paths.empty()) return false;

    std::unordered_map<uint32_t, std::vector<CoveredRun>> by_group;
    by_group.reserve(hp.paths.size() * 4 + 1);

    for (size_t i = 0; i < hp.paths.size(); ++i) {
        std::vector<CoveredRun> runs = covered_runs_for_path_(
            hp.paths[i],
            static_cast<uint32_t>(i),
            node_index
        );

        if (runs.empty()) return false;

        for (const CoveredRun& run : runs) {
            by_group[run.group_id].push_back(run);
        }
    }

    uint32_t best_group = UINT32_MAX;
    uint64_t best_clip = 0;
    std::vector<CoveredRun> best_runs;

    for (const auto& kv : by_group) {
        const std::vector<CoveredRun>& runs = kv.second;

        std::vector<CoveredRun> per_path(hp.paths.size());
        std::vector<uint8_t> seen(hp.paths.size(), 0);

        for (const CoveredRun& run : runs) {
            if (run.path_id >= per_path.size()) continue;

            CoveredRun& old = per_path[run.path_id];
            const bool better = !seen[run.path_id] || run.len > old.len || (run.len == old.len && run.beg < old.beg);

            if (better) {
                per_path[run.path_id] = run;
                seen[run.path_id] = 1;
            }
        }

        bool all_paths_covered = true;
        for (uint8_t x : seen) {
            if (!x) {
                all_paths_covered = false;
                break;
            }
        }

        if (!all_paths_covered) continue;

        uint64_t clip = 0;
        for (const CoveredRun& run : per_path) {
            clip += run.len;
        }

        if (clip > best_clip || (clip == best_clip && kv.first < best_group)) {
            best_group = kv.first;
            best_clip = clip;
            best_runs = std::move(per_path);
        }
    }

    if (best_group == UINT32_MAX) return false;

    plan.group_id = best_group;
    plan.runs = std::move(best_runs);

    return true;
}

bool HomologousPathDeduplicator::build_trimmed_piece_(
    const HomologousPath& hp,
    const TrimPlan& plan,
    bool keep_prefix,
    HomologousPath& out
) const {
    if (!plan.valid()) return false;
    if (plan.runs.size() != hp.paths.size()) return false;

    out = HomologousPath();
    out.type = hp.type;
    out.min_hi_similarity = hp.min_hi_similarity;
    out.max_hi_similarity = hp.max_hi_similarity;
    out.path_guided = hp.path_guided;

    for (size_t i = 0; i < hp.paths.size(); ++i) {
        const std::vector<Vertex>& old_path = hp.paths[i];
        const CoveredRun& run = plan.runs[i];

        if (run.beg > old_path.size()) return false;
        if (run.end() > old_path.size()) return false;

        std::vector<Vertex> new_path;

        if (keep_prefix) {
            if (run.beg == 0) return false;
            new_path.assign(old_path.begin(), old_path.begin() + run.beg);
        } else {
            if (run.end() >= old_path.size()) return false;
            new_path.assign(old_path.begin() + run.end(), old_path.end());
        }

        if (!append_trimmed_path_(hp, i, new_path, out)) {
            return false;
        }
    }

    return passes_thresholds_(out);
}

bool HomologousPathDeduplicator::append_trimmed_path_(
    const HomologousPath& old_hp,
    size_t i,
    const std::vector<Vertex>& path,
    HomologousPath& out
) const {
    if (path.empty()) return false;

    const HomologousParam& params = params_for_(old_hp.type);
    Vertex source = path.front();
    Vertex start = path.front();
    std::vector<Vertex> len_path = path;

    if (old_hp.type == Type::SameSource && i < old_hp.sources.size()) {
        source = old_hp.sources[i];

        if (path.front() == source) {
            if (path.size() < 2) return false;
            start = path[1];
            len_path.assign(path.begin() + 1, path.end());
        } else {
            start = path.front();
        }
    }

    const uint64_t bp = path_bp_(len_path);

    if (bp < params.min_path_bp) return false;

    out.sources.push_back(source);
    out.starts.push_back(start);
    out.ends.push_back(path.back());
    out.lens.push_back(bp);
    out.paths.push_back(path);

    return true;
}

bool HomologousPathDeduplicator::passes_thresholds_(
    const HomologousPath& hp
) const {
    if (hp.paths.size() < 2) return false;
    if (hp.sources.size() != hp.paths.size()) return false;
    if (hp.starts.size() != hp.paths.size()) return false;
    if (hp.ends.size() != hp.paths.size()) return false;
    if (hp.lens.size() != hp.paths.size()) return false;

    const HomologousParam& params = params_for_(hp.type);

    for (size_t i = 0; i < hp.paths.size(); ++i) {
        if (hp.lens[i] < params.min_path_bp) return false;
    }

    return true;
}

void HomologousPathDeduplicator::index_group_(
    const HomologousPath& hp,
    uint32_t group_id,
    std::unordered_map<uint32_t, std::vector<NodeHit>>& node_index
) const {
    for (uint32_t path_id = 0; path_id < hp.paths.size(); ++path_id) {
        const std::vector<Vertex>& path = hp.paths[path_id];

        for (uint32_t pos = 0; pos < path.size(); ++pos) {
            const uint32_t sid = path[pos].segment_id();
            node_index[sid].push_back(NodeHit{group_id, path_id, pos});
        }
    }
}
/* ================================================================================================================
 *                                          HOMOLOGOUS PATH DEDUPLICATION END
 * ================================================================================================================ */


/* ================================================================================================================
 *                                         PATH CLUSTERING START
 * ================================================================================================================ */
PathClusterer::PathClusterer(
    const GfaGraph& graph,
    double min_similarity,
    uint32_t kmer,
    uint32_t window
)
    : graph_(graph),
      min_similarity_(min_similarity),
      kmer_(kmer),
      window_(window)
{}

std::vector<uint32_t> PathClusterer::cluster(
    const std::vector<std::vector<uint32_t>>& paths,
    Type type
) const
{
    std::vector<uint32_t> labels(paths.size(), UINT32_MAX);
    if (paths.empty()) return labels;

    std::vector<std::string> sequences(paths.size());
    std::vector<uint64_t> lengths(paths.size(), 0);
    std::vector<uint8_t> can_cluster(paths.size(), 0);
    const uint64_t min_cluster_bp = uint64_t(std::max(kmer_, window_)) * 10;

    for (size_t i = 0; i < paths.size(); ++i) {
        const std::vector<uint32_t> trimmed_path = GfaBubbleFinder::trim_path_cluster_anchors(paths[i], type);
        sequences[i] = graph_.get_path_sequence(trimmed_path);
        lengths[i] = sequences[i].size();
        can_cluster[i] = lengths[i] >= min_cluster_bp;
    }

    minimizerdna::Options options;
    options.k = kmer_;
    options.w = window_;
    options.canonical = false;

    const minimizerdna::MinimizerBuilder builder(options);
    std::vector<minimizerdna::Sketch> sketches(paths.size());

    for (size_t i = 0; i < paths.size(); ++i) {
        if (can_cluster[i]) {
            sketches[i] = builder.build(sequences[i]);
        }
    }

    std::vector<size_t> order(paths.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (lengths[a] != lengths[b]) return lengths[a] > lengths[b];
        return a < b;
    });

    std::vector<size_t> representatives;
    representatives.reserve(paths.size());
    for (size_t i : order) {
        if (!can_cluster[i]) {
            labels[i] = static_cast<uint32_t>(representatives.size());
            representatives.push_back(i);
            continue;
        }

        for (size_t cluster = 0; cluster < representatives.size(); ++cluster) {
            const size_t rep = representatives[cluster];
            if (!can_cluster[rep]) continue;

            if (sketches[i].jaccard(sketches[rep]) >= min_similarity_) {
                labels[i] = static_cast<uint32_t>(cluster);
                break;
            }
        }

        if (labels[i] == UINT32_MAX) {
            labels[i] = static_cast<uint32_t>(representatives.size());
            representatives.push_back(i);
        }
    }

    return labels;
}

void GfaBubbleFinder::label_path_clusters()
{
    const size_t total = bubbles_.size() + homologous_paths_.size();

    log_stream() << "Clustering bubble and homologous paths ...\n";
    log_stream() << "  - Bubble path groups: " << bubbles_.size() << "\n";
    log_stream() << "  - Homologous path groups: " << homologous_paths_.size() << "\n";

    if (total == 0) {
        log_stream() << "Finished clustering bubble and homologous paths\n\n";
        return;
    }

    PathClusterer clusterer(graph_, params_.path_sim, mm_opt_.k, mm_opt_.w);
    ThreadPool pool(std::max<uint32_t>(1, params_.threads));
    ProgressTracker prog(total);

    std::vector<std::future<std::vector<uint32_t>>> bubble_futs;
    bubble_futs.reserve(bubbles_.size());

    for (Bubble& bubble : bubbles_) {
        bubble_futs.push_back(pool.submit([&clusterer, &bubble, &prog]() {
            std::vector<uint32_t> clusters = clusterer.cluster(bubble.get_paths(), Type::Normal);
            prog.hit();
            return clusters;
        }));
    }

    std::vector<std::future<std::vector<uint32_t>>> homo_futs;
    homo_futs.reserve(homologous_paths_.size());

    for (HomologousPath& hp : homologous_paths_) {
        homo_futs.push_back(pool.submit([&clusterer, &hp, &prog]() {
            if (hp.path_guided && hp.type == Type::DiffSource) {
                std::unordered_set<uint32_t> owners;
                bool duplicate_owner = false;

                for (Vertex source : hp.sources) {
                    if (!owners.insert(source.vertex_id()).second) {
                        duplicate_owner = true;
                        break;
                    }
                }

                // Every occurrence must remain a representative when one source has multiple contig paths.
                if (duplicate_owner) {
                    std::vector<uint32_t> clusters(hp.paths.size());
                    std::iota(clusters.begin(), clusters.end(), 0);
                    prog.hit();
                    return clusters;
                }
            }

            std::vector<std::vector<uint32_t>> paths;
            paths.reserve(hp.paths.size());

            for (const std::vector<Vertex>& path : hp.paths) {
                std::vector<uint32_t> vertices;
                vertices.reserve(path.size());

                for (Vertex v : path) {
                    vertices.push_back(v.vertex_id());
                }

                paths.push_back(std::move(vertices));
            }

            std::vector<uint32_t> clusters = clusterer.cluster(paths, hp.type);
            prog.hit();
            return clusters;
        }));
    }

    for (size_t i = 0; i < bubble_futs.size(); ++i) {
        bubbles_[i].set_path_clusters(bubble_futs[i].get());
    }

    for (size_t i = 0; i < homo_futs.size(); ++i) {
        homologous_paths_[i].path_clusters = homo_futs[i].get();
    }

    prog.finish();
    pool.stop();

    log_stream() << "Finished clustering bubble and homologous paths\n\n";
}
/* ================================================================================================================
 *                                         PATH CLUSTERING END
 * ================================================================================================================ */


std::vector<uint32_t> GfaBubbleFinder::collect_orientation_conflict_segments() const
{
    const uint32_t node_count = static_cast<uint32_t>(graph_.getNumNodes());
    constexpr uint8_t ORIENTATION_MASK = 3u;
    constexpr uint8_t CONFLICT = 4u;
    constexpr uint8_t BOUNDARY = 8u;
    constexpr uint8_t REPORTED = 16u;

    std::vector<uint8_t> state(node_count, 0);
    std::vector<uint32_t> first_pos(node_count, UINT32_MAX);

    std::vector<uint32_t> out;
    std::vector<uint32_t> touched;
    std::vector<uint32_t> boundaries;
    std::vector<uint32_t> path_segments;
    std::vector<int32_t> repeat_delta;

    out.reserve(1024);
    touched.reserve(256);
    boundaries.reserve(16);
    path_segments.reserve(256);

    auto mark_boundary = [&](uint32_t seg_id) {
        if (seg_id >= node_count || (state[seg_id] & BOUNDARY)) return;
        state[seg_id] |= BOUNDARY;
        boundaries.push_back(seg_id);
    };

    auto clear_boundaries = [&]() {
        for (uint32_t seg_id : boundaries) state[seg_id] &= ~BOUNDARY;
        boundaries.clear();
    };

    auto add_segment = [&](uint32_t seg_id) {
        if (seg_id >= node_count) return;
        if (graph_.getNodeDeleted(seg_id)) return;
        if (state[seg_id] & REPORTED) return;

        state[seg_id] |= REPORTED;
        out.push_back(seg_id);
    };

    auto vertex_id_of = [](const auto& item) -> uint32_t {
        using Item = std::decay_t<decltype(item)>;
        if constexpr (std::is_same_v<Item, Vertex>) return item.vertex_id();
        else return static_cast<uint32_t>(item);
    };

    auto scan_paths = [&](const auto& paths) {
        for (uint32_t seg_id : touched) {
            state[seg_id] &= REPORTED | BOUNDARY;
        }
        touched.clear();

        for (const auto& path : paths) {
            path_segments.clear();
            repeat_delta.assign(path.size() + 1, 0);

            for (size_t i = 0; i < path.size(); ++i) {
                const uint32_t vertex_id = vertex_id_of(path[i]);
                const uint32_t seg_id = Vertex::get_segment_id(vertex_id);

                if (seg_id >= node_count) continue;
                if (state[seg_id] & BOUNDARY) continue;
                if (graph_.getNodeDeleted(seg_id)) continue;

                if ((state[seg_id] & ORIENTATION_MASK) == 0) touched.push_back(seg_id);

                state[seg_id] |= Vertex::get_is_reverse(vertex_id) ? 2u : 1u;

                if (first_pos[seg_id] == UINT32_MAX) {
                    first_pos[seg_id] = static_cast<uint32_t>(i);
                    path_segments.push_back(seg_id);
                } else {
                    ++repeat_delta[first_pos[seg_id]];
                    --repeat_delta[i + 1];
                }
            }

            int32_t repeat_depth = 0;
            for (size_t i = 0; i < path.size(); ++i) {
                repeat_depth += repeat_delta[i];
                if (repeat_depth == 0) continue;

                const uint32_t vertex_id = vertex_id_of(path[i]);
                const uint32_t seg_id = Vertex::get_segment_id(vertex_id);

                if (seg_id >= node_count) continue;
                if (state[seg_id] & BOUNDARY) continue;
                if (graph_.getNodeDeleted(seg_id)) continue;
                state[seg_id] |= CONFLICT;
            }

            for (uint32_t seg_id : path_segments) first_pos[seg_id] = UINT32_MAX;
        }

        bool has_conflict = false;

        for (uint32_t seg_id : touched) {
            if ((state[seg_id] & ORIENTATION_MASK) == ORIENTATION_MASK) state[seg_id] |= CONFLICT;
            if (state[seg_id] & CONFLICT) has_conflict = true;
        }

        return has_conflict;
    };

    // ------------------------------------------------ Bubbles ------------------------------------------------
    for (const Bubble& bubble : bubbles_) {
        const std::vector<std::vector<uint32_t>>& paths = bubble.get_paths();
        if (paths.size() < 2) continue;

        mark_boundary(Vertex::get_segment_id(bubble.get_source()));
        mark_boundary(Vertex::get_segment_id(bubble.get_sink()));

        const bool has_conflict = scan_paths(paths);
        clear_boundaries();
        if (!has_conflict) continue;

        for (uint32_t seg_id : touched) {
            add_segment(seg_id);
        }
    }

    // ------------------------------------------------ Same-source homologous paths ------------------------------------------------
    for (const HomologousPath& hp : homologous_paths_) {
        if (hp.type == Type::DiffSource) continue;
        if (hp.paths.size() < 2) continue;

        for (Vertex v : hp.starts) mark_boundary(v.segment_id());
        for (Vertex v : hp.ends) mark_boundary(v.segment_id());

        const bool has_conflict = scan_paths(hp.paths);
        clear_boundaries();
        if (!has_conflict) continue;

        for (uint32_t seg_id : touched) {
            if (!(state[seg_id] & CONFLICT)) continue;
            add_segment(seg_id);
        }
    }

    return out;
}


/* ================================================================================================================
 *                                         BUBBLE WRITER START
 * ================================================================================================================ */
BubbleWriter::BubbleWriter(
    const GfaGraph& graph,
    const std::vector<Bubble>& bubbles,
    const std::vector<HomologousPath>& homologous_paths
)
    : graph_(graph), bubbles_(bubbles), homologous_paths_(homologous_paths)
{}

void BubbleWriter::save_gfa(
    const std::string& output_file,
    const uint32_t min_len,
    const uint32_t min_num,
    bool write_seq,
    const std::string& command_line
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

    oss << "H\tVN:Z:1.0\n";
    oss << "H\tTS:Z:liftasm\n";
    if (!command_line.empty()) {
        std::string cl = command_line;
        for (char& c : cl) {
            if (c == '\t' || c == '\n' || c == '\r') c = ' ';
        }
        oss << "H\tCL:Z:" << cl << "\n";
    }
    oss << "H\tCO:Z:  ID:i   record id\n";
    oss << "H\tCO:Z:  TP:Z   record type (bubble | homo)\n";
    oss << "H\tCO:Z:  SR:Z   source node(s)\n";
    oss << "H\tCO:Z:  SK:Z   sink node(s)\n";
    oss << "H\tCO:Z:  TY:Z   subtype\n";
    oss << "H\tCO:Z:  LN:i   total length (bp) of inner nodes\n";
    oss << "H\tCO:Z:  NC:i   number of unique inner nodes\n";
    oss << "H\tCO:Z:  CG:i   path cluster id\n";
    oss << "H\tCO:Z:  CX:i   complex node flag\n";
    oss << "H\tCO:Z:  VT:Z   variant category from member nodes\n";
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
        const std::string variant_types = collect_variant_types_(inner_segs);

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
            << "\tNC:i:" << inner_nodes;
        if (!variant_types.empty()) oss << "\tVT:Z:" << variant_types;
        oss << "\n";
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
        const std::string variant_types = collect_variant_types_(inner_segs);

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
            << "\tNC:i:" << inner_nodes;
        if (!variant_types.empty()) oss << "\tVT:Z:" << variant_types;
        oss << "\n";
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
        oss << "\tLN:i:" << graph_.getNodeLength(sid);
        if (graph_.getNodeComplex(sid)) oss << "\tCX:i:1";
        const std::string variant_types = node_variant_types_(sid);
        if (!variant_types.empty()) oss << "\tVT:Z:" << variant_types;
        oss << "\n";

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

            const std::string seg_list_str = path_name_u32_(vertex_list);
            const auto overlaps = graph_.get_path_overlaps(vertex_list);
            const std::string ovlp_str = overlap_string_(overlaps);
            const std::string seq = graph_.get_path_sequence(vertex_list);

            uint64_t path_len = 0;
            if (!seq.empty() && seq!="*") {
                path_len = seq.size();
            } else {
                path_len = path_span_u32_(vertex_list, overlaps, bb.get_source(), bb.get_sink());
            }

            oss << "P\t"
                << "b" << cur_id << ".p" << pid << "\t"
                << seg_list_str << "\t"
                << ovlp_str << "\t"
                << "ID:i:" << cur_id << "\t"
                << "TP:Z:bubble\t"
                << "PI:i:" << pid << "\t"
                << "CG:i:" << (pid < bb.get_path_clusters().size() ? bb.get_path_clusters()[pid] : pid) << "\t"
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

            const std::string seg_list_str = path_name_vtx_(p);
            const auto overlaps = graph_.get_path_overlaps(p);
            const std::string ovlp_str = overlap_string_(overlaps);
            const std::string seq = graph_.get_path_sequence(p);

            uint64_t path_len = 0;
            if (!seq.empty() && seq!="*") {
                path_len = seq.size();
            } else {
                path_len = path_span_vtx_(p, overlaps);
            }

            oss << "P\t"
                << "h" << cur_id << ".p" << pid << "\t"
                << seg_list_str << "\t"
                << ovlp_str << "\t"
                << "ID:i:" << cur_id << "\t"
                << "TP:Z:homo\t"
                << "PI:i:" << pid << "\t"
                << "CG:i:" << (pid < hp.get_path_clusters().size() ? hp.get_path_clusters()[pid] : pid) << "\t"
                << "LN:i:" << path_len << "\t"
                << "SEQ:Z:" << (seq.empty() || !write_seq ? "*" : seq)
                << "\n";

            saver.save(oss.str()); oss.str(""); oss.clear();
        }
    }

    log_stream() << "  - Total records written: " << nb + nh << "\n\n";
}

// ------------------------------------------------ Contig P-line genotype index ------------------------------------------------
namespace {

static bool is_component_path_name_(std::string_view name)
{
    constexpr std::string_view component = "component";
    if (name.starts_with(component)) {
        const size_t first = component.size();
        if (first < name.size() && std::all_of(name.begin() + first, name.end(), [](unsigned char c) {
                return std::isdigit(c);
            })) return true;
    }

    const size_t marker = name.rfind(".component");
    if (marker == std::string_view::npos) return false;
    const size_t first = marker + component.size() + 1;
    return first < name.size() &&
        std::all_of(name.begin() + first, name.end(), [](unsigned char c) {
            return std::isdigit(c);
        });
}

// Match complete bubble alleles against real contig P-lines. Component paths
// describe graph backbones and must not contribute sample genotypes.
class ContigPathIndex {
public:
    ContigPathIndex(const GfaGraph& graph, const std::vector<Bubble>& bubbles)
    {
        anchor_vertices_.reserve(bubbles.size() * 2);
        for (const Bubble& bubble : bubbles) {
            if (bubble.get_source() == UINT32_MAX || bubble.get_sink() == UINT32_MAX) continue;
            anchor_vertices_.insert(bubble.get_source());
            anchor_vertices_.insert(bubble.get_sink() ^ 1u);
        }
        occurrences_.reserve(anchor_vertices_.size());

        std::unordered_map<std::string, uint32_t> sample_name_to_id;
        sample_name_to_id.reserve(graph.getPaths().size() + graph.getSampleNames().size());
        paths_.reserve(graph.getPaths().size());

        auto intern_sample = [&](const std::string& name) {
            const auto [sample, inserted] = sample_name_to_id.emplace(
                name, static_cast<uint32_t>(sample_names_.size())
            );
            if (inserted) sample_names_.push_back(name);
            return sample->second;
        };

        for (const GfaPath& path : graph.getPaths()) {
            if (path.segments.empty() || path.segments.size() > UINT32_MAX ||
                is_component_path_name_(path.name)) continue;

            bool valid = true;
            for (const PathSegment& segment : path.segments) {
                if (segment.node_id >= graph.getNumNodes() || segment.node_id > UINT32_MAX / 2 ||
                    graph.getNodeDeleted(segment.node_id)) {
                    valid = false;
                    break;
                }
            }
            if (!valid || paths_.size() >= UINT32_MAX) continue;

            const std::string& sample_name = path.sample.empty() ? path.name : path.sample;
            const uint32_t sample = intern_sample(sample_name);

            const uint32_t path_id = static_cast<uint32_t>(paths_.size());
            paths_.push_back({&path, sample});
            for (uint32_t position = 0; position < path.segments.size(); ++position) {
                const uint32_t vertex = vertex_at_(path, position);
                if (anchor_vertices_.contains(vertex)) {
                    occurrences_[vertex].push_back({path_id, position});
                }
            }
        }

        graph_sample_to_path_sample_.assign(graph.getSampleNames().size(), UINT32_MAX);
        for (uint32_t sample = 0; sample < graph.getSampleNames().size(); ++sample) {
            if (is_haplotype_sample(graph, sample)) {
                graph_sample_to_path_sample_[sample] = intern_sample(graph.getSampleNames()[sample]);
            }
        }
    }

    bool empty() const { return paths_.empty(); }
    size_t path_count() const { return paths_.size(); }
    const std::vector<std::string>& sample_names() const { return sample_names_; }

    std::vector<uint32_t> map_graph_samples(const std::vector<uint32_t>& graph_samples) const
    {
        std::vector<uint32_t> out;
        out.reserve(graph_samples.size());
        for (uint32_t sample : graph_samples) {
            if (sample < graph_sample_to_path_sample_.size() &&
                graph_sample_to_path_sample_[sample] != UINT32_MAX) {
                out.push_back(graph_sample_to_path_sample_[sample]);
            }
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    std::vector<std::vector<uint32_t>> match_samples(const Bubble& bubble) const
    {
        const auto& alleles = bubble.get_paths();
        std::vector<std::vector<uint32_t>> out(alleles.size());

        for (uint32_t allele_id = 0; allele_id < alleles.size(); ++allele_id) {
            const std::vector<uint32_t>& allele = alleles[allele_id];
            if (allele.empty() || allele.front() != bubble.get_source() ||
                allele.back() != bubble.get_sink()) continue;

            auto add_sample = [&](const IndexedPath& path) {
                std::vector<uint32_t>& samples = out[allele_id];
                if (std::find(samples.begin(), samples.end(), path.sample_id) == samples.end()) {
                    samples.push_back(path.sample_id);
                }
            };
            find_matches_(allele, false, add_sample);
            find_matches_(allele, true, add_sample);
        }
        return out;
    }

private:
    struct Occurrence {
        uint32_t path_id{0};
        uint32_t position{0};
    };

    struct IndexedPath {
        const GfaPath* path{nullptr};
        uint32_t sample_id{UINT32_MAX};
    };

    static uint32_t vertex_at_(const GfaPath& path, uint32_t position)
    {
        const PathSegment& segment = path.segments[position];
        return Vertex::make_vertex(static_cast<uint32_t>(segment.node_id), segment.is_reverse);
    }

    template<class Callback>
    void find_matches_(const std::vector<uint32_t>& allele, bool reverse, Callback&& callback) const
    {
        const uint32_t first = reverse ? allele.back() ^ 1u : allele.front();
        const auto found = occurrences_.find(first);
        if (found == occurrences_.end()) return;

        for (const Occurrence& occurrence : found->second) {
            const IndexedPath& indexed = paths_[occurrence.path_id];
            const GfaPath& path = *indexed.path;
            if (allele.size() > path.segments.size() - occurrence.position) continue;

            bool equal = true;
            for (size_t i = 0; i < allele.size(); ++i) {
                const uint32_t expected = reverse ? allele[allele.size() - i - 1] ^ 1u : allele[i];
                if (vertex_at_(path, occurrence.position + static_cast<uint32_t>(i)) != expected) {
                    equal = false;
                    break;
                }
            }
            if (equal) callback(indexed);
        }
    }

    std::vector<std::string> sample_names_;
    std::vector<IndexedPath> paths_;
    std::unordered_set<uint32_t> anchor_vertices_;
    std::unordered_map<uint32_t, std::vector<Occurrence>> occurrences_;
    std::vector<uint32_t> graph_sample_to_path_sample_;
};

// ------------------------------------------------ Primary path coordinate index ------------------------------------------------
class PrimaryPathIndex {
public:
    struct Match {
        uint32_t path_id{UINT32_MAX};
        uint32_t allele_id{UINT32_MAX};
        uint32_t begin{0};
        uint32_t end{0};
        uint64_t source_pos{0};
        uint64_t sink_pos{0};
        bool reverse{false};
    };

    PrimaryPathIndex(
        const GfaGraph& graph,
        const std::vector<ReferencePath>& reference_paths,
        const std::vector<Bubble>& bubbles
    )
        : graph_(graph)
    {
        anchor_vertices_.reserve(bubbles.size() * 2);
        for (const Bubble& bubble : bubbles) {
            if (bubble.get_source() == UINT32_MAX || bubble.get_sink() == UINT32_MAX) continue;
            anchor_vertices_.insert(bubble.get_source());
            anchor_vertices_.insert(bubble.get_sink() ^ 1u);
        }
        occurrences_.reserve(anchor_vertices_.size());
        paths_.reserve(reference_paths.size());
        for (const ReferencePath& path : reference_paths) add_path_(path);
    }

    const std::string& name(uint32_t path_id) const { return paths_[path_id].path->name; }
    uint64_t length(uint32_t path_id) const { return paths_[path_id].length; }
    size_t size() const { return paths_.size(); }

    bool find_unique(const Bubble& bubble, Match& match, bool& multiple) const
    {
        match = {};
        multiple = false;
        bool found = false;

        const auto record = [&](const Occurrence& occurrence, uint32_t allele_id, uint32_t path_size, bool reverse) {
            const IndexedPath& path = paths_[occurrence.path_id];
            const uint32_t last = occurrence.position + path_size - 1;

            Match candidate;
            candidate.path_id = occurrence.path_id;
            candidate.allele_id = allele_id;
            candidate.begin = occurrence.position;
            candidate.end = last;
            candidate.reverse = reverse;

            if (!reverse) {
                candidate.source_pos = path.begins[occurrence.position]
                    + graph_.getNodeLength(Vertex::get_segment_id(path.path->vertices[occurrence.position])) - 1;
                candidate.sink_pos = path.begins[last];
            } else {
                candidate.sink_pos = path.begins[occurrence.position]
                    + graph_.getNodeLength(Vertex::get_segment_id(path.path->vertices[occurrence.position])) - 1;
                candidate.source_pos = path.begins[last];
            }

            if (!found) {
                match = candidate;
                found = true;
                return true;
            }

            // Multiple identical alleles at the same primary interval are one placement.
            if (match.path_id == candidate.path_id && match.begin == candidate.begin && match.end == candidate.end) return true;
            multiple = true;
            return false;
        };

        const auto& alleles = bubble.get_paths();
        for (uint32_t allele_id = 0; allele_id < alleles.size() && !multiple; ++allele_id) {
            const std::vector<uint32_t>& allele = alleles[allele_id];
            if (allele.empty() || allele.front() != bubble.get_source() || allele.back() != bubble.get_sink()) continue;

            find_matches_(allele, false, [&](const Occurrence& occurrence) {
                return record(occurrence, allele_id, static_cast<uint32_t>(allele.size()), false);
            });
            if (multiple) break;

            find_matches_(allele, true, [&](const Occurrence& occurrence) {
                return record(occurrence, allele_id, static_cast<uint32_t>(allele.size()), true);
            });
        }
        return found && !multiple;
    }

private:
    static constexpr uint64_t HASH_BASE_ = 0x9e3779b185ebca87ULL;

    struct Occurrence {
        uint32_t path_id{0};
        uint32_t position{0};
    };

    struct IndexedPath {
        const ReferencePath* path{nullptr};
        std::vector<uint64_t> begins;
        std::vector<uint64_t> hash;
        uint64_t length{0};
    };

    void add_path_(const ReferencePath& path)
    {
        IndexedPath indexed;
        indexed.path = &path;
        indexed.begins.resize(path.vertices.size());
        indexed.hash.resize(path.vertices.size() + 1, 0);

        const uint64_t vertex_count = graph_.getNumNodes() * 2;
        for (uint32_t vertex : path.vertices) {
            const uint32_t segment = Vertex::get_segment_id(vertex);
            if (vertex >= vertex_count || graph_.getNodeDeleted(segment) || graph_.getNodeLength(segment) == 0) {
                return;
            }
        }

        const std::vector<uint32_t> overlaps = graph_.get_path_overlaps(path.vertices);
        uint64_t assembled = 0;

        for (size_t i = 0; i < path.vertices.size(); ++i) {
            const uint32_t vertex = path.vertices[i];
            const uint64_t node_length = graph_.getNodeLength(Vertex::get_segment_id(vertex));
            const uint64_t overlap = i == 0 ? 0 : std::min<uint64_t>(overlaps[i - 1], node_length);
            if (overlap > assembled) return;

            indexed.begins[i] = assembled - overlap;
            indexed.hash[i + 1] = indexed.hash[i] * HASH_BASE_ + uint64_t(vertex) + 1;
            assembled = indexed.begins[i] + node_length;
        }
        indexed.length = assembled;

        const uint32_t path_id = static_cast<uint32_t>(paths_.size());
        paths_.push_back(std::move(indexed));
        for (uint32_t position = 0; position < path.vertices.size(); ++position) {
            const uint32_t vertex = path.vertices[position];
            if (anchor_vertices_.contains(vertex)) occurrences_[vertex].push_back({path_id, position});
        }
    }

    template<class Callback>
    void find_matches_(const std::vector<uint32_t>& allele, bool reverse, Callback&& callback) const
    {
        const uint32_t first = reverse ? (allele.back() ^ 1u) : allele.front();
        const auto occurrence_it = occurrences_.find(first);
        if (occurrence_it == occurrences_.end()) return;

        uint64_t expected_hash = 0;
        uint64_t hash_power = 1;
        if (!reverse) {
            for (uint32_t vertex : allele) {
                expected_hash = expected_hash * HASH_BASE_ + uint64_t(vertex) + 1;
                hash_power *= HASH_BASE_;
            }
        } else {
            for (auto it = allele.rbegin(); it != allele.rend(); ++it) {
                expected_hash = expected_hash * HASH_BASE_ + uint64_t(*it ^ 1u) + 1;
                hash_power *= HASH_BASE_;
            }
        }

        const size_t size = allele.size();
        for (const Occurrence& occurrence : occurrence_it->second) {
            const IndexedPath& path = paths_[occurrence.path_id];
            if (size > path.path->vertices.size() - occurrence.position) continue;

            const uint64_t observed_hash = path.hash[occurrence.position + size]
                - path.hash[occurrence.position] * hash_power;
            if (observed_hash != expected_hash) continue;

            bool equal = true;
            for (size_t i = 0; i < size; ++i) {
                const uint32_t expected = reverse ? (allele[size - i - 1] ^ 1u) : allele[i];
                if (path.path->vertices[occurrence.position + i] != expected) {
                    equal = false;
                    break;
                }
            }
            if (equal && !callback(occurrence)) return;
        }
    }

    const GfaGraph& graph_;
    std::vector<IndexedPath> paths_;
    std::unordered_set<uint32_t> anchor_vertices_;
    std::unordered_map<uint32_t, std::vector<Occurrence>> occurrences_;
};

} // namespace


std::vector<MosaicSite> BubbleWriter::save_vcf(
    const std::string& output_prefix,
    const std::string& paf_file,
    const std::string& ref_file,
    uint32_t min_mapq,
    uint32_t min_aln_len,
    const std::vector<ReferencePath>* reference_paths,
    double max_mosaic_freq
) const
{
    // ------------------------------------------------ Initialization ------------------------------------------------
    const bool use_primary_paths = reference_paths != nullptr;
    const bool use_liftover = !use_primary_paths && !paf_file.empty();
    const bool use_ref = !ref_file.empty();
    const bool use_external_ref_sequence = use_liftover && use_ref;
    const bool use_ref_sequence = use_primary_paths || use_external_ref_sequence;

    std::optional<PrimaryPathIndex> primary_index;
    if (use_primary_paths) primary_index.emplace(graph_, *reference_paths, bubbles_);

    const ContigPathIndex contig_paths(graph_, bubbles_);
    const bool use_p_genotypes = !contig_paths.empty();

    log_stream() << "Writing bubble VCF to " << output_prefix << ".bubbles.vcf ...\n";
    if (use_primary_paths) {
        log_stream() << "  - Using " << primary_index->size() << " primary contig paths as reference\n";
        if (primary_index->size() != reference_paths->size()) {
            warning_stream() << "  ! Skipped " << reference_paths->size() - primary_index->size()
                             << " invalid primary contig path(s)\n";
        }
    } else if (!use_liftover) {
        log_stream() << "  - No PAF provided; using bubble positions from GFA\n";
    }
    if (!use_primary_paths && !use_ref) {
        log_stream() << "  - No reference provided; REF will not be extracted\n";
    }
    if (!use_primary_paths && !use_liftover && use_ref) {
        log_stream() << "  - No PAF provided; REF will not be extracted\n";
    }
    if (use_p_genotypes) {
        log_stream() << "  - Using " << contig_paths.path_count() << " contig P-lines from "
                     << contig_paths.sample_names().size() << " samples for genotypes\n";
    }

    std::unordered_map<std::string, std::vector<liftover::Pafinfo>> pafs;
    std::unordered_map<std::string, std::vector<uint32_t>> paf_ends;
    std::unordered_map<std::string, std::vector<uint32_t>> paf_idx_ends;
    if (use_liftover) {
        paf_reader(paf_file, min_mapq, min_aln_len, pafs, paf_ends, paf_idx_ends);
    }

    liftover::SeqDB db;
    if (use_external_ref_sequence && !db.load(ref_file)) {
        error_stream() << "Load ref failed: " << ref_file << "\n";
        std::exit(1);
    }

    SAVE saver(output_prefix + ".bubbles.vcf");
    std::ostringstream oss;

    const std::vector<std::string>& sample_names = use_p_genotypes ?
        contig_paths.sample_names() : graph_.getSampleNames();


    // ------------------------------------------------ VCF Output ------------------------------------------------
    oss << "##fileformat=VCFv4.3\n";
    oss << "##source=liftasm bubble\n";
    if (use_primary_paths) {
        for (uint32_t path_id = 0; path_id < primary_index->size(); ++path_id) {
            oss << "##contig=<ID=" << primary_index->name(path_id)
                << ",length=" << primary_index->length(path_id) << ">\n";
        }
    } else {
        for (const auto& name : db.names) {
            const std::string* seq = db.get(name);
            if (!seq) continue;
            oss << "##contig=<ID=" << name << ",length=" << seq->size() << ">\n";
        }
    }
    oss << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position (1-based, inclusive)\">\n";
    oss << "##INFO=<ID=SRC,Number=1,Type=String,Description=\"Source node; reserved label separators are shown as |\">\n";
    oss << "##INFO=<ID=SNK,Number=1,Type=String,Description=\"Sink node; reserved label separators are shown as |\">\n";
    oss << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of paths\">\n";
    oss << "##INFO=<ID=LN,Number=1,Type=Integer,Description=\"Total inner bp\">\n";
    oss << "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of inner nodes\">\n";
    oss << "##INFO=<ID=PN,Number=.,Type=String,Description=\"Path names in genotype order; reserved label separators are shown as |\">\n";
    oss << "##INFO=<ID=VT,Number=.,Type=String,Description=\"Variant categories from bubble nodes\">\n";
    oss << "##INFO=<ID=MT,Number=R,Type=String,Description=\"Mutation type for REF and ALT alleles from tissue frequency\">\n";
    oss << "##INFO=<ID=TF,Number=R,Type=Float,Description=\"Fraction of callable tissues carrying REF and ALT alleles\">\n";
    oss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    oss << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    if (sample_names.empty()) {
        oss << "\t" << output_prefix;
    } else {
        for (const auto& s : sample_names) oss << "\t" << s;
    }
    oss << "\n";
    saver.save(oss.str());
    oss.str(""); oss.clear();

    std::vector<VcfRecord> vcf_records;
    vcf_records.reserve(bubbles_.size());
    uint64_t primary_unmatched = 0;
    uint64_t primary_ambiguous = 0;

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
        const std::string variant_types = collect_variant_types_(inner_segs);

        std::string chrom;
        uint32_t ref_src_pos = 0;
        uint32_t ref_sink_pos = 0;
        bool bubble_need_rc_to_ref = false;
        uint32_t primary_allele = UINT32_MAX;

        BoundaryWindow src_bw, sink_bw;
        if (use_primary_paths) {
            PrimaryPathIndex::Match match;
            bool multiple = false;
            if (!primary_index->find_unique(bb, match, multiple)) {
                multiple ? ++primary_ambiguous : ++primary_unmatched;
                continue;
            }
            if (match.source_pos > UINT32_MAX || match.sink_pos > UINT32_MAX) continue;

            chrom = primary_index->name(match.path_id);
            ref_src_pos = static_cast<uint32_t>(match.source_pos);
            ref_sink_pos = static_cast<uint32_t>(match.sink_pos);
            bubble_need_rc_to_ref = match.reverse;
            primary_allele = match.allele_id;
        } else if (!use_liftover) {
            if (!boundary_window_(src_name, true, src_len, src_bw)) continue;
            if (!boundary_window_(sink_name, false, sink_len, sink_bw)) continue;
            chrom = src_bw.chrom;
            ref_src_pos = src_bw.boundary;
            ref_sink_pos = sink_bw.boundary;
            bubble_need_rc_to_ref = (ref_src_pos > ref_sink_pos);
        } else {
            if (!boundary_window_(src_name, true, src_len, src_bw)) continue;
            if (!boundary_window_(sink_name, false, sink_len, sink_bw)) continue;
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
                if (!map_boundary_(sh, src_bw, sref)) continue;

                for (const auto& th : sink_hits) {
                    if (sh.tname != th.tname) continue;
                    if (sh.strand != th.strand) continue;

                    uint32_t tref = 0;
                    if (!map_boundary_(th, sink_bw, tref)) continue;

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
        std::vector<std::vector<uint32_t>> path_sample_ids;
        std::vector<uint32_t> bubble_path_ids;

        std::vector<std::vector<uint32_t>> allele_sample_ids = use_p_genotypes ?
            contig_paths.match_samples(bb) :
            std::vector<std::vector<uint32_t>>(bb.get_paths().size());

        for (uint32_t path_id = 0; path_id < bb.get_paths().size(); ++path_id) {
            const auto& path = bb.get_paths()[path_id];
            if (!use_p_genotypes) {
                allele_sample_ids[path_id] = path_sample_ids_(path, src_seg, sink_seg);
            }

            std::vector<uint32_t> vcf_samples;
            for (uint32_t vertex : path) {
                const uint32_t sid = Vertex::get_segment_id(vertex);
                if (sid == src_seg || sid == sink_seg || !node_is_vcf_variant(graph_, sid)) continue;
                for (uint32_t sample : graph_.getNodeSampleIds(sid)) vcf_samples.push_back(sample);
            }

            std::sort(vcf_samples.begin(), vcf_samples.end());
            vcf_samples.erase(std::unique(vcf_samples.begin(), vcf_samples.end()), vcf_samples.end());

            if (!vcf_samples.empty()) {
                if (use_p_genotypes) vcf_samples = contig_paths.map_graph_samples(vcf_samples);
                auto& samples = allele_sample_ids[path_id];
                samples.insert(samples.end(), vcf_samples.begin(), vcf_samples.end());
                std::sort(samples.begin(), samples.end());
                samples.erase(std::unique(samples.begin(), samples.end()), samples.end());
            }
        }

        path_names.reserve(bb.get_paths().size());
        allele_seqs.reserve(bb.get_paths().size());
        has_allele_seq.reserve(bb.get_paths().size());
        path_sample_ids.reserve(bb.get_paths().size());
        bubble_path_ids.reserve(bb.get_paths().size());

        for (uint32_t path_id = 0; path_id < bb.get_paths().size(); ++path_id) {
            const auto& p = bb.get_paths()[path_id];
            if (p.empty()) continue;

            path_names.push_back(path_name_u32_(p, src_seg, sink_seg));
            path_sample_ids.push_back(allele_sample_ids[path_id]);
            bubble_path_ids.push_back(path_id);

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
            rec.vt = variant_types;
            rec.path_names = std::move(path_names);

            if (!have_any_allele_seq) {
                rec.ref = ".";
                rec.alt = ".";
                rec.gt_tokens.assign(rec.path_names.size(), ".");
                rec.sample_gt_tokens = sample_gt_tokens_(rec.gt_tokens, path_sample_ids, sample_names.size());
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
                rec.alt = join_strings_(alt_alleles, ",");
            }

            rec.sample_gt_tokens = sample_gt_tokens_(rec.gt_tokens, path_sample_ids, sample_names.size());
            vcf_records.push_back(std::move(rec));
            continue;
        }

        VcfRecord rec;
        rec.chrom = chrom;
        rec.pos = beg + 1;
        rec.end = end;
        rec.src = src_name;
        rec.snk = sink_name;
        rec.ns = static_cast<uint32_t>(bb.get_paths().size());
        rec.ln = inner_bp;
        rec.nc = inner_nodes;
        rec.vt = variant_types;
        rec.path_names = std::move(path_names);

        if (!have_any_allele_seq) {
            rec.ref = ".";
            rec.alt = ".";
            rec.gt_tokens.assign(rec.path_names.size(), ".");
            rec.sample_gt_tokens = sample_gt_tokens_(rec.gt_tokens, path_sample_ids, sample_names.size());
            vcf_records.push_back(std::move(rec));
            continue;
        }

        const std::string* ref_ptr = nullptr;
        std::string ref_seq;

        if (use_primary_paths) {
            const auto selected = std::find(bubble_path_ids.begin(), bubble_path_ids.end(), primary_allele);
            if (selected == bubble_path_ids.end()) continue;

            const size_t selected_index = static_cast<size_t>(selected - bubble_path_ids.begin());
            if (!has_allele_seq[selected_index]) continue;
            ref_seq = allele_seqs[selected_index];
        } else {
            ref_ptr = db.get(chrom);
            if (!ref_ptr || end > ref_ptr->size()) continue;
            ref_seq = ref_ptr->substr(beg, end - beg);
        }
        if (ref_seq.empty()) continue;

        uint32_t out_pos1 = beg + 1;

        std::vector<std::string> called_alleles;
        called_alleles.reserve(allele_seqs.size());

        for (size_t i = 0; i < allele_seqs.size(); ++i) {
            if (!has_allele_seq[i]) continue;
            called_alleles.push_back(allele_seqs[i]);
        }

        while (!called_alleles.empty() && !first_base_all_same_(ref_seq, called_alleles)) {
            if (use_primary_paths) break;
            if (beg == 0) break;
            --beg;
            out_pos1 = beg + 1;
            ref_seq.insert(ref_seq.begin(), (*ref_ptr)[beg]);
            for (auto& a : called_alleles) {
                a.insert(a.begin(), (*ref_ptr)[beg]);
            }
        }
        if (use_primary_paths && !first_base_all_same_(ref_seq, called_alleles)) continue;

        while (can_trim_prefix_(ref_seq, called_alleles)) {
            ref_seq.erase(ref_seq.begin());
            for (auto& a : called_alleles) a.erase(a.begin());
            ++out_pos1;
        }

        while (can_trim_suffix_(ref_seq, called_alleles)) {
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
        rec.alt = join_strings_(alt_alleles, ",");
        rec.sample_gt_tokens = sample_gt_tokens_(rec.gt_tokens, path_sample_ids, sample_names.size());
        vcf_records.push_back(std::move(rec));
    }

    // Nested bubbles can normalize to the same VCF allele. Keep the larger
    // bubble as the representative and merge its sample evidence below.
    std::sort(vcf_records.begin(), vcf_records.end(), [](const VcfRecord& a, const VcfRecord& b) {
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        if (a.pos != b.pos) return a.pos < b.pos;
        if (a.ref != b.ref) return a.ref < b.ref;
        if (a.alt != b.alt) return a.alt < b.alt;
        if (a.ln != b.ln) return a.ln > b.ln;
        if (a.nc != b.nc) return a.nc > b.nc;
        if (a.ns != b.ns) return a.ns > b.ns;
        if (a.src != b.src) return a.src < b.src;
        return a.snk < b.snk;
    });

    const size_t before_merge = vcf_records.size();
    size_t kept = 0;
    for (size_t i = 0; i < vcf_records.size(); ++i) {
        if (kept > 0 && same_vcf_allele_(vcf_records[kept - 1], vcf_records[i])) {
            // A complete P-line match is already a haploid call. Keep the
            // larger representative instead of combining nested-bubble calls.
            if (!use_p_genotypes) {
                merge_sample_genotypes_(vcf_records[kept - 1], vcf_records[i]);
            }
            continue;
        }
        if (kept != i) vcf_records[kept] = std::move(vcf_records[i]);
        ++kept;
    }
    vcf_records.resize(kept);

    BubbleAnnotation annotation(sample_names, max_mosaic_freq);
    std::vector<MosaicSite> mosaic_sites;
    for (size_t i = 0; i < vcf_records.size(); ++i) {
        normalize_vcf_record_(vcf_records, i);
        annotation.annotate(vcf_records[i]);
        if (vcf_records[i].mosaic && vcf_records[i].pos > 0) {
            mosaic_sites.push_back({
                vcf_records[i].chrom,
                vcf_records[i].pos - 1,
                vcf_records[i].end
            });
        }
        saver.save(vcf_record_string_(vcf_records[i]));
    }

    if (use_primary_paths && primary_unmatched > 0) {
        warning_stream() << "  ! Skipped " << primary_unmatched
                         << " bubble(s) without an exact primary-path allele\n";
    }
    if (use_primary_paths && primary_ambiguous > 0) {
        warning_stream() << "  ! Skipped " << primary_ambiguous
                         << " bubble(s) with multiple exact primary-path placements\n";
    }
    if (before_merge > vcf_records.size()) {
        log_stream() << "  - Duplicate VCF alleles merged: "
                     << before_merge - vcf_records.size() << "\n";
    }
    log_stream() << "  - Total VCF records written: " << vcf_records.size() << "\n\n";
    return mosaic_sites;
}

std::string BubbleWriter::node_variant_types_(uint32_t sid) const
{
    const GfaNode* node = graph_.getNode(sid);
    if (!node || graph_.getNodeDeleted(sid)) return {};

    const uint8_t* tag = GfaAuxParser::get_aux_data(node->aux.aux_data, "VT");
    if (!tag || *tag != 'Z') return {};
    return reinterpret_cast<const char*>(tag + 1);
}

std::string BubbleWriter::collect_variant_types_(const std::unordered_set<uint32_t>& segments) const
{
    std::vector<std::string> values;
    for (uint32_t sid : segments) {
        const std::string csv = node_variant_types_(sid);
        size_t beg = 0;
        while (beg < csv.size()) {
            const size_t end = csv.find(',', beg);
            const size_t len = end == std::string::npos ? csv.size() - beg : end - beg;
            if (len > 0) values.emplace_back(csv.substr(beg, len));
            if (end == std::string::npos) break;
            beg = end + 1;
        }
    }

    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return join_strings_(values, ",");
}

std::string BubbleWriter::join_strings_(const std::vector<std::string>& values, std::string_view separator)
{
    if (values.empty()) return {};

    size_t size = separator.size() * (values.size() - 1);
    for (const std::string& value : values) size += value.size();

    std::string out;
    out.reserve(size);
    for (const std::string& value : values) {
        if (!out.empty()) out.append(separator);
        out += value;
    }
    return out;
}

std::string BubbleWriter::path_name_u32_(
    const std::vector<uint32_t>& path,
    uint32_t skip_a,
    uint32_t skip_b
) const {
    std::string out;
    out.reserve(path.size() * 16);
    for (uint32_t vertex : path) {
        const uint32_t sid = Vertex::get_segment_id(vertex);
        if (sid == skip_a || sid == skip_b) continue;
        out.push_back(Vertex::get_is_reverse(vertex) ? '<' : '>');
        out += graph_.getNodeName(sid);
    }
    return out.empty() ? "*" : out;
}

std::string BubbleWriter::path_name_vtx_(const std::vector<Vertex>& path) const
{
    std::string out;
    out.reserve(path.size() * 16);
    for (Vertex vertex : path) {
        out.push_back(vertex.is_reverse() ? '<' : '>');
        out += graph_.getNodeName(vertex.segment_id());
    }
    return out.empty() ? "*" : out;
}

std::string BubbleWriter::overlap_string_(const std::vector<uint32_t>& overlaps)
{
    if (overlaps.empty()) return "*";

    std::string out;
    out.reserve(overlaps.size() * 6);
    for (uint32_t overlap : overlaps) {
        if (!out.empty()) out.push_back(',');
        out += std::to_string(overlap);
        out.push_back('M');
    }
    return out;
}

uint64_t BubbleWriter::path_span_u32_(
    const std::vector<uint32_t>& path,
    const std::vector<uint32_t>& overlaps,
    uint32_t skip_a,
    uint32_t skip_b
) const {
    uint64_t total = 0;
    for (uint32_t vertex : path) {
        const uint32_t sid = Vertex::get_segment_id(vertex);
        if (vertex != skip_a && vertex != skip_b) total += graph_.getNodeLength(sid);
    }
    uint64_t overlap_bp = 0;
    for (uint32_t overlap : overlaps) overlap_bp += overlap;
    return total >= overlap_bp ? total - overlap_bp : 0;
}

uint64_t BubbleWriter::path_span_vtx_(
    const std::vector<Vertex>& path,
    const std::vector<uint32_t>& overlaps
) const {
    uint64_t total = 0;
    for (Vertex vertex : path) total += graph_.getNodeLength(vertex.segment_id());
    uint64_t overlap_bp = 0;
    for (uint32_t overlap : overlaps) overlap_bp += overlap;
    return total >= overlap_bp ? total - overlap_bp : 0;
}

bool BubbleWriter::boundary_window_(
    const std::string& name_with_sign,
    bool is_source,
    uint64_t plain_len,
    BoundaryWindow& out
) const {
    static constexpr uint32_t FLANK = 100;
    out = {};
    if (name_with_sign.empty()) return false;

    std::string base = name_with_sign;
    bool reverse = false;
    if (base.back() == '+' || base.back() == '-') {
        reverse = base.back() == '-';
        base.pop_back();
    }
    if (base.empty()) return false;

    if (base.find(':') == std::string::npos) {
        if (plain_len == 0) return false;
        out.chrom = base;
        out.win_end = static_cast<uint32_t>(std::min<uint64_t>(plain_len, UINT32_MAX));
        if (out.win_end == 0) return false;
        out.boundary = is_source ? (reverse ? 0 : out.win_end - 1) : (reverse ? out.win_end - 1 : 0);
        out.boundary_is_left = out.boundary == 0;
        return true;
    }

    gfaName parser(';');
    const auto pieces = parser.parse_composite_with_dir(parser.force_name_dir(base, reverse));
    if (pieces.empty()) return false;

    const auto& piece = is_source ? pieces.back() : pieces.front();
    if (piece.hi <= piece.lo) return false;
    out.chrom = piece.root;

    const bool take_head = is_source ? piece.rev : !piece.rev;
    if (take_head) {
        out.win_beg = static_cast<uint32_t>(piece.lo);
        out.win_end = static_cast<uint32_t>(std::min<uint64_t>(piece.hi, piece.lo + FLANK));
        out.boundary = out.win_beg;
        out.boundary_is_left = true;
    } else {
        out.win_beg = static_cast<uint32_t>(piece.hi > FLANK ? std::max<uint64_t>(piece.lo, piece.hi - FLANK) : piece.lo);
        out.win_end = static_cast<uint32_t>(piece.hi);
        if (out.win_end <= out.win_beg) return false;
        out.boundary = out.win_end - 1;
        out.boundary_is_left = false;
    }
    return out.win_end > out.win_beg;
}

bool BubbleWriter::map_boundary_(
    const liftover::LIFTresult& hit,
    const BoundaryWindow& window,
    uint32_t& ref_pos
) {
    if (hit.tend <= hit.tbeg) return false;
    if (window.boundary_is_left) {
        ref_pos = hit.strand == '+' ? static_cast<uint32_t>(hit.tbeg) : static_cast<uint32_t>(hit.tend - 1);
    } else {
        ref_pos = hit.strand == '+' ? static_cast<uint32_t>(hit.tend - 1) : static_cast<uint32_t>(hit.tbeg);
    }
    return true;
}

bool BubbleWriter::first_base_all_same_(const std::string& ref, const std::vector<std::string>& alleles)
{
    if (ref.empty()) return false;
    for (const std::string& allele : alleles) {
        if (allele.empty() || allele.front() != ref.front()) return false;
    }
    return true;
}

bool BubbleWriter::can_trim_prefix_(const std::string& ref, const std::vector<std::string>& alleles)
{
    if (ref.size() <= 1 || alleles.empty()) return false;
    for (const std::string& allele : alleles) {
        if (allele.size() <= 1 || allele.front() != ref.front()) return false;
    }
    return true;
}

bool BubbleWriter::can_trim_suffix_(const std::string& ref, const std::vector<std::string>& alleles)
{
    if (ref.size() <= 1 || alleles.empty()) return false;
    for (const std::string& allele : alleles) {
        if (allele.size() <= 1 || allele.back() != ref.back()) return false;
    }
    return true;
}

std::string BubbleWriter::path_root_(const std::string& path_name) const
{
    if (path_name.empty() || path_name == "*") return {};

    std::string name = path_name;
    if (name.front() == '<' || name.front() == '>') name.erase(name.begin());
    const size_t cut = name.find_first_of("<>,");
    if (cut != std::string::npos) name.resize(cut);

    BoundaryWindow window;
    return boundary_window_(name, true, 1, window) ? window.chrom : name;
}

std::vector<uint32_t> BubbleWriter::path_sample_ids_(
    const std::vector<uint32_t>& path,
    uint32_t src_seg,
    uint32_t sink_seg
) const {
    std::vector<uint32_t> ids;
    for (uint32_t vertex : path) {
        const uint32_t sid = Vertex::get_segment_id(vertex);
        if (sid == src_seg || sid == sink_seg) continue;
        graph_.merge_sample_ids_into(ids, graph_.getNodeSampleIds(sid));
    }
    return ids;
}

std::vector<std::string> BubbleWriter::sample_gt_tokens_(
    const std::vector<std::string>& path_gt_tokens,
    const std::vector<std::vector<uint32_t>>& path_sample_ids,
    size_t sample_count
) const {
    std::vector<std::vector<std::string>> by_sample(sample_count);
    const size_t n = std::min(path_gt_tokens.size(), path_sample_ids.size());

    for (size_t i = 0; i < n; ++i) {
        const std::string& gt = path_gt_tokens[i];
        if (gt.empty() || gt == ".") continue;
        for (uint32_t sid : path_sample_ids[i]) {
            if (sid >= by_sample.size()) continue;
            std::vector<std::string>& values = by_sample[sid];
            if (std::find(values.begin(), values.end(), gt) == values.end()) values.push_back(gt);
        }
    }

    std::vector<std::string> out;
    out.reserve(by_sample.size());
    for (std::vector<std::string>& values : by_sample) {
        if (values.empty()) {
            out.emplace_back(".");
            continue;
        }
        std::sort(values.begin(), values.end(), [](const std::string& a, const std::string& b) {
            if (a == ".") return false;
            if (b == ".") return true;
            return std::stoi(a) < std::stoi(b);
        });
        out.push_back(join_strings_(values, "/"));
    }
    return out;
}

bool BubbleWriter::same_vcf_allele_(const VcfRecord& a, const VcfRecord& b)
{
    return a.chrom == b.chrom && a.pos == b.pos && a.ref == b.ref && a.alt == b.alt;
}

std::string BubbleWriter::merge_genotypes_(std::string_view a, std::string_view b)
{
    std::vector<uint32_t> alleles;
    auto collect = [&alleles](std::string_view genotype) {
        uint32_t allele = 0;
        bool have_allele = false;
        for (size_t i = 0; i <= genotype.size(); ++i) {
            const char c = i < genotype.size() ? genotype[i] : '/';
            if (c >= '0' && c <= '9') {
                allele = allele * 10 + static_cast<uint32_t>(c - '0');
                have_allele = true;
            } else if (c == '/' || c == '|') {
                if (have_allele) alleles.push_back(allele);
                allele = 0;
                have_allele = false;
            }
        }
    };
    collect(a);
    collect(b);
    if (alleles.empty()) return ".";

    std::sort(alleles.begin(), alleles.end());
    alleles.erase(std::unique(alleles.begin(), alleles.end()), alleles.end());
    std::string merged;
    for (uint32_t allele : alleles) {
        if (!merged.empty()) merged.push_back('/');
        merged += std::to_string(allele);
    }
    return merged;
}

void BubbleWriter::merge_sample_genotypes_(VcfRecord& target, const VcfRecord& duplicate)
{
    if (target.sample_gt_tokens.size() < duplicate.sample_gt_tokens.size()) {
        target.sample_gt_tokens.resize(duplicate.sample_gt_tokens.size(), ".");
    }
    for (size_t i = 0; i < duplicate.sample_gt_tokens.size(); ++i) {
        target.sample_gt_tokens[i] = merge_genotypes_(
            target.sample_gt_tokens[i], duplicate.sample_gt_tokens[i]
        );
    }
}

// Keep graph labels readable while replacing characters that delimit VCF INFO values.
std::string BubbleWriter::vcf_info_text_(std::string_view value, bool replace_comma)
{
    std::string out;
    out.reserve(value.size());
    for (const unsigned char c : value) {
        const bool separator = c == ';' || c == '=' || c == '%' ||
            (replace_comma && c == ',') || c == '\r' || c == '\n' || c == '\t';
        out.push_back(separator ? '|' : static_cast<char>(c));
    }
    return out;
}

std::string BubbleWriter::vcf_record_string_(const VcfRecord& record)
{
    const std::vector<std::string>& gt = record.sample_gt_tokens.empty() ? record.gt_tokens : record.sample_gt_tokens;
    std::vector<std::string> path_names;
    path_names.reserve(record.path_names.size());
    for (const std::string& name : record.path_names) {
        path_names.push_back(vcf_info_text_(name));
    }

    std::ostringstream out;
    out << record.chrom << '\t' << record.pos << '\t' << record.id << '\t'
        << record.ref << '\t' << record.alt << "\t.\tPASS\t"
        << "END=" << record.end
        << ";SRC=" << vcf_info_text_(record.src)
        << ";SNK=" << vcf_info_text_(record.snk)
        << ";NS=" << record.ns
        << ";LN=" << record.ln
        << ";NC=" << record.nc
        << ";PN=" << join_strings_(path_names, ",");
    if (!record.vt.empty()) out << ";VT=" << vcf_info_text_(record.vt, false);
    if (!record.mt.empty()) out << ";MT=" << record.mt;
    if (!record.tf.empty()) out << ";TF=" << record.tf;
    out << "\tGT\t" << join_strings_(gt, "\t") << '\n';
    return out.str();
}

void BubbleWriter::normalize_vcf_record_(std::vector<VcfRecord>& records, size_t index) const
{
    VcfRecord& record = records[index];
    record.id = std::to_string(index);
    const size_t n = record.path_names.size();
    if (n <= 1 || record.gt_tokens.size() != n) return;

    std::vector<std::string> roots(n);
    std::vector<size_t> missing_indices;
    for (size_t i = 0; i < n; ++i) {
        roots[i] = path_root_(record.path_names[i]);
        if (roots[i].empty()) missing_indices.push_back(i);
    }

    std::vector<std::string> previous_roots;
    if (index > 0) {
        for (const std::string& name : records[index - 1].path_names) {
            std::string root = path_root_(name);
            if (!root.empty() && std::find(previous_roots.begin(), previous_roots.end(), root) == previous_roots.end()) {
                previous_roots.push_back(std::move(root));
            }
        }
    }

    if (!previous_roots.empty() && !missing_indices.empty()) {
        std::unordered_set<std::string> used(roots.begin(), roots.end());
        std::vector<std::string> missing_roots;
        for (const std::string& root : previous_roots) {
            if (!used.contains(root)) missing_roots.push_back(root);
        }
        const size_t count = std::min(missing_indices.size(), missing_roots.size());
        for (size_t i = 0; i < count; ++i) roots[missing_indices[i]] = missing_roots[i];
    }
    for (size_t i : missing_indices) if (roots[i].empty()) roots[i] = "~";

    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return roots[a] != roots[b] ? roots[a] < roots[b] : a < b;
    });

    bool unchanged = true;
    for (size_t i = 0; i < n; ++i) unchanged &= order[i] == i;
    if (unchanged) return;

    std::vector<std::string> names;
    std::vector<std::string> genotypes;
    names.reserve(n);
    genotypes.reserve(n);
    for (size_t i : order) {
        names.push_back(std::move(record.path_names[i]));
        genotypes.push_back(std::move(record.gt_tokens[i]));
    }
    record.path_names = std::move(names);
    record.gt_tokens = std::move(genotypes);
}
/* ================================================================================================================
 *                                         BUBBLE WRITER END
 * ================================================================================================================ */

} // namespace GfaBubble
