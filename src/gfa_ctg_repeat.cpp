#include "gfa_ctg_collapser.hpp"

#include "ThreadPool.hpp"
#include "logger.hpp"
#include "progress_tracker.hpp"

#include <algorithm>
#include <cctype>
#include <deque>
#include <future>
#include <map>
#include <numeric>
#include <unordered_set>

namespace {

struct VertexPathHash {
    size_t operator()(const std::vector<uint32_t>& path) const noexcept {
        size_t hash = path.size();
        for (uint32_t vertex : path) {
            hash ^= static_cast<size_t>(vertex) + 0x9e3779b9u + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

uint8_t nucleotide(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

char complement(char c) {
    switch (c) {
        case 'A': return 'T'; case 'a': return 't';
        case 'C': return 'G'; case 'c': return 'g';
        case 'G': return 'C'; case 'g': return 'c';
        case 'T': return 'A'; case 't': return 'a';
        default: return c;
    }
}

void append_oriented(std::string& out, const GfaNode& node, bool reverse) {
    if (!reverse) {
        out += node.sequence;
        return;
    }

    const size_t old_size = out.size();
    out.resize(old_size + node.sequence.size());
    for (size_t i = 0; i < node.sequence.size(); ++i) {
        out[old_size + i] = complement(node.sequence[node.sequence.size() - i - 1]);
    }
}

std::vector<uint64_t> build_tandem_repeat_mask(
    const std::string& sequence,
    uint32_t min_len,
    uint32_t max_period,
    uint32_t max_mismatch
) {
    const uint32_t length = static_cast<uint32_t>(sequence.size());
    std::vector<uint64_t> mask((length + 63) >> 6, 0);
    if (min_len == 0 || max_period == 0 || length < min_len) return mask;

    std::vector<std::pair<uint32_t, uint32_t>> intervals;
    const uint32_t max_p = std::min(max_period, length >> 1);

    for (uint32_t period = 1; period <= max_p; ++period) {
        std::vector<std::pair<uint32_t, uint32_t>> period_intervals;
        std::deque<uint32_t> mismatches;
        uint32_t begin = 0;

        for (uint32_t i = 0; i + period < length; ++i) {
            const uint32_t right = i + period;
            const uint8_t left_base = nucleotide(sequence[i]);
            const uint8_t right_base = nucleotide(sequence[right]);

            // Unknown bases split a repeat instead of becoming an artificial motif.
            if (right_base > 3) {
                begin = right + 1;
                mismatches.clear();
                continue;
            }
            if (left_base > 3) {
                begin = std::max(begin, i + 1);
                while (!mismatches.empty() && mismatches.front() < begin) mismatches.pop_front();
                continue;
            }

            if (left_base != right_base) mismatches.push_back(i);
            while (mismatches.size() > max_mismatch) {
                begin = mismatches.front() + 1;
                mismatches.pop_front();
            }

            const uint32_t end = right + 1;
            if (end - begin < min_len) continue;

            if (!period_intervals.empty() && begin <= period_intervals.back().second) {
                period_intervals.back().second = end;
            } else {
                period_intervals.emplace_back(begin, end);
            }
        }

        intervals.insert(intervals.end(), period_intervals.begin(), period_intervals.end());
    }

    if (intervals.empty()) return mask;
    std::sort(intervals.begin(), intervals.end());

    std::vector<std::pair<uint32_t, uint32_t>> merged;
    merged.reserve(intervals.size());
    for (const auto& interval : intervals) {
        if (!merged.empty() && interval.first <= merged.back().second) {
            merged.back().second = std::max(merged.back().second, interval.second);
        } else {
            merged.push_back(interval);
        }
    }

    for (const auto& interval : merged) {
        for (uint32_t pos = interval.first; pos < interval.second; ++pos) {
            mask[pos >> 6] |= uint64_t(1) << (pos & 63);
        }
    }
    return mask;
}

bool range_is_repeat(const std::vector<uint64_t>& mask, uint32_t begin, uint32_t end) {
    if (begin >= end) return false;
    for (uint32_t pos = begin; pos < end; ++pos) {
        if ((mask[pos >> 6] & (uint64_t(1) << (pos & 63))) == 0) return false;
    }
    return true;
}

bool interval_free(const std::map<uint32_t, uint32_t>& intervals, uint32_t begin, uint32_t end) {
    auto next = intervals.lower_bound(begin);
    if (next != intervals.end() && next->first < end) return false;
    return next == intervals.begin() || std::prev(next)->second <= begin;
}

uint64_t edge_key(uint32_t source, uint32_t sink) {
    return (uint64_t(source) << 32) | sink;
}

uint64_t canonical_edge(uint32_t source, uint32_t sink) {
    return std::min(edge_key(source, sink), edge_key(sink ^ 1u, source ^ 1u));
}

} // namespace

CtgRepeatNormalizer::CtgRepeatNormalizer(
    GfaCtgCollapser& graph,
    const CollapseOpts& options,
    BubbleOpts bubble_options
)
    : graph_(graph),
      min_repeat_len_(options.repeat_mask_min_len),
      max_repeat_period_(options.repeat_mask_max_period),
      max_repeat_mismatch_(options.repeat_mask_max_mismatch),
      bubble_options_(std::move(bubble_options))
{}

bool CtgRepeatNormalizer::is_component_path_(const std::string& name) {
    constexpr std::string_view global_prefix = "component";
    if (name.size() > global_prefix.size() && name.starts_with(global_prefix) &&
        std::all_of(name.begin() + global_prefix.size(), name.end(), [](unsigned char c) { return std::isdigit(c); })) {
        return true;
    }

    const size_t marker = name.rfind(".component");
    return marker != std::string::npos && marker + 10 < name.size() &&
        std::all_of(name.begin() + marker + 10, name.end(), [](unsigned char c) { return std::isdigit(c); });
}

uint32_t CtgRepeatNormalizer::path_vertex_(const PathSegment& segment) {
    return Vertex::make_vertex(static_cast<uint32_t>(segment.node_id), segment.is_reverse);
}

CtgRepeatNormalizer::PathAnnotation CtgRepeatNormalizer::annotate_path_(size_t path_id) const {
    PathAnnotation annotation;
    if (path_id >= graph_.paths_.size()) return annotation;

    const GfaPath& path = graph_.paths_[path_id];
    annotation.repeat.assign(path.segments.size(), 0);
    if (path.segments.empty()) return annotation;

    uint64_t total_length = 0;
    for (size_t i = 0; i < path.segments.size(); ++i) {
        const PathSegment& segment = path.segments[i];
        const uint32_t sid = static_cast<uint32_t>(segment.node_id);
        if (sid >= graph_.nodes_.size()) return annotation;
        const GfaNode& node = graph_.nodes_[sid];
        if (node.deleted || node.sequence.empty() || node.sequence == "*" || node.sequence.size() != node.length) return annotation;
        if (i > 0 && graph_.get_edge_ow(path_vertex_(path.segments[i - 1]), path_vertex_(segment)) != 0) return annotation;
        total_length += node.length;
    }
    if (total_length > UINT32_MAX) return annotation;

    std::string sequence;
    sequence.reserve(static_cast<size_t>(total_length));
    std::vector<uint32_t> offsets(path.segments.size() + 1, 0);

    for (size_t i = 0; i < path.segments.size(); ++i) {
        const PathSegment& segment = path.segments[i];
        append_oriented(sequence, graph_.nodes_[segment.node_id], segment.is_reverse);
        offsets[i + 1] = static_cast<uint32_t>(sequence.size());
    }

    const std::vector<uint64_t> mask = build_tandem_repeat_mask(
        sequence, min_repeat_len_, max_repeat_period_, max_repeat_mismatch_
    );
    for (size_t i = 0; i < path.segments.size(); ++i) {
        annotation.repeat[i] = range_is_repeat(mask, offsets[i], offsets[i + 1]);
    }
    annotation.valid = true;
    return annotation;
}

void CtgRepeatNormalizer::annotate_paths_() {
    log_stream() << "Annotating tandem repeats in contig paths ...\n";
    annotations_.assign(graph_.paths_.size(), PathAnnotation{});

    std::vector<size_t> path_ids;
    path_ids.reserve(graph_.paths_.size());
    uint64_t max_path_length = 0;
    for (size_t path_id = 0; path_id < graph_.paths_.size(); ++path_id) {
        const GfaPath& path = graph_.paths_[path_id];
        if (is_component_path_(path.name) || path.segments.empty()) continue;

        uint64_t length = 0;
        for (const PathSegment& segment : path.segments) {
            if (segment.node_id < graph_.nodes_.size()) length += graph_.nodes_[segment.node_id].length;
        }
        max_path_length = std::max(max_path_length, length);
        path_ids.push_back(path_id);
    }
    if (path_ids.empty()) {
        log_stream() << "  - Contig paths annotated: 0\n\n";
        return;
    }

    // Bound temporary sequence and bit-mask memory while still processing contigs in parallel.
    constexpr uint64_t memory_budget = 512ULL << 20;
    const uint64_t bytes_per_worker = max_path_length + ((max_path_length + 7) >> 3) + 1;
    const uint32_t memory_workers = static_cast<uint32_t>(std::max<uint64_t>(1, memory_budget / bytes_per_worker));
    const uint32_t threads = std::max<uint32_t>(1, std::min<uint32_t>(bubble_options_.threads, memory_workers));

    ThreadPool pool(threads);
    ProgressTracker progress(path_ids.size());
    std::vector<std::future<std::pair<size_t, PathAnnotation>>> futures;
    futures.reserve(path_ids.size());

    for (size_t path_id : path_ids) {
        futures.emplace_back(pool.submit([this, path_id, &progress]() {
            PathAnnotation annotation = annotate_path_(path_id);
            progress.hit();
            return std::make_pair(path_id, std::move(annotation));
        }));
    }
    size_t annotated = 0;
    for (auto& future : futures) {
        auto result = future.get();
        annotated += result.second.valid;
        annotations_[result.first] = std::move(result.second);
    }
    pool.stop();
    progress.finish();
    log_stream() << "  - Contig paths annotated: " << annotated << "\n\n";
}

void CtgRepeatNormalizer::index_bubble_starts_(const std::vector<GfaBubble::Bubble>& bubbles) {
    std::unordered_set<uint32_t> wanted;
    wanted.reserve(bubbles.size() * 2 + 1);
    for (const GfaBubble::Bubble& bubble : bubbles) {
        wanted.insert(bubble.get_source());
        wanted.insert(bubble.get_sink() ^ 1u);
    }

    starts_.clear();
    starts_.reserve(wanted.size() * 2 + 1);
    path_occurrences_.assign(graph_.nodes_.size(), 0);

    for (uint32_t path_id = 0; path_id < graph_.paths_.size(); ++path_id) {
        const GfaPath& path = graph_.paths_[path_id];
        for (uint32_t pos = 0; pos < path.segments.size(); ++pos) {
            const uint32_t sid = static_cast<uint32_t>(path.segments[pos].node_id);
            if (sid >= graph_.nodes_.size() || graph_.nodes_[sid].deleted) continue;
            ++path_occurrences_[sid];

            const uint32_t vertex = path_vertex_(path.segments[pos]);
            if (wanted.contains(vertex)) starts_[vertex].push_back({path_id, pos});
        }
    }
}

bool CtgRepeatNormalizer::build_candidate_(
    const GfaBubble::Bubble& bubble,
    Candidate& candidate
) const {
    struct Slice {
        uint32_t path{0};
        uint32_t begin{0}, end{0};
        bool reverse{false};
        bool real{false};
    };
    struct Branch {
        std::vector<uint32_t> vertices;
        std::vector<Slice> slices;
        uint64_t length{0};
        uint32_t real_occurrences{0};
        bool real_repeat{true};
    };

    const uint32_t source = bubble.get_source();
    const uint32_t sink = bubble.get_sink();
    if (source == UINT32_MAX || sink == UINT32_MAX) return false;

    std::vector<Branch> branches;
    std::unordered_map<std::vector<uint32_t>, size_t, VertexPathHash> branch_index;
    std::unordered_set<uint32_t> allowed_vertices;
    std::unordered_set<uint64_t> represented_edges;
    size_t max_steps = 0;

    for (const std::vector<uint32_t>& path : bubble.get_paths()) {
        if (path.size() < 3 || path.front() != source || path.back() != sink) return false;
        if (branch_index.contains(path)) continue;

        Branch branch;
        branch.vertices = path;
        for (size_t i = 1; i + 1 < path.size(); ++i) {
            const uint32_t sid = Vertex::get_segment_id(path[i]);
            if (sid >= graph_.nodes_.size() || graph_.nodes_[sid].deleted) return false;
            branch.length += graph_.nodes_[sid].length;
        }
        if (branch.length == 0) return false;

        const size_t id = branches.size();
        branch_index.emplace(branch.vertices, id);
        max_steps = std::max(max_steps, path.size());
        for (uint32_t vertex : path) allowed_vertices.insert(vertex);
        for (size_t i = 1; i < path.size(); ++i) {
            represented_edges.insert((uint64_t(path[i - 1]) << 32) | path[i]);
        }
        branches.push_back(std::move(branch));
    }
    if (branches.size() < 2) return false;

    // Every local graph edge must be represented by an enumerated bubble path.
    for (uint32_t vertex : allowed_vertices) {
        if (vertex != sink) {
            for (const GfaArc* arc : graph_.getArcsFromVertex(vertex)) {
                if (!arc || arc->get_del()) continue;
                const uint64_t edge = (uint64_t(vertex) << 32) | arc->get_target_vertex_id();
                if (!represented_edges.contains(edge)) return false;
            }
        }
        if (vertex != source) {
            for (const GfaArc* arc : graph_.getArcsToVertex(vertex)) {
                if (!arc || arc->get_del()) continue;
                const uint64_t edge = (uint64_t(arc->get_source_vertex_id()) << 32) | vertex;
                if (!represented_edges.contains(edge)) return false;
            }
        }
    }

    auto collect = [&](uint32_t start, uint32_t target, bool reverse) {
        const auto found = starts_.find(start);
        if (found == starts_.end()) return;

        for (const Occurrence& occurrence : found->second) {
            if (occurrence.path >= graph_.paths_.size()) continue;
            const GfaPath& path = graph_.paths_[occurrence.path];
            if (occurrence.pos >= path.segments.size() || path_vertex_(path.segments[occurrence.pos]) != start) continue;

            std::vector<uint32_t> raw;
            raw.reserve(max_steps);
            raw.push_back(start);
            uint32_t end = occurrence.pos;

            for (uint32_t pos = occurrence.pos + 1; pos < path.segments.size() && raw.size() < max_steps; ++pos) {
                const uint32_t vertex = path_vertex_(path.segments[pos]);
                if (vertex == target) {
                    raw.push_back(vertex);
                    end = pos;
                    break;
                }
                if (!allowed_vertices.contains(reverse ? vertex ^ 1u : vertex)) break;
                raw.push_back(vertex);
            }
            if (end == occurrence.pos) continue;

            std::vector<uint32_t> normalized;
            if (!reverse) {
                normalized = raw;
            } else {
                normalized.reserve(raw.size());
                for (auto it = raw.rbegin(); it != raw.rend(); ++it) normalized.push_back(*it ^ 1u);
            }

            const auto branch_found = branch_index.find(normalized);
            if (branch_found == branch_index.end()) continue;
            Branch& branch = branches[branch_found->second];

            const bool duplicate = std::any_of(branch.slices.begin(), branch.slices.end(), [&](const Slice& slice) {
                return slice.path == occurrence.path && slice.begin == occurrence.pos && slice.end == end;
            });
            if (duplicate) continue;

            const bool real = !is_component_path_(path.name);
            branch.slices.push_back({occurrence.path, occurrence.pos, end, reverse, real});
            if (!real) continue;

            ++branch.real_occurrences;
            if (occurrence.path >= annotations_.size() || !annotations_[occurrence.path].valid) {
                branch.real_repeat = false;
                continue;
            }
            const std::vector<uint8_t>& repeat = annotations_[occurrence.path].repeat;
            for (uint32_t pos = occurrence.pos + 1; pos < end; ++pos) {
                if (pos >= repeat.size() || !repeat[pos]) {
                    branch.real_repeat = false;
                    break;
                }
            }
        }
    };

    collect(source, sink, false);
    collect(sink ^ 1u, source ^ 1u, true);

    for (const Branch& branch : branches) {
        if (branch.real_occurrences == 0 || !branch.real_repeat) return false;
    }

    size_t template_id = 0;
    for (size_t i = 1; i < branches.size(); ++i) {
        if (branches[i].length > branches[template_id].length ||
            (branches[i].length == branches[template_id].length && branches[i].vertices < branches[template_id].vertices)) {
            template_id = i;
        }
    }

    const Branch& template_branch = branches[template_id];
    candidate = Candidate{};
    candidate.length = template_branch.length;

    std::unordered_set<uint32_t> template_nodes;
    for (size_t i = 1; i + 1 < template_branch.vertices.size(); ++i) {
        template_nodes.insert(Vertex::get_segment_id(template_branch.vertices[i]));
    }
    candidate.template_nodes.assign(template_nodes.begin(), template_nodes.end());

    std::unordered_set<uint32_t> removed_nodes;
    std::unordered_map<uint32_t, uint32_t> covered_occurrences;

    auto replacement = [&](bool reverse) {
        std::vector<PathSegment> result;
        result.reserve(template_branch.vertices.size() - 2);
        if (!reverse) {
            for (size_t i = 1; i + 1 < template_branch.vertices.size(); ++i) {
                const uint32_t vertex = template_branch.vertices[i];
                result.push_back({Vertex::get_segment_id(vertex), Vertex::get_is_reverse(vertex)});
            }
        } else {
            for (size_t i = template_branch.vertices.size() - 1; i-- > 1;) {
                const uint32_t vertex = template_branch.vertices[i] ^ 1u;
                result.push_back({Vertex::get_segment_id(vertex), Vertex::get_is_reverse(vertex)});
            }
        }
        return result;
    };

    for (size_t branch_id = 0; branch_id < branches.size(); ++branch_id) {
        if (branch_id == template_id) continue;
        const Branch& branch = branches[branch_id];

        for (size_t i = 1; i + 1 < branch.vertices.size(); ++i) {
            const uint32_t sid = Vertex::get_segment_id(branch.vertices[i]);
            if (!template_nodes.contains(sid)) removed_nodes.insert(sid);
            graph_.merge_sample_ids_into(candidate.sample_ids, graph_.nodes_[sid].sample_ids);
        }
        for (const Slice& slice : branch.slices) {
            candidate.edits.push_back({slice.path, slice.begin + 1, slice.end, replacement(slice.reverse)});
            const GfaPath& path = graph_.paths_[slice.path];
            for (uint32_t pos = slice.begin + 1; pos < slice.end; ++pos) {
                ++covered_occurrences[static_cast<uint32_t>(path.segments[pos].node_id)];
            }
        }
    }
    if (candidate.edits.empty() || removed_nodes.empty()) return false;

    for (uint32_t sid : removed_nodes) {
        if (sid >= path_occurrences_.size() || covered_occurrences[sid] != path_occurrences_[sid]) return false;
    }

    for (size_t branch_id = 0; branch_id < branches.size(); ++branch_id) {
        if (branch_id == template_id) continue;
        const std::vector<uint32_t>& path = branches[branch_id].vertices;
        for (size_t i = 1; i < path.size(); ++i) {
            const uint32_t source = Vertex::get_segment_id(path[i - 1]);
            const uint32_t sink = Vertex::get_segment_id(path[i]);
            if (!removed_nodes.contains(source) && !removed_nodes.contains(sink)) {
                candidate.alternative_edges.push_back(canonical_edge(path[i - 1], path[i]));
            }
        }
    }

    candidate.removed_nodes.assign(removed_nodes.begin(), removed_nodes.end());
    std::sort(candidate.sample_ids.begin(), candidate.sample_ids.end());
    candidate.sample_ids.erase(std::unique(candidate.sample_ids.begin(), candidate.sample_ids.end()), candidate.sample_ids.end());
    std::sort(candidate.alternative_edges.begin(), candidate.alternative_edges.end());
    candidate.alternative_edges.erase(
        std::unique(candidate.alternative_edges.begin(), candidate.alternative_edges.end()),
        candidate.alternative_edges.end()
    );
    std::sort(candidate.removed_nodes.begin(), candidate.removed_nodes.end());
    std::sort(candidate.template_nodes.begin(), candidate.template_nodes.end());
    return true;
}

size_t CtgRepeatNormalizer::apply_candidates_(std::vector<Candidate> candidates) {
    log_stream() << "Collapsing repeat-only bubbles ...\n";
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
        if (a.length != b.length) return a.length < b.length;
        return a.edits.size() < b.edits.size();
    });

    std::vector<std::map<uint32_t, uint32_t>> used_intervals(graph_.paths_.size());
    std::unordered_set<uint32_t> removed_nodes;
    std::unordered_set<uint32_t> protected_nodes;
    std::vector<Candidate> accepted;
    accepted.reserve(candidates.size());

    for (Candidate& candidate : candidates) {
        bool conflict = false;
        for (const Edit& edit : candidate.edits) {
            if (edit.path >= used_intervals.size() || !interval_free(used_intervals[edit.path], edit.begin, edit.end)) {
                conflict = true;
                break;
            }
        }
        if (!conflict) {
            for (uint32_t sid : candidate.removed_nodes) {
                if (removed_nodes.contains(sid) || protected_nodes.contains(sid)) {
                    conflict = true;
                    break;
                }
            }
        }
        if (!conflict) {
            for (uint32_t sid : candidate.template_nodes) {
                if (removed_nodes.contains(sid)) {
                    conflict = true;
                    break;
                }
            }
        }
        if (conflict) continue;

        for (const Edit& edit : candidate.edits) used_intervals[edit.path].emplace(edit.begin, edit.end);
        removed_nodes.insert(candidate.removed_nodes.begin(), candidate.removed_nodes.end());
        protected_nodes.insert(candidate.template_nodes.begin(), candidate.template_nodes.end());
        accepted.push_back(std::move(candidate));
    }
    if (accepted.empty()) {
        log_stream() << "  - Bubbles collapsed: 0\n\n";
        return 0;
    }

    std::vector<std::vector<Edit>> path_edits(graph_.paths_.size());
    std::unordered_set<uint64_t> alternative_edges;
    size_t alternative_edge_count = 0;
    for (const Candidate& candidate : accepted) {
        alternative_edge_count += candidate.alternative_edges.size();
    }
    alternative_edges.reserve(alternative_edge_count);
    for (Candidate& candidate : accepted) {
        for (uint32_t sid : candidate.template_nodes) {
            graph_.merge_sample_ids_into(graph_.nodes_[sid].sample_ids, candidate.sample_ids);
        }
        for (uint64_t edge : candidate.alternative_edges) {
            alternative_edges.insert(edge);
        }
        for (Edit& edit : candidate.edits) path_edits[edit.path].push_back(std::move(edit));
    }

    for (size_t path_id = 0; path_id < graph_.paths_.size(); ++path_id) {
        std::vector<Edit>& edits = path_edits[path_id];
        if (edits.empty()) continue;
        std::sort(edits.begin(), edits.end(), [](const Edit& a, const Edit& b) { return a.begin < b.begin; });

        GfaPath& path = graph_.paths_[path_id];
        size_t final_size = path.segments.size();
        for (const Edit& edit : edits) {
            final_size -= edit.end - edit.begin;
            final_size += edit.replacement.size();
        }

        std::vector<PathSegment> rewritten;
        rewritten.reserve(final_size);
        uint32_t cursor = 0;
        for (const Edit& edit : edits) {
            rewritten.insert(rewritten.end(), path.segments.begin() + cursor, path.segments.begin() + edit.begin);
            rewritten.insert(rewritten.end(), edit.replacement.begin(), edit.replacement.end());
            cursor = edit.end;
        }
        rewritten.insert(rewritten.end(), path.segments.begin() + cursor, path.segments.end());
        path.segments.swap(rewritten);
    }

    for (const GfaPath& path : graph_.paths_) {
        for (size_t i = 1; i < path.segments.size(); ++i) {
            alternative_edges.erase(canonical_edge(
                path_vertex_(path.segments[i - 1]),
                path_vertex_(path.segments[i])
            ));
        }
    }

    size_t removed_arcs = 0;
    for (GfaArc& arc : graph_.arcs_) {
        const uint64_t edge = canonical_edge(arc.get_source_vertex_id(), arc.get_target_vertex_id());
        if (!arc.get_del() && alternative_edges.contains(edge)) {
            arc.set_del(true);
            ++removed_arcs;
        }
    }

    for (uint32_t sid : removed_nodes) graph_.delete_segment(sid);
    log_stream() << "  - Bubbles collapsed: " << accepted.size() << '\n';
    log_stream() << "  - Alternative nodes removed: " << removed_nodes.size() << '\n';
    log_stream() << "  - Alternative arcs removed: " << removed_arcs << "\n\n";

    graph_.rebuild_after_edits();
    graph_.merge_linear_chains();
    graph_.rebuild_component_paths();
    return accepted.size();
}

size_t CtgRepeatNormalizer::run() {
    if (min_repeat_len_ == 0 || max_repeat_period_ == 0 || graph_.paths_.empty()) {
        log_stream() << "  - Skipped: repeat detection is disabled or no P-lines are available\n\n";
        return 0;
    }

    graph_.build_nodes_connectivity_index(false);
    graph_.build_vertex_topological_index(false);
    bubble_options_.keep_nested = false;
    bubble_options_.skip_comp = false;
    bubble_options_.threads = std::max<uint32_t>(1, bubble_options_.threads);

    GfaBubble::GfaBubbleFinder finder(graph_, bubble_options_);
    finder.find_bubbles();
    const std::vector<GfaBubble::Bubble>& bubbles = finder.get_bubbles();
    if (bubbles.empty()) {
        log_stream() << "  - Bubbles collapsed: 0\n\n";
        return 0;
    }

    annotate_paths_();
    index_bubble_starts_(bubbles);

    log_stream() << "Selecting bubbles fully covered by tandem repeats ...\n";
    std::vector<Candidate> candidates;
    candidates.reserve(bubbles.size());
    ProgressTracker progress(bubbles.size());
    for (const GfaBubble::Bubble& bubble : bubbles) {
        Candidate candidate;
        if (build_candidate_(bubble, candidate)) candidates.push_back(std::move(candidate));
        progress.hit();
    }
    progress.finish();

    log_stream() << "  - Repeat-only bubbles selected: " << candidates.size() << "\n\n";
    return apply_candidates_(std::move(candidates));
}

size_t GfaCtgCollapser::normalize_repeat_paths(
    const CollapseOpts& options,
    BubbleOpts bubble_options
) {
    log_stream() << "Detecting tandem repeats in contig paths and collapsing repeat-induced bubbles ...\n\n";
    return CtgRepeatNormalizer(*this, options, std::move(bubble_options)).run();
}
