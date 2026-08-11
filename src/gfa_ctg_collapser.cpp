#include "../include/gfa_ctg_collapser.hpp"

#include "../include/ThreadPool.hpp"
#include "../include/gfa_name.hpp"
#include "../include/gfa_name_check.hpp"
#include "../include/logger.hpp"
#include "../include/progress_tracker.hpp"
#include "../include/save.hpp"

#include <algorithm>
#include <cctype>
#include <deque>
#include <future>
#include <iomanip>
#include <map>
#include <sstream>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

extern "C" {
#include "minimap.h"
}

namespace {

using IntervalIndex = std::unordered_map<uint32_t, std::map<uint32_t, uint32_t>>;
using ScopedIntervalIndex = std::unordered_map<uint64_t, std::map<uint32_t, uint32_t>>;

bool is_global_component_path(const std::string& name) {
    constexpr std::string_view prefix = "component";
    return name.size() > prefix.size() && name.starts_with(prefix) &&
        std::all_of(name.begin() + prefix.size(), name.end(), [](unsigned char c) {
            return std::isdigit(c);
        });
}

bool is_component_path(const std::string& name) {
    if (is_global_component_path(name)) return true;
    const size_t marker = name.rfind(".component");
    if (marker == std::string::npos || marker + 10 >= name.size()) return false;
    return std::all_of(name.begin() + marker + 10, name.end(), [](unsigned char c) {
        return std::isdigit(c);
    });
}

struct CtgPathEdge {
    uint32_t next{UINT32_MAX};
    bool enabled{true};
};

using CtgPathGraph = std::unordered_map<uint32_t, std::vector<CtgPathEdge>>;

bool interval_free(const IntervalIndex& index, uint32_t sid, uint32_t beg, uint32_t end) {
    const auto found = index.find(sid);
    if (found == index.end()) return true;

    const auto& intervals = found->second;
    auto next = intervals.lower_bound(beg);
    if (next != intervals.end() && next->first < end) return false;
    return next == intervals.begin() || std::prev(next)->second <= beg;
}

void add_interval(IntervalIndex& index, uint32_t sid, uint32_t beg, uint32_t end) {
    index[sid].emplace(beg, end);
}

bool interval_free(const ScopedIntervalIndex& index, uint64_t key, uint32_t beg, uint32_t end) {
    const auto found = index.find(key);
    if (found == index.end()) return true;

    const auto& intervals = found->second;
    auto next = intervals.lower_bound(beg);
    if (next != intervals.end() && next->first < end) return false;
    return next == intervals.begin() || std::prev(next)->second <= beg;
}

void add_interval(ScopedIntervalIndex& index, uint64_t key, uint32_t beg, uint32_t end) {
    index[key].emplace(beg, end);
}

uint64_t interval_union_length(std::vector<std::pair<uint32_t, uint32_t>> intervals) {
    if (intervals.empty()) return 0;
    std::sort(intervals.begin(), intervals.end());

    uint64_t total = 0;
    uint32_t beg = intervals.front().first;
    uint32_t end = intervals.front().second;
    for (size_t i = 1; i < intervals.size(); ++i) {
        if (intervals[i].first > end) {
            total += end - beg;
            beg = intervals[i].first;
            end = intervals[i].second;
        } else {
            end = std::max(end, intervals[i].second);
        }
    }
    return total + end - beg;
}

uint32_t component_root(std::vector<uint32_t>& groups, uint32_t id) {
    uint32_t root = id;
    while (groups[root] != root) root = groups[root];
    while (groups[id] != id) {
        const uint32_t next = groups[id];
        groups[id] = root;
        id = next;
    }
    return root;
}

void merge_component_groups(std::vector<uint32_t>& groups, uint32_t a, uint32_t b) {
    a = component_root(groups, a);
    b = component_root(groups, b);
    if (a != b) groups[b] = a;
}

uint32_t left_boundary(const IntervalIndex& index, uint32_t sid, uint32_t pos) {
    const auto found = index.find(sid);
    if (found == index.end()) return 0;

    const auto& intervals = found->second;
    auto next = intervals.lower_bound(pos);
    return next == intervals.begin() ? 0 : std::prev(next)->second;
}

uint32_t right_boundary(const IntervalIndex& index, uint32_t sid, uint32_t pos, uint32_t length) {
    const auto found = index.find(sid);
    if (found == index.end()) return length;

    const auto next = found->second.lower_bound(pos);
    return next == found->second.end() ? length : next->first;
}

uint32_t path_vertex(const PathSegment& segment) {
    return Vertex::make_vertex(static_cast<uint32_t>(segment.node_id), segment.is_reverse);
}

void add_path_edge(CtgPathGraph& graph, uint32_t from, uint32_t to) {
    auto& edges = graph[from];
    for (const CtgPathEdge& edge : edges) {
        if (edge.next == to) return;
    }
    edges.push_back({to, true});
}

uint64_t path_length(const GfaPath& path, const std::vector<GfaNode>& nodes) {
    uint64_t length = 0;
    for (const PathSegment& segment : path.segments) {
        if (segment.node_id < nodes.size()) length += nodes[segment.node_id].length;
    }
    return length;
}

GfaPath longest_component_path(
    CtgPathGraph& graph,
    const std::vector<GfaNode>& nodes,
    uint64_t& removed_cycle_edges
) {
    std::vector<uint32_t> vertices;
    vertices.reserve(graph.size());
    for (auto& item : graph) {
        vertices.push_back(item.first);
        std::sort(item.second.begin(), item.second.end(), [](const CtgPathEdge& a, const CtgPathEdge& b) {
            return a.next < b.next;
        });
    }
    std::sort(vertices.begin(), vertices.end());

    struct DfsFrame { uint32_t vertex; size_t edge_index; };
    std::unordered_map<uint32_t, uint8_t> color;
    color.reserve(vertices.size() * 2 + 1);
    std::vector<uint32_t> postorder;
    postorder.reserve(vertices.size());
    std::vector<DfsFrame> stack;

    removed_cycle_edges = 0;
    for (uint32_t root : vertices) {
        if (color[root] != 0) continue;
        color[root] = 1;
        stack.push_back({root, 0});

        while (!stack.empty()) {
            DfsFrame& frame = stack.back();
            auto& edges = graph[frame.vertex];
            if (frame.edge_index == edges.size()) {
                color[frame.vertex] = 2;
                postorder.push_back(frame.vertex);
                stack.pop_back();
                continue;
            }

            CtgPathEdge& edge = edges[frame.edge_index++];
            const uint8_t next_color = color[edge.next];
            if (next_color == 1) {
                edge.enabled = false;
                ++removed_cycle_edges;
            } else if (next_color == 0) {
                color[edge.next] = 1;
                stack.push_back({edge.next, 0});
            }
        }
    }

    std::unordered_map<uint32_t, uint64_t> best_length;
    best_length.reserve(vertices.size() * 2 + 1);
    uint32_t start = UINT32_MAX;
    uint64_t longest = 0;
    for (uint32_t vertex : postorder) {
        const uint32_t sid = Vertex::get_segment_id(vertex);
        const uint64_t node_length = sid < nodes.size() ? nodes[sid].length : 0;
        uint64_t score = node_length;
        for (const CtgPathEdge& edge : graph[vertex]) {
            if (edge.enabled) score = std::max(score, node_length + best_length[edge.next]);
        }
        best_length[vertex] = score;
        if (start == UINT32_MAX || score > longest || (score == longest && vertex < start)) {
            start = vertex;
            longest = score;
        }
    }

    GfaPath path;
    if (start == UINT32_MAX) return path;

    std::unordered_set<uint32_t> visited;
    visited.reserve(vertices.size() / 2 + 1);
    uint32_t current = start;
    while (current != UINT32_MAX) {
        const uint32_t sid = Vertex::get_segment_id(current);
        if (sid >= nodes.size() || nodes[sid].deleted || !visited.insert(sid).second) break;
        path.segments.push_back({sid, Vertex::get_is_reverse(current)});

        uint32_t next = UINT32_MAX;
        uint64_t next_length = 0;
        for (const CtgPathEdge& edge : graph[current]) {
            const uint32_t next_sid = Vertex::get_segment_id(edge.next);
            if (!edge.enabled || visited.count(next_sid)) continue;
            const uint64_t score = best_length[edge.next];
            if (next == UINT32_MAX || score > next_length || (score == next_length && edge.next < next)) {
                next = edge.next;
                next_length = score;
            }
        }
        current = next;
    }
    return path;
}

} // namespace

std::string GfaCtgCollapser::file_prefix_(const std::string& path) {
    const size_t slash = path.find_last_of("/\\");
    std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
    if (name.ends_with(".gz")) name.resize(name.size() - 3);
    if (name.ends_with(".gfa")) name.resize(name.size() - 4);
    return name.empty() ? "sample" : name;
}

uint32_t GfaCtgCollapser::alignment_span_(const GfaAlignment& alignment, uint32_t segment_len) {
    if (alignment.read_end_pos <= alignment.read_start_pos || alignment.position_on_unitig >= segment_len) return 0;
    return std::min(alignment.read_end_pos - alignment.read_start_pos, segment_len - alignment.position_on_unitig);
}

bool GfaCtgCollapser::intervals_overlap_(uint32_t a_beg, uint32_t a_end, uint32_t b_beg, uint32_t b_end) {
    return a_beg < b_end && b_beg < a_end;
}

void GfaCtgCollapser::load_ctgs(
    const std::vector<std::string>& hap1_files,
    const std::vector<std::string>& hap2_files,
    const std::vector<std::string>& sample_names,
    uint32_t min_contig_len
) {
    std::vector<std::string> files;
    std::vector<std::string> names;
    std::vector<uint32_t> duplicate_groups;
    files.reserve(hap1_files.size() * 2);
    names.reserve(hap1_files.size() * 2);
    duplicate_groups.reserve(hap1_files.size() * 2);

    std::vector<std::string> samples = sample_names;
    if (samples.empty()) {
        samples.reserve(hap1_files.size());
        for (const std::string& file : hap1_files) samples.push_back(file_prefix_(file));
    }
    sample_names_ = samples;

    for (size_t i = 0; i < hap1_files.size(); ++i) {
        files.push_back(hap1_files[i]);
        names.push_back(samples[i] + ".hap1");
        duplicate_groups.push_back(static_cast<uint32_t>(i));
        files.push_back(hap2_files[i]);
        names.push_back(samples[i] + ".hap2");
        duplicate_groups.push_back(static_cast<uint32_t>(i));
    }

    load_from_GFA(files, names, /*read_links=*/false, duplicate_groups);
    paths_.clear();
    filter_short_contigs_(min_contig_len);
    index_origins_(samples);
    initialize_backbones_();
    CtgCollapseDebugger::inputs(hap1_files, hap2_files, samples);
}

void GfaCtgCollapser::filter_short_contigs_(uint32_t min_contig_len) {
    if (min_contig_len == 0) return;

    uint64_t removed_bp = 0;
    uint64_t removed_contigs = 0;
    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        if (nodes_[sid].deleted || nodes_[sid].length >= min_contig_len) continue;
        removed_bp += nodes_[sid].length;
        ++removed_contigs;
        delete_segment(sid);
    }

    const size_t old_alignments = alignments_.size();
    alignments_.erase(
        std::remove_if(alignments_.begin(), alignments_.end(), [this](const GfaAlignment& alignment) {
            return alignment.unitig_node_id >= nodes_.size() || nodes_[alignment.unitig_node_id].deleted;
        }),
        alignments_.end()
    );
    total_alignments_ = alignments_.size();

    log_stream() << "Filtering short contigs (min_length=" << min_contig_len << ") ...\n";
    log_stream() << "  - Contigs removed: " << removed_contigs << "\n";
    log_stream() << "  - Sequence removed: " << std::fixed << std::setprecision(2) << static_cast<double>(removed_bp) / 1'000'000.0 << " Mb\n";
    log_stream() << "  - Alignments removed: " << (old_alignments - alignments_.size()) << "\n\n";
}

void GfaCtgCollapser::index_origins_(const std::vector<std::string>& sample_names) {
    origins_.assign(nodes_.size(), SegmentOrigin{});
    std::vector<SegmentOrigin> by_sample_id(sample_id_to_name_.size());

    for (uint32_t sample = 0; sample < sample_names.size(); ++sample) {
        const auto hap1_it = sample_name_to_id_.find(sample_names[sample] + ".hap1");
        const auto hap2_it = sample_name_to_id_.find(sample_names[sample] + ".hap2");
        if (hap1_it == sample_name_to_id_.end() || hap2_it == sample_name_to_id_.end()) continue;
        by_sample_id[hap1_it->second] = {sample, 1};
        by_sample_id[hap2_it->second] = {sample, 2};
    }

    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        for (uint32_t sample_id : nodes_[sid].sample_ids) {
            if (sample_id < by_sample_id.size() && by_sample_id[sample_id].hap != 0) {
                origins_[sid] = by_sample_id[sample_id];
                break;
            }
        }
    }
}

void GfaCtgCollapser::namespace_segment_names_(const std::string& suffix) {
    gfa_name_check::FileRenamePlan plan(5, suffix);
    name_to_id_map_.clear();
    name_to_id_map_.reserve(nodes_.size() * 2 + 1);

    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        GfaNode& node = nodes_[sid];
        node.name = plan.register_segment_name(node.name, name_to_id_map_);
        name_to_id_map_.emplace(node.name, sid);
    }
    plan.print_log("");
}

void GfaCtgCollapser::initialize_backbones_() {
    backbones_.clear();
    backbones_.reserve(nodes_.size());

    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        if (nodes_[sid].deleted || origins_[sid].hap == 0) continue;
        const SegmentOrigin& origin = origins_[sid];
        const std::string& sample = sample_names_[origin.sample];
        backbones_.push_back({sample + "." + nodes_[sid].name, sid});
    }
}

std::vector<GfaCtgCollapser::Anchor> GfaCtgCollapser::collect_unique_shared_anchors_() const {
    log_stream() << "Collecting unique shared-read anchors ...\n";

    std::unordered_map<std::string, ReadOccurrence> reads;
    reads.reserve(alignments_.size() * 2 + 1);

    for (const GfaAlignment& alignment : alignments_) {
        if (alignment.unitig_node_id >= origins_.size()) continue;
        const SegmentOrigin& origin = origins_[alignment.unitig_node_id];
        if (origin.hap == 0) continue;

        const std::string key = std::to_string(origin.sample) + '\t' + alignment.read_name;
        ReadOccurrence& occurrence = reads[key];
        if (origin.hap == 1) {
            if (occurrence.hap1_count++ == 0) occurrence.hap1 = &alignment;
        } else {
            if (occurrence.hap2_count++ == 0) occurrence.hap2 = &alignment;
        }
    }

    std::vector<Anchor> anchors;
    anchors.reserve(reads.size());
    for (const auto& item : reads) {
        const ReadOccurrence& occurrence = item.second;
        if (occurrence.hap1_count != 1 || occurrence.hap2_count != 1) continue;

        const GfaAlignment& hap1 = *occurrence.hap1;
        const GfaAlignment& hap2 = *occurrence.hap2;
        const uint32_t a_sid = static_cast<uint32_t>(hap1.unitig_node_id);
        const uint32_t b_sid = static_cast<uint32_t>(hap2.unitig_node_id);
        const uint32_t a_span = alignment_span_(hap1, nodes_[a_sid].length);
        const uint32_t b_span = alignment_span_(hap2, nodes_[b_sid].length);
        if (a_span == 0 || b_span == 0) continue;

        anchors.push_back({
            a_sid, b_sid,
            hap1.position_on_unitig,
            hap1.position_on_unitig + a_span,
            hap2.position_on_unitig,
            hap2.position_on_unitig + b_span,
            hap1.read_is_rev != hap2.read_is_rev,
            &hap1.read_name
        });
    }

    log_stream() << "  - Unique anchors: " << anchors.size() << "\n\n";

    return anchors;
}

std::vector<GfaCtgCollapser::AnchorBlock> GfaCtgCollapser::build_anchor_blocks_(std::vector<Anchor> anchors) const {
    log_stream() << "Building and filtering anchor blocks ...\n";

    std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b) {
        return std::tie(a.a_sid, a.b_sid, a.reverse, a.a_beg, a.b_beg) < std::tie(b.a_sid, b.b_sid, b.reverse, b.a_beg, b.b_beg);
    });

    std::vector<AnchorBlock> blocks;
    size_t group_beg = 0;
    while (group_beg < anchors.size()) {
        size_t group_end = group_beg + 1;
        while (group_end < anchors.size() && anchors[group_end].a_sid == anchors[group_beg].a_sid && anchors[group_end].b_sid == anchors[group_beg].b_sid && anchors[group_end].reverse == anchors[group_beg].reverse) {
            ++group_end;
        }

        std::vector<uint32_t> tails;
        std::vector<size_t> tail_index;
        std::vector<size_t> previous(group_end - group_beg, SIZE_MAX);
        for (size_t i = group_beg; i < group_end; ++i) {
            const uint32_t key = anchors[i].reverse ? UINT32_MAX - anchors[i].b_end : anchors[i].b_beg;
            const size_t pos = std::lower_bound(tails.begin(), tails.end(), key) - tails.begin();
            if (pos > 0) previous[i - group_beg] = tail_index[pos - 1];
            if (pos == tails.size()) {
                tails.push_back(key);
                tail_index.push_back(i - group_beg);
            } else {
                tails[pos] = key;
                tail_index[pos] = i - group_beg;
            }
        }

        std::vector<size_t> chain;
        if (!tail_index.empty()) {
            for (size_t i = tail_index.back(); i != SIZE_MAX; i = previous[i]) chain.push_back(group_beg + i);
            std::reverse(chain.begin(), chain.end());
        }

        std::vector<uint8_t> selected(group_end - group_beg, 0);
        for (size_t index : chain) selected[index - group_beg] = 1;
        for (size_t i = group_beg; i < group_end; ++i) {
            CtgCollapseDebugger::anchor_decision(
                *this, anchors[i], selected[i - group_beg] ? "keep" : "drop"
            );
        }

        for (size_t index : chain) {
            const Anchor& anchor = anchors[index];
            if (!blocks.empty()) {
                AnchorBlock& block = blocks.back();
                const bool same_group = block.a_sid == anchor.a_sid && block.b_sid == anchor.b_sid && block.reverse == anchor.reverse;
                const bool touches_both = intervals_overlap_(block.a_beg, block.a_end, anchor.a_beg, anchor.a_end) && intervals_overlap_(block.b_beg, block.b_end, anchor.b_beg, anchor.b_end);
                if (same_group && touches_both) {
                    block.a_beg = std::min(block.a_beg, anchor.a_beg);
                    block.a_end = std::max(block.a_end, anchor.a_end);
                    block.b_beg = std::min(block.b_beg, anchor.b_beg);
                    block.b_end = std::max(block.b_end, anchor.b_end);
                    ++block.anchors;
                    continue;
                }
            }
            blocks.push_back({anchor.a_sid, anchor.b_sid, anchor.a_beg, anchor.a_end, anchor.b_beg, anchor.b_end, 1, anchor.reverse});
        }
        group_beg = group_end;
    }

    log_stream() << "  - Merged blocks: " << blocks.size() << "\n\n";
    CtgCollapseDebugger::blocks(*this, blocks, "monotonic-merged");

    return blocks;
}

std::vector<GfaCtgCollapser::AnchorBlock> GfaCtgCollapser::select_nonconflicting_blocks_(
    std::vector<AnchorBlock> blocks,
    uint32_t short_contig_len
) const {
    log_stream() << "Selecting non-conflicting anchor blocks ...\n";

    struct PairScore {
        uint64_t anchors{0};
        uint64_t span{0};
    };
    std::map<std::pair<uint32_t, uint32_t>, PairScore> scores;
    for (const AnchorBlock& block : blocks) {
        PairScore& score = scores[{block.a_sid, block.b_sid}];
        score.anchors += block.anchors;
        score.span += uint64_t(block.a_end - block.a_beg) + (block.b_end - block.b_beg);
    }

    std::unordered_map<uint32_t, uint32_t> best_b;
    std::unordered_map<uint32_t, uint32_t> best_a;
    for (const auto& item : scores) {
        const uint32_t a_sid = item.first.first;
        const uint32_t b_sid = item.first.second;
        const PairScore& score = item.second;

        const auto a_it = best_b.find(a_sid);
        if (a_it == best_b.end() ||
            std::tie(score.anchors, score.span, b_sid) >
            std::tie(scores.at({a_sid, a_it->second}).anchors, scores.at({a_sid, a_it->second}).span, a_it->second)) {
            best_b[a_sid] = b_sid;
        }

        const auto b_it = best_a.find(b_sid);
        if (b_it == best_a.end() ||
            std::tie(score.anchors, score.span, a_sid) >
            std::tie(scores.at({b_it->second, b_sid}).anchors, scores.at({b_it->second, b_sid}).span, b_it->second)) {
            best_a[b_sid] = a_sid;
        }
    }

    for (const auto& item : scores) {
        const uint32_t a_sid = item.first.first;
        const uint32_t b_sid = item.first.second;
        CtgCollapseDebugger::pair_score(
            *this, a_sid, b_sid, item.second.anchors, item.second.span,
            best_b[a_sid] == b_sid || best_a[b_sid] == a_sid
        );
    }

    struct PairCandidate {
        uint32_t a_sid{0};
        uint32_t b_sid{0};
        uint64_t anchors{0};
        uint64_t span{0};
    };
    std::vector<PairCandidate> candidates;
    candidates.reserve(scores.size());
    for (const auto& item : scores) {
        const uint32_t a_sid = item.first.first;
        const uint32_t b_sid = item.first.second;
        if (best_b[a_sid] != b_sid && best_a[b_sid] != a_sid) continue;
        candidates.push_back({a_sid, b_sid, item.second.anchors, item.second.span});
    }
    std::sort(candidates.begin(), candidates.end(), [](const PairCandidate& a, const PairCandidate& b) {
        if (a.anchors != b.anchors) return a.anchors > b.anchors;
        if (a.span != b.span) return a.span > b.span;
        return std::tie(a.a_sid, a.b_sid) < std::tie(b.a_sid, b.b_sid);
    });

    std::unordered_map<uint32_t, uint32_t> short_partner;
    std::unordered_set<uint64_t> accepted_pairs;
    accepted_pairs.reserve(candidates.size() * 2 + 1);
    for (const PairCandidate& pair : candidates) {
        const bool a_short = short_contig_len > 0 && nodes_[pair.a_sid].length <= short_contig_len;
        const bool b_short = short_contig_len > 0 && nodes_[pair.b_sid].length <= short_contig_len;
        const auto a_it = short_partner.find(pair.a_sid);
        const auto b_it = short_partner.find(pair.b_sid);
        if ((a_short && a_it != short_partner.end() && a_it->second != pair.b_sid) || (b_short && b_it != short_partner.end() && b_it->second != pair.a_sid)) continue;

        accepted_pairs.insert((uint64_t(pair.a_sid) << 32) | pair.b_sid);
        if (a_short) short_partner[pair.a_sid] = pair.b_sid;
        if (b_short) short_partner[pair.b_sid] = pair.a_sid;
    }

    blocks.erase(std::remove_if(blocks.begin(), blocks.end(), [&](const AnchorBlock& block) {
        const bool one_way_best = best_b[block.a_sid] == block.b_sid || best_a[block.b_sid] == block.a_sid;
        const bool accepted = accepted_pairs.count((uint64_t(block.a_sid) << 32) | block.b_sid) != 0;
        const bool remove = !one_way_best || !accepted;
        if (!one_way_best) CtgCollapseDebugger::block_decision(*this, block, "drop:not-best-from-either-side");
        else if (!accepted) CtgCollapseDebugger::block_decision(*this, block, "drop:short-contig-multiple-partners");
        return remove;
    }), blocks.end());

    std::sort(blocks.begin(), blocks.end(), [](const AnchorBlock& a, const AnchorBlock& b) {
        const uint64_t a_span = uint64_t(a.a_end - a.a_beg) + (a.b_end - a.b_beg);
        const uint64_t b_span = uint64_t(b.a_end - b.a_beg) + (b.b_end - b.b_beg);
        if (a.anchors != b.anchors) return a.anchors > b.anchors;
        if (a_span != b_span) return a_span > b_span;
        return std::tie(a.a_sid, a.a_beg, a.b_sid, a.b_beg) < std::tie(b.a_sid, b.a_beg, b.b_sid, b.b_beg);
    });

    IntervalIndex hap1_used, hap2_used;
    std::vector<AnchorBlock> kept;
    kept.reserve(blocks.size());
    for (const AnchorBlock& block : blocks) {
        if (!interval_free(hap1_used, block.a_sid, block.a_beg, block.a_end) || !interval_free(hap2_used, block.b_sid, block.b_beg, block.b_end)) {
            CtgCollapseDebugger::block_decision(*this, block, "drop:coordinate-conflict");
            continue;
        }
        add_interval(hap1_used, block.a_sid, block.a_beg, block.a_end);
        add_interval(hap2_used, block.b_sid, block.b_beg, block.b_end);
        kept.push_back(block);
        CtgCollapseDebugger::block_decision(*this, block, "keep");
    }

    std::sort(kept.begin(), kept.end(), [](const AnchorBlock& a, const AnchorBlock& b) {
        return std::tie(a.a_sid, a.a_beg, a.b_sid, a.b_beg) < std::tie(b.a_sid, b.a_beg, b.b_sid, b.b_beg);
    });

    log_stream() << "  - Selected blocks: " << blocks.size() << "\n\n";
    CtgCollapseDebugger::blocks(*this, kept, "selected");

    return kept;
}

std::vector<GfaCtgCollapser::AnchorBlock> GfaCtgCollapser::filter_spanning_blocks_(
    std::vector<AnchorBlock> blocks,
    double min_coverage
) const {
    if (blocks.empty() || min_coverage <= 0.0) return blocks;

    struct Bounds {
        uint32_t a_beg{UINT32_MAX}, a_end{0};
        uint32_t b_beg{UINT32_MAX}, b_end{0};
    };
    using Key = std::tuple<uint32_t, uint32_t, bool>;
    std::map<Key, Bounds> bounds;
    for (const AnchorBlock& block : blocks) {
        Bounds& span = bounds[{block.a_sid, block.b_sid, block.reverse}];
        span.a_beg = std::min(span.a_beg, block.a_beg);
        span.a_end = std::max(span.a_end, block.a_end);
        span.b_beg = std::min(span.b_beg, block.b_beg);
        span.b_end = std::max(span.b_end, block.b_end);
    }

    std::unordered_set<uint64_t> supported;
    supported.reserve(bounds.size() * 2 + 1);
    for (const auto& item : bounds) {
        const auto& [a_sid, b_sid, reverse] = item.first;
        const Bounds& span = item.second;
        const double a_cov = nodes_[a_sid].length == 0 ? 0.0 : static_cast<double>(span.a_end - span.a_beg) / nodes_[a_sid].length;
        const double b_cov = nodes_[b_sid].length == 0 ? 0.0 : static_cast<double>(span.b_end - span.b_beg) / nodes_[b_sid].length;
        if (std::max(a_cov, b_cov) >= min_coverage) {
            supported.insert((uint64_t(a_sid) << 33) | (uint64_t(b_sid) << 1) | reverse);
        }
    }

    blocks.erase(std::remove_if(blocks.begin(), blocks.end(), [&](const AnchorBlock& block) {
        const uint64_t key = (uint64_t(block.a_sid) << 33) | (uint64_t(block.b_sid) << 1) | block.reverse;
        const bool remove = !supported.contains(key);
        if (remove) CtgCollapseDebugger::block_decision(*this, block, "drop:insufficient-contig-span");
        return remove;
    }), blocks.end());

    log_stream() << "Filtering anchor pairs by contig span coverage ...\n";
    log_stream() << "  - Span-supported blocks: " << blocks.size() << "\n\n";
    return blocks;
}

void GfaCtgCollapser::align_unanchored_contigs_(
    const std::vector<AnchorBlock>& blocks,
    uint32_t short_contig_len,
    const std::string& prefix
) {
    std::vector<uint8_t> anchored(nodes_.size(), 0);
    for (const AnchorBlock& block : blocks) {
        anchored[block.a_sid] = 1;
        anchored[block.b_sid] = 1;
    }

    std::vector<ComponentBackbone> contigs;
    size_t hap1_count = 0;
    for (uint8_t hap = 1; hap <= 2; ++hap) {
        for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
            if (nodes_[sid].deleted || anchored[sid] || origins_[sid].hap != hap) continue;
            const std::string sequence = get_oriented_sequence(Vertex(Vertex::make_vertex(sid, false)));
            if (sequence.empty() || sequence == "*") continue;

            ComponentBackbone contig;
            contig.name = nodes_[sid].name;
            contig.sequence = sequence;
            contig.pieces.push_back({Vertex::make_vertex(sid, false), 0, nodes_[sid].length});
            contig.sample = hap - 1;
            contigs.push_back(std::move(contig));
            if (hap == 1) ++hap1_count;
        }
    }
    if (hap1_count == 0 || hap1_count == contigs.size()) return;
    for (uint32_t i = 0; i < contigs.size(); ++i) contigs[i].id = i;

    log_stream() << "Aligning contigs without shared-read anchors ...\n";
    log_stream() << "  - Unanchored hap1 contigs: " << hap1_count << "\n";
    log_stream() << "  - Unanchored hap2 contigs: " << contigs.size() - hap1_count << "\n\n";
    std::vector<BubbleAlignment> alignments = align_component_backbones_(
        contigs, {"hap1", "hap2"}, prefix, short_contig_len
    );
    if (alignments.empty() || min_component_coverage_ <= 0.0) return;

    struct Bounds {
        uint32_t a_beg{UINT32_MAX}, a_end{0};
        uint32_t b_beg{UINT32_MAX}, b_end{0};
    };
    using Key = std::tuple<uint32_t, uint32_t, bool>;
    std::map<Key, Bounds> bounds;
    for (const BubbleAlignment& alignment : alignments) {
        const uint32_t a_sid = Vertex::get_segment_id(alignment.v_a);
        const uint32_t b_sid = Vertex::get_segment_id(alignment.v_b);
        const bool reverse = Vertex::get_is_reverse(alignment.v_a) != Vertex::get_is_reverse(alignment.v_b);
        Bounds& span = bounds[{a_sid, b_sid, reverse}];
        span.a_beg = std::min(span.a_beg, alignment.beg_a);
        span.a_end = std::max(span.a_end, alignment.end_a);
        span.b_beg = std::min(span.b_beg, alignment.beg_b);
        span.b_end = std::max(span.b_end, alignment.end_b);
    }

    alignments.erase(std::remove_if(alignments.begin(), alignments.end(), [&](const BubbleAlignment& alignment) {
        const uint32_t a_sid = Vertex::get_segment_id(alignment.v_a);
        const uint32_t b_sid = Vertex::get_segment_id(alignment.v_b);
        const bool reverse = Vertex::get_is_reverse(alignment.v_a) != Vertex::get_is_reverse(alignment.v_b);
        const Bounds& span = bounds.at({a_sid, b_sid, reverse});
        const double a_cov = nodes_[a_sid].length == 0 ? 0.0 : static_cast<double>(span.a_end - span.a_beg) / nodes_[a_sid].length;
        const double b_cov = nodes_[b_sid].length == 0 ? 0.0 : static_cast<double>(span.b_end - span.b_beg) / nodes_[b_sid].length;
        return std::max(a_cov, b_cov) < min_component_coverage_;
    }), alignments.end());

    bubble_aligns_.insert(
        bubble_aligns_.end(),
        std::make_move_iterator(alignments.begin()),
        std::make_move_iterator(alignments.end())
    );

    log_stream() << "Span-supported unanchored alignments: " << alignments.size() << "\n\n";
    return;
}

std::vector<GfaCtgCollapser::AnchorBlock> GfaCtgCollapser::build_alignment_regions_(
    const std::vector<AnchorBlock>& anchors,
    bool anchor_only
) const {
    log_stream() << "Preparing haplotype alignment regions ...\n";

    if (anchor_only || anchors.empty()) {
        log_stream() << "  - Alignment regions: " << anchors.size() << "\n\n";
        CtgCollapseDebugger::blocks(*this, anchors, "alignment-regions");
        return anchors;
    }

    using GroupKey = std::tuple<uint32_t, uint32_t, bool>;
    std::map<GroupKey, std::vector<AnchorBlock>> groups;
    std::unordered_map<uint64_t, uint8_t> orientations;
    IntervalIndex a_anchors, b_anchors;

    std::vector<AnchorBlock> regions = anchors;
    regions.reserve(anchors.size() * 3 + 2);
    for (const AnchorBlock& anchor : anchors) {
        groups[{anchor.a_sid, anchor.b_sid, anchor.reverse}].push_back(anchor);
        orientations[(uint64_t(anchor.a_sid) << 32) | anchor.b_sid] |= anchor.reverse ? 2u : 1u;
        add_interval(a_anchors, anchor.a_sid, anchor.a_beg, anchor.a_end);
        add_interval(b_anchors, anchor.b_sid, anchor.b_beg, anchor.b_end);
    }

    auto add_region = [&](AnchorBlock region, const char* decision) {
        if (region.a_beg >= region.a_end || region.b_beg >= region.b_end) return;
        if (!interval_free(a_anchors, region.a_sid, region.a_beg, region.a_end) || !interval_free(b_anchors, region.b_sid, region.b_beg, region.b_end)) {
            CtgCollapseDebugger::block_decision(*this, region, "drop:anchor-conflict");
            return;
        }
        regions.push_back(region);
        CtgCollapseDebugger::block_decision(*this, region, decision);
    };

    for (auto& item : groups) {
        auto& blocks = item.second;
        std::sort(blocks.begin(), blocks.end(), [](const AnchorBlock& a, const AnchorBlock& b) {
            return a.a_beg < b.a_beg;
        });

        const AnchorBlock& first = blocks.front();
        const AnchorBlock& last = blocks.back();
        const uint64_t pair_key = (uint64_t(first.a_sid) << 32) | first.b_sid;
        const bool one_orientation = orientations[pair_key] == (first.reverse ? 2u : 1u);
        if (!one_orientation) continue;

        for (size_t i = 1; i < blocks.size(); ++i) {
            const AnchorBlock& left = blocks[i - 1];
            const AnchorBlock& right = blocks[i];
            AnchorBlock gap{left.a_sid, left.b_sid, left.a_end, right.a_beg, 0, 0, 0, left.reverse};
            if (left.reverse) {
                gap.b_beg = right.b_end;
                gap.b_end = left.b_beg;
            } else {
                gap.b_beg = left.b_end;
                gap.b_end = right.b_beg;
            }
            add_region(gap, "keep:inter-anchor");
        }

        AnchorBlock prefix{
            first.a_sid, first.b_sid,
            left_boundary(a_anchors, first.a_sid, first.a_beg), first.a_beg,
            0, 0, 0, first.reverse
        };
        AnchorBlock suffix{
            last.a_sid, last.b_sid,
            last.a_end, right_boundary(a_anchors, last.a_sid, last.a_end, nodes_[last.a_sid].length),
            0, 0, 0, last.reverse
        };
        if (first.reverse) {
            prefix.b_beg = first.b_end;
            prefix.b_end = right_boundary(b_anchors, first.b_sid, first.b_end, nodes_[first.b_sid].length);
            suffix.b_beg = left_boundary(b_anchors, last.b_sid, last.b_beg);
            suffix.b_end = last.b_beg;
        } else {
            prefix.b_beg = left_boundary(b_anchors, first.b_sid, first.b_beg);
            prefix.b_end = first.b_beg;
            suffix.b_beg = last.b_end;
            suffix.b_end = right_boundary(b_anchors, last.b_sid, last.b_end, nodes_[last.b_sid].length);
        }
        add_region(prefix, "keep:prefix");
        add_region(suffix, "keep:suffix");
    }

    std::sort(regions.begin(), regions.end(), [](const AnchorBlock& a, const AnchorBlock& b) {
        return std::tie(a.a_sid, a.a_beg, a.b_sid, a.b_beg) < std::tie(b.a_sid, b.a_beg, b.b_sid, b.b_beg);
    });

    log_stream() << "  - Alignment regions: " << regions.size() << "\n\n";
    CtgCollapseDebugger::blocks(*this, regions, "alignment-regions");

    return regions;
}

std::vector<GfaCtgCollapser::BubbleAlignment> GfaCtgCollapser::align_block_(const AnchorBlock& block) {
    std::string a_seq = slice_seq_or_star_(nodes_, block.a_sid, block.a_beg, block.a_end, false);
    std::string b_seq = slice_seq_or_star_(nodes_, block.b_sid, block.b_beg, block.b_end, block.reverse);
    if (a_seq.empty() || b_seq.empty() || a_seq == "*" || b_seq == "*") {
        CtgCollapseDebugger::alignment(*this, block, 0);
        return {};
    }

    std::vector<BubbleAlignment> alignments = align_and_pack_(
        Vertex::make_vertex(block.a_sid, false),
        Vertex::make_vertex(block.b_sid, block.reverse),
        nodes_[block.a_sid].name,
        nodes_[block.b_sid].name,
        a_seq, b_seq,
        block.a_beg, block.a_end,
        block.b_beg, block.b_end
    );
    CtgCollapseDebugger::alignment(*this, block, alignments.size());
    return alignments;
}

void GfaCtgCollapser::align_regions_(const std::vector<AnchorBlock>& regions) {
    log_stream() << "Aligning haplotype regions ...\n";
    ThreadPool pool(alignOpts_.threads);
    std::vector<std::future<std::vector<BubbleAlignment>>> futures;
    futures.reserve(regions.size());
    ProgressTracker progress(regions.size());
    for (const AnchorBlock& block : regions) {
        futures.emplace_back(pool.submit([this, block, &progress]() {
            std::vector<BubbleAlignment> result = align_block_(block);
            progress.hit();
            return result;
        }));
    }

    for (auto& future : futures) {
        std::vector<BubbleAlignment> result = future.get();
        bubble_aligns_.insert(bubble_aligns_.end(), std::make_move_iterator(result.begin()), std::make_move_iterator(result.end()));
    }
    progress.finish();
    pool.stop();
    log_stream() << "  - Alignment regions: " << regions.size() << "\n";
    log_stream() << "  - Alignments: " << bubble_aligns_.size() << "\n\n";
}

std::vector<GfaPath> GfaCtgCollapser::build_component_paths_(
    const std::vector<size_t>& path_ids,
    const std::string& name_prefix,
    size_t* component_count,
    uint64_t* removed_cycle_edges
) const {
    const auto& component_ids = get_nodes_connectivity_index();

    std::map<uint32_t, std::vector<size_t>> groups;
    for (size_t path_id : path_ids) {
        if (path_id >= paths_.size() || paths_[path_id].segments.empty()) continue;
        const uint64_t sid = paths_[path_id].segments.front().node_id;
        if (sid >= component_ids.size() || component_ids[sid] == UINT32_MAX) continue;
        groups[component_ids[sid]].push_back(path_id);
    }

    std::vector<GfaPath> component_paths;
    component_paths.reserve(groups.size());
    uint64_t removed = 0;
    for (const auto& item : groups) {
        const uint32_t component = item.first;
        const std::vector<size_t>& members = item.second;

        CtgPathGraph graph;
        for (size_t path_idx : members) {
            const GfaPath& path = paths_[path_idx];
            for (const PathSegment& segment : path.segments) {
                const uint32_t vertex = path_vertex(segment);
                graph.try_emplace(vertex);
                graph.try_emplace(vertex ^ 1u);
            }
            for (size_t i = 1; i < path.segments.size(); ++i) {
                const uint32_t from = path_vertex(path.segments[i - 1]);
                const uint32_t to = path_vertex(path.segments[i]);
                add_path_edge(graph, from, to);
                add_path_edge(graph, to ^ 1u, from ^ 1u);
            }
        }

        uint64_t component_cycle_edges = 0;
        GfaPath path = longest_component_path(graph, nodes_, component_cycle_edges);
        if (path.segments.empty()) continue;
        removed += component_cycle_edges;
        path.name = name_prefix + "component" + std::to_string(component);
        component_paths.push_back(std::move(path));
    }
    if (component_count) *component_count = groups.size();
    if (removed_cycle_edges) *removed_cycle_edges = removed;
    return component_paths;
}

void GfaCtgCollapser::append_component_paths_() {
    if (paths_.empty()) return;

    build_nodes_connectivity_index(false);
    log_stream() << "Building longest paths for connected components ...\n";

    std::vector<size_t> path_ids(paths_.size());
    for (size_t i = 0; i < path_ids.size(); ++i) path_ids[i] = i;
    const std::string name_prefix = sample_names_.size() == 1 ? sample_names_.front() + "." : "";
    size_t component_count = 0;
    uint64_t removed_cycle_edges = 0;
    std::vector<GfaPath> component_paths = build_component_paths_(
        path_ids, name_prefix, &component_count, &removed_cycle_edges
    );

    log_stream() << "  - Components with backbone paths: " << component_count << "\n";
    log_stream() << "  - Cycle edges removed: " << removed_cycle_edges << "\n";
    log_stream() << "  - Longest paths added: " << component_paths.size() << "\n\n";

    paths_.insert(paths_.end(), std::make_move_iterator(component_paths.begin()), std::make_move_iterator(component_paths.end()));
}

void GfaCtgCollapser::set_component_filter(double min_coverage, double end_fraction) noexcept {
    min_component_coverage_ = min_coverage;
    component_end_fraction_ = end_fraction;
}

void GfaCtgCollapser::rebuild_component_paths() {
    if (paths_.empty()) return;

    std::vector<GfaPath> kept_paths;
    kept_paths.reserve(paths_.size());
    size_t sample_component_paths = 0;
    for (GfaPath& path : paths_) {
        if (is_global_component_path(path.name)) continue;
        if (path.name.find(".component") != std::string::npos) ++sample_component_paths;
        kept_paths.push_back(std::move(path));
    }
    if (kept_paths.empty()) return;

    paths_ = std::move(kept_paths);
    build_nodes_connectivity_index(false);
    log_stream() << "Rebuilding longest paths for connected components ...\n";

    std::vector<size_t> all_paths(paths_.size());
    for (size_t i = 0; i < all_paths.size(); ++i) all_paths[i] = i;
    std::vector<GfaPath> global = build_component_paths_(all_paths, "");
    const size_t global_count = global.size();
    paths_.insert(paths_.end(), std::make_move_iterator(global.begin()), std::make_move_iterator(global.end()));

    log_stream() << "  - Sample backbone paths preserved: " << sample_component_paths << "\n";
    log_stream() << "  - Global backbone paths: " << global_count << "\n\n";
}

void GfaCtgCollapser::rebuild_backbone_paths_(const SegReplace::Expander& expander) {
    paths_.clear();
    paths_.reserve(backbones_.size() * 2);

    for (const Backbone& backbone : backbones_) {
        if (backbone.sid >= cuts_.size()) continue;
        GfaPath path;
        path.name = backbone.name;

        const auto& cuts = cuts_[backbone.sid].v;
        for (size_t i = 1; i < cuts.size(); ++i) {
            const SegReplace::Seg piece = SegReplace::Interval::pack(backbone.sid, cuts[i - 1], cuts[i], false);
            const SegReplace::Expansion expansion = expander.query(piece);
            for (SegReplace::Seg interval : expansion) {
                const uint32_t sid = materialized_segment_(interval);
                if (sid == UINT32_MAX || sid >= nodes_.size() || nodes_[sid].deleted) continue;
                const bool reverse = SegReplace::Interval::is_reverse(interval);
                if (!path.segments.empty() && path.segments.back().node_id == sid && path.segments.back().is_reverse == reverse) continue;
                path.segments.push_back({sid, reverse});
            }
        }
        if (!path.segments.empty()) {
            paths_.push_back(std::move(path));
        }
    }
}

std::vector<GfaCtgCollapser::ComponentBackbone> GfaCtgCollapser::component_backbones_(
    const std::vector<std::string>& sample_names,
    const std::unordered_map<std::string, uint32_t>* queries
) const {
    std::vector<uint32_t> sample_ids(sample_names.size(), UINT32_MAX);
    for (size_t i = 0; i < sample_names.size(); ++i) {
        const auto found = sample_name_to_id_.find(sample_names[i]);
        if (found != sample_name_to_id_.end()) sample_ids[i] = found->second;
    }

    std::vector<ComponentBackbone> result;
    const std::vector<uint32_t>& component_ids = get_nodes_connectivity_index();
    for (const GfaPath& path : paths_) {
        const bool query_path = queries && queries->contains(path.name);
        bool component_path = is_global_component_path(path.name);
        if (!queries) {
            for (const std::string& sample_name : sample_names) {
                component_path = component_path || path.name.starts_with(sample_name + ".component");
            }
        }
        if ((!component_path && !query_path) || path.segments.empty()) continue;
        const uint32_t first_sid = static_cast<uint32_t>(path.segments.front().node_id);
        if (first_sid >= nodes_.size()) continue;

        uint32_t sample = UINT32_MAX;
        if (!queries) {
            for (size_t i = 0; i < sample_names.size(); ++i) {
                if (path.name.starts_with(sample_names[i] + ".component")) {
                    sample = static_cast<uint32_t>(i);
                    break;
                }
            }
            if (sample == UINT32_MAX) {
                for (uint32_t id : nodes_[first_sid].sample_ids) {
                    const auto found = std::find(sample_ids.begin(), sample_ids.end(), id);
                    if (found != sample_ids.end()) {
                        sample = static_cast<uint32_t>(found - sample_ids.begin());
                        break;
                    }
                }
            }
            if (sample == UINT32_MAX) continue;
        }

        ComponentBackbone backbone;
        if (queries) {
            if (!query_path && (first_sid >= component_ids.size() || component_ids[first_sid] == UINT32_MAX)) continue;
            backbone.name = path.name;
            backbone.component = query_path ? queries->at(path.name) : component_ids[first_sid];
        } else {
            const std::string sample_prefix = sample_names[sample] + ".";
            backbone.name = path.name.starts_with(sample_prefix) ? path.name : sample_prefix + path.name;
            backbone.sample = sample;
        }
        for (const PathSegment& segment : path.segments) {
            const uint32_t sid = static_cast<uint32_t>(segment.node_id);
            if (sid >= nodes_.size() || nodes_[sid].deleted) continue;
            const std::string seq = get_oriented_sequence(Vertex(Vertex::make_vertex(sid, segment.is_reverse)));
            if (seq.empty() || seq == "*") {
                backbone.sequence.clear();
                backbone.pieces.clear();
                break;
            }
            const uint32_t beg = static_cast<uint32_t>(backbone.sequence.size());
            backbone.sequence += seq;
            backbone.pieces.push_back({Vertex::make_vertex(sid, segment.is_reverse), beg, static_cast<uint32_t>(backbone.sequence.size())});
        }
        if (!backbone.pieces.empty()) result.push_back(std::move(backbone));
    }
    for (uint32_t i = 0; i < result.size(); ++i) result[i].id = i;
    return result;
}

std::vector<GfaCtgCollapser::BubbleAlignment> GfaCtgCollapser::split_component_alignment_(
    const ComponentBackbone& ref,
    const ComponentBackbone& query,
    bool query_reverse,
    uint32_t ref_beg,
    uint32_t query_beg,
    uint8_t mapq,
    const std::vector<CIGAR::COp>& input_ops
) const {
    std::vector<PathPiece> query_pieces;
    if (!query_reverse) {
        query_pieces = query.pieces;
    } else {
        query_pieces.reserve(query.pieces.size());
        uint32_t pos = 0;
        for (auto it = query.pieces.rbegin(); it != query.pieces.rend(); ++it) {
            const uint32_t len = it->end - it->beg;
            query_pieces.push_back({it->vertex ^ 1u, pos, pos + len});
            pos += len;
        }
    }

    std::vector<BubbleAlignment> result;
    size_t ri = 0, qi = 0;
    uint32_t rp = ref_beg, qp = query_beg;
    while (ri < ref.pieces.size() && rp >= ref.pieces[ri].end) ++ri;
    while (qi < query_pieces.size() && qp >= query_pieces[qi].end) ++qi;

    size_t start_ri = ri, start_qi = qi;
    uint32_t start_rp = rp, start_qp = qp;
    std::vector<CIGAR::COp> piece_ops;

    auto flush = [&]() {
        if (piece_ops.empty() || start_ri >= ref.pieces.size() || start_qi >= query_pieces.size()) return;
        const PathPiece& a = ref.pieces[start_ri];
        const PathPiece& b = query_pieces[start_qi];
        uint32_t a0 = start_rp - a.beg, a1 = rp - a.beg;
        uint32_t b0 = start_qp - b.beg, b1 = qp - b.beg;
        const uint32_t a_sid = Vertex::get_segment_id(a.vertex);
        const uint32_t b_sid = Vertex::get_segment_id(b.vertex);
        if (a0 >= a1 || b0 >= b1 || a_sid == b_sid) { piece_ops.clear(); return; }
        if (Vertex::get_is_reverse(a.vertex)) {
            const uint32_t len = nodes_[a_sid].length;
            const uint32_t x = len - a1; a1 = len - a0; a0 = x;
        }
        if (Vertex::get_is_reverse(b.vertex)) {
            const uint32_t len = nodes_[b_sid].length;
            const uint32_t x = len - b1; b1 = len - b0; b0 = x;
        }
        result.emplace_back(nodes_[a_sid].name, nodes_[b_sid].name, a0, a1, b0, b1, a.vertex, b.vertex, mapq, std::move(piece_ops));
        piece_ops.clear();
    };

    for (const CIGAR::COp& op : input_ops) {
        uint32_t left = op.len;
        const bool consume_ref = path_pair_split::consumes_a(op.op);
        const bool consume_query = path_pair_split::consumes_b(op.op);
        while (left && ri < ref.pieces.size() && qi < query_pieces.size()) {
            uint32_t step = left;
            if (consume_ref) step = std::min(step, ref.pieces[ri].end - rp);
            if (consume_query) step = std::min(step, query_pieces[qi].end - qp);
            if (step == 0) {
                flush();
                while (ri < ref.pieces.size() && rp >= ref.pieces[ri].end) ++ri;
                while (qi < query_pieces.size() && qp >= query_pieces[qi].end) ++qi;
                start_ri = ri; start_qi = qi; start_rp = rp; start_qp = qp;
                continue;
            }
            if (!piece_ops.empty() && piece_ops.back().op == op.op) piece_ops.back().len += step;
            else piece_ops.emplace_back(step, op.op);
            if (consume_ref) rp += step;
            if (consume_query) qp += step;
            left -= step;
            if ((consume_ref && rp == ref.pieces[ri].end) || (consume_query && qp == query_pieces[qi].end)) {
                flush();
                while (ri < ref.pieces.size() && rp >= ref.pieces[ri].end) ++ri;
                while (qi < query_pieces.size() && qp >= query_pieces[qi].end) ++qi;
                start_ri = ri; start_qi = qi; start_rp = rp; start_qp = qp;
            }
        }
    }
    flush();
    return result;
}

void GfaCtgCollapser::collapse_ctgs(
    const std::vector<std::string>& vcf_files,
    const std::string& prefix,
    double min_anchor_coverage,
    uint32_t short_contig_len,
    bool anchor_only
) {
    log_stream() << "Collapsing haplotype contigs from unique shared-read anchors ...\n\n";

    bubble_aligns_.clear();
    cuts_.clear();
    rulemap_.clear();

    Stats variant_stats;
    if (!vcf_files.empty()) {
        variant_stats = read_vcfs_(vcf_files);
        validate_variant_order_();
        create_alt_nodes_("V");
        variant_stats.alleles = alt_nodes_.size();
    }

    initialize_cuts_();
    std::vector<Anchor> anchors = collect_unique_shared_anchors_();
    std::vector<AnchorBlock> blocks = build_anchor_blocks_(std::move(anchors));
    blocks = select_nonconflicting_blocks_(std::move(blocks), short_contig_len);
    blocks = filter_spanning_blocks_(std::move(blocks), min_anchor_coverage);
    const std::vector<AnchorBlock> regions = build_alignment_regions_(blocks, anchor_only);
    align_regions_(regions);
    align_unanchored_contigs_(blocks, short_contig_len, prefix);
    if (!alt_nodes_.empty()) inject_variant_alignments_();

    const auto groups = build_align_groups_();
    build_propagate_prune_cuts_(groups);
    if (DEBUG_ENABLED) print_cuts(cuts_);
    build_rulemap_(groups);

    SegReplace::Expander expander = build_SegReplace_();
    std::vector<GfaPath> backbone_paths;
    backbone_paths.reserve(backbones_.size());
    for (const auto& backbone : backbones_) {
        if (backbone.sid >= nodes_.size() || nodes_[backbone.sid].deleted) continue;
        GfaPath path;
        path.name = backbone.name;
        path.segments.push_back({backbone.sid, false});
        backbone_paths.push_back(std::move(path));
    }
    GfaPathCycleGuard(*this).filter_rules(backbone_paths, expander);
    expander.save_map(prefix + ".collapse.map");
    expand_and_rewire_edges_(expander);
    if (!alt_nodes_.empty()) connect_alt_nodes_(expander);
    rebuild_backbone_paths_(expander);
    remove_unused_nodes_();
    append_component_paths_();
    CtgCollapseDebugger::paths(*this);
    if (!alt_nodes_.empty()) rename_alt_nodes_();

    std::string message = "Finished haplotype contig collapse";
    if (!vcf_files.empty()) message += " with " + std::to_string(variant_stats.alleles) + " VCF allele(s)";
    log_stream() << message << ".\n\n";
}

std::vector<std::string> GfaCtgCollapser::collapse_samples(
    const std::vector<std::string>& hap1_files,
    const std::vector<std::string>& hap2_files,
    const std::vector<std::string>& sample_names,
    const std::vector<std::string>& vcf_files,
    const std::string& prefix,
    double min_anchor_coverage,
    uint32_t short_contig_len,
    uint32_t min_contig_len,
    bool anchor_only,
    const std::string& command_line
) const {
    std::vector<std::string> names = sample_names;
    if (names.empty()) {
        names.reserve(hap1_files.size());
        for (const std::string& file : hap1_files) names.push_back(file_prefix_(file));
    }

    std::vector<std::string> all_files;
    std::vector<uint32_t> sample_groups;
    all_files.reserve(hap1_files.size() * 2);
    sample_groups.reserve(hap1_files.size() * 2);
    for (size_t i = 0; i < hap1_files.size(); ++i) {
        all_files.push_back(hap1_files[i]);
        sample_groups.push_back(static_cast<uint32_t>(i));
        all_files.push_back(hap2_files[i]);
        sample_groups.push_back(static_cast<uint32_t>(i));
    }
    const bool namespace_samples = input_segment_names_need_namespace_(all_files, sample_groups);
    if (namespace_samples) {
        warning_stream() << "Duplicate segment names found across contig samples; " << "adding a sample suffix to every segment name\n\n";
    }

    std::unordered_set<std::string> used_names;
    std::vector<std::string> outputs;
    std::vector<std::string> output_names;
    outputs.reserve(hap1_files.size());
    output_names.reserve(hap1_files.size());

    for (size_t i = 0; i < hap1_files.size(); ++i) {
        const std::string base_name = names[i].empty() ? "sample" + std::to_string(i + 1) : names[i];
        std::string unique_name = base_name;
        for (uint32_t suffix = 1; !used_names.insert(unique_name).second; ++suffix) {
            unique_name = base_name + "_" + std::to_string(suffix);
        }
        if (unique_name != base_name) {
            warning_stream() << "  ! Duplicate sample name '" << base_name << "'; using '" << unique_name << "'\n";
        }

        const std::string sample_prefix = prefix + "." + unique_name;
        const std::vector<std::string> sample_vcf = vcf_files.empty() ? std::vector<std::string>{} : std::vector<std::string>{vcf_files[i]};

        log_stream() << "Processing contig sample " << (i + 1) << '/' << hap1_files.size() << " (" << unique_name << ") ...\n\n";
        GfaCtgCollapser sample(*this);
        sample.load_ctgs({hap1_files[i]}, {hap2_files[i]}, {unique_name}, min_contig_len);
        if (namespace_samples) sample.namespace_segment_names_("_" + std::to_string(i + 1));
        sample.collapse_ctgs(
            sample_vcf, sample_prefix, min_anchor_coverage,
            short_contig_len, anchor_only
        );
        sample.save_to_disk(sample_prefix + ".collapse.gfa", true, false, true, command_line);
        sample.save_to_disk(sample_prefix + ".collapse.noseq.gfa", true, false, false, command_line);
        outputs.push_back(sample_prefix + ".collapse.gfa");
        output_names.push_back(unique_name);
    }

    if (outputs.size() > 1) {
        GfaCtgCollapser combined(*this);
        combined.collapse_backbone_samples_(outputs, output_names, prefix, short_contig_len, command_line);
        return {prefix + ".iter0.collapse.gfa"};
    }
    return outputs;
}

std::vector<std::vector<GfaCtgCollapser::ComponentHit>> GfaCtgCollapser::select_component_hits_(
    std::vector<std::vector<ComponentHit>> hits,
    const std::vector<ComponentBackbone>& backbones,
    uint32_t short_contig_len,
    std::vector<uint32_t>& component_groups,
    bool place_queries
) const {
    struct PairGroup {
        uint32_t qid{0};
        uint32_t rid{0};
        uint32_t query_beg{UINT32_MAX}, query_end{0};
        uint32_t raw_query_beg{UINT32_MAX}, raw_query_end{0};
        uint32_t ref_beg{UINT32_MAX}, ref_end{0};
        uint64_t span{0};
        uint64_t matches{0};
        uint64_t block_len{0};
        bool short_component{false};
        bool strong{false};
        std::vector<ComponentHit> hits;
    };

    std::vector<PairGroup> pairs;
    for (uint32_t qid = 0; qid < hits.size() && qid < backbones.size(); ++qid) {
        const uint32_t query_len = static_cast<uint32_t>(backbones[qid].sequence.size());
        std::map<std::pair<uint32_t, bool>, std::vector<ComponentHit>> pair_hits;
        for (ComponentHit& hit : hits[qid]) {
            if (hit.rid < backbones.size()) pair_hits[{hit.rid, hit.reverse}].push_back(std::move(hit));
        }

        for (auto& item : pair_hits) {
            std::vector<ComponentHit>& candidates = item.second;
            const uint32_t rid = item.first.first;
            const uint64_t ref_len = backbones[rid].sequence.size();
            const bool short_query = !place_queries && short_contig_len > 0 && query_len <= short_contig_len;
            const bool short_reference = !place_queries && short_contig_len > 0 && ref_len <= short_contig_len;

            if (short_query || short_reference) {
                std::vector<ComponentHit> eligible;
                eligible.reserve(candidates.size());
                for (ComponentHit& candidate : candidates) {
                    const uint64_t aligned_len = std::min(
                        candidate.ref_end - candidate.ref_beg,
                        candidate.query_end - candidate.query_beg
                    );
                    if (short_query && static_cast<double>(aligned_len) / query_len < MIN_JACCARD_FOR_ALIGN_) continue;
                    if (short_reference && static_cast<double>(aligned_len) / ref_len < MIN_JACCARD_FOR_ALIGN_) continue;
                    eligible.push_back(std::move(candidate));
                }
                candidates = std::move(eligible);
                if (candidates.empty()) continue;
            }

            std::vector<std::pair<uint32_t, uint32_t>> query_intervals, ref_intervals;
            query_intervals.reserve(candidates.size());
            ref_intervals.reserve(candidates.size());
            for (const ComponentHit& candidate : candidates) {
                query_intervals.emplace_back(candidate.raw_query_beg, candidate.raw_query_end);
                ref_intervals.emplace_back(candidate.ref_beg, candidate.ref_end);
            }
            const uint64_t support_span = std::min(
                interval_union_length(std::move(query_intervals)),
                interval_union_length(std::move(ref_intervals))
            );

            std::sort(candidates.begin(), candidates.end(), [short_query, short_reference](const ComponentHit& a, const ComponentHit& b) {
                const uint32_t a_len = std::min(a.ref_end - a.ref_beg, a.query_end - a.query_beg);
                const uint32_t b_len = std::min(b.ref_end - b.ref_beg, b.query_end - b.query_beg);
                if (short_query || short_reference) {
                    const uint64_t lhs = static_cast<uint64_t>(a.matches) * b.block_len;
                    const uint64_t rhs = static_cast<uint64_t>(b.matches) * a.block_len;
                    if (lhs != rhs) return lhs > rhs;
                }
                if (a_len != b_len) return a_len > b_len;
                if (a.matches != b.matches) return a.matches > b.matches;
                if (a.mapq != b.mapq) return a.mapq > b.mapq;
                return std::tie(a.ref_beg, a.query_beg) < std::tie(b.ref_beg, b.query_beg);
            });

            PairGroup group;
            group.qid = qid;
            group.rid = rid;
            group.short_component = short_query || short_reference;
            for (ComponentHit& candidate : candidates) {
                bool accepted = true;
                for (const ComponentHit& selected : group.hits) {
                    const bool query_overlap = intervals_overlap_(
                        candidate.query_beg, candidate.query_end, selected.query_beg, selected.query_end
                    );
                    const bool ref_overlap = intervals_overlap_(
                        candidate.ref_beg, candidate.ref_end, selected.ref_beg, selected.ref_end
                    );
                    if (query_overlap || ref_overlap) {
                        MmWfaHit hit{
                            candidate.ref_beg, candidate.ref_end,
                            candidate.query_beg, candidate.query_end,
                            candidate.mapq, candidate.cigar
                        };
                        const MmWfaHit occupied{
                            selected.ref_beg, selected.ref_end,
                            selected.query_beg, selected.query_end,
                            selected.mapq, selected.cigar
                        };
                        if (!trim_small_overlap_(hit, occupied)) {
                            accepted = false;
                            break;
                        }
                        candidate.ref_beg = hit.r_beg;
                        candidate.ref_end = hit.r_end;
                        candidate.query_beg = hit.q_beg;
                        candidate.query_end = hit.q_end;
                        candidate.cigar = std::move(hit.cigar);
                        candidate.ops = CIGAR::parse(candidate.cigar);
                        candidate.matches = CIGAR::match_len(candidate.cigar);
                        candidate.block_len = std::max(
                            candidate.ref_end - candidate.ref_beg,
                            candidate.query_end - candidate.query_beg
                        );
                        if (candidate.reverse) {
                            candidate.raw_query_beg = query_len - candidate.query_end;
                            candidate.raw_query_end = query_len - candidate.query_beg;
                        } else {
                            candidate.raw_query_beg = candidate.query_beg;
                            candidate.raw_query_end = candidate.query_end;
                        }
                    }
                    if ((candidate.ref_beg < selected.ref_beg) != (candidate.query_beg < selected.query_beg)) {
                        accepted = false;
                        break;
                    }
                }
                if (!accepted || candidate.ops.empty()) continue;
                group.query_beg = std::min(group.query_beg, candidate.query_beg);
                group.query_end = std::max(group.query_end, candidate.query_end);
                group.raw_query_beg = std::min(group.raw_query_beg, candidate.raw_query_beg);
                group.raw_query_end = std::max(group.raw_query_end, candidate.raw_query_end);
                group.ref_beg = std::min(group.ref_beg, candidate.ref_beg);
                group.ref_end = std::max(group.ref_end, candidate.ref_end);
                group.matches += candidate.matches;
                group.block_len += candidate.block_len;
                group.hits.push_back(std::move(candidate));
                if (short_query || short_reference) break;
            }
            group.span = support_span;

            if (place_queries) {
                // A standalone ms component may occupy only a small part of its target; cover the query, not the reference.
                std::vector<std::pair<uint32_t, uint32_t>> intervals;
                intervals.reserve(group.hits.size());
                for (const ComponentHit& hit : group.hits) {
                    intervals.emplace_back(hit.raw_query_beg, hit.raw_query_end);
                }
                group.span = interval_union_length(std::move(intervals));
                group.strong = query_len > 0 &&
                    static_cast<double>(group.span) / query_len >= MIN_JACCARD_FOR_ALIGN_;
                if (group.strong) pairs.push_back(std::move(group));
                continue;
            }

            const uint64_t shorter = std::min<uint64_t>(query_len, backbones[group.rid].sequence.size());
            const uint64_t query_tolerance = component_end_fraction_ > 0.0 ? std::max<uint64_t>(10'000, query_len * component_end_fraction_) : 0;
            const uint64_t ref_tolerance = component_end_fraction_ > 0.0 ? std::max<uint64_t>(10'000, ref_len * component_end_fraction_) : 0;
            const bool opposite_ends = component_end_fraction_ > 0.0 && (
                (group.query_beg <= query_tolerance && ref_len - group.ref_end <= ref_tolerance) ||
                (query_len - group.query_end <= query_tolerance && group.ref_beg <= ref_tolerance)
            );
            const double coverage = shorter > 0 ? static_cast<double>(group.span) / shorter : 0.0;
            group.strong = short_query || short_reference || coverage >= min_component_coverage_ || opposite_ends;
            if (!group.hits.empty()) pairs.push_back(std::move(group));
        }
    }

    std::sort(pairs.begin(), pairs.end(), [](const PairGroup& a, const PairGroup& b) {
        if (a.strong != b.strong) return a.strong > b.strong;
        if (a.short_component != b.short_component) return a.short_component < b.short_component;
        if (a.short_component && b.short_component) {
            const uint64_t lhs = a.matches * b.block_len;
            const uint64_t rhs = b.matches * a.block_len;
            if (lhs != rhs) return lhs > rhs;
        }
        if (a.span != b.span) return a.span > b.span;
        if (a.matches != b.matches) return a.matches > b.matches;
        return std::tie(a.qid, a.rid) < std::tie(b.qid, b.rid);
    });

    if (place_queries) {
        // Relocate only when every complete placement agrees on one target component.
        std::vector<std::vector<ComponentHit>> selected(backbones.size());
        std::vector<uint32_t> target(backbones.size(), UINT32_MAX);
        std::vector<uint8_t> ambiguous(backbones.size(), 0);
        for (PairGroup& pair : pairs) {
            if (!pair.strong) continue;
            const uint32_t component = backbones[pair.rid].component;
            if (target[pair.qid] == UINT32_MAX) {
                target[pair.qid] = component;
                selected[pair.qid] = std::move(pair.hits);
            } else if (target[pair.qid] != component) {
                ambiguous[pair.qid] = 1;
            }
        }
        for (uint32_t qid = 0; qid < selected.size(); ++qid) {
            if (!ambiguous[qid]) continue;
            warning_stream() << "  ! Ambiguous ms placement: " << backbones[qid].name << '\n';
            selected[qid].clear();
        }
        return selected;
    }

    for (const PairGroup& pair : pairs) {
        if (pair.strong && !pair.short_component) {
            merge_component_groups(component_groups, pair.qid, pair.rid);
        }
    }

    struct RootScore {
        uint64_t span{0};
        uint64_t matches{0};
        uint64_t partner_len{0};
    };
    std::vector<std::unordered_map<uint32_t, RootScore>> short_scores(backbones.size());
    for (const PairGroup& pair : pairs) {
        if (!pair.strong || !pair.short_component) continue;
        const bool short_query = short_contig_len > 0 &&
            backbones[pair.qid].sequence.size() <= short_contig_len;
        const bool short_reference = short_contig_len > 0 &&
            backbones[pair.rid].sequence.size() <= short_contig_len;
        if (short_query) {
            RootScore& score = short_scores[pair.qid][component_root(component_groups, pair.rid)];
            score.span += pair.span;
            score.matches += pair.matches;
            score.partner_len = std::max<uint64_t>(score.partner_len, backbones[pair.rid].sequence.size());
        }
        if (short_reference) {
            RootScore& score = short_scores[pair.rid][component_root(component_groups, pair.qid)];
            score.span += pair.span;
            score.matches += pair.matches;
            score.partner_len = std::max<uint64_t>(score.partner_len, backbones[pair.qid].sequence.size());
        }
    }

    std::vector<uint32_t> short_group(backbones.size(), UINT32_MAX);
    static constexpr uint64_t SIMILAR_SUPPORT_PERCENT = 95;
    for (uint32_t component = 0; component < short_scores.size(); ++component) {
        RootScore best;
        for (const auto& item : short_scores[component]) {
            const uint64_t larger_span = std::max(item.second.span, best.span);
            const uint64_t smaller_span = std::min(item.second.span, best.span);
            const bool similar_span =
                larger_span == 0 || smaller_span * 100 >= larger_span * SIMILAR_SUPPORT_PERCENT;
            const bool better =
                short_group[component] == UINT32_MAX ||
                (!similar_span && item.second.span > best.span) ||
                (similar_span && item.second.partner_len > best.partner_len) ||
                (similar_span && item.second.partner_len == best.partner_len && item.second.span > best.span) ||
                (item.second.span == best.span && item.second.partner_len == best.partner_len &&
                 item.second.matches > best.matches) ||
                (item.second.span == best.span && item.second.partner_len == best.partner_len &&
                 item.second.matches == best.matches && item.first < short_group[component]);
            if (better) {
                short_group[component] = item.first;
                best = item.second;
            }
        }
    }

    for (const PairGroup& pair : pairs) {
        if (!pair.strong || !pair.short_component) continue;
        const bool short_query = short_group[pair.qid] != UINT32_MAX;
        const bool short_reference = short_group[pair.rid] != UINT32_MAX;
        if (short_query && component_root(component_groups, short_group[pair.qid]) != component_root(component_groups, pair.rid)) continue;
        if (short_reference && component_root(component_groups, short_group[pair.rid]) != component_root(component_groups, pair.qid)) continue;
        merge_component_groups(component_groups, pair.qid, pair.rid);
    }

    std::vector<uint32_t> group_sizes(backbones.size(), 0);
    for (uint32_t component = 0; component < backbones.size(); ++component) {
        ++group_sizes[component_root(component_groups, component)];
    }
    std::vector<std::vector<uint32_t>> orphan_targets(backbones.size());
    for (const PairGroup& pair : pairs) {
        orphan_targets[pair.qid].push_back(pair.rid);
        orphan_targets[pair.rid].push_back(pair.qid);
    }
    for (uint32_t component = 0; component < backbones.size(); ++component) {
        const uint32_t root = component_root(component_groups, component);
        if (group_sizes[root] != 1 || orphan_targets[component].empty()) continue;

        uint32_t target_root = UINT32_MAX;
        bool ambiguous = false;
        for (uint32_t target : orphan_targets[component]) {
            const uint32_t candidate_root = component_root(component_groups, target);
            if (candidate_root == root) continue;
            if (target_root == UINT32_MAX) target_root = candidate_root;
            else if (target_root != candidate_root) {
                ambiguous = true;
                break;
            }
        }
        if (!ambiguous && target_root != UINT32_MAX) {
            merge_component_groups(component_groups, target_root, component);
            group_sizes[target_root] += group_sizes[root];
            group_sizes[root] = 0;
        }
    }

    std::vector<std::vector<ComponentHit>> selected(backbones.size());
    ScopedIntervalIndex query_occupied, ref_occupied;
    for (PairGroup& pair : pairs) {
        if (component_root(component_groups, pair.qid) != component_root(component_groups, pair.rid)) continue;
        const uint64_t pair_key = (static_cast<uint64_t>(pair.qid) << 32) | pair.rid;
        for (ComponentHit& hit : pair.hits) {
            if (!interval_free(query_occupied, pair_key, hit.raw_query_beg, hit.raw_query_end) || !interval_free(ref_occupied, pair_key, hit.ref_beg, hit.ref_end)) continue;
            add_interval(query_occupied, pair_key, hit.raw_query_beg, hit.raw_query_end);
            add_interval(ref_occupied, pair_key, hit.ref_beg, hit.ref_end);
            selected[pair.qid].push_back(std::move(hit));
        }
    }

    for (std::vector<ComponentHit>& query_hits : selected) {
        std::sort(query_hits.begin(), query_hits.end(), [](const ComponentHit& a, const ComponentHit& b) {
            return std::tie(a.query_beg, a.query_end, a.rid) < std::tie(b.query_beg, b.query_end, b.rid);
        });
    }
    return selected;
}

uint32_t GfaCtgCollapser::best_reference_sample_(
    const std::vector<ComponentBackbone>& backbones,
    uint64_t& reference_n50
) {
    uint32_t sample_count = 0;
    for (const ComponentBackbone& backbone : backbones) {
        sample_count = std::max(sample_count, backbone.sample + 1);
    }

    std::vector<std::vector<uint64_t>> lengths(sample_count);
    std::vector<uint64_t> totals(sample_count, 0);
    for (const ComponentBackbone& backbone : backbones) {
        const uint64_t length = backbone.sequence.size();
        lengths[backbone.sample].push_back(length);
        totals[backbone.sample] += length;
    }

    uint32_t best_sample = 0;
    uint64_t best_n50 = 0;
    for (uint32_t sample = 0; sample < sample_count; ++sample) {
        std::sort(lengths[sample].begin(), lengths[sample].end(), std::greater<uint64_t>());
        const uint64_t half = (totals[sample] + 1) / 2;
        uint64_t accumulated = 0;
        uint64_t n50 = 0;
        for (uint64_t length : lengths[sample]) {
            accumulated += length;
            if (accumulated >= half) {
                n50 = length;
                break;
            }
        }
        if (n50 > best_n50 || (n50 == best_n50 && totals[sample] > totals[best_sample])) {
            best_sample = sample;
            best_n50 = n50;
        }
    }
    reference_n50 = best_n50;
    return best_sample;
}

std::vector<GfaCtgCollapser::BubbleAlignment> GfaCtgCollapser::align_component_backbones_(
    const std::vector<ComponentBackbone>& backbones,
    const std::vector<std::string>& sample_names,
    const std::string& prefix,
    uint32_t short_contig_len,
    const std::unordered_map<std::string, uint32_t>* query_components,
    std::unordered_set<std::string>* placed_queries
) {
    if (backbones.empty()) return {};
    if (placed_queries) placed_queries->clear();

    mm_idxopt_t index_options;
    mm_mapopt_t map_options;
    mm_set_opt(nullptr, &index_options, &map_options);
    if (mm_set_opt(MM2_PRESET_.c_str(), &index_options, &map_options) < 0) {
        error_stream() << "Invalid mm2 preset: " << MM2_PRESET_ << "\n";
        std::exit(1);
    }
    index_options.k = static_cast<short>(chainOpts_.k);
    index_options.w = static_cast<short>(chainOpts_.w);
    map_options.flag |= MM_F_CIGAR | MM_F_EQX;
    map_options.best_n = static_cast<short>(anchorOpts_.max_kept);
    map_options.zdrop = extendOpts_.dyn_zdrop;

    SAVE paf(prefix + ".paf");
    SAVE paf_all(prefix + ".all.paf");
    std::vector<BubbleAlignment> alignments;

    std::vector<uint32_t> component_groups(backbones.size());
    for (uint32_t i = 0; i < component_groups.size(); ++i) component_groups[i] = i;
    std::vector<std::vector<ComponentHit>> all_hits(backbones.size());
    uint64_t component_hits = 0;
    log_stream() << (query_components ? "Aligning misassembly contigs to component backbones ...\n" : "Aligning sample component backbones ...\n");

    struct MapResult {
        std::vector<ComponentHit> hits;
    };

    std::vector<const ComponentBackbone*> references;
    std::vector<const ComponentBackbone*> queries;
    if (query_components) {
        for (const ComponentBackbone& backbone : backbones) {
            if (query_components->contains(backbone.name)) queries.push_back(&backbone);
            else references.push_back(&backbone);
        }
        log_stream() << "  - Component backbones: " << references.size() << '\n';
        log_stream() << "  - Misassembly queries: " << queries.size() << '\n';
    } else {
        uint32_t max_sample = 0;
        for (const ComponentBackbone& backbone : backbones) max_sample = std::max(max_sample, backbone.sample);
        std::vector<std::vector<const ComponentBackbone*>> samples(max_sample + 1);
        for (const ComponentBackbone& backbone : backbones) samples[backbone.sample].push_back(&backbone);

        uint64_t reference_n50 = 0;
        const uint32_t reference_sample = best_reference_sample_(backbones, reference_n50);
        references = samples[reference_sample];
        queries.reserve(backbones.size() - references.size());
        for (uint32_t sample = 0; sample <= max_sample; ++sample) {
            if (sample != reference_sample) {
                queries.insert(queries.end(), samples[sample].begin(), samples[sample].end());
            }
        }
        std::ostringstream n50_text;
        n50_text << std::fixed << std::setprecision(2) << static_cast<double>(reference_n50) / 1'000'000.0;
        const std::string reference_name = reference_sample < sample_names.size() ? sample_names[reference_sample] : "sample" + std::to_string(reference_sample + 1);
        log_stream() << "  - Reference sample: " << reference_name << " (component N50: " << n50_text.str() << " Mb)\n";
    }
    if (references.empty() || queries.empty()) {
        log_stream() << "  - Component alignments: 0\n  - Node alignments: 0\n\n";
        return {};
    }

    std::vector<const char*> ref_sequences, ref_names;
    ref_sequences.reserve(references.size());
    ref_names.reserve(references.size());
    for (const ComponentBackbone* ref : references) {
        ref_sequences.push_back(ref->sequence.c_str());
        ref_names.push_back(ref->name.c_str());
    }
    const int hpc = (index_options.flag & MM_I_HPC) ? 1 : 0;
    mm_idx_t* index = mm_idx_str(
        index_options.w, index_options.k, hpc, index_options.bucket_bits,
        static_cast<int>(references.size()), ref_sequences.data(), ref_names.data()
    );
    if (!index) {
        error_stream() << "Failed to build minimap2 index for component backbones\n";
        std::exit(1);
    }
    mm_mapopt_t round_options = map_options;
    mm_mapopt_update(&round_options, index);

    ThreadPool pool(alignOpts_.threads);
    struct PendingMap {
        const ComponentBackbone* query = nullptr;
        std::future<MapResult> future;
    };
    std::deque<PendingMap> futures;
    const size_t max_pending = std::max<size_t>(1, static_cast<size_t>(alignOpts_.threads) * 2);
    ProgressTracker progress(queries.size());

    for (const ComponentBackbone* query : queries) {
        futures.push_back({query, pool.submit([this, index, round_options, &references, query, query_components, &progress]() mutable {
                MapResult result;
                mm_tbuf_t* buffer = mm_tbuf_init();
                int count = 0;
                mm_reg1_t* regs = mm_map(
                    index, static_cast<int>(query->sequence.size()), query->sequence.c_str(),
                    &count, buffer, &round_options, query->name.c_str()
                );
                result.hits.reserve(count);
                for (int i = 0; i < count; ++i) {
                    const mm_reg1_t& hit = regs[i];
                    if (!hit.p || hit.rid < 0 || static_cast<size_t>(hit.rid) >= references.size()) continue;
                    if (hit.mapq < MIN_MAPQ_ || hit.parent != hit.id) continue;  // filter based on mapping quality
                    const mm_extra_t* extra = static_cast<const mm_extra_t*>(hit.p);
                    if (!extra || extra->n_cigar == 0) continue;

                    ComponentHit candidate;
                    candidate.rid = static_cast<uint32_t>(hit.rid);
                    candidate.ref_beg = static_cast<uint32_t>(hit.rs);
                    candidate.ref_end = static_cast<uint32_t>(hit.re);
                    candidate.raw_query_beg = static_cast<uint32_t>(hit.qs);
                    candidate.raw_query_end = static_cast<uint32_t>(hit.qe);
                    candidate.reverse = hit.rev;
                    candidate.query_beg = hit.rev ? static_cast<uint32_t>(query->sequence.size() - hit.qe) : candidate.raw_query_beg;
                    candidate.query_end = hit.rev ? static_cast<uint32_t>(query->sequence.size() - hit.qs) : candidate.raw_query_end;
                    candidate.matches = static_cast<uint32_t>(hit.mlen);
                    candidate.block_len = static_cast<uint32_t>(hit.blen);
                    candidate.mapq = hit.mapq;
                    candidate.cigar = CIGAR::to_string(extra->cigar, extra->n_cigar);
                    candidate.ops = CIGAR::parse(candidate.cigar);
                    if (candidate.ops.empty() || CIGAR::match_ratio(candidate.cigar) < MIN_MATCH_RATIO_) continue;  // filter based on match ratio

                    const ComponentBackbone& ref = *references[candidate.rid];
                    if (query_components && ref.component == query->component) continue;
                    const uint32_t span = std::max(
                        candidate.ref_end - candidate.ref_beg,
                        candidate.query_end - candidate.query_beg
                    );
                    const uint64_t shorter = std::min(ref.sequence.size(), query->sequence.size());
                    if (shorter > 0 && static_cast<double>(span) / shorter < MIN_ALI_RATIO_) continue;  // filter based on alignment ratio
                    candidate.rid = ref.id;
                    result.hits.push_back(std::move(candidate));
                }
                if (regs) {
                    for (int i = 0; i < count; ++i) std::free(regs[i].p);
                    std::free(regs);
                }
                mm_tbuf_destroy(buffer);
                progress.hit();
                return result;
        })});

        if (futures.size() >= max_pending) {
            PendingMap pending = std::move(futures.front());
            futures.pop_front();
            MapResult result = pending.future.get();
            std::vector<ComponentHit>& query_hits = all_hits[pending.query->id];
            query_hits.insert(
                query_hits.end(),
                std::make_move_iterator(result.hits.begin()),
                std::make_move_iterator(result.hits.end())
            );
        }
    }
    while (!futures.empty()) {
        PendingMap pending = std::move(futures.front());
        futures.pop_front();
        MapResult result = pending.future.get();
        std::vector<ComponentHit>& query_hits = all_hits[pending.query->id];
        query_hits.insert(
            query_hits.end(),
            std::make_move_iterator(result.hits.begin()),
            std::make_move_iterator(result.hits.end())
        );
    }
    pool.stop();
    mm_idx_destroy(index);
    progress.finish();

    for (uint32_t qid = 0; qid < all_hits.size(); ++qid) {
        for (const ComponentHit& hit : all_hits[qid]) {
            const ComponentBackbone& query = backbones[qid];
            const ComponentBackbone& ref = backbones[hit.rid];
            std::ostringstream line;
            line << query.name << '\t' << query.sequence.size() << '\t'
                 << hit.raw_query_beg << '\t' << hit.raw_query_end << '\t'
                 << (hit.reverse ? '-' : '+') << '\t' << ref.name << '\t' << ref.sequence.size() << '\t'
                 << hit.ref_beg << '\t' << hit.ref_end << '\t' << hit.matches << '\t' << hit.block_len << '\t'
                 << static_cast<unsigned>(hit.mapq) << "\ttp:A:P\tcg:Z:" << hit.cigar << '\n';
            paf_all.save(line.str());
        }
    }

    std::vector<std::vector<ComponentHit>> selected = select_component_hits_(
        std::move(all_hits), backbones, short_contig_len,
        component_groups, query_components != nullptr
    );
    for (uint32_t qid = 0; qid < selected.size(); ++qid) {
        const ComponentBackbone& query = backbones[qid];
        for (ComponentHit& hit : selected[qid]) {
            const ComponentBackbone& ref = backbones[hit.rid];
            std::vector<BubbleAlignment> split = split_component_alignment_(
                ref, query, hit.reverse, hit.ref_beg, hit.query_beg, hit.mapq, hit.ops
            );
            const bool usable = !split.empty();
            alignments.insert(
                alignments.end(), std::make_move_iterator(split.begin()), std::make_move_iterator(split.end())
            );
            if (usable && placed_queries) placed_queries->insert(query.name);
            ++component_hits;
            std::ostringstream line;
            line << query.name << '\t' << query.sequence.size() << '\t'
                 << hit.raw_query_beg << '\t' << hit.raw_query_end << '\t'
                 << (hit.reverse ? '-' : '+') << '\t' << ref.name << '\t' << ref.sequence.size() << '\t'
                 << hit.ref_beg << '\t' << hit.ref_end << '\t' << hit.matches << '\t' << hit.block_len << '\t'
                 << static_cast<unsigned>(hit.mapq) << "\ttp:A:P\tcg:Z:" << hit.cigar << '\n';
            paf.save(line.str());
        }
    }

    log_stream() << "  - Component alignments: " << component_hits << "\n";
    log_stream() << "  - Node alignments: " << alignments.size() << "\n\n";
    return alignments;
}

void GfaCtgCollapser::apply_component_alignments_(const std::string& prefix) {
    initialize_cuts_();
    dedup_aligns_();
    const auto groups = build_align_groups_();
    build_propagate_prune_cuts_(groups);
    build_rulemap_(groups);
    SegReplace::Expander expander = build_SegReplace_();
    GfaPathCycleGuard(*this).filter_path_rules(expander);
    expander.save_map(prefix + ".collapse.map");
    expand_and_rewire_edges_(expander);
    rewrite_paths_(expander);
    remove_unused_nodes_();
    append_component_paths_();
}

std::unordered_set<std::string> GfaCtgCollapser::collapse_misassemblies(
    const std::unordered_map<std::string, uint32_t>& source_components,
    const std::string& prefix
) {
    bool have_component_paths = false;
    for (const GfaPath& path : paths_) {
        if (is_global_component_path(path.name)) {
            have_component_paths = true;
            break;
        }
    }
    if (!have_component_paths) rebuild_component_paths();

    const std::vector<uint32_t>& component_ids = get_nodes_connectivity_index();
    std::unordered_map<uint32_t, std::string> component_paths;
    for (const GfaPath& path : paths_) {
        if (!is_global_component_path(path.name) || path.segments.empty()) continue;
        const uint32_t sid = static_cast<uint32_t>(path.segments.front().node_id);
        if (sid < component_ids.size()) component_paths[component_ids[sid]] = path.name;
    }

    std::unordered_map<std::string, uint32_t> query_components;
    std::unordered_map<std::string, std::vector<std::string>> query_members;
    for (const auto& source : source_components) {
        const auto path = std::find_if(paths_.begin(), paths_.end(), [&](const GfaPath& item) {
            return item.name == source.first;
        });
        if (path == paths_.end() || path->segments.empty()) continue;
        const uint32_t sid = static_cast<uint32_t>(path->segments.front().node_id);
        if (sid >= component_ids.size()) continue;
        const auto component = component_paths.find(component_ids[sid]);
        if (component == component_paths.end()) continue;
        query_components[component->second] = source.second;
        query_members[component->second].push_back(source.first);
    }

    std::unordered_set<std::string> placements;
    const std::vector<ComponentBackbone> backbones = component_backbones_({}, &query_components);
    bubble_aligns_ = align_component_backbones_(
        backbones, {}, prefix, 0, &query_components, &placements
    );

    // Existing component paths were mapping references; rebuild them from rewritten sample paths.
    paths_.erase(
        std::remove_if(paths_.begin(), paths_.end(), [](const GfaPath& path) {
            return is_component_path(path.name);
        }),
        paths_.end()
    );
    if (bubble_aligns_.empty()) {
        append_component_paths_();
        return {};
    }
    apply_component_alignments_(prefix);

    std::unordered_set<std::string> integrated;
    for (const std::string& component : placements) {
        for (const std::string& name : query_members[component]) {
            const auto path = std::find_if(paths_.begin(), paths_.end(), [&](const GfaPath& item) {
                return item.name == name;
            });
            if (path == paths_.end()) continue;
            for (const PathSegment& segment : path->segments) {
                std::vector<std::string> roots;
                gfaName::collect_roots_from_name(nodes_[segment.node_id].name, roots);
                if (std::any_of(roots.begin(), roots.end(), [&](const std::string& root) {
                    return root != name;
                })) {
                    integrated.insert(name);
                    break;
                }
            }
        }
    }
    log_stream() << "  - Misassembly contigs integrated: " << integrated.size() << "\n\n";
    return integrated;
}

void GfaCtgCollapser::collapse_backbone_samples_(
    const std::vector<std::string>& files,
    const std::vector<std::string>& sample_names,
    const std::string& prefix,
    uint32_t short_contig_len,
    const std::string& command_line
) {
    log_stream() << "Loading sample graphs for cross-sample component collapse ...\n\n";
    load_from_GFA(files);
    const std::vector<ComponentBackbone> backbones = component_backbones_(sample_names);
    log_stream() << "Samples: " << sample_names.size() << "\n";
    log_stream() << "Component backbones: " << backbones.size() << "\n\n";
    std::string file_prefix = prefix + ".iter0";
    bubble_aligns_ = align_component_backbones_(backbones, sample_names, file_prefix, short_contig_len);

    if (bubble_aligns_.empty()) {
        warning_stream() << "  ! No valid cross-sample component alignments; writing the combined graph unchanged\n";
        append_component_paths_();
        save_to_disk(file_prefix + ".collapse.gfa", true, false, true, command_line);
        save_to_disk(file_prefix + ".collapse.noseq.gfa", true, false, false, command_line);
        return;
    }

    apply_component_alignments_(file_prefix);
    save_to_disk(file_prefix + ".collapse.gfa", true, false, true, command_line);
    save_to_disk(file_prefix + ".collapse.noseq.gfa", true, false, false, command_line);
}


/* ================================================================================================================
 *                                      CTG COLLAPSE DEBUG START
 * ================================================================================================================ */
void CtgCollapseDebugger::inputs(
    const std::vector<std::string>& hap1_files,
    const std::vector<std::string>& hap2_files,
    const std::vector<std::string>& sample_names
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Input pairs: " << hap1_files.size() << '\n';
    for (size_t i = 0; i < hap1_files.size(); ++i) {
        debug_stream() << "  - sample=" << sample_names[i] << " c1=" << hap1_files[i] << " c2=" << hap2_files[i] << '\n';
    }
    debug_stream() << "\n";
}

void CtgCollapseDebugger::anchor_decision(
    const GfaCtgCollapser& collapser,
    const GfaCtgCollapser::Anchor& anchor,
    const char* decision
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Anchor " << decision
                   << ": read=" << (anchor.read_name ? *anchor.read_name : "?")
                   << " c1=" << collapser.nodes_[anchor.a_sid].name << ':' << anchor.a_beg << '-' << anchor.a_end
                   << " c2=" << collapser.nodes_[anchor.b_sid].name << ':' << anchor.b_beg << '-' << anchor.b_end
                   << " orientation=" << (anchor.reverse ? "reverse" : "forward") << '\n';
}

void CtgCollapseDebugger::blocks(
    const GfaCtgCollapser& collapser,
    const std::vector<GfaCtgCollapser::AnchorBlock>& blocks,
    const char* stage
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Anchor blocks (" << stage << "): " << blocks.size() << '\n';
    for (const auto& block : blocks) {
        debug_stream() << "  - c1=" << collapser.nodes_[block.a_sid].name << ':' << block.a_beg << '-' << block.a_end
                       << " c2=" << collapser.nodes_[block.b_sid].name << ':' << block.b_beg << '-' << block.b_end
                       << " orientation=" << (block.reverse ? "reverse" : "forward")
                       << " anchors=" << block.anchors << '\n';
    }
    debug_stream() << "\n";
}

void CtgCollapseDebugger::pair_score(
    const GfaCtgCollapser& collapser,
    uint32_t a_sid,
    uint32_t b_sid,
    uint64_t anchors,
    uint64_t span,
    bool selected
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Contig pair: c1=" << collapser.nodes_[a_sid].name
                   << " c2=" << collapser.nodes_[b_sid].name
                   << " anchors=" << anchors
                   << " span=" << span
                   << " selected_by_either_side=" << (selected ? "yes" : "no") << '\n';
}

void CtgCollapseDebugger::block_decision(
    const GfaCtgCollapser& collapser,
    const GfaCtgCollapser::AnchorBlock& block,
    const char* decision
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Block " << decision
                   << ": c1=" << collapser.nodes_[block.a_sid].name << ':' << block.a_beg << '-' << block.a_end
                   << " c2=" << collapser.nodes_[block.b_sid].name << ':' << block.b_beg << '-' << block.b_end
                   << " orientation=" << (block.reverse ? "reverse" : "forward")
                   << " anchors=" << block.anchors << '\n';
}

void CtgCollapseDebugger::alignment(
    const GfaCtgCollapser& collapser,
    const GfaCtgCollapser::AnchorBlock& block,
    size_t alignment_count
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Block alignment: c1=" << collapser.nodes_[block.a_sid].name
                   << ':' << block.a_beg << '-' << block.a_end
                   << " c2=" << collapser.nodes_[block.b_sid].name
                   << ':' << block.b_beg << '-' << block.b_end
                   << " orientation=" << (block.reverse ? "reverse" : "forward")
                   << " results=" << alignment_count << '\n';
}

void CtgCollapseDebugger::paths(const GfaCtgCollapser& collapser) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "GFA paths: " << collapser.paths_.size() << '\n';
    for (const GfaPath& path : collapser.paths_) {
        std::ostringstream out;
        out << "  - " << path.name << ':';
        for (const PathSegment& segment : path.segments) {
            out << ' ' << collapser.nodes_[segment.node_id].name << (segment.is_reverse ? '-' : '+');
        }
        debug_stream() << out.str() << '\n';
    }
    debug_stream() << "\n";
}
/* ================================================================================================================
 *                                       CTG COLLAPSE DEBUG END
 * ================================================================================================================ */
