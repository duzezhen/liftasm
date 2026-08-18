#include "gfa_gapfill.hpp"
#include "gfa_bubble.hpp"
#include "gfa_gapfill_logger.hpp"

#include "../include/ThreadPool.hpp"
#include "logger.hpp"
#include "progress_tracker.hpp"
#include "save.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <numeric>
#include <queue>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

extern "C" {
#include "minimap.h"
}

namespace {

// Detect cycles while selecting globally compatible gap-fill connections.
class DisjointSet {
public:
    explicit DisjointSet(size_t n) : parent_(n), rank_(n, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    uint32_t find(uint32_t x) {
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }

    bool join(uint32_t a, uint32_t b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (rank_[a] < rank_[b]) std::swap(a, b);
        parent_[b] = a;
        if (rank_[a] == rank_[b]) ++rank_[a];
        return true;
    }

private:
    std::vector<uint32_t> parent_;
    std::vector<uint8_t> rank_;
};

// Track haplotype intervals already assigned on each bridge contig.
using PhaseIntervals = std::map<uint32_t, uint32_t>;

struct LocalPhaseAssignments {
    std::array<PhaseIntervals, 2> hap;
};

// Format paired left/right metrics in the final TSV report.
std::string metric_pair(double left, double right) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2);
    if (left < 0.0) out << '.';
    else out << left;
    out << '/';
    if (right < 0.0) out << '.';
    else out << right;
    return out.str();
}

// Detect and merge local phase intervals during connection selection.
bool overlaps_phase_interval(
    const PhaseIntervals& intervals,
    uint32_t begin,
    uint32_t end
) {
    auto next = intervals.upper_bound(end);
    if (next == intervals.begin()) return false;
    return std::prev(next)->second >= begin;
}

void add_phase_interval(PhaseIntervals& intervals, uint32_t begin, uint32_t end) {
    auto next = intervals.lower_bound(begin);
    if (next != intervals.begin()) {
        auto previous = std::prev(next);
        if (previous->second >= begin) {
            begin = previous->first;
            end = std::max(end, previous->second);
            next = intervals.erase(previous);
        }
    }
    while (next != intervals.end() && next->first <= end) {
        end = std::max(end, next->second);
        next = intervals.erase(next);
    }
    intervals.emplace(begin, end);
}

// Normalize sample and contig names used in output files and GFA segments.
std::string safe_name(std::string name) {
    for (char& c : name) {
        const bool valid = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '.' || c == '_' || c == '-';
        if (!valid) c = '_';
    }
    return name.empty() ? "sample" : name;
}

std::string numbered_contig_name(uint8_t hap, const char* kind, size_t number) {
    std::ostringstream name;
    name << 'h' << static_cast<uint32_t>(hap) << kind << std::setw(6) << std::setfill('0') << number << 'l';
    return name.str();
}

std::string compact_contig_name(std::string name) {
    size_t marker = name.find(".hap1.");
    if (marker == std::string::npos) marker = name.find(".hap2.");
    if (marker != std::string::npos) name.replace(marker, 6, ".");
    return name;
}

std::string output_contig_name(const std::string& name, const std::string& sample) {
    const std::string compact = compact_contig_name(name);
    const std::string prefix = sample + '.';
    return compact.compare(0, prefix.size(), prefix) == 0 ?
        compact.substr(prefix.size()) : compact;
}

// Shared nodes of a contig pair
struct PlotRelation {
    uint32_t a{0}, b{0};   // contig indexes
    uint64_t shared_bp{0}; // total shared base pairs
    std::vector<std::pair<uint64_t, uint64_t>> anchors;  // shared graph-node offsets on contig a and b
};

int64_t layout_median(std::vector<int64_t> values) {
    if (values.empty()) return 0;
    const size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    return values[mid];
}

uint64_t layout_deviation(const std::vector<int64_t>& values, int64_t center) {
    uint64_t sum = 0;
    for (int64_t value : values) sum += static_cast<uint64_t>(std::llabs(value - center));
    return sum;
}

}  // namespace


// Configure and run the complete gap-fill pipeline.
GfaGapfill::GfaGapfill(Params params) : params_(params) {}

void GfaGapfill::set_alignment_options(
    const opt::ChainOpts& chain,
    const opt::AnchorOpts& anchor,
    const opt::ExtendOpts& extend,
    const opt::AlignOpts& align,
    const std::string& preset,
    uint32_t max_occ,
    double min_match,
    double min_ali_ratio,
    uint8_t min_mapq,
    double boundary_coverage
) {
    mm2_.k = chain.k;
    mm2_.w = chain.w;
    mm2_.best_n = anchor.max_kept;
    mm2_.zdrop = extend.dyn_zdrop;
    mm2_.preset = preset;
    mm2_.max_occ = max_occ;
    mm2_.min_match = min_match;
    mm2_.min_ali_ratio = min_ali_ratio;
    mm2_.min_mapq = min_mapq;
    mm2_.boundary_coverage = boundary_coverage;
    params_.threads = std::max<uint32_t>(1, align.threads);
}


// ================================================= Clear state =================================================
void GfaGapfill::clear_state_() {
    prepared_candidates_.clear();
    prepared_selected_.clear();
    candidates_prepared_ = false;
    fragments_.clear();
    sample_chains_.clear();
    records_.clear();
    relocations_.clear();
    relocations_.reserve(misassemblies_.size());
    for (const MisassemblyContig& ms : misassemblies_) {
        RelocationRecord record;
        record.sample = ms.sample;
        record.hap = ms.hap;
        record.source = ms.source;
        record.side = ms.side;
        record.source_component = ms.source_component;
        record.status = relocated_misassemblies_.count(ms.name) ? "relocated" : "unmapped";
        record.relocated_contig = ms.name;
        relocations_.push_back(std::move(record));
    }
}

size_t GfaGapfill::gapfill(const std::vector<GfaBubble::Bubble>& bubbles) {
    log_stream() << "Gap-filling sample contigs with cross-sample bridging paths ...\n\n";
    std::vector<Candidate> candidates, selected;
    if (candidates_prepared_) {
        candidates = std::move(prepared_candidates_);
        selected = std::move(prepared_selected_);
        candidates_prepared_ = false;
    } else {
        clear_state_();
        index_bubble_nodes_(bubbles);
        build_fragments_();
        build_graph_order_();
        candidates = build_candidates_();
        selected = select_candidates_(candidates);
        refine_boundaries_(selected, candidates);
    }

    mark_used_relocations_(selected);
    build_sample_chains_(selected);
    mark_redundant_fragments_(selected);

    write_html_(candidates);
    print_summary_(selected.size());

    return selected.size();
}

std::vector<GfaGapfill::MisassemblyContig> GfaGapfill::prepare_misassemblies(
    const std::vector<GfaBubble::Bubble>& bubbles
) {
    log_stream() << "Detecting terminal misassemblies before gap filling ...\n\n";
    clear_state_();
    index_bubble_nodes_(bubbles);
    build_fragments_();
    build_graph_order_();

    std::vector<Candidate> candidates = build_candidates_(true);
    check_unplaced_sequence_(candidates);
    std::vector<MisassemblyContig> result = split_misassemblies_(candidates);
    if (result.empty()) {
        size_t kept = 0;
        for (size_t i = 0; i < candidates.size(); ++i) {
            const Candidate& candidate = candidates[i];
            const uint64_t fill_bp = candidate.right_boundary.bridge_cut_bp -
                candidate.left_boundary.bridge_cut_bp;
            if ((!candidate.target_overlaps && candidate.target_distance > params_.max_gap_bp) ||
                fill_bp > params_.max_gap_bp) continue;
            if (kept != i) candidates[kept] = std::move(candidates[i]);
            ++kept;
        }
        candidates.resize(kept);
        std::vector<Candidate> selected = select_candidates_(candidates);
        refine_boundaries_(selected, candidates);
        prepared_candidates_ = std::move(candidates);
        prepared_selected_ = std::move(selected);
        candidates_prepared_ = true;
    }
    return result;
}

void GfaGapfill::set_misassemblies(
    const std::vector<MisassemblyContig>& records,
    const std::unordered_set<std::string>& relocated
) {
    misassemblies_ = records;
    relocated_misassemblies_ = relocated;
    misassembly_index_.clear();
    misassembly_index_.reserve(records.size());
    for (uint32_t i = 0; i < records.size(); ++i) {
        misassembly_index_.emplace(records[i].name, i);
    }
}


// ================================================= Contig preparation and graph layout =================================================
void GfaGapfill::index_bubble_nodes_(const std::vector<GfaBubble::Bubble>& bubbles) {
    log_stream() << "Indexing bubble nodes used as phase evidence ...\n";
    bubble_nodes_.assign(getNumNodes(), 0);

    std::vector<uint8_t> boundary(getNumNodes(), 0);
    for (const GfaBubble::Bubble& bubble : bubbles) {
        const uint32_t source = Vertex::get_segment_id(bubble.get_source());
        const uint32_t sink = Vertex::get_segment_id(bubble.get_sink());
        if (source < boundary.size()) boundary[source] = 1;
        if (sink < boundary.size()) boundary[sink] = 1;
    }

    std::vector<uint32_t> path_count(getNumNodes(), 0);
    std::vector<uint32_t> last_path(getNumNodes(), 0);
    std::vector<uint32_t> touched;
    uint32_t path_id = 1;

    for (const GfaBubble::Bubble& bubble : bubbles) {
        const auto& paths = bubble.get_paths();
        if (paths.size() < 2) continue;
        const uint32_t source = Vertex::get_segment_id(bubble.get_source());
        const uint32_t sink = Vertex::get_segment_id(bubble.get_sink());
        bool eligible = true;
        touched.clear();
        for (const std::vector<uint32_t>& path : paths) {
            uint64_t path_bp = 0;
            for (uint32_t vertex : path) {
                const uint32_t segment = Vertex::get_segment_id(vertex);
                if (segment >= bubble_nodes_.size() || segment == source || segment == sink || last_path[segment] == path_id) continue;
                last_path[segment] = path_id;
                path_bp += nodes_[segment].length;
                if (!boundary[segment] && path_count[segment]++ == 0) {
                    touched.push_back(segment);
                }
            }
            if (path_bp > params_.phase_path_len) eligible = false;
            ++path_id;
        }
        for (uint32_t segment : touched) {
            if (eligible && path_count[segment] < paths.size()) {
                bubble_nodes_[segment] = 1;
            }
            path_count[segment] = 0;
        }
    }

    uint64_t bubble_bp = 0;
    size_t bubble_node_count = 0;
    for (uint32_t segment = 0; segment < bubble_nodes_.size(); ++segment) {
        if (!bubble_nodes_[segment]) continue;
        ++bubble_node_count;
        bubble_bp += nodes_[segment].length;
    }
    log_stream() << "  - Bubbles indexed: " << bubbles.size() << '\n';
    log_stream() << "  - Bubble-interior nodes: " << bubble_node_count << '\n';
    log_stream() << "  - Bubble-interior length: " << bubble_bp << " bp\n\n";
}

void GfaGapfill::build_fragments_() {
    log_stream() << "Indexing sample contig paths ...\n";
    fragments_.reserve(paths_.size());

    // ------------------------------------------------ Read path names ------------------------------------------------
    for (size_t path_id = 0; path_id < paths_.size(); ++path_id) {
        GfaPath& path = paths_[path_id];
        Fragment fragment;
        if (!parse_path_name_(path.name, fragment.sample, fragment.hap)) {
            GfaGapfillDebugger::path(path.name, "drop:unrecognized-sample");
            continue;
        }
        // ------------------------------------------------ Normalize and check paths ------------------------------------------------
        path.name = compact_contig_name(path.name);
        if (path.segments.empty()) {
            GfaGapfillDebugger::path(path.name, "drop:empty");
            continue;
        }

        // ------------------------------------------------ Validate path nodes ------------------------------------------------
        fragment.path_id = path_id;
        std::unordered_set<uint32_t> seen;
        fragment.vertices.reserve(path.segments.size());
        bool valid = true;
        for (const PathSegment& segment : path.segments) {
            const uint32_t sid = static_cast<uint32_t>(segment.node_id);
            if (sid >= nodes_.size() || nodes_[sid].deleted) {
                GfaGapfillDebugger::path(path.name, "drop:invalid-segment");
                valid = false;
                break;
            }
            if (!seen.insert(sid).second) {
                GfaGapfillDebugger::path(path.name, "drop:repeated-segment node=" + nodes_[sid].name);
                valid = false;
                break;
            }

            const uint32_t component = connectivity_index_[sid];
            if (fragment.component == UINT32_MAX) fragment.component = component;
            else if (fragment.component != component) {
                GfaGapfillDebugger::path(path.name, "drop:multiple-components");
                valid = false;
                break;
            }
            fragment.vertices.push_back(Vertex::make_vertex(sid, segment.is_reverse));
            fragment.length += nodes_[sid].length;
        }
        if (!valid) continue;

        // ------------------------------------------------ Choose path orientation ------------------------------------------------
        uint32_t forward_drops = 0;
        for (size_t i = 1; i < fragment.vertices.size(); ++i) {
            if (vertex_topo_rank(Vertex(fragment.vertices[i - 1])) > vertex_topo_rank(Vertex(fragment.vertices[i]))) ++forward_drops;
        }
        std::vector<uint32_t> reverse(fragment.vertices.rbegin(), fragment.vertices.rend());
        for (uint32_t& vertex : reverse) vertex ^= 1u;
        uint32_t reverse_drops = 0;
        for (size_t i = 1; i < reverse.size(); ++i) {
            if (vertex_topo_rank(Vertex(reverse[i - 1])) > vertex_topo_rank(Vertex(reverse[i]))) ++reverse_drops;
        }

        const uint32_t forward_start = vertex_topo_rank(Vertex(fragment.vertices.front()));
        const uint32_t reverse_start = vertex_topo_rank(Vertex(reverse.front()));
        if (reverse_drops < forward_drops || (reverse_drops == forward_drops && reverse_start < forward_start)) {
            fragment.vertices = std::move(reverse);
            fragment.reverse = !fragment.reverse;
        }

        // ------------------------------------------------ Build fragment indexes ------------------------------------------------
        rebuild_fragment_index_(fragment);
        fragment.eligible = fragment.length >= params_.min_contig_bp;

        // ------------------------------------------------ Attach relocation records ------------------------------------------------
        const auto relocation = misassembly_index_.find(path.name);
        if (relocation != misassembly_index_.end()) {
            fragment.eligible = true;
            fragment.relocation = static_cast<int32_t>(relocation->second);
            RelocationRecord& record = relocations_[relocation->second];
            if (record.status == "relocated") {
                record.target_component = fragment.component;
                record.target_backbone = "component" + std::to_string(fragment.component);
            }
        }

        // ------------------------------------------------ Save fragment ------------------------------------------------
        GfaGapfillDebugger::path(
            path.name, fragment.eligible ? "keep" : "exclude:short",
            fragment.sample, fragment.hap, fragment.component,
            fragment.vertices.size(), fragment.length,
            nodes_[Vertex::get_segment_id(fragment.vertices.front())].name,
            Vertex::get_is_reverse(fragment.vertices.front()),
            nodes_[Vertex::get_segment_id(fragment.vertices.back())].name,
            Vertex::get_is_reverse(fragment.vertices.back())
        );
        fragments_.push_back(std::move(fragment));
    }
    log_stream() << "  - Contig paths indexed: " << fragments_.size() << "\n\n";
}

void GfaGapfill::build_graph_order_() {
    log_stream() << "Ordering sample contigs from shared graph nodes ...\n";

    // ------------------------------------------------ Collect layout contigs ------------------------------------------------
    std::vector<LayoutContig> layout;
    layout.reserve(fragments_.size());
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible) continue;
        LayoutContig contig;
        contig.id = i;
        contig.component = fragment.component;
        contig.hap = fragment.hap;
        contig.bp = fragment.length;
        contig.name = paths_[fragment.path_id].name;
        contig.sample = fragment.sample;
        // ------------------------------------------------ Build shared-node anchors ------------------------------------------------
        contig.anchors.reserve(fragment.vertices.size());
        for (uint32_t j = 0; j < fragment.vertices.size(); ++j) {
            const uint32_t segment = Vertex::get_segment_id(fragment.vertices[j]);
            contig.anchors.push_back({
                segment, nodes_[segment].length,
                fragment.path_bp[j] + nodes_[segment].length / 2
            });
        }
        layout.push_back(std::move(contig));
    }

    // ------------------------------------------------ Place contigs on graph coordinates ------------------------------------------------
    place_contigs_(layout);

    // ------------------------------------------------ Apply layout and orientation ------------------------------------------------
    for (const LayoutContig& contig : layout) {
        Fragment& fragment = fragments_[contig.id];
        fragment.layout_start = contig.start;
        if (!contig.reverse) continue;
        std::reverse(fragment.vertices.begin(), fragment.vertices.end());
        for (uint32_t& vertex : fragment.vertices) vertex ^= 1u;
        fragment.reverse = !fragment.reverse;
        rebuild_fragment_index_(fragment);
    }

    // ------------------------------------------------ Sort and report fragments ------------------------------------------------
    std::sort(fragments_.begin(), fragments_.end(), [](const Fragment& a, const Fragment& b) {
        if (a.sample != b.sample) return a.sample < b.sample;
        if (a.component != b.component) return a.component < b.component;
        if (a.layout_start != b.layout_start) return a.layout_start < b.layout_start;
        return a.path_id < b.path_id;
    });
    if (DEBUG_ENABLED) {
        std::map<std::pair<std::string, uint32_t>, size_t> groups;
        for (const Fragment& fragment : fragments_) ++groups[{fragment.sample, fragment.component}];
        for (const auto& group : groups) {
            GfaGapfillDebugger::component(group.first.first, group.first.second, group.second);
        }
    }
    log_stream() << "  - Contigs ordered from shared graph nodes: " << layout.size() << "\n\n";
}

bool GfaGapfill::parse_path_name_(const std::string& name, std::string& sample, uint8_t& hap) {
    const size_t h1 = name.find(".hap1.");
    const size_t h2 = name.find(".hap2.");
    size_t marker = h1 != std::string::npos ? h1 : h2;
    if (marker != std::string::npos) {
        if (marker == 0) return false;
        sample = name.substr(0, marker);
        hap = h1 != std::string::npos ? 1 : 2;
        return true;
    }

    const size_t t1 = name.find(".h1tg");
    const size_t t2 = name.find(".h2tg");
    const size_t m1 = name.find(".h1ms");
    const size_t m2 = name.find(".h2ms");
    marker = t1 != std::string::npos ? t1 : (t2 != std::string::npos ? t2 : (m1 != std::string::npos ? m1 : m2));
    if (marker == std::string::npos || marker == 0) return false;
    sample = name.substr(0, marker);
    hap = t1 != std::string::npos || m1 != std::string::npos ? 1 : 2;
    return true;
}

void GfaGapfill::rebuild_fragment_index_(Fragment& fragment) const {
    fragment.vertex_index.clear();
    fragment.vertex_index.reserve(fragment.vertices.size());
    fragment.bubble_positions.clear();
    fragment.rank_drops.assign(fragment.vertices.size(), 0);
    fragment.path_bp.assign(fragment.vertices.size() + 1, 0);
    fragment.bubble_bp.assign(fragment.vertices.size() + 1, 0);
    for (uint32_t i = 0; i < fragment.vertices.size(); ++i) {
        const uint32_t vertex = fragment.vertices[i];
        const uint32_t segment = Vertex::get_segment_id(vertex);
        const uint32_t rank = vertex_topo_rank(Vertex(vertex));
        fragment.vertex_index.emplace_back(vertex, i);
        fragment.path_bp[i + 1] = fragment.path_bp[i] + nodes_[segment].length;
        fragment.bubble_bp[i + 1] = fragment.bubble_bp[i] + (segment < bubble_nodes_.size() && bubble_nodes_[segment] ? nodes_[segment].length : 0);
        if (segment < bubble_nodes_.size() && bubble_nodes_[segment]) {
            fragment.bubble_positions.push_back(i);
        }
        if (i > 0) {
            fragment.rank_drops[i] = fragment.rank_drops[i - 1];
            if (vertex_topo_rank(Vertex(fragment.vertices[i - 1])) > rank) {
                ++fragment.rank_drops[i];
            }
        }
    }
    std::sort(fragment.vertex_index.begin(), fragment.vertex_index.end());
    fragment.length = fragment.path_bp.back();
}

void GfaGapfill::place_contigs_(std::vector<LayoutContig>& contigs) {
    // ------------------------------------------------ Index shared nodes ------------------------------------------------
    struct Occurrence {
        uint32_t contig;  // contig index
        uint64_t offset;  // node offset on contig
        uint32_t length;  // node length
    };

    //  Input: contig -> node
    // Output:   node -> contig
    std::unordered_map<uint64_t, std::vector<Occurrence>> occurrences;
    for (uint32_t i = 0; i < contigs.size(); ++i) {
        for (const auto& anchor : contigs[i].anchors) {
            const uint64_t key = (static_cast<uint64_t>(contigs[i].component) << 32) | anchor.segment;
            occurrences[key].push_back({i, anchor.offset, anchor.length});
        }
    }

    // ------------------------------------------------ Build contig relations ------------------------------------------------
    std::unordered_map<uint64_t, PlotRelation> relation_map;
    for (const auto& item : occurrences) {
        const auto& found = item.second;
        for (size_t i = 0; i < found.size(); ++i) {
            for (size_t j = i + 1; j < found.size(); ++j) {
                if (found[i].contig == found[j].contig) continue;
                uint32_t a = found[i].contig, b = found[j].contig;
                uint64_t ao = found[i].offset, bo = found[j].offset;
                if (a > b) { std::swap(a, b); std::swap(ao, bo); }
                const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;  // key = (contig_a << 32) | contig_b
                PlotRelation& relation = relation_map[key];
                relation.a = a;
                relation.b = b;
                relation.shared_bp += std::min(found[i].length, found[j].length);
                relation.anchors.emplace_back(ao, bo);
            }
        }
    }

    // ------------------------------------------------ Build adjacency list ------------------------------------------------
    std::vector<PlotRelation> relations;
    relations.reserve(relation_map.size());
    for (auto& item : relation_map) relations.push_back(std::move(item.second));
    std::vector<std::vector<uint32_t>> adjacency(contigs.size());
    for (uint32_t i = 0; i < relations.size(); ++i) {
        adjacency[relations[i].a].push_back(i);
        adjacency[relations[i].b].push_back(i);
    }

    // ------------------------------------------------ Place each graph component ------------------------------------------------
    std::map<uint32_t, std::vector<uint32_t>> components;
    for (uint32_t i = 0; i < contigs.size(); ++i) components[contigs[i].component].push_back(i);
    std::vector<uint8_t> placed(contigs.size(), 0);

    for (const auto& component : components) {
        int64_t cursor = 0;
        size_t remaining = component.second.size();
        while (remaining > 0) {
            // ------------------------------------------------ Start a new layout chain ------------------------------------------------
            uint32_t root = UINT32_MAX;
            for (uint32_t index : component.second) {
                if (!placed[index] && (root == UINT32_MAX || contigs[index].bp > contigs[root].bp)) root = index;
            }
            contigs[root].start = cursor;
            contigs[root].reverse = false;
            placed[root] = 1;
            --remaining;

            using QueueItem = std::tuple<uint64_t, uint32_t, uint32_t, uint32_t>;  // (shared_bp, relation_id, from, to)
            // ------------------------------------------------ Expand by strongest shared nodes ------------------------------------------------
            std::priority_queue<QueueItem> queue;
            for (uint32_t relation : adjacency[root]) {
                const PlotRelation& edge = relations[relation];
                queue.emplace(edge.shared_bp, relation, root, edge.a == root ? edge.b : edge.a);
            }

            while (!queue.empty()) {
                const auto [weight, relation_id, from, to] = queue.top();
                queue.pop();
                (void)weight;
                if (!placed[from] || placed[to]) continue;
                const PlotRelation& relation = relations[relation_id];
                std::vector<int64_t> forward, reverse;
                forward.reserve(relation.anchors.size());
                reverse.reserve(relation.anchors.size());
                for (const auto& anchor : relation.anchors) {
                    const uint64_t from_offset = relation.a == from ? anchor.first : anchor.second;
                    const uint64_t to_offset = relation.a == to ? anchor.first : anchor.second;
                    const LayoutContig& source = contigs[from];
                    const int64_t global = source.start + static_cast<int64_t>(
                        source.reverse ? source.bp - from_offset : from_offset
                    );
                    forward.push_back(global - static_cast<int64_t>(to_offset));
                    reverse.push_back(global - static_cast<int64_t>(contigs[to].bp - to_offset));
                }
                // ------------------------------------------------ Choose orientation and position ------------------------------------------------
                const int64_t forward_start = layout_median(forward);
                const int64_t reverse_start = layout_median(reverse);
                contigs[to].reverse = layout_deviation(reverse, reverse_start) < layout_deviation(forward, forward_start);
                contigs[to].start = contigs[to].reverse ? reverse_start : forward_start;
                placed[to] = 1;
                --remaining;
                for (uint32_t next_relation : adjacency[to]) {
                    const PlotRelation& edge = relations[next_relation];
                    const uint32_t next = edge.a == to ? edge.b : edge.a;
                    if (!placed[next]) queue.emplace(edge.shared_bp, next_relation, to, next);
                }
            }

            int64_t end = cursor;
            for (uint32_t index : component.second) {
                if (placed[index]) end = std::max(end, contigs[index].start + static_cast<int64_t>(contigs[index].bp));
            }
            // Leave space for unconnected contigs.
            cursor = end + std::max<int64_t>(1'000'000, end / 100);
        }

        // ------------------------------------------------ Normalize component start ------------------------------------------------
        int64_t begin = contigs[component.second.front()].start;
        for (uint32_t index : component.second) begin = std::min(begin, contigs[index].start);
        for (uint32_t index : component.second) contigs[index].start -= begin;
    }
}


// ================================================= Shared-node candidate evidence =================================================
std::vector<GfaGapfill::Candidate> GfaGapfill::build_candidates_(
    bool include_long_gaps
) const {
    // ------------------------------------------------ Group contigs by sample ------------------------------------------------
    std::vector<std::string> group_keys(fragments_.size());
    std::map<std::string, std::vector<uint32_t>> fragments_by_sample;
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible || fragment.redundant) continue;
        fragments_by_sample[fragment.sample].push_back(i);
        group_keys[i] = fragment.sample + '\t' + std::to_string(fragment.component);
    }

    // ------------------------------------------------ Collect node positions ------------------------------------------------
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> node_paths;  // (segment_id, contig_id, position_in_fragment)
    size_t path_nodes = 0;
    for (const Fragment& fragment : fragments_) {
        if (!fragment.eligible) continue;
        path_nodes += fragment.vertices.size();
    }
    node_paths.reserve(path_nodes);
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        if (!fragments_[i].eligible) continue;
        for (uint32_t pos = 0; pos < fragments_[i].vertices.size(); ++pos) {
            node_paths.emplace_back(
                Vertex::get_segment_id(fragments_[i].vertices[pos]), i, pos
            );
        }
    }
    std::sort(node_paths.begin(), node_paths.end());

    // ------------------------------------------------ Build overlap index ------------------------------------------------
    OverlapIndex overlaps(fragments_.size());
    for (size_t begin = 0; begin < node_paths.size();) {
        size_t end = begin + 1;
        while (end < node_paths.size() && std::get<0>(node_paths[end]) == std::get<0>(node_paths[begin])) {
            ++end;
        }
        const uint32_t sid = std::get<0>(node_paths[begin]);
        for (size_t i = begin; i < end; ++i) {
            for (size_t j = i + 1; j < end; ++j) {
                const uint32_t a = std::get<1>(node_paths[i]);
                const uint32_t b = std::get<1>(node_paths[j]);
                const uint32_t a_pos = std::get<2>(node_paths[i]);
                const uint32_t b_pos = std::get<2>(node_paths[j]);
                PathOverlap& ab = overlaps[a][b];
                PathOverlap& ba = overlaps[b][a];
                ab.bp += nodes_[sid].length;
                ba.bp += nodes_[sid].length;
                ++ab.nodes;
                ++ba.nodes;
                ab.a_begin = std::min(ab.a_begin, a_pos);
                ab.a_end = std::max(ab.a_end, a_pos);
                ab.b_begin = std::min(ab.b_begin, b_pos);
                ab.b_end = std::max(ab.b_end, b_pos);
                ba.a_begin = std::min(ba.a_begin, b_pos);
                ba.a_end = std::max(ba.a_end, b_pos);
                ba.b_begin = std::min(ba.b_begin, a_pos);
                ba.b_end = std::max(ba.b_end, a_pos);
            }
        }
        begin = end;
    }
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>().swap(node_paths);

    // ------------------------------------------------ Dump overlap details ------------------------------------------------
    if (DEBUG_ENABLED) {
        const GfaGapfillBoundary boundary(*this);
        for (uint32_t a = 0; a < overlaps.size(); ++a) {
            for (const auto& item : overlaps[a]) {
                const uint32_t b = item.first;
                if (a >= b) continue;
                const OverlapSupport support = overlap_support_(fragments_[a], fragments_[b], item.second);
                const PathOverlap& overlap = item.second;
                const Boundary left_boundary = boundary.inspect(
                    a, overlap.a_begin, b, overlap.b_begin,
                    false, true
                );
                const Boundary right_boundary = boundary.inspect(
                    a, overlap.a_end, b, overlap.b_end,
                    true, true
                );
                const GfaGapfillDebugger::BoundarySupport left_support_debug{
                    left_boundary.minimizer,
                    left_boundary.phase.similarity,
                    left_boundary.phase.target_bp,
                    left_boundary.phase.bridge_bp
                };
                const GfaGapfillDebugger::BoundarySupport right_support_debug{
                    right_boundary.minimizer,
                    right_boundary.phase.similarity,
                    right_boundary.phase.target_bp,
                    right_boundary.phase.bridge_bp
                };
                GfaGapfillDebugger::overlap(
                    paths_[fragments_[a].path_id].name,
                    paths_[fragments_[b].path_id].name,
                    support.shared_bp, support.overlap_bp,
                    support.similarity,
                    fragments_[a].length == 0 ? 0.0 : static_cast<double>(support.shared_bp) / fragments_[a].length,
                    fragments_[b].length == 0 ? 0.0 : static_cast<double>(support.shared_bp) / fragments_[b].length,
                    left_support_debug, right_support_debug
                );
            }
        }
    }

    // ------------------------------------------------ Build sample groups ------------------------------------------------
    std::vector<std::vector<uint32_t>> sample_groups;
    sample_groups.reserve(fragments_by_sample.size());
    for (auto& item : fragments_by_sample) sample_groups.push_back(std::move(item.second));

    // ------------------------------------------------ Search samples in parallel ------------------------------------------------
    log_stream() << "Building gap-fill candidates by sample ...\n";
    ThreadPool pool(params_.threads);
    std::vector<std::future<CandidateBatch>> futures;
    futures.reserve(sample_groups.size());
    ProgressTracker progress(sample_groups.size());

    for (size_t sample_id = 0; sample_id < sample_groups.size(); ++sample_id) {
        futures.emplace_back(pool.submit([&, sample_id] {
            CandidateBatch batch = build_sample_candidates_(
                sample_groups[sample_id], overlaps, group_keys, include_long_gaps
            );
            progress.hit();
            return batch;
        }));
    }

    // ------------------------------------------------ Merge sample results ------------------------------------------------
    std::vector<Candidate> candidates;
    std::vector<CandidateBatch> results(futures.size());
    for (size_t i = 0; i < futures.size(); ++i) {
        results[i] = futures[i].get();
    }
    progress.finish();
    pool.stop();

    for (size_t i = 0; i < results.size(); ++i) {
        CandidateBatch& batch = results[i];
        const std::string& sample = fragments_[sample_groups[i].front()].sample;
        GfaGapfillDebugger::search(
            sample, sample_groups[i].size(), batch.tested, batch.low_support,
            batch.wrong_group, batch.anchor_failed, batch.candidates.size()
        );
        candidates.insert(
            candidates.end(),
            std::make_move_iterator(batch.candidates.begin()),
            std::make_move_iterator(batch.candidates.end())
        );
    }

    log_stream() << "  - Candidate connections: " << candidates.size() << "\n\n";
    return candidates;
}

GfaGapfill::CandidateBatch GfaGapfill::build_sample_candidates_(
    const std::vector<uint32_t>& sample_fragments,
    const OverlapIndex& overlaps,
    const std::vector<std::string>& group_keys,
    bool include_long_gaps
) const {
    // ------------------------------------------------ Prepare sample search ------------------------------------------------
    CandidateBatch batch;
    const GfaGapfillBoundary boundary(*this);

    // ------------------------------------------------ Try each left contig ------------------------------------------------
    for (uint32_t left_id : sample_fragments) {
        const Fragment& left = fragments_[left_id];
        if (left.length < params_.min_contig_bp) continue;

        // ------------------------------------------------ Find bridge and right contigs ------------------------------------------------
        // Keep all bridge paths for each right contig.
        /*
         S.left -> T.bridge1 -> S.right
         S.left -> T.bridge2 -> S.right
         S.left -> U.bridge3 -> S.right
        */
        std::unordered_map<uint32_t, std::vector<Candidate>> evidence_by_right;
        // Enumerate bridge paths that connect this left contig to a later contig.
        for (const auto& bridge_item : overlaps[left_id]) {
            const uint32_t bridge_id = bridge_item.first;
            const Fragment& bridge = fragments_[bridge_id];
            if (bridge.sample == left.sample) continue;

            // The left side must share enough nodes with the bridge.
            const OverlapSupport left_support = overlap_support_(
                left, bridge, bridge_item.second
            );
            if (left_support.shared_bp < params_.min_overlap_bp) {
                ++batch.low_support;
                continue;
            }

            // Find the right contig through this bridge.
            for (const auto& right_item : overlaps[bridge_id]) {
                const uint32_t right_id = right_item.first;
                const Fragment& right = fragments_[right_id];
                ++batch.tested;

                // The right side must also have enough shared nodes.
                const OverlapSupport right_support = overlap_support_(
                    bridge, right, right_item.second
                );
                if (right_support.shared_bp < params_.min_overlap_bp) {
                    ++batch.low_support;
                    continue;
                }
                // Left and right must belong to the same target path group.
                if (right_id == left_id || group_keys[right_id] != group_keys[left_id] || !ordered_pair_(left, right, include_long_gaps)) {
                    ++batch.wrong_group;
                    continue;
                }

                // The bridge must cover the left-to-right layout range.
                const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
                const int64_t bridge_end = bridge.layout_start + static_cast<int64_t>(bridge.length);
                if (bridge.layout_start > left_end || bridge_end < right.layout_start) {
                    ++batch.wrong_group;
                    continue;
                }

                // ------------------------------------------------ Build one candidate ------------------------------------------------
                // Save this possible left -> bridge -> right connection.
                Candidate candidate;
                candidate.left = left_id;
                candidate.right = right_id;
                candidate.bridge = bridge_id;
                const int64_t target_gap = right.layout_start - left_end;
                candidate.target_overlaps = target_gap < 0;
                candidate.target_distance = static_cast<uint64_t>(target_gap < 0 ? -target_gap : target_gap);
                candidate.left_shared = left_support.shared_bp;
                candidate.right_shared = right_support.shared_bp;
                // Get graph boundaries before the later mm2 refinement.
                if (!boundary.prepare(left_id, right_id, bridge_id, candidate.left_boundary, candidate.right_boundary)) {
                    ++batch.anchor_failed;
                    continue;
                }
                const uint64_t fill_bp = candidate.right_boundary.bridge_cut_bp -
                    candidate.left_boundary.bridge_cut_bp;
                if (!include_long_gaps && fill_bp > params_.max_gap_bp) {
                    ++batch.wrong_group;
                    continue;
                }
                evidence_by_right[right_id].push_back(candidate);
            }
        }

        // ------------------------------------------------ Compare bridge haplotypes ------------------------------------------------
        for (auto& item : evidence_by_right) {
            // check if the left and right contigs are spanned by a homologous contig in this sample
            const bool homolog_span = homolog_spans_(left_id, item.first, overlaps);

            // Keep the strongest spanning contig for each bridge haplotype.
            std::map<std::string, std::array<size_t, 2>> best_phase;
            for (size_t i = 0; i < item.second.size(); ++i) {
                const Candidate& candidate = item.second[i];
                const Fragment& bridge = fragments_[candidate.bridge];
                if (bridge.hap == 0 || bridge.hap > 2) continue;
                auto inserted = best_phase.emplace(
                    bridge.sample, std::array<size_t, 2>{{SIZE_MAX, SIZE_MAX}}
                );
                size_t& best = inserted.first->second[bridge.hap - 1];
                if (best == SIZE_MAX || raw_phase_score_(candidate) > raw_phase_score_(item.second[best])) {
                    best = i;
                }
            }

            for (Candidate& candidate : item.second) {
                const Fragment& bridge = fragments_[candidate.bridge];
                candidate.phase_score = 0.5;
                if (bridge.hap > 0 && bridge.hap <= 2) {
                    const size_t alternative_id = best_phase.at(bridge.sample)[2 - bridge.hap];
                    if (alternative_id != SIZE_MAX) {
                        const Candidate& alternative = item.second[alternative_id];
                        const double phase[2] = {
                            candidate.left_boundary.phase.similarity,
                            candidate.right_boundary.phase.similarity
                        };
                        const double other[2] = {
                            alternative.left_boundary.phase.similarity,
                            alternative.right_boundary.phase.similarity
                        };
                        double relative = 1.0;
                        uint32_t compared = 0;
                        for (uint32_t side = 0; side < 2; ++side) {
                            if (phase[side] < 0.0 || other[side] < 0.0) continue;
                            const double total = phase[side] + other[side];
                            relative *= total > 0.0 ? phase[side] / total : 0.5;
                            ++compared;
                        }
                        if (compared > 0) {
                            candidate.phase_score = std::pow(relative, 1.0 / compared);
                        }
                    }
                }
            }

            // ------------------------------------------------ Add sample support ------------------------------------------------
            update_sample_support_(item.second);

            // ------------------------------------------------ Keep the best bridge ------------------------------------------------
            // Boundary alignment is reserved for these coarse candidates instead of every bridge.
            size_t best = 0;
            for (Candidate& candidate : item.second) {
                candidate.homolog_span = homolog_span;
            }
            for (size_t i = 1; i < item.second.size(); ++i) {
                if (candidate_better_(item.second[i], item.second[best])) best = i;
            }
            batch.candidates.push_back(std::move(item.second[best]));
        }
    }
    return batch;
}

bool GfaGapfill::ordered_pair_(
    const Fragment& left,
    const Fragment& right,
    bool include_long_gaps
) const {
    if (left.sample != right.sample || left.component != right.component || left.length < params_.min_contig_bp || right.length < params_.min_contig_bp) return false;
    if (left.layout_start > right.layout_start) return false;
    const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
    const int64_t right_end = right.layout_start + static_cast<int64_t>(right.length);
    if (right.layout_start >= left_end) {
        return include_long_gaps || static_cast<uint64_t>(right.layout_start - left_end) <= params_.max_gap_bp;
    }
    const uint64_t overlap = static_cast<uint64_t>(std::min(left_end, right_end) - right.layout_start);
    const double fraction = static_cast<double>(overlap) / std::min(left.length, right.length);
    return fraction <= params_.max_target_overlap;
}

GfaGapfill::OverlapSupport GfaGapfill::overlap_support_(
    const Fragment& a,
    const Fragment& b,
    const PathOverlap& overlap
) const {
    OverlapSupport support;
    if (overlap.nodes == 0) return support;

    const auto path_region = [](const Fragment& fragment, uint32_t beg, uint32_t end) {
        if (beg > end || end >= fragment.vertices.size()) return std::pair<uint32_t, uint64_t>{0, 0};
        return std::pair<uint32_t, uint64_t>{
            end - beg + 1,
            fragment.path_bp[end + 1] - fragment.path_bp[beg]
        };
    };

    const auto a_region = path_region(a, overlap.a_begin, overlap.a_end);
    const auto b_region = path_region(b, overlap.b_begin, overlap.b_end);
    const uint32_t total_nodes = a_region.first + b_region.first;
    if (total_nodes < overlap.nodes) return support;

    support.shared_bp = overlap.bp;
    support.overlap_bp = std::min(a_region.second, b_region.second);
    support.similarity = support.overlap_bp == 0 ? 0.0 : std::min(1.0, static_cast<double>(support.shared_bp) / support.overlap_bp);
    return support;
}

bool GfaGapfill::homolog_spans_(
    uint32_t left_id,
    uint32_t right_id,
    const OverlapIndex& overlaps
) const {
    const Fragment& left = fragments_[left_id];
    const Fragment& right = fragments_[right_id];
    const GfaGapfillBoundary boundary(*this);
    for (const auto& item : overlaps[left_id]) {
        const uint32_t bridge_id = item.first;
        const Fragment& bridge = fragments_[bridge_id];
        if (bridge.sample != left.sample || bridge.hap == left.hap || bridge.component != left.component) continue;

        const auto right_it = overlaps[bridge_id].find(right_id);
        if (right_it == overlaps[bridge_id].end()) continue;
        const OverlapSupport left_support = overlap_support_(left, bridge, item.second);
        const OverlapSupport right_support = overlap_support_(bridge, right, right_it->second);
        if (left_support.shared_bp < params_.min_overlap_bp || right_support.shared_bp < params_.min_overlap_bp || left_support.similarity < params_.min_similarity || right_support.similarity < params_.min_similarity) continue;

        const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
        const int64_t bridge_end = bridge.layout_start + static_cast<int64_t>(bridge.length);
        if (bridge.layout_start > left_end || bridge_end < right.layout_start) continue;
        if (boundary.spans(left_id, right_id, bridge_id)) return true;
    }
    return false;
}

double GfaGapfill::raw_phase_score_(const Candidate& candidate) {
    const double left = candidate.left_boundary.phase.similarity;
    const double right = candidate.right_boundary.phase.similarity;
    if (left >= 0.0 && right >= 0.0) return std::sqrt(left * right);
    return -1.0;
}

void GfaGapfill::update_sample_support_(std::vector<Candidate>& evidence) const {
    // Take the part of each bridge that actually fills the gap.
    std::vector<std::vector<uint32_t>> fill_paths;
    fill_paths.reserve(evidence.size());
    for (size_t i = 0; i < evidence.size(); ++i) {
        const Candidate& witness = evidence[i];
        const Fragment& bridge = fragments_[witness.bridge];
        const uint32_t begin = witness.left_boundary.bridge_pos + 1;
        const uint32_t end = witness.right_boundary.bridge_pos;
        fill_paths.emplace_back(
            bridge.vertices.begin() + begin,
            bridge.vertices.begin() + end
        );
    }

    // Group similar bridge paths together.
    const GfaBubble::PathClusterer clusterer(
        *this, params_.path_difference, mm2_.k, mm2_.w
    );
    const std::vector<uint32_t> labels = clusterer.cluster(
        fill_paths, GfaBubble::Type::Tip
    );

    // Count each sample once per distinct bridge allele. Homozygous
    // haplotypes therefore contribute one vote, heterozygous ones one each.
    std::map<std::string, std::vector<size_t>> by_sample;
    for (size_t i = 0; i < evidence.size(); ++i) {
        by_sample[fragments_[evidence[i].bridge].sample].push_back(i);
    }
    std::unordered_map<uint32_t, uint32_t> support_by_label;
    for (const auto& item : by_sample) {
        std::unordered_set<uint32_t> sample_labels;
        for (size_t index : item.second) sample_labels.insert(labels[index]);
        for (uint32_t label : sample_labels) ++support_by_label[label];
    }

    // Write support rate and update the final candidate confidence.
    for (size_t i = 0; i < evidence.size(); ++i) {
        Candidate& candidate = evidence[i];
        candidate.sample_support = support_by_label[labels[i]];
        candidate.spanning_samples = static_cast<uint32_t>(by_sample.size());
        candidate.sample_score = candidate.spanning_samples == 0 ? 0.0 : static_cast<double>(candidate.sample_support) / candidate.spanning_samples;
        update_confidence_(candidate);
    }
}


// ================================================= Candidate selection =================================================
std::vector<GfaGapfill::Candidate> GfaGapfill::select_candidates_(
    std::vector<Candidate>& candidates
) {
    log_stream() << "Selecting unambiguous gap-fill connections ...\n";

    // ------------------------------------------------ Sort candidates ------------------------------------------------
    std::sort(candidates.begin(), candidates.end(), connection_better_);  // Put the most useful connections first.

    // ------------------------------------------------ Prepare selection state ------------------------------------------------
    std::vector<Candidate> selected;
    // Track used contig ends, bridge phase intervals and connected components.
    std::vector<int32_t> incoming(fragments_.size(), -1), outgoing(fragments_.size(), -1);
    std::vector<std::unordered_map<std::string, LocalPhaseAssignments>> bridge_phase(fragments_.size());
    DisjointSet sets(fragments_.size());

    // ------------------------------------------------ Check each candidate ------------------------------------------------
    for (Candidate& candidate : candidates) {
        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        const std::string& left_name = paths_[fragments_[candidate.left].path_id].name;
        const std::string& right_name = paths_[fragments_[candidate.right].path_id].name;
        const std::string& bridge_name = paths_[fragments_[candidate.bridge].path_id].name;

        // ------------------------------------------------ Check overlap and confidence ------------------------------------------------
        if (candidate.target_overlaps && static_cast<double>(candidate.target_distance) / std::min(left.length, right.length) > params_.max_target_overlap) {
            candidate.status = Candidate::COORDINATE_CONFLICT;
            GfaGapfillDebugger::selection(
                left_name, right_name, bridge_name,
                "drop:long-overlap", candidate.confidence
            );
            continue;
        }
        if (candidate.confidence < params_.min_confidence && !candidate.homolog_span) {
            candidate.status = Candidate::LOW_CONFIDENCE;
            GfaGapfillDebugger::selection(
                left_name, right_name, bridge_name,
                "drop:low-confidence", candidate.confidence
            );
            continue;
        }
        if (outgoing[candidate.left] >= 0 || incoming[candidate.right] >= 0) {
            candidate.status = Candidate::USED_END;
            GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "drop:used-end", candidate.confidence);
            continue;
        }

        // ------------------------------------------------ Check phase conflict ------------------------------------------------
        const uint8_t candidate_hap = left.hap;
        LocalPhaseAssignments& local_phase = bridge_phase[candidate.bridge][left.sample];
        // Do not reuse the same bridge interval for the opposite haplotype.
        if (candidate_hap != 0 && overlaps_phase_interval(
                local_phase.hap[2 - candidate_hap],
                candidate.left_boundary.bridge_pos,
                candidate.right_boundary.bridge_pos
            )
        ) {
            candidate.status = Candidate::PHASE_CONFLICT;
            GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "drop:phase-conflict", candidate.confidence);
            continue;
        }

        // ------------------------------------------------ Check target coordinates ------------------------------------------------
        uint32_t left_entry = 0;
        if (incoming[candidate.left] >= 0) {
            left_entry = selected[incoming[candidate.left]].right_boundary.target_pos;
        }
        uint32_t right_exit = static_cast<uint32_t>(fragments_[candidate.right].vertices.size() - 1);
        if (outgoing[candidate.right] >= 0) {
            right_exit = selected[outgoing[candidate.right]].left_boundary.target_pos;
        }
        // Keep the new cuts inside the still-unused target ranges.
        if (left_entry > candidate.left_boundary.target_pos || candidate.right_boundary.target_pos > right_exit) {
            candidate.status = Candidate::COORDINATE_CONFLICT;
            GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "drop:coordinate-conflict", candidate.confidence);
            continue;
        }

        // ------------------------------------------------ Check cycle ------------------------------------------------
        if (sets.find(candidate.left) == sets.find(candidate.right)) {
            candidate.status = Candidate::CYCLE;
            GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "drop:cycle", candidate.confidence);
            continue;
        }

        // ------------------------------------------------ Keep candidate ------------------------------------------------
        sets.join(candidate.left, candidate.right);
        const int32_t index = static_cast<int32_t>(selected.size());
        selected.push_back(candidate);
        candidate.status = Candidate::KEPT;
        selected.back().status = Candidate::KEPT;
        outgoing[candidate.left] = index;
        incoming[candidate.right] = index;
        if (candidate_hap != 0) {
            add_phase_interval(
                local_phase.hap[candidate_hap - 1],
                candidate.left_boundary.bridge_pos,
                candidate.right_boundary.bridge_pos
            );
        }
        GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "keep", candidate.confidence);
    }
    log_stream() << "  - Connections selected: " << selected.size() << "\n\n";
    return selected;
}

bool GfaGapfill::candidate_better_(const Candidate& a, const Candidate& b) {
    if (a.confidence != b.confidence) return a.confidence > b.confidence;
    if (a.sample_score != b.sample_score) return a.sample_score > b.sample_score;
    if (a.phase_score != b.phase_score) return a.phase_score > b.phase_score;
    const uint32_t a_phase_boundaries = (a.left_boundary.phase.similarity >= 0.0) + (a.right_boundary.phase.similarity >= 0.0);
    const uint32_t b_phase_boundaries = (b.left_boundary.phase.similarity >= 0.0) + (b.right_boundary.phase.similarity >= 0.0);
    if (a_phase_boundaries != b_phase_boundaries) return a_phase_boundaries > b_phase_boundaries;
    const double a_raw_phase = raw_phase_score_(a);
    const double b_raw_phase = raw_phase_score_(b);
    if (a_raw_phase != b_raw_phase) return a_raw_phase > b_raw_phase;
    const uint64_t as = a.left_shared + a.right_shared;
    const uint64_t bs = b.left_shared + b.right_shared;
    if (as != bs) return as > bs;
    const uint32_t alen = a.right_boundary.bridge_pos - a.left_boundary.bridge_pos;
    const uint32_t blen = b.right_boundary.bridge_pos - b.left_boundary.bridge_pos;
    if (alen != blen) return alen < blen;
    return a.bridge < b.bridge;
}

bool GfaGapfill::connection_better_(const Candidate& a, const Candidate& b) {
    // A complete contig on the other haplotype fixes the target pairing before
    // local gap evidence is used to resolve the remaining ambiguous ends.
    if (a.homolog_span != b.homolog_span) return a.homolog_span;
    if (a.target_overlaps != b.target_overlaps) return !a.target_overlaps;
    if (a.target_distance != b.target_distance) return a.target_distance < b.target_distance;
    return candidate_better_(a, b);
}


// ================================================= Refine left/bridge and bridge/right cut positions with local alignment =================================================
void GfaGapfill::refine_boundaries_(
    std::vector<Candidate>& selected,
    std::vector<Candidate>& candidates
) const {
    log_stream() << "Refining selected gap boundaries from local alignments ...\n";
    const GfaGapfillBoundary boundary(*this);
    ThreadPool pool(params_.threads);
    std::vector<std::future<bool>> futures;
    futures.reserve(selected.size());
    ProgressTracker progress(selected.size());

    // Only selected graph connections need base-level cut refinement.
    for (Candidate& candidate : selected) {
        Candidate* task = &candidate;
        futures.emplace_back(pool.submit([&boundary, task, &progress]() {
            const bool refined = boundary.refine(
                task->left, task->right, task->bridge,
                task->left_boundary, task->right_boundary
            );
            progress.hit();
            return refined;
        }));
    }

    size_t refined = 0;
    // Finish all alignments before final end assignment.
    for (size_t i = 0; i < selected.size(); ++i) {
        if (futures[i].get()) ++refined;
    }

    // Refresh scores and write the cached evidence once.
    for (Candidate& candidate : selected) {
        update_confidence_(candidate);
        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        const Fragment& bridge = fragments_[candidate.bridge];
        GfaGapfillDebugger::refined_boundary(
            {paths_[left.path_id].name, left.length, left.reverse},
            candidate.left_boundary.target_cut_bp,
            {paths_[right.path_id].name, right.length, right.reverse},
            candidate.right_boundary.target_cut_bp,
            {paths_[bridge.path_id].name, bridge.length, bridge.reverse},
            candidate.left_boundary.bridge_cut_bp,
            candidate.right_boundary.bridge_cut_bp,
            {candidate.left_boundary.minimizer,
             candidate.left_boundary.phase.similarity,
             candidate.left_boundary.phase.target_bp,
             candidate.left_boundary.phase.bridge_bp},
            {candidate.right_boundary.minimizer,
             candidate.right_boundary.phase.similarity,
             candidate.right_boundary.phase.target_bp,
             candidate.right_boundary.phase.bridge_bp},
            candidate.left_boundary.identity,
            candidate.right_boundary.identity
        );
        GfaGapfillDebugger::candidate(
            paths_[fragments_[candidate.left].path_id].name,
            paths_[fragments_[candidate.right].path_id].name,
            paths_[fragments_[candidate.bridge].path_id].name,
            fragments_[candidate.left].component,
            candidate.left_boundary.phase.similarity,
            candidate.right_boundary.phase.similarity,
            candidate.left_boundary.phase.target_bp,
            candidate.left_boundary.phase.bridge_bp,
            candidate.right_boundary.phase.target_bp,
            candidate.right_boundary.phase.bridge_bp,
            candidate.left_boundary.minimizer,
            candidate.right_boundary.minimizer,
            candidate.homolog_span,
            candidate.sample_support, candidate.spanning_samples,
            candidate.sample_score, candidate.phase_score,
            candidate.alignment_score, candidate.confidence
        );
    }

    // Keep the HTML/report copy of selected candidates in sync.
    size_t selected_id = 0;
    for (Candidate& candidate : candidates) {
        if (candidate.status == Candidate::KEPT) {
            candidate = selected[selected_id++];
        }
    }
    progress.finish();
    pool.stop();
    log_stream() << "  - Boundaries refined: " << refined << "\n\n";
}

void GfaGapfill::update_confidence_(Candidate& candidate) {
    const double left = candidate.left_boundary.identity >= 0.0 && candidate.left_boundary.coverage >= 0.0 ?
        std::sqrt(candidate.left_boundary.identity * candidate.left_boundary.coverage) :
        (candidate.left_boundary.minimizer >= 0.0 ? candidate.left_boundary.minimizer : 0.5);
    const double right = candidate.right_boundary.identity >= 0.0 && candidate.right_boundary.coverage >= 0.0 ?
        std::sqrt(candidate.right_boundary.identity * candidate.right_boundary.coverage) :
        (candidate.right_boundary.minimizer >= 0.0 ? candidate.right_boundary.minimizer : 0.5);
    candidate.alignment_score = std::min(left, right);
    candidate.confidence =
        std::pow(std::clamp(candidate.sample_score, 0.0, 1.0), 0.5) *
        std::pow(std::clamp(candidate.phase_score, 0.0, 1.0), 0.1) *
        std::pow(std::clamp(candidate.alignment_score, 0.0, 1.0), 0.4);
}


// ================================================= Misassembly detection and relocation =================================================
void GfaGapfill::check_unplaced_sequence_(std::vector<Candidate>& candidates) const {
    log_stream() << "Checking long unplaced contig ends ...\n";

    // ------------------------------------------------ Bit encoding ------------------------------------------------
    // Left/right end
    constexpr uint32_t RIGHT_END_CODE = 0;
    constexpr uint32_t LEFT_END_CODE = 1;
    constexpr uint32_t END_KEY_SHIFT = 1;
    // left/right sequence checks
    constexpr uint8_t CHECK_LEFT = 1u << 0;
    constexpr uint8_t CHECK_RIGHT = 1u << 1;
    // left, then right
    constexpr uint32_t TERMINAL_SLOTS = 2;
    constexpr uint32_t LEFT_TERMINAL_OFFSET = 0;
    constexpr uint32_t RIGHT_TERMINAL_OFFSET = 1;

    // ------------------------------------------------ Helper Functions ------------------------------------------------
    const auto check_candidate = [this](  // Compare the unplaced left/right ends with the bridge sequence.
        Candidate& candidate, bool check_left, bool check_right
    ) {
        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        const Fragment& bridge = fragments_[candidate.bridge];

        if (check_left) {
            const uint64_t target_begin = left.path_bp[candidate.left_boundary.junction_target_pos + 1];
            const uint64_t bridge_begin = bridge.path_bp[candidate.left_boundary.junction_bridge_pos + 1];
            const uint64_t compare_bp = std::min(left.length - target_begin, bridge.length - bridge_begin);
            const std::string unplaced = fragment_subsequence_(left, target_begin, target_begin + compare_bp);
            const std::string filled = fragment_subsequence_(bridge, bridge_begin, bridge_begin + compare_bp);
            candidate.left_unplaced_similarity = mm2_similarity_(filled, unplaced);
            candidate.left_misassembly = candidate.left_unplaced_similarity >= 0.0 && candidate.left_unplaced_similarity < params_.misassembly_similarity;
        }
        if (check_right) {
            const uint64_t target_end = right.path_bp[candidate.right_boundary.junction_target_pos];
            const uint64_t bridge_end = bridge.path_bp[candidate.right_boundary.junction_bridge_pos];
            const uint64_t compare_bp = std::min(target_end, bridge_end);
            const std::string unplaced = fragment_subsequence_(right, target_end - compare_bp, target_end);
            const std::string filled = fragment_subsequence_(bridge, bridge_end - compare_bp, bridge_end);
            candidate.right_unplaced_similarity = mm2_similarity_(filled, unplaced);
            candidate.right_misassembly = candidate.right_unplaced_similarity >= 0.0 && candidate.right_unplaced_similarity < params_.misassembly_similarity;
        }
    };


    // ------------------------------------------------ Keep one candidate per end ------------------------------------------------
    std::unordered_map<uint64_t, uint32_t> best;  // Key = (fragment_id << 1) | side, Value = candidate index; left side = 1, right side = 0
    for (uint32_t i = 0; i < candidates.size(); ++i) {
        Candidate& candidate = candidates[i];
        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        candidate.left_unplaced = left.length - left.path_bp[candidate.left_boundary.junction_target_pos + 1];
        candidate.right_unplaced = right.path_bp[candidate.right_boundary.junction_target_pos];
        const uint64_t keys[2] = {
            (static_cast<uint64_t>(candidate.left) << END_KEY_SHIFT) | LEFT_END_CODE,
            (static_cast<uint64_t>(candidate.right) << END_KEY_SHIFT) | RIGHT_END_CODE
        };
        const uint64_t lengths[2] = {candidate.left_unplaced, candidate.right_unplaced};
        for (uint32_t side = 0; side < 2; ++side) {
            if (lengths[side] <= params_.misassembly_check_bp) continue;
            const auto found = best.find(keys[side]);
            if (found == best.end() || candidate_better_(candidate, candidates[found->second])) {
                best[keys[side]] = i;
            }
        }
    }

    // ------------------------------------------------ Group fragments by haplotype ------------------------------------------------
    std::map<std::pair<std::string, uint8_t>, std::vector<uint32_t>> haplotypes;  // Sample name + haplotype -> fragment IDs
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        haplotypes[{fragments_[i].sample, fragments_[i].hap}].push_back(i);
    }
    if (haplotypes.size() <= 2) {
        log_stream() << "  - Misassembly check skipped: only " << haplotypes.size() << " haplotypes\n\n";
        return;
    }

    // ------------------------------------------------ Find shared terminal nodes ------------------------------------------------
    // rule: <=4 haplotypes -> support 1 is checked;
    //        >4 haplotypes -> support 1/2 is checked.
    const uint32_t max_rare_support = haplotypes.size() <= 4 ? 1 : 2;
    std::unordered_map<uint32_t, std::vector<uint32_t>> terminal_nodes;  // Segment ID -> candidate IDs
    std::vector<uint64_t> terminal_bp(candidates.size() * 2, 0);
    for (const auto& item : best) {
        const uint32_t candidate_id = item.second;
        const bool left_side = (item.first & LEFT_END_CODE) == LEFT_END_CODE;
        const Candidate& candidate = candidates[candidate_id];
        const Fragment& fragment = fragments_[left_side ? candidate.left : candidate.right];
        const uint32_t begin = left_side ? candidate.left_boundary.junction_target_pos + 1 : 0;
        const uint32_t end = left_side ? fragment.vertices.size() : candidate.right_boundary.junction_target_pos;
        const uint32_t terminal_id = candidate_id * TERMINAL_SLOTS + (left_side ? LEFT_TERMINAL_OFFSET : RIGHT_TERMINAL_OFFSET);
        terminal_bp[terminal_id] = left_side ? candidate.left_unplaced : candidate.right_unplaced;
        for (uint32_t pos = begin; pos < end; ++pos) {
            terminal_nodes[Vertex::get_segment_id(fragment.vertices[pos])].push_back(terminal_id);
        }
    }

    // ------------------------------------------------ Count haplotype support ------------------------------------------------
    std::vector<uint32_t> haplotype_support(terminal_bp.size(), 0);
    for (const auto& haplotype : haplotypes) {
        std::unordered_set<uint32_t> seen_nodes;
        std::unordered_map<uint32_t, uint64_t> shared_bp;
        for (uint32_t fragment_id : haplotype.second) {
            for (uint32_t vertex : fragments_[fragment_id].vertices) {
                const uint32_t segment = Vertex::get_segment_id(vertex);
                const auto found = terminal_nodes.find(segment);
                if (found == terminal_nodes.end() || !seen_nodes.insert(segment).second) continue;
                for (uint32_t terminal_id : found->second) {
                    shared_bp[terminal_id] += nodes_[segment].length;
                }
            }
        }
        for (const auto& support : shared_bp) {
            if (support.second >= terminal_bp[support.first] * params_.misassembly_similarity) {
                ++haplotype_support[support.first];
            }
        }
    }

    // ------------------------------------------------ Keep supported ends and check the rest ------------------------------------------------
    std::vector<uint8_t> sides(candidates.size(), 0);  // Bit 0 = unplaced, Bit 1 = left side, Bit 2 = right side, Bit 3 = both sides
    size_t supported = 0, checked = 0;
    for (const auto& item : best) {
        const uint32_t candidate_id = item.second;
        const bool left_side = item.first & 1;
        const uint32_t terminal_id = candidate_id * 2 + !left_side;
        if (haplotype_support[terminal_id] > max_rare_support) {
            const Candidate& candidate = candidates[candidate_id];
            const Fragment& fragment = fragments_[left_side ? candidate.left : candidate.right];
            log_stream() << "  - Preserved " << (left_side ? "right end of " : "left end of ") << paths_[fragment.path_id].name << ": supported by " << haplotype_support[terminal_id] << " haplotypes\n";
            ++supported;
            continue;
        }
        sides[candidate_id] |= left_side ? CHECK_LEFT : CHECK_RIGHT;
        ++checked;
    }
    std::vector<uint32_t> work;
    work.reserve(best.size());
    for (uint32_t i = 0; i < sides.size(); ++i) {
        if (sides[i]) work.push_back(i);
    }

    // ------------------------------------------------ Compare selected ends with bridge sequence ------------------------------------------------
    if (!work.empty()) {
        ThreadPool pool(std::min<size_t>(params_.threads, work.size()));
        std::vector<std::future<void>> futures;
        futures.reserve(work.size());
        for (uint32_t i : work) {
            futures.emplace_back(pool.submit([&, i] {
                check_candidate(
                    candidates[i], sides[i] & CHECK_LEFT, sides[i] & CHECK_RIGHT
                );
            }));
        }
        ProgressTracker progress(work.size());
        for (auto& future : futures) {
            future.get();
            progress.hit();
        }
        progress.finish();
        pool.stop();
    }

    // ------------------------------------------------ Report possible misassemblies ------------------------------------------------
    size_t possible = 0;
    for (uint32_t i : work) {
        const Candidate& candidate = candidates[i];
        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        const Fragment& bridge = fragments_[candidate.bridge];
        if ((sides[i] & CHECK_LEFT) && candidate.left_unplaced_similarity < 0.0) {
            warning_stream() << "  ! Preserving unchecked left end of " << paths_[left.path_id].name << ": sequence unavailable\n";
        }
        if ((sides[i] & CHECK_RIGHT) && candidate.right_unplaced_similarity < 0.0) {
            warning_stream() << "  ! Preserving unchecked right end of " << paths_[right.path_id].name << ": sequence unavailable\n";
        }
        possible += candidate.left_misassembly;
        possible += candidate.right_misassembly;
        if (candidate.left_misassembly) {
            uint64_t begin = left.path_bp[candidate.left_boundary.junction_target_pos + 1];
            uint64_t end = left.length;
            if (left.reverse) {
                end = left.length - begin;
                begin = 0;
            }
            log_stream() << "  - Possible misassembly: " << paths_[left.path_id].name
                         << ':' << begin << '-' << end
                         << " | compared with " << paths_[bridge.path_id].name
                         << " | similarity=" << std::fixed << std::setprecision(4)
                         << candidate.left_unplaced_similarity
            << " | haplotypes=" << haplotype_support[i * TERMINAL_SLOTS + LEFT_TERMINAL_OFFSET]
                         << '/' << haplotypes.size() << '\n';
        }
        if (candidate.right_misassembly) {
            uint64_t begin = 0;
            uint64_t end = right.path_bp[candidate.right_boundary.junction_target_pos];
            if (right.reverse) {
                begin = right.length - end;
                end = right.length;
            }
            log_stream() << "  - Possible misassembly: " << paths_[right.path_id].name
                         << ':' << begin << '-' << end
                         << " | compared with " << paths_[bridge.path_id].name
                         << " | similarity=" << std::fixed << std::setprecision(4)
                         << candidate.right_unplaced_similarity
                         << " | haplotypes=" << haplotype_support[i * TERMINAL_SLOTS + RIGHT_TERMINAL_OFFSET]
                         << '/' << haplotypes.size() << '\n';
        }
    }
    log_stream() << "  - Long ends checked: " << checked << '\n';
    log_stream() << "  - Supported ends preserved: " << supported << '\n';
    log_stream() << "  - Possible misassemblies: " << possible << "\n\n";
}

std::vector<GfaGapfill::MisassemblyContig> GfaGapfill::split_misassemblies_(
    const std::vector<Candidate>& candidates
) {
    // ------------------------------------------------ Helper Functions ------------------------------------------------
    const auto misassembly_contig_name = [](
        const std::string& sample, uint8_t hap, size_t number
    ) {
        std::ostringstream name;
        name << sample << ".h" << static_cast<uint32_t>(hap) << "ms" << std::setw(6) << std::setfill('0') << number << 'l';
        return name.str();
    };

    // ------------------------------------------------ Bit encoding ------------------------------------------------
    enum SplitEnd : uint64_t {
        SPLIT_RIGHT_END = 0,
        SPLIT_LEFT_END = 1
    };
    constexpr uint32_t SPLIT_END_BITS = 1;
    constexpr uint32_t VERTEX_REVERSE_BIT = 1u;
    constexpr SplitEnd CANDIDATE_SPLIT_ENDS[2] = {SPLIT_RIGHT_END, SPLIT_LEFT_END};

    // ------------------------------------------------ Collect misassembly ends ------------------------------------------------
    struct Split {
        uint32_t source{UINT32_MAX}, begin{0}, end{0};  // Fragment ID and vertex range to split
        const char* side{nullptr};  // left/right/whole
        std::string sequence;
    };

    std::vector<Split> splits;
    std::unordered_set<uint64_t> seen;
    for (const Candidate& candidate : candidates) {
        const uint32_t sources[2] = {candidate.left, candidate.right};
        const uint32_t begins[2] = {candidate.left_boundary.junction_target_pos + 1, 0};
        const uint32_t ends[2] = {
            static_cast<uint32_t>(fragments_[candidate.left].vertices.size()),
            candidate.right_boundary.junction_target_pos
        };
        const bool flagged[2] = {candidate.left_misassembly, candidate.right_misassembly};
        for (uint32_t side = 0; side < 2; ++side) {
            if (!flagged[side] || begins[side] >= ends[side]) continue;
            const uint64_t key = (static_cast<uint64_t>(sources[side]) << SPLIT_END_BITS) | CANDIDATE_SPLIT_ENDS[side];
            if (!seen.insert(key).second) continue;
            splits.push_back({
                sources[side], begins[side], ends[side],
                CANDIDATE_SPLIT_ENDS[side] == SPLIT_RIGHT_END ? "right" : "left", {}
            });
        }
    }

    // ------------------------------------------------ Find matching terminal blocks ------------------------------------------------
    // Terminal block can occur on more than one misassembled haplotype. Split every contig where that block is also terminal.
    const size_t seed_count = splits.size();
    for (size_t seed_id = 0; seed_id < seed_count; ++seed_id) {
        const Split& seed = splits[seed_id];
        const Fragment& source = fragments_[seed.source];
        std::unordered_set<uint32_t> terminal_nodes;  // segment IDs
        uint64_t terminal_bp = 0;
        for (uint32_t pos = seed.begin; pos < seed.end; ++pos) {
            const uint32_t sid = Vertex::get_segment_id(source.vertices[pos]);
            if (terminal_nodes.insert(sid).second) terminal_bp += nodes_[sid].length;
        }
        if (terminal_bp == 0) continue;

        for (uint32_t fragment_id = 0; fragment_id < fragments_.size(); ++fragment_id) {
            if (fragment_id == seed.source) continue;
            const Fragment& fragment = fragments_[fragment_id];
            uint32_t first = UINT32_MAX, last = 0;
            uint64_t shared_bp = 0;
            std::unordered_set<uint32_t> shared_nodes;  // segment IDs
            for (uint32_t pos = 0; pos < fragment.vertices.size(); ++pos) {
                const uint32_t sid = Vertex::get_segment_id(fragment.vertices[pos]);
                if (!terminal_nodes.contains(sid) || !shared_nodes.insert(sid).second) continue;
                first = std::min(first, pos);
                last = std::max(last, pos);
                shared_bp += nodes_[sid].length;
            }
            if (first == UINT32_MAX || shared_bp < terminal_bp * params_.misassembly_similarity) continue;

            const uint64_t prefix_bp = fragment.path_bp[first];
            const uint64_t suffix_bp = fragment.length - fragment.path_bp[last + 1];
            const bool prefix = prefix_bp <= params_.end_skip_bp;
            const bool suffix = suffix_bp <= params_.end_skip_bp;
            if (!prefix && !suffix) continue;
            const uint64_t fragment_key = static_cast<uint64_t>(fragment_id) << SPLIT_END_BITS;
            const uint64_t right_end_key = fragment_key | SPLIT_RIGHT_END;
            const uint64_t left_end_key = fragment_key | SPLIT_LEFT_END;
            if (prefix && suffix) {
                if (seen.contains(left_end_key) || seen.contains(right_end_key)) continue;
                seen.insert(left_end_key);
                seen.insert(right_end_key);
                splits.push_back({fragment_id, 0, static_cast<uint32_t>(fragment.vertices.size()), "whole", {}});
            } else if (prefix) {
                if (!seen.insert(left_end_key).second) continue;
                splits.push_back({fragment_id, 0, last + 1, "left", {}});
            } else {
                if (!seen.insert(right_end_key).second) continue;
                splits.push_back({fragment_id, first, static_cast<uint32_t>(fragment.vertices.size()), "right", {}});
            }
        }
    }

    // ------------------------------------------------ Build split sequences ------------------------------------------------
    std::vector<uint8_t> split_nodes(nodes_.size(), false);  // Nodes needed to be deleted after all paths have been rewritten
    size_t valid_splits = 0;
    for (size_t i = 0; i < splits.size(); ++i) {
        Split& split = splits[i];
        split.sequence = fragment_sequence_(fragments_[split.source], split.begin, split.end);
        if (split.sequence.empty()) {
            warning_stream() << "  ! Cannot split misassembly without sequence: " << paths_[fragments_[split.source].path_id].name << '\n';
            continue;
        }
        for (uint32_t pos = split.begin; pos < split.end; ++pos) {
            split_nodes[Vertex::get_segment_id(fragments_[split.source].vertices[pos])] = true;
        }
        if (valid_splits != i) splits[valid_splits] = std::move(split);
        ++valid_splits;
    }
    splits.resize(valid_splits);

    // ------------------------------------------------ Find retained source ranges ------------------------------------------------
    std::vector<uint32_t> retained_begin(fragments_.size(), 0);
    std::vector<uint32_t> retained_end;
    retained_end.reserve(fragments_.size());
    for (const Fragment& fragment : fragments_) retained_end.push_back(fragment.vertices.size());
    for (const Split& split : splits) {
        if (split.begin == 0) retained_begin[split.source] = std::max(retained_begin[split.source], split.end);
        if (split.end == fragments_[split.source].vertices.size()) {
            retained_end[split.source] = std::min(retained_end[split.source], split.begin);
        }
    }

    // ------------------------------------------------ Cut fully contig from backbone ------------------------------------------------
    std::vector<uint8_t> valid_source(fragments_.size(), true);
    for (const Split& split : splits) {
        if (retained_begin[split.source] >= retained_end[split.source]) valid_source[split.source] = false;
    }

    // ------------------------------------------------ Create single misassembly paths ------------------------------------------------
    std::unordered_set<std::string> names;
    names.reserve(nodes_.size() + paths_.size());
    for (const GfaNode& node : nodes_) names.insert(node.name);
    for (const GfaPath& path : paths_) names.insert(path.name);

    std::vector<MisassemblyContig> records;
    records.reserve(splits.size());
    size_t number = 0;
    for (const Split& split : splits) {
        const Fragment& source = fragments_[split.source];
        std::string name;
        do {
            name = misassembly_contig_name(source.sample, source.hap, ++number);
        } while (!names.insert(name).second);

        std::vector<uint32_t> sample_ids;
        const auto sample = sample_name_to_id_.find(source.sample);
        if (sample != sample_name_to_id_.end()) sample_ids.push_back(sample->second);
        // Create a new segment
        const uint32_t segment = add_segment(name, split.sequence, true, sample_ids);
        // Create a new path
        GfaPath path;
        path.name = name;
        path.segments.push_back({segment, false});
        paths_.push_back(std::move(path));
        records.push_back({
            name, source.sample, paths_[source.path_id].name,
            split.side, source.hap, source.component
        });
    }

    // ------------------------------------------------ Trim original paths ------------------------------------------------
    size_t trimmed = 0;
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        if (retained_begin[i] == 0 && retained_end[i] == fragments_[i].vertices.size()) continue;
        const Fragment& source = fragments_[i];
        if (!valid_source[i]) {
            paths_[source.path_id].segments.clear();
            ++trimmed;
            continue;
        }
        std::vector<uint32_t> retained(
            source.vertices.begin() + retained_begin[i],
            source.vertices.begin() + retained_end[i]
        );
        if (source.reverse) {
            std::reverse(retained.begin(), retained.end());
            for (uint32_t& vertex : retained) vertex ^= VERTEX_REVERSE_BIT;
        }
        GfaPath& path = paths_[source.path_id];
        path.segments.clear();
        path.segments.reserve(retained.size());
        for (uint32_t vertex : retained) {
            path.segments.push_back({
                Vertex::get_segment_id(vertex), Vertex::get_is_reverse(vertex)
            });
        }
        ++trimmed;
    }
    paths_.erase(
        std::remove_if(paths_.begin(), paths_.end(), [](const GfaPath& path) {
            return path.segments.empty();
        }),
        paths_.end()
    );

    // ------------------------------------------------ Remove unused split nodes ------------------------------------------------
    // Delete old split nodes
    std::vector<uint8_t> path_nodes(nodes_.size(), false);
    for (const GfaPath& path : paths_) {
        std::string sample;
        uint8_t hap = 0;
        if (!parse_path_name_(path.name, sample, hap)) continue;
        for (const PathSegment& segment : path.segments) {
            if (segment.node_id < path_nodes.size()) path_nodes[segment.node_id] = true;
        }
    }
    size_t removed = 0;
    for (uint32_t sid = 0; sid < split_nodes.size(); ++sid) {
        if (split_nodes[sid] && !path_nodes[sid] && delete_segment(sid)) ++removed;
    }
    if (removed > 0) {
        paths_.erase(
            std::remove_if(paths_.begin(), paths_.end(), [&](const GfaPath& path) {
                return std::any_of(path.segments.begin(), path.segments.end(), [&](const PathSegment& segment) {
                    return segment.node_id < nodes_.size() && nodes_[segment.node_id].deleted;
                });
            }),
            paths_.end()
        );
        rebuild_after_edits();
    }

    log_stream() << "Splitting terminal misassemblies into single contigs ...\n";
    log_stream() << "  - Misassembly contigs created: " << records.size() << '\n';
    log_stream() << "  - Misassembly nodes removed: " << removed << "\n\n";
    return records;
}

void GfaGapfill::mark_used_relocations_(const std::vector<Candidate>& selected) {
    for (const Candidate& candidate : selected) {
        if (candidate.right_boundary.bridge_pos <= candidate.left_boundary.bridge_pos + 1) continue;
        const uint32_t ids[3] = {candidate.left, candidate.right, candidate.bridge};
        for (uint32_t id : ids) {
            const int32_t relocation = fragments_[id].relocation;
            if (relocation >= 0) relocations_[relocation].used_for_gap = true;
        }
    }
}


// ================================================= Assign haplotypes =================================================
void GfaGapfill::build_sample_chains_(const std::vector<Candidate>& selected) {
    log_stream() << "Building gap-filled sample paths ...\n\n";

    // ------------------------------------------------ Build selected connection indexes ------------------------------------------------
    std::vector<int32_t> incoming(fragments_.size(), -1), outgoing(fragments_.size(), -1);
    for (size_t i = 0; i < selected.size(); ++i) {
        incoming[selected[i].right] = static_cast<int32_t>(i);
        outgoing[selected[i].left] = static_cast<int32_t>(i);
    }

    // ------------------------------------------------ Build chains from selected connections ------------------------------------------------
    std::map<std::pair<std::string, uint32_t>, std::vector<Chain>> groups;  // Sample name + component -> chains
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        if (incoming[i] >= 0) continue;

        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible || fragment.redundant) continue;  // Short fragment or redundant fragment which is already included in a gap-filled chain
        Chain chain = build_chain_(i, incoming, outgoing, selected);
        if (chain.vertices.empty()) continue;
        // Example: sampleA + component3 -> [chain1, chain2, ...]
        groups[{fragment.sample, fragment.component}].push_back(std::move(chain));
    }

    // ------------------------------------------------ Compare overlapping chains based on shared nodes ------------------------------------------------
    for (auto& group : groups) {
        std::vector<Chain>& chains = group.second;
        // ------------------------------------------------ Find shared graph nodes ------------------------------------------------
        std::vector<std::pair<uint32_t, uint32_t>> node_chains;  // (segment_id, chain_id)
        size_t total_nodes = 0;
        for (const Chain& chain : chains) total_nodes += chain.vertices.size();
        node_chains.reserve(total_nodes);
        for (uint32_t i = 0; i < chains.size(); ++i) {
            for (uint32_t vertex : chains[i].vertices) {
                node_chains.emplace_back(Vertex::get_segment_id(vertex), i);
            }
        }
        std::sort(node_chains.begin(), node_chains.end());

        const uint64_t min_shared = static_cast<uint64_t>(std::ceil(params_.min_overlap_bp * params_.min_similarity));

        std::unordered_map<uint64_t, uint64_t> shared_bp;  // Key = (chain_a << 32) | chain_b, Value = shared base pairs
        for (size_t begin = 0; begin < node_chains.size();) {
            size_t end = begin + 1;
            while (end < node_chains.size() && node_chains[end].first == node_chains[begin].first) {
                ++end;
            }
            for (size_t i = begin; i < end; ++i) {
                if (i > begin && node_chains[i].second == node_chains[i - 1].second) continue;
                for (size_t j = i + 1; j < end; ++j) {
                    if (node_chains[j].second == node_chains[j - 1].second) continue;
                    const uint32_t a = node_chains[i].second;
                    const uint32_t b = node_chains[j].second;
                    const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;
                    shared_bp[key] += nodes_[node_chains[begin].first].length;
                }
            }
            begin = end;
        }

        // ------------------------------------------------ Compare overlapping chains based on layout ------------------------------------------------
        std::vector<std::pair<int64_t, uint32_t>> chain_order;  // (chain_begin, chain_id)
        std::vector<std::pair<int64_t, int64_t>> chain_bounds(chains.size());  // (chain_begin, chain_end)
        chain_order.reserve(chains.size());
        for (uint32_t i = 0; i < chains.size(); ++i) {
            int64_t begin = std::numeric_limits<int64_t>::max();
            int64_t end = std::numeric_limits<int64_t>::min();
            for (uint32_t fragment_id : chains[i].source_fragments) {
                const Fragment& fragment = fragments_[fragment_id];
                begin = std::min(begin, fragment.layout_start);
                end = std::max(end, fragment.layout_start + static_cast<int64_t>(fragment.length));
            }
            chain_bounds[i] = {begin, end};
            chain_order.emplace_back(begin, i);
        }
        std::sort(chain_order.begin(), chain_order.end());
        std::vector<int32_t> best_partner(chains.size(), -1);
        std::vector<uint64_t> best_overlap(chains.size(), 0);
        for (size_t i = 0; i < chain_order.size(); ++i) {
            const uint32_t a = chain_order[i].second;
            const bool mixed_a = chains[a].hap_bp[0] > 0 && chains[a].hap_bp[1] > 0;
            for (size_t j = i + 1; j < chain_order.size(); ++j) {
                const uint32_t b = chain_order[j].second;
                if (chain_bounds[b].first >= chain_bounds[a].second) break;  // no overlap
                const bool mixed_b = chains[b].hap_bp[0] > 0 && chains[b].hap_bp[1] > 0;
                if (!mixed_a && !mixed_b) continue;  // both chains are already assigned to a single haplotype
                const int64_t overlap = std::min(chain_bounds[a].second, chain_bounds[b].second) - chain_bounds[b].first;
                const int64_t shorter = std::min(chain_bounds[a].second - chain_bounds[a].first, chain_bounds[b].second - chain_bounds[b].first);
                if (overlap < static_cast<int64_t>(min_shared) || overlap * 2 < shorter) continue;
                const uint64_t overlap_bp = static_cast<uint64_t>(overlap);
                if (overlap_bp > best_overlap[a]) {
                    best_overlap[a] = overlap_bp;
                    best_partner[a] = static_cast<int32_t>(b);
                }
                if (overlap_bp > best_overlap[b]) {
                    best_overlap[b] = overlap_bp;
                    best_partner[b] = static_cast<int32_t>(a);
                }
            }
        }
        for (uint32_t a = 0; a < chains.size(); ++a) {
            if (best_partner[a] < 0) continue;
            const uint32_t b = static_cast<uint32_t>(best_partner[a]);
            const uint32_t first = std::min(a, b);
            const uint32_t second = std::max(a, b);
            const uint64_t key = (static_cast<uint64_t>(first) << 32) | second;
            shared_bp[key] = std::max(shared_bp[key], best_overlap[a]);
        }

        // ------------------------------------------------ Build chain conflict graph ------------------------------------------------
        std::vector<std::vector<uint32_t>> conflicts(chains.size());  // Chain ID -> conflicting chain IDs
        for (const auto& item : shared_bp) {
            if (item.second < min_shared) continue;
            const uint32_t a = static_cast<uint32_t>(item.first >> 32);
            const uint32_t b = static_cast<uint32_t>(item.first);
            conflicts[a].push_back(b);
            conflicts[b].push_back(a);
        }

        // ------------------------------------------------ Assign chain haplotypes ------------------------------------------------
        std::vector<int8_t> color(chains.size(), -1);
        std::array<uint64_t, 2> component_bp{0, 0};
        for (uint32_t root = 0; root < chains.size(); ++root) {
            if (color[root] >= 0) continue;
            std::queue<uint32_t> queue;
            std::vector<uint32_t> block;
            bool conflict = false;
            color[root] = 0;
            queue.push(root);
            while (!queue.empty()) {
                const uint32_t u = queue.front();
                queue.pop();
                block.push_back(u);
                for (uint32_t v : conflicts[u]) {
                    if (color[v] < 0) {
                        color[v] = color[u] ^ 1;
                        queue.push(v);
                    } else if (color[v] == color[u]) {
                        conflict = true;
                    }
                }
            }

            uint64_t normal = 0;
            uint64_t swapped = 0;
            uint64_t color_bp[2] = {0, 0};
            for (uint32_t i : block) {
                const Chain& chain = chains[i];
                const uint64_t length = chain_length_(chain);
                color_bp[color[i]] += length;
                if (color[i] == 0) {
                    normal += chain.hap_bp[0];
                    swapped += chain.hap_bp[1];
                } else {
                    normal += chain.hap_bp[1];
                    swapped += chain.hap_bp[0];
                }
            }
            bool swap = swapped > normal;
            if (swapped == normal) {
                const uint64_t normal_hap1 = component_bp[0] + color_bp[0];
                const uint64_t normal_hap2 = component_bp[1] + color_bp[1];
                const uint64_t swapped_hap1 = component_bp[0] + color_bp[1];
                const uint64_t swapped_hap2 = component_bp[1] + color_bp[0];
                const uint64_t normal_diff = normal_hap1 > normal_hap2 ? normal_hap1 - normal_hap2 : normal_hap2 - normal_hap1;
                const uint64_t swapped_diff = swapped_hap1 > swapped_hap2 ? swapped_hap1 - swapped_hap2 : swapped_hap2 - swapped_hap1;
                swap = swapped_diff < normal_diff;
            }
            for (uint32_t i : block) {
                const uint8_t assigned = static_cast<uint8_t>(color[i]) ^ swap;
                Chain& chain = chains[i];
                chain.hap = assigned + 1;
                component_bp[assigned] += chain_length_(chain);
                GfaGapfillDebugger::assignment(
                    chain.sample, chain.component, chain.hap,
                    chain.hap_bp[0], chain.hap_bp[1], conflict
                );
                sample_chains_[chain.sample][assigned].push_back(std::move(chain));
            }
        }
    }

    // ------------------------------------------------ Keep short single fragments ------------------------------------------------
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (fragment.eligible || fragment.relocation >= 0 || fragment.redundant) continue;
        Chain chain;
        chain.sample = fragment.sample;
        chain.hap = fragment.hap;
        chain.component = fragment.component;
        chain.source_fragments.push_back(i);
        chain.vertices.push_back(fragment.vertices.front());
        chain.parts.push_back({i, 0, fragment.length, false});
        sample_chains_[fragment.sample][fragment.hap - 1].push_back(std::move(chain));
    }
}

GfaGapfill::Chain GfaGapfill::build_chain_(
    uint32_t start,
    const std::vector<int32_t>& incoming,
    const std::vector<int32_t>& outgoing,
    const std::vector<Candidate>& selected
) const {
    Chain chain;
    chain.sample = fragments_[start].sample;
    chain.hap = fragments_[start].hap;
    chain.component = fragments_[start].component;

    uint32_t current = start;

    while (current != UINT32_MAX) {
        const Fragment& fragment = fragments_[current];
        chain.source_fragments.push_back(current);
        const uint32_t begin = incoming[current] < 0 ? 0 : selected[incoming[current]].right_boundary.target_pos;
        const uint32_t end = outgoing[current] < 0 ? static_cast<uint32_t>(fragment.vertices.size() - 1) : selected[outgoing[current]].left_boundary.target_pos;
        const uint64_t begin_bp = incoming[current] < 0 ? 0 : selected[incoming[current]].right_boundary.target_cut_bp;
        const uint64_t end_bp = outgoing[current] < 0 ? fragment.length : selected[outgoing[current]].left_boundary.target_cut_bp;
        if (begin_bp < end_bp) {
            chain.parts.push_back({current, begin_bp, end_bp, false});
            chain.hap_bp[fragment.hap - 1] += end_bp - begin_bp;
        }
        for (uint32_t i = begin; i <= end; ++i) {
            chain.vertices.push_back(fragment.vertices[i]);
        }
        if (outgoing[current] < 0) break;

        const Candidate& edge = selected[outgoing[current]];
        const Fragment& bridge = fragments_[edge.bridge];
        if (edge.left_boundary.bridge_cut_bp <= edge.right_boundary.bridge_cut_bp) {
            chain.parts.push_back({
                edge.bridge, edge.left_boundary.bridge_cut_bp,
                edge.right_boundary.bridge_cut_bp, true
            });
        }
        if (edge.left_boundary.bridge_cut_bp < edge.right_boundary.bridge_cut_bp) {
            for (uint32_t i = edge.left_boundary.bridge_pos + 1; i < edge.right_boundary.bridge_pos; ++i) {
                chain.vertices.push_back(bridge.vertices[i]);
            }
        }
        chain.gaps.push_back(edge);
        current = edge.right;
    }
    return chain;
}


// ================================================= Deduplication =================================================
void GfaGapfill::mark_redundant_fragments_(const std::vector<Candidate>& selected) {
    log_stream() << "Checking contigs covered by filled gaps ...\n";

    // ------------------------------------------------ Group selected gaps ------------------------------------------------
    std::vector<uint8_t> used(fragments_.size(), false);
    std::vector<uint8_t> assigned_hap(fragments_.size(), 0);  // 0 = unassigned, 1 = hap1, 2 = hap2
    for (const auto& sample : sample_chains_) {
        for (const auto& haplotype : sample.second) {
            for (const Chain& chain : haplotype) {
                for (uint32_t fragment_id : chain.source_fragments) {
                    assigned_hap[fragment_id] = chain.hap;
                }
            }
        }
    }

    using Group = std::tuple<std::string, uint32_t, uint8_t>;  // sample, component, hap
    std::map<Group, std::vector<size_t>> gaps_by_group;
    for (size_t candidate_id = 0; candidate_id < selected.size(); ++candidate_id) {
        const Candidate& candidate = selected[candidate_id];
        used[candidate.left] = true;
        used[candidate.right] = true;
        used[candidate.bridge] = true;
        const Fragment& left = fragments_[candidate.left];
        /**************************** may have bugs ********************************/
        // Coverage for left/right with different haplotypes uses the final chain hap. (2026-08-11 v0.1.3-r21)
        const uint8_t hap = assigned_hap[candidate.left];
        if (hap > 0) gaps_by_group[{left.sample, left.component, hap}].push_back(candidate_id);
        /**************************** may have bugs ********************************/
    }

    // ------------------------------------------------ Check unused fragments ------------------------------------------------
    std::vector<std::string> gap_sequences(selected.size());
    size_t redundant = 0;
    uint64_t redundant_bp = 0;
    for (uint32_t fragment_id = 0; fragment_id < fragments_.size(); ++fragment_id) {
        Fragment& fragment = fragments_[fragment_id];
        if (!fragment.eligible || used[fragment_id]) continue;
        const uint8_t hap = assigned_hap[fragment_id];
        if (hap == 0) continue;
        const auto group = gaps_by_group.find({fragment.sample, fragment.component, hap});
        if (group == gaps_by_group.end()) continue;

        // ------------------------------------------------ Check layout coverage ------------------------------------------------
        const int64_t fragment_begin = fragment.layout_start;
        const int64_t fragment_end = fragment.layout_start + static_cast<int64_t>(fragment.length);
        std::string sequence;
        bool sequence_loaded = false;
        for (size_t candidate_id : group->second) {
            const Candidate& candidate = selected[candidate_id];
            const Fragment& bridge = fragments_[candidate.bridge];
            if (candidate.right_boundary.bridge_pos <= candidate.left_boundary.bridge_pos + 1) continue;

            const int64_t gap_begin = bridge.layout_start + static_cast<int64_t>(candidate.left_boundary.bridge_cut_bp);
            const int64_t gap_end = bridge.layout_start + static_cast<int64_t>(candidate.right_boundary.bridge_cut_bp);
            const int64_t overlap_begin = std::max(fragment_begin, gap_begin);
            const int64_t overlap_end = std::min(fragment_end, gap_end);
            if (overlap_end <= overlap_begin || static_cast<uint64_t>(overlap_end - overlap_begin) * 10 < fragment.length * 8) continue;

            // ------------------------------------------------ Compare covered sequences ------------------------------------------------
            std::string& gap_sequence = gap_sequences[candidate_id];
            if (gap_sequence.empty()) {
                gap_sequence = fragment_subsequence_(
                    bridge, candidate.left_boundary.bridge_cut_bp,
                    candidate.right_boundary.bridge_cut_bp
                );
            }
            if (!sequence_loaded) {
                sequence = fragment_sequence_(fragment, 0, static_cast<uint32_t>(fragment.vertices.size()));
                sequence_loaded = true;
            }
            const double similarity = mm2_similarity_(gap_sequence, sequence);
            if (similarity < params_.dedup_similarity) continue;

            // ------------------------------------------------ Mark redundant fragment ------------------------------------------------
            fragment.redundant = true;
            ++redundant;
            redundant_bp += fragment.length;
            log_stream() << "  - Covered contig: " << paths_[fragment.path_id].name << " by " << paths_[bridge.path_id].name << " (similarity=" << std::fixed << std::setprecision(4) << similarity << ")\n";
            break;
        }
    }

    // ------------------------------------------------ Remove redundant chains ------------------------------------------------
    for (auto& sample : sample_chains_) {
        for (auto& haplotype : sample.second) {
            haplotype.erase(
                std::remove_if(haplotype.begin(), haplotype.end(), [&](const Chain& chain) {
                    return chain.source_fragments.size() == 1 &&
                        fragments_[chain.source_fragments.front()].redundant;
                }),
                haplotype.end()
            );
        }
    }
    log_stream() << "  - Redundant contigs removed: " << redundant << " (" << redundant_bp << " bp)\n\n";
}

double GfaGapfill::mm2_similarity_(
    const std::string& reference,
    const std::string& query
) const {
    if (reference.empty() || query.empty()) return -1.0;

    mm_idxopt_t index_options;
    mm_mapopt_t map_options;
    mm_set_opt(nullptr, &index_options, &map_options);
    if (mm_set_opt(mm2_.preset.c_str(), &index_options, &map_options) < 0) return -1.0;
    index_options.k = static_cast<short>(mm2_.k);
    index_options.w = static_cast<short>(mm2_.w);
    map_options.flag |= MM_F_CIGAR | MM_F_EQX;
    map_options.best_n = static_cast<short>(mm2_.best_n);
    map_options.zdrop = mm2_.zdrop;

    const char* reference_sequences[1] = {reference.c_str()};
    const char* reference_names[1] = {"gap"};
    const int hpc = (index_options.flag & MM_I_HPC) ? 1 : 0;
    mm_idx_t* index = mm_idx_str(
        index_options.w, index_options.k, hpc, index_options.bucket_bits,
        1, reference_sequences, reference_names
    );
    if (!index) return -1.0;
    mm_mapopt_update(&map_options, index);

    mm_tbuf_t* buffer = mm_tbuf_init();
    if (!buffer) {
        mm_idx_destroy(index);
        return -1.0;
    }
    int count = 0;
    mm_reg1_t* hits = mm_map(
        index, static_cast<int>(query.size()), query.c_str(),
        &count, buffer, &map_options, "unplaced"
    );

    uint32_t best_matches = 0;
    for (int i = 0; i < count; ++i) {
        const mm_reg1_t& hit = hits[i];
        if (!hit.p || hit.rev || hit.parent != hit.id || hit.mapq < mm2_.min_mapq || hit.blen <= 0) continue;
        const double match_ratio = static_cast<double>(hit.mlen) / hit.blen;
        const uint64_t aligned_span = std::max<int64_t>(hit.re - hit.rs, hit.qe - hit.qs);
        const size_t shorter = std::min(reference.size(), query.size());
        const double alignment_ratio = shorter == 0 ? 0.0 : static_cast<double>(aligned_span) / shorter;
        if (match_ratio < mm2_.min_match || alignment_ratio < mm2_.min_ali_ratio) continue;
        best_matches = std::max(best_matches, static_cast<uint32_t>(hit.mlen));
    }
    if (hits) {
        for (int i = 0; i < count; ++i) std::free(hits[i].p);
        std::free(hits);
    }
    mm_tbuf_destroy(buffer);
    mm_idx_destroy(index);

    return std::min(1.0, static_cast<double>(best_matches) / query.size());
}


// ================================================= Select the primary chains =================================================
GfaGapfill::PrimarySelection GfaGapfill::select_primary_chains_() const {
    log_stream() << "Selecting primary haplotypes by component ...\n";

    // ------------------------------------------------ Group chains by component and sample ------------------------------------------------
    using HapChains = std::array<std::vector<const Chain*>, 2>;  // haplotype -> chains
    std::map<uint32_t, std::map<std::string, HapChains>> components;  // component -> sample -> haplotype -> chains
    std::unordered_map<const Chain*, uint64_t> chain_lengths;  // Chain pointer -> chain length
    ChainRefs unmapped_misassemblies;
    for (const auto& sample : sample_chains_) {
        for (uint8_t hap = 0; hap < 2; ++hap) {
            for (const Chain& chain : sample.second[hap]) {
                chain_lengths.emplace(&chain, chain_length_(chain));
                bool unmapped = false;
                for (uint32_t id : chain.source_fragments) {
                    const int32_t relocation = fragments_[id].relocation;
                    if (relocation >= 0 && !relocations_[relocation].used_for_gap && relocations_[relocation].status == "unmapped") {
                        unmapped = true;
                        break;
                    }
                }
                if (unmapped) unmapped_misassemblies[hap].push_back(&chain);
                else components[chain.component][sample.first][hap].push_back(&chain);
            }
        }
    }

    // ------------------------------------------------ Estimate component sizes ------------------------------------------------
    std::map<uint32_t, uint64_t> component_sizes;  // keey = component, value = total base pairs
    for (const auto& component : components) {
        uint64_t component_bp = 0;
        for (const auto& sample : component.second) {
            for (uint8_t hap = 0; hap < 2; ++hap) {
                uint64_t total = 0;
                for (const Chain* chain : sample.second[hap]) {
                    total += chain_lengths.at(chain);
                }
                component_bp = std::max(component_bp, total);
            }
        }
        component_sizes.emplace(component.first, component_bp);
    }

    // ------------------------------------------------ Select the best sample per component ------------------------------------------------
    std::map<uint32_t, HapChains> primary_components;  // key = component, value = haplotype -> chains
    for (const auto& component : components) {
        const uint64_t component_bp = component_sizes.at(component.first);

        const HapChains* best = nullptr;
        std::string best_sample;
        std::tuple<bool, uint64_t, uint64_t, uint64_t, uint64_t> best_score;  // complete_pair, sum_n50, min_n50, total_bp, contigs
        uint64_t best_n50[2] = {0, 0};
        size_t best_contigs = 0;

        for (const auto& sample : component.second) {
            std::vector<uint64_t> lengths[2];
            uint64_t total = 0;
            size_t contigs = 0;
            for (uint8_t hap = 0; hap < 2; ++hap) {
                lengths[hap].reserve(sample.second[hap].size());
                for (const Chain* chain : sample.second[hap]) {
                    const uint64_t length = chain_lengths.at(chain);
                    lengths[hap].push_back(length);
                    total += length;
                    ++contigs;
                }
            }
            const uint64_t hap_n50[2] = {
                ng50_(std::move(lengths[0]), component_bp),
                ng50_(std::move(lengths[1]), component_bp)
            };
            const bool complete_pair = hap_n50[0] > 0 && hap_n50[1] > 0;
            const auto score = std::make_tuple(
                complete_pair,
                hap_n50[0] + hap_n50[1],
                std::min(hap_n50[0], hap_n50[1]),
                total,
                UINT64_MAX - static_cast<uint64_t>(contigs)
            );
            GfaGapfillDebugger::primary_candidate(component.first, sample.first, component_bp, hap_n50[0], hap_n50[1]);
            const bool better = best == nullptr || score > best_score || (score == best_score && sample.first < best_sample);
            if (!better) continue;
            best = &sample.second;
            best_sample = sample.first;
            best_score = score;
            best_n50[0] = hap_n50[0];
            best_n50[1] = hap_n50[1];
            best_contigs = contigs;
        }
        if (best == nullptr) continue;
        primary_components[component.first] = *best;
        log_stream() << "  - Primary component " << component.first << ": " << best_sample
                     << " (component_bp=" << component_bp
                     << ", hap1_NG50=" << best_n50[0]
                     << ", hap2_NG50=" << best_n50[1]
                     << ", mean_NG50=" << (best_n50[0] + best_n50[1]) / 2
                     << ", contigs=" << best_contigs << ")\n";
    }
    std::cerr << '\n';

    // ------------------------------------------------ Remove small and unmapped chains ------------------------------------------------
    PrimarySelection selection;
    selection.low_quality = unmatched_small_primary_chains_(primary_components, chain_lengths);
    for (uint8_t hap = 0; hap < 2; ++hap) {
        selection.low_quality[hap].insert(
            selection.low_quality[hap].end(),
            unmapped_misassemblies[hap].begin(), unmapped_misassemblies[hap].end()
        );
    }
    // ------------------------------------------------ Drop tiny primary components ------------------------------------------------
    for (auto component = primary_components.begin(); component != primary_components.end();) {
        uint64_t hap_bp[2]{0, 0};
        for (uint8_t hap = 0; hap < 2; ++hap) {
            for (const Chain* chain : component->second[hap]) hap_bp[hap] += chain_lengths.at(chain);
        }
        if (std::max(hap_bp[0], hap_bp[1]) <= params_.dedup_component_bp) {
            component = primary_components.erase(component);
        } else {
            ++component;
        }
    }
    // ------------------------------------------------ Order components by haplotype balance ------------------------------------------------
    std::vector<std::pair<uint64_t, uint32_t>> component_order;
    component_order.reserve(primary_components.size());
    for (const auto& component : primary_components) {
        uint64_t bp[2]{0, 0};
        for (uint8_t hap = 0; hap < 2; ++hap) {
            for (const Chain* chain : component.second[hap]) bp[hap] += chain_lengths.at(chain);
        }
        const uint64_t difference = bp[0] > bp[1] ? bp[0] - bp[1] : bp[1] - bp[0];
        component_order.emplace_back(UINT64_MAX - difference, component.first);
    }
    std::sort(component_order.begin(), component_order.end());

    // ------------------------------------------------ Align haplotype labels across components ------------------------------------------------
    uint64_t primary_bp[2]{0, 0};
    for (const auto& ordered_component : component_order) {
        const uint32_t component_id = ordered_component.second;
        HapChains& pair = primary_components[component_id];
        uint64_t bp[2]{0, 0};
        for (uint8_t hap = 0; hap < 2; ++hap) {
            for (const Chain* chain : pair[hap]) bp[hap] += chain_lengths.at(chain);
        }
        const uint64_t normal_a = primary_bp[0] + bp[0];
        const uint64_t normal_b = primary_bp[1] + bp[1];
        const uint64_t swapped_a = primary_bp[0] + bp[1];
        const uint64_t swapped_b = primary_bp[1] + bp[0];
        const uint64_t normal_diff = normal_a > normal_b ? normal_a - normal_b : normal_b - normal_a;
        const uint64_t swapped_diff = swapped_a > swapped_b ? swapped_a - swapped_b : swapped_b - swapped_a;
        if (swapped_diff < normal_diff) {
            std::swap(pair[0], pair[1]);
            std::swap(bp[0], bp[1]);
        }
        primary_bp[0] += bp[0];
        primary_bp[1] += bp[1];
    }
    // ------------------------------------------------ Collect final primary chains ------------------------------------------------
    for (const auto& component : primary_components) {
        selection.primary[0].insert(selection.primary[0].end(), component.second[0].begin(), component.second[0].end());
        selection.primary[1].insert(selection.primary[1].end(), component.second[1].begin(), component.second[1].end());
    }
    return selection;
}

// Test whether a small component is almost completely represented by a larger haplotype.
bool sequence_contained(
    const mm_idx_t* index,
    const mm_mapopt_t& map_options,
    const std::string& sequence,
    uint8_t min_mapq,
    double min_match,
    double min_similarity
) {
    mm_tbuf_t* buffer = mm_tbuf_init();
    if (!buffer) return false;
    int count = 0;
    mm_reg1_t* hits = mm_map(
        index, static_cast<int>(sequence.size()), sequence.c_str(),
        &count, buffer, &map_options, "small_component"
    );
    bool contained = false;
    for (int i = 0; i < count && !contained; ++i) {
        const mm_reg1_t& hit = hits[i];
        if (!hit.p || hit.parent != hit.id || hit.mapq < min_mapq || hit.blen <= 0) continue;
        const double match_ratio = static_cast<double>(hit.mlen) / hit.blen;
        const double contained_ratio = static_cast<double>(hit.mlen) / sequence.size();
        contained = match_ratio >= min_match && contained_ratio >= min_similarity;
    }
    if (hits) {
        for (int i = 0; i < count; ++i) std::free(hits[i].p);
        std::free(hits);
    }
    mm_tbuf_destroy(buffer);
    return contained;
}

GfaGapfill::ChainRefs GfaGapfill::unmatched_small_primary_chains_(
    const std::map<uint32_t, ChainRefs>& components,
    const std::unordered_map<const Chain*, uint64_t>& chain_lengths
) const {
    // ------------------------------------------------ Split small and large components ------------------------------------------------
    struct Query {
        const Chain* chain;
        std::string sequence;
    };

    ChainRefs small_chains;
    ChainRefs large;
    for (const auto& component : components) {
        uint64_t hap_bp[2]{0, 0};
        for (uint8_t hap = 0; hap < 2; ++hap) {
            for (const Chain* chain : component.second[hap]) hap_bp[hap] += chain_lengths.at(chain);
        }
        const bool is_small = std::max(hap_bp[0], hap_bp[1]) <= params_.dedup_component_bp;
        for (uint8_t hap = 0; hap < 2; ++hap) {
            auto& destination = is_small ? small_chains[hap] : large[hap];
            destination.insert(destination.end(), component.second[hap].begin(), component.second[hap].end());
        }
    }

    // ------------------------------------------------ Prepare mm2 options ------------------------------------------------
    log_stream() << "Aligning small primary components to larger haplotypes ...\n";
    mm_idxopt_t index_options;
    mm_mapopt_t default_map_options;
    mm_set_opt(nullptr, &index_options, &default_map_options);
    ChainRefs low_quality;
    if (mm_set_opt(mm2_.preset.c_str(), &index_options, &default_map_options) < 0) return small_chains;
    index_options.k = static_cast<short>(mm2_.k);
    index_options.w = static_cast<short>(mm2_.w);

    // ------------------------------------------------ Check each haplotype separately ------------------------------------------------
    for (uint8_t hap = 0; hap < 2; ++hap) {
        // ------------------------------------------------ references ------------------------------------------------
        std::vector<std::string> reference_sequences;
        std::vector<std::string> reference_names;
        reference_sequences.reserve(large[hap].size());
        reference_names.reserve(large[hap].size());
        for (size_t i = 0; i < large[hap].size(); ++i) {
            std::string sequence = get_path_sequence(large[hap][i]->vertices);
            if (sequence == "*") continue;
            reference_sequences.push_back(std::move(sequence));
            reference_names.push_back("large_" + std::to_string(i));
        }
        std::vector<const char*> sequences;
        std::vector<const char*> names;
        sequences.reserve(reference_sequences.size());
        names.reserve(reference_names.size());
        for (size_t i = 0; i < reference_sequences.size(); ++i) {
            sequences.push_back(reference_sequences[i].c_str());
            names.push_back(reference_names[i].c_str());
        }

        // ------------------------------------------------ queries ------------------------------------------------
        std::vector<Query> queries;
        queries.reserve(small_chains[hap].size());
        for (const Chain* chain : small_chains[hap]) {
            std::string sequence = get_path_sequence(chain->vertices);
            if (sequence == "*") sequence.clear();
            queries.push_back({chain, std::move(sequence)});
        }
        if (queries.empty()) continue;
        if (sequences.empty()) {
            low_quality[hap] = small_chains[hap];
            continue;
        }

        // ------------------------------------------------ Index ------------------------------------------------
        const int hpc = (index_options.flag & MM_I_HPC) ? 1 : 0;
        mm_idx_t* index = mm_idx_str(
            index_options.w, index_options.k, hpc, index_options.bucket_bits,
            static_cast<int>(sequences.size()), sequences.data(), names.data()
        );
        if (!index) {
            low_quality[hap] = small_chains[hap];
            continue;
        }
        mm_mapopt_t map_options = default_map_options;
        map_options.flag |= MM_F_CIGAR | MM_F_EQX;
        map_options.best_n = static_cast<short>(mm2_.best_n);
        map_options.zdrop = mm2_.zdrop;
        mm_mapopt_update(&map_options, index);

        // ------------------------------------------------ Map ------------------------------------------------
        ThreadPool pool(params_.threads);
        std::vector<std::future<bool>> futures;
        futures.reserve(queries.size());
        for (const Query& query : queries) {
            if (query.sequence.empty()) {
                futures.emplace_back();
                continue;
            }
            futures.emplace_back(pool.submit(
                sequence_contained, index, std::cref(map_options),
                std::cref(query.sequence), mm2_.min_mapq,
                mm2_.min_match, params_.dedup_similarity
            ));
        }
        // ------------------------------------------------ Collect low-quality chains ------------------------------------------------
        ProgressTracker progress(queries.size());
        for (size_t i = 0; i < futures.size(); ++i) {
            if (!futures[i].valid() || !futures[i].get()) low_quality[hap].push_back(queries[i].chain);
            progress.hit();
        }
        progress.finish();
        mm_idx_destroy(index);
        log_stream() << "  - hap" << static_cast<uint32_t>(hap + 1) << ": checked=" << queries.size() << ", contained=" << queries.size() - low_quality[hap].size() << ", low_quality=" << low_quality[hap].size() << '\n';
    }
    std::cerr << '\n';
    return low_quality;
}


// ================================================= Sequence, chain, and output helpers =================================================
void GfaGapfill::print_summary_(size_t selected) const {
    size_t path_count = 0;
    for (const auto& item : sample_chains_) {
        path_count += item.second[0].size() + item.second[1].size();
    }
    log_stream() << "Gap-fill summary:\n";
    log_stream() << "  - Samples: " << sample_chains_.size() << "\n";
    log_stream() << "  - Contig paths: " << fragments_.size() << "\n";
    log_stream() << "  - Confident gap bridges: " << selected << "\n";
    log_stream() << "  - Sample paths written: " << path_count << "\n\n";
}

// Select representative haplotypes, assemble sequence and write final GFA and reports.
void GfaGapfill::save_samples(
    const std::string& prefix,
    const std::string& command_line
) {
    log_stream() << "Writing gap-filled sample graphs ...\n";
    std::unordered_set<std::string> filenames;
    for (const auto& item : sample_chains_) {
        std::string name = safe_name(item.first);
        const std::string base = name;
        size_t suffix = 1;
        while (!filenames.insert(name).second) name = base + "_" + std::to_string(suffix++);
        std::vector<const Chain*> hap1;
        std::vector<const Chain*> hap2;
        hap1.reserve(item.second[0].size());
        hap2.reserve(item.second[1].size());
        for (const Chain& chain : item.second[0]) hap1.push_back(&chain);
        for (const Chain& chain : item.second[1]) hap2.push_back(&chain);
        save_haplotype_(item.first, 1, hap1, prefix + "." + name + ".hap1", command_line);
        save_haplotype_(item.first, 2, hap2, prefix + "." + name + ".hap2", command_line);
    }
    std::cerr << '\n';

    const PrimarySelection primary = select_primary_chains_();
    log_stream() << "Writing primary haplotype graphs ...\n";
    save_haplotype_("primary", 1, primary.primary[0], prefix + ".primary.hap1", command_line, true);
    save_haplotype_("primary", 2, primary.primary[1], prefix + ".primary.hap2", command_line, true);
    std::cerr << '\n';

    log_stream() << "Writing low-quality primary graphs ...\n";
    save_haplotype_("primary.lowQ", 1, primary.low_quality[0], prefix + ".primary.hap1.lowQ", command_line, true, "lq");
    save_haplotype_("primary.lowQ", 2, primary.low_quality[1], prefix + ".primary.hap2.lowQ", command_line, true, "lq");

    SAVE report(prefix + ".gapfill.tsv");
    report.save(
        "sample\thap\tcomponent\tcontig\tleft\tright\tbridge\tgap_interval\t"
        "phase\talign_identity\talign_coverage\tsample_support\tconfidence\n"
    );
    for (const GapRecord& record : records_) {
        std::ostringstream line;
        line << record.sample << '\t' << static_cast<uint32_t>(record.hap) << '\t'
             << record.component << '\t'
             << record.contig << '\t' << record.left << '\t' << record.right << '\t'
             << record.bridge << '\t'
             << record.gap_beg << '-' << record.gap_end << '\t'
             << std::fixed << std::setprecision(2);
        line << metric_pair(record.phase_left, record.phase_right) << '\t'
             << metric_pair(record.identity_left, record.identity_right) << '\t'
             << metric_pair(record.coverage_left, record.coverage_right) << '\t'
             << record.sample_support << '/' << record.spanning_samples << '\t'
             << record.confidence << '\n';
        report.save(line.str());
    }

    SAVE relocation_report(prefix + ".gapfill.relocations.tsv");
    relocation_report.save(
        "sample\thap\tsource_component\tsource_contig\tside\tstatus\trelocated_contig\t"
        "target_component\ttarget_backbone\ttarget_interval\tstrand\tp\tused_for_gap\n"
    );
    for (const RelocationRecord& record : relocations_) {
        std::ostringstream line;
        line << record.sample << '\t' << static_cast<uint32_t>(record.hap) << '\t'
             << record.source_component << '\t' << record.source << '\t' << record.side << '\t'
             << record.status << '\t'
             << (record.relocated_contig.empty() ? "." : record.relocated_contig) << '\t';
        if (record.target_component == UINT32_MAX) line << ".\t.\t.\t.\t";
        else line << record.target_component << '\t' << record.target_backbone << '\t' << record.target_beg << '-' << record.target_end << '\t' << record.strand << '\t';
        line << std::fixed << std::setprecision(2) << record.probability << '\t' << (record.used_for_gap ? 1 : 0) << '\n';
        relocation_report.save(line.str());
    }
    log_stream() << "  - Samples written: " << sample_chains_.size() << "\n\n";
}

void GfaGapfill::save_haplotype_(
    const std::string& label,
    uint8_t hap,
    const std::vector<const Chain*>& chains,
    const std::string& prefix,
    const std::string& command_line,
    bool primary,
    const char* primary_kind
) {
    // ------------------------------------------------ Open files ------------------------------------------------
    const std::string suffix = ".gapfill";
    SAVE with_seq(prefix + suffix + ".gfa");
    SAVE no_seq(prefix + suffix + ".noseq.gfa");

    // ------------------------------------------------ GFA header ------------------------------------------------
    std::string header = "H\tVN:Z:1.0\nH\tTS:Z:liftasm\n";
    if (!command_line.empty()) {
        std::string command = command_line;
        for (char& c : command) {
            if (c == '\t' || c == '\n' || c == '\r') c = ' ';
        }
        header += "H\tCL:Z:" + command + '\n';
    }
    header += "H\tCO:Z:  LN:i   segment length\n";
    header += "H\tCO:Z:  SN:Z   node source/sample name\n";
    header += "H\tCO:Z:  CP:i   graph component\n";
    header += "H\tCO:Z:  GF:Z   gap-filled intervals (0-based, half-open)\n";
    header += "H\tCO:Z:  SC:Z   target pieces as contig@start-end+/- (source coordinates; 0-based, half-open)\n";
    header += "H\tCO:Z:  GS:Z   gap-source pieces in GF order, using the same coordinate format\n";
    header += "H\tCO:Z:  MS:Z   source contig of an unused misassembly\n";
    header += "H\tCO:Z:  ME:Z   source end of an unused misassembly\n";
    header += "H\tCO:Z:  MR:Z   reason the misassembly was not used\n";
    with_seq.save(header);
    no_seq.save(header);

    // ------------------------------------------------ Reserve names of unchanged chains ------------------------------------------------
    std::unordered_set<std::string> output_names;
    for (const Chain* chain_ptr : chains) {
        const Chain& chain = *chain_ptr;
        if (!chain.gaps.empty() || chain.source_fragments.empty()) continue;
        const Fragment& source = fragments_[chain.source_fragments.front()];
        output_names.insert(output_contig_name(paths_[source.path_id].name, chain.sample));
    }

    // ------------------------------------------------ Process chains ------------------------------------------------
    size_t filled_number = 0;
    size_t primary_number = 0;
    for (const Chain* chain_ptr : chains) {
        const Chain& chain = *chain_ptr;
        if (chain.source_fragments.empty()) continue;
        const bool gap_filled = !chain.gaps.empty();
        const Fragment& first_source = fragments_[chain.source_fragments.front()];

        // ------------------------------------------------ Choose output name ------------------------------------------------
        std::string name = output_contig_name(paths_[first_source.path_id].name, chain.sample);
        if (primary) {
            name = numbered_contig_name(hap, primary_kind, ++primary_number);
        } else if (gap_filled) {
            do {
                name = numbered_contig_name(hap, "gf", ++filled_number);
            } while (!output_names.insert(name).second);
        }

        // ------------------------------------------------ Choose graph path ------------------------------------------------
        std::vector<uint32_t> original_vertices;
        const std::vector<uint32_t>* assembled_vertices = &chain.vertices;
        if (!gap_filled) {
            const GfaPath& path = paths_[first_source.path_id];
            original_vertices.reserve(path.segments.size());
            for (const PathSegment& segment : path.segments) {
                original_vertices.push_back(Vertex::make_vertex(
                    static_cast<uint32_t>(segment.node_id), segment.is_reverse
                ));
            }
            assembled_vertices = &original_vertices;
        }

        // ------------------------------------------------ Prepare sequence state ------------------------------------------------
        uint64_t assembled_length = 0;
        bool sequence_available = true;
        std::string sequence;
        std::vector<std::pair<uint64_t, uint64_t>> intervals;

        // ------------------------------------------------ Build chain sequence ------------------------------------------------
        if (gap_filled) {
            sequence.reserve(chain_length_(chain));
            for (const Chain::Part& chain_part : chain.parts) {
                const uint64_t part_length = chain_part.end - chain_part.begin;
                if (chain_part.filled) {
                    intervals.emplace_back(assembled_length, assembled_length + part_length);
                }
                if (chain_part.begin < chain_part.end) {
                    std::string part = fragment_subsequence_(
                        fragments_[chain_part.fragment], chain_part.begin, chain_part.end
                    );
                    if (part.empty()) sequence_available = false;
                    for (char& base : part) {
                        base = chain_part.filled ?
                            static_cast<char>(std::tolower(static_cast<unsigned char>(base))) :
                            static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
                    }
                    if (sequence_available) sequence += part;
                }
                assembled_length += part_length;
            }
        } else {
            // ------------------------------------------------ Build unchanged sequence ------------------------------------------------
            for (size_t i = 0; i < assembled_vertices->size(); ++i) {
                const uint32_t vertex = (*assembled_vertices)[i];
                const uint32_t sid = Vertex::get_segment_id(vertex);
                const uint32_t overlap = i == 0 ? 0 : get_edge_ow((*assembled_vertices)[i - 1], vertex);
                const uint32_t added = nodes_[sid].length > overlap ? nodes_[sid].length - overlap : 0;
                std::string part = get_oriented_sequence(Vertex(vertex), overlap);
                if (part.size() != added) sequence_available = false;
                for (char& base : part) {
                    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
                }
                if (sequence_available) sequence += part;
                assembled_length += added;
            }
        }
        if (!sequence_available) sequence.clear();

        // ------------------------------------------------ Build segment tags ------------------------------------------------
        std::ostringstream tags;
        tags << "\tLN:i:" << assembled_length << "\tSN:Z:" << chain.sample
             << "\tCP:i:" << chain.component;
        if (!intervals.empty()) {
            tags << "\tGF:Z:";
            for (size_t i = 0; i < intervals.size(); ++i) {
                if (i) tags << ',';
                tags << intervals[i].first << '-' << intervals[i].second;
            }
        }
        if (gap_filled) {
            tags << "\tSC:Z:";
            bool first = true;
            for (const Chain::Part& part : chain.parts) {
                if (part.filled) continue;
                if (!first) tags << ',';
                first = false;
                const Fragment& source = fragments_[part.fragment];
                uint64_t begin = part.begin, end = part.end;
                if (source.reverse) {
                    begin = source.length - part.end;
                    end = source.length - part.begin;
                }
                tags << output_contig_name(paths_[source.path_id].name, chain.sample)
                     << '@' << begin << '-' << end << (source.reverse ? '-' : '+');
            }
            tags << "\tGS:Z:";
            first = true;
            for (const Chain::Part& part : chain.parts) {
                if (!part.filled) continue;
                if (!first) tags << ',';
                first = false;
                const Fragment& source = fragments_[part.fragment];
                uint64_t begin = part.begin, end = part.end;
                if (source.reverse) {
                    begin = source.length - part.end;
                    end = source.length - part.begin;
                }
                tags << paths_[source.path_id].name
                     << '@' << begin << '-' << end << (source.reverse ? '-' : '+');
            }
        }

        // ------------------------------------------------ Write GFA segment records ------------------------------------------------
        with_seq.save("S\t" + name + '\t');
        with_seq.save(sequence_available ? sequence : "*");
        with_seq.save(tags.str() + '\n');
        no_seq.save("S\t" + name + "\t*" + tags.str() + '\n');

        // ------------------------------------------------ Add gap report records ------------------------------------------------
        for (size_t i = 0; !primary && i < intervals.size() && i < chain.gaps.size(); ++i) {
            const Candidate& candidate = chain.gaps[i];
            records_.push_back({
                chain.sample,
                name,
                paths_[fragments_[candidate.left].path_id].name,
                paths_[fragments_[candidate.right].path_id].name,
                paths_[fragments_[candidate.bridge].path_id].name,
                hap,
                chain.component,
                intervals[i].first,
                intervals[i].second,
                candidate.left_boundary.phase.similarity,
                candidate.right_boundary.phase.similarity,
                candidate.left_boundary.identity,
                candidate.left_boundary.coverage,
                candidate.right_boundary.identity,
                candidate.right_boundary.coverage,
                candidate.sample_support,
                candidate.spanning_samples,
                candidate.confidence
            });
        }
        GfaGapfillDebugger::chain(
            name, hap, chain.component, chain.gaps.size() + 1,
            intervals.size(), assembled_length
        );
    }

    log_stream() << "  - " << label << ".hap" << static_cast<uint32_t>(hap) << ": " << chains.size() << " contigs\n";
}


// ================================================= helpers =================================================
std::string GfaGapfill::fragment_sequence_(
    const Fragment& fragment,
    uint32_t begin,
    uint32_t end
) const {
    if (begin >= end || end > fragment.vertices.size()) return {};

    std::string sequence;
    sequence.reserve(fragment.path_bp[end] - fragment.path_bp[begin]);
    for (uint32_t i = begin; i < end; ++i) {
        const uint32_t overlap = i == begin ? 0 : get_edge_ow(fragment.vertices[i - 1], fragment.vertices[i]);
        std::string part = get_oriented_sequence(Vertex(fragment.vertices[i]), overlap);
        if (part.empty()) return {};
        sequence += part;
    }
    return sequence;
}

std::string GfaGapfill::fragment_subsequence_(
    const Fragment& fragment,
    uint64_t begin,
    uint64_t end
) const {
    if (begin >= end || end > fragment.length) return {};

    const auto first = std::upper_bound(
        fragment.path_bp.begin(), fragment.path_bp.end(), begin
    );
    uint32_t pos = first == fragment.path_bp.begin() ? 0 : static_cast<uint32_t>(first - fragment.path_bp.begin() - 1);
    std::string sequence;
    sequence.reserve(end - begin);
    while (pos < fragment.vertices.size() && fragment.path_bp[pos] < end) {
        const uint64_t node_begin = fragment.path_bp[pos];
        const uint64_t node_end = fragment.path_bp[pos + 1];
        std::string part = get_oriented_sequence(Vertex(fragment.vertices[pos]), 0);
        if (part.empty() || part.size() != node_end - node_begin) return {};
        const uint64_t from = std::max(begin, node_begin) - node_begin;
        const uint64_t to = std::min(end, node_end) - node_begin;
        sequence.append(part, from, to - from);
        ++pos;
    }
    return sequence.size() == end - begin ? sequence : std::string{};
}

uint64_t GfaGapfill::chain_length_(const Chain& chain) const {
    if (!chain.parts.empty()) {
        uint64_t length = 0;
        for (const Chain::Part& part : chain.parts) length += part.end - part.begin;
        return length;
    }
    uint64_t length = 0;
    for (size_t i = 0; i < chain.vertices.size(); ++i) {
        const uint32_t node_length = nodes_[Vertex::get_segment_id(chain.vertices[i])].length;
        const uint32_t overlap = i == 0 ? 0 : get_edge_ow(chain.vertices[i - 1], chain.vertices[i]);
        length += node_length > overlap ? node_length - overlap : 0;
    }
    return length;
}

uint64_t GfaGapfill::ng50_(std::vector<uint64_t> lengths, uint64_t component_bp) {
    if (component_bp == 0) return 0;
    std::sort(lengths.begin(), lengths.end(), std::greater<uint64_t>());
    const uint64_t threshold = component_bp / 2 + component_bp % 2;
    uint64_t covered = 0;
    for (uint64_t length : lengths) {
        covered += length;
        if (covered >= threshold) return length;
    }
    return 0;
}


/* =================================================================================================================
 * GAP BOUNDARY START
 * ================================================================================================================= */
// Build graph anchors and local evidence once, then refine selected cuts with mm2
bool GfaGapfillBoundary::prepare(
    uint32_t left_id, uint32_t right_id, uint32_t bridge_id,
    Boundary& left_boundary, Boundary& right_boundary
) const {
    if (!find_anchors_(left_id, right_id, bridge_id, left_boundary, right_boundary)) return false;
    const GfaGapfill::Fragment& left = gapfill_.fragments_[left_id];
    const GfaGapfill::Fragment& right = gapfill_.fragments_[right_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];

    // Shared-node cuts are safe fallbacks when local alignment is unavailable.
    left_boundary.target_cut_bp = left.path_bp[left_boundary.target_pos + 1];
    left_boundary.bridge_cut_bp = bridge.path_bp[left_boundary.bridge_pos + 1];
    right_boundary.bridge_cut_bp = bridge.path_bp[right_boundary.bridge_pos];
    right_boundary.target_cut_bp = right.path_bp[right_boundary.target_pos];

    // Cache all local evidence and alignment coordinates once per candidate.
    left_boundary.phase = phase_similarity_(
        left_id, left_boundary.target_pos,
        bridge_id, left_boundary.bridge_pos, true
    );
    right_boundary.phase = phase_similarity_(
        right_id, right_boundary.target_pos,
        bridge_id, right_boundary.bridge_pos, false
    );
    left_boundary.minimizer = minimizer_similarity_(
        left_id, left_boundary.target_pos,
        bridge_id, left_boundary.bridge_pos, true, 1
    );
    right_boundary.minimizer = minimizer_similarity_(
        right_id, right_boundary.target_pos,
        bridge_id, right_boundary.bridge_pos, false, 1
    );

    prepare_alignment_window_(left_id, bridge_id, left_boundary, right_boundary, true, left_boundary);
    prepare_alignment_window_(right_id, bridge_id, left_boundary, right_boundary, false, right_boundary);
    return true;
}

bool GfaGapfillBoundary::spans(
    uint32_t left, uint32_t right, uint32_t bridge
) const {
    Boundary left_boundary, right_boundary;
    return find_anchors_(left, right, bridge, left_boundary, right_boundary);
}

GfaGapfillBoundary::Boundary GfaGapfillBoundary::inspect(
    uint32_t target,
    uint32_t target_pos,
    uint32_t bridge,
    uint32_t bridge_pos,
    bool before,
    bool expand_minimizer
) const {
    Boundary boundary;
    boundary.target_pos = target_pos;
    boundary.bridge_pos = bridge_pos;
    boundary.phase = phase_similarity_(
        target, target_pos, bridge, bridge_pos, before
    );
    boundary.minimizer = minimizer_similarity_(
        target, target_pos, bridge, bridge_pos, before,
        expand_minimizer ? boundary.phase.windows : 1
    );
    return boundary;
}

bool GfaGapfillBoundary::refine(
    uint32_t left_id, uint32_t right_id, uint32_t bridge_id,
    Boundary& left_boundary, Boundary& right_boundary
) const {
    const GfaGapfill::Fragment& left = gapfill_.fragments_[left_id];
    const GfaGapfill::Fragment& right = gapfill_.fragments_[right_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];
    std::array<Alignment, 2> alignments;
    bool mapped_left = false;
    bool mapped_right = false;
    if (left_boundary.alignable) {
        mapped_left = align_(
            gapfill_.fragment_subsequence_(
                bridge, left_boundary.bridge_begin, left_boundary.bridge_end
            ),
            gapfill_.fragment_subsequence_(
                left, left_boundary.target_begin, left_boundary.target_end
            ),
            true, alignments[0]
        );
    }
    if (right_boundary.alignable) {
        mapped_right = align_(
            gapfill_.fragment_subsequence_(
                bridge, right_boundary.bridge_begin, right_boundary.bridge_end
            ),
            gapfill_.fragment_subsequence_(
                right, right_boundary.target_begin, right_boundary.target_end
            ),
            false, alignments[1]
        );
    }
    if (left_boundary.alignable && alignments[0].hits > 0) {
        left_boundary.identity = alignments[0].identity;
        left_boundary.coverage = alignments[0].coverage;
    }
    if (right_boundary.alignable && alignments[1].hits > 0) {
        right_boundary.identity = alignments[1].identity;
        right_boundary.coverage = alignments[1].coverage;
    }

    const GfaGapfillDebugger::Alignment left_result{
        left_boundary.target_begin, left_boundary.target_end,
        left_boundary.bridge_begin, left_boundary.bridge_end,
        alignments[0].query_cut, alignments[0].reference_cut,
        alignments[0].query_begin, alignments[0].query_end,
        alignments[0].reference_begin, alignments[0].reference_end,
        alignments[0].hits, alignments[0].mapq,
        alignments[0].identity, alignments[0].coverage,
        alignments[0].accepted, alignments[0].reverse
    };
    const GfaGapfillDebugger::Alignment right_result{
        right_boundary.target_begin, right_boundary.target_end,
        right_boundary.bridge_begin, right_boundary.bridge_end,
        alignments[1].query_cut, alignments[1].reference_cut,
        alignments[1].query_begin, alignments[1].query_end,
        alignments[1].reference_begin, alignments[1].reference_end,
        alignments[1].hits, alignments[1].mapq,
        alignments[1].identity, alignments[1].coverage,
        alignments[1].accepted, alignments[1].reverse
    };
    GfaGapfillDebugger::boundary_alignment(
        {gapfill_.paths_[left.path_id].name, left.length, left.reverse},
        {gapfill_.paths_[right.path_id].name, right.length, right.reverse},
        {gapfill_.paths_[bridge.path_id].name, bridge.length, bridge.reverse},
        left_result,
        right_result
    );

    /**************************** may have bugs ********************************/
    /**
     * @date 2026-08-13
     * @version v0.1.3-r21
     * 
     * @note
     *   left_cut = std::max(old_left_cut, refined_left_cut);
     *   right_cut = std::min(old_right_cut, refined_right_cut);
     */
    uint64_t left_cut = left_boundary.target_cut_bp;
    uint64_t bridge_left_cut = left_boundary.bridge_cut_bp;
    uint64_t bridge_right_cut = right_boundary.bridge_cut_bp;
    uint64_t right_cut = right_boundary.target_cut_bp;
    if (mapped_left) {
        left_cut = left_boundary.target_begin + alignments[0].query_cut;
        bridge_left_cut = left_boundary.bridge_begin + alignments[0].reference_cut;
    }
    if (mapped_right) {
        bridge_right_cut = right_boundary.bridge_begin + alignments[1].reference_cut;
        right_cut = right_boundary.target_begin + alignments[1].query_cut;
    }
    if (bridge_left_cut > bridge_right_cut) {
        if (!mapped_left || !mapped_right || !splice_overlap_(
                alignments[0], alignments[1], left_boundary, right_boundary,
                left_cut, right_cut, bridge_left_cut, bridge_right_cut
            )) return false;
    }
    left_boundary.target_cut_bp = left_cut;
    left_boundary.bridge_cut_bp = bridge_left_cut;
    right_boundary.bridge_cut_bp = bridge_right_cut;
    right_boundary.target_cut_bp = right_cut;
    /**************************** may have bugs ********************************/

    return mapped_left || mapped_right;
}

// Anchor discovery keeps graph order and rejects ambiguous or non-unique bridge intervals.
bool GfaGapfillBoundary::find_anchors_(
    uint32_t left_id, uint32_t right_id, uint32_t bridge_id,
    Boundary& left_boundary, Boundary& right_boundary
) const {
    const GfaGapfill::Fragment& left = gapfill_.fragments_[left_id];
    const GfaGapfill::Fragment& right = gapfill_.fragments_[right_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];
    GfaGapfill::Candidate candidate;

    // ------------------------------------------------ Check direct target overlap ------------------------------------------------
    // Reject target pairs whose shared nodes imply a large or incorrectly ordered overlap.
    uint64_t left_sum = 0, right_sum = 0, direct_bp = 0;
    uint32_t direct_matches = 0;
    for (uint32_t left_pos = 0; left_pos < left.vertices.size(); ++left_pos) {
        const uint32_t right_pos = find_position_(right_id, left.vertices[left_pos]);
        if (right_pos == UINT32_MAX) continue;
        left_sum += left_pos;
        right_sum += right_pos;
        direct_bp += gapfill_.nodes_[Vertex::get_segment_id(left.vertices[left_pos])].length;
        ++direct_matches;
    }
    if (direct_matches > 0) {
        const uint64_t shorter_bp = std::min(left.length, right.length);
        if (shorter_bp == 0 || static_cast<double>(direct_bp) / shorter_bp > gapfill_.params_.max_target_overlap) {
            return false;
        }
        if (2 * left_sum + direct_matches <= direct_matches * left.vertices.size() || 2 * right_sum + direct_matches >= direct_matches * right.vertices.size()) {
            return false;
        }
    }

    // ------------------------------------------------ Collect shared nodes ------------------------------------------------
    std::vector<std::pair<uint32_t, uint32_t>> left_matches;  // {bridge_pos, left_pos}
    std::vector<std::pair<uint32_t, uint32_t>> right_matches;  // {bridge_pos, right_pos}}
    left_matches.reserve(std::min(left.vertices.size(), bridge.vertices.size()));
    right_matches.reserve(std::min(right.vertices.size(), bridge.vertices.size()));

    for (uint32_t bridge_pos = 0; bridge_pos < bridge.vertices.size(); ++bridge_pos) {
        const uint32_t vertex = bridge.vertices[bridge_pos];
        const uint32_t left_pos = find_position_(left_id, vertex);
        if (left_pos != UINT32_MAX) left_matches.emplace_back(bridge_pos, left_pos);
        const uint32_t right_pos = find_position_(right_id, vertex);
        if (right_pos != UINT32_MAX) right_matches.emplace_back(bridge_pos, right_pos);
    }
    if (left_matches.empty() || right_matches.empty()) return false;

    // ------------------------------------------------ Choose the best junctions ------------------------------------------------
    std::vector<uint32_t> suffix(right_matches.size());
    suffix.back() = static_cast<uint32_t>(right_matches.size() - 1);
    for (size_t i = right_matches.size() - 1; i > 0; --i) {
        const uint32_t current = static_cast<uint32_t>(i - 1);
        const uint32_t best = suffix[i];
        suffix[current] = right_matches[current].second < right_matches[best].second ? current : best;
    }

    // Trend to select the junctions that retain the most target sequence, and in case of ties, the smallest bridge span
    uint64_t best_retained = 0;
    uint32_t best_span = UINT32_MAX;
    bool found = false;
    for (const auto& left_match : left_matches) {
        const auto it = std::upper_bound(
            right_matches.begin(), right_matches.end(),
            std::pair<uint32_t, uint32_t>{left_match.first, UINT32_MAX}
        );
        if (it == right_matches.end()) continue;

        const size_t first = static_cast<size_t>(it - right_matches.begin());
        const auto& right_match = right_matches[suffix[first]];
        if (bridge.rank_drops[right_match.first] != bridge.rank_drops[left_match.first]) continue;  // left and right matches must be on the same phase region
        const uint64_t retained = left.path_bp[left_match.second + 1] + right.length - right.path_bp[right_match.second];
        const uint32_t span = right_match.first - left_match.first;
        if (!found || retained > best_retained || (retained == best_retained && span < best_span)) {
            candidate.left_boundary.target_pos = left_match.second;
            candidate.left_boundary.bridge_pos = left_match.first;
            candidate.right_boundary.bridge_pos = right_match.first;
            candidate.right_boundary.target_pos = right_match.second;
            best_retained = retained;
            best_span = span;
            found = true;
        }
    }
    if (!found) return false;

    // ------------------------------------------------ Try a shared crossover junction ------------------------------------------------
    GfaGapfill::Candidate normal = candidate;
    normal.left_boundary.junction_target_pos = normal.left_boundary.target_pos;
    normal.left_boundary.junction_bridge_pos = normal.left_boundary.bridge_pos;
    normal.right_boundary.junction_bridge_pos = normal.right_boundary.bridge_pos;
    normal.right_boundary.junction_target_pos = normal.right_boundary.target_pos;

    GfaGapfill::Candidate crossover = normal;
    uint64_t crossover_retained = 0;
    bool have_crossover = false;
    // A shared three-way node is preferred when the two target contigs overlap.
    for (const auto& match : left_matches) {
        const uint32_t vertex = bridge.vertices[match.first];
        const uint32_t right_pos = find_position_(right_id, vertex);
        if (right_pos == UINT32_MAX) continue;
        const uint64_t retained = left.path_bp[match.second] + right.length - right.path_bp[right_pos];
        const bool prefer_right_intersection = right.length > left.length && match.first < crossover.left_boundary.bridge_pos;
        const bool prefer_left_intersection = left.length > right.length && match.first > crossover.left_boundary.bridge_pos;
        if (!have_crossover || retained > crossover_retained || (retained == crossover_retained && (prefer_right_intersection || prefer_left_intersection))) {
            crossover.left_boundary.target_pos = match.second;
            crossover.left_boundary.bridge_pos = match.first;
            crossover.right_boundary.bridge_pos = match.first;
            crossover.right_boundary.target_pos = right_pos;
            crossover.left_boundary.junction_target_pos = match.second;
            crossover.left_boundary.junction_bridge_pos = match.first;
            crossover.right_boundary.junction_bridge_pos = match.first;
            crossover.right_boundary.junction_target_pos = right_pos;
            crossover_retained = retained;
            have_crossover = true;
        }
    }

    // ------------------------------------------------ Move boundaries to the trusted regions ------------------------------------------------
    if (have_crossover && adjust_boundaries_(left_id, right_id, bridge_id, left_matches, right_matches, crossover.left_boundary, crossover.right_boundary)) {
        candidate = crossover;
    } else {
        candidate = normal;
        if (!adjust_boundaries_(left_id, right_id, bridge_id, left_matches, right_matches, candidate.left_boundary, candidate.right_boundary)) return false;
    }

    if (!is_cycle_free_(left_id, right_id, bridge_id, candidate.left_boundary, candidate.right_boundary)) return false;

    // ------------------------------------------------ Move junctions outside the shared overlap ------------------------------------------------
    // Misassembly checks start outside the complete target/bridge overlap: after the last shared node on the left and before the first on the right.
    const uint32_t left_bridge_rank = bridge.rank_drops[candidate.left_boundary.junction_bridge_pos];
    for (const auto& match : left_matches) {
        if (bridge.rank_drops[match.first] != left_bridge_rank || match.first < candidate.left_boundary.junction_bridge_pos || match.second < candidate.left_boundary.junction_target_pos) continue;
        if (match.second > candidate.left_boundary.junction_target_pos || (match.second == candidate.left_boundary.junction_target_pos && match.first > candidate.left_boundary.junction_bridge_pos)) {
            candidate.left_boundary.junction_bridge_pos = match.first;
            candidate.left_boundary.junction_target_pos = match.second;
        }
    }

    const uint32_t right_bridge_rank = bridge.rank_drops[candidate.right_boundary.junction_bridge_pos];
    for (const auto& match : right_matches) {
        if (bridge.rank_drops[match.first] != right_bridge_rank || match.first > candidate.right_boundary.junction_bridge_pos || match.second > candidate.right_boundary.junction_target_pos) continue;
        if (match.second < candidate.right_boundary.junction_target_pos || (match.second == candidate.right_boundary.junction_target_pos && match.first < candidate.right_boundary.junction_bridge_pos)) {
            candidate.right_boundary.junction_bridge_pos = match.first;
            candidate.right_boundary.junction_target_pos = match.second;
        }
    }

    // ------------------------------------------------ Return the selected boundaries ------------------------------------------------
    left_boundary = candidate.left_boundary;
    right_boundary = candidate.right_boundary;
    return true;
}

bool GfaGapfillBoundary::adjust_boundaries_(
    uint32_t left_id, uint32_t right_id, uint32_t bridge_id,
    const std::vector<std::pair<uint32_t, uint32_t>>& left_matches,
    const std::vector<std::pair<uint32_t, uint32_t>>& right_matches,
    Boundary& left_boundary, Boundary& right_boundary
) const {
    const GfaGapfill::Fragment& left = gapfill_.fragments_[left_id];
    const GfaGapfill::Fragment& right = gapfill_.fragments_[right_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];
    // Move each junction past the unreliable contig end to the nearest shared node.
    const uint64_t left_junction = left.path_bp[left_boundary.junction_target_pos + 1];
    const uint64_t left_trusted_end = left.length > gapfill_.params_.end_skip_bp ? left.length - gapfill_.params_.end_skip_bp : 0;
    const uint64_t left_goal = std::min(left_junction, left_trusted_end);
    uint64_t best_left_distance = UINT64_MAX;
    uint32_t best_left_bridge = left_boundary.junction_bridge_pos;
    uint32_t best_left_pos = left_boundary.junction_target_pos;
    bool have_left = false;
    for (const auto& match : left_matches) {
        if (match.first > left_boundary.junction_bridge_pos || match.second > left_boundary.junction_target_pos) continue;  // in the left
        if (bridge.rank_drops[match.first] != bridge.rank_drops[left_boundary.junction_bridge_pos]) continue;  // need to be on the same phase region
        const uint64_t position = left.path_bp[match.second + 1];
        if (position > left_goal) continue;
        const uint64_t distance = left_goal - position;
        if (distance < best_left_distance) {
            best_left_distance = distance;
            best_left_bridge = match.first;
            best_left_pos = match.second;
            have_left = true;
        }
    }

    const uint64_t right_junction = right.path_bp[right_boundary.junction_target_pos];
    const uint64_t right_goal = std::max(right_junction, gapfill_.params_.end_skip_bp);
    uint64_t best_right_distance = UINT64_MAX;
    uint32_t best_right_bridge = right_boundary.junction_bridge_pos;
    uint32_t best_right_pos = right_boundary.junction_target_pos;
    bool have_right = false;
    for (const auto& match : right_matches) {
        if (match.first < right_boundary.junction_bridge_pos || match.second < right_boundary.junction_target_pos) continue;  // in the right
        if (bridge.rank_drops[match.first] != bridge.rank_drops[right_boundary.junction_bridge_pos]) continue;  // need to be on the same phase region
        const uint64_t position = right.path_bp[match.second];
        if (position < right_goal) continue;
        const uint64_t distance = position - right_goal;
        if (distance < best_right_distance) {
            best_right_distance = distance;
            best_right_bridge = match.first;
            best_right_pos = match.second;
            have_right = true;
        }
    }

    if (!have_left || !have_right) return false;
    left_boundary.target_pos = best_left_pos;
    left_boundary.bridge_pos = best_left_bridge;
    right_boundary.bridge_pos = best_right_bridge;
    right_boundary.target_pos = best_right_pos;
    return right_boundary.bridge_pos > left_boundary.bridge_pos;
}

bool GfaGapfillBoundary::is_cycle_free_(
    uint32_t left_id, uint32_t right_id, uint32_t bridge_id,
    const Boundary& left_boundary, const Boundary& right_boundary
) const {
    const GfaGapfill::Fragment& left = gapfill_.fragments_[left_id];
    const GfaGapfill::Fragment& right = gapfill_.fragments_[right_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];
    std::unordered_set<uint32_t> seen;
    for (uint32_t i = 0; i <= left_boundary.target_pos; ++i) {
        if (!seen.insert(Vertex::get_segment_id(left.vertices[i])).second) return false;
    }
    for (uint32_t i = left_boundary.bridge_pos + 1; i < right_boundary.bridge_pos; ++i) {
        if (!seen.insert(Vertex::get_segment_id(bridge.vertices[i])).second) return false;
    }
    for (uint32_t i = right_boundary.target_pos; i < right.vertices.size(); ++i) {
        if (!seen.insert(Vertex::get_segment_id(right.vertices[i])).second) return false;
    }
    return true;
}

// Local phase and minimizer evidence uses the trusted sequence inside each gap boundary.
uint64_t GfaGapfillBoundary::phase_available_(
    uint32_t fragment_id,
    uint32_t pos,
    bool before
) const {
    const GfaGapfill::Fragment& fragment = gapfill_.fragments_[fragment_id];
    if (before) {
        return std::min(fragment.path_bp[pos + 1], fragment.length > gapfill_.params_.end_skip_bp ? fragment.length - gapfill_.params_.end_skip_bp : 0);
    }
    const uint64_t begin = std::min(std::max(fragment.path_bp[pos], gapfill_.params_.end_skip_bp), fragment.length);
    return fragment.length - begin;
}

std::pair<uint32_t, uint32_t> GfaGapfillBoundary::phase_region_(
    uint32_t fragment_id,
    uint32_t pos,
    bool before,
    uint64_t windows
) const {
    const GfaGapfill::Fragment& fragment = gapfill_.fragments_[fragment_id];
    const uint64_t window = gapfill_.params_.phase_window_bp;
    const uint64_t available = phase_available_(fragment_id, pos, before);
    const uint64_t max_windows = available / window + (available % window != 0);
    const uint64_t span = windows >= max_windows ? available : window * windows;
    uint64_t begin_bp = 0;
    uint64_t end_bp = 0;
    if (before) {
        end_bp = std::min(
            fragment.path_bp[pos + 1],
            fragment.length > gapfill_.params_.end_skip_bp ? fragment.length - gapfill_.params_.end_skip_bp : 0
        );
        begin_bp = end_bp - span;
    } else {
        begin_bp = std::max(fragment.path_bp[pos], gapfill_.params_.end_skip_bp);
        begin_bp = std::min(begin_bp, fragment.length);
        end_bp = begin_bp + span;
    }
    if (end_bp <= begin_bp) return {0, 0};

    const auto begin_it = std::lower_bound(fragment.path_bp.begin(), fragment.path_bp.end(), begin_bp);
    const auto end_it = std::upper_bound(fragment.path_bp.begin(), fragment.path_bp.end(), end_bp);
    const uint32_t begin = std::min(
        static_cast<uint32_t>(begin_it - fragment.path_bp.begin()),
        static_cast<uint32_t>(fragment.vertices.size())
    );
    const uint32_t end = end_it == fragment.path_bp.begin() ? 0 : std::min(
        static_cast<uint32_t>(end_it - fragment.path_bp.begin() - 1),
        static_cast<uint32_t>(fragment.vertices.size())
    );
    return end > begin ? std::make_pair(begin, end) : std::make_pair(0u, 0u);
}

GfaGapfillBoundary::PhaseSupport GfaGapfillBoundary::phase_similarity_(
    uint32_t target_id,
    uint32_t target_pos,
    uint32_t bridge_id,
    uint32_t bridge_pos,
    bool before
) const {
    const GfaGapfill::Fragment& target = gapfill_.fragments_[target_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];
    PhaseSupport support;
    const uint64_t window = gapfill_.params_.phase_window_bp;
    const uint64_t target_available = phase_available_(target_id, target_pos, before);
    const uint64_t bridge_available = phase_available_(bridge_id, bridge_pos, before);
    const uint64_t target_windows = target_available / window + (target_available % window != 0);
    const uint64_t bridge_windows = bridge_available / window + (bridge_available % window != 0);
    uint64_t windows = 1;
    support.windows = windows;
    auto target_region = phase_region_(target_id, target_pos, before, windows);
    auto bridge_region = phase_region_(bridge_id, bridge_pos, before, windows);
    uint64_t target_bp = target.bubble_bp[target_region.second] - target.bubble_bp[target_region.first];
    uint64_t bridge_bp = bridge.bubble_bp[bridge_region.second] - bridge.bubble_bp[bridge_region.first];
    // Expand by whole phase windows until both paths contain enough bubble evidence.
    while (target_bp < gapfill_.params_.phase_min_bp || bridge_bp < gapfill_.params_.phase_min_bp) {
        if ((target_bp < gapfill_.params_.phase_min_bp && windows >= target_windows) || (bridge_bp < gapfill_.params_.phase_min_bp && windows >= bridge_windows)) {
            support.target_bp = target_bp;
            support.bridge_bp = bridge_bp;
            support.windows = windows;
            return support;
        }
        ++windows;
        support.windows = windows;
        target_region = phase_region_(target_id, target_pos, before, windows);
        bridge_region = phase_region_(bridge_id, bridge_pos, before, windows);
        target_bp = target.bubble_bp[target_region.second] - target.bubble_bp[target_region.first];
        bridge_bp = bridge.bubble_bp[bridge_region.second] - bridge.bubble_bp[bridge_region.first];
    }
    support.target_bp = target_bp;
    support.bridge_bp = bridge_bp;

    const GfaGapfill::Fragment* query = &target;
    const GfaGapfill::Fragment* reference = &bridge;
    auto query_region = target_region;
    auto reference_region = bridge_region;
    const auto target_begin = std::lower_bound(
        target.bubble_positions.begin(), target.bubble_positions.end(), target_region.first
    );
    const auto target_end = std::lower_bound(
        target.bubble_positions.begin(), target.bubble_positions.end(), target_region.second
    );
    const auto bridge_begin = std::lower_bound(
        bridge.bubble_positions.begin(), bridge.bubble_positions.end(), bridge_region.first
    );
    const auto bridge_end = std::lower_bound(
        bridge.bubble_positions.begin(), bridge.bubble_positions.end(), bridge_region.second
    );
    if (bridge_end - bridge_begin < target_end - target_begin) {
        query = &bridge;
        reference = &target;
        query_region = bridge_region;
        reference_region = target_region;
    }

    const auto query_begin = std::lower_bound(
        query->bubble_positions.begin(), query->bubble_positions.end(), query_region.first
    );
    const auto query_end = std::lower_bound(
        query->bubble_positions.begin(), query->bubble_positions.end(), query_region.second
    );
    for (auto it = query_begin; it != query_end; ++it) {
        const uint32_t pos = *it;
        const uint32_t reference_id = reference == &target ? target_id : bridge_id;
        const uint32_t reference_pos = find_position_(reference_id, query->vertices[pos]);
        if (reference_pos >= reference_region.first && reference_pos < reference_region.second) {
            support.shared_bp += gapfill_.nodes_[Vertex::get_segment_id(query->vertices[pos])].length;
        }
    }
    const uint64_t compared_bp = std::min(target_bp, bridge_bp);
    support.similarity = std::min(1.0, static_cast<double>(support.shared_bp) / compared_bp);
    return support;
}

double GfaGapfillBoundary::minimizer_similarity_(
    uint32_t target,
    uint32_t target_pos,
    uint32_t bridge,
    uint32_t bridge_pos,
    bool before,
    uint64_t windows
) const {
    const GfaGapfill::Fragment& target_fragment = gapfill_.fragments_[target];
    const GfaGapfill::Fragment& bridge_fragment = gapfill_.fragments_[bridge];
    const auto target_region = phase_region_(target, target_pos, before, windows);
    const auto bridge_region = phase_region_(bridge, bridge_pos, before, windows);
    const std::string target_sequence = gapfill_.fragment_sequence_(target_fragment, target_region.first, target_region.second);
    const std::string bridge_sequence = gapfill_.fragment_sequence_(bridge_fragment, bridge_region.first, bridge_region.second);
    if (target_sequence.empty() || bridge_sequence.empty()) return -1.0;

    minimizerdna::Options options;
    options.k = gapfill_.mm2_.k;
    options.w = gapfill_.mm2_.w;
    const minimizerdna::MinimizerBuilder builder(options);
    const minimizerdna::Sketch target_sketch = builder.build(target_sequence);
    const minimizerdna::Sketch bridge_sketch = builder.build(bridge_sequence);
    if (target_sketch.empty() || bridge_sketch.empty()) return -1.0;
    return target_sketch.jaccard(bridge_sketch);
}

// Alignment windows reuse cached graph anchors and only update cuts when minimap2 passes both filters.
bool GfaGapfillBoundary::prepare_alignment_window_(
    uint32_t target_id, uint32_t bridge_id,
    const Boundary& left_boundary,
    const Boundary& right_boundary,
    bool before,
    Boundary& boundary
) const {
    const GfaGapfill::Fragment& target = gapfill_.fragments_[target_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];
    if (target.vertices.empty() || bridge.vertices.empty()) return false;
    uint64_t anchor_begin = 0;
    uint64_t anchor_end = 0;
    if (before) {
        anchor_end = target.length > gapfill_.params_.end_skip_bp ? target.length - gapfill_.params_.end_skip_bp : 0;
        anchor_begin = anchor_end > gapfill_.params_.phase_window_bp ? anchor_end - gapfill_.params_.phase_window_bp : 0;
    } else {
        anchor_begin = std::min(target.length, gapfill_.params_.end_skip_bp);
        anchor_end = anchor_begin + std::min(gapfill_.params_.phase_window_bp, target.length - anchor_begin);
    }
    if (anchor_begin >= anchor_end) return false;

    const uint32_t first = std::min(
        static_cast<uint32_t>(std::upper_bound(target.path_bp.begin(), target.path_bp.end(), anchor_begin) - target.path_bp.begin() - 1),
        static_cast<uint32_t>(target.vertices.size() - 1)
    );
    const uint32_t last = std::min(
        static_cast<uint32_t>(std::lower_bound(target.path_bp.begin(), target.path_bp.end(), anchor_end) - target.path_bp.begin()),
        static_cast<uint32_t>(target.vertices.size() - 1)
    );
    // Anchor the outer edge of the target phase window on the same bridge rank.
    uint32_t target_anchor = UINT32_MAX;
    uint32_t bridge_anchor = UINT32_MAX;

    if (before) {
        for (uint32_t pos = first; pos <= last; ++pos) {
            if (!shared_anchor_(target_id, bridge_id, left_boundary, right_boundary, pos, true, bridge_anchor)) continue;
            target_anchor = pos;
            break;
        }
        for (uint32_t pos = first; target_anchor == UINT32_MAX && pos > 0; --pos) {
            if (shared_anchor_(target_id, bridge_id, left_boundary, right_boundary, pos - 1, true, bridge_anchor)) {
                target_anchor = pos - 1;
            }
        }
    } else {
        for (uint32_t pos = last + 1; target_anchor == UINT32_MAX && pos > first; --pos) {
            if (shared_anchor_(target_id, bridge_id, left_boundary, right_boundary, pos - 1, false, bridge_anchor)) {
                target_anchor = pos - 1;
            }
        }
        for (uint32_t pos = last + 1; target_anchor == UINT32_MAX && pos < target.vertices.size(); ++pos) {
            if (shared_anchor_(target_id, bridge_id, left_boundary, right_boundary, pos, false, bridge_anchor)) {
                target_anchor = pos;
            }
        }
    }
    if (target_anchor == UINT32_MAX) return false;

    // Anchor in trusted sequence, but include the complete contig end in the alignment.
    if (before) {
        boundary.target_begin = target.path_bp[target_anchor];
        boundary.target_end = target.length;
        boundary.bridge_begin = bridge.path_bp[bridge_anchor];
        const uint64_t query_length = boundary.target_end - boundary.target_begin;
        const uint64_t reference_length = query_length + (query_length + 9) / 10;
        boundary.bridge_end = boundary.bridge_begin + std::min(reference_length, bridge.length - boundary.bridge_begin);
    } else {
        boundary.target_begin = 0;
        boundary.target_end = target.path_bp[target_anchor + 1];
        boundary.bridge_end = bridge.path_bp[bridge_anchor + 1];
        const uint64_t query_length = boundary.target_end - boundary.target_begin;
        const uint64_t reference_length = query_length + (query_length + 9) / 10;
        boundary.bridge_begin = boundary.bridge_end > reference_length ? boundary.bridge_end - reference_length : 0;
    }
    boundary.alignable = boundary.target_begin < boundary.target_end && boundary.bridge_begin < boundary.bridge_end;
    return boundary.alignable;
}

bool GfaGapfillBoundary::shared_anchor_(
    uint32_t target_id, uint32_t bridge_id,
    const Boundary& left_boundary,
    const Boundary& right_boundary,
    uint32_t target_pos,
    bool before,
    uint32_t& bridge_pos
) const {
    const GfaGapfill::Fragment& target = gapfill_.fragments_[target_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];
    bridge_pos = find_position_(bridge_id, target.vertices[target_pos]);
    if (bridge_pos == UINT32_MAX) return false;
    const uint32_t rank = before ? left_boundary.bridge_pos : right_boundary.bridge_pos;
    const bool ordered = before ? bridge_pos < right_boundary.bridge_pos : bridge_pos > left_boundary.bridge_pos;
    return ordered && bridge.rank_drops[bridge_pos] == bridge.rank_drops[rank];
}

bool GfaGapfillBoundary::align_(
    const std::string& reference,
    const std::string& query,
    bool before,
    Alignment& alignment
) const {
    if (reference.empty() || query.empty()) return false;

    mm_idxopt_t index_options;
    mm_mapopt_t map_options;
    mm_set_opt(nullptr, &index_options, &map_options);
    if (mm_set_opt(gapfill_.mm2_.preset.c_str(), &index_options, &map_options) < 0) return false;
    index_options.k = static_cast<short>(gapfill_.mm2_.k);
    index_options.w = static_cast<short>(gapfill_.mm2_.w);
    map_options.best_n = static_cast<short>(gapfill_.mm2_.best_n);
    map_options.mid_occ = static_cast<int32_t>(gapfill_.mm2_.max_occ);
    map_options.zdrop = gapfill_.mm2_.zdrop;
    map_options.flag |= MM_F_CIGAR;
    map_options.flag &= ~MM_F_EQX;

    const char* reference_sequence = reference.c_str();
    const char* reference_name = "bridge";
    mm_idx_t* index = mm_idx_str(
        index_options.w, index_options.k,
        (index_options.flag & MM_I_HPC) ? 1 : 0,
        index_options.bucket_bits, 1,
        &reference_sequence, &reference_name
    );
    if (!index) return false;
    mm_mapopt_update(&map_options, index);
    mm_tbuf_t* buffer = mm_tbuf_init();
    if (!buffer) {
        mm_idx_destroy(index);
        return false;
    }

    int count = 0;
    mm_reg1_t* hits = mm_map(
        index, static_cast<int>(query.size()), query.c_str(),
        &count, buffer, &map_options, "boundary"
    );
    std::vector<AlignmentHit> candidates;
    candidates.reserve(count);
    for (int i = 0; i < count; ++i) {
        const mm_reg1_t& hit = hits[i];
        AlignmentHit candidate;
        if (hit.p) {
            candidate.cigar = hit.p->cigar;
            candidate.cigar_size = hit.p->n_cigar;
        }
        candidate.query_begin = hit.qs;
        candidate.query_end = hit.qe;
        candidate.reference_begin = hit.rs;
        candidate.reference_end = hit.re;
        candidate.matches = hit.mlen > 0 ? static_cast<uint64_t>(hit.mlen) : 0;
        candidate.columns = hit.blen > 0 ? static_cast<uint64_t>(hit.blen) : 0;
        candidate.mapq = hit.mapq;
        candidate.reverse = hit.rev;
        candidates.push_back(candidate);
    }
    const double coverage = filter_and_chain_hits_(
        candidates, query.size(), reference.size(), before, alignment
    );
    alignment.accepted = alignment.cut_mapped && coverage >= gapfill_.mm2_.boundary_coverage;
    if (hits) {
        for (int i = 0; i < count; ++i) std::free(hits[i].p);
        std::free(hits);
    }
    mm_tbuf_destroy(buffer);
    mm_idx_destroy(index);
    return alignment.accepted;
}

double GfaGapfillBoundary::filter_and_chain_hits_(
    std::vector<AlignmentHit>& hits,
    uint64_t query_length,
    uint64_t reference_length,
    bool before,
    Alignment& alignment
) const {
    const uint64_t shorter_length = std::min(query_length, reference_length);
    if (shorter_length == 0) return -1.0;

    // Filter each independent hit before rebuilding one forward collinear chain.
    size_t retained = 0;
    for (AlignmentHit& hit : hits) {
        if (!hit.cigar || hit.cigar_size == 0 || hit.reverse || hit.mapq < gapfill_.mm2_.min_mapq || hit.query_begin >= hit.query_end || hit.reference_begin >= hit.reference_end) continue;

        const double match_ratio = hit.columns == 0 ? 0.0 : static_cast<double>(hit.matches) / hit.columns;
        const uint64_t aligned_span = std::min(
            hit.query_end - hit.query_begin,
            hit.reference_end - hit.reference_begin
        );
        const double alignment_ratio = static_cast<double>(aligned_span) / shorter_length;
        if (match_ratio < gapfill_.mm2_.min_match || alignment_ratio < gapfill_.mm2_.min_ali_ratio) continue;
        hits[retained++] = hit;
    }
    hits.resize(retained);
    std::sort(hits.begin(), hits.end());

    // The longest hit is the anchor; add only non-overlapping collinear hits.
    std::map<uint32_t, const AlignmentHit*> chain;
    for (const AlignmentHit& hit : hits) {
        const auto next = chain.lower_bound(hit.query_begin);
        if (next != chain.end() && (hit.query_end > next->second->query_begin || hit.reference_end > next->second->reference_begin)) continue;
        if (next != chain.begin()) {
            const AlignmentHit& previous = *std::prev(next)->second;
            if (previous.query_end > hit.query_begin || previous.reference_end > hit.reference_begin) continue;
        }
        chain.emplace(hit.query_begin, &hit);
    }
    if (chain.empty()) return -1.0;

    uint64_t matches = 0;
    uint64_t columns = 0;
    uint64_t covered = 0;
    alignment.mapq = UINT8_MAX;
    for (const auto& item : chain) {
        const AlignmentHit& hit = *item.second;
        matches += hit.matches;
        columns += hit.columns;
        covered += hit.query_end - hit.query_begin;
        alignment.mapq = std::min(alignment.mapq, hit.mapq);

        uint64_t query = hit.query_begin;
        uint64_t reference = hit.reference_begin;
        for (uint32_t i = 0; i < hit.cigar_size; ++i) {
            const uint32_t length = hit.cigar[i] >> 4;
            const uint32_t operation = hit.cigar[i] & 0xf;
            if (operation == MM_CIGAR_MATCH) {
                alignment.matches.push_back({
                    static_cast<uint32_t>(query), static_cast<uint32_t>(reference), length
                });
                if (before || !alignment.cut_mapped) {
                    alignment.query_cut = before ? query + length : query;
                    alignment.reference_cut = before ? reference + length : reference;
                    alignment.cut_mapped = true;
                }
            }
            if (operation == MM_CIGAR_MATCH || operation == MM_CIGAR_INS) query += length;
            if (operation == MM_CIGAR_MATCH || operation == MM_CIGAR_DEL || operation == MM_CIGAR_N_SKIP) reference += length;
        }
    }
    const AlignmentHit& first = *chain.begin()->second;
    const AlignmentHit& last = *chain.rbegin()->second;
    alignment.query_begin = first.query_begin;
    alignment.query_end = last.query_end;
    alignment.reference_begin = first.reference_begin;
    alignment.reference_end = last.reference_end;
    alignment.hits = static_cast<uint32_t>(chain.size());
    alignment.identity = columns == 0 ? 0.0 : static_cast<double>(matches) / columns;
    alignment.coverage = std::min(1.0, static_cast<double>(covered) / query_length);
    return alignment.coverage;
}

bool GfaGapfillBoundary::splice_overlap_(
    const Alignment& left,
    const Alignment& right,
    const Boundary& left_boundary,
    const Boundary& right_boundary,
    uint64_t& left_cut,
    uint64_t& right_cut,
    uint64_t& bridge_left_cut,
    uint64_t& bridge_right_cut
) const {
    const uint64_t overlap_begin = right_boundary.bridge_begin + right.reference_cut;
    const uint64_t overlap_end = left_boundary.bridge_begin + left.reference_cut;
    if (overlap_begin >= overlap_end || left.matches.empty() || right.matches.empty()) return false;

    // Find the longest reference interval covered by M blocks from both targets.
    size_t i = 0, j = 0;
    uint64_t best_span = 0;
    const Alignment::MatchBlock* best_left = nullptr;
    const Alignment::MatchBlock* best_right = nullptr;
    uint64_t best_begin = 0, best_end = 0;
    while (i < left.matches.size() && j < right.matches.size()) {
        const Alignment::MatchBlock& a = left.matches[i];
        const Alignment::MatchBlock& b = right.matches[j];
        const uint64_t a_begin = left_boundary.bridge_begin + a.reference_begin;
        const uint64_t b_begin = right_boundary.bridge_begin + b.reference_begin;
        const uint64_t a_end = a_begin + a.length;
        const uint64_t b_end = b_begin + b.length;
        const uint64_t begin = std::max({a_begin, b_begin, overlap_begin});
        const uint64_t end = std::min({a_end, b_end, overlap_end});
        if (begin < end && (!best_left || end - begin > best_span)) {
            best_span = end - begin;
            best_left = &a;
            best_right = &b;
            best_begin = begin;
            best_end = end;
        }
        if (a_end <= b_end) ++i;
        if (b_end <= a_end) ++j;
    }
    if (best_left) {
        bridge_left_cut = best_begin + (best_end - best_begin) / 2;
        bridge_right_cut = bridge_left_cut;
        left_cut = left_boundary.target_begin + best_left->query_begin + bridge_left_cut - (left_boundary.bridge_begin + best_left->reference_begin);
        right_cut = right_boundary.target_begin + best_right->query_begin + bridge_right_cut - (right_boundary.bridge_begin + best_right->reference_begin);
        return left_cut <= left_boundary.target_end && right_cut <= right_boundary.target_end;
    }

    // No common M: use the nearest ordered M ends and fill their interval from the bridge.
    size_t left_index = 0;
    const Alignment::MatchBlock* nearest_left = nullptr;
    const Alignment::MatchBlock* nearest_right = nullptr;
    const Alignment::MatchBlock* latest_left = nullptr;
    uint64_t nearest_distance = UINT64_MAX;
    uint64_t nearest_loss = UINT64_MAX;
    uint64_t nearest_left_end = 0, nearest_right_begin = 0, latest_left_end = 0;
    for (const Alignment::MatchBlock& b : right.matches) {
        const uint64_t b_begin = right_boundary.bridge_begin + b.reference_begin;
        while (left_index < left.matches.size()) {
            const Alignment::MatchBlock& a = left.matches[left_index];
            const uint64_t a_end = left_boundary.bridge_begin + a.reference_begin + a.length;
            if (a_end > b_begin) break;
            if (!latest_left || a_end > latest_left_end) {
                latest_left = &a;
                latest_left_end = a_end;
            }
            ++left_index;
        }
        if (!latest_left) continue;
        const uint64_t distance = b_begin - latest_left_end;
        const uint64_t left_query_end = latest_left->query_begin + latest_left->length;
        const uint64_t loss = (left.query_cut > left_query_end ? left.query_cut - left_query_end : 0) + (b.query_begin > right.query_cut ? b.query_begin - right.query_cut : 0);
        if (distance < nearest_distance || (distance == nearest_distance && loss < nearest_loss)) {
            nearest_distance = distance;
            nearest_loss = loss;
            nearest_left = latest_left;
            nearest_right = &b;
            nearest_left_end = latest_left_end;
            nearest_right_begin = b_begin;
        }
    }
    if (!nearest_left) return false;
    if (nearest_right_begin - nearest_left_end > gapfill_.params_.max_gap_bp) return false;

    bridge_left_cut = nearest_left_end;
    bridge_right_cut = nearest_right_begin;
    left_cut = left_boundary.target_begin + nearest_left->query_begin + nearest_left->length;
    right_cut = right_boundary.target_begin + nearest_right->query_begin;
    return left_cut <= left_boundary.target_end && right_cut <= right_boundary.target_end;
}

uint32_t GfaGapfillBoundary::find_position_(
    uint32_t fragment,
    uint32_t vertex
) const {
    const auto& index = gapfill_.fragments_[fragment].vertex_index;
    const auto it = std::lower_bound(
        index.begin(), index.end(), std::pair<uint32_t, uint32_t>{vertex, 0}
    );
    return it != index.end() && it->first == vertex ? it->second : UINT32_MAX;
}

/* =================================================================================================================
 * GAP BOUNDARY END
 * ================================================================================================================= */
