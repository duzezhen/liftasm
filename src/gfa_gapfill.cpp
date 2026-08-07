#include "gfa_gapfill.hpp"
#include "gfa_bubble.hpp"

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

using PhaseIntervals = std::map<uint32_t, uint32_t>;

struct LocalPhaseAssignments {
    std::array<PhaseIntervals, 2> hap;
};

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

std::string safe_name(std::string name) {
    for (char& c : name) {
        const bool valid = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '.' || c == '_' || c == '-';
        if (!valid) c = '_';
    }
    return name.empty() ? "sample" : name;
}

std::string numbered_contig_name(uint8_t hap, const char* kind, size_t number) {
    std::ostringstream name;
    name << 'h' << static_cast<uint32_t>(hap) << kind
         << std::setw(6) << std::setfill('0') << number << 'l';
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

double phase_odds(double similarity) {
    constexpr double epsilon = 1e-6;
    similarity = std::clamp(similarity, epsilon, 1.0 - epsilon);
    return similarity / (1.0 - similarity);
}

std::string json_string(const std::string& value) {
    std::string out;
    out.reserve(value.size() + 2);
    out.push_back('"');
    for (unsigned char c : value) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    const char hex[] = "0123456789abcdef";
                    out += "\\u00";
                    out.push_back(hex[c >> 4]);
                    out.push_back(hex[c & 15]);
                } else {
                    out.push_back(static_cast<char>(c));
                }
        }
    }
    out.push_back('"');
    return out;
}

bool closest_match_cut(
    const mm_reg1_t& hit,
    bool left_boundary,
    uint64_t& query_cut,
    uint64_t& reference_cut
) {
    if (!hit.p || hit.p->n_cigar == 0) return false;

    uint64_t query = hit.qs;
    uint64_t reference = hit.rs;
    bool found = false;
    for (uint32_t i = 0; i < hit.p->n_cigar; ++i) {
        const uint32_t length = hit.p->cigar[i] >> 4;
        const uint32_t operation = hit.p->cigar[i] & 0xf;
        if (operation == MM_CIGAR_MATCH) {
            if (!left_boundary) {
                query_cut = query;
                reference_cut = reference;
                return true;
            }
            query_cut = query + length;
            reference_cut = reference + length;
            found = true;
        }
        if (operation == MM_CIGAR_MATCH || operation == MM_CIGAR_INS) query += length;
        if (operation == MM_CIGAR_MATCH || operation == MM_CIGAR_DEL ||
            operation == MM_CIGAR_N_SKIP) reference += length;
    }
    return found;
}

std::string metric(double value) {
    if (value < 0.0) return "NA";
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

std::string boundary_table_header() {
    std::ostringstream out;
    out << "  - " << std::left << std::setw(9) << "Boundary"
        << std::right << std::setw(8) << "Phase"
        << std::setw(18) << "Bubble bp"
        << std::setw(12) << "Minimizer"
        << std::setw(12) << "Alignment";
    return out.str();
}

std::string boundary_table_row(
    const char* label,
    const GfaGapfillDebugger::BoundarySupport& support,
    double alignment
) {
    std::ostringstream bubble;
    bubble << support.first_bubble_bp << '/' << support.second_bubble_bp;
    std::ostringstream out;
    out << "  - " << std::left << std::setw(9) << label
        << std::right << std::setw(8) << metric(support.phase)
        << std::setw(18) << bubble.str()
        << std::setw(12) << metric(support.minimizer)
        << std::setw(12) << metric(alignment);
    return out.str();
}

void place_contigs(std::vector<GfaGapfillPlotter::Contig>& contigs);

}  // namespace

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
    double boundary_identity,
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
    mm2_.boundary_identity = boundary_identity;
    mm2_.boundary_coverage = boundary_coverage;
    params_.threads = std::max<uint32_t>(1, align.threads);
}

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
    marker = t1 != std::string::npos ? t1 : t2;
    if (marker == std::string::npos || marker == 0) return false;
    sample = name.substr(0, marker);
    hap = t1 != std::string::npos ? 1 : 2;
    return true;
}

uint32_t GfaGapfill::find_position_(const Fragment& fragment, uint32_t vertex) {
    const auto it = std::lower_bound(
        fragment.vertex_index.begin(), fragment.vertex_index.end(),
        std::pair<uint32_t, uint32_t>{vertex, 0}
    );
    return it != fragment.vertex_index.end() && it->first == vertex ? it->second : UINT32_MAX;
}

std::pair<uint32_t, uint64_t> GfaGapfill::path_region_(
    const Fragment& fragment,
    uint32_t beg,
    uint32_t end
) {
    if (beg > end || end >= fragment.vertices.size()) return {0, 0};
    return {end - beg + 1, fragment.path_bp[end + 1] - fragment.path_bp[beg]};
}

bool GfaGapfill::candidate_better_(const Candidate& a, const Candidate& b) {
    if (a.phase_consistent != b.phase_consistent) return a.phase_consistent;
    if (a.phase_score != b.phase_score) return a.phase_score > b.phase_score;
    if (a.probability != b.probability) return a.probability > b.probability;
    if (a.score != b.score) return a.score > b.score;
    const uint64_t as = a.left_shared + a.right_shared;
    const uint64_t bs = b.left_shared + b.right_shared;
    if (as != bs) return as > bs;
    const uint32_t alen = a.bridge_right - a.bridge_left;
    const uint32_t blen = b.bridge_right - b.bridge_left;
    if (alen != blen) return alen < blen;
    return a.bridge < b.bridge;
}

bool GfaGapfill::connection_better_(const Candidate& a, const Candidate& b) {
    if (a.target_distance != b.target_distance) return a.target_distance < b.target_distance;
    if (a.target_overlaps != b.target_overlaps) return !a.target_overlaps;
    return candidate_better_(a, b);
}

std::string GfaGapfill::vertex_name_(uint32_t vertex) const {
    const Vertex v(vertex);
    return nodes_[v.segment_id()].name + (v.is_reverse() ? "-" : "+");
}

bool GfaGapfill::build_fragment_(size_t path_id, Fragment& fragment) {
    GfaPath& path = paths_[path_id];
    if (!parse_path_name_(path.name, fragment.sample, fragment.hap)) {
        GfaGapfillDebugger::path(path.name, "drop:unrecognized-sample");
        return false;
    }
    path.name = compact_contig_name(path.name);
    if (path.segments.empty()) {
        GfaGapfillDebugger::path(path.name, "drop:empty");
        return false;
    }

    fragment.path_id = path_id;
    std::unordered_set<uint32_t> seen;
    fragment.vertices.reserve(path.segments.size());

    for (const PathSegment& segment : path.segments) {
        const uint32_t sid = static_cast<uint32_t>(segment.node_id);
        if (sid >= nodes_.size() || nodes_[sid].deleted) {
            GfaGapfillDebugger::path(path.name, "drop:invalid-segment");
            return false;
        }
        if (!seen.insert(sid).second) {
            GfaGapfillDebugger::path(path.name, "drop:repeated-segment node=" + nodes_[sid].name);
            return false;
        }

        const uint32_t component = connectivity_index_[sid];
        if (fragment.component == UINT32_MAX) fragment.component = component;
        else if (fragment.component != component) {
            GfaGapfillDebugger::path(path.name, "drop:multiple-components");
            return false;
        }

        fragment.vertices.push_back(Vertex::make_vertex(sid, segment.is_reverse));
        fragment.length += nodes_[sid].length;
    }

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

    rebuild_fragment_index_(fragment);
    fragment.eligible = fragment.length >= params_.min_contig_bp;
    GfaGapfillDebugger::path(
        path.name, fragment.eligible ? "keep" : "exclude:short", fragment.sample, fragment.hap, fragment.component,
        fragment.vertices.size(), fragment.length,
        vertex_name_(fragment.vertices.front()), vertex_name_(fragment.vertices.back())
    );
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

void GfaGapfill::build_fragments_() {
    log_stream() << "Indexing sample contig paths ...\n";
    fragments_.reserve(paths_.size());
    for (size_t path_id = 0; path_id < paths_.size(); ++path_id) {
        Fragment fragment;
        if (build_fragment_(path_id, fragment)) fragments_.push_back(std::move(fragment));
    }
    log_stream() << "  - Contig paths indexed: " << fragments_.size() << "\n\n";
}

void GfaGapfill::build_graph_order_() {
    log_stream() << "Ordering sample contigs from shared graph nodes ...\n";
    std::vector<GfaGapfillPlotter::Contig> layout;
    layout.reserve(fragments_.size());
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible) continue;
        GfaGapfillPlotter::Contig contig;
        contig.id = i;
        contig.component = fragment.component;
        contig.hap = fragment.hap;
        contig.bp = fragment.length;
        contig.name = paths_[fragment.path_id].name;
        contig.sample = fragment.sample;
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
    place_contigs(layout);

    for (const GfaGapfillPlotter::Contig& contig : layout) {
        Fragment& fragment = fragments_[contig.id];
        fragment.layout_start = contig.start;
        if (!contig.reverse) continue;
        std::reverse(fragment.vertices.begin(), fragment.vertices.end());
        for (uint32_t& vertex : fragment.vertices) vertex ^= 1u;
        fragment.reverse = !fragment.reverse;
        rebuild_fragment_index_(fragment);
    }

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

bool GfaGapfill::ordered_pair_(const Fragment& left, const Fragment& right) const {
    if (left.sample != right.sample || left.component != right.component ||
        left.length < params_.min_contig_bp || right.length < params_.min_contig_bp) return false;
    if (left.layout_start > right.layout_start) return false;
    const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
    const int64_t right_end = right.layout_start + static_cast<int64_t>(right.length);
    if (right.layout_start >= left_end) {
        if (left.relocation >= 0 || right.relocation >= 0) return true;
        return static_cast<uint64_t>(right.layout_start - left_end) <= params_.max_gap_bp;
    }
    const uint64_t overlap = static_cast<uint64_t>(
        std::min(left_end, right_end) - right.layout_start
    );
    return static_cast<double>(overlap) / std::min(left.length, right.length) <= params_.max_target_overlap;
}

GfaGapfill::OverlapSupport GfaGapfill::overlap_support_(
    const Fragment& a,
    const Fragment& b,
    const PathOverlap& overlap
) const {
    OverlapSupport support;
    if (overlap.nodes == 0) return support;

    const auto a_region = path_region_(a, overlap.a_begin, overlap.a_end);
    const auto b_region = path_region_(b, overlap.b_begin, overlap.b_end);
    const uint32_t total_nodes = a_region.first + b_region.first;
    if (total_nodes < overlap.nodes) return support;

    support.shared_bp = overlap.bp;
    support.overlap_bp = std::min(a_region.second, b_region.second);
    support.similarity = support.overlap_bp == 0 ? 0.0 : std::min(
        1.0, static_cast<double>(support.shared_bp) / support.overlap_bp
    );
    const uint64_t a_bubble = a.bubble_bp[overlap.a_end + 1] - a.bubble_bp[overlap.a_begin];
    const uint64_t b_bubble = b.bubble_bp[overlap.b_end + 1] - b.bubble_bp[overlap.b_begin];
    support.bubble_shared_bp = std::min({overlap.bubble_bp, a_bubble, b_bubble});
    support.bubble_overlap_bp = std::min(a_bubble, b_bubble);
    support.phase_similarity = support.bubble_overlap_bp == 0 ? 1.0 :
        static_cast<double>(support.bubble_shared_bp) / support.bubble_overlap_bp;
    return support;
}

bool GfaGapfill::gap_boundary_supported_(
    const Fragment& a,
    const Fragment& b,
    const PathOverlap& overlap,
    bool right_boundary
) const {
    const uint32_t a_pos = right_boundary ? overlap.a_end : overlap.a_begin;
    const uint32_t b_pos = right_boundary ? overlap.b_end : overlap.b_begin;
    return boundary_similarity_(
        a, a_pos, b, b_pos, right_boundary, 1
    ) >= params_.min_similarity;
}

std::string GfaGapfill::boundary_sequence_(
    const Fragment& fragment,
    uint32_t pos,
    bool before,
    uint64_t windows
) const {
    const auto region = boundary_phase_region_(fragment, pos, before, windows);
    return fragment_sequence_(fragment, region.first, region.second);
}

double GfaGapfill::boundary_similarity_(
    const Fragment& target,
    uint32_t target_pos,
    const Fragment& bridge,
    uint32_t bridge_pos,
    bool before,
    uint64_t windows
) const {
    const std::string target_sequence = boundary_sequence_(target, target_pos, before, windows);
    const std::string bridge_sequence = boundary_sequence_(bridge, bridge_pos, before, windows);
    if (target_sequence.empty() || bridge_sequence.empty()) return -1.0;

    minimizerdna::Options options;
    options.k = mm2_.k;
    options.w = mm2_.w;
    const minimizerdna::MinimizerBuilder builder(options);
    const minimizerdna::Sketch target_sketch = builder.build(target_sequence);
    const minimizerdna::Sketch bridge_sketch = builder.build(bridge_sequence);
    if (target_sketch.empty() || bridge_sketch.empty()) return -1.0;
    return target_sketch.jaccard(bridge_sketch);
}

std::pair<uint32_t, uint32_t> GfaGapfill::boundary_phase_region_(
    const Fragment& fragment,
    uint32_t pos,
    bool before,
    uint64_t windows
) const {
    const uint64_t window = params_.phase_window_bp;
    const uint64_t available = boundary_phase_available_(fragment, pos, before);
    const uint64_t max_windows = available / window + (available % window != 0);
    const uint64_t span = windows >= max_windows ? available : window * windows;
    uint64_t begin_bp = 0;
    uint64_t end_bp = 0;
    if (before) {
        end_bp = std::min(
            fragment.path_bp[pos + 1],
            fragment.length > params_.phase_skip_bp ?
                fragment.length - params_.phase_skip_bp : 0
        );
        begin_bp = end_bp - span;
    } else {
        begin_bp = std::max(fragment.path_bp[pos], params_.phase_skip_bp);
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

std::pair<uint64_t, uint64_t> GfaGapfill::boundary_alignment_interval_(
    const Fragment& fragment,
    bool before
) const {
    uint64_t begin_bp = 0;
    uint64_t end_bp = 0;
    if (before) {
        end_bp = fragment.length > params_.phase_skip_bp ?
            fragment.length - params_.phase_skip_bp : 0;
        begin_bp = end_bp > params_.phase_window_bp ?
            end_bp - params_.phase_window_bp : 0;
    } else {
        begin_bp = std::min(fragment.length, params_.phase_skip_bp);
        end_bp = begin_bp + std::min(
            params_.phase_window_bp, fragment.length - begin_bp
        );
    }
    return end_bp > begin_bp ?
        std::make_pair(begin_bp, end_bp) : std::pair<uint64_t, uint64_t>{0, 0};
}

uint64_t GfaGapfill::boundary_phase_available_(
    const Fragment& fragment,
    uint32_t pos,
    bool before
) const {
    if (before) {
        return std::min(
            fragment.path_bp[pos + 1],
            fragment.length > params_.phase_skip_bp ?
                fragment.length - params_.phase_skip_bp : 0
        );
    }
    const uint64_t begin = std::min(
        std::max(fragment.path_bp[pos], params_.phase_skip_bp),
        fragment.length
    );
    return fragment.length - begin;
}

GfaGapfill::BoundaryPhaseSupport GfaGapfill::boundary_phase_similarity_(
    const Fragment& target,
    uint32_t target_pos,
    const Fragment& bridge,
    uint32_t bridge_pos,
    bool before
) const {
    BoundaryPhaseSupport support;
    const uint64_t window = params_.phase_window_bp;
    const uint64_t target_available = boundary_phase_available_(target, target_pos, before);
    const uint64_t bridge_available = boundary_phase_available_(bridge, bridge_pos, before);
    const uint64_t target_windows = target_available / window + (target_available % window != 0);
    const uint64_t bridge_windows = bridge_available / window + (bridge_available % window != 0);
    uint64_t windows = 1;
    support.windows = windows;
    auto target_region = boundary_phase_region_(target, target_pos, before, windows);
    auto bridge_region = boundary_phase_region_(bridge, bridge_pos, before, windows);
    uint64_t target_bp = target.bubble_bp[target_region.second] - target.bubble_bp[target_region.first];
    uint64_t bridge_bp = bridge.bubble_bp[bridge_region.second] - bridge.bubble_bp[bridge_region.first];
    while (target_bp < params_.phase_min_bp || bridge_bp < params_.phase_min_bp) {
        if ((target_bp < params_.phase_min_bp && windows >= target_windows) || (bridge_bp < params_.phase_min_bp && windows >= bridge_windows)) {
            support.target_bp = target_bp;
            support.bridge_bp = bridge_bp;
            support.windows = windows;
            return support;
        }
        ++windows;
        support.windows = windows;
        target_region = boundary_phase_region_(target, target_pos, before, windows);
        bridge_region = boundary_phase_region_(bridge, bridge_pos, before, windows);
        target_bp = target.bubble_bp[target_region.second] - target.bubble_bp[target_region.first];
        bridge_bp = bridge.bubble_bp[bridge_region.second] - bridge.bubble_bp[bridge_region.first];
    }
    support.target_bp = target_bp;
    support.bridge_bp = bridge_bp;

    const Fragment* query = &target;
    const Fragment* reference = &bridge;
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
        const uint32_t reference_pos = find_position_(*reference, query->vertices[pos]);
        if (reference_pos >= reference_region.first && reference_pos < reference_region.second) {
            support.shared_bp += nodes_[Vertex::get_segment_id(query->vertices[pos])].length;
        }
    }
    const uint64_t compared_bp = std::min(target_bp, bridge_bp);
    support.similarity = std::min(
        1.0, static_cast<double>(support.shared_bp) / compared_bp
    );
    return support;
}

bool GfaGapfill::find_anchors_(
    const Fragment& left,
    const Fragment& right,
    const Fragment& bridge,
    Candidate& candidate
) const {
    uint64_t left_sum = 0, right_sum = 0, direct_bp = 0;
    uint32_t direct_matches = 0;
    for (uint32_t left_pos = 0; left_pos < left.vertices.size(); ++left_pos) {
        const uint32_t right_pos = find_position_(right, left.vertices[left_pos]);
        if (right_pos == UINT32_MAX) continue;
        left_sum += left_pos;
        right_sum += right_pos;
        direct_bp += nodes_[Vertex::get_segment_id(left.vertices[left_pos])].length;
        ++direct_matches;
    }
    if (direct_matches > 0) {
        const uint64_t shorter_bp = std::min(left.length, right.length);
        if (shorter_bp == 0 || static_cast<double>(direct_bp) / shorter_bp > params_.max_target_overlap) {
            return false;
        }
        if (2 * left_sum + direct_matches <= direct_matches * left.vertices.size() || 2 * right_sum + direct_matches >= direct_matches * right.vertices.size()) {
            return false;
        }
    }

    std::vector<std::pair<uint32_t, uint32_t>> left_matches;
    std::vector<std::pair<uint32_t, uint32_t>> right_matches;
    left_matches.reserve(std::min(left.vertices.size(), bridge.vertices.size()));
    right_matches.reserve(std::min(right.vertices.size(), bridge.vertices.size()));

    for (uint32_t bridge_pos = 0; bridge_pos < bridge.vertices.size(); ++bridge_pos) {
        const uint32_t vertex = bridge.vertices[bridge_pos];
        const uint32_t left_pos = find_position_(left, vertex);
        if (left_pos != UINT32_MAX) left_matches.emplace_back(bridge_pos, left_pos);
        const uint32_t right_pos = find_position_(right, vertex);
        if (right_pos != UINT32_MAX) right_matches.emplace_back(bridge_pos, right_pos);
    }
    if (left_matches.empty() || right_matches.empty()) return false;

    std::vector<uint32_t> suffix(right_matches.size());
    suffix.back() = static_cast<uint32_t>(right_matches.size() - 1);
    for (size_t i = right_matches.size() - 1; i > 0; --i) {
        const uint32_t current = static_cast<uint32_t>(i - 1);
        const uint32_t best = suffix[i];
        suffix[current] = right_matches[current].second < right_matches[best].second ? current : best;
    }

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
        if (bridge.rank_drops[right_match.first] != bridge.rank_drops[left_match.first]) continue;
        const uint64_t retained = left.path_bp[left_match.second + 1] + right.length - right.path_bp[right_match.second];
        const uint32_t span = right_match.first - left_match.first;
        if (!found || retained > best_retained || (retained == best_retained && span < best_span)) {
            candidate.left_pos = left_match.second;
            candidate.bridge_left = left_match.first;
            candidate.bridge_right = right_match.first;
            candidate.right_pos = right_match.second;
            best_retained = retained;
            best_span = span;
            found = true;
        }
    }
    if (!found) return false;

    Candidate normal = candidate;
    normal.junction_left_pos = normal.left_pos;
    normal.junction_bridge_left = normal.bridge_left;
    normal.junction_bridge_right = normal.bridge_right;
    normal.junction_right_pos = normal.right_pos;

    Candidate crossover = normal;
    uint64_t crossover_retained = 0;
    bool have_crossover = false;
    for (const auto& match : left_matches) {
        const uint32_t vertex = bridge.vertices[match.first];
        const uint32_t right_pos = find_position_(right, vertex);
        if (right_pos == UINT32_MAX) continue;
        const uint64_t retained = left.path_bp[match.second] + right.length - right.path_bp[right_pos];
        const bool prefer_right_intersection = right.length > left.length && match.first < crossover.bridge_left;
        const bool prefer_left_intersection = left.length > right.length && match.first > crossover.bridge_left;
        if (!have_crossover || retained > crossover_retained || (retained == crossover_retained && (prefer_right_intersection || prefer_left_intersection))) {
            crossover.left_pos = match.second;
            crossover.bridge_left = match.first;
            crossover.bridge_right = match.first;
            crossover.right_pos = right_pos;
            crossover.junction_left_pos = match.second;
            crossover.junction_bridge_left = match.first;
            crossover.junction_bridge_right = match.first;
            crossover.junction_right_pos = right_pos;
            crossover_retained = retained;
            have_crossover = true;
        }
    }

    if (have_crossover && expand_boundary_(
            left, right, bridge, left_matches, right_matches, crossover)) {
        candidate = crossover;
    } else {
        candidate = normal;
        if (!expand_boundary_(left, right, bridge, left_matches, right_matches, candidate)) return false;
    }

    return candidate_path_unique_(left, right, bridge, candidate);
}

bool GfaGapfill::candidate_path_unique_(
    const Fragment& left,
    const Fragment& right,
    const Fragment& bridge,
    const Candidate& candidate
) const {
    std::unordered_set<uint32_t> seen;
    for (uint32_t i = 0; i <= candidate.left_pos; ++i) {
        if (!seen.insert(Vertex::get_segment_id(left.vertices[i])).second) return false;
    }
    for (uint32_t i = candidate.bridge_left + 1; i < candidate.bridge_right; ++i) {
        if (!seen.insert(Vertex::get_segment_id(bridge.vertices[i])).second) return false;
    }
    for (uint32_t i = candidate.right_pos; i < right.vertices.size(); ++i) {
        if (!seen.insert(Vertex::get_segment_id(right.vertices[i])).second) return false;
    }
    return true;
}

bool GfaGapfill::homolog_spans_(
    uint32_t left_id,
    uint32_t right_id,
    const OverlapIndex& overlaps
) const {
    const Fragment& left = fragments_[left_id];
    const Fragment& right = fragments_[right_id];
    for (const auto& item : overlaps[left_id]) {
        const uint32_t bridge_id = item.first;
        const Fragment& bridge = fragments_[bridge_id];
        if (bridge.sample != left.sample || bridge.hap == left.hap || bridge.component != left.component) continue;

        const auto right_it = overlaps[bridge_id].find(right_id);
        if (right_it == overlaps[bridge_id].end()) continue;
        const OverlapSupport left_support = overlap_support_(left, bridge, item.second);
        const OverlapSupport right_support = overlap_support_(bridge, right, right_it->second);
        if (
            left_support.overlap_bp < params_.min_overlap_bp || right_support.overlap_bp < params_.min_overlap_bp ||
            left_support.similarity < params_.min_similarity || right_support.similarity < params_.min_similarity
        ) continue;

        const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
        const int64_t bridge_end = bridge.layout_start + static_cast<int64_t>(bridge.length);
        if (bridge.layout_start > left_end || bridge_end < right.layout_start) continue;

        Candidate witness;
        if (find_anchors_(left, right, bridge, witness)) return true;
    }
    return false;
}

bool GfaGapfill::expand_boundary_(
    const Fragment& left,
    const Fragment& right,
    const Fragment& bridge,
    const std::vector<std::pair<uint32_t, uint32_t>>& left_matches,
    const std::vector<std::pair<uint32_t, uint32_t>>& right_matches,
    Candidate& candidate
) const {
    const uint64_t left_junction = left.path_bp[candidate.junction_left_pos + 1];
    const uint64_t left_trusted_end = left.length > params_.phase_skip_bp ?
        left.length - params_.phase_skip_bp : 0;
    const uint64_t left_goal = std::min(left_junction, left_trusted_end);
    uint64_t best_left_distance = UINT64_MAX;
    uint32_t best_left_bridge = candidate.junction_bridge_left;
    uint32_t best_left_pos = candidate.junction_left_pos;
    bool have_left = false;
    for (const auto& match : left_matches) {
        if (match.first > candidate.junction_bridge_left || match.second > candidate.junction_left_pos) continue;
        if (bridge.rank_drops[match.first] != bridge.rank_drops[candidate.junction_bridge_left]) continue;
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

    const uint64_t right_junction = right.path_bp[candidate.junction_right_pos];
    const uint64_t right_goal = std::max(right_junction, params_.phase_skip_bp);
    uint64_t best_right_distance = UINT64_MAX;
    uint32_t best_right_bridge = candidate.junction_bridge_right;
    uint32_t best_right_pos = candidate.junction_right_pos;
    bool have_right = false;
    for (const auto& match : right_matches) {
        if (match.first < candidate.junction_bridge_right || match.second < candidate.junction_right_pos) continue;
        if (bridge.rank_drops[match.first] != bridge.rank_drops[candidate.junction_bridge_right]) continue;
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
    candidate.left_pos = best_left_pos;
    candidate.bridge_left = best_left_bridge;
    candidate.bridge_right = best_right_bridge;
    candidate.right_pos = best_right_pos;
    return candidate.bridge_right > candidate.bridge_left;
}

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
    uint32_t pos = first == fragment.path_bp.begin() ? 0 :
        static_cast<uint32_t>(first - fragment.path_bp.begin() - 1);
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

bool GfaGapfill::align_boundary_(
    const std::string& reference,
    const std::string& query,
    bool before,
    BoundaryAlignment& alignment
) const {
    if (reference.empty() || query.empty()) return false;

    mm_idxopt_t index_options;
    mm_mapopt_t map_options;
    mm_set_opt(nullptr, &index_options, &map_options);
    if (mm_set_opt(mm2_.preset.c_str(), &index_options, &map_options) < 0) return false;
    index_options.k = static_cast<short>(mm2_.k);
    index_options.w = static_cast<short>(mm2_.w);
    map_options.best_n = 1;
    map_options.mid_occ = static_cast<int32_t>(mm2_.max_occ);
    map_options.zdrop = mm2_.zdrop;
    map_options.flag |= MM_F_CIGAR;

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
    alignment.hits = static_cast<uint32_t>(count);
    const mm_reg1_t* best = nullptr;
    const mm_reg1_t* raw_best = nullptr;
    uint64_t best_distance = UINT64_MAX;
    for (int i = 0; i < count; ++i) {
        const mm_reg1_t& hit = hits[i];
        if (hit.blen <= 0 || hit.parent != hit.id) continue;
        if (!raw_best || hit.mlen > raw_best->mlen) raw_best = &hit;
        if (hit.rev) continue;
        const double identity = static_cast<double>(hit.mlen) / hit.blen;
        const double coverage = static_cast<double>(hit.qe - hit.qs) / query.size();
        if (identity < mm2_.boundary_identity || coverage < mm2_.boundary_coverage) continue;

        const uint64_t distance = before ?
            static_cast<uint64_t>(hit.qs) + hit.rs :
            query.size() - hit.qe + reference.size() - hit.re;
        if (!best || distance < best_distance ||
            (distance == best_distance && hit.mlen > best->mlen)) {
            best = &hit;
            best_distance = distance;
        }
    }

    const mm_reg1_t* reported = best ? best : raw_best;
    if (reported) {
        alignment.query_begin = reported->qs;
        alignment.query_end = reported->qe;
        alignment.reference_begin = reported->rs;
        alignment.reference_end = reported->re;
        alignment.mapq = reported->mapq;
        alignment.identity = static_cast<double>(reported->mlen) / reported->blen;
        alignment.coverage = static_cast<double>(reported->qe - reported->qs) / query.size();
        alignment.reverse = reported->rev;
        if (!reported->rev) {
            closest_match_cut(
                *reported, before, alignment.query_cut, alignment.reference_cut
            );
        }
    }
    if (best) {
        alignment.cut_mapped = closest_match_cut(
            *best, before, alignment.query_cut, alignment.reference_cut
        );
        alignment.accepted = alignment.cut_mapped;
    }
    if (hits) {
        for (int i = 0; i < count; ++i) std::free(hits[i].p);
        std::free(hits);
    }
    mm_tbuf_destroy(buffer);
    mm_idx_destroy(index);
    return alignment.accepted;
}

bool GfaGapfill::shared_boundary_anchor_(
    const Fragment& target,
    const Fragment& bridge,
    const Candidate& candidate,
    uint32_t target_pos,
    bool before,
    uint32_t& bridge_pos
) const {
    bridge_pos = find_position_(bridge, target.vertices[target_pos]);
    if (bridge_pos == UINT32_MAX) return false;
    const uint32_t rank = before ? candidate.bridge_left : candidate.bridge_right;
    const bool ordered = before ?
        bridge_pos < candidate.bridge_right : bridge_pos > candidate.bridge_left;
    return ordered && bridge.rank_drops[bridge_pos] == bridge.rank_drops[rank];
}

bool GfaGapfill::boundary_alignment_window_(
    const Fragment& target,
    const Fragment& bridge,
    const Candidate& candidate,
    bool before,
    BoundaryWindow& window
) const {
    if (target.vertices.empty() || bridge.vertices.empty()) return false;
    const auto interval = boundary_alignment_interval_(target, before);
    if (interval.first >= interval.second) return false;

    const uint32_t first = std::min(
        static_cast<uint32_t>(std::upper_bound(
            target.path_bp.begin(), target.path_bp.end(), interval.first
        ) - target.path_bp.begin() - 1),
        static_cast<uint32_t>(target.vertices.size() - 1)
    );
    const uint32_t last = std::min(
        static_cast<uint32_t>(std::lower_bound(
            target.path_bp.begin(), target.path_bp.end(), interval.second
        ) - target.path_bp.begin()),
        static_cast<uint32_t>(target.vertices.size() - 1)
    );
    uint32_t target_anchor = UINT32_MAX;
    uint32_t bridge_anchor = UINT32_MAX;

    if (before) {
        for (uint32_t pos = first; pos <= last; ++pos) {
            if (!shared_boundary_anchor_(
                    target, bridge, candidate, pos, true, bridge_anchor)) continue;
            target_anchor = pos;
            break;
        }
        for (uint32_t pos = first; target_anchor == UINT32_MAX && pos > 0; --pos) {
            if (shared_boundary_anchor_(
                    target, bridge, candidate, pos - 1, true, bridge_anchor)) {
                target_anchor = pos - 1;
            }
        }
    } else {
        for (uint32_t pos = last + 1; target_anchor == UINT32_MAX && pos > first; --pos) {
            if (shared_boundary_anchor_(
                    target, bridge, candidate, pos - 1, false, bridge_anchor)) {
                target_anchor = pos - 1;
            }
        }
        for (uint32_t pos = last + 1;
             target_anchor == UINT32_MAX && pos < target.vertices.size(); ++pos) {
            if (shared_boundary_anchor_(
                    target, bridge, candidate, pos, false, bridge_anchor)) {
                target_anchor = pos;
            }
        }
    }
    if (target_anchor == UINT32_MAX) return false;

    if (before) {
        window.target_begin = target.path_bp[target_anchor];
        window.target_end = interval.second;
        window.bridge_begin = bridge.path_bp[bridge_anchor];
        const uint64_t query_length = window.target_end - window.target_begin;
        const uint64_t reference_length = query_length + (query_length + 9) / 10;
        window.bridge_end = window.bridge_begin + std::min(
            reference_length, bridge.length - window.bridge_begin
        );
    } else {
        window.target_begin = interval.first;
        window.target_end = target.path_bp[target_anchor + 1];
        window.bridge_end = bridge.path_bp[bridge_anchor + 1];
        const uint64_t query_length = window.target_end - window.target_begin;
        const uint64_t reference_length = query_length + (query_length + 9) / 10;
        window.bridge_begin = window.bridge_end > reference_length ?
            window.bridge_end - reference_length : 0;
    }
    return window.target_begin < window.target_end &&
           window.bridge_begin < window.bridge_end;
}

bool GfaGapfill::refine_boundary_(Candidate& candidate) const {
    const Fragment& left = fragments_[candidate.left];
    const Fragment& right = fragments_[candidate.right];
    const Fragment& bridge = fragments_[candidate.bridge];
    candidate.left_cut_bp = left.path_bp[candidate.left_pos + 1];
    candidate.bridge_left_cut_bp = bridge.path_bp[candidate.bridge_left + 1];
    candidate.bridge_right_cut_bp = bridge.path_bp[candidate.bridge_right];
    candidate.right_cut_bp = right.path_bp[candidate.right_pos];

    std::array<BoundaryWindow, 2> windows;
    const bool have_left = boundary_alignment_window_(
        left, bridge, candidate, true, windows[0]
    );
    const bool have_right = boundary_alignment_window_(
        right, bridge, candidate, false, windows[1]
    );
    std::array<BoundaryAlignment, 2> alignments;
    bool mapped_left = false;
    bool mapped_right = false;
    if (have_left) {
        mapped_left = align_boundary_(
            fragment_subsequence_(bridge, windows[0].bridge_begin, windows[0].bridge_end),
            fragment_subsequence_(left, windows[0].target_begin, windows[0].target_end),
            true, alignments[0]
        );
    }
    if (have_right) {
        mapped_right = align_boundary_(
            fragment_subsequence_(bridge, windows[1].bridge_begin, windows[1].bridge_end),
            fragment_subsequence_(right, windows[1].target_begin, windows[1].target_end),
            false, alignments[1]
        );
    }
    if (mapped_left) {
        candidate.left_alignment_similarity = alignments[0].identity;
    }
    if (mapped_right) {
        candidate.right_alignment_similarity = alignments[1].identity;
    }

    const GfaGapfillDebugger::Alignment left_result{
        windows[0].target_begin, windows[0].target_end,
        windows[0].bridge_begin, windows[0].bridge_end,
        alignments[0].query_cut, alignments[0].reference_cut,
        alignments[0].query_begin, alignments[0].query_end,
        alignments[0].reference_begin, alignments[0].reference_end,
        alignments[0].hits, alignments[0].mapq,
        alignments[0].identity, alignments[0].coverage,
        alignments[0].accepted, alignments[0].reverse
    };
    const GfaGapfillDebugger::Alignment right_result{
        windows[1].target_begin, windows[1].target_end,
        windows[1].bridge_begin, windows[1].bridge_end,
        alignments[1].query_cut, alignments[1].reference_cut,
        alignments[1].query_begin, alignments[1].query_end,
        alignments[1].reference_begin, alignments[1].reference_end,
        alignments[1].hits, alignments[1].mapq,
        alignments[1].identity, alignments[1].coverage,
        alignments[1].accepted, alignments[1].reverse
    };
    GfaGapfillDebugger::boundary_alignment(
        paths_[left.path_id].name,
        paths_[right.path_id].name,
        paths_[bridge.path_id].name,
        left_result,
        right_result
    );

    uint64_t left_cut = candidate.left_cut_bp;
    uint64_t bridge_left_cut = candidate.bridge_left_cut_bp;
    uint64_t bridge_right_cut = candidate.bridge_right_cut_bp;
    uint64_t right_cut = candidate.right_cut_bp;
    if (mapped_left) {
        left_cut = windows[0].target_begin + alignments[0].query_cut;
        bridge_left_cut = windows[0].bridge_begin + alignments[0].reference_cut;
    }
    if (mapped_right) {
        bridge_right_cut = windows[1].bridge_begin + alignments[1].reference_cut;
        right_cut = windows[1].target_begin + alignments[1].query_cut;
    }
    if (bridge_left_cut >= bridge_right_cut) return false;
    candidate.left_cut_bp = left_cut;
    candidate.bridge_left_cut_bp = bridge_left_cut;
    candidate.bridge_right_cut_bp = bridge_right_cut;
    candidate.right_cut_bp = right_cut;
    return mapped_left || mapped_right;
}

void GfaGapfill::refine_boundaries_(
    std::vector<Candidate>& selected,
    std::vector<Candidate>& candidates
) const {
    log_stream() << "Refining selected gap boundaries from local alignments ...\n";
    ThreadPool pool(params_.threads);
    std::vector<std::future<bool>> futures;
    futures.reserve(selected.size());
    ProgressTracker progress(selected.size());
    for (Candidate& candidate : selected) {
        Candidate* task = &candidate;
        futures.emplace_back(pool.submit([this, task, &progress]() {
            const bool refined = refine_boundary_(*task);
            progress.hit();
            return refined;
        }));
    }

    size_t refined = 0;
    for (size_t i = 0; i < selected.size(); ++i) {
        const bool changed = futures[i].get();
        const Candidate& candidate = selected[i];
        for (Candidate& item : candidates) {
            if (item.status != Candidate::KEPT || item.left != candidate.left || item.right != candidate.right || item.bridge != candidate.bridge) continue;
            item.left_alignment_similarity = candidate.left_alignment_similarity;
            item.right_alignment_similarity = candidate.right_alignment_similarity;
            item.left_cut_bp = candidate.left_cut_bp;
            item.bridge_left_cut_bp = candidate.bridge_left_cut_bp;
            item.bridge_right_cut_bp = candidate.bridge_right_cut_bp;
            item.right_cut_bp = candidate.right_cut_bp;
            if (changed) {
                item.left_pos = candidate.left_pos;
                item.bridge_left = candidate.bridge_left;
                item.bridge_right = candidate.bridge_right;
                item.right_pos = candidate.right_pos;
                item.left_phase = candidate.left_phase;
                item.right_phase = candidate.right_phase;
                item.left_boundary_similarity = candidate.left_boundary_similarity;
                item.right_boundary_similarity = candidate.right_boundary_similarity;
            }
            break;
        }
        if (changed) ++refined;
        GfaGapfillDebugger::refined_boundary(
            paths_[fragments_[candidate.left].path_id].name,
            candidate.left_cut_bp,
            paths_[fragments_[candidate.right].path_id].name,
            candidate.right_cut_bp, fragments_[candidate.right].length,
            paths_[fragments_[candidate.bridge].path_id].name,
            candidate.bridge_left_cut_bp, candidate.bridge_right_cut_bp,
            {candidate.left_boundary_similarity, candidate.left_phase.similarity,
             candidate.left_phase.target_bp, candidate.left_phase.bridge_bp},
            {candidate.right_boundary_similarity, candidate.right_phase.similarity,
             candidate.right_phase.target_bp, candidate.right_phase.bridge_bp},
            candidate.left_alignment_similarity,
            candidate.right_alignment_similarity
        );
    }
    progress.finish();
    pool.stop();
    log_stream() << "  - Boundaries refined: " << refined << "\n\n";
}

void GfaGapfill::update_support_(
    Candidate& candidate,
    const std::vector<Candidate>& evidence
) const {
    std::vector<std::vector<uint32_t>> fill_paths;
    fill_paths.reserve(evidence.size());
    size_t representative = 0;
    for (size_t i = 0; i < evidence.size(); ++i) {
        const Candidate& witness = evidence[i];
        const Fragment& bridge = fragments_[witness.bridge];
        const uint32_t begin = witness.bridge_left + 1;
        const uint32_t end = witness.bridge_right;
        fill_paths.emplace_back(
            bridge.vertices.begin() + begin,
            bridge.vertices.begin() + end
        );
        if (witness.bridge == candidate.bridge && witness.bridge_left == candidate.bridge_left && witness.bridge_right == candidate.bridge_right) {
            representative = i;
        }
    }

    const GfaBubble::PathClusterer clusterer(
        *this, params_.path_difference, mm2_.k, mm2_.w
    );
    const std::vector<uint32_t> labels = clusterer.cluster(
        fill_paths, GfaBubble::Type::Tip
    );
    const uint32_t target_label = labels[representative];

    std::map<std::string, std::vector<size_t>> by_sample;
    for (size_t i = 0; i < evidence.size(); ++i) {
        by_sample[fragments_[evidence[i].bridge].sample].push_back(i);
    }
    candidate.sample_support = 0;
    candidate.informative_samples = 0;
    candidate.spanning_samples = by_sample.size();

    for (const auto& item : by_sample) {
        const std::vector<size_t>& alleles = item.second;
        uint32_t label = labels[alleles.front()];
        bool homozygous = true;
        for (size_t index : alleles) {
            if (labels[index] != label) {
                homozygous = false;
                break;
            }
        }

        size_t chosen = alleles.front();
        if (!homozygous) {
            bool phased = false;
            bool tied = false;
            for (size_t index : alleles) {
                if (!evidence[index].phase_consistent) continue;
                if (!phased || evidence[index].phase_score > evidence[chosen].phase_score) {
                    chosen = index;
                    phased = true;
                    tied = false;
                } else if (evidence[index].phase_score == evidence[chosen].phase_score && labels[index] != labels[chosen]) {
                    tied = true;
                }
            }
            if (!phased || tied) continue;
            label = labels[chosen];
        }

        ++candidate.informative_samples;
        if (label == target_label) ++candidate.sample_support;
    }

    candidate.probability = candidate.informative_samples == 0 ? 0.0 :
        (1.0 + candidate.sample_support) / (2.0 + candidate.informative_samples);
}

void GfaGapfill::check_unplaced_sequence_(std::vector<Candidate>& candidates, bool retry) const {
    log_stream() << (retry ? "Rechecking" : "Checking") << " long unplaced contig ends ...\n";
    for (Candidate& candidate : candidates) check_unplaced_candidate_(candidate);
    std::cerr << '\n';
}

void GfaGapfill::check_unplaced_candidate_(Candidate& candidate) const {
    const Fragment& left = fragments_[candidate.left];
    const Fragment& right = fragments_[candidate.right];
    const Fragment& bridge = fragments_[candidate.bridge];

    candidate.left_unplaced = left.length - left.path_bp[candidate.junction_left_pos + 1];
    candidate.right_unplaced = right.path_bp[candidate.junction_right_pos];
    if (candidate.left_unplaced > params_.misassembly_check_bp) {
        const std::string unplaced = fragment_sequence_(
            left, candidate.junction_left_pos + 1,
            static_cast<uint32_t>(left.vertices.size())
        );
        const std::string filled = fragment_sequence_(
            bridge, candidate.junction_bridge_left + 1, candidate.bridge_right
        );
        candidate.left_unplaced_similarity = mm2_similarity_(filled, unplaced);
        if (candidate.left_unplaced_similarity < 0.0) {
            warning_stream() << "  ! Preserving unchecked left end of " << paths_[left.path_id].name << ": sequence unavailable\n";
        }
        candidate.left_misassembly =
            candidate.left_unplaced_similarity < params_.misassembly_similarity;
    }
    if (candidate.right_unplaced > params_.misassembly_check_bp) {
        const std::string unplaced = fragment_sequence_(
            right, 0, candidate.junction_right_pos
        );
        const std::string filled = fragment_sequence_(
            bridge, candidate.bridge_left + 1, candidate.junction_bridge_right
        );
        candidate.right_unplaced_similarity = mm2_similarity_(filled, unplaced);
        if (candidate.right_unplaced_similarity < 0.0) {
            warning_stream() << "  ! Preserving unchecked right end of " << paths_[right.path_id].name << ": sequence unavailable\n";
        }
        candidate.right_misassembly =
            candidate.right_unplaced_similarity < params_.misassembly_similarity;
    }

    if (candidate.left_misassembly || candidate.right_misassembly) {
        std::ostringstream message;
        message << "  - Possible misassembly: " << paths_[left.path_id].name << " -> " << paths_[right.path_id].name << " via " << paths_[bridge.path_id].name;
        if (candidate.left_misassembly) {
            message << " left_unplaced=" << candidate.left_unplaced << " left_similarity=" << std::fixed << std::setprecision(4) << candidate.left_unplaced_similarity;
        }
        if (candidate.right_misassembly) {
            message << " right_unplaced=" << candidate.right_unplaced << " right_similarity=" << std::fixed << std::setprecision(4) << candidate.right_unplaced_similarity;
        }
        log_stream() << message.str() << '\n';
    }
}

bool GfaGapfill::relocate_misassemblies_(
    const std::vector<Candidate>& selected,
    std::unordered_set<uint32_t>& affected_components
) {
    log_stream() << "Relocating misassembled contig ends ...\n";

    struct Proposal {
        uint32_t source{UINT32_MAX};
        uint32_t begin{0}, end{0};
        uint32_t record{UINT32_MAX};
        uint32_t backbone{UINT32_MAX};
        std::vector<uint32_t> vertices;
        int64_t layout_start{0};
        bool accepted{false};
    };
    struct Backbone {
        uint32_t component{UINT32_MAX};
        uint32_t fragment{UINT32_MAX};
        std::string sequence;
        std::vector<uint64_t> offsets;
    };

    std::vector<Proposal> proposals;
    std::unordered_set<uint64_t> seen;
    for (const Candidate& candidate : selected) {
        const uint32_t fragment_ids[2] = {candidate.left, candidate.right};
        const uint32_t begins[2] = {
            candidate.junction_left_pos + 1, 0
        };
        const uint32_t ends[2] = {
            static_cast<uint32_t>(fragments_[candidate.left].vertices.size()),
            candidate.junction_right_pos
        };
        const bool flagged[2] = {
            candidate.left_misassembly, candidate.right_misassembly
        };
        for (uint32_t side = 0; side < 2; ++side) {
            if (!flagged[side]) continue;
            const uint64_t key = (static_cast<uint64_t>(fragment_ids[side]) << 1) | side;
            if (!seen.insert(key).second) continue;

            const Fragment& source = fragments_[fragment_ids[side]];
            RelocationRecord record;
            record.sample = source.sample;
            record.hap = source.hap;
            record.source = paths_[source.path_id].name;
            record.side = side == 0 ? "right" : "left";
            record.source_component = source.component;
            record.status = "unmapped";
            record.unplaced_vertices.assign(
                source.vertices.begin() + begins[side],
                source.vertices.begin() + ends[side]
            );
            relocations_.push_back(std::move(record));
            proposals.push_back({
                fragment_ids[side], begins[side], ends[side],
                static_cast<uint32_t>(relocations_.size() - 1)
            });
        }
    }
    if (proposals.empty()) {
        log_stream() << "  - Relocated misassembled ends: " << 0 << '\n';
        log_stream() << "  - Components scheduled for gap-fill retry: " << 0 << "\n\n";
        return false;
    }

    std::map<uint32_t, uint32_t> best_fragment;
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible) continue;
        const auto found = best_fragment.find(fragment.component);
        if (found != best_fragment.end() && fragments_[found->second].length >= fragment.length) continue;
        bool sequence_available = true;
        for (uint32_t vertex : fragment.vertices) {
            const std::string& sequence = nodes_[Vertex::get_segment_id(vertex)].sequence;
            if (sequence.empty() || sequence == "*") {
                sequence_available = false;
                break;
            }
        }
        if (sequence_available) best_fragment[fragment.component] = i;
    }

    std::vector<Backbone> backbones;
    backbones.reserve(best_fragment.size());
    for (const auto& item : best_fragment) {
        const Fragment& fragment = fragments_[item.second];
        Backbone backbone;
        backbone.component = item.first;
        backbone.fragment = item.second;
        backbone.sequence = fragment_sequence_(fragment, 0, fragment.vertices.size());
        backbone.offsets.reserve(fragment.vertices.size() + 1);
        backbone.offsets.push_back(0);
        for (uint32_t i = 0; i < fragment.vertices.size(); ++i) {
            const uint32_t overlap = i == 0 ? 0 : get_edge_ow(fragment.vertices[i - 1], fragment.vertices[i]);
            const uint32_t length = nodes_[Vertex::get_segment_id(fragment.vertices[i])].length;
            backbone.offsets.push_back(backbone.offsets.back() + (length > overlap ? length - overlap : 0));
        }
        backbones.push_back(std::move(backbone));
    }

    std::unordered_set<uint32_t> source_components;
    for (const Proposal& proposal : proposals) {
        source_components.insert(fragments_[proposal.source].component);
    }
    for (uint32_t source_component : source_components) {
        std::vector<uint32_t> reference_ids;
        std::vector<const char*> sequences;
        std::vector<const char*> names;
        for (uint32_t i = 0; i < backbones.size(); ++i) {
            if (backbones[i].component == source_component) continue;
            reference_ids.push_back(i);
            sequences.push_back(backbones[i].sequence.c_str());
            names.push_back(paths_[fragments_[backbones[i].fragment].path_id].name.c_str());
        }
        if (sequences.empty()) continue;

        mm_idxopt_t index_options;
        mm_mapopt_t map_options;
        mm_set_opt(nullptr, &index_options, &map_options);
        if (mm_set_opt(mm2_.preset.c_str(), &index_options, &map_options) < 0) continue;
        index_options.k = static_cast<short>(mm2_.k);
        index_options.w = static_cast<short>(mm2_.w);
        map_options.flag |= MM_F_CIGAR | MM_F_EQX;
        map_options.best_n = static_cast<short>(mm2_.best_n);
        map_options.zdrop = mm2_.zdrop;
        const int hpc = (index_options.flag & MM_I_HPC) ? 1 : 0;
        mm_idx_t* index = mm_idx_str(
            index_options.w, index_options.k, hpc, index_options.bucket_bits,
            static_cast<int>(sequences.size()), sequences.data(), names.data()
        );
        if (!index) continue;
        mm_mapopt_update(&map_options, index);
        mm_tbuf_t* buffer = mm_tbuf_init();
        if (!buffer) {
            mm_idx_destroy(index);
            continue;
        }

        for (Proposal& proposal : proposals) {
            const Fragment& source = fragments_[proposal.source];
            if (source.component != source_component) continue;
            RelocationRecord& record = relocations_[proposal.record];
            const std::string query = fragment_sequence_(source, proposal.begin, proposal.end);
            if (query.empty()) {
                record.status = "sequence-unavailable";
                continue;
            }

            int count = 0;
            mm_reg1_t* hits = mm_map(
                index, static_cast<int>(query.size()), query.c_str(),
                &count, buffer, &map_options, record.source.c_str()
            );
            const mm_reg1_t* best = nullptr;
            for (int i = 0; i < count; ++i) {
                const mm_reg1_t& hit = hits[i];
                if (!hit.p || hit.parent != hit.id || hit.mapq < mm2_.min_mapq || hit.blen <= 0 || hit.rid < 0 || static_cast<size_t>(hit.rid) >= reference_ids.size()) continue;
                const double match_ratio = static_cast<double>(hit.mlen) / hit.blen;
                const Backbone& reference = backbones[reference_ids[hit.rid]];
                const uint64_t shorter = std::min(query.size(), reference.sequence.size());
                const uint32_t aligned_span = std::max(hit.re - hit.rs, hit.qe - hit.qs);
                const double aligned_ratio = shorter == 0 ? 0.0 :
                    static_cast<double>(aligned_span) / shorter;
                if (match_ratio < mm2_.min_match || aligned_ratio < mm2_.min_ali_ratio) continue;
                if (!best || hit.mlen > best->mlen ||
                    (hit.mlen == best->mlen && hit.mapq > best->mapq)) best = &hit;
            }

            if (best) {
                const Backbone& backbone = backbones[reference_ids[best->rid]];
                const Fragment& target = fragments_[backbone.fragment];
                const auto first_it = std::upper_bound(
                    backbone.offsets.begin(), backbone.offsets.end(), best->rs
                );
                const auto last_it = std::lower_bound(
                    backbone.offsets.begin(), backbone.offsets.end(), best->re
                );
                uint32_t first = first_it == backbone.offsets.begin() ? 0 : static_cast<uint32_t>(first_it - backbone.offsets.begin() - 1);
                uint32_t last = static_cast<uint32_t>(last_it - backbone.offsets.begin());
                first = std::min(first, static_cast<uint32_t>(target.vertices.size()));
                last = std::min(last, static_cast<uint32_t>(target.vertices.size()));
                if (first < last) {
                    proposal.backbone = backbone.fragment;
                    proposal.vertices.assign(
                        target.vertices.begin() + first, target.vertices.begin() + last
                    );
                    proposal.layout_start = target.layout_start +
                        static_cast<int64_t>(backbone.offsets[first]);
                    proposal.accepted = true;
                    record.status = "relocated";
                    record.target_component = backbone.component;
                    record.target_backbone = paths_[target.path_id].name;
                    record.target_beg = best->rs;
                    record.target_end = best->re;
                    record.strand = best->rev ? "-" : "+";
                    const double identity = static_cast<double>(best->mlen) / best->blen;
                    const uint64_t shorter = std::min(query.size(), backbone.sequence.size());
                    const uint32_t span = std::max(best->re - best->rs, best->qe - best->qs);
                    const double coverage = shorter == 0 ? 0.0 : std::min(1.0, static_cast<double>(span) / shorter);
                    const double mapping_confidence = 1.0 - std::pow(10.0, -static_cast<double>(best->mapq) / 10.0);
                    record.probability = identity * coverage * mapping_confidence;
                }
            }
            if (hits) {
                for (int i = 0; i < count; ++i) std::free(hits[i].p);
                std::free(hits);
            }
        }
        mm_tbuf_destroy(buffer);
        mm_idx_destroy(index);
    }

    std::vector<uint32_t> retained_begin(fragments_.size(), 0);
    std::vector<uint32_t> retained_end;
    retained_end.reserve(fragments_.size());
    for (const Fragment& fragment : fragments_) {
        retained_end.push_back(static_cast<uint32_t>(fragment.vertices.size()));
    }
    for (const Proposal& proposal : proposals) {
        if (!proposal.accepted) continue;
        if (proposal.begin == 0) retained_begin[proposal.source] =
            std::max(retained_begin[proposal.source], proposal.end);
        if (proposal.end == fragments_[proposal.source].vertices.size()) {
            retained_end[proposal.source] = std::min(retained_end[proposal.source], proposal.begin);
        }
    }
    for (Proposal& proposal : proposals) {
        if (!proposal.accepted) continue;
        if (retained_begin[proposal.source] >= retained_end[proposal.source]) {
            retained_begin[proposal.source] = 0;
            retained_end[proposal.source] =
                static_cast<uint32_t>(fragments_[proposal.source].vertices.size());
            for (Proposal& other : proposals) {
                if (other.source != proposal.source) continue;
                other.accepted = false;
                relocations_[other.record].status = "conflicting-trims";
            }
        }
    }

    std::unordered_set<std::string> path_names;
    for (const GfaPath& path : paths_) path_names.insert(path.name);
    uint32_t relocated_count = 0;
    for (Proposal& proposal : proposals) {
        if (!proposal.accepted) continue;
        Fragment& source = fragments_[proposal.source];
        const RelocationRecord& record = relocations_[proposal.record];
        std::string name = paths_[source.path_id].name + ".reloc" + (proposal.begin == 0 ? "L" : "R");
        const std::string base = name;
        while (!path_names.insert(name).second) {
            name = base + std::to_string(++relocated_count);
        }
        relocations_[proposal.record].relocated_contig = name;

        GfaPath path;
        path.name = name;
        path.segments.reserve(proposal.vertices.size());
        for (uint32_t vertex : proposal.vertices) {
            path.segments.push_back({
                Vertex::get_segment_id(vertex), Vertex::get_is_reverse(vertex)
            });
        }
        paths_.push_back(std::move(path));

        Fragment relocated;
        relocated.path_id = paths_.size() - 1;
        relocated.sample = source.sample;
        relocated.component = record.target_component;
        relocated.hap = source.hap;
        relocated.layout_start = proposal.layout_start;
        relocated.vertices = std::move(proposal.vertices);
        relocated.relocation = static_cast<int32_t>(proposal.record);
        rebuild_fragment_index_(relocated);
        relocated.eligible = relocated.length >= params_.min_contig_bp;
        const uint32_t source_component = source.component;
        fragments_.push_back(std::move(relocated));
        affected_components.insert(source_component);
        affected_components.insert(record.target_component);
    }

    const size_t original_fragments = retained_end.size();
    for (uint32_t i = 0; i < original_fragments; ++i) {
        if (retained_begin[i] == 0 && retained_end[i] == fragments_[i].vertices.size()) continue;
        Fragment& source = fragments_[i];
        const uint64_t removed_prefix = source.path_bp[retained_begin[i]];
        std::vector<uint32_t> retained(
            source.vertices.begin() + retained_begin[i],
            source.vertices.begin() + retained_end[i]
        );
        source.vertices = retained;
        source.layout_start += static_cast<int64_t>(removed_prefix);
        rebuild_fragment_index_(source);

        if (source.reverse) {
            std::reverse(retained.begin(), retained.end());
            for (uint32_t& vertex : retained) vertex ^= 1u;
        }
        GfaPath& path = paths_[source.path_id];
        path.segments.clear();
        path.segments.reserve(retained.size());
        for (uint32_t vertex : retained) {
            path.segments.push_back({
                Vertex::get_segment_id(vertex), Vertex::get_is_reverse(vertex)
            });
        }
    }

    if (affected_components.empty()) {
        log_stream() << "  - Relocated misassembled ends: " << 0 << '\n';
        log_stream() << "  - Components scheduled for gap-fill retry: " << 0 << "\n\n";
        return false;
    }

    size_t accepted = 0;
    for (const Proposal& proposal : proposals) accepted += proposal.accepted;
    log_stream() << "  - Relocated misassembled ends: " << accepted << '\n';
    log_stream() << "  - Components scheduled for gap-fill retry: " << affected_components.size() << "\n\n";
    return true;
}

void GfaGapfill::mark_used_relocations_(const std::vector<Candidate>& selected) {
    for (const Candidate& candidate : selected) {
        if (candidate.bridge_right <= candidate.bridge_left + 1) continue;
        const uint32_t ids[3] = {candidate.left, candidate.right, candidate.bridge};
        for (uint32_t id : ids) {
            const int32_t relocation = fragments_[id].relocation;
            if (relocation >= 0) relocations_[relocation].used_for_gap = true;
        }
    }
}

void GfaGapfill::mark_redundant_fragments_(const std::vector<Candidate>& selected) {
    log_stream() << "Checking contigs covered by filled gaps ...\n";
    std::vector<uint8_t> used(fragments_.size(), 0);
    using Group = std::tuple<std::string, uint32_t, uint8_t>;
    std::map<Group, std::vector<size_t>> gaps_by_group;
    for (size_t candidate_id = 0; candidate_id < selected.size(); ++candidate_id) {
        const Candidate& candidate = selected[candidate_id];
        used[candidate.left] = 1;
        used[candidate.right] = 1;
        used[candidate.bridge] = 1;
        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        if (left.hap == right.hap) {
            gaps_by_group[{left.sample, left.component, left.hap}].push_back(candidate_id);
        }
    }

    std::vector<std::string> gap_sequences(selected.size());
    size_t redundant = 0;
    uint64_t redundant_bp = 0;
    for (uint32_t fragment_id = 0; fragment_id < fragments_.size(); ++fragment_id) {
        Fragment& fragment = fragments_[fragment_id];
        if (!fragment.eligible || used[fragment_id]) continue;
        const auto group = gaps_by_group.find({fragment.sample, fragment.component, fragment.hap});
        if (group == gaps_by_group.end()) continue;

        const int64_t fragment_begin = fragment.layout_start;
        const int64_t fragment_end = fragment.layout_start + static_cast<int64_t>(fragment.length);
        std::string sequence;
        bool sequence_loaded = false;
        for (size_t candidate_id : group->second) {
            const Candidate& candidate = selected[candidate_id];
            const Fragment& bridge = fragments_[candidate.bridge];
            if (candidate.bridge_right <= candidate.bridge_left + 1) continue;

            const int64_t gap_begin = bridge.layout_start +
                static_cast<int64_t>(candidate.bridge_left_cut_bp);
            const int64_t gap_end = bridge.layout_start +
                static_cast<int64_t>(candidate.bridge_right_cut_bp);
            const int64_t overlap_begin = std::max(fragment_begin, gap_begin);
            const int64_t overlap_end = std::min(fragment_end, gap_end);
            if (overlap_end <= overlap_begin || static_cast<uint64_t>(overlap_end - overlap_begin) * 10 < fragment.length * 8) continue;

            std::string& gap_sequence = gap_sequences[candidate_id];
            if (gap_sequence.empty()) {
                gap_sequence = fragment_subsequence_(
                    bridge, candidate.bridge_left_cut_bp,
                    candidate.bridge_right_cut_bp
                );
            }
            if (!sequence_loaded) {
                sequence = fragment_sequence_(
                    fragment, 0, static_cast<uint32_t>(fragment.vertices.size())
                );
                sequence_loaded = true;
            }
            const double similarity = mm2_similarity_(gap_sequence, sequence);
            if (similarity < params_.dedup_similarity) continue;

            fragment.redundant = true;
            ++redundant;
            redundant_bp += fragment.length;
            log_stream() << "  - Covered contig: " << paths_[fragment.path_id].name << " by " << paths_[bridge.path_id].name << " (similarity=" << std::fixed << std::setprecision(4) << similarity << ")\n";
            break;
        }
    }
    log_stream() << "  - Redundant contigs removed: " << redundant << " (" << redundant_bp << " bp)\n\n";
}

std::vector<GfaGapfill::Candidate> GfaGapfill::build_candidates_(
    const std::unordered_set<uint32_t>* components
) const {
    std::vector<std::string> group_keys(fragments_.size());
    std::map<std::string, std::vector<uint32_t>> fragments_by_sample;
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible || fragment.redundant) continue;
        if (components && components->find(fragment.component) == components->end()) continue;
        fragments_by_sample[fragment.sample].push_back(i);
        group_keys[i] = fragment.sample + '\t' + std::to_string(fragment.component);
    }

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> node_paths;
    size_t path_nodes = 0;
    for (const Fragment& fragment : fragments_) {
        if (!fragment.eligible) continue;
        if (!components || components->find(fragment.component) != components->end()) {
            path_nodes += fragment.vertices.size();
        }
    }
    node_paths.reserve(path_nodes);
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        if (!fragments_[i].eligible) continue;
        if (components && components->find(fragments_[i].component) == components->end()) continue;
        for (uint32_t pos = 0; pos < fragments_[i].vertices.size(); ++pos) {
            node_paths.emplace_back(
                Vertex::get_segment_id(fragments_[i].vertices[pos]), i, pos
            );
        }
    }
    std::sort(node_paths.begin(), node_paths.end());

    OverlapIndex overlaps(fragments_.size());
    for (size_t begin = 0; begin < node_paths.size();) {
        size_t end = begin + 1;
        while (end < node_paths.size() &&
               std::get<0>(node_paths[end]) == std::get<0>(node_paths[begin])) {
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
                if (sid < bubble_nodes_.size() && bubble_nodes_[sid]) {
                    ab.bubble_bp += nodes_[sid].length;
                    ba.bubble_bp += nodes_[sid].length;
                }
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

    if (DEBUG_ENABLED) {
        for (uint32_t a = 0; a < overlaps.size(); ++a) {
            for (const auto& item : overlaps[a]) {
                const uint32_t b = item.first;
                if (a >= b) continue;
                const OverlapSupport support = overlap_support_(fragments_[a], fragments_[b], item.second);
                const PathOverlap& overlap = item.second;
                const BoundaryPhaseSupport left_phase = boundary_phase_similarity_(
                    fragments_[a], overlap.a_begin,
                    fragments_[b], overlap.b_begin,
                    false
                );
                const BoundaryPhaseSupport right_phase = boundary_phase_similarity_(
                    fragments_[a], overlap.a_end,
                    fragments_[b], overlap.b_end,
                    true
                );
                const GfaGapfillDebugger::BoundarySupport left_boundary{
                    boundary_similarity_(
                        fragments_[a], overlap.a_begin,
                        fragments_[b], overlap.b_begin,
                        false, left_phase.windows
                    ),
                    left_phase.similarity,
                    left_phase.target_bp,
                    left_phase.bridge_bp
                };
                const GfaGapfillDebugger::BoundarySupport right_boundary{
                    boundary_similarity_(
                        fragments_[a], overlap.a_end,
                        fragments_[b], overlap.b_end,
                        true, right_phase.windows
                    ),
                    right_phase.similarity,
                    right_phase.target_bp,
                    right_phase.bridge_bp
                };
                GfaGapfillDebugger::overlap(
                    paths_[fragments_[a].path_id].name,
                    paths_[fragments_[b].path_id].name,
                    support.shared_bp, support.overlap_bp,
                    support.similarity,
                    fragments_[a].length == 0 ? 0.0 : static_cast<double>(support.shared_bp) / fragments_[a].length,
                    fragments_[b].length == 0 ? 0.0 : static_cast<double>(support.shared_bp) / fragments_[b].length,
                    left_boundary, right_boundary
                );
            }
        }
    }

    struct CandidateBatch {
        std::vector<Candidate> candidates;
        uint64_t tested{0}, low_support{0}, wrong_group{0}, anchor_failed{0};
    };

    std::vector<std::vector<uint32_t>> sample_groups;
    sample_groups.reserve(fragments_by_sample.size());
    for (auto& item : fragments_by_sample) sample_groups.push_back(std::move(item.second));

    log_stream() << "Building gap-fill candidates by sample ...\n";
    ThreadPool pool(params_.threads);
    std::vector<std::future<CandidateBatch>> futures;
    futures.reserve(sample_groups.size());
    ProgressTracker progress(sample_groups.size());

    for (size_t sample_id = 0; sample_id < sample_groups.size(); ++sample_id) {
        futures.emplace_back(pool.submit([&, sample_id]() {
            CandidateBatch batch;
            for (uint32_t left_id : sample_groups[sample_id]) {
                const Fragment& left = fragments_[left_id];
                if (left.length < params_.min_contig_bp) continue;
                std::unordered_map<uint32_t, Candidate> by_right;
                std::unordered_map<uint32_t, std::vector<Candidate>> evidence_by_right;

                for (const auto& bridge_item : overlaps[left_id]) {
                    const uint32_t bridge_id = bridge_item.first;
                    const Fragment& bridge = fragments_[bridge_id];
                    if (bridge.sample == left.sample) continue;
                    const OverlapSupport left_support = overlap_support_(
                        left, bridge, bridge_item.second
                    );
                    const double left_fraction = static_cast<double>(left_support.shared_bp) /
                        std::min(left.length, bridge.length);
                    if (left_support.overlap_bp < params_.min_overlap_bp ||
                        (left_support.similarity < params_.min_similarity &&
                         !gap_boundary_supported_(left, bridge, bridge_item.second, true))) {
                        ++batch.low_support;
                        continue;
                    }

                    for (const auto& right_item : overlaps[bridge_id]) {
                        const uint32_t right_id = right_item.first;
                        ++batch.tested;
                        const Fragment& right = fragments_[right_id];
                        const OverlapSupport right_support = overlap_support_(
                            bridge, right, right_item.second
                        );
                        const double right_fraction = static_cast<double>(right_support.shared_bp) /
                            std::min(right.length, bridge.length);
                        if (right_support.overlap_bp < params_.min_overlap_bp ||
                            (right_support.similarity < params_.min_similarity &&
                             !gap_boundary_supported_(bridge, right, right_item.second, false))) {
                            ++batch.low_support;
                            continue;
                        }
                        if (right_id == left_id || group_keys[right_id] != group_keys[left_id] ||
                            !ordered_pair_(left, right)) {
                            ++batch.wrong_group;
                            continue;
                        }
                        const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
                        const int64_t bridge_end = bridge.layout_start + static_cast<int64_t>(bridge.length);
                        if (bridge.layout_start > left_end || bridge_end < right.layout_start) {
                            ++batch.wrong_group;
                            continue;
                        }

                        Candidate candidate;
                        candidate.left = left_id;
                        candidate.right = right_id;
                        candidate.bridge = bridge_id;
                        const int64_t target_gap = right.layout_start - left_end;
                        candidate.target_overlaps = target_gap < 0;
                        candidate.target_distance = static_cast<uint64_t>(
                            target_gap < 0 ? -target_gap : target_gap
                        );
                        candidate.left_shared = left_support.shared_bp;
                        candidate.right_shared = right_support.shared_bp;
                        candidate.left_similarity = left_support.phase_similarity;
                        candidate.right_similarity = right_support.phase_similarity;
                        if (!find_anchors_(left, right, bridge, candidate)) {
                            ++batch.anchor_failed;
                            continue;
                        }

                        candidate.left_phase = boundary_phase_similarity_(
                            left, candidate.left_pos, bridge, candidate.bridge_left, true
                        );
                        candidate.right_phase = boundary_phase_similarity_(
                            right, candidate.right_pos, bridge, candidate.bridge_right, false
                        );
                        const double left_boundary = boundary_similarity_(
                            left, candidate.left_pos, bridge, candidate.bridge_left,
                            true, candidate.left_phase.windows
                        );
                        const double right_boundary = boundary_similarity_(
                            right, candidate.right_pos, bridge, candidate.bridge_right,
                            false, candidate.right_phase.windows
                        );
                        candidate.left_boundary_similarity = left_boundary;
                        candidate.right_boundary_similarity = right_boundary;
                        if (candidate.left_phase.similarity >= 0.0) {
                            candidate.left_similarity = candidate.left_phase.similarity;
                        } else if (left_boundary >= 0.0) {
                            candidate.left_similarity = left_boundary;
                        }
                        if (candidate.right_phase.similarity >= 0.0) {
                            candidate.right_similarity = candidate.right_phase.similarity;
                        } else if (right_boundary >= 0.0) {
                            candidate.right_similarity = right_boundary;
                        }

                        candidate.score = std::sqrt(
                            candidate.left_similarity * candidate.right_similarity *
                            left_fraction * right_fraction
                        );
                        candidate.phase_score =
                            phase_odds(candidate.left_similarity) *
                            phase_odds(candidate.right_similarity) *
                            std::sqrt(left_fraction * right_fraction);
                        evidence_by_right[right_id].push_back(candidate);
                    }
                }

                for (auto& item : evidence_by_right) {
                    const bool homolog_span = homolog_spans_(left_id, item.first, overlaps);
                    std::map<std::string, std::pair<double, double>> best_phase;
                    for (const Candidate& candidate : item.second) {
                        auto& best = best_phase[fragments_[candidate.bridge].sample];
                        best.first = std::max(best.first, candidate.left_similarity);
                        best.second = std::max(best.second, candidate.right_similarity);
                    }

                    Candidate best;
                    for (Candidate& candidate : item.second) {
                        const auto& phase = best_phase[fragments_[candidate.bridge].sample];
                        candidate.phase_consistent = candidate.left_similarity + 1e-12 >= phase.first && candidate.right_similarity + 1e-12 >= phase.second;
                        if (best.left == UINT32_MAX || candidate_better_(candidate, best)) {
                            best = candidate;
                        }
                    }
                    update_support_(best, item.second);
                    best.homolog_span = homolog_span;
                    by_right[item.first] = std::move(best);
                }

                for (auto& item : by_right) {
                    batch.candidates.push_back(std::move(item.second));
                }
            }
            progress.hit();
            return batch;
        }));
    }

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

    for (const Candidate& candidate : candidates) {
        GfaGapfillDebugger::candidate(
            paths_[fragments_[candidate.left].path_id].name,
            paths_[fragments_[candidate.right].path_id].name,
            paths_[fragments_[candidate.bridge].path_id].name,
            fragments_[candidate.left].component,
            candidate.left_phase.similarity, candidate.right_phase.similarity,
            candidate.left_phase.target_bp, candidate.left_phase.bridge_bp,
            candidate.right_phase.target_bp, candidate.right_phase.bridge_bp,
            candidate.left_boundary_similarity, candidate.right_boundary_similarity,
            candidate.phase_consistent, candidate.homolog_span,
            candidate.sample_support, candidate.informative_samples,
            candidate.spanning_samples,
            candidate.probability
        );
    }
    log_stream() << "  - Candidate connections: " << candidates.size() << "\n\n";
    return candidates;
}

std::vector<GfaGapfill::Candidate> GfaGapfill::select_candidates_(
    std::vector<Candidate>& candidates
) {
    log_stream() << "Selecting unambiguous gap-fill connections ...\n";
    std::sort(candidates.begin(), candidates.end(), connection_better_);

    std::vector<Candidate> selected;
    std::vector<int32_t> incoming(fragments_.size(), -1), outgoing(fragments_.size(), -1);
    std::vector<std::unordered_map<std::string, LocalPhaseAssignments>> bridge_phase(fragments_.size());
    DisjointSet sets(fragments_.size());
    for (Candidate& candidate : candidates) {
        const std::string& left_name = paths_[fragments_[candidate.left].path_id].name;
        const std::string& right_name = paths_[fragments_[candidate.right].path_id].name;
        const std::string& bridge_name = paths_[fragments_[candidate.bridge].path_id].name;
        if (outgoing[candidate.left] >= 0 || incoming[candidate.right] >= 0) {
            candidate.status = Candidate::USED_END;
            GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "drop:used-end", candidate.probability);
            continue;
        }

        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        const uint8_t candidate_hap = left.hap;
        LocalPhaseAssignments& local_phase = bridge_phase[candidate.bridge][left.sample];
        if (candidate_hap != 0 && overlaps_phase_interval(
                local_phase.hap[2 - candidate_hap],
                candidate.bridge_left,
                candidate.bridge_right
            )) {
            candidate.status = Candidate::PHASE_CONFLICT;
            GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "drop:phase-conflict", candidate.probability);
            continue;
        }

        uint32_t left_entry = 0;
        if (incoming[candidate.left] >= 0) {
            left_entry = selected[incoming[candidate.left]].right_pos;
        }
        uint32_t right_exit = static_cast<uint32_t>(fragments_[candidate.right].vertices.size() - 1);
        if (outgoing[candidate.right] >= 0) {
            right_exit = selected[outgoing[candidate.right]].left_pos;
        }
        if (left_entry > candidate.left_pos || candidate.right_pos > right_exit) {
            candidate.status = Candidate::COORDINATE_CONFLICT;
            GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "drop:coordinate-conflict", candidate.probability);
            continue;
        }

        if (sets.find(candidate.left) == sets.find(candidate.right)) {
            candidate.status = Candidate::CYCLE;
            GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "drop:cycle", candidate.probability);
            continue;
        }

        if (candidate.probability < params_.min_probability && !candidate.homolog_span) {
            candidate.status = Candidate::LOW_PROBABILITY;
            GfaGapfillDebugger::selection(
                left_name, right_name, bridge_name,
                "drop:low-probability", candidate.probability
            );
            continue;
        }

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
                candidate.bridge_left,
                candidate.bridge_right
            );
        }
        GfaGapfillDebugger::selection(left_name, right_name, bridge_name, "keep", candidate.probability);
    }
    log_stream() << "  - Connections selected: " << selected.size() << "\n\n";
    return selected;
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
        if (fragment.hap == 1) ++chain.hap1_votes;
        else ++chain.hap2_votes;
        const uint32_t begin = incoming[current] < 0 ? 0 : selected[incoming[current]].right_pos;
        const uint32_t end = outgoing[current] < 0 ?
            static_cast<uint32_t>(fragment.vertices.size() - 1) :
            selected[outgoing[current]].left_pos;
        const uint64_t begin_bp = incoming[current] < 0 ? 0 :
            selected[incoming[current]].right_cut_bp;
        const uint64_t end_bp = outgoing[current] < 0 ? fragment.length :
            selected[outgoing[current]].left_cut_bp;
        if (begin_bp < end_bp) {
            chain.parts.push_back({current, begin_bp, end_bp, false});
        }
        for (uint32_t i = begin; i <= end; ++i) {
            chain.vertices.push_back(fragment.vertices[i]);
        }
        if (outgoing[current] < 0) break;

        const Candidate& edge = selected[outgoing[current]];
        const Fragment& bridge = fragments_[edge.bridge];
        if (edge.bridge_left_cut_bp < edge.bridge_right_cut_bp) {
            chain.parts.push_back({
                edge.bridge, edge.bridge_left_cut_bp,
                edge.bridge_right_cut_bp, true
            });
        }
        for (uint32_t i = edge.bridge_left + 1; i < edge.bridge_right; ++i) {
            chain.vertices.push_back(bridge.vertices[i]);
        }
        chain.gaps.push_back(edge);
        current = edge.right;
    }
    return chain;
}

void GfaGapfill::build_sample_chains_(const std::vector<Candidate>& selected) {
    log_stream() << "Building gap-filled sample paths ...\n\n";
    std::vector<int32_t> incoming(fragments_.size(), -1), outgoing(fragments_.size(), -1);
    for (size_t i = 0; i < selected.size(); ++i) {
        incoming[selected[i].right] = static_cast<int32_t>(i);
        outgoing[selected[i].left] = static_cast<int32_t>(i);
    }

    std::map<std::pair<std::string, uint32_t>, std::vector<Chain>> groups;
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        if (incoming[i] >= 0) continue;

        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible || fragment.redundant) continue;
        Chain chain = build_chain_(i, incoming, outgoing, selected);
        if (chain.vertices.empty()) continue;
        groups[{fragment.sample, fragment.component}].push_back(std::move(chain));
    }

    for (auto& group : groups) {
        std::vector<Chain>& chains = group.second;
        std::vector<std::pair<uint32_t, uint32_t>> node_chains;
        size_t total_nodes = 0;
        for (const Chain& chain : chains) total_nodes += chain.vertices.size();
        node_chains.reserve(total_nodes);
        for (uint32_t i = 0; i < chains.size(); ++i) {
            for (uint32_t vertex : chains[i].vertices) {
                node_chains.emplace_back(Vertex::get_segment_id(vertex), i);
            }
        }
        std::sort(node_chains.begin(), node_chains.end());

        const uint64_t min_shared = static_cast<uint64_t>(
            std::ceil(params_.min_overlap_bp * params_.min_similarity)
        );
        std::unordered_map<uint64_t, uint64_t> shared_bp;
        for (size_t begin = 0; begin < node_chains.size();) {
            size_t end = begin + 1;
            while (end < node_chains.size() &&
                   node_chains[end].first == node_chains[begin].first) {
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

        std::vector<std::pair<int64_t, uint32_t>> chain_order;
        std::vector<std::pair<int64_t, int64_t>> chain_bounds(chains.size());
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
            const bool mixed_a = chains[a].hap1_votes > 0 && chains[a].hap2_votes > 0;
            for (size_t j = i + 1; j < chain_order.size(); ++j) {
                const uint32_t b = chain_order[j].second;
                if (chain_bounds[b].first >= chain_bounds[a].second) break;
                const bool mixed_b = chains[b].hap1_votes > 0 && chains[b].hap2_votes > 0;
                if (!mixed_a && !mixed_b) continue;
                const int64_t overlap = std::min(chain_bounds[a].second, chain_bounds[b].second) -
                    chain_bounds[b].first;
                const int64_t shorter = std::min(
                    chain_bounds[a].second - chain_bounds[a].first,
                    chain_bounds[b].second - chain_bounds[b].first
                );
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

        std::vector<std::vector<uint32_t>> conflicts(chains.size());
        for (const auto& item : shared_bp) {
            if (item.second < min_shared) continue;
            const uint32_t a = static_cast<uint32_t>(item.first >> 32);
            const uint32_t b = static_cast<uint32_t>(item.first);
            conflicts[a].push_back(b);
            conflicts[b].push_back(a);
        }

        std::vector<int8_t> color(chains.size(), -1);
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
            for (uint32_t i : block) {
                if (color[i] == 0) {
                    normal += chains[i].hap1_votes;
                    swapped += chains[i].hap2_votes;
                } else {
                    normal += chains[i].hap2_votes;
                    swapped += chains[i].hap1_votes;
                }
            }
            const bool swap = swapped > normal;
            for (uint32_t i : block) {
                const uint8_t assigned = static_cast<uint8_t>(color[i]) ^ swap;
                Chain& chain = chains[i];
                chain.hap = assigned + 1;
                GfaGapfillDebugger::assignment(
                    chain.sample, chain.component, chain.hap,
                    chain.hap1_votes, chain.hap2_votes, conflict
                );
                sample_chains_[chain.sample][assigned].push_back(std::move(chain));
            }
        }
    }

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

size_t GfaGapfill::gapfill(const std::vector<GfaBubble::Bubble>& bubbles) {
    log_stream() << "Gap-filling sample contigs with cross-sample bridging paths ...\n\n";
    fragments_.clear();
    sample_chains_.clear();
    records_.clear();
    relocations_.clear();
    index_bubble_nodes_(bubbles);
    build_fragments_();
    build_graph_order_();
    std::vector<Candidate> candidates = build_candidates_();
    std::vector<Candidate> selected = select_candidates_(candidates);
    check_unplaced_sequence_(selected);

    std::unordered_set<uint32_t> affected_components;
    if (relocate_misassemblies_(selected, affected_components)) {
        std::vector<Candidate> unchanged_candidates;
        std::vector<Candidate> unchanged_selected;
        unchanged_candidates.reserve(candidates.size());
        unchanged_selected.reserve(selected.size());
        for (const Candidate& candidate : candidates) {
            if (affected_components.find(fragments_[candidate.left].component) == affected_components.end()) {
                unchanged_candidates.push_back(candidate);
            }
        }
        for (const Candidate& candidate : selected) {
            if (affected_components.find(fragments_[candidate.left].component) == affected_components.end()) {
                unchanged_selected.push_back(candidate);
            }
        }

        std::vector<Candidate> retry_candidates = build_candidates_(&affected_components);
        std::vector<Candidate> retry_selected = select_candidates_(retry_candidates);
        check_unplaced_sequence_(retry_selected, true);
        unchanged_candidates.insert(unchanged_candidates.end(), retry_candidates.begin(), retry_candidates.end());
        unchanged_selected.insert(unchanged_selected.end(), retry_selected.begin(), retry_selected.end());
        candidates = std::move(unchanged_candidates);
        selected = std::move(unchanged_selected);
    }
    refine_boundaries_(selected, candidates);
    mark_used_relocations_(selected);
    mark_redundant_fragments_(selected);
    write_html_(candidates);
    build_sample_chains_(selected);
    print_summary_(selected.size());
    return selected.size();
}

void GfaGapfill::write_html_(const std::vector<Candidate>& candidates) const {
    if (params_.html_file.empty()) return;

    std::vector<GfaGapfillPlotter::Contig> contigs;
    contigs.reserve(fragments_.size());
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible) continue;
        GfaGapfillPlotter::Contig contig;
        contig.id = i;
        contig.component = fragment.component;
        contig.hap = fragment.hap;
        contig.bp = fragment.length;
        contig.name = paths_[fragment.path_id].name;
        contig.sample = fragment.sample;
        contig.first = vertex_name_(fragment.vertices.front());
        contig.last = vertex_name_(fragment.vertices.back());
        contig.anchors.reserve(fragment.vertices.size());
        for (uint32_t j = 0; j < fragment.vertices.size(); ++j) {
            const uint32_t segment = Vertex::get_segment_id(fragment.vertices[j]);
            contig.anchors.push_back({
                segment,
                nodes_[segment].length,
                fragment.path_bp[j] + nodes_[segment].length / 2
            });
        }
        contigs.push_back(std::move(contig));
    }

    std::vector<GfaGapfillPlotter::Connection> connections;
    connections.reserve(candidates.size());
    for (const Candidate& candidate : candidates) {
        GfaGapfillPlotter::Connection connection;
        connection.component = fragments_[candidate.left].component;
        connection.left = candidate.left;
        connection.right = candidate.right;
        connection.bridge = candidate.bridge;
        connection.phase_left = candidate.left_phase.similarity;
        connection.phase_right = candidate.right_phase.similarity;
        connection.phase_left_target_bp = candidate.left_phase.target_bp;
        connection.phase_left_bridge_bp = candidate.left_phase.bridge_bp;
        connection.phase_right_target_bp = candidate.right_phase.target_bp;
        connection.phase_right_bridge_bp = candidate.right_phase.bridge_bp;
        connection.boundary_left = candidate.left_boundary_similarity;
        connection.boundary_right = candidate.right_boundary_similarity;
        connection.alignment_left = candidate.left_alignment_similarity;
        connection.alignment_right = candidate.right_alignment_similarity;
        connection.probability = candidate.probability;
        connection.homolog_span = candidate.homolog_span;
        connection.sample_support = candidate.sample_support;
        connection.informative_samples = candidate.informative_samples;
        connection.spanning_samples = candidate.spanning_samples;
        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        const Fragment& bridge = fragments_[candidate.bridge];
        const bool final_boundary = candidate.status == Candidate::KEPT;
        connection.left_cut_bp = final_boundary ? candidate.left_cut_bp : left.path_bp[candidate.left_pos + 1];
        connection.right_cut_bp = final_boundary ? candidate.right_cut_bp : right.path_bp[candidate.right_pos];
        connection.bridge_left_bp = final_boundary ? candidate.bridge_left_cut_bp : bridge.path_bp[candidate.bridge_left + 1];
        connection.bridge_right_bp = final_boundary ? candidate.bridge_right_cut_bp : bridge.path_bp[candidate.bridge_right];
        switch (candidate.status) {
            case Candidate::KEPT: connection.status = "keep"; break;
            case Candidate::LOW_PROBABILITY: connection.status = "drop:low-probability"; break;
            case Candidate::USED_END: connection.status = "drop:used-end"; break;
            case Candidate::PHASE_CONFLICT: connection.status = "drop:phase-conflict"; break;
            case Candidate::COORDINATE_CONFLICT: connection.status = "drop:coordinate-conflict"; break;
            case Candidate::CYCLE: connection.status = "drop:cycle"; break;
            default: connection.status = "pending"; break;
        }
        connections.push_back(std::move(connection));
    }

    GfaGapfillPlotter::save(params_.html_file, contigs, connections);
}

std::string GfaGapfill::vertices_sequence_(const std::vector<uint32_t>& vertices) const {
    size_t reserve = 0;
    for (uint32_t vertex : vertices) {
        const std::string& sequence = nodes_[Vertex::get_segment_id(vertex)].sequence;
        if (sequence.empty() || sequence == "*") return {};
        reserve += sequence.size();
    }

    std::string sequence;
    sequence.reserve(reserve);
    for (size_t i = 0; i < vertices.size(); ++i) {
        const uint32_t overlap = i == 0 ? 0 : get_edge_ow(vertices[i - 1], vertices[i]);
        std::string part = get_oriented_sequence(Vertex(vertices[i]), overlap);
        const uint32_t length = nodes_[Vertex::get_segment_id(vertices[i])].length;
        if (part.size() != (length > overlap ? length - overlap : 0)) return {};
        for (char& base : part) {
            base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
        }
        sequence += part;
    }
    return sequence;
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

    log_stream() << "Aligning small primary components to larger haplotypes ...\n";
    mm_idxopt_t index_options;
    mm_mapopt_t default_map_options;
    mm_set_opt(nullptr, &index_options, &default_map_options);
    ChainRefs low_quality;
    if (mm_set_opt(mm2_.preset.c_str(), &index_options, &default_map_options) < 0) return small_chains;
    index_options.k = static_cast<short>(mm2_.k);
    index_options.w = static_cast<short>(mm2_.w);

    for (uint8_t hap = 0; hap < 2; ++hap) {
        std::vector<std::string> reference_sequences;
        std::vector<std::string> reference_names;
        reference_sequences.reserve(large[hap].size());
        reference_names.reserve(large[hap].size());
        for (size_t i = 0; i < large[hap].size(); ++i) {
            std::string sequence = vertices_sequence_(large[hap][i]->vertices);
            if (sequence.empty()) continue;
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

        std::vector<Query> queries;
        queries.reserve(small_chains[hap].size());
        for (const Chain* chain : small_chains[hap]) {
            queries.push_back({chain, vertices_sequence_(chain->vertices)});
        }
        if (queries.empty()) continue;
        if (sequences.empty()) {
            low_quality[hap] = small_chains[hap];
            continue;
        }

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

GfaGapfill::PrimarySelection GfaGapfill::select_primary_chains_() const {
    log_stream() << "Selecting primary haplotypes by component ...\n";
    using HapChains = std::array<std::vector<const Chain*>, 2>;
    std::map<uint32_t, std::map<std::string, HapChains>> components;
    std::unordered_map<const Chain*, uint64_t> chain_lengths;
    for (const auto& sample : sample_chains_) {
        for (uint8_t hap = 0; hap < 2; ++hap) {
            for (const Chain& chain : sample.second[hap]) {
                components[chain.component][sample.first][hap].push_back(&chain);
                chain_lengths.emplace(&chain, chain_length_(chain));
            }
        }
    }

    std::map<uint32_t, uint64_t> component_sizes;
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

    std::map<uint32_t, HapChains> primary_components;
    for (const auto& component : components) {
        const uint64_t component_bp = component_sizes.at(component.first);

        const HapChains* best = nullptr;
        std::string best_sample;
        std::tuple<bool, uint64_t, uint64_t, uint64_t, uint64_t> best_score;
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
            GfaGapfillDebugger::primary_candidate(
                component.first, sample.first, component_bp, hap_n50[0], hap_n50[1]
            );
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

    PrimarySelection selection;
    selection.low_quality = unmatched_small_primary_chains_(primary_components, chain_lengths);
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
    for (const auto& component : primary_components) {
        selection.primary[0].insert(selection.primary[0].end(), component.second[0].begin(), component.second[0].end());
        selection.primary[1].insert(selection.primary[1].end(), component.second[1].begin(), component.second[1].end());
    }
    return selection;
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
    const std::string suffix = ".gapfill";
    SAVE with_seq(prefix + suffix + ".gfa");
    SAVE no_seq(prefix + suffix + ".noseq.gfa");
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
    header += "H\tCO:Z:  SC:Z   target contigs in assembly order (+/- orientation)\n";
    header += "H\tCO:Z:  GS:Z   gap-source contigs in GF interval order (+/- orientation)\n";
    header += "H\tCO:Z:  MS:Z   source contig of an unused misassembly\n";
    header += "H\tCO:Z:  ME:Z   source end of an unused misassembly\n";
    header += "H\tCO:Z:  MR:Z   reason the misassembly was not used\n";
    with_seq.save(header);
    no_seq.save(header);

    std::unordered_set<std::string> output_names;
    for (const Chain* chain_ptr : chains) {
        const Chain& chain = *chain_ptr;
        if (!chain.gaps.empty() || chain.source_fragments.empty()) continue;
        const Fragment& source = fragments_[chain.source_fragments.front()];
        output_names.insert(output_contig_name(paths_[source.path_id].name, chain.sample));
    }

    size_t filled_number = 0;
    size_t primary_number = 0;
    for (const Chain* chain_ptr : chains) {
        const Chain& chain = *chain_ptr;
        if (chain.source_fragments.empty()) continue;
        const bool gap_filled = !chain.gaps.empty();
        const Fragment& first_source = fragments_[chain.source_fragments.front()];

        std::string name = output_contig_name(paths_[first_source.path_id].name, chain.sample);
        if (primary) {
            name = numbered_contig_name(hap, primary_kind, ++primary_number);
        } else if (gap_filled) {
            do {
                name = numbered_contig_name(hap, "gf", ++filled_number);
            } while (!output_names.insert(name).second);
        }

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

        uint64_t assembled_length = 0;
        bool sequence_available = true;
        std::string sequence;
        bool in_gap = false;
        uint64_t gap_begin = 0;
        std::vector<std::pair<uint64_t, uint64_t>> intervals;

        if (gap_filled) {
            sequence.reserve(chain_length_(chain));
            for (const Chain::Part& chain_part : chain.parts) {
                if (chain_part.filled && !in_gap) {
                    gap_begin = assembled_length;
                    in_gap = true;
                } else if (!chain_part.filled && in_gap) {
                    intervals.emplace_back(gap_begin, assembled_length);
                    in_gap = false;
                }
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
                assembled_length += chain_part.end - chain_part.begin;
            }
        } else {
            for (size_t i = 0; i < assembled_vertices->size(); ++i) {
                const uint32_t vertex = (*assembled_vertices)[i];
                const uint32_t sid = Vertex::get_segment_id(vertex);
                const uint32_t overlap = i == 0 ? 0 :
                    get_edge_ow((*assembled_vertices)[i - 1], vertex);
                const uint32_t added = nodes_[sid].length > overlap ?
                    nodes_[sid].length - overlap : 0;
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
        if (in_gap) intervals.emplace_back(gap_begin, assembled_length);

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
            for (size_t i = 0; i < chain.source_fragments.size(); ++i) {
                if (i) tags << ',';
                const Fragment& source = fragments_[chain.source_fragments[i]];
                tags << output_contig_name(paths_[source.path_id].name, chain.sample) << (source.reverse ? '-' : '+');
            }
            tags << "\tGS:Z:";
            for (size_t i = 0; i < chain.gaps.size(); ++i) {
                if (i) tags << ',';
                const Fragment& source = fragments_[chain.gaps[i].bridge];
                tags << paths_[source.path_id].name << (source.reverse ? '-' : '+');
            }
        }

        with_seq.save("S\t" + name + '\t');
        with_seq.save(sequence_available ? sequence : "*");
        with_seq.save(tags.str() + '\n');
        no_seq.save("S\t" + name + "\t*" + tags.str() + '\n');

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
                candidate.left_phase.similarity,
                candidate.right_phase.similarity,
                candidate.left_alignment_similarity,
                candidate.right_alignment_similarity,
                candidate.probability
            });
        }
        GfaGapfillDebugger::chain(
            name, hap, chain.component, chain.gaps.size() + 1,
            intervals.size(), assembled_length
        );
    }

    size_t unplaced_number = 0;
    size_t unplaced_count = 0;
    for (const RelocationRecord& record : relocations_) {
        if (primary || record.sample != label || record.hap != hap || record.used_for_gap ||
            record.unplaced_vertices.empty()) continue;

        std::string name;
        do {
            name = numbered_contig_name(hap, "ms", ++unplaced_number);
        } while (!output_names.insert(name).second);

        uint64_t length = 0;
        for (size_t i = 0; i < record.unplaced_vertices.size(); ++i) {
            const uint32_t vertex = record.unplaced_vertices[i];
            const uint32_t overlap = i == 0 ? 0 :
                get_edge_ow(record.unplaced_vertices[i - 1], vertex);
            const uint32_t node_length = nodes_[Vertex::get_segment_id(vertex)].length;
            length += node_length > overlap ? node_length - overlap : 0;
        }
        const std::string sequence = vertices_sequence_(record.unplaced_vertices);
        std::ostringstream tags;
        tags << "\tLN:i:" << length << "\tSN:Z:" << label
             << "\tCP:i:" << record.source_component
             << "\tMS:Z:" << record.source << "\tME:Z:" << record.side
             << "\tMR:Z:" << record.status;
        with_seq.save("S\t" + name + '\t' + (sequence.empty() ? "*" : sequence) + tags.str() + '\n');
        no_seq.save("S\t" + name + "\t*" + tags.str() + '\n');
        ++unplaced_count;
    }

    log_stream() << "  - " << label << ".hap" << static_cast<uint32_t>(hap) << ": " << chains.size() + unplaced_count << " contigs\n";
}

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
        "phase_left\tphase_right\talign_left\talign_right\tp\n"
    );
    for (const GapRecord& record : records_) {
        std::ostringstream line;
        line << record.sample << '\t' << static_cast<uint32_t>(record.hap) << '\t'
             << record.component << '\t'
             << record.contig << '\t' << record.left << '\t' << record.right << '\t'
             << record.bridge << '\t'
             << record.gap_beg << '-' << record.gap_end << '\t'
             << std::fixed << std::setprecision(2);
        if (record.phase_left < 0.0) line << ".\t";
        else line << record.phase_left << '\t';
        if (record.phase_right < 0.0) line << ".\t";
        else line << record.phase_right << '\t';
        if (record.alignment_left < 0.0) line << ".\t";
        else line << record.alignment_left << '\t';
        if (record.alignment_right < 0.0) line << ".\t";
        else line << record.alignment_right << '\t';
        line << record.probability << '\n';
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

void GfaGapfillDebugger::path(
    const std::string& name,
    const std::string& decision,
    const std::string& sample,
    uint8_t hap,
    uint32_t component,
    size_t nodes,
    uint64_t bp,
    const std::string& first,
    const std::string& last
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Path: " << name << '\n';
    debug_stream() << "  - Decision: " << decision << '\n';
    if (decision != "keep") return;
    debug_stream() << "  - Sample: " << sample << '\n';
    debug_stream() << "  - Component: " << component << " hap=" << static_cast<uint32_t>(hap) << '\n';
    debug_stream() << "  - Size: nodes=" << nodes << " bp=" << bp << '\n';
    debug_stream() << "  - Endpoints: " << first << " -> " << last << '\n';
}

void GfaGapfillDebugger::component(const std::string& sample, uint32_t component, size_t paths) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Component: " << component << '\n';
    debug_stream() << "  - Sample: " << sample << '\n';
    debug_stream() << "  - Contig paths: " << paths << '\n';
}

void GfaGapfillDebugger::overlap(
    const std::string& a,
    const std::string& b,
    uint64_t shared_bp,
    uint64_t overlap_bp,
    double similarity,
    double a_fraction,
    double b_fraction,
    const BoundarySupport& left,
    const BoundarySupport& right
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Overlap support:\n";
    debug_stream() << "  - Paths: " << a << " <-> " << b << '\n';
    debug_stream() << "  - Whole-node bp: shared=" << shared_bp << " overlap=" << overlap_bp << " similarity=" << std::fixed << std::setprecision(3) << similarity << '\n';
    debug_stream() << "  - Contig overlap fractions: " << std::fixed << std::setprecision(3) << a_fraction << ", " << b_fraction << '\n';
    std::ostringstream left_line;
    left_line << std::fixed << std::setprecision(3) << "  - Left overlap boundary: minimizer=";
    if (left.minimizer < 0.0) left_line << "NA";
    else left_line << left.minimizer;
    left_line << " phase=";
    if (left.phase < 0.0) left_line << "NA";
    else left_line << left.phase;
    left_line << " bubble_bp=" << left.first_bubble_bp << '/' << left.second_bubble_bp;
    debug_stream() << left_line.str() << '\n';

    std::ostringstream right_line;
    right_line << std::fixed << std::setprecision(3) << "  - Right overlap boundary: minimizer=";
    if (right.minimizer < 0.0) right_line << "NA";
    else right_line << right.minimizer;
    right_line << " phase=";
    if (right.phase < 0.0) right_line << "NA";
    else right_line << right.phase;
    right_line << " bubble_bp=" << right.first_bubble_bp << '/' << right.second_bubble_bp;
    debug_stream() << right_line.str() << '\n';
}

void GfaGapfillDebugger::search(
    const std::string& sample,
    size_t paths,
    uint64_t tested,
    uint64_t low_support,
    uint64_t wrong_group,
    uint64_t anchor_failed,
    size_t candidates
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Candidate search: " << sample << '\n';
    debug_stream() << "  - Contig paths: " << paths << '\n';
    debug_stream() << "  - Pairs tested: " << tested << '\n';
    debug_stream() << "  - Rejected: low_support=" << low_support << " wrong_group=" << wrong_group << " anchor_failed=" << anchor_failed << '\n';
    debug_stream() << "  - Candidates: " << candidates << '\n';
}

void GfaGapfillDebugger::candidate(
    const std::string& left,
    const std::string& right,
    const std::string& bridge,
    uint32_t component,
    double left_phase,
    double right_phase,
    uint64_t left_target_bp,
    uint64_t left_bridge_bp,
    uint64_t right_target_bp,
    uint64_t right_bridge_bp,
    double left_boundary,
    double right_boundary,
    bool phase_consistent,
    bool homolog_span,
    uint32_t sample_support,
    uint32_t informative_samples,
    uint32_t spanning_samples,
    double probability
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Candidate:\n";
    debug_stream() << "  - Component: " << component << '\n';
    debug_stream() << "  - Left: " << left << '\n';
    debug_stream() << "  - Right: " << right << '\n';
    debug_stream() << "  - Bridge: " << bridge << '\n';
    const BoundarySupport left_support{
        left_boundary, left_phase, left_target_bp, left_bridge_bp
    };
    const BoundarySupport right_support{
        right_boundary, right_phase, right_target_bp, right_bridge_bp
    };
    debug_stream() << boundary_table_header() << '\n';
    debug_stream() << boundary_table_row("Left", left_support, -1.0) << '\n';
    debug_stream() << boundary_table_row("Right", right_support, -1.0) << '\n';
    debug_stream() << "  - Best phase partner: " << (phase_consistent ? "yes" : "no") << '\n';
    debug_stream() << "  - Other hap spans gap: " << (homolog_span ? "yes" : "no") << '\n';
    debug_stream() << "  - Sample support: " << sample_support << '/' << informative_samples << " spanning=" << spanning_samples << " p=" << std::fixed << std::setprecision(4) << probability << '\n';
}

void GfaGapfillDebugger::refined_boundary(
    const std::string& left, uint64_t left_cut,
    const std::string& right, uint64_t right_cut, uint64_t right_length,
    const std::string& bridge, uint64_t bridge_begin, uint64_t bridge_end,
    const BoundarySupport& left_support,
    const BoundarySupport& right_support,
    double left_alignment,
    double right_alignment
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Final gap boundary:\n";
    debug_stream() << "  - Left: " << left << ":0-" << left_cut << '\n';
    debug_stream() << "  - Right: " << right << ':' << right_cut << '-' << right_length << '\n';
    debug_stream() << "  - Bridge: " << bridge << ':' << bridge_begin << '-' << bridge_end << '\n';
    debug_stream() << boundary_table_header() << '\n';
    debug_stream() << boundary_table_row("Left", left_support, left_alignment) << '\n';
    debug_stream() << boundary_table_row("Right", right_support, right_alignment) << '\n';
}

void GfaGapfillDebugger::boundary_alignment(
    const std::string& left,
    const std::string& right,
    const std::string& bridge,
    const Alignment& left_alignment,
    const Alignment& right_alignment
) {
    if (!DEBUG_ENABLED) return;

    debug_stream() << "Boundary alignment: " << left << " -> " << right << " via " << bridge << '\n';
    const Alignment* alignments[2] = {&left_alignment, &right_alignment};
    const char* labels[2] = {"Left", "Right"};
    const std::string* names[2] = {&left, &right};
    for (uint32_t side = 0; side < 2; ++side) {
        const Alignment& alignment = *alignments[side];
        debug_stream() << "  - " << labels[side] << ": " << *names[side] << ':'
                       << alignment.target_begin << '-' << alignment.target_end
                       << " vs " << bridge << ':' << alignment.bridge_begin
                       << '-' << alignment.bridge_end
                       << " query=" << alignment.target_begin + alignment.query_begin
                       << '-' << alignment.target_begin + alignment.query_end
                       << " ref=" << alignment.bridge_begin + alignment.reference_begin
                       << '-' << alignment.bridge_begin + alignment.reference_end
                       << (alignment.reverse ? '-' : '+')
                       << " hits=" << alignment.hits
                       << " mapq=" << static_cast<uint32_t>(alignment.mapq)
                       << " identity=" << std::fixed << std::setprecision(3) << alignment.identity
                       << " coverage=" << alignment.coverage
                       << " " << (alignment.accepted ? "keep" : "drop")
                       << std::endl;
    }
}

void GfaGapfillDebugger::selection(
    const std::string& left,
    const std::string& right,
    const std::string& bridge,
    const std::string& decision,
    double probability
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Selection: " << decision << '\n';
    debug_stream() << "  - Probability: " << std::fixed << std::setprecision(2) << probability << '\n';
    debug_stream() << "  - Left: " << left << '\n';
    debug_stream() << "  - Right: " << right << '\n';
    debug_stream() << "  - Bridge: " << bridge << '\n';
}

void GfaGapfillDebugger::assignment(
    const std::string& sample,
    uint32_t component,
    uint8_t hap,
    uint32_t hap1_votes,
    uint32_t hap2_votes,
    bool conflict
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Haplotype assignment:\n";
    debug_stream() << "  - Sample: " << sample << '\n';
    debug_stream() << "  - Component: " << component << '\n';
    debug_stream() << "  - Input votes: hap1=" << hap1_votes << " hap2=" << hap2_votes << '\n';
    debug_stream() << "  - Output haplotype: " << static_cast<uint32_t>(hap) << '\n';
    if (conflict) {
        debug_stream() << "  - Warning: more than two overlapping contigs\n";
    }
}

void GfaGapfillDebugger::primary_candidate(
    uint32_t component,
    const std::string& sample,
    uint64_t component_bp,
    uint64_t hap1_ng50,
    uint64_t hap2_ng50
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Primary candidate: component=" << component
                   << " sample=" << sample
                   << " component_bp=" << component_bp
                   << " hap1_NG50=" << hap1_ng50
                   << " hap2_NG50=" << hap2_ng50
                   << " mean_NG50=" << (hap1_ng50 + hap2_ng50) / 2
                   << '\n';
}

void GfaGapfillDebugger::chain(
    const std::string& name,
    uint8_t hap,
    uint32_t component,
    size_t contigs,
    size_t gaps,
    uint64_t bp
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Output contig: " << name << '\n';
    debug_stream() << "  - Component: " << component << " hap=" << static_cast<uint32_t>(hap) << '\n';
    debug_stream() << "  - Sources: contigs=" << contigs << " gaps=" << gaps << '\n';
    debug_stream() << "  - Length: " << bp << " bp\n";
}

namespace {

struct PlotRelation {
    uint32_t a{0}, b{0};
    uint64_t shared_bp{0};
    std::vector<std::pair<uint64_t, uint64_t>> anchors;
};

int64_t median(std::vector<int64_t> values) {
    if (values.empty()) return 0;
    const size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    return values[mid];
}

uint64_t deviation(const std::vector<int64_t>& values, int64_t center) {
    uint64_t sum = 0;
    for (int64_t value : values) sum += static_cast<uint64_t>(std::llabs(value - center));
    return sum;
}

void place_contigs(std::vector<GfaGapfillPlotter::Contig>& contigs) {
    using Contig = GfaGapfillPlotter::Contig;
    struct Occurrence { uint32_t contig; uint64_t offset; uint32_t length; };

    std::unordered_map<uint64_t, std::vector<Occurrence>> occurrences;
    for (uint32_t i = 0; i < contigs.size(); ++i) {
        for (const auto& anchor : contigs[i].anchors) {
            const uint64_t key = (static_cast<uint64_t>(contigs[i].component) << 32) | anchor.segment;
            occurrences[key].push_back({i, anchor.offset, anchor.length});
        }
    }

    std::unordered_map<uint64_t, PlotRelation> relation_map;
    for (const auto& item : occurrences) {
        const auto& found = item.second;
        for (size_t i = 0; i < found.size(); ++i) {
            for (size_t j = i + 1; j < found.size(); ++j) {
                if (found[i].contig == found[j].contig) continue;
                uint32_t a = found[i].contig, b = found[j].contig;
                uint64_t ao = found[i].offset, bo = found[j].offset;
                if (a > b) { std::swap(a, b); std::swap(ao, bo); }
                const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;
                PlotRelation& relation = relation_map[key];
                relation.a = a;
                relation.b = b;
                relation.shared_bp += std::min(found[i].length, found[j].length);
                relation.anchors.emplace_back(ao, bo);
            }
        }
    }

    std::vector<PlotRelation> relations;
    relations.reserve(relation_map.size());
    for (auto& item : relation_map) relations.push_back(std::move(item.second));
    std::vector<std::vector<uint32_t>> adjacency(contigs.size());
    for (uint32_t i = 0; i < relations.size(); ++i) {
        adjacency[relations[i].a].push_back(i);
        adjacency[relations[i].b].push_back(i);
    }

    std::map<uint32_t, std::vector<uint32_t>> components;
    for (uint32_t i = 0; i < contigs.size(); ++i) components[contigs[i].component].push_back(i);
    std::vector<uint8_t> placed(contigs.size(), 0);

    for (const auto& component : components) {
        int64_t cursor = 0;
        size_t remaining = component.second.size();
        while (remaining > 0) {
            uint32_t root = UINT32_MAX;
            for (uint32_t index : component.second) {
                if (!placed[index] && (root == UINT32_MAX || contigs[index].bp > contigs[root].bp)) root = index;
            }
            contigs[root].start = cursor;
            contigs[root].reverse = false;
            placed[root] = 1;
            --remaining;

            using QueueItem = std::tuple<uint64_t, uint32_t, uint32_t, uint32_t>;
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
                    const Contig& source = contigs[from];
                    const int64_t global = source.start + static_cast<int64_t>(
                        source.reverse ? source.bp - from_offset : from_offset
                    );
                    forward.push_back(global - static_cast<int64_t>(to_offset));
                    reverse.push_back(global - static_cast<int64_t>(contigs[to].bp - to_offset));
                }
                const int64_t forward_start = median(forward);
                const int64_t reverse_start = median(reverse);
                contigs[to].reverse = deviation(reverse, reverse_start) < deviation(forward, forward_start);
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
            cursor = end + std::max<int64_t>(1'000'000, end / 100);
        }

        int64_t begin = contigs[component.second.front()].start;
        for (uint32_t index : component.second) begin = std::min(begin, contigs[index].start);
        for (uint32_t index : component.second) contigs[index].start -= begin;
    }
}

}  // namespace

void GfaGapfillPlotter::save(
    const std::string& file,
    const std::vector<Contig>& contigs,
    const std::vector<Connection>& connections
) {
    log_stream() << "Writing interactive gap-fill candidate graph ...\n";
    std::vector<Contig> layout = contigs;
    place_contigs(layout);
    SAVE out(file);
    out.save(R"HTML(<!doctype html>
<html><head><meta charset="utf-8"><title>liftasm gap-fill candidates</title>
<style>
*{box-sizing:border-box}body{margin:0;font:14px system-ui,sans-serif;color:#182338;background:#f6f8fc}
header{height:54px;padding:10px 18px;background:linear-gradient(105deg,#122039,#29466d);color:white;font-size:18px;font-weight:650;letter-spacing:.2px;box-shadow:0 2px 8px #10192a33}
header small{margin-left:10px;color:#c9d7ea;font-size:12px;font-weight:450}
main{display:grid;grid-template-columns:400px 1fr;height:calc(100vh - 54px)}
aside{display:flex;min-height:0;flex-direction:column;border-right:1px solid #d8dee9;background:#fff;overflow:hidden;box-shadow:2px 0 8px #1720330a;z-index:2}.controls{padding:14px 16px 7px;border-bottom:1px solid #e4e9f0;background:#fbfcfe}
label{display:block;margin:0 0 5px;color:#536174;font-size:12px;font-weight:700;text-transform:uppercase}
select,input[type=search]{width:100%;padding:7px;margin-bottom:12px;border:1px solid #b8c1cf;border-radius:5px;background:white}input[type=range]{width:100%;margin:2px 0 14px}
#selection{padding:12px 16px;border-bottom:1px solid #e1e6ee;background:#fff;flex:none}#selection.focused{background:linear-gradient(135deg,#fff8f4,#fff)}#selection.focused label{color:#c94622}#detail{max-height:34vh;overflow:auto;white-space:pre-wrap;padding:11px 12px;border:1px solid #dfe5ed;border-radius:8px;background:#f4f6f9;line-height:1.48;overflow-wrap:anywhere}#selection.focused #detail{border-color:#f0a48e;background:#fff;box-shadow:0 4px 14px #c9462215}
.coordinates{margin-top:9px}.evidence{width:100%;margin:10px 0 8px;border-collapse:collapse;font-size:11px;white-space:nowrap}.evidence th,.evidence td{padding:5px 4px;border-bottom:1px solid #e7eaf0;text-align:right}.evidence th:first-child{text-align:left}.evidence thead th{color:#607087;font-weight:700}.support{color:#536174;font-size:12px}
#candidatePanel{min-height:0;flex:1;overflow:auto;padding:11px 16px 14px}#summary{margin:0 0 10px;color:#536174}.tools{display:flex;gap:7px;margin:-4px 0 7px}.tool{padding:6px 11px;border:1px solid #bcc6d4;border-radius:6px;background:#fff;color:#34435a;cursor:pointer}.tool:hover{background:#edf4ff}
.candidate{width:100%;padding:10px 11px;margin:0 0 8px;text-align:left;border:1px solid #dde3ec;border-left:4px solid #aab4c3;border-radius:8px;background:white;cursor:pointer;transition:background .14s,border-color .14s,box-shadow .14s,transform .14s}.candidate.keep{border-left-color:#159447}.candidate:hover{background:#edf4ff;transform:translateY(-1px);box-shadow:0 4px 10px #1b35551a}.candidate.active{border-color:#e4572e;border-left:7px solid #e4572e;background:linear-gradient(110deg,#fff0e9,#fff);box-shadow:0 0 0 2px #e4572e26,0 7px 18px #b83f1e24;transform:none}.candidate.active .route{color:#b53618}.route{display:block;font-weight:700;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.via,.metrics{display:block;margin-top:3px;color:#657286;font-size:11px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
#view{overflow:auto;position:relative;background:#fff;scroll-behavior:smooth}svg{min-height:100%;background:#fff}.lane{fill:#f8faff}.lane.alt{fill:#fff}.row{stroke:#e2e7ef;stroke-width:1}.axis{font-size:10px;fill:#77859a}.sample{font-size:12px;font-weight:700;fill:#34435a}.edge{fill:none;stroke:#9aa5b5;stroke-width:1.5;opacity:.20;cursor:pointer;transition:opacity .18s,stroke-width .18s}.edge.keep{stroke:#159447;stroke-width:2.5;opacity:.82}.edge.related{opacity:.42}.edge.active{stroke:#e4572e;stroke-width:4;opacity:1}.edge.dim{opacity:.025}.bridge-span{stroke-width:7;stroke-linecap:round;opacity:.32}.bridge-span.keep{opacity:.9}.anchor{stroke:white;stroke-width:1.5;fill:#7d8999;opacity:.7}.anchor.keep{fill:#159447;opacity:1}.anchor.active{fill:#e4572e}
.node{transition:opacity .18s}.node rect{fill:#edf2f8;fill-opacity:.82;stroke:#8492a6;stroke-width:1}.node.h1 rect{fill:#d9eaff;stroke:#3274c5}.node.h2 rect{fill:#ffe2e5;stroke:#d25763}.node.active rect{fill-opacity:1;stroke:#e4572e;stroke-width:3;filter:drop-shadow(0 3px 5px #e4572e55)}.node text{font-size:10px;pointer-events:none;fill:#1d2a3d}
.legend{display:flex;gap:12px;margin:2px 0 12px;color:#536174}.dot:before{content:'';display:inline-block;width:12px;height:3px;margin-right:5px;vertical-align:middle;background:#9aa5b5}.dot.keep:before{background:#159447}
#tooltip{position:fixed;display:none;max-width:420px;padding:7px 9px;border-radius:5px;background:#17243a;color:white;font-size:12px;line-height:1.4;white-space:pre-line;pointer-events:none;z-index:5;box-shadow:0 3px 12px #10182855}
</style></head><body><header>liftasm gap-fill <small>contig → bridge interval → contig</small></header><main><aside><div class="controls">
<label for="component">Component</label><select id="component"></select>
<label for="status">Connections</label><select id="status"><option value="all">All candidates</option><option value="keep">Kept only</option><option value="drop">Rejected only</option></select>
<label for="query">Find candidate</label><input id="query" type="search" placeholder="contig or sample">
<label for="zoom">Horizontal zoom: <span id="zoomValue">1x</span></label><input id="zoom" type="range" min="1" max="20" value="1" step="1">
<div class="tools"><button class="tool" id="fit" type="button">Fit view</button><button class="tool" id="clear" type="button">Clear selection</button></div>
<div class="legend"><span class="dot keep">kept</span><span class="dot">rejected</span></div></div>
<section id="selection"><label>Selected candidate</label><div id="detail">Click a connection or candidate.</div></section>
<section id="candidatePanel"><div id="summary"></div><div id="list"></div></section>
</aside><div id="view"><svg id="graph"></svg></div></main><div id="tooltip"></div><script>
const N=)HTML");

    out.save("[");
    for (size_t i = 0; i < layout.size(); ++i) {
        const Contig& n = layout[i];
        if (i) out.save(",");
        std::ostringstream line;
        line << "{i:" << n.id << ",c:" << n.component << ",h:" << static_cast<uint32_t>(n.hap)
             << ",b:" << n.bp << ",x:" << n.start << ",rv:" << (n.reverse ? 1 : 0)
             << ",n:" << json_string(n.name)
             << ",s:" << json_string(n.sample) << ",a:" << json_string(n.first)
             << ",z:" << json_string(n.last) << "}";
        out.save(line.str());
    }
    out.save("];\nconst E=[");
    for (size_t i = 0; i < connections.size(); ++i) {
        const Connection& e = connections[i];
        if (i) out.save(",");
        std::ostringstream line;
        line << std::setprecision(8)
             << "{i:" << i << ",c:" << e.component << ",l:" << e.left
             << ",r:" << e.right << ",b:" << e.bridge
             << ",pl:" << e.phase_left << ",prh:" << e.phase_right
             << ",plt:" << e.phase_left_target_bp << ",plb:" << e.phase_left_bridge_bp
             << ",prt:" << e.phase_right_target_bp << ",prb:" << e.phase_right_bridge_bp
             << ",ml:" << e.boundary_left << ",mr:" << e.boundary_right
             << ",al:" << e.alignment_left << ",ar:" << e.alignment_right
             << ",pr:" << e.probability << ",ss:" << e.sample_support
             << ",ns:" << e.informative_samples
             << ",sp:" << e.spanning_samples
             << ",hs:" << (e.homolog_span ? 1 : 0)
             << ",lp:" << e.left_cut_bp << ",rp:" << e.right_cut_bp
             << ",bl:" << e.bridge_left_bp << ",br:" << e.bridge_right_bp
             << ",st:" << json_string(e.status) << "}";
        out.save(line.str());
    }
    out.save(R"HTML(];
const $=function(id){return document.getElementById(id)},svg=$('graph'),view=$('view'),comp=$('component'),filter=$('status'),query=$('query'),zoom=$('zoom'),zoomValue=$('zoomValue'),list=$('list'),detail=$('detail'),selection=$('selection'),summary=$('summary'),tooltip=$('tooltip');
const ns='http://www.w3.org/2000/svg',byId=new Map(N.map(function(n){return[n.i,n]}));
function S(tag,attrs,text){const x=document.createElementNS(ns,tag);Object.keys(attrs||{}).forEach(function(k){x.setAttribute(k,attrs[k])});if(text!==undefined)x.textContent=text;return x}
function visible(e){const status=filter.value==='all'||(filter.value==='keep'&&e.st==='keep')||(filter.value==='drop'&&e.st!=='keep');if(!status)return false;const q=query.value.trim().toLowerCase();if(!q)return true;const a=byId.get(e.l),b=byId.get(e.r),c=byId.get(e.b);return[a.n,b.n,c.n,a.s,b.s,c.s,e.st].join(' ').toLowerCase().includes(q)}
function fmt(n){return n.toLocaleString()}
function val(n){return n<0?'NA':n.toFixed(2)}
function esc(s){return s.replace(/[&<>"']/g,function(c){return{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]})}
function phase(e){const values=[e.pl,e.prh].filter(function(x){return x>=0});return values.length?values.reduce(function(a,b){return a+b},0)/values.length:-1}
function tipText(e){return byId.get(e.l).n+' → '+byId.get(e.r).n+'\nvia '+byId.get(e.b).n+'\nphase '+val(e.pl)+' / '+val(e.prh)+' · support '+e.ss+'/'+e.ns+' · p '+e.pr.toFixed(2)}
function bindTip(mark,e){mark.onpointermove=function(event){tooltip.style.display='block';tooltip.style.left=Math.min(innerWidth-430,event.clientX+14)+'px';tooltip.style.top=Math.min(innerHeight-80,event.clientY+14)+'px';tooltip.textContent=tipText(e)};mark.onpointerleave=function(){tooltip.style.display='none'}}
function describe(e){const l=byId.get(e.l),r=byId.get(e.r),b=byId.get(e.b);return '<div><b>Status:</b> '+esc(e.st)+' &nbsp; <b>Component:</b> '+e.c+'</div><div class="coordinates"><b>Left:</b> '+esc(l.n)+':0-'+fmt(e.lp)+'<br><b>Right:</b> '+esc(r.n)+':'+fmt(e.rp)+'-'+fmt(r.b)+'<br><b>Bridge:</b> '+esc(b.n)+':'+fmt(e.bl)+'-'+fmt(e.br)+'</div><table class="evidence"><thead><tr><th>Boundary</th><th>Phase</th><th>Bubble bp</th><th>Minimizer</th><th>Alignment</th></tr></thead><tbody><tr><th>Left</th><td>'+val(e.pl)+'</td><td>'+fmt(e.plt)+'/'+fmt(e.plb)+'</td><td>'+val(e.ml)+'</td><td>'+val(e.al)+'</td></tr><tr><th>Right</th><td>'+val(e.prh)+'</td><td>'+fmt(e.prt)+'/'+fmt(e.prb)+'</td><td>'+val(e.mr)+'</td><td>'+val(e.ar)+'</td></tr></tbody></table><div class="support">Other hap spans gap: '+(e.hs?'yes':'no')+'<br>Sample support: '+e.ss+' / '+e.ns+'<br>Spanning samples: '+e.sp+'<br>Posterior p: '+e.pr.toFixed(2)+'</div>'}
let selectedEdge=null;
function centerEdge(id){const marks=svg.querySelectorAll('[data-e="'+id+'"]');let left=Infinity,right=-Infinity,top=Infinity,bottom=-Infinity;marks.forEach(function(mark){const box=mark.getBBox();left=Math.min(left,box.x);right=Math.max(right,box.x+box.width);top=Math.min(top,box.y);bottom=Math.max(bottom,box.y+box.height)});if(left!==Infinity)view.scrollTo({left:Math.max(0,(left+right-view.clientWidth)/2),top:Math.max(0,(top+bottom-view.clientHeight)/2),behavior:'smooth'})}
function selectEdge(id){selectedEdge=selectedEdge===id?null:id;render();if(selectedEdge!==null)requestAnimationFrame(function(){centerEdge(id)})}
function anchorPosition(n,offset){return n.x+(n.rv?n.b-offset:offset)}
function anchorX(n,offset,scale,left,origin){return left+(anchorPosition(n,offset)-origin)*scale}
function render(){
  const cid=Number(comp.value),allNodes=N.filter(function(n){return n.c===cid}),allEdges=E.filter(function(e){return e.c===cid&&visible(e)}),focus=selectedEdge===null?null:E[selectedEdge];
  if(focus&&(focus.c!==cid||!allEdges.some(function(e){return e.i===focus.i}))){selectedEdge=null;return render()}
  const focused=focus!==null,nodeIds=focused?new Set([focus.l,focus.r,focus.b]):null;
  const nodes=focused?allNodes.filter(function(n){return nodeIds.has(n.i)}):allNodes;
  const edges=focused?[focus]:allEdges;
  const rows=[...new Set(nodes.map(function(n){return n.s+'\t'+n.h}))].sort(function(a,b){const x=a.split('\t'),y=b.split('\t');return x[0].localeCompare(y[0])||Number(x[1])-Number(y[1])});
  let origin=0,end=Math.max(1,...nodes.map(function(n){return n.x+n.b}));
  if(focused){const positions=[anchorPosition(byId.get(focus.l),focus.lp),anchorPosition(byId.get(focus.b),focus.bl),anchorPosition(byId.get(focus.b),focus.br),anchorPosition(byId.get(focus.r),focus.rp)],range=Math.max(1,Math.max(...positions)-Math.min(...positions)),padding=Math.max(500000,range*.16);origin=Math.max(0,Math.min(...positions)-padding);end=Math.max(...positions)+padding}
  const left=170,rowHeight=focused?86:58,top=58,span=Math.max(1,end-origin),factor=Number(zoom.value);
  const plotWidth=Math.max(1000,$('view').clientWidth-left-40)*factor,scale=plotWidth/span,width=left+plotWidth+35,height=Math.max(300,top+rows.length*rowHeight+45),pos=new Map();
  svg.replaceChildren();svg.setAttribute('viewBox','0 0 '+width+' '+height);svg.setAttribute('width',width);svg.setAttribute('height',height);zoomValue.textContent=factor+'x';
  rows.forEach(function(row,index){const parts=row.split('\t'),y=top+index*rowHeight;svg.appendChild(S('rect',{x:0,y:y-rowHeight/2,width:width,height:rowHeight,class:'lane '+(index%2?'alt':'')}));svg.appendChild(S('line',{x1:left,y1:y,x2:width-20,y2:y,class:'row'}));svg.appendChild(S('text',{x:14,y:y+4,class:'sample'},parts[0]+'  hap'+parts[1]));pos.set(row,y)});
  for(let tick=0;tick<=5;++tick){const bp=origin+span*tick/5,x=left+plotWidth*tick/5;svg.appendChild(S('line',{x1:x,y1:25,x2:x,y2:height-20,class:'row'}));svg.appendChild(S('text',{x:x,y:18,'text-anchor':'middle',class:'axis'},(bp/1000000).toFixed(span>=10000000?0:1)+' Mb'))}
  edges.forEach(function(e){const a=byId.get(e.l),m=byId.get(e.b),z=byId.get(e.r),active=focused?' active':'';if(!a||!m||!z)return;const points=[{x:anchorX(a,e.lp,scale,left,origin),y:pos.get(a.s+'\t'+a.h)},{x:anchorX(m,e.bl,scale,left,origin),y:pos.get(m.s+'\t'+m.h)},{x:anchorX(m,e.br,scale,left,origin),y:pos.get(m.s+'\t'+m.h)},{x:anchorX(z,e.rp,scale,left,origin),y:pos.get(z.s+'\t'+z.h)}];[[points[0],points[1]],[points[2],points[3]]].forEach(function(pair){const middle=(pair[0].y+pair[1].y)/2,path=S('path',{d:'M'+pair[0].x+','+pair[0].y+' C'+pair[0].x+','+middle+' '+pair[1].x+','+middle+' '+pair[1].x+','+pair[1].y,class:'edge '+(e.st==='keep'?'keep':'')+active,'data-e':e.i});path.onclick=function(){selectEdge(e.i)};bindTip(path,e);svg.appendChild(path)})});
  nodes.sort(function(a,b){return a.s.localeCompare(b.s)||a.h-b.h||a.x-b.x}).forEach(function(n){const y=pos.get(n.s+'\t'+n.h),nodeBegin=Math.max(n.x,origin),nodeEnd=Math.min(n.x+n.b,origin+span);if(nodeEnd<=nodeBegin)return;const x=left+(nodeBegin-origin)*scale,w=Math.max(1,(nodeEnd-nodeBegin)*scale),g=S('g',{class:'node h'+n.h+(focused?' active':''),'data-n':n.i});g.appendChild(S('rect',{x:x,y:y-(focused?16:11),width:w,height:focused?32:22,rx:3}));if(w>55)g.appendChild(S('text',{x:x+w/2,y:y+4,'text-anchor':'middle'},n.n.length>34?n.n.slice(0,31)+'...':n.n));g.appendChild(S('title',{},n.n+'\ncomponent: '+(n.x/1000000).toFixed(3)+'-'+((n.x+n.b)/1000000).toFixed(3)+' Mb\nlength: '+fmt(n.b)+' bp\norientation: '+(n.rv?'reverse':'forward')+'\n'+n.a+' -> '+n.z));g.onclick=function(event){event.stopPropagation();const edge=edges.find(function(e){return e.l===n.i||e.r===n.i||e.b===n.i});if(edge)selectEdge(edge.i)};svg.appendChild(g)});
  edges.forEach(function(e){const m=byId.get(e.b),y=pos.get(m.s+'\t'+m.h),x1=anchorX(m,e.bl,scale,left,origin),x2=anchorX(m,e.br,scale,left,origin),kind=(e.st==='keep'?'keep ':'')+(focused?'active':'');const bridgeMark=S('line',{x1:x1,y1:y,x2:x2,y2:y,class:'edge bridge-span '+kind,'data-e':e.i});bridgeMark.onclick=function(){selectEdge(e.i)};bindTip(bridgeMark,e);svg.appendChild(bridgeMark);[x1,x2].forEach(function(x){const dot=S('circle',{cx:x,cy:y,r:focused?6:4,class:'edge anchor '+kind,'data-e':e.i});dot.onclick=function(){selectEdge(e.i)};bindTip(dot,e);svg.appendChild(dot)})});
  list.replaceChildren();edges.slice().sort(function(a,b){return(a.st==='keep'?0:1)-(b.st==='keep'?0:1)||byId.get(a.l).x-byId.get(b.l).x}).forEach(function(e){const button=document.createElement('button'),route=document.createElement('span'),via=document.createElement('span'),metrics=document.createElement('span'),p=phase(e);button.className='candidate '+(e.st==='keep'?'keep ':'')+(focused?'active':'');button.dataset.e=e.i;route.className='route';route.textContent=byId.get(e.l).n+' → '+byId.get(e.r).n;via.className='via';via.textContent='via '+byId.get(e.b).n;metrics.className='metrics';metrics.textContent=e.st+' · phase '+val(p)+' · support '+e.ss+'/'+e.ns+' · p '+e.pr.toFixed(2);button.append(route,via,metrics);button.onclick=function(){selectEdge(e.i)};list.appendChild(button)});
  selection.classList.toggle('focused',focused);if(focused)detail.innerHTML=describe(focus);else detail.textContent='Click a connection or candidate.';const kept=edges.filter(function(e){return e.st==='keep'}).length;summary.textContent=focused?'Focused connection · '+nodes.length+' contigs':nodes.length+' contigs · '+kept+' kept · '+(edges.length-kept)+' rejected';view.scrollTo({left:0,top:0})
}
[...new Set(N.map(function(n){return n.c}))].sort(function(a,b){return a-b}).forEach(function(c){const o=document.createElement('option');o.value=c;o.textContent='component '+c;comp.appendChild(o)});comp.onchange=function(){selectedEdge=null;render()};filter.onchange=function(){selectedEdge=null;render()};query.oninput=function(){selectedEdge=null;render()};$('fit').onclick=function(){zoom.value=1;render()};$('clear').onclick=function(){selectedEdge=null;render()};let zoomFrame=0;zoom.oninput=function(){cancelAnimationFrame(zoomFrame);zoomFrame=requestAnimationFrame(render)};svg.onclick=function(event){if(!event.target.closest('[data-e],[data-n]')){selectedEdge=null;render()}};document.onkeydown=function(event){if(event.key==='Escape')$('clear').click()};render();
</script></body></html>)HTML");
    log_stream() << "  - HTML: " << file << "\n\n";
}
