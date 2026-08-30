#include "gfa_gapfill.hpp"
#include "gfa_bubble.hpp"
#include "gfa_gapfill_logger.hpp"
#include "gfa_walker.hpp"

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

bool is_backbone_name(const std::string& name) {
    constexpr char prefix[] = "component";
    constexpr size_t prefix_length = sizeof(prefix) - 1;
    return name.size() > prefix_length &&
        name.compare(0, prefix_length, prefix) == 0 &&
        std::all_of(name.begin() + prefix_length, name.end(), [](unsigned char c) {
            return std::isdigit(c);
        });
}

// Shared nodes of a contig pair
struct PlotRelation {
    uint32_t a{0}, b{0};   // contig indexes
    uint64_t shared_bp{0}; // total shared base pairs
    std::vector<std::pair<uint64_t, uint64_t>> anchors;  // shared graph-node offsets on contig a and b
};

int64_t layout_median(std::vector<int64_t>& values) {
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
GfaGapfill::GfaGapfill(GapfillOpts options) : params_(std::move(options)) {
    // Gapfill uses node boundaries from a deoverlapped graph.
    set_forbid_overlap(true);
    full_phase_samples_.insert(
        params_.full_phase_samples.begin(), params_.full_phase_samples.end()
    );
    mm2_.k = params_.alignment_options.chain.k;
    mm2_.w = params_.alignment_options.chain.w;
    mm2_.best_n = params_.alignment_options.anchor.max_kept;
    mm2_.zdrop = params_.alignment_options.extend.dyn_zdrop;
    mm2_.minimizer_freq = params_.alignment_options.chain.mid_occ_frac;
    mm2_.max_occ = params_.alignment_options.chain.max_mid_occ;
    mm2_.preset = params_.mm2_preset;
    mm2_.min_match = params_.min_match;
    mm2_.min_ali_ratio = params_.min_ali_ratio;
    mm2_.min_mapq = params_.min_mapq;
    params_.threads = std::max<uint32_t>(1, params_.threads);
}


// ================================================= Clear state =================================================
void GfaGapfill::clear_state_() {
    fragments_.clear();
    backbone_bp_.clear();
    source_ids_.clear();
    source_scope_.clear();
    source_transitions_.clear();
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
    log_stream() << "Gap-filling sample contigs with source-aware graph walks ...\n\n";
    clear_state_();
    index_bubble_nodes_(bubbles);
    build_fragments_();
    build_graph_order_();
    index_source_transitions_();
    std::vector<Candidate> candidates = build_candidates_();

    std::vector<Candidate> selected = select_and_walk_(candidates);

    mark_used_relocations_(selected);
    build_sample_chains_(selected);
    mark_redundant_fragments_(selected);

    write_html_(candidates);
    print_summary_(selected.size());

    return selected.size();
}

std::vector<GfaGapfill::MisassemblyContig> GfaGapfill::prepare_misassemblies() {
    log_stream() << "Detecting terminal misassemblies before gap filling ...\n\n";
    clear_state_();
    bubble_nodes_.assign(getNumNodes(), 0);
    build_fragments_();
    return GfaGapfillMisassembly(*this).run();
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

    // ------------------------------------------------ Index haplotype SN labels ------------------------------------------------
    {
        const std::vector<std::string>& source_names = getSampleNames();
        for (uint32_t id = 0; id < source_names.size(); ++id) {
            source_ids_.emplace(source_names[id], id);
            if (source_names[id].find(".hap1") != std::string::npos || source_names[id].find(".hap2") != std::string::npos) {
                source_scope_.insert(id);
            }
        }
    }

    // ------------------------------------------------ Read path names ------------------------------------------------
    for (size_t path_id = 0; path_id < paths_.size(); ++path_id) {
        GfaPath& path = paths_[path_id];
        Fragment fragment;
        if (!parse_path_name_(path.name, fragment.sample, fragment.hap)) {
            GfaGapfillDebugger::path(path.name, "drop:unrecognized-sample");
            continue;
        }

        // P names are authoritative local haplotype labels. Register them even
        // when the input S lines do not carry matching SN tags.
        const std::string source_name = fragment.sample + ".hap" + std::to_string(fragment.hap);
        if (!source_ids_.count(source_name)) {
            const uint32_t source_id = static_cast<uint32_t>(sample_id_to_name_.size());
            source_ids_.emplace(source_name, source_id);
            source_scope_.insert(source_id);
            sample_name_to_id_.emplace(source_name, source_id);
            sample_id_to_name_.push_back(source_name);
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
        fragment.eligible = fragment.length >= params_.min_contig;

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

    // ------------------------------------------------ Validate full-phase sample names ------------------------------------------------
    if (!full_phase_samples_.empty()) {
        std::unordered_set<std::string> samples;
        for (const Fragment& fragment : fragments_) samples.insert(fragment.sample);
        for (const std::string& sample : full_phase_samples_) {
            if (!samples.count(sample)) {
                warning_stream() << "  ! Unknown --full_phase sample: " << sample << '\n';
            }
        }
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
        contig.bp = fragment.length;
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
    const size_t backbone_placed = place_contigs_on_backbones_(layout);
    if (backbone_placed < layout.size()) place_contigs_(layout);

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

void GfaGapfill::index_source_transitions_() {
    log_stream() << "Indexing haplotype-supported graph transitions ...\n";

    // ------------------------------------------------ Collect oriented P-path transitions ------------------------------------------------
    uint64_t path_transitions = 0;
    for (const Fragment& fragment : fragments_) {
        if (fragment.hap < 1 || fragment.hap > 2 || fragment.vertices.size() < 2) continue;
        const auto source = source_ids_.find(fragment.sample + ".hap" + std::to_string(fragment.hap));
        if (source == source_ids_.end()) continue;

        path_transitions += fragment.vertices.size() - 1;
        for (size_t i = 1; i < fragment.vertices.size(); ++i) {
            const uint32_t from = fragment.vertices[i - 1];
            const uint32_t to = fragment.vertices[i];
            source_transitions_[(static_cast<uint64_t>(from) << 32) | to].push_back(source->second);
            source_transitions_[(static_cast<uint64_t>(to ^ 1u) << 32) | (from ^ 1u)].push_back(source->second);
        }
    }

    // ------------------------------------------------ Deduplicate source support ------------------------------------------------
    for (auto& transition : source_transitions_) {
        std::vector<uint32_t>& sources = transition.second;
        std::sort(sources.begin(), sources.end());
        sources.erase(std::unique(sources.begin(), sources.end()), sources.end());
    }

    log_stream() << "  - P-path transitions: " << path_transitions << '\n';
    log_stream() << "  - Oriented transitions indexed: " << source_transitions_.size() << "\n\n";
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
    fragment.path_bp.assign(fragment.vertices.size() + 1, 0);
    fragment.bubble_bp.assign(fragment.vertices.size() + 1, 0);
    for (uint32_t i = 0; i < fragment.vertices.size(); ++i) {
        const uint32_t vertex = fragment.vertices[i];
        const uint32_t segment = Vertex::get_segment_id(vertex);
        fragment.vertex_index.emplace_back(vertex, i);
        fragment.path_bp[i + 1] = fragment.path_bp[i] + nodes_[segment].length;
        fragment.bubble_bp[i + 1] = fragment.bubble_bp[i] + (segment < bubble_nodes_.size() && bubble_nodes_[segment] ? nodes_[segment].length : 0);
        if (segment < bubble_nodes_.size() && bubble_nodes_[segment]) {
            fragment.bubble_positions.push_back(i);
        }
    }
    std::sort(fragment.vertex_index.begin(), fragment.vertex_index.end());
    fragment.length = fragment.path_bp.back();
}

size_t GfaGapfill::place_contigs_on_backbones_(std::vector<LayoutContig>& contigs) {
    struct BackbonePath {
        size_t path_id{0};
        uint64_t bp{0};
    };
    struct BackboneNode {
        uint64_t offset{0};
        uint32_t component{UINT32_MAX};
        bool unique{false};
    };

    // ------------------------------------------------ Select one backbone per graph component ------------------------------------------------
    std::unordered_map<uint32_t, BackbonePath> backbones;
    for (size_t path_id = 0; path_id < paths_.size(); ++path_id) {
        const GfaPath& path = paths_[path_id];
        if (!is_backbone_name(path.name) || path.segments.empty()) continue;

        uint32_t component = UINT32_MAX;
        uint64_t bp = 0;
        bool valid = true;
        for (const PathSegment& segment : path.segments) {
            const uint32_t sid = static_cast<uint32_t>(segment.node_id);
            if (sid >= nodes_.size() || nodes_[sid].deleted) {
                valid = false;
                break;
            }
            const uint32_t node_component = connectivity_index_[sid];
            if (component == UINT32_MAX) component = node_component;
            else if (component != node_component) {
                valid = false;
                break;
            }
            bp += nodes_[sid].length;
        }
        if (!valid) continue;

        auto found = backbones.find(component);
        if (found == backbones.end() || bp > found->second.bp) {
            backbones[component] = {path_id, bp};
        }
    }
    if (backbones.empty()) return 0;

    // ------------------------------------------------ Index unambiguous backbone node coordinates ------------------------------------------------
    std::vector<BackboneNode> backbone_nodes(nodes_.size());
    for (const auto& item : backbones) {
        const uint32_t component = item.first;
        const BackbonePath& backbone = item.second;
        const GfaPath& path = paths_[backbone.path_id];
        uint64_t offset = 0;
        for (const PathSegment& segment : path.segments) {
            const uint32_t sid = static_cast<uint32_t>(segment.node_id);
            BackboneNode& node = backbone_nodes[sid];
            if (node.component == UINT32_MAX) {
                node.offset = offset + nodes_[sid].length / 2;
                node.component = component;
                node.unique = true;
            } else {
                node.unique = false;
            }
            offset += nodes_[sid].length;
        }
        backbone_bp_[component] = backbone.bp;
    }

    // ------------------------------------------------ Project each contig directly onto its backbone ------------------------------------------------
    size_t placed = 0;
    for (LayoutContig& contig : contigs) {
        if (!backbones.count(contig.component)) continue;

        std::vector<int64_t> forward, reverse;
        forward.reserve(contig.anchors.size());
        reverse.reserve(contig.anchors.size());
        for (const LayoutAnchor& anchor : contig.anchors) {
            const BackboneNode& node = backbone_nodes[anchor.segment];
            if (!node.unique || node.component != contig.component) continue;
            const int64_t backbone_offset = static_cast<int64_t>(node.offset);
            forward.push_back(backbone_offset - static_cast<int64_t>(anchor.offset));
            reverse.push_back(backbone_offset - static_cast<int64_t>(contig.bp - anchor.offset));
        }
        if (forward.empty()) continue;

        const int64_t forward_start = layout_median(forward);
        const int64_t reverse_start = layout_median(reverse);
        contig.reverse = layout_deviation(reverse, reverse_start) <
            layout_deviation(forward, forward_start);
        contig.start = contig.reverse ? reverse_start : forward_start;
        contig.placed = true;
        ++placed;
    }

    log_stream() << "  - Backbone paths indexed: " << backbones.size() << '\n';
    log_stream() << "  - Contigs placed on backbone: " << placed << '\n';
    log_stream() << "  - Contigs requiring graph fallback: " << contigs.size() - placed << '\n';
    return placed;
}

void GfaGapfill::place_contigs_(std::vector<LayoutContig>& contigs) {
    // ------------------------------------------------ Index shared nodes ------------------------------------------------
    struct Occurrence {
        uint32_t contig;  // contig index
        uint64_t offset;  // node offset on contig
        uint32_t length;  // node length
    };

    // Only nodes touched by an unplaced contig can contribute to fallback.
    std::unordered_set<uint64_t> fallback_nodes;
    for (const LayoutContig& contig : contigs) {
        if (contig.placed) continue;
        for (const LayoutAnchor& anchor : contig.anchors) {
            fallback_nodes.insert((static_cast<uint64_t>(contig.component) << 32) | anchor.segment);
        }
    }

    //  Input: contig -> node
    // Output:   node -> contig
    std::unordered_map<uint64_t, std::vector<Occurrence>> occurrences;
    occurrences.reserve(fallback_nodes.size());
    for (uint32_t i = 0; i < contigs.size(); ++i) {
        for (const auto& anchor : contigs[i].anchors) {
            const uint64_t key = (static_cast<uint64_t>(contigs[i].component) << 32) | anchor.segment;
            if (!fallback_nodes.count(key)) continue;
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
                if (contigs[a].placed && contigs[b].placed) continue;
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
    for (uint32_t i = 0; i < contigs.size(); ++i) placed[i] = contigs[i].placed;

    using QueueItem = std::tuple<uint64_t, uint32_t, uint32_t, uint32_t>;  // (shared_bp, relation_id, from, to)

    for (const auto& component : components) {
        int64_t cursor = 0;
        size_t remaining = 0;
        bool backbone_seeded = false;
        for (uint32_t index : component.second) {
            if (placed[index]) {
                backbone_seeded = true;
                cursor = std::max(cursor, contigs[index].start + static_cast<int64_t>(contigs[index].bp));
            } else {
                ++remaining;
            }
        }

        // ------------------------------------------------ Expand from all fixed backbone placements ------------------------------------------------
        std::priority_queue<QueueItem> queue;
        auto queue_from = [&](uint32_t from) {
            for (uint32_t relation : adjacency[from]) {
                const PlotRelation& edge = relations[relation];
                const uint32_t to = edge.a == from ? edge.b : edge.a;
                if (!placed[to]) queue.emplace(edge.shared_bp, relation, from, to);
            }
        };
        auto expand = [&]() {
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
                const int64_t forward_start = layout_median(forward);
                const int64_t reverse_start = layout_median(reverse);
                contigs[to].reverse = layout_deviation(reverse, reverse_start) <
                    layout_deviation(forward, forward_start);
                contigs[to].start = contigs[to].reverse ? reverse_start : forward_start;
                contigs[to].placed = true;
                placed[to] = 1;
                --remaining;
                queue_from(to);
            }
        };

        if (backbone_seeded) {
            for (uint32_t index : component.second) {
                if (placed[index]) queue_from(index);
            }
            expand();
        } else {
            cursor = 0;
        }

        bool have_layout = backbone_seeded;
        while (remaining > 0) {
            if (have_layout) {
                int64_t end = cursor;
                for (uint32_t index : component.second) {
                    if (placed[index]) end = std::max(
                        end, contigs[index].start + static_cast<int64_t>(contigs[index].bp)
                    );
                }
                cursor = end + std::max<int64_t>(1'000'000, end / 100);
            }

            // ------------------------------------------------ Start a new layout chain ------------------------------------------------
            uint32_t root = UINT32_MAX;
            for (uint32_t index : component.second) {
                if (!placed[index] && (root == UINT32_MAX || contigs[index].bp > contigs[root].bp)) root = index;
            }
            contigs[root].start = cursor;
            contigs[root].reverse = false;
            contigs[root].placed = true;
            placed[root] = 1;
            --remaining;
            queue_from(root);
            expand();
            have_layout = true;
        }

        // Preserve absolute backbone coordinates; normalize reference-free layouts.
        if (!backbone_seeded) {
            int64_t begin = contigs[component.second.front()].start;
            for (uint32_t index : component.second) begin = std::min(begin, contigs[index].start);
            for (uint32_t index : component.second) contigs[index].start -= begin;
        }
    }
}


// ================================================= Shared-node candidate evidence =================================================
std::vector<GfaGapfill::Candidate> GfaGapfill::build_candidates_(
    bool include_long_gaps
) const {
    // ------------------------------------------------ Group contigs by sample ------------------------------------------------
    std::map<std::string, std::vector<uint32_t>> fragments_by_sample;
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (!fragment.eligible || fragment.redundant) continue;
        fragments_by_sample[fragment.sample].push_back(i);
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
                    false
                );
                const Boundary right_boundary = boundary.inspect(
                    a, overlap.a_end, b, overlap.b_end,
                    true
                );
                const GfaGapfillDebugger::BoundarySupport left_support_debug{
                    left_boundary.phase.similarity,
                    left_boundary.phase.target_bp,
                    left_boundary.phase.bridge_bp
                };
                const GfaGapfillDebugger::BoundarySupport right_support_debug{
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
                sample_groups[sample_id], overlaps, include_long_gaps
            );
            progress.hit();
            return batch;
        }));
    }

    // ------------------------------------------------ Merge sample results ------------------------------------------------
    std::vector<Candidate> candidates, phase_edges;
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
            batch.wrong_group, batch.anchor_failed,
            batch.candidates.size() + batch.phase_edges.size()
        );
        candidates.insert(
            candidates.end(),
            std::make_move_iterator(batch.candidates.begin()),
            std::make_move_iterator(batch.candidates.end())
        );
        phase_edges.insert(
            phase_edges.end(),
            std::make_move_iterator(batch.phase_edges.begin()),
            std::make_move_iterator(batch.phase_edges.end())
        );
    }

    // Misassembly detection inspects every valid spanning connection. Formal
    // gap filling uses cross-haplotype edges only as assignment evidence.
    if (include_long_gaps) {
        candidates.insert(
            candidates.end(),
            std::make_move_iterator(phase_edges.begin()),
            std::make_move_iterator(phase_edges.end())
        );
    } else {
        group_full_phase_candidates_(candidates, phase_edges);
    }

    log_stream() << "  - Candidate connections: " << candidates.size() << "\n\n";
    return candidates;
}

GfaGapfill::CandidateBatch GfaGapfill::build_sample_candidates_(
    const std::vector<uint32_t>& sample_fragments,
    const OverlapIndex& overlaps,
    bool include_long_gaps
) const {
    // ------------------------------------------------ Prepare sample search ------------------------------------------------
    CandidateBatch batch;
    const GfaGapfillBoundary boundary(*this);

    // ------------------------------------------------ Try each left contig ------------------------------------------------
    for (uint32_t left_id : sample_fragments) {
        const Fragment& left = fragments_[left_id];
        if (left.length < params_.min_contig) continue;
        const bool target_full_phase = full_phase_samples_.count(left.sample);

        // ------------------------------------------------ Find bridge and right contigs ------------------------------------------------
        // Keep all bridge paths for each right contig.
        /*
         S.left -> T.bridge1 -> S.right
         S.left -> T.bridge2 -> S.right
         S.left -> U.bridge3 -> S.right
        */
        std::map<uint32_t, std::vector<Candidate>> evidence_by_right;
        // Enumerate bridge paths that connect this left contig to a later contig.
        for (const auto& bridge_item : overlaps[left_id]) {
            const uint32_t bridge_id = bridge_item.first;
            const Fragment& bridge = fragments_[bridge_id];
            if (bridge.sample == left.sample) continue;

            /**************************** may have bugs ********************************/
            if (target_full_phase && !full_phase_samples_.count(bridge.sample)) continue;  // Skip non-full-phase bridges if the left contig is full-phase
            /**************************** may have bugs ********************************/

            // The left side must share enough nodes with the bridge.
            const OverlapSupport left_support = overlap_support_(
                left, bridge, bridge_item.second
            );
            if (left_support.shared_bp < params_.min_overlap) {
                ++batch.low_support;
                continue;
            }

            // Find the right contig through this bridge.
            for (const auto& right_item : overlaps[bridge_id]) {
                const uint32_t right_id = right_item.first;
                const Fragment& right = fragments_[right_id];
                ++batch.tested;

                // Left and right must belong to the same target path group.
                if (right_id == left_id || right.sample != left.sample || right.component != left.component ||
                    !ordered_pair_(left, right, include_long_gaps)) {
                    ++batch.wrong_group;
                    continue;
                }
                if (!target_full_phase && right.hap != left.hap) {
                    ++batch.wrong_group;
                    continue;
                }
                if (right.hap == left.hap && skips_backbone_(left_id, right_id)) {
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

                // The right side must also have enough shared nodes.
                const OverlapSupport right_support = overlap_support_(
                    bridge, right, right_item.second
                );
                if (right_support.shared_bp < params_.min_overlap) {
                    ++batch.low_support;
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
                // Boundary fixes the trusted shared nodes used by the later graph walk.
                if (!boundary.prepare(left_id, right_id, bridge_id, candidate.left_boundary, candidate.right_boundary)) {
                    ++batch.anchor_failed;
                    continue;
                }
                const uint64_t fill_bp = candidate.right_boundary.bridge_cut_bp -
                    candidate.left_boundary.bridge_cut_bp;
                if (!include_long_gaps && fill_bp > params_.max_gap) {
                    ++batch.wrong_group;
                    continue;
                }
                evidence_by_right[right_id].push_back(candidate);
            }
        }

        // ------------------------------------------------ Compare bridge haplotypes ------------------------------------------------
        for (auto& item : evidence_by_right) {
            std::sort(item.second.begin(), item.second.end(), [&](const Candidate& a, const Candidate& b) {
                const Fragment& ab = fragments_[a.bridge];
                const Fragment& bb = fragments_[b.bridge];
                return std::tie(ab.sample, ab.hap, a.bridge) < std::tie(bb.sample, bb.hap, b.bridge);
            });

            // check if the left and right contigs are spanned by a homologous contig in this sample
            const bool homolog_span = fragments_[item.first].hap == left.hap &&
                homolog_spans_(left_id, item.first, overlaps);

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

            // ------------------------------------------------ Keep the required spanning evidence ------------------------------------------------
            for (Candidate& candidate : item.second) {
                candidate.homolog_span = homolog_span;
            }

            if (target_full_phase) {
                // Reciprocal target connections must be supported by the two
                // phased haplotypes of one donor sample. Keep one bridge for
                // each donor haplotype until those pairs are formed.
                std::map<std::pair<std::string, uint8_t>, size_t> best;
                for (size_t i = 0; i < item.second.size(); ++i) {
                    const Fragment& bridge = fragments_[item.second[i].bridge];
                    if (bridge.hap == 0 || bridge.hap > 2) continue;
                    const auto key = std::make_pair(bridge.sample, bridge.hap);
                    const auto found = best.find(key);
                    if (found == best.end() || candidate_better_(item.second[i], item.second[found->second])) {
                        best[key] = i;
                    }
                }
                for (const auto& selected : best) {
                    Candidate candidate = std::move(item.second[selected.second]);
                    if (fragments_[candidate.left].hap == fragments_[candidate.right].hap) {
                        batch.candidates.push_back(std::move(candidate));
                    } else {
                        batch.phase_edges.push_back(std::move(candidate));
                    }
                }
            } else {
                size_t best = 0;
                for (size_t i = 1; i < item.second.size(); ++i) {
                    if (candidate_better_(item.second[i], item.second[best])) best = i;
                }
                batch.candidates.push_back(std::move(item.second[best]));
            }
        }
    }
    return batch;
}

bool GfaGapfill::ordered_pair_(
    const Fragment& left,
    const Fragment& right,
    bool include_long_gaps
) const {
    if (left.sample != right.sample || left.component != right.component || left.length < params_.min_contig || right.length < params_.min_contig) return false;
    if (left.layout_start > right.layout_start) return false;
    const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
    const int64_t right_end = right.layout_start + static_cast<int64_t>(right.length);
    if (right.layout_start >= left_end) {
        return include_long_gaps || static_cast<uint64_t>(right.layout_start - left_end) <= params_.max_gap;
    }
    const uint64_t overlap = static_cast<uint64_t>(std::min(left_end, right_end) - right.layout_start);
    const double fraction = static_cast<double>(overlap) / std::min(left.length, right.length);
    return fraction <= params_.max_overlap;
}

bool GfaGapfill::skips_backbone_(uint32_t left_id, uint32_t right_id) const {
    const Fragment& left = fragments_[left_id];
    const Fragment& right = fragments_[right_id];
    if (left.hap < 1 || left.hap > 2) return false;

    const int64_t gap_begin = left.layout_start + static_cast<int64_t>(left.length);
    const int64_t gap_end = right.layout_start;
    if (gap_begin >= gap_end) return false;

    // fragments_ is already sorted by sample, component and layout. A long
    // contig that occupies the open gap must be connected first. Contigs
    // nested inside either endpoint do not block this connection.
    const uint32_t begin = std::min(left_id, right_id) + 1;
    const uint32_t end = std::max(left_id, right_id);
    for (uint32_t id = begin; id < end; ++id) {
        const Fragment& middle = fragments_[id];
        if (middle.sample != left.sample || middle.component != left.component ||
            middle.hap != left.hap || middle.length < params_.min_contig ||
            !middle.eligible || middle.redundant) continue;
        const int64_t middle_begin = middle.layout_start;
        const int64_t middle_end = middle_begin + static_cast<int64_t>(middle.length);
        if (std::max(gap_begin, middle_begin) < std::min(gap_end, middle_end)) return true;
    }
    return false;
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
        if (left_support.shared_bp < params_.min_overlap || right_support.shared_bp < params_.min_overlap || left_support.similarity < params_.min_similarity || right_support.similarity < params_.min_similarity) continue;

        const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
        const int64_t bridge_end = bridge.layout_start + static_cast<int64_t>(bridge.length);
        if (bridge.layout_start > left_end || bridge_end < right.layout_start) continue;
        if (!boundary.spans(left_id, right_id, bridge_id)) continue;
        return true;
    }
    return false;
}

double GfaGapfill::raw_phase_score_(const Candidate& candidate) {
    const double left = candidate.left_boundary.phase.similarity;
    const double right = candidate.right_boundary.phase.similarity;
    if (left >= 0.0 && right >= 0.0) return std::sqrt(left * right);
    return -1.0;
}

void GfaGapfill::group_full_phase_candidates_(
    std::vector<Candidate>& candidates,
    const std::vector<Candidate>& phase_edges
) const {
    // ------------------------------------------------ Collect haplotype-internal gaps ------------------------------------------------
    struct Gap {
        uint32_t left{UINT32_MAX}, right{UINT32_MAX};
        uint8_t hap{0};
        int64_t begin{0}, end{0};
        bool homolog_span{false};
        std::vector<size_t> edges;
    };
    std::vector<Gap> gaps;
    std::map<std::pair<uint32_t, uint32_t>, size_t> gap_by_ends;
    std::map<std::pair<uint32_t, uint32_t>, std::vector<size_t>> direct_by_ends;
    std::map<std::pair<uint32_t, uint32_t>, std::vector<size_t>> cross_by_ends;
    std::map<std::pair<std::string, uint32_t>, std::vector<size_t>> gaps_by_component;

    for (size_t id = 0; id < candidates.size(); ++id) {
        const Candidate& candidate = candidates[id];
        const Fragment& left = fragments_[candidate.left];
        const Fragment& right = fragments_[candidate.right];
        if (!full_phase_samples_.count(left.sample) || left.hap != right.hap ||
            left.hap < 1 || left.hap > 2) continue;
        const uint64_t bridge_bp = candidate.right_boundary.bridge_cut_bp -
            candidate.left_boundary.bridge_cut_bp;
        if ((!candidate.target_overlaps && candidate.target_distance > params_.max_gap) ||
            bridge_bp > params_.max_gap) continue;

        const auto key = std::make_pair(candidate.left, candidate.right);
        direct_by_ends[key].push_back(id);
        auto found = gap_by_ends.find(key);
        if (found == gap_by_ends.end()) {
            const int64_t left_end = left.layout_start + static_cast<int64_t>(left.length);
            const int64_t right_begin = right.layout_start;
            const size_t gap_id = gaps.size();
            gaps.push_back({
                candidate.left, candidate.right, left.hap,
                std::min(left_end, right_begin), std::max(left_end, right_begin),
                candidate.homolog_span, {id}
            });
            gap_by_ends.emplace(key, gap_id);
            gaps_by_component[{left.sample, left.component}].push_back(gap_id);
        } else {
            Gap& gap = gaps[found->second];
            gap.homolog_span = gap.homolog_span || candidate.homolog_span;
            gap.edges.push_back(id);
        }
    }
    for (size_t id = 0; id < phase_edges.size(); ++id) {
        cross_by_ends[{phase_edges[id].left, phase_edges[id].right}].push_back(id);
    }
    if (gaps.empty()) return;

    // ------------------------------------------------ Form overlapping diploid loci ------------------------------------------------
    struct GapPair {
        size_t hap1{SIZE_MAX}, hap2{SIZE_MAX};
        uint32_t locus{UINT32_MAX};
    };
    std::vector<size_t> parent(gaps.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto root = [&](size_t id) {
        while (parent[id] != id) {
            parent[id] = parent[parent[id]];
            id = parent[id];
        }
        return id;
    };
    const auto join = [&](size_t a, size_t b) {
        a = root(a);
        b = root(b);
        if (a != b) parent[b] = a;
    };
    std::vector<uint8_t> paired(gaps.size(), 0);
    for (auto& group : gaps_by_component) {
        std::vector<size_t>& ids = group.second;
        std::sort(ids.begin(), ids.end(), [&](size_t a, size_t b) {
            return std::tie(gaps[a].begin, gaps[a].end, gaps[a].hap, a) <
                std::tie(gaps[b].begin, gaps[b].end, gaps[b].hap, b);
        });
        for (size_t i = 0; i < ids.size(); ++i) {
            const Gap& first = gaps[ids[i]];
            if (first.homolog_span) continue;
            for (size_t j = i + 1; j < ids.size() && gaps[ids[j]].begin < first.end; ++j) {
                const Gap& second = gaps[ids[j]];
                if (second.homolog_span || first.hap == second.hap ||
                    std::max(first.begin, second.begin) >= std::min(first.end, second.end)) continue;

                std::array<uint32_t, 4> ends{first.left, first.right, second.left, second.right};
                std::sort(ends.begin(), ends.end());
                if (std::adjacent_find(ends.begin(), ends.end()) != ends.end()) continue;

                join(ids[i], ids[j]);
                paired[ids[i]] = paired[ids[j]] = 1;
            }
        }
    }

    // A simple two-gap locus has one unambiguous diploid decision. If several
    // gaps overlap transitively, keep their same-haplotype candidates instead
    // of discarding valid connections or guessing a larger matching.
    std::map<size_t, std::vector<size_t>> gaps_by_locus;
    for (size_t id = 0; id < gaps.size(); ++id) {
        if (paired[id]) gaps_by_locus[root(id)].push_back(id);
    }
    std::vector<GapPair> gap_pairs;
    uint32_t next_locus = 0;
    for (const auto& item : gaps_by_locus) {
        if (item.second.size() != 2) continue;
        const Gap& first = gaps[item.second[0]];
        const Gap& second = gaps[item.second[1]];
        if (first.hap == second.hap) continue;
        gap_pairs.push_back({
            first.hap == 1 ? item.second[0] : item.second[1],
            first.hap == 2 ? item.second[0] : item.second[1],
            next_locus++
        });
    }
    if (gap_pairs.empty()) return;

    // ------------------------------------------------ Enumerate straight and crossed assignments ------------------------------------------------
    struct PairOption {
        Candidate first, second;
        uint32_t locus{UINT32_MAX};
        uint32_t phase_evidence{0};
        double phase{0.5}, relative_phase{0.0};
        uint64_t shared{0}, distance{0};
        bool crossed{false};
    };
    std::vector<PairOption> options;
    std::vector<uint8_t> in_locus(candidates.size(), 0);

    const auto add_options = [&](
        const std::vector<Candidate>& source,
        const std::vector<size_t>* first_ids,
        const std::vector<size_t>* second_ids,
        uint32_t locus,
        uint64_t distance,
        bool crossed
    ) {
        if (first_ids == nullptr || second_ids == nullptr) return;
        for (size_t first_id : *first_ids) {
            const Candidate& first = source[first_id];
            const Fragment& first_bridge = fragments_[first.bridge];
            for (size_t second_id : *second_ids) {
                const Candidate& second = source[second_id];
                const Fragment& second_bridge = fragments_[second.bridge];
                if (first_bridge.sample != second_bridge.sample ||
                    first_bridge.hap < 1 || first_bridge.hap > 2 ||
                    second_bridge.hap < 1 || second_bridge.hap > 2 ||
                    first_bridge.hap == second_bridge.hap) continue;

                const double similarities[4] = {
                    first.left_boundary.phase.similarity,
                    first.right_boundary.phase.similarity,
                    second.left_boundary.phase.similarity,
                    second.right_boundary.phase.similarity
                };
                uint32_t evidence = 0;
                double phase = 0.0;
                for (double similarity : similarities) {
                    if (similarity < 0.0) continue;
                    phase += similarity;
                    ++evidence;
                }
                options.push_back({
                    first, second, locus, evidence,
                    evidence ? phase / evidence : 0.5,
                    std::sqrt(first.phase_score * second.phase_score),
                    first.left_shared + first.right_shared + second.left_shared + second.right_shared,
                    distance, crossed
                });
            }
        }
    };
    const auto edge_ids = [](const auto& index, uint32_t left, uint32_t right) {
        const auto found = index.find({left, right});
        return found == index.end() ? nullptr : &found->second;
    };

    for (const GapPair& pair : gap_pairs) {
        const Gap& hap1 = gaps[pair.hap1];
        const Gap& hap2 = gaps[pair.hap2];
        const uint64_t distance = static_cast<uint64_t>(
            std::llabs(hap1.begin - hap2.begin) + std::llabs(hap1.end - hap2.end)
        );

        const size_t option_begin = options.size();
        add_options(
            candidates,
            edge_ids(direct_by_ends, hap1.left, hap1.right),
            edge_ids(direct_by_ends, hap2.left, hap2.right),
            pair.locus, distance, false
        );
        add_options(
            phase_edges,
            edge_ids(cross_by_ends, hap1.left, hap2.right),
            edge_ids(cross_by_ends, hap2.left, hap1.right),
            pair.locus, distance, true
        );
        if (options.size() != option_begin) {
            for (size_t id : hap1.edges) in_locus[id] = 1;
            for (size_t id : hap2.edges) in_locus[id] = 1;
        }
    }

    // More observed boundaries are stronger than neutral missing values. The
    // phase similarity then chooses straight versus crossed assignments.
    std::sort(options.begin(), options.end(), [](const PairOption& a, const PairOption& b) {
        if (a.locus != b.locus) return a.locus < b.locus;
        if (a.phase_evidence != b.phase_evidence) return a.phase_evidence > b.phase_evidence;
        if (a.phase != b.phase) return a.phase > b.phase;
        if (a.relative_phase != b.relative_phase) return a.relative_phase > b.relative_phase;
        if (a.shared != b.shared) return a.shared > b.shared;
        if (a.crossed != b.crossed) return !a.crossed;
        if (a.distance != b.distance) return a.distance < b.distance;
        return std::tie(a.first.left, a.first.right, a.first.bridge,
                        a.second.left, a.second.right, a.second.bridge) <
            std::tie(b.first.left, b.first.right, b.first.bridge,
                     b.second.left, b.second.right, b.second.bridge);
    });

    // ------------------------------------------------ Replace locus edges with atomic alternatives ------------------------------------------------
    std::vector<Candidate> grouped;
    grouped.reserve(candidates.size() + 2 * options.size());
    for (size_t id = 0; id < candidates.size(); ++id) {
        if (!in_locus[id]) grouped.push_back(std::move(candidates[id]));
    }
    uint32_t next_group = 0;
    for (PairOption& option : options) {
        option.first.phase_group = next_group;
        option.second.phase_group = next_group++;
        option.first.phase_locus = option.locus;
        option.second.phase_locus = option.locus;
        grouped.push_back(std::move(option.first));
        grouped.push_back(std::move(option.second));
    }
    candidates.swap(grouped);
}


// ================================================= Candidate selection =================================================
std::vector<GfaGapfill::Candidate> GfaGapfill::select_and_walk_(
    std::vector<Candidate>& candidates
) {
    // Candidate indexes stay stable throughout selection and graph walking.
    std::vector<size_t> selected_ids;
    size_t failed = 0;
    do {
        selected_ids = select_candidates_(candidates);
        failed = resolve_fill_paths_(candidates, selected_ids);
        if (failed) {
            log_stream() << "Retrying after " << failed << " graph walk failure(s) ...\n";
        }
    } while (failed);

    std::vector<Candidate> selected;
    selected.reserve(selected_ids.size());
    for (size_t id : selected_ids) selected.push_back(candidates[id]);
    return selected;
}

std::vector<size_t> GfaGapfill::select_candidates_(
    std::vector<Candidate>& candidates
) {
    log_stream() << "Selecting unambiguous gap-fill connections ...\n";

    struct Unit {
        size_t first{SIZE_MAX};
        size_t second{SIZE_MAX};
    };
    std::vector<Unit> units;
    std::map<uint32_t, std::vector<size_t>> pairs;
    for (size_t id = 0; id < candidates.size(); ++id) {
        Candidate& candidate = candidates[id];
        candidate.walk_hap = 0;
        const bool crossover = candidate.left_boundary.bridge_pos == candidate.right_boundary.bridge_pos;
        candidate.right_boundary.target_cut_bp = fragments_[candidate.right].path_bp[
            candidate.right_boundary.target_pos + static_cast<uint32_t>(crossover)
        ];
        candidate.fill.reset();
        if (candidate.status != Candidate::NO_PATH && candidate.status != Candidate::PAIR_CONFLICT) {
            candidate.status = Candidate::PENDING;
        }
        if (candidate.status != Candidate::PENDING) continue;
        if (candidate.phase_group == UINT32_MAX) {
            units.push_back({id, SIZE_MAX});
            continue;
        }
        pairs[candidate.phase_group].push_back(id);
    }
    for (auto& item : pairs) {
        if (item.second.size() == 2 &&
            candidates[item.second[0]].phase_locus == candidates[item.second[1]].phase_locus) {
            units.push_back({item.second[0], item.second[1]});
        } else {
            for (size_t id : item.second) {
                candidates[id].status = Candidate::PAIR_CONFLICT;
            }
        }
    }

    // ------------------------------------------------ Sort complete connection units ------------------------------------------------
    const auto weakest = [&](const Unit& unit) {
        if (unit.second == SIZE_MAX) return unit.first;
        return candidate_better_(candidates[unit.first], candidates[unit.second]) ?
            unit.second : unit.first;
    };
    const auto stable_key = [&](const Unit& unit) {
        const Candidate& first_candidate = candidates[unit.first];
        std::array<uint32_t, 6> key{
            first_candidate.left, first_candidate.right, first_candidate.bridge,
            UINT32_MAX, UINT32_MAX, UINT32_MAX
        };
        if (unit.second != SIZE_MAX) {
            const Candidate& second_candidate = candidates[unit.second];
            std::array<uint32_t, 3> first{key[0], key[1], key[2]};
            std::array<uint32_t, 3> second{
                second_candidate.left, second_candidate.right, second_candidate.bridge
            };
            if (second < first) std::swap(first, second);
            std::copy(first.begin(), first.end(), key.begin());
            std::copy(second.begin(), second.end(), key.begin() + 3);
        }
        return key;
    };
    const auto boundary_phase = [&](const Unit& unit) {
        double total = 0.0;
        uint32_t evidence = 0, boundaries = 0;
        const size_t edges[2] = {unit.first, unit.second};
        for (size_t id : edges) {
            if (id == SIZE_MAX) continue;
            const Candidate& edge = candidates[id];
            const double values[2] = {
                edge.left_boundary.phase.similarity,
                edge.right_boundary.phase.similarity
            };
            for (double value : values) {
                total += value >= 0.0 ? value : 0.5;
                evidence += value >= 0.0;
                ++boundaries;
            }
        }
        return std::pair<double, uint32_t>{total / boundaries, evidence};
    };
    std::sort(units.begin(), units.end(), [&](const Unit& a, const Unit& b) {
        const Fragment& af = fragments_[candidates[a.first].left];
        const Fragment& bf = fragments_[candidates[b.first].left];
        if (af.sample != bf.sample) return af.sample < bf.sample;
        if (af.component != bf.component) return af.component < bf.component;

        const size_t aw = weakest(a);
        const size_t bw = weakest(b);
        // Ordinary samples keep the original r22 ordering. Full-phase
        // samples use one fixed unit ordering for single and paired edges.
        if (!full_phase_samples_.count(af.sample)) {
            if (connection_better_(candidates[aw], candidates[bw])) return true;
            if (connection_better_(candidates[bw], candidates[aw])) return false;
            return stable_key(a) < stable_key(b);
        }

        const bool a_homolog = candidates[a.first].homolog_span ||
            (a.second != SIZE_MAX && candidates[a.second].homolog_span);
        const bool b_homolog = candidates[b.first].homolog_span ||
            (b.second != SIZE_MAX && candidates[b.second].homolog_span);
        if (a_homolog != b_homolog) return a_homolog;

        const auto a_boundary_phase = boundary_phase(a);
        const auto b_boundary_phase = boundary_phase(b);
        if (a_boundary_phase.second != b_boundary_phase.second) {
            return a_boundary_phase.second > b_boundary_phase.second;
        }
        if (a_boundary_phase.first != b_boundary_phase.first) {
            return a_boundary_phase.first > b_boundary_phase.first;
        }
        const double a_relative_phase = a.second != SIZE_MAX ?
            std::sqrt(candidates[a.first].phase_score * candidates[a.second].phase_score) :
            candidates[a.first].phase_score;
        const double b_relative_phase = b.second != SIZE_MAX ?
            std::sqrt(candidates[b.first].phase_score * candidates[b.second].phase_score) :
            candidates[b.first].phase_score;
        if (a_relative_phase != b_relative_phase) return a_relative_phase > b_relative_phase;

        const bool a_crossed = fragments_[candidates[a.first].left].hap !=
            fragments_[candidates[a.first].right].hap;
        const bool b_crossed = fragments_[candidates[b.first].left].hap !=
            fragments_[candidates[b.first].right].hap;
        if (a_crossed != b_crossed) return !a_crossed;

        if (connection_better_(candidates[aw], candidates[bw])) return true;
        if (connection_better_(candidates[bw], candidates[aw])) return false;
        return stable_key(a) < stable_key(b);
    });

    // ------------------------------------------------ Prepare selection state ------------------------------------------------
    std::vector<size_t> selected;
    // Track used contig ends and connected components.
    std::vector<size_t> incoming(fragments_.size(), SIZE_MAX), outgoing(fragments_.size(), SIZE_MAX);
    DisjointSet sets(fragments_.size());
    std::unordered_set<uint32_t> used_phase_loci;

    // ------------------------------------------------ Check each connection unit ------------------------------------------------
    for (Unit& unit : units) {
        const size_t edges[2] = {unit.first, unit.second};
        const size_t count = unit.second != SIZE_MAX ? 2 : 1;
        Candidate::Status rejected = Candidate::PENDING;
        const char* decision = "keep";

        const uint32_t phase_locus = candidates[unit.first].phase_locus;
        if (phase_locus != UINT32_MAX && used_phase_loci.count(phase_locus)) {
            rejected = Candidate::COORDINATE_CONFLICT;
            decision = "drop:phase-locus";
        }

        // All ends and target intervals in a full-phase pair are checked before either edge is committed.
        for (size_t i = 0; i < count && rejected == Candidate::PENDING; ++i) {
            Candidate& edge = candidates[edges[i]];
            const Fragment& left = fragments_[edge.left];
            const Fragment& right = fragments_[edge.right];
            if (edge.target_overlaps && static_cast<double>(edge.target_distance) /
                    std::min(left.length, right.length) > params_.max_overlap) {
                rejected = Candidate::COORDINATE_CONFLICT;
                decision = "drop:long-overlap";
            } else if (outgoing[edge.left] != SIZE_MAX || incoming[edge.right] != SIZE_MAX) {
                rejected = Candidate::USED_END;
                decision = "drop:used-end";
            } else {
                const uint32_t left_entry = incoming[edge.left] == SIZE_MAX ? 0 :
                    candidates[incoming[edge.left]].right_boundary.target_pos;
                const uint32_t right_exit = outgoing[edge.right] == SIZE_MAX ?
                    static_cast<uint32_t>(right.vertices.size() - 1) :
                    candidates[outgoing[edge.right]].left_boundary.target_pos;
                if (left_entry > edge.left_boundary.target_pos ||
                    edge.right_boundary.target_pos > right_exit) {
                    rejected = Candidate::COORDINATE_CONFLICT;
                    decision = "drop:coordinate-conflict";
                }
            }
        }

        // Adding either edge must not close a cycle; paired edges are tested atomically.
        if (rejected == Candidate::PENDING) {
            const Candidate& first = candidates[unit.first];
            const uint32_t a = sets.find(first.left);
            const uint32_t b = sets.find(first.right);
            if (a == b) {
                rejected = Candidate::CYCLE;
                decision = "drop:cycle";
            } else if (unit.second != SIZE_MAX) {
                const Candidate& second = candidates[unit.second];
                const uint32_t c = sets.find(second.left);
                const uint32_t d = sets.find(second.right);
                if (c == d || ((c == a || c == b) && (d == a || d == b))) {
                    rejected = Candidate::CYCLE;
                    decision = "drop:cycle";
                }
            }
        }

        for (size_t i = 0; i < count; ++i) {
            const size_t id = edges[i];
            Candidate& edge = candidates[id];
            if (rejected != Candidate::PENDING) {
                edge.status = rejected;
            } else {
                edge.status = Candidate::KEPT;
                selected.push_back(id);
                outgoing[edge.left] = id;
                incoming[edge.right] = id;
                sets.join(edge.left, edge.right);
            }
            GfaGapfillDebugger::selection(
                paths_[fragments_[edge.left].path_id].name,
                paths_[fragments_[edge.right].path_id].name,
                paths_[fragments_[edge.bridge].path_id].name,
                decision
            );
        }
        if (rejected == Candidate::PENDING && phase_locus != UINT32_MAX) {
            used_phase_loci.insert(phase_locus);
        }
    }

    // Propagate the output haplotype from each chain start.  After a
    // full-phase exchange this follows the assembled chain, not the original
    // h1/h2 label of each downstream contig.
    for (size_t fragment = 0; fragment < fragments_.size(); ++fragment) {
        if (incoming[fragment] != SIZE_MAX || outgoing[fragment] == SIZE_MAX) continue;
        const uint8_t hap = fragments_[fragment].hap;
        if (hap < 1 || hap > 2) continue;
        size_t current = fragment;
        while (outgoing[current] != SIZE_MAX) {
            const size_t id = outgoing[current];
            candidates[id].walk_hap = hap;
            current = candidates[id].right;
        }
    }
    log_stream() << "  - Connections selected: " << selected.size() << "\n\n";
    return selected;
}

bool GfaGapfill::candidate_better_(const Candidate& a, const Candidate& b) {
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
    return std::tie(a.left, a.right, a.bridge,
                    a.left_boundary.target_pos, a.right_boundary.target_pos) <
        std::tie(b.left, b.right, b.bridge,
                 b.left_boundary.target_pos, b.right_boundary.target_pos);
}

bool GfaGapfill::connection_better_(const Candidate& a, const Candidate& b) {
    // A complete contig on the other haplotype fixes the target pairing before
    // local gap evidence is used to resolve the remaining ambiguous ends.
    if (a.homolog_span != b.homolog_span) return a.homolog_span;
    if (a.target_overlaps != b.target_overlaps) return !a.target_overlaps;
    if (a.target_distance != b.target_distance) return a.target_distance < b.target_distance;
    return candidate_better_(a, b);
}


// ================================================= Node-level gap filling =================================================
size_t GfaGapfill::resolve_fill_paths_(
    std::vector<Candidate>& candidates,
    std::vector<size_t>& selected
) const {
    log_stream() << "Walking graph nodes between selected Boundary anchors ...\n";

    // ------------------------------------------------ Stable sample order ------------------------------------------------
    std::sort(selected.begin(), selected.end(), [&](size_t a_id, size_t b_id) {
        const Candidate& a = candidates[a_id];
        const Candidate& b = candidates[b_id];
        const Fragment& al = fragments_[a.left];
        const Fragment& bl = fragments_[b.left];
        if (al.sample != bl.sample) return al.sample < bl.sample;

        if (full_phase_samples_.count(al.sample)) {
            if (al.component != bl.component) return al.component < bl.component;
            if (a.phase_group != b.phase_group) return a.phase_group < b.phase_group;
            if (al.hap != bl.hap) return al.hap < bl.hap;
        } else {
            const uint8_t ah = al.hap == 1 ? 0 : (al.hap == 2 ? 1 : 2);
            const uint8_t bh = bl.hap == 1 ? 0 : (bl.hap == 2 ? 1 : 2);
            if (ah != bh) return ah < bh;
            if (al.component != bl.component) return al.component < bl.component;
        }
        return std::tie(al.layout_start, a.left, a.right, a_id) <
            std::tie(bl.layout_start, b.left, b.right, b_id);
    });

    // ------------------------------------------------ Group work and original haplotype nodes by sample ------------------------------------------------
    struct SampleWork {
        size_t begin{0}, end{0};
        std::vector<uint32_t> fragments;
    };
    std::vector<SampleWork> work;
    std::unordered_map<std::string, size_t> work_by_sample;
    for (size_t begin = 0; begin < selected.size();) {
        const std::string& sample = fragments_[candidates[selected[begin]].left].sample;
        size_t end = begin + 1;
        while (end < selected.size() && fragments_[candidates[selected[end]].left].sample == sample) ++end;
        work_by_sample.emplace(sample, work.size());
        work.push_back({begin, end, {}});
        begin = end;
    }
    for (uint32_t fragment_id = 0; fragment_id < fragments_.size(); ++fragment_id) {
        const Fragment& fragment = fragments_[fragment_id];
        if (fragment.hap < 1 || fragment.hap > 2) continue;
        const auto found = work_by_sample.find(fragment.sample);
        if (found != work_by_sample.end()) work[found->second].fragments.push_back(fragment_id);
    }

    // ------------------------------------------------ Resolve samples in parallel ------------------------------------------------
    size_t resolved = 0;
    if (!work.empty()) {
        ThreadPool pool(std::min<size_t>(params_.threads, work.size()));
        ProgressTracker progress(work.size());
        std::vector<std::future<size_t>> futures;
        futures.reserve(work.size());
        for (size_t work_id = 0; work_id < work.size(); ++work_id) {
            futures.emplace_back(pool.submit([&, work_id] {
                SampleWork& item = work[work_id];
                const size_t count = resolve_sample_fill_paths_(
                    candidates, selected, item.begin, item.end, item.fragments
                );
                progress.hit();
                return count;
            }));
        }
        for (std::future<size_t>& future : futures) resolved += future.get();
        progress.finish();
        pool.stop();
    }

    // ------------------------------------------------ Report resolved walks serially ------------------------------------------------
    for (size_t id : selected) {
        const Candidate& candidate = candidates[id];
        if (!candidate.fill) continue;
        const Fragment& left = fragments_[candidate.left];
        GfaGapfillDebugger::candidate(
            paths_[left.path_id].name,
            paths_[fragments_[candidate.right].path_id].name,
            paths_[fragments_[candidate.bridge].path_id].name,
            left.component,
            candidate.left_boundary.phase.similarity,
            candidate.right_boundary.phase.similarity,
            candidate.left_boundary.phase.target_bp,
            candidate.left_boundary.phase.bridge_bp,
            candidate.right_boundary.phase.target_bp,
            candidate.right_boundary.phase.bridge_bp,
            candidate.homolog_span,
            candidate.phase_score
        );
        std::vector<std::string> walk_vertices;
        if (DEBUG_ENABLED) {
            walk_vertices.reserve(candidate.fill->vertices.size());
            for (uint32_t vertex : candidate.fill->vertices) {
                walk_vertices.push_back(
                    nodes_[Vertex::get_segment_id(vertex)].name +
                    (Vertex::get_is_reverse(vertex) ? "-" : "+")
                );
            }
        }
        GfaGapfillDebugger::node_walk(
            paths_[left.path_id].name,
            paths_[fragments_[candidate.right].path_id].name,
            paths_[fragments_[candidate.bridge].path_id].name,
            candidate.fill->length,
            candidate.fill->left_cut,
            candidate.fill->right_cut,
            walk_vertices
        );
    }

    const size_t failed = selected.size() - resolved;
    log_stream() << "  - Node walks resolved: " << resolved << '\n';
    log_stream() << "  - Node walks rejected: " << failed << "\n\n";
    return failed;
}

size_t GfaGapfill::resolve_sample_fill_paths_(
    std::vector<Candidate>& candidates,
    const std::vector<size_t>& selected,
    size_t begin,
    size_t end,
    const std::vector<uint32_t>& sample_fragments
) const {
    // Original haplotype nodes are the initial soft masks. Filled internal
    // nodes are added only after a single edge or complete atomic pair succeeds.
    std::vector<uint8_t> used(nodes_.size(), 0); // bit 0: hap1, bit 1: hap2
    for (uint32_t fragment_id : sample_fragments) {
        const Fragment& fragment = fragments_[fragment_id];
        const uint8_t bit = static_cast<uint8_t>(1u << (fragment.hap - 1));
        for (uint32_t vertex : fragment.vertices) {
            used[Vertex::get_segment_id(vertex)] |= bit;
        }
    }

    const auto add_internal = [](const FillPath& fill, std::vector<uint8_t>& used, uint8_t bit) {
        for (size_t i = 1; i + 1 < fill.vertices.size(); ++i) {
            used[Vertex::get_segment_id(fill.vertices[i])] |= bit;
        }
    };
    size_t resolved = 0;

    // ------------------------------------------------ Resolve single connections and atomic pairs ------------------------------------------------
    for (size_t i = begin; i < end;) {
        const Candidate& current = candidates[selected[i]];
        size_t mate = SIZE_MAX;
        if (current.phase_group != UINT32_MAX && i + 1 < end &&
            candidates[selected[i + 1]].phase_group == current.phase_group) {
            mate = i + 1;
        }

        const size_t count = mate == SIZE_MAX ? 1 : 2;
        size_t edges[2] = {selected[i], mate == SIZE_MAX ? SIZE_MAX : selected[mate]};
        if (count == 2) {
            const bool second_better = candidate_better_(candidates[edges[1]], candidates[edges[0]]);
            const bool tied = !second_better && !candidate_better_(candidates[edges[0]], candidates[edges[1]]);
            if (second_better || (tied && fragments_[candidates[edges[1]].left].hap < fragments_[candidates[edges[0]].left].hap)) {
                std::swap(edges[0], edges[1]);
            }
        }

        Candidate first = candidates[edges[0]];
        const uint8_t first_hap = first.walk_hap ? first.walk_hap : fragments_[first.left].hap;
        bool ok = first_hap >= 1 && first_hap <= 2;
        if (ok) {
            const uint8_t blocked_bit = static_cast<uint8_t>(1u << (2 - first_hap));
            ok = resolve_fill_path_(first, used, blocked_bit);
        }

        Candidate second;
        if (ok && count == 2) {
            second = candidates[edges[1]];
            const uint8_t second_hap = second.walk_hap ? second.walk_hap : fragments_[second.left].hap;
            ok = second_hap >= 1 && second_hap <= 2;
            if (ok) {
                const uint8_t blocked_bit = static_cast<uint8_t>(1u << (2 - second_hap));
                std::vector<uint32_t> temporary;
                for (size_t pos = 1; pos + 1 < first.fill->vertices.size(); ++pos) {
                    const uint32_t segment = Vertex::get_segment_id(first.fill->vertices[pos]);
                    if (!(used[segment] & blocked_bit)) {
                        used[segment] |= blocked_bit;
                        temporary.push_back(segment);
                    }
                }
                ok = resolve_fill_path_(second, used, blocked_bit);
                for (uint32_t segment : temporary) used[segment] &= static_cast<uint8_t>(~blocked_bit);
            }
        }

        if (!ok) {
            for (size_t edge = 0; edge < count; ++edge) {
                candidates[edges[edge]].fill.reset();
                candidates[edges[edge]].status = Candidate::NO_PATH;
            }
        } else {
            candidates[edges[0]] = std::move(first);
            if (count == 2) candidates[edges[1]] = std::move(second);
            for (size_t edge = 0; edge < count; ++edge) {
                Candidate& candidate = candidates[edges[edge]];
                const uint8_t hap = candidate.walk_hap ? candidate.walk_hap : fragments_[candidate.left].hap;
                add_internal(*candidate.fill, used, static_cast<uint8_t>(1u << (hap - 1)));
                ++resolved;
            }
        }

        i += count;
    }
    return resolved;
}

bool GfaGapfill::resolve_fill_path_(
    Candidate& candidate,
    const std::vector<uint8_t>& soft_blocked,
    uint8_t soft_block_mask
) const {
    const Fragment& left = fragments_[candidate.left];
    const Fragment& right = fragments_[candidate.right];
    const Fragment& bridge = fragments_[candidate.bridge];
    const uint32_t source = left.vertices[candidate.left_boundary.target_pos];
    const uint32_t sink = right.vertices[candidate.right_boundary.target_pos];

    // ------------------------------------------------ Configure source-aware walk ------------------------------------------------
    GfaSourceWalker::Options options;
    const std::string bridge_source = bridge.sample + ".hap" + std::to_string(bridge.hap);
    const auto named_source = source_ids_.find(bridge_source);
    if (named_source != source_ids_.end()) options.preferred_source = named_source->second;

    // Use the dominant local SN carried by the real spanning bridge. This is
    // robust when P-path sample names and node provenance names are different.
    std::vector<double> bridge_source_weight(getSampleNames().size(), 0.0);
    for (uint32_t pos = candidate.left_boundary.bridge_pos + 1; pos < candidate.right_boundary.bridge_pos; ++pos) {
        const uint32_t segment = Vertex::get_segment_id(bridge.vertices[pos]);
        std::vector<uint32_t> sources;
        for (uint32_t source_id : getNodeSampleIds(segment)) {
            if (source_scope_.empty() || source_scope_.count(source_id)) sources.push_back(source_id);
        }
        if (sources.empty()) continue;
        const double share = static_cast<double>(
            bridge.path_bp[pos + 1] - bridge.path_bp[pos]
        ) / sources.size();
        for (uint32_t source_id : sources) bridge_source_weight[source_id] += share;
    }
    for (uint32_t source_id = 0; source_id < bridge_source_weight.size(); ++source_id) {
        if (bridge_source_weight[source_id] > 0.0 && (options.preferred_source == UINT32_MAX || bridge_source_weight[source_id] > bridge_source_weight[options.preferred_source])) {
            options.preferred_source = source_id;
        }
    }
    options.source_scope = source_scope_.empty() ? nullptr : &source_scope_;
    options.transition_sources = &source_transitions_;
    options.soft_blocked = soft_blocked.empty() ? nullptr : &soft_blocked;
    options.soft_block_mask = soft_block_mask;
    options.max_depth = params_.max_depth;
    options.max_bp = params_.max_gap;
    options.max_states = params_.DFS_guard;

    const GfaSourceWalker::Result result = GfaSourceWalker(*this).walk(source, sink, options);
    if (!result || result.truncated) return false;

    // ------------------------------------------------ Build node walk sequence ------------------------------------------------
    auto fill = std::make_shared<FillPath>();
    fill->vertices = result.vertices;
    for (uint32_t vertex : result.vertices) {
        fill->length += nodes_[Vertex::get_segment_id(vertex)].length;
    }
    fill->sequence = get_path_sequence(result.vertices);
    if (fill->sequence == "*") fill->sequence.clear();

    // Source belongs to the left target; sink belongs to the right target.
    if (source == sink) {
        fill->left_cut = fill->right_cut = fill->length;
        candidate.right_boundary.target_cut_bp = right.path_bp[candidate.right_boundary.target_pos + 1];
    } else {
        fill->left_cut = nodes_[Vertex::get_segment_id(source)].length;
        fill->right_cut = fill->length - nodes_[Vertex::get_segment_id(sink)].length;
    }
    if (fill->left_cut > fill->right_cut || fill->right_cut - fill->left_cut > params_.max_gap || fill->right_cut > fill->length || candidate.right_boundary.target_cut_bp > right.length) return false;

    candidate.fill = std::move(fill);
    return true;
}


// ================================================= Misassembly detection and relocation =================================================
void GfaGapfill::mark_used_relocations_(const std::vector<Candidate>& selected) {
    for (const Candidate& candidate : selected) {
        const uint32_t ids[2] = {candidate.left, candidate.right};
        for (uint32_t id : ids) {
            const int32_t relocation = fragments_[id].relocation;
            if (relocation >= 0) relocations_[relocation].used_for_gap = true;
        }
        if (!candidate.fill) continue;
        std::unordered_set<uint32_t> fill_segments;
        for (size_t i = 1; i + 1 < candidate.fill->vertices.size(); ++i) {
            fill_segments.insert(Vertex::get_segment_id(candidate.fill->vertices[i]));
        }
        for (const Fragment& fragment : fragments_) {
            if (fragment.relocation < 0) continue;
            for (uint32_t vertex : fragment.vertices) {
                if (fill_segments.count(Vertex::get_segment_id(vertex))) {
                    relocations_[fragment.relocation].used_for_gap = true;
                    break;
                }
            }
        }
    }
}


// ================================================= Assign haplotypes =================================================
void GfaGapfill::build_sample_chains_(const std::vector<Candidate>& selected) {
    log_stream() << "Building gap-filled sample paths ...\n\n";

    // ------------------------------------------------ Index selected connections ------------------------------------------------
    std::vector<int32_t> incoming(fragments_.size(), -1), outgoing(fragments_.size(), -1);
    for (size_t i = 0; i < selected.size(); ++i) {
        incoming[selected[i].right] = static_cast<int32_t>(i);
        outgoing[selected[i].left] = static_cast<int32_t>(i);
    }

    // ------------------------------------------------ Build chains in their selected left-end phase ------------------------------------------------
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (incoming[i] >= 0 || !fragment.eligible || fragment.redundant ||
            fragment.hap == 0 || fragment.hap > 2) continue;
        Chain chain = build_chain_(i, incoming, outgoing, selected);
        if (!chain.vertices.empty()) {
            sample_chains_[chain.sample][chain.hap - 1].push_back(std::move(chain));
        }
    }

    // ------------------------------------------------ Keep short single fragments ------------------------------------------------
    for (uint32_t i = 0; i < fragments_.size(); ++i) {
        const Fragment& fragment = fragments_[i];
        if (fragment.eligible || fragment.relocation >= 0 || fragment.redundant ||
            fragment.vertices.empty() || fragment.hap == 0 || fragment.hap > 2) continue;
        Chain chain;
        chain.sample = fragment.sample;
        chain.hap = fragment.hap;
        chain.component = fragment.component;
        chain.source_fragments.push_back(i);
        chain.vertices = fragment.vertices;
        chain.parts.push_back({i, UINT32_MAX, 0, fragment.length, false});
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
            chain.parts.push_back({current, UINT32_MAX, begin_bp, end_bp, false});
        }
        for (uint32_t i = begin; i <= end; ++i) {
            if (!chain.vertices.empty() && chain.vertices.back() == fragment.vertices[i]) continue;
            chain.vertices.push_back(fragment.vertices[i]);
        }
        if (outgoing[current] < 0) break;

        const Candidate& edge = selected[outgoing[current]];
        if (!edge.fill) break;
        const uint32_t gap = static_cast<uint32_t>(chain.gaps.size());
        chain.gaps.push_back(edge);
        chain.parts.push_back({
            UINT32_MAX, gap, edge.fill->left_cut, edge.fill->right_cut, true
        });
        for (size_t i = 1; i + 1 < edge.fill->vertices.size(); ++i) {
            if (!chain.vertices.empty() && chain.vertices.back() == edge.fill->vertices[i]) continue;
            chain.vertices.push_back(edge.fill->vertices[i]);
        }
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
            if (!candidate.fill || candidate.fill->left_cut >= candidate.fill->right_cut) continue;
            if (candidate.fill->sequence.empty()) continue;

            const int64_t gap_begin = fragments_[candidate.left].layout_start +
                static_cast<int64_t>(candidate.left_boundary.target_cut_bp);
            const int64_t gap_end = fragments_[candidate.right].layout_start +
                static_cast<int64_t>(candidate.right_boundary.target_cut_bp);
            const int64_t overlap_begin = std::max(fragment_begin, gap_begin);
            const int64_t overlap_end = std::min(fragment_end, gap_end);
            if (overlap_end <= overlap_begin || static_cast<uint64_t>(overlap_end - overlap_begin) * 10 < fragment.length * 8) continue;

            // ------------------------------------------------ Compare covered sequences ------------------------------------------------
            std::string& gap_sequence = gap_sequences[candidate_id];
            if (gap_sequence.empty()) {
                gap_sequence = candidate.fill->sequence.substr(
                    candidate.fill->left_cut,
                    candidate.fill->right_cut - candidate.fill->left_cut
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
            log_stream() << "  - Covered contig: " << paths_[fragment.path_id].name
                         << " by bridge " << paths_[fragments_[candidate.bridge].path_id].name
                         << " (similarity=" << std::fixed << std::setprecision(4)
                         << similarity << ")\n";
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
    map_options.mid_occ_frac = mm2_.minimizer_freq;
    map_options.min_mid_occ = std::min(map_options.min_mid_occ, mm2_.max_occ - 1);
    map_options.max_mid_occ = mm2_.max_occ;
    map_options.mid_occ = 0;

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
        if (std::max(hap_bp[0], hap_bp[1]) <= params_.dedup_component) {
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
        const bool is_small = std::max(hap_bp[0], hap_bp[1]) <= params_.dedup_component;
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
    default_map_options.mid_occ_frac = mm2_.minimizer_freq;
    default_map_options.min_mid_occ = std::min(default_map_options.min_mid_occ, mm2_.max_occ - 1);
    default_map_options.max_mid_occ = mm2_.max_occ;
    default_map_options.mid_occ = 0;

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
    log_stream() << "  - Node-walk gaps filled: " << selected << "\n";
    log_stream() << "  - Sample paths written: " << path_count << "\n\n";
}

// Select representative haplotypes, assemble sequence and write final GFA and reports.
GfaGapfill::PrimaryPaths GfaGapfill::save_samples(
    const std::string& prefix,
    const std::string& command_line
) {
    PrimaryPaths primary_paths;
    log_stream() << "Writing gap-filled sample graphs ...\n";
    std::unordered_set<std::string> filenames{"primary"};
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
    save_haplotype_("primary", 1, primary.primary[0], prefix + ".primary.hap1", command_line, true, "pr", &primary_paths[0]);
    save_haplotype_("primary", 2, primary.primary[1], prefix + ".primary.hap2", command_line, true, "pr", &primary_paths[1]);
    std::cerr << '\n';

    log_stream() << "Writing low-quality primary graphs ...\n";
    save_haplotype_("primary.lowQ", 1, primary.low_quality[0], prefix + ".primary.hap1.lowQ", command_line, true, "lq");
    save_haplotype_("primary.lowQ", 2, primary.low_quality[1], prefix + ".primary.hap2.lowQ", command_line, true, "lq");

    SAVE report(prefix + ".gapfill.tsv");
    report.save(
        "sample\thap\tcomponent\tcontig\tleft\tright\tbridge\tgap_interval\tphase\n"
    );
    for (const GapRecord& record : records_) {
        std::ostringstream line;
        line << record.sample << '\t' << static_cast<uint32_t>(record.hap) << '\t'
             << record.component << '\t'
             << record.contig << '\t' << record.left << '\t' << record.right << '\t'
             << record.bridge << '\t'
             << record.gap_beg << '-' << record.gap_end << '\t'
             << std::fixed << std::setprecision(2);
        line << metric_pair(record.phase_left, record.phase_right) << '\n';
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
    return primary_paths;
}

void GfaGapfill::save_haplotype_(
    const std::string& label,
    uint8_t hap,
    const std::vector<const Chain*>& chains,
    const std::string& prefix,
    const std::string& command_line,
    bool primary,
    const char* primary_kind,
    std::vector<GfaBubble::ReferencePath>* reference_paths
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
    header += "H\tCO:Z:  SC:Z   target pieces as contig:start-end+/- (source coordinates; 0-based, half-open)\n";
    header += "H\tCO:Z:  GS:Z   bridge-source:start-end+/- pieces in GF order (walk coordinates; 0-based, half-open)\n";
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
                    std::string part;
                    if (chain_part.filled) {
                        const Candidate& gap = chain.gaps[chain_part.gap];
                        if (gap.fill && chain_part.end <= gap.fill->sequence.size()) {
                            part = gap.fill->sequence.substr(
                                chain_part.begin, chain_part.end - chain_part.begin
                            );
                        }
                    } else {
                        part = fragment_subsequence_(
                            fragments_[chain_part.fragment], chain_part.begin, chain_part.end
                        );
                    }
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

        // ------------------------------------------------ Export the exact primary node path ------------------------------------------------
        if (reference_paths) {
            uint64_t path_length = 0;
            for (size_t i = 0; i < assembled_vertices->size(); ++i) {
                const uint32_t vertex = (*assembled_vertices)[i];
                const uint32_t node_length = nodes_[Vertex::get_segment_id(vertex)].length;
                const uint32_t overlap = i == 0 ? 0 : get_edge_ow((*assembled_vertices)[i - 1], vertex);
                path_length += node_length > overlap ? node_length - overlap : 0;
            }
            if (path_length == assembled_length) {
                reference_paths->push_back({name, *assembled_vertices});
            } else {
                warning_stream() << "  ! Primary path length does not match output contig: "
                                 << name << " (path=" << path_length
                                 << ", output=" << assembled_length << ")\n";
            }
        }

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
                     << ':' << begin << '-' << end << (source.reverse ? '-' : '+');
            }
            tags << "\tGS:Z:";
            first = true;
            for (const Chain::Part& part : chain.parts) {
                if (!part.filled) continue;
                if (!first) tags << ',';
                first = false;
                const Candidate& gap = chain.gaps[part.gap];
                const Fragment& bridge = fragments_[gap.bridge];
                tags << bridge.sample << ".hap" << static_cast<uint32_t>(bridge.hap)
                     << ':' << part.begin << '-' << part.end << '+';
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
                candidate.right_boundary.phase.similarity
            });
        }
        GfaGapfillDebugger::chain(
            name, hap, chain.component, chain.gaps.size() + 1,
            intervals.size(), assembled_length
        );
    }

    log_stream() << "  - " << label << ".hap" << static_cast<uint32_t>(hap) << ": " << chains.size() << " contigs\n";
}


// ================================================= Sequence helpers =================================================
std::string GfaGapfill::fragment_sequence_(
    const Fragment& fragment,
    uint32_t begin,
    uint32_t end
) const {
    if (begin >= end || end > fragment.vertices.size()) return {};

    std::string sequence;
    sequence.reserve(fragment.path_bp[end] - fragment.path_bp[begin]);
    for (uint32_t i = begin; i < end; ++i) {
        const uint32_t overlap = i == begin ? 0 :
            get_edge_ow(fragment.vertices[i - 1], fragment.vertices[i]);
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
// Build ordered graph anchors and cache the local evidence once.
bool GfaGapfillBoundary::prepare(
    uint32_t left_id, uint32_t right_id, uint32_t bridge_id,
    Boundary& left_boundary, Boundary& right_boundary
) const {
    if (!find_anchors_(left_id, right_id, bridge_id, left_boundary, right_boundary)) return false;
    const GfaGapfill::Fragment& left = gapfill_.fragments_[left_id];
    const GfaGapfill::Fragment& right = gapfill_.fragments_[right_id];
    const GfaGapfill::Fragment& bridge = gapfill_.fragments_[bridge_id];

    // Shared-node cuts are the final node-level gap boundaries. A three-way
    // crossover owns its shared node on the left and starts the right target
    // immediately after it, so the node appears exactly once.
    left_boundary.target_cut_bp = left.path_bp[left_boundary.target_pos + 1];
    left_boundary.bridge_cut_bp = bridge.path_bp[left_boundary.bridge_pos + 1];
    if (left_boundary.bridge_pos == right_boundary.bridge_pos) {
        right_boundary.bridge_cut_bp = left_boundary.bridge_cut_bp;
        right_boundary.target_cut_bp = right.path_bp[right_boundary.target_pos + 1];
    } else {
        right_boundary.bridge_cut_bp = bridge.path_bp[right_boundary.bridge_pos];
        right_boundary.target_cut_bp = right.path_bp[right_boundary.target_pos];
    }

    // Cache local phase evidence once per candidate.
    left_boundary.phase = phase_similarity_(
        left_id, left_boundary.target_pos,
        bridge_id, left_boundary.bridge_pos, true
    );
    right_boundary.phase = phase_similarity_(
        right_id, right_boundary.target_pos,
        bridge_id, right_boundary.bridge_pos, false
    );
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
    bool before
) const {
    Boundary boundary;
    boundary.target_pos = target_pos;
    boundary.bridge_pos = bridge_pos;
    boundary.phase = phase_similarity_(
        target, target_pos, bridge, bridge_pos, before
    );
    return boundary;
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
        if (shorter_bp == 0 || static_cast<double>(direct_bp) / shorter_bp > gapfill_.params_.max_overlap) {
            return false;
        }
        if (2 * left_sum + direct_matches <= direct_matches * left.vertices.size() || 2 * right_sum + direct_matches >= direct_matches * right.vertices.size()) {
            return false;
        }
    }

    // ------------------------------------------------ Collect shared nodes ------------------------------------------------
    std::vector<std::pair<uint32_t, uint32_t>> left_matches;  // {bridge_pos, left_pos}
    std::vector<std::pair<uint32_t, uint32_t>> right_matches;  // {bridge_pos, right_pos}
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

    // ------------------------------------------------ Find an ordered source/sink pair ------------------------------------------------
    // Without a three-way crossover, retain as much target sequence as
    // possible while keeping source before sink on the spanning bridge.
    GfaGapfill::Candidate ordered;
    bool have_ordered = false;
    uint64_t ordered_retained = 0;
    uint32_t ordered_span = UINT32_MAX;
    std::vector<uint32_t> suffix(right_matches.size());
    suffix.back() = static_cast<uint32_t>(right_matches.size() - 1);
    for (size_t i = right_matches.size() - 1; i > 0; --i) {
        const uint32_t current = static_cast<uint32_t>(i - 1);
        const uint32_t best = suffix[i];
        suffix[current] = right_matches[current].second < right_matches[best].second ? current : best;
    }
    for (const auto& left_match : left_matches) {
        const auto it = std::upper_bound(
            right_matches.begin(), right_matches.end(),
            std::pair<uint32_t, uint32_t>{left_match.first, UINT32_MAX}
        );
        if (it == right_matches.end()) continue;

        const auto& right_match = right_matches[suffix[it - right_matches.begin()]];
        const uint64_t retained = left.path_bp[left_match.second + 1] +
            right.length - right.path_bp[right_match.second];
        const uint32_t span = right_match.first - left_match.first;
        if (!have_ordered || retained > ordered_retained ||
            (retained == ordered_retained && span < ordered_span) ||
            (retained == ordered_retained && span == ordered_span &&
             std::tie(left_match.first, right_match.first) <
             std::tie(ordered.left_boundary.bridge_pos, ordered.right_boundary.bridge_pos))) {
            ordered.left_boundary.target_pos = left_match.second;
            ordered.left_boundary.bridge_pos = left_match.first;
            ordered.right_boundary.bridge_pos = right_match.first;
            ordered.right_boundary.target_pos = right_match.second;
            ordered_retained = retained;
            ordered_span = span;
            have_ordered = true;
        }
    }

    // ------------------------------------------------ Find a three-way crossover ------------------------------------------------
    GfaGapfill::Candidate crossover;
    bool have_crossover = false;
    uint64_t crossover_retained = 0;
    for (const auto& match : left_matches) {
        const uint32_t right_pos = find_position_(right_id, bridge.vertices[match.first]);
        if (right_pos == UINT32_MAX) continue;

        const uint64_t retained = left.path_bp[match.second + 1] +
            right.length - right.path_bp[right_pos + 1];
        if (!have_crossover || retained > crossover_retained ||
            (retained == crossover_retained && match.first < crossover.left_boundary.bridge_pos)) {
            crossover.left_boundary.target_pos = match.second;
            crossover.left_boundary.bridge_pos = match.first;
            crossover.right_boundary.bridge_pos = match.first;
            crossover.right_boundary.target_pos = right_pos;
            crossover_retained = retained;
            have_crossover = true;
        }
    }

    // ------------------------------------------------ Select and validate the boundary ------------------------------------------------
    if (have_crossover) {
        candidate = crossover;
    } else if (have_ordered) {
        candidate = ordered;
    } else {
        return false;
    }

    if (!is_cycle_free_(left_id, right_id, bridge_id, candidate.left_boundary, candidate.right_boundary)) {
        const bool selected_crossover = candidate.left_boundary.bridge_pos == candidate.right_boundary.bridge_pos;
        const GfaGapfill::Candidate* fallback = selected_crossover && have_ordered ? &ordered :
            (!selected_crossover && have_crossover ? &crossover : nullptr);
        if (!fallback || !is_cycle_free_(left_id, right_id, bridge_id, fallback->left_boundary, fallback->right_boundary)) {
            return false;
        }
        candidate = *fallback;
    }

    left_boundary = candidate.left_boundary;
    right_boundary = candidate.right_boundary;
    return true;
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
    const bool crossover = left_boundary.bridge_pos == right_boundary.bridge_pos;
    const uint32_t right_begin = right_boundary.target_pos + static_cast<uint32_t>(crossover);
    for (uint32_t i = right_begin; i < right.vertices.size(); ++i) {
        if (!seen.insert(Vertex::get_segment_id(right.vertices[i])).second) return false;
    }
    return true;
}

std::pair<uint32_t, uint32_t> GfaGapfillBoundary::phase_region_(
    uint32_t fragment_id,
    uint32_t pos,
    bool before
) const {
    const GfaGapfill::Fragment& fragment = gapfill_.fragments_[fragment_id];
    const uint64_t window = gapfill_.params_.phase_win;
    uint64_t begin_bp = 0;
    uint64_t end_bp = 0;
    if (before) {
        end_bp = fragment.path_bp[pos + 1];
        begin_bp = end_bp > window ? end_bp - window : 0;
    } else {
        begin_bp = fragment.path_bp[pos];
        end_bp = std::min(fragment.length, begin_bp + window);
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
    const auto target_region = phase_region_(target_id, target_pos, before);
    const auto bridge_region = phase_region_(bridge_id, bridge_pos, before);
    const uint64_t target_bp = target.bubble_bp[target_region.second] - target.bubble_bp[target_region.first];
    const uint64_t bridge_bp = bridge.bubble_bp[bridge_region.second] - bridge.bubble_bp[bridge_region.first];
    support.target_bp = target_bp;
    support.bridge_bp = bridge_bp;
    if (target_bp == 0 || bridge_bp == 0) return support;

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
