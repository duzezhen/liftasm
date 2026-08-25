#include "../include/gfa_clean.hpp"
#include "../include/progress_tracker.hpp"

#include <algorithm>
#include <sstream>
#include <utility>

GfaCleaner::GfaCleaner(CleanOpts params) : params_(std::move(params)) {}

void GfaCleaner::build_ctg_read_index_(
    const std::vector<std::string>& ctg_files,
    ReadIndex& read_index
) {
    GfaGraph ctg_graph;
    ctg_graph.load_from_GFA(ctg_files);
    ctg_graph.build_nodes_connectivity_index();

    log_stream() << "Building contig and read connectivity ...\n";

    const std::vector<GfaAlignment>& alignments = ctg_graph.getAlignments();
    read_placements_.clear();
    read_placements_.reserve(alignments.size());
    read_index.clear();
    read_index.reserve(alignments.size());

    uint32_t next_read_id = 0;
    ProgressTracker progress(alignments.size());

    for (const GfaAlignment& alignment : alignments) {
        progress.hit();

        const uint32_t contig_id = static_cast<uint32_t>(alignment.unitig_node_id);
        const uint32_t component = ctg_graph.node_connectivity_id(Vertex(contig_id, false));
        if (component == UINT32_MAX) continue;

        auto [it, inserted] = read_index.try_emplace(alignment.read_name);
        if (inserted) it->second.id = next_read_id++;
        it->second.components.push_back(component);

        read_placements_.push_back({
            component,
            contig_id,
            ctg_graph.getNodeName(contig_id),
            alignment.position_on_unitig,
            alignment.read_is_rev,
            alignment.read_name,
            alignment.read_start_pos,
            alignment.read_end_pos,
            alignment.optional_tags
        });
    }
    progress.finish();

    for (auto& [read, info] : read_index) {
        (void)read;
        std::sort(info.components.begin(), info.components.end());
        info.components.erase(
            std::unique(info.components.begin(), info.components.end()),
            info.components.end()
        );
    }

    std::sort(read_placements_.begin(), read_placements_.end(),
        [](const ReadPlacement& a, const ReadPlacement& b) {
            if (a.component != b.component) return a.component < b.component;
            if (a.contig_id != b.contig_id) return a.contig_id < b.contig_id;
            if (a.contig_position != b.contig_position) return a.contig_position < b.contig_position;
            if (a.read_start != b.read_start) return a.read_start < b.read_start;
            return a.read_name < b.read_name;
        }
    );

    uint32_t component_count = 0;
    for (uint32_t component : ctg_graph.get_nodes_connectivity_index()) {
        if (component != UINT32_MAX) component_count = std::max(component_count, component + 1);
    }

    log_stream() << "  - Contig components: " << component_count << "\n";
    log_stream() << "  - Read placements: " << read_placements_.size() << "\n";
    log_stream() << "  - Unique reads: " << read_index.size() << "\n\n";
}

void GfaCleaner::save_read_placements(const std::string& output_file) const {
    log_stream() << "Writing ordered read placements to " << output_file << " ...\n";

    SAVE saver(output_file);
    saver.save("component\tcontig\tcontig_pos\tstrand\tread\tread_start\tread_end\toptional\n");

    std::ostringstream out;
    for (const ReadPlacement& placement : read_placements_) {
        out << placement.component << '\t'
            << placement.contig_name << '\t'
            << placement.contig_position << '\t'
            << (placement.is_reverse ? '-' : '+') << '\t'
            << placement.read_name << '\t'
            << placement.read_start << '\t'
            << placement.read_end << '\t'
            << placement.optional_tags << '\n';
        saver.save(out.str());
        out.str("");
        out.clear();
    }

    log_stream() << "  - Read placements written: " << read_placements_.size() << "\n\n";
}

std::vector<GfaCleaner::NodeEvidence> GfaCleaner::build_utg_evidence_(
    const ReadIndex& read_index
) const {
    struct Assignment {
        uint32_t node;
        uint32_t label;
        uint32_t read;
    };

    log_stream() << "Projecting contig connectivity onto unitigs ...\n";

    std::vector<Assignment> assignments;
    assignments.reserve(alignments_.size());

    for (const GfaAlignment& alignment : alignments_) {
        const auto it = read_index.find(alignment.read_name);
        if (it == read_index.end() || alignment.unitig_node_id >= nodes_.size()) continue;

        for (uint32_t component : it->second.components) {
            assignments.push_back({
                static_cast<uint32_t>(alignment.unitig_node_id),
                component,
                it->second.id
            });
        }
    }

    std::sort(assignments.begin(), assignments.end(), [](const Assignment& a, const Assignment& b) {
        if (a.node != b.node) return a.node < b.node;
        if (a.label != b.label) return a.label < b.label;
        return a.read < b.read;
    });
    assignments.erase(
        std::unique(assignments.begin(), assignments.end(), [](const Assignment& a, const Assignment& b) {
            return a.node == b.node && a.label == b.label && a.read == b.read;
        }),
        assignments.end()
    );

    std::vector<NodeEvidence> evidence(nodes_.size());
    for (size_t i = 0; i < assignments.size();) {
        NodeEvidence& node = evidence[assignments[i].node];
        const uint32_t node_id = assignments[i].node;

        while (i < assignments.size() && assignments[i].node == node_id) {
            const uint32_t label = assignments[i].label;
            uint32_t count = 0;
            while (i < assignments.size() &&
                   assignments[i].node == node_id &&
                   assignments[i].label == label) {
                ++count;
                ++i;
            }

            node.labels.push_back({label, count});
            node.total += count;
            if (count > node.dominant_count) {
                node.dominant_label = label;
                node.dominant_count = count;
            }
        }
    }

    size_t labeled_nodes = 0;
    for (const NodeEvidence& node : evidence) labeled_nodes += node.total != 0;

    log_stream() << "  - Matched read-component assignments: " << assignments.size() << "\n";
    log_stream() << "  - Unitigs labeled: " << labeled_nodes << " / " << nodes_.size() << "\n\n";
    return evidence;
}

bool GfaCleaner::resolved_different_(const NodeEvidence& a, const NodeEvidence& b) const {
    if (a.total < params_.min_reads || b.total < params_.min_reads) return false;

    const double purity_a = static_cast<double>(a.dominant_count) / a.total;
    const double purity_b = static_cast<double>(b.dominant_count) / b.total;
    if (purity_a < params_.min_purity || purity_b < params_.min_purity) return false;
    return a.dominant_label != b.dominant_label;
}

uint64_t GfaCleaner::component_pair_key_(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | b;
}

GfaCleaner::ComponentOverlapIndex GfaCleaner::build_component_overlap_index_(
    const std::vector<NodeEvidence>& evidence,
    const std::unordered_set<uint64_t>& candidate_pairs
) const {
    ComponentOverlapIndex index;
    index.node_counts.reserve(candidate_pairs.size() * 2 + 1);
    index.shared_nodes.reserve(candidate_pairs.size());

    std::unordered_map<uint32_t, std::vector<uint32_t>> higher_neighbors;
    higher_neighbors.reserve(candidate_pairs.size() * 2 + 1);
    for (uint64_t key : candidate_pairs) {
        higher_neighbors[static_cast<uint32_t>(key >> 32)].push_back(static_cast<uint32_t>(key));
    }
    for (auto& [component, neighbors] : higher_neighbors) {
        (void)component;
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    for (const NodeEvidence& node : evidence) {
        for (const LabelCount& label : node.labels) ++index.node_counts[label.label];

        for (const LabelCount& label : node.labels) {
            const auto neighbor_it = higher_neighbors.find(label.label);
            if (neighbor_it == higher_neighbors.end()) continue;

            for (uint32_t other : neighbor_it->second) {
                const auto found = std::lower_bound(
                    node.labels.begin(), node.labels.end(), other,
                    [](const LabelCount& item, uint32_t value) { return item.label < value; }
                );
                if (found != node.labels.end() && found->label == other) {
                    ++index.shared_nodes[component_pair_key_(label.label, other)];
                }
            }
        }
    }

    return index;
}

double GfaCleaner::component_overlap_(
    const ComponentOverlapIndex& index,
    uint32_t a,
    uint32_t b
) {
    const auto count_a = index.node_counts.find(a);
    const auto count_b = index.node_counts.find(b);
    if (count_a == index.node_counts.end() || count_b == index.node_counts.end()) return 0.0;

    const uint32_t denominator = std::min(count_a->second, count_b->second);
    if (denominator == 0) return 0.0;

    const auto shared = index.shared_nodes.find(component_pair_key_(a, b));
    const uint32_t shared_count = shared == index.shared_nodes.end() ? 0 : shared->second;
    return static_cast<double>(shared_count) / denominator;
}

size_t GfaCleaner::clean(const std::vector<std::string>& ctg_files) {
    ReadIndex read_index;
    build_ctg_read_index_(ctg_files, read_index);
    std::vector<NodeEvidence> evidence = build_utg_evidence_(read_index);

    log_stream() << "Removing weak component-conflicting links ...\n";

    std::unordered_set<uint64_t> candidate_pairs;
    candidate_pairs.reserve(arcs_.size() / 16 + 1);

    for (const GfaArc& arc : arcs_) {
        if (arc.get_del() || arc.get_comp()) continue;
        const uint32_t source = arc.get_source_segment_id();
        const uint32_t target = arc.get_target_segment_id();
        if (source == target || source >= evidence.size() || target >= evidence.size()) continue;
        if (!resolved_different_(evidence[source], evidence[target])) continue;

        candidate_pairs.insert(component_pair_key_(
            evidence[source].dominant_label,
            evidence[target].dominant_label
        ));
    }

    const ComponentOverlapIndex overlap_index = build_component_overlap_index_(evidence, candidate_pairs);

    std::unordered_set<uint64_t> weak_pairs;
    weak_pairs.reserve(candidate_pairs.size());
    for (uint64_t key : candidate_pairs) {
        const uint32_t a = static_cast<uint32_t>(key >> 32);
        const uint32_t b = static_cast<uint32_t>(key);
        if (component_overlap_(overlap_index, a, b) < params_.min_comp_overlap) {
            weak_pairs.insert(key);
        }
    }

    std::unordered_set<uint64_t> rejected_links;
    rejected_links.reserve(arcs_.size() / 32 + 1);

    ProgressTracker progress(arcs_.size());
    for (const GfaArc& arc : arcs_) {
        progress.hit();
        if (arc.get_del() || arc.get_comp()) continue;

        const uint32_t source = arc.get_source_segment_id();
        const uint32_t target = arc.get_target_segment_id();
        if (source == target || source >= evidence.size() || target >= evidence.size()) continue;
        if (!resolved_different_(evidence[source], evidence[target])) continue;

        const uint32_t label_a = evidence[source].dominant_label;
        const uint32_t label_b = evidence[target].dominant_label;
        if (weak_pairs.count(component_pair_key_(label_a, label_b))) {
            rejected_links.insert(arc.get_link_id());
        }
    }
    progress.finish();

    if (!rejected_links.empty()) {
        for (GfaArc& arc : arcs_) {
            if (rejected_links.count(arc.get_link_id())) arc.set_del(true);
        }
        cleanup();
        countLinks();
        calculateDeg();
    }

    log_stream() << "  - Component pairs evaluated: " << candidate_pairs.size() << "\n";
    log_stream() << "  - Low-overlap component pairs: " << weak_pairs.size() << "\n";
    log_stream() << "  - Links removed: " << rejected_links.size() << "\n";
    log_stream() << "  - Links retained: " << total_uniq_links_ << "\n\n";
    return rejected_links.size();
}
