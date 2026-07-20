#include "../include/gfa_split.hpp"

#include <algorithm>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "../include/logger.hpp"
#include "../include/save.hpp"

bool GfaSplitter::save_component_to_disk_(
    const std::string& filename,
    uint32_t component_id,
    bool write_paths,
    bool write_align,
    bool write_seq
) const {
    auto in_comp = [&](uint32_t sid) -> bool {
        if (sid >= nodes_.size() || nodes_[sid].deleted) return false;
        return node_connectivity_id(Vertex(sid, false)) == component_id;
    };

    SAVE saver(filename);
    std::ostringstream oss;

    bool wrote_node = false;

    // 1. S-lines
    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        if (!in_comp(sid)) continue;

        const auto& n = nodes_[sid];

        oss << "S\t" << n.name << '\t';
        if (n.sequence.empty() || n.sequence == "*" || !write_seq) oss << '*';
        else oss << n.sequence;

        oss << "\tLN:i:" << n.length;

        const std::string sn = sample_ids_to_csv_(n.sample_ids);
        if (!sn.empty()) oss << "\tSN:Z:" << sn;

        oss << GfaAuxParser::write_aux_as_text(n.aux.aux_data);
        oss << '\n';

        saver.save(oss.str());
        oss.str("");
        oss.clear();

        wrote_node = true;
    }

    if (!wrote_node) return false;

    // 2. L-lines
    for (const auto& a : arcs_) {
        if (a.get_del() || a.get_comp()) continue;

        const uint32_t v = a.get_source_vertex_id();
        const uint32_t w = a.get_target_vertex_id();

        const uint32_t vseg = NodeHandle::get_segment_id(v);
        const uint32_t wseg = NodeHandle::get_segment_id(w);

        if (!in_comp(vseg) || !in_comp(wseg)) continue;

        const bool vrev = NodeHandle::get_is_reverse(v);
        const bool wrev = NodeHandle::get_is_reverse(w);

        oss << "L\t" << nodes_[vseg].name << '\t' << (vrev ? '-' : '+') << '\t'
            << nodes_[wseg].name << '\t' << (wrev ? '-' : '+') << '\t'
            << format_overlap_field(a.ov, a.ow);

        const uint64_t lid = a.get_link_id();
        if (lid < link_aux_.size()) {
            oss << GfaAuxParser::write_aux_as_text(link_aux_[lid].aux_data);
        }

        oss << '\n';

        saver.save(oss.str());
        oss.str("");
        oss.clear();
    }

    // 3. P-lines
    if (write_paths) {
        for (const auto& p : paths_) {
            if (p.segments.empty()) continue;

            bool ok = true;
            for (const auto& s : p.segments) {
                if (!in_comp(s.node_id)) {
                    ok = false;
                    break;
                }
            }
            if (!ok) continue;

            oss << "P\t" << p.name << '\t';

            for (size_t i = 0; i < p.segments.size(); ++i) {
                const auto& s = p.segments[i];
                if (i) oss << ',';
                oss << nodes_[s.node_id].name << (s.is_reverse ? '-' : '+');
            }

            oss << "\t*\n";

            saver.save(oss.str());
            oss.str("");
            oss.clear();
        }
    }

    // 4. A-lines
    if (write_align) {
        for (const auto& a : alignments_) {
            if (!in_comp(a.unitig_node_id)) continue;

            oss << "A\t"
                << nodes_[a.unitig_node_id].name << '\t'
                << a.position_on_unitig << '\t'
                << (a.read_is_rev ? '-' : '+') << '\t'
                << a.read_name << '\t'
                << a.read_start_pos << '\t'
                << a.read_end_pos;

            if (!a.optional_tags.empty()) oss << '\t' << a.optional_tags;
            oss << '\n';

            saver.save(oss.str());
            oss.str("");
            oss.clear();
        }
    }

    return true;
}

uint32_t GfaSplitter::split_by_components(
    const std::string& prefix,
    bool write_paths,
    bool write_align,
    bool skip_comp
) {
    build_nodes_connectivity_index(skip_comp);

    log_stream() << "Splitting graph into components and writing to disk ...\n";

    std::unordered_set<uint32_t> comps;
    comps.reserve(nodes_.size());

    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        if (nodes_[sid].deleted) continue;

        const uint32_t cid = node_connectivity_id(Vertex(sid, false));
        if (cid != UINT32_MAX) comps.insert(cid);
    }

    std::vector<uint32_t> comp_ids(comps.begin(), comps.end());
    std::sort(comp_ids.begin(), comp_ids.end());

    uint32_t written = 0;

    for (uint32_t cid : comp_ids) {
        const std::string out_gfa = prefix + ".comp" + std::to_string(cid) + ".gfa";
        const std::string out_noseq = prefix + ".comp" + std::to_string(cid) + ".noseq.gfa";

        const bool ok = save_component_to_disk_(
            out_gfa,
            cid,
            write_paths,
            write_align,
            true
        );

        if (ok) {
            save_component_to_disk_(
                out_noseq,
                cid,
                write_paths,
                write_align,
                false
            );
            ++written;
        }
    }

    log_stream() << "  - Components written: " << written << " / " << comp_ids.size() << "\n\n";

    return written;
}