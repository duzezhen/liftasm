#include "../include/gfa_parser.hpp"
#include "../include/gfa_walker.hpp"
#include "../include/gfa_name.hpp"
#include "../include/gfa_name_check.hpp"
#include "../include/progress_tracker.hpp"

#include <iomanip>
#include <cstring>
#include <limits>
#include <iterator>
#include <stack>
#include <queue>
#include <sstream>
#include <unordered_set>
#include <algorithm>

using namespace wfa; // WFAlignerGapAffine


// ------------------------------------------------ Helper functions ------------------------------------------------
static std::string gfa_sample_name_from_path_(const std::string& path) {
    size_t p = path.find_last_of("/\\");
    std::string s = (p == std::string::npos) ? path : path.substr(p + 1);

    auto strip = [&](const std::string& suf) {
        if (s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0) {
            s.resize(s.size() - suf.size());
            return true;
        }
        return false;
    };

    if (!strip(".gfa.gz")) {
        strip(".gfa");
        strip(".gz");
    }

    // If the file is a collapse or deoverlap graph, don't add anything tag to the nodes.
    if (s.find(".collapse") != std::string::npos ||
        s.find(".deoverlap") != std::string::npos) {
        return "";
    }

    static const char* hifiasm_suffixes[] = {
        ".asm.bp.hap1.p_ctg.noseq",
        ".asm.bp.hap2.p_ctg.noseq",
        ".asm.bp.hap1.p_ctg",
        ".asm.bp.hap2.p_ctg",

        ".asm.bp.p_ctg.noseq",
        ".asm.bp.p_ctg",

        ".asm.bp.p_utg.noseq",
        ".asm.bp.p_utg",

        ".asm.bp.r_utg.noseq",
        ".asm.bp.r_utg",

        ".asm.bp.hap1.p_ctg.lowQ",
        ".asm.bp.hap2.p_ctg.lowQ",
        ".asm.bp.p_ctg.lowQ",
        ".asm.bp.p_utg.lowQ",
        ".asm.bp.r_utg.lowQ"
    };

    for (auto suf : hifiasm_suffixes) {
        if (strip(suf)) break;
    }

    return s.empty() ? "gfa" : s;
}

bool GfaGraph::input_segment_names_need_namespace_(
    const std::vector<std::string>& filenames,
    const std::vector<uint32_t>& duplicate_groups
) {
    const bool grouped = !duplicate_groups.empty();
    std::unordered_map<std::string, uint32_t> first_file;
    std::unordered_set<std::string> grouped_names;
    bool found_cross_file_duplicate = false;

    for (uint32_t file_id = 0; file_id < filenames.size(); ++file_id) {
        std::unordered_set<std::string> local_names;
        GzChunkReader reader(filenames[file_id]);
        std::string line;
        while (reader.readLine(line)) {
            if (line.empty() || line[0] != 'S') continue;
            std::stringstream in(line);
            char type = 0;
            std::string name;
            in >> type >> name;
            if (type != 'S' || name.empty()) continue;

            if (!local_names.insert(name).second) {
                error_stream() << "Duplicate segment name '" << name << "' within GFA file '" << filenames[file_id] << "'\n";
                std::exit(1);
            }

            if (grouped) {
                const std::string key = std::to_string(duplicate_groups[file_id]) + '\t' + name;
                if (!grouped_names.insert(key).second) {
                    error_stream() << "Duplicate segment name '" << name << "' within contig sample " << (duplicate_groups[file_id] + 1) << "\n";
                    std::exit(1);
                }
            }

            const auto inserted = first_file.emplace(name, file_id);
            if (!inserted.second && inserted.first->second != file_id) {
                found_cross_file_duplicate = true;
            }
        }
        reader.close();
    }
    return found_cross_file_duplicate;
}


GfaGraph::GfaGraph()  = default;
GfaGraph::~GfaGraph() = default;

// Load GFA file(s) into the graph structure
void GfaGraph::load_from_GFA(
    const std::vector<std::string>& filenames,
    const std::vector<std::string>& sample_names,
    bool read_links,
    const std::vector<uint32_t>& duplicate_groups
) {
    if (filenames.empty()) {
        error_stream() << "No GFA input file provided.\n";
        std::exit(1);
    }

    if (!sample_names.empty() && sample_names.size() != filenames.size()) {
        error_stream() << "-n/--name must have the same number of values as -g/--gfa\n";
        std::exit(1);
    }
    if (!duplicate_groups.empty() && duplicate_groups.size() != filenames.size()) {
        error_stream() << "Internal error: duplicate-name groups must match the number of GFA files\n";
        std::exit(1);
    }

    const bool has_cross_input_duplicate = filenames.size() > 1 && input_segment_names_need_namespace_(filenames, duplicate_groups);
    const bool namespace_all = has_cross_input_duplicate;
    if (has_cross_input_duplicate) {
        warning_stream() << "Duplicate segment names found across input " << (duplicate_groups.empty() ? "files" : "samples") << "; adding a namespace suffix to every segment name\n\n";
    }

    std::vector<gfa_name_check::FileRenamePlan> plans;
    plans.reserve(filenames.size());

    // pass-1: register S names file by file, and check for duplicates and invalid names
    for (size_t file_id = 0; file_id < filenames.size(); ++file_id) {
        const auto& filename = filenames[file_id];
        const uint32_t namespace_id = duplicate_groups.empty() ? file_id : duplicate_groups[file_id];
        const std::string suffix = namespace_all ? "_" + std::to_string(namespace_id + 1) : "";
        gfa_name_check::FileRenamePlan plan(/*max_log_examples=*/5, suffix);

        GzChunkReader zr(filename);
        std::string line;
        while (zr.readLine(line)) {
            if (line.empty() || line[0] != 'S') continue;
            if (line.back() == '\n' || line.back() == '\r') line.pop_back();

            std::stringstream ss(line);
            char c;
            std::string raw_name;
            ss >> c >> raw_name;
            if (raw_name.empty()) continue;

            const std::string fixed_name = plan.register_segment_name(raw_name, name_to_id_map_, filename);

            if (name_to_id_map_.insert({fixed_name, total_segments_}).second) {
                ++total_segments_;
            }
        }
        zr.close();

        plans.emplace_back(std::move(plan));
    }

    nodes_.resize(total_segments_);

    // pass-2: rewrite names with the corresponding file plan, then parse lines to build the graph
    for (std::size_t i = 0; i < filenames.size(); ++i) {
        const auto& filename = filenames[i];
        const auto& plan = plans[i];

        const uint32_t source_id = static_cast<uint32_t>(i);
        const std::string sample_name = sample_names.empty() ? gfa_sample_name_from_path_(filename) : sample_names[i];

        log_stream() << "Loading genome graph from '" << filename << "' (SN Tag: " << (sample_name.empty() ? "N/A" : sample_name) << ") ...\n";

        GzChunkReader zr(filename);
        std::string line;
        while (zr.readLine(line)) {
            if (line.empty()) continue;
            if (line.back() == '\n' || line.back() == '\r') line.pop_back();

            line = plan.rewrite_gfa_line(line);

            std::stringstream ss(line);
            char type;
            ss >> type;
            ss.clear();
            ss.str(line);

            switch (type) {
                case 'S': parseSLine(ss, sample_name); break;
                case 'L': if (read_links) parseLLine(ss); break;
                case 'P': parsePLine(ss); break;
                case 'A': parseALine(ss); break;
                default:  break;
            }
        }
        zr.close();

        plan.print_log(filename);
    }

    if (forbid_overlap_ && has_overlap_) {
        error_stream() << "Overlap-style GFA is not supported. Please check the input file.\n";
        std::exit(1);
    }

    // Finalize the graph after loading
    finalize_();

    // print graph stats
    print_graph_stats();

    return;
}


uint32_t GfaGraph::get_or_add_segment(const std::string& name, bool* is_new) {
    auto it = name_to_id_map_.find(name);
    if (it != name_to_id_map_.end()) {
        if (is_new) *is_new = false;
        return it->second;
    }

    uint32_t id = static_cast<uint32_t>(total_segments_);
    name_to_id_map_[name] = id;
    nodes_.resize(id + 1);
    nodes_[id].name = name;
    ++total_segments_;

    if (is_new) *is_new = true;
    return id;
}

std::string GfaGraph::slice_seq_or_star_(
    const std::vector<GfaNode>& nodes, uint32_t seg_id, uint32_t beg, uint32_t end, bool is_rev
) {
    const auto& n = nodes[seg_id];
    if (n.sequence.empty() || n.sequence == "*" || end > n.sequence.size() || beg > end) {
        return std::string("*");
    }
    std::string sub = n.sequence.substr(beg, end - beg);
    if (is_rev) sub = seqUtils::revcomp(sub);
    return sub;
}

std::string GfaGraph::get_oriented_sequence(Vertex vertex, uint32_t ow) const {
    uint32_t seg_id = vertex.segment_id();
    bool rev = vertex.is_reverse();
    const std::string& seq = nodes_[seg_id].sequence;
    if (seq == "*") return std::string();
    if (ow >= seq.size()) return std::string();
    if (rev) {
        std::string rc = seqUtils::revcomp(seq);
        return rc.substr(ow);
    }
    return seq.substr(ow);
}

std::string GfaGraph::get_path_sequence(const std::vector<uint32_t>& path_vtx_ids) const {
    if (path_vtx_ids.empty()) return "*";

    size_t total_len = 0;
    for (uint32_t vtx_id : path_vtx_ids) {
        uint32_t seg_id = NodeHandle::get_segment_id(vtx_id);
        total_len += getNodeLength(seg_id);
    }

    std::string out;
    out.reserve(total_len);

    bool first = true;
    uint32_t prev_vtx_id = 0;

    for (uint32_t vtx_id : path_vtx_ids) {
        std::string s = get_oriented_sequence(Vertex(vtx_id));
        if (s.empty()) return "*";

        if (first) {
            out += s;
            first = false;
        } else {
            uint32_t ovlp = 0;
            auto arcs = getArcsFromVertex(prev_vtx_id);
            for (const auto* arc : arcs) {
                if (arc && arc->w == vtx_id) {
                    ovlp = arc->ow > 0 ? static_cast<uint32_t>(arc->ow) : 0u;
                    break;
                }
            }
            if (ovlp >= s.size()) out += "";
            else out += s.substr(ovlp);
        }

        prev_vtx_id = vtx_id;
    }

    return out;
}

std::string GfaGraph::get_path_sequence(const std::vector<Vertex>& path_vtxs) const {
    if (path_vtxs.empty()) return "*";

    size_t total_len = 0;
    for (const Vertex& vtx : path_vtxs) {
        uint32_t seg_id = vtx.segment_id();
        total_len += getNodeLength(seg_id);
    }

    std::string out;
    out.reserve(total_len);

    bool first = true;
    Vertex prev_vtx = Vertex(UINT32_MAX);

    for (const Vertex& vtx : path_vtxs) {
        std::string s = get_oriented_sequence(vtx);
        if (s.empty()) return "*";

        if (first) {
            out += s;
            first = false;
        } else {
            uint32_t ovlp = 0;
            auto arcs = getArcsFromVertex(prev_vtx.vertex_id());
            for (const auto* arc : arcs) {
                if (arc && arc->w == vtx.vertex_id()) {
                    ovlp = arc->ow > 0 ? static_cast<uint32_t>(arc->ow) : 0u;
                    break;
                }
            }
            if (ovlp >= s.size()) out += "";
            else out += s.substr(ovlp);
        }

        prev_vtx = vtx;
    }

    return out;
}

std::vector<uint32_t> GfaGraph::get_path_overlaps(const std::vector<uint32_t>& path_vtx_ids) const {
    std::vector<uint32_t> overlaps;
    if (path_vtx_ids.size() <= 1) return overlaps;

    overlaps.reserve(path_vtx_ids.size() - 1);

    for (size_t i = 1; i < path_vtx_ids.size(); ++i) {
        const uint32_t prev_v = path_vtx_ids[i - 1];
        const uint32_t curr_v = path_vtx_ids[i];

        uint32_t ovlp = 0;
        const auto& arcs = getArcsFromVertex(prev_v);
        for (const auto* arc : arcs) {
            if (arc && !arc->get_del() && arc->w == curr_v) {
                ovlp = (arc->ow > 0 ? static_cast<uint32_t>(arc->ow) : 0u);
                break;
            }
        }

        overlaps.push_back(ovlp);
    }

    return overlaps;
}

std::vector<uint32_t> GfaGraph::get_path_overlaps(const std::vector<Vertex>& path_vtxs) const {
    std::vector<uint32_t> overlaps;
    if (path_vtxs.size() <= 1) return overlaps;

    overlaps.reserve(path_vtxs.size() - 1);

    for (size_t i = 1; i < path_vtxs.size(); ++i) {
        const uint32_t prev_v = path_vtxs[i - 1].vertex_id();
        const uint32_t curr_v = path_vtxs[i].vertex_id();

        uint32_t ovlp = 0;
        const auto& arcs = getArcsFromVertex(prev_v);
        for (const auto* arc : arcs) {
            if (arc && !arc->get_del() && arc->w == curr_v) {
                ovlp = (arc->ow > 0 ? static_cast<uint32_t>(arc->ow) : 0u);
                break;
            }
        }

        overlaps.push_back(ovlp);
    }

    return overlaps;
}

uint32_t GfaGraph::get_edge_ow(uint32_t from_vtx_id, uint32_t to_vtx_id) const {
    const auto& arcs = getArcsFromVertex(from_vtx_id);
    for (const auto* arc : arcs) {
        if (!arc || arc->get_del()) continue;
        if (arc->w != to_vtx_id) continue;
        return (arc->ow > 0 ? static_cast<uint32_t>(arc->ow) : 0u);
    }
    return 0u;
}
uint32_t GfaGraph::get_edge_ow(Vertex from_vtx, Vertex to_vtx) const {
    const auto& arcs = getArcsFromVertex(from_vtx.vertex_id());
    for (const auto* arc : arcs) {
        if (!arc || arc->get_del()) continue;
        if (arc->w != to_vtx.vertex_id()) continue;
        return (arc->ow > 0 ? static_cast<uint32_t>(arc->ow) : 0u);
    }
    return 0u;
}


GfaArc* GfaGraph::add_arc(
    uint32_t v, uint32_t w, int32_t ov, int32_t ow,
    int64_t  link_id, bool comp_flag
) {
    /* ---------- 1. Ensure arcs_ and link_aux_ capacity are synchronized ---------- */
    if (arcs_.size() == arcs_.capacity()) {               // hit capacity →  Expand 2×, starting at 16
        size_t old_cap = arcs_.capacity();
        size_t new_cap = old_cap ? old_cap << 1 : 16;
        arcs_.reserve(new_cap);
        link_aux_.reserve(new_cap);                       // Two vectors synchronized reserve
    }

    /* ---------- 2. First insert an empty link_aux slot to make "subscript = link_id" always true ---------- */
    link_aux_.emplace_back();

    /* The index of the current arc (also the default link_id) */
    size_t idx = arcs_.size();

    /* ---------- 3. Construct arc ---------- */
    GfaArc arc{};
    arc.v_lv = static_cast<uint64_t>(v) << 32;            // low 32bit (lv) -> 0
    arc.w    = w;
    arc.ov   = ov;
    arc.ow   = ow;
    arc.rank = -1;
    arc.set_comp(comp_flag);
    arc.set_del(false);
    arc.set_strong(false);

    uint64_t real_link_id = (link_id >= 0) ? static_cast<uint64_t>(link_id) : static_cast<uint64_t>(idx);
    arc.set_link_id(real_link_id);

    /* If this is a complementary arc (link_id >= 0), inherit the rank of the main arc */
    if (link_id >= 0 && link_id < static_cast<int64_t>(arcs_.size()))
        arc.rank = arcs_[real_link_id].rank;

    /* ---------- 4. push to arcs_ ---------- */
    arcs_.push_back(arc);

    return &arcs_.back();
}

/*------------------------------------------------------------*/
/*                      S-line parser                         */
/*------------------------------------------------------------*/
bool GfaGraph::parseSLine(std::stringstream& ss, const std::string& sample_name) {
    // ------------------------------------------------ S  Name  seq ------------------------------------------------
    char c; std::string name, seq;
    ss >> c >> name >> seq;
    uint32_t len = 0; bool ln_tag_found = false;
    std::string rest; std::getline(ss, rest);

    // ------------------------------------------------ AUX ------------------------------------------------
    auto aux_vec = GfaAuxParser::parse_aux_string(rest);

    // ------------------------------------------------ SN tag ------------------------------------------------
    std::vector<uint32_t> sample_ids;
    const uint8_t* pSN = GfaAuxParser::get_aux_data(aux_vec, "SN");
    if (pSN && *pSN == 'Z') {
        sample_ids = parse_sample_csv_(reinterpret_cast<const char*>(pSN + 1));
        has_node_source_tag_ = true;
        GfaAuxParser::delete_aux_tag(aux_vec, pSN);
    }

    uint32_t file_sample_id = intern_sample_(sample_name);
    if (file_sample_id != UINT32_MAX) {
        sample_ids.push_back(file_sample_id);
        has_node_source_tag_ = true;
    }

    std::sort(sample_ids.begin(), sample_ids.end());
    sample_ids.erase(std::unique(sample_ids.begin(), sample_ids.end()), sample_ids.end());

    // ------------------------------------------------ LN tag ------------------------------------------------
    const uint8_t* pLN = GfaAuxParser::get_aux_data(aux_vec, "LN");
    if (pLN && *pLN == 'i') {
        int32_t tmp = 0;
        std::memcpy(&tmp, pLN + 1, sizeof(tmp));
        if (tmp < 0) tmp = 0;
        len = static_cast<uint32_t>(tmp);
        ln_tag_found = true;
        GfaAuxParser::delete_aux_tag(aux_vec, pLN);
    }
    if (!ln_tag_found) {
        len = (seq == "*") ? 0u : static_cast<uint32_t>(seq.length());
    }

    // ------------------------------------------------ SR tag ------------------------------------------------
    const uint8_t* pSR = GfaAuxParser::get_aux_data(aux_vec, "SR");
    if (pSR && *pSR == 'i') {
        int32_t rank = 0;
        std::memcpy(&rank, pSR + 1, sizeof(rank));
        if (rank > 0 && static_cast<uint32_t>(rank) > max_rank_) {
            max_rank_ = static_cast<uint32_t>(rank);
        }
        GfaAuxParser::delete_aux_tag(aux_vec, pSR);
    }

    // ------------------------------------------------ CX tag (complex) ------------------------------------------------
    bool is_complex = false;
    const uint8_t* pCX = GfaAuxParser::get_aux_data(aux_vec, "CX");
    if (pCX && *pCX == 'i') {
        int32_t x = 0;
        std::memcpy(&x, pCX + 1, sizeof(x));
        is_complex = (x != 0);
        GfaAuxParser::delete_aux_tag(aux_vec, pCX);
    }

    // ------------------------------------------------ Node index ------------------------------------------------
    uint32_t id = get_or_add_segment(name);
    GfaNode& n = nodes_[id];
    n.name = name;
    std::transform(seq.begin(), seq.end(), seq.begin(), [](unsigned char c) { return std::toupper(c); });  // 2026-04-16
    n.sequence = seq;
    n.length = len;
    merge_sample_ids_into(n.sample_ids, sample_ids);
    n.is_complex = n.is_complex || is_complex;
    GfaAux incoming_aux;
    incoming_aux.aux_data = std::move(aux_vec);
    merge_variant_tags_into(incoming_aux, n.aux);
    n.aux = std::move(incoming_aux);
    total_segment_length_ += len;

    return true;
}

/*------------------------------------------------------------*/
/*                      L-line parser                         */
/*------------------------------------------------------------*/
bool GfaGraph::parseLLine(std::stringstream& ss) {
    char c;
    std::string v_name, w_name;
    char v_ori, w_ori;
    std::string ovlap_field;
    ss >> c >> v_name >> v_ori >> w_name >> w_ori >> ovlap_field;
    std::string rest; std::getline(ss, rest);

    uint32_t v_seg = get_or_add_segment(v_name);
    bool v_rev = (v_ori == '-');
    uint32_t v_vertex = (v_seg << 1) | (v_rev ? 1 : 0);

    uint32_t w_seg = get_or_add_segment(w_name);
    bool w_rev = (w_ori == '-');
    uint32_t w_vertex = (w_seg << 1) | (w_rev ? 1 : 0);

    /* parse overlap like original C code */
    int32_t ov = INT32_MAX, ow = INT32_MAX;
    if (ovlap_field == "*") {
        ov = ow = 0;
    } else if (!ovlap_field.empty() && ovlap_field[0] == ':') {
        /* ":OW" */
        try { ow = std::stoi(ovlap_field.substr(1)); } 
        catch (...) { ow = INT32_MAX; }
    } else if (!ovlap_field.empty() && std::isdigit(static_cast<unsigned char>(ovlap_field[0]))) {
        size_t pos = 0;
        while (pos < ovlap_field.size() && std::isdigit(static_cast<unsigned char>(ovlap_field[pos]))) ++pos;
        if (pos == ovlap_field.size()) {
            /* pure number */
            ov = ow = static_cast<int32_t>(std::stoi(ovlap_field));
        } else if (ovlap_field[pos] == ':') {
            /* "OV:OW" */
            try {
                ov = std::stoi(ovlap_field.substr(0, pos));
                ow = std::stoi(ovlap_field.substr(pos + 1));
            } catch (...) {
                ov = ow = 0;
            }
        } else {
            /* assumed CIGAR */
            ov = ow = 0;
            size_t i = 0;
            while (i < ovlap_field.size()) {
                size_t j = i;
                while (j < ovlap_field.size() && std::isdigit(static_cast<unsigned char>(ovlap_field[j]))) ++j;
                if (j == ovlap_field.size()) break;
                long val = std::strtol(ovlap_field.c_str() + i, nullptr, 10);
                char op = ovlap_field[j];
                if (op == 'M' || op == 'D' || op == 'N') ov += val;
                if (op == 'M' || op == 'I' || op == 'S') ow += val;
                i = j + 1;
            }
        }
    } else {
        ov = ow = 0;
    }

    if ((ov != INT32_MAX && ov > 0) || (ow != INT32_MAX && ow > 0)) {
        has_overlap_ = true;
    }

    /* add arc pair */
    GfaArc* prim = add_arc(v_vertex, w_vertex, ov, ow, -1, false);

    /* parse aux for link */
    auto aux_vec = GfaAuxParser::parse_aux_string(rest);
    GfaAux& laux = link_aux_[prim->get_link_id()];
    laux.aux_data = std::move(aux_vec);

    /* SR */
    const uint8_t* pSR = GfaAuxParser::get_aux_data(laux.aux_data, "SR");
    if (pSR && *pSR == 'i') {
        prim->rank = *reinterpret_cast<const int32_t*>(pSR + 1);
        GfaAuxParser::delete_aux_tag(laux.aux_data, pSR);
    }
    /* L1 L2 adjust lengths */
    const uint8_t* pL1 = GfaAuxParser::get_aux_data(laux.aux_data, "L1");
    if (pL1 && *pL1 == 'i') {
        int32_t l1 = *reinterpret_cast<const int32_t*>(pL1 + 1);
        if (ov != INT32_MAX) {
            uint32_t new_len = static_cast<uint32_t>(std::max<int64_t>(0, static_cast<int64_t>(ov) + l1));

            if (!nodes_[v_seg].sequence.empty() && nodes_[v_seg].sequence != "*") {
                uint32_t seq_len = static_cast<uint32_t>(nodes_[v_seg].sequence.size());
                if (new_len > seq_len) new_len = seq_len;
            }

            nodes_[v_seg].length = std::max<uint32_t>(nodes_[v_seg].length, new_len);

            uint32_t ov_u = (ov < 0 ? 0u : static_cast<uint32_t>(ov));
            if (ov_u > nodes_[v_seg].length) ov_u = nodes_[v_seg].length;
            prim->set_source_segment_len(nodes_[v_seg].length - ov_u);
        }
        GfaAuxParser::delete_aux_tag(laux.aux_data, pL1);
    }
    const uint8_t* pL2 = GfaAuxParser::get_aux_data(laux.aux_data, "L2");
    if (pL2 && *pL2 == 'i') {
        int32_t l2 = *reinterpret_cast<const int32_t*>(pL2 + 1);
        if (ow != INT32_MAX) {
            uint32_t new_len = static_cast<uint32_t>(std::max<int64_t>(0, static_cast<int64_t>(ow) + l2));

            if (!nodes_[w_seg].sequence.empty() && nodes_[w_seg].sequence != "*") {
                uint32_t seq_len = static_cast<uint32_t>(nodes_[w_seg].sequence.size());
                if (new_len > seq_len) new_len = seq_len;
            }

            nodes_[w_seg].length = std::max<uint32_t>(nodes_[w_seg].length, new_len);
        }
        GfaAuxParser::delete_aux_tag(laux.aux_data, pL2);
    }

    return true;
}

/*------------------------------------------------------------*/
/*                      P-line parser                         */
/*------------------------------------------------------------*/
bool GfaGraph::parsePLine(std::stringstream& ss) {
    char c; std::string pname, segs, cigar;
    ss >> c >> pname >> segs >> cigar;
    if (pname.empty() || segs.empty()) return false;

    GfaPath p; p.name = pname;
    std::stringstream segs_ss(segs);
    std::string tok;
    while (std::getline(segs_ss, tok, ',')) {
        if (tok.size() < 2) continue;
        std::string seg_name = tok.substr(0, tok.size() - 1);
        char ori = tok.back();
        auto it = name_to_id_map_.find(seg_name);
        if (it == name_to_id_map_.end()) {
            warning_stream() << "  ! Path references undefined segment '" << seg_name << "'\n";
            continue;
        }
        p.segments.push_back({it->second, ori=='-'});
    }
    paths_.push_back(std::move(p));
    return true;
}

/*------------------------------------------------------------*/
/*                      A-line parser                         */
/*------------------------------------------------------------*/
bool GfaGraph::parseALine(std::stringstream& ss) {
    char c; std::string unitig; uint32_t pos; char strand; std::string rname;
    uint32_t rstart, rend;
    ss >> c >> unitig >> pos >> strand >> rname >> rstart >> rend;
    std::string rest; std::getline(ss, rest);

    auto it = name_to_id_map_.find(unitig);
    if (it == name_to_id_map_.end()) {
        warning_stream() << "  ! A-line refers undefined unitig '" << unitig << "'\n";
        return false;
    }
    GfaAlignment a;
    a.unitig_node_id = it->second;
    a.position_on_unitig = pos;
    a.read_is_rev = (strand == '-');
    a.read_name = rname;
    a.read_start_pos = rstart;
    a.read_end_pos = rend;
    a.optional_tags = rest.empty()? "" : rest.substr(1); /* drop leading space */

    alignments_.push_back(std::move(a));
    ++total_alignments_;
    return true;
}

/*============================================================*/
/*           ------------  finalize_  ------------            */
/*============================================================*/
void GfaGraph::finalize_() {
    fixNoSeg();
    arcSort();
    arcIndex();
    fixSemiArc();
    fixSymmAdd();
    fixArcLen();
    cleanup();
    countLinks();
    calculateDeg();
}

/*------------------  step1: mark empty seg ------------------*/
void GfaGraph::fixNoSeg() {
    bool found = false;
    for (auto& n : nodes_) {
        if (!n.deleted && n.length == 0) {
            n.deleted = true;
            found = true;
            warning_stream() << "  ! Delete zero-length segment '" << n.name << "'\n";
        }
    }

    if (!found) return;
    for (GfaArc& a : arcs_) {
        const uint32_t v = a.get_source_segment_id();
        const uint32_t w = a.get_target_segment_id();
        if ((v < nodes_.size() && nodes_[v].deleted) || (w < nodes_.size() && nodes_[w].deleted)) {
            a.set_del(true);
        }
    }
}

/*------------------  step2: sort arcs   ---------------------*/
void GfaGraph::arcSort() {
    std::sort(arcs_.begin(), arcs_.end());
}

/*------------------  step3: build index -------------------- */
void GfaGraph::arcIndex() {
    uint64_t max_v = total_segments_ ? (total_segments_ * 2 - 1) : 0;
    arc_indexs_.assign(max_v + 1, 0);

    uint64_t cur_v = ~0ULL, start = 0; uint32_t cnt = 0;
    for (uint64_t i = 0; i < arcs_.size(); ++i) {
        uint32_t v = arcs_[i].get_source_vertex_id();
        if (v != cur_v) {
            if (cur_v != ~0ULL)
                arc_indexs_[cur_v] = (start << 32) | cnt;
            cur_v = v; start = i; cnt = 1;
        } else ++cnt;
    }
    if (cur_v != ~0ULL) arc_indexs_[cur_v] = (start << 32) | cnt;
}

/*------------------  step4: fill semi arcs ------------------*/
void GfaGraph::fixSemiArc() {
    std::unordered_map<uint32_t,std::vector<size_t>> outmap;
    for (size_t i = 0; i < arcs_.size(); ++i)
        outmap[arcs_[i].get_source_vertex_id()].push_back(i);

    for (size_t idx = 0; idx < arcs_.size(); ++idx) {
        GfaArc& a = arcs_[idx];
        if (a.get_del()) continue;
        bool needOv = (a.ov == INT32_MAX);
        bool needOw = (a.ow == INT32_MAX);
        if (!needOv && !needOw) continue;

        uint32_t wcomp = a.get_target_vertex_id() ^ 1;
        size_t found = SIZE_MAX;
        int multi = 0;

        auto it = outmap.find(wcomp);
        if (it != outmap.end()) {
            for (size_t j : it->second) {
                const GfaArc& b = arcs_[j];
                if (b.get_del()) continue;
                if (b.get_target_vertex_id() == (a.get_source_vertex_id() ^ 1)) {
                    ++multi;
                    found = j;
                }
            }
        }
        if (multi == 1) {
            const GfaArc& b = arcs_[found];
            if (needOw) a.ow = b.ov;
            if (needOv) a.ov = b.ow;
        } else {
            a.set_del(true);
            warning_stream() << "  ! arc " << nodes_[NodeHandle::get_segment_id(a.get_source_vertex_id())].name
                << " -> " << nodes_[NodeHandle::get_segment_id(a.get_target_vertex_id())].name
                << " removed (missing overlap)\n";
        }
    }
}

/*------------------  step5: add complement ------------------*/
void GfaGraph::fixSymmAdd() {
    auto findArc = [&](uint32_t v, uint32_t w, int32_t ov, int32_t ow) -> size_t {
        if (v >= arc_indexs_.size()) return SIZE_MAX;
        uint64_t pk  = arc_indexs_[v];
        uint64_t st  = pk >> 32;
        uint32_t cnt = (uint32_t)pk;
        for (uint32_t i = 0; i < cnt; ++i) {
            const GfaArc& a = arcs_[st + i];
            // if (a.get_del() || a.get_comp()) continue;
            if (a.get_del()) continue;
            if (a.w == w && a.ov == ow && a.ow == ov)
                return st + i;
        }
        return SIZE_MAX;
    };

    std::vector<GfaArc> addlist;

    for (const GfaArc& a : arcs_) {
        if (a.get_del() || a.get_comp()) continue;
        uint32_t v=a.get_source_vertex_id(), w=a.get_target_vertex_id();
        uint32_t w2=w^1, v2=v^1;

        size_t idx = findArc(w2, v2, a.ov, a.ow);
        if (idx == SIZE_MAX) {
            GfaArc b = a;
            b.v_lv = (static_cast<uint64_t>(w2)<<32) | nodes_[w2>>1].length;
            b.w = v2;
            std::swap(b.ov,b.ow);
            b.set_comp(true);
            addlist.push_back(b);
        } else {
            GfaArc& b = arcs_[idx];
            b.set_comp(true);
            b.set_link_id(a.get_link_id());
        }
    }
    if (!addlist.empty()) {
        arcs_.insert(arcs_.end(), addlist.begin(), addlist.end());
        arcSort();
        arcIndex();
    }
}

/*------------------  step6: adjust lengths ------------------*/
void GfaGraph::fixArcLen() {
    for (GfaArc& a : arcs_) {
        if (a.get_del()) continue;

        uint32_t sseg = NodeHandle::get_segment_id(a.get_source_vertex_id());
        uint32_t tseg = NodeHandle::get_segment_id(a.get_target_vertex_id());

        if (nodes_[sseg].deleted || nodes_[tseg].deleted) {
            a.set_del(true);
            continue;
        }

        uint32_t slen = (!nodes_[sseg].sequence.empty() && nodes_[sseg].sequence != "*")
                        ? static_cast<uint32_t>(nodes_[sseg].sequence.size())
                        : nodes_[sseg].length;

        uint32_t tlen = (!nodes_[tseg].sequence.empty() && nodes_[tseg].sequence != "*")
                        ? static_cast<uint32_t>(nodes_[tseg].sequence.size())
                        : nodes_[tseg].length;

        if (a.ov == INT32_MAX) a.ov = 0;
        if (a.ow == INT32_MAX) a.ow = 0;

        if (a.ov < 0) a.ov = 0;
        if (a.ow < 0) a.ow = 0;

        if (a.ov > static_cast<int32_t>(slen)) {
            warning_stream() << "  ! ov longer than source segment '" << nodes_[sseg].name << "': " << a.ov << " > " << slen << ", clamp to " << slen << '\n';
            a.ov = static_cast<int32_t>(slen);
            a.ow = a.ov;
        }

        if (a.ow > static_cast<int32_t>(tlen)) {
            warning_stream() << "  ! ow longer than target segment '" << nodes_[tseg].name << "': " << a.ow << " > " << tlen << ", clamp to " << tlen << '\n';
            a.ow = static_cast<int32_t>(tlen);
            a.ov = a.ow;
        }

        a.set_source_segment_len(slen - static_cast<uint32_t>(a.ov));
    }
}

/*------------------  step7: purge deleted -------------------*/
void GfaGraph::cleanup() {
    arcs_.erase(std::remove_if(arcs_.begin(), arcs_.end(), [](const GfaArc& a){ return a.get_del(); }), arcs_.end());
    arcSort();
    arcIndex();
}

/*------------------  step8: Count links (unique L-lines): one per forward / reverse pair -------------------*/
void GfaGraph::countLinks() {
    total_uniq_links_ = 0;
    for (const auto& arc : arcs_) {
        if (!arc.get_comp() && !arc.get_del())
            ++total_uniq_links_;
    }
}

/*------------------  step9: calculate degree -----------------*/
void GfaGraph::calculateDeg() {
    total_deg_ = 0;
    max_deg_ = 0;
    for (size_t i = 0; i < arc_indexs_.size(); ++i) {
        uint32_t count = static_cast<uint32_t>(arc_indexs_[i]);  // low 32 bits are the count
        total_deg_ += count;
        if (count > max_deg_) max_deg_ = count;
    }
    return;
}

/*============================================================*/
/*                  getters & summary                         */
/*============================================================*/
uint64_t GfaGraph::getNodeInternalId(const std::string& name) const {
    auto it = name_to_id_map_.find(name);
    if (it == name_to_id_map_.end())
        throw std::runtime_error("segment '" + name + "' not found");
    return it->second;
}
std::vector<size_t> GfaGraph::getArcsIdxFromVertex(uint32_t v, bool skip_self) const {
    std::vector<size_t> outidx;
    if (v >= arc_indexs_.size()) return outidx;

    uint64_t pk = arc_indexs_[v];
    uint64_t st = pk >> 32;
    uint32_t cnt = (uint32_t)pk;

    outidx.reserve(cnt);
    for (uint32_t i = 0; i < cnt; ++i) {
        size_t idx = (size_t)(st + i);
        const auto& a = arcs_[idx];
        if (a.get_del()) continue;
        if (skip_self && a.get_target_vertex_id() == v) continue;  // skip self-loop
        outidx.push_back(idx);
    }
    return outidx;
}
std::vector<const GfaArc*> GfaGraph::getArcsFromVertex(uint32_t v, bool skip_self) const {
    std::vector<const GfaArc*> out;
    if (v >= arc_indexs_.size()) return out;

    uint64_t pk = arc_indexs_[v];
    uint64_t st = pk >> 32;
    uint32_t cnt = (uint32_t)pk;

    out.reserve(cnt);
    for (uint32_t i = 0; i < cnt; ++i) {
        const auto& a = arcs_[st + i];
        if (a.get_del()) continue;
        if (skip_self && a.get_target_vertex_id() == v) continue;  // skip self-loop
        out.push_back(&a);
    }

    return out;
}
std::vector<size_t> GfaGraph::getArcsIdxToVertex(uint32_t v, bool skip_self) const {
    std::vector<size_t> v_inidx;

    uint32_t vr = v ^ 1;
    std::vector<const GfaArc*> vr_out = getArcsFromVertex(vr, skip_self);
    v_inidx.reserve(vr_out.size());

    for (const GfaArc* erev : vr_out) {
        if (!erev) continue;
        uint32_t erev_w = erev->get_target_vertex_id();
        uint32_t erev_w_rev = erev_w ^ 1;
        for (size_t idx : getArcsIdxFromVertex(erev_w_rev, skip_self)) {
            const GfaArc& e = arcs_[idx];
            if (e.get_del() || e.get_target_vertex_id() != v) continue;
            v_inidx.push_back(idx);
        }
    }
    return v_inidx;
}
std::vector<const GfaArc*> GfaGraph::getArcsToVertex(uint32_t v, bool skip_self) const {
    std::vector<const GfaArc*> in;

    uint32_t vr = v ^ 1;
    std::vector<const GfaArc*> vr_out = getArcsFromVertex(vr, skip_self);
    in.reserve(vr_out.size());

    for (const GfaArc* erev : vr_out) {
        if (!erev) continue;
        uint32_t erev_w = erev->get_target_vertex_id();
        uint32_t erev_w_rev = erev_w ^ 1;
        for (const GfaArc* a : getArcsFromVertex(erev_w_rev, skip_self)) {
            if (a->get_del() || a->get_target_vertex_id() != v) continue;
            in.push_back(a);
        }
    }

    return in;
}

std::vector<std::vector<uint32_t>> GfaGraph::enumerate_paths_greedy_DFS(
    const uint32_t src, const uint32_t sink, const std::unordered_set<uint32_t>& region_set,
    const uint32_t max_depth, const uint32_t max_paths,
    const bool skip_comp, bool& hit_limits, const uint64_t DFS_guard, const uint32_t stall_round_limit
) const {
    GfaWalker walker(*this);
    return walker.enumerate_paths_greedy_DFS(src, sink, region_set, max_depth, max_paths, skip_comp, hit_limits, DFS_guard, stall_round_limit);
};

std::vector<std::vector<uint32_t>> GfaGraph::open_walk(
    const uint32_t src,
    const std::unordered_set<uint32_t>& region_set,
    const uint32_t max_depth,
    const bool skip_comp,
    const uint64_t DFS_guard,
    const uint64_t walk_bp,
    const std::unordered_set<uint32_t>* blocked_seg
) const {
    GfaWalker walker(*this);
    return walker.open_walk(src, region_set, max_depth, skip_comp, DFS_guard, walk_bp, blocked_seg);
}


bool GfaGraph::build_nodes_connectivity_index(bool skip_comp)
{
    log_stream() << "Building nodes connectivity index ...\n";

    std::vector<uint32_t>().swap(connectivity_index_);  // clear and free memory
    std::vector<uint32_t>().swap(connectivity_sizes_);

    const uint32_t nseg = static_cast<uint32_t>(getNumNodes());
    connectivity_index_.assign(nseg, UINT32_MAX);
    if (nseg == 0) return false;

    std::vector<uint32_t> parent(nseg, UINT32_MAX);
    std::vector<uint8_t> rank(nseg, 0);
    for (uint32_t s = 0; s < nseg; ++s) {
        if (!getNodeDeleted(s)) parent[s] = s;
    }

    auto find_root = [&](uint32_t x) -> uint32_t {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](uint32_t a, uint32_t b) {
        uint32_t ra = find_root(a);
        uint32_t rb = find_root(b);
        if (ra == rb) return;
        if (rank[ra] < rank[rb]) std::swap(ra, rb);
        parent[rb] = ra;
        if (rank[ra] == rank[rb]) ++rank[ra];
    };

    for (const GfaArc& a : getAllArcs()) {
        if (a.get_del() || (skip_comp && a.get_comp())) continue;
        const uint32_t s = a.get_source_segment_id();
        const uint32_t t = a.get_target_segment_id();
        if (s >= nseg || t >= nseg) continue;
        if (parent[s] == UINT32_MAX || parent[t] == UINT32_MAX) continue;
        unite(s, t);
    }

    std::vector<uint32_t> root_to_comp(nseg, UINT32_MAX);
    uint32_t ncomp = 0;
    for (uint32_t s = 0; s < nseg; ++s) {
        if (parent[s] == UINT32_MAX) continue;
        const uint32_t r = find_root(s);
        if (root_to_comp[r] == UINT32_MAX) {
            root_to_comp[r] = ncomp++;
        }
        connectivity_index_[s] = root_to_comp[r];
    }

    connectivity_sizes_.assign(ncomp, 0);
    for (uint32_t cid : connectivity_index_) {
        if (cid != UINT32_MAX) ++connectivity_sizes_[cid];
    }

    log_stream() << "  - Nodes connectivity components indexed: " << ncomp << "\n" << "\n";

    return true;
}

uint32_t GfaGraph::node_connectivity_id(Vertex v) const
{
    if (connectivity_index_.empty()) {
        warning_stream() << "  ! Connectivity index not built yet\n";
        return UINT32_MAX;
    }
    const uint32_t sid = v.segment_id();
    if (sid >= connectivity_index_.size()) return UINT32_MAX;
    return connectivity_index_[sid];
}

uint32_t GfaGraph::node_connectivity_size(Vertex v) const
{
    const uint32_t cid = node_connectivity_id(v);
    return cid < connectivity_sizes_.size() ? connectivity_sizes_[cid] : 0;
}

bool GfaGraph::nodes_connected(Vertex a, Vertex b) const
{
    const uint32_t ca = node_connectivity_id(a);
    const uint32_t cb = node_connectivity_id(b);
    return (ca != UINT32_MAX && cb != UINT32_MAX && ca == cb);
}

bool GfaGraph::build_vertex_topological_index(bool skip_comp)
{
    log_stream() << "Building nodes topological index ...\n";

    const uint32_t nseg = static_cast<uint32_t>(getNumNodes());
    const uint32_t V = nseg * 2;

    topo_index_.clear();
    topo_index_.resize(V);

    auto valid_vertex = [&](uint32_t v) -> bool {
        const uint32_t sid = NodeHandle::get_segment_id(v);
        return sid < nseg && !getNodeDeleted(sid);
    };

    std::vector<uint8_t> self_loop(V, 0);
    std::vector<uint32_t> disc(V, UINT32_MAX);
    std::vector<uint32_t> low(V, UINT32_MAX);
    std::vector<uint8_t> on_stack(V, 0);
    std::vector<uint32_t> stack;
    stack.reserve(V);

    std::vector<uint32_t> comp_id(V, UINT32_MAX);
    std::vector<uint32_t> comp_size;
    std::vector<uint32_t> comp_min_v;

    uint32_t timer = 0;

    std::function<void(uint32_t)> dfs = [&](uint32_t u) {
        disc[u] = low[u] = timer++;
        stack.push_back(u);
        on_stack[u] = 1;

        for (const GfaArc* a : getArcsFromVertex(u)) {
            if (!a || a->get_del() || (skip_comp && a->get_comp())) continue;

            const uint32_t w = a->get_target_vertex_id();
            if (!valid_vertex(w)) continue;

            if (w == u) {
                self_loop[u] = 1;
                continue;
            }

            if (disc[w] == UINT32_MAX) {
                dfs(w);
                low[u] = std::min(low[u], low[w]);
            } else if (on_stack[w]) {
                low[u] = std::min(low[u], disc[w]);
            }
        }

        if (low[u] != disc[u]) return;

        const uint32_t cid = static_cast<uint32_t>(comp_size.size());
        uint32_t sz = 0;
        uint32_t min_v = UINT32_MAX;

        while (true) {
            uint32_t x = stack.back();
            stack.pop_back();
            on_stack[x] = 0;

            comp_id[x] = cid;
            ++sz;
            min_v = std::min(min_v, x);

            if (x == u) break;
        }

        comp_size.push_back(sz);
        comp_min_v.push_back(min_v);
    };

    for (uint32_t v = 0; v < V; ++v) {
        if (!valid_vertex(v)) continue;
        if (disc[v] == UINT32_MAX) dfs(v);
    }

    const uint32_t C = static_cast<uint32_t>(comp_size.size());

    std::vector<std::vector<uint32_t>> dag(C);
    std::vector<uint32_t> indeg(C, 0);

    for (uint32_t u = 0; u < V; ++u) {
        if (!valid_vertex(u)) continue;

        const uint32_t cu = comp_id[u];
        if (cu == UINT32_MAX) continue;

        for (const GfaArc* a : getArcsFromVertex(u)) {
            if (!a || a->get_del() || (skip_comp && a->get_comp())) continue;

            const uint32_t w = a->get_target_vertex_id();
            if (!valid_vertex(w)) continue;

            const uint32_t cw = comp_id[w];
            if (cw == UINT32_MAX || cu == cw) continue;

            dag[cu].push_back(cw);
        }
    }

    for (uint32_t c = 0; c < C; ++c) {
        auto& xs = dag[c];
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());

        for (uint32_t y : xs) {
            ++indeg[y];
        }
    }

    auto cmp = [&](uint32_t a, uint32_t b) {
        return comp_min_v[a] > comp_min_v[b];
    };

    std::priority_queue<uint32_t, std::vector<uint32_t>, decltype(cmp)> q(cmp);

    for (uint32_t c = 0; c < C; ++c) {
        if (indeg[c] == 0) q.push(c);
    }

    std::vector<uint32_t> comp_rank(C, UINT32_MAX);
    uint32_t rank = 0;

    while (!q.empty()) {
        const uint32_t c = q.top();
        q.pop();

        comp_rank[c] = rank++;

        for (uint32_t y : dag[c]) {
            if (--indeg[y] == 0) q.push(y);
        }
    }

    for (uint32_t v = 0; v < V; ++v) {
        if (!valid_vertex(v)) continue;

        const uint32_t c = comp_id[v];
        topo_index_[v].scc_id = c;
        topo_index_[v].topo_rank = comp_rank[c];
        topo_index_[v].scc_size = comp_size[c];
        topo_index_[v].in_cycle = (comp_size[c] > 1) || self_loop[v];
    }

    log_stream() << "  - SCCs: " << C << "\n";
    log_stream() << "  - Node topological ranks indexed\n\n";

    if (DEBUG_ENABLED) {
        std::vector<uint32_t> vertices;
        vertices.reserve(V);

        for (uint32_t v = 0; v < V; ++v) {
            if (valid_vertex(v)) vertices.push_back(v);
        }

        std::sort(vertices.begin(), vertices.end(), [&](uint32_t a, uint32_t b) {
            if (topo_index_[a].topo_rank != topo_index_[b].topo_rank)
                return topo_index_[a].topo_rank < topo_index_[b].topo_rank;
            return a < b;
        });

        for (uint32_t v : vertices) {
            const uint32_t sid = Vertex::get_segment_id(v);
            const std::string seg_name = getNodeName(sid);
            const char strand = Vertex::get_is_reverse(v) ? '-' : '+';

            debug_stream()
                << "  - rank=" << topo_index_[v].topo_rank
                << " name=" << seg_name << strand
                << " scc_id=" << topo_index_[v].scc_id
                << " scc_size=" << topo_index_[v].scc_size
                << " in_cycle=" << static_cast<int>(topo_index_[v].in_cycle)
                << "\n";
        }

        debug_stream() << "\n";
    }

    return topo_index_.size() == V;
}

uint32_t GfaGraph::vertex_topo_rank(Vertex v) const
{
    if (topo_index_.empty()) {
        warning_stream() << "  ! Topological index not built yet\n";
        return UINT32_MAX;
    }
    const uint32_t vid = v.vertex_id();
    if (vid >= topo_index_.size()) return UINT32_MAX;
    return topo_index_[vid].topo_rank;
}

 bool GfaGraph::vertex_in_cycle(Vertex v) const
{
    if (topo_index_.empty()) {
        warning_stream() << "  ! Topological index not built yet\n";
        return false;
    }
    const uint32_t vid = v.vertex_id();
    if (vid >= topo_index_.size()) return false;
    return topo_index_[vid].in_cycle;
}

const GfaTopoIndex& GfaGraph::vertex_topo_info(Vertex v) const
{
    static const GfaTopoIndex empty;
    if (topo_index_.empty()) {
        warning_stream() << "  ! Topological index not built yet\n";
        return empty;
    }
    const uint32_t vid = v.vertex_id();
    if (vid >= topo_index_.size()) return empty;
    return topo_index_[vid];
}

void GfaGraph::detect_complex_regions(
    uint16_t branch_degree,
    uint16_t hub_degree,
    uint32_t min_nodes,
    uint32_t min_branches,
    bool skip_comp
) {
    log_stream() << "Detecting complex graph regions ...\n";

    if (nodes_.empty()) return;

    const uint32_t segment_count = static_cast<uint32_t>(nodes_.size());
    std::vector<std::vector<uint32_t>> adjacent(segment_count);
    std::vector<uint8_t> self_loop(segment_count, 0);

    for (const GfaArc& arc : arcs_) {
        if (arc.get_del() || (skip_comp && arc.get_comp())) continue;

        const uint32_t v = arc.get_source_vertex_id();
        const uint32_t w = arc.get_target_vertex_id();
        const uint32_t from = Vertex::get_segment_id(v);
        const uint32_t to = Vertex::get_segment_id(w);
        if (from >= segment_count || to >= segment_count) continue;
        if (nodes_[from].deleted || nodes_[to].deleted) continue;

        if (from == to) {
            self_loop[from] = 1;
        } else {
            adjacent[from].push_back(to);
            adjacent[to].push_back(from);
        }
    }

    for (std::vector<uint32_t>& xs : adjacent) {
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    }

    auto mark = [&](uint32_t sid, uint64_t& count) {
        if (sid >= segment_count || nodes_[sid].deleted || nodes_[sid].is_complex) return;
        nodes_[sid].is_complex = true;
        ++count;
    };

    uint64_t cycle_nodes = 0;
    uint64_t tangle_blocks = 0;
    uint64_t tangle_nodes = 0;

    for (uint32_t sid = 0; sid < segment_count; ++sid) {
        if (self_loop[sid]) mark(sid, cycle_nodes);
    }

    {
        struct DfsFrame { uint32_t node; size_t next; };

        std::vector<uint32_t> discover(segment_count, UINT32_MAX);
        std::vector<uint32_t> low(segment_count, UINT32_MAX);
        std::vector<uint32_t> parent(segment_count, UINT32_MAX);
        std::vector<uint64_t> edge_stack;
        std::vector<DfsFrame> dfs;

        std::vector<uint32_t> block_degree(segment_count, 0);
        std::vector<uint32_t> block_stamp(segment_count, 0);
        std::vector<uint32_t> block_index(segment_count, 0);
        std::vector<uint32_t> block_nodes;
        uint32_t time = 0;
        uint32_t stamp = 0;

        auto emit_block = [&](uint32_t stop_from, uint32_t stop_to) {
            if (++stamp == 0) {
                std::fill(block_stamp.begin(), block_stamp.end(), 0);
                stamp = 1;
            }

            block_nodes.clear();
            uint32_t block_edges = 0;
            while (!edge_stack.empty()) {
                const uint64_t packed = edge_stack.back();
                edge_stack.pop_back();
                const uint32_t a = static_cast<uint32_t>(packed >> 32);
                const uint32_t b = static_cast<uint32_t>(packed);

                if (block_stamp[a] != stamp) {
                    block_stamp[a] = stamp;
                    block_degree[a] = 0;
                    block_nodes.push_back(a);
                }
                if (block_stamp[b] != stamp) {
                    block_stamp[b] = stamp;
                    block_degree[b] = 0;
                    block_nodes.push_back(b);
                }
                ++block_degree[a];
                ++block_degree[b];
                ++block_edges;
                if (a == stop_from && b == stop_to) break;
            }

            for (uint32_t i = 0; i < block_nodes.size(); ++i) {
                const uint32_t v = block_nodes[i];
                block_index[v] = i;
            }

            const uint32_t n = static_cast<uint32_t>(block_nodes.size());
            std::vector<uint8_t> visited(n, 0);
            std::vector<uint8_t> local_seen(n, 0);
            std::vector<uint32_t> work;
            std::vector<uint32_t> group;
            std::vector<uint32_t> local_nodes;
            std::vector<uint32_t> complex_nodes;
            std::vector<uint8_t> selected(n, 0);
            uint32_t complex_groups = 0;

            for (uint32_t seed = 0; seed < n; ++seed) {
                if (visited[seed] || block_degree[block_nodes[seed]] < branch_degree) continue;

                visited[seed] = 1;
                work.assign(1, seed);
                group.clear();
                bool group_has_hub = false;

                while (!work.empty()) {
                    const uint32_t i = work.back();
                    work.pop_back();
                    group.push_back(i);
                    group_has_hub = group_has_hub || block_degree[block_nodes[i]] >= hub_degree;

                    for (uint32_t w : adjacent[block_nodes[i]]) {
                        if (block_stamp[w] != stamp) continue;
                        const uint32_t j = block_index[w];
                        if (!visited[j] && block_degree[w] >= branch_degree) {
                            visited[j] = 1;
                            work.push_back(j);
                        }
                    }
                }

                local_nodes.clear();
                for (uint32_t i : group) {
                    if (!local_seen[i]) {
                        local_seen[i] = 1;
                        local_nodes.push_back(i);
                    }
                    for (uint32_t w : adjacent[block_nodes[i]]) {
                        if (block_stamp[w] != stamp) continue;
                        const uint32_t j = block_index[w];
                        if (!local_seen[j]) {
                            local_seen[j] = 1;
                            local_nodes.push_back(j);
                        }
                    }
                }

                const bool complex_group = group.size() >= min_branches &&
                    (group_has_hub || local_nodes.size() >= min_nodes);
                if (complex_group) {
                    ++complex_groups;
                    for (uint32_t i : group) {
                        selected[i] = 1;
                        complex_nodes.push_back(block_nodes[i]);
                    }
                }
                for (uint32_t i : local_nodes) local_seen[i] = 0;
            }

            // Recover simple paths internal to the tangle, but do not extend along arms
            // that touch the complex core at only one end.
            std::vector<uint8_t> simple_seen(n, 0);
            std::vector<uint8_t> touch_seen(n, 0);
            std::vector<uint32_t> touches;
            for (uint32_t seed = 0; seed < n; ++seed) {
                if (simple_seen[seed] || block_degree[block_nodes[seed]] >= branch_degree) continue;

                simple_seen[seed] = 1;
                work.assign(1, seed);
                group.clear();
                touches.clear();

                while (!work.empty()) {
                    const uint32_t i = work.back();
                    work.pop_back();
                    group.push_back(i);

                    for (uint32_t w : adjacent[block_nodes[i]]) {
                        if (block_stamp[w] != stamp) continue;
                        const uint32_t j = block_index[w];
                        if (block_degree[w] < branch_degree) {
                            if (!simple_seen[j]) {
                                simple_seen[j] = 1;
                                work.push_back(j);
                            }
                        } else if (selected[j] && !touch_seen[j]) {
                            touch_seen[j] = 1;
                            touches.push_back(j);
                        }
                    }
                }

                if (touches.size() >= 2) {
                    for (uint32_t i : group) {
                        selected[i] = 1;
                        complex_nodes.push_back(block_nodes[i]);
                    }
                }
                for (uint32_t i : touches) touch_seen[i] = 0;
            }

            const bool complex = !complex_nodes.empty();

            if (DEBUG_ENABLED) {
                std::vector<uint32_t> boundaries;
                for (uint32_t v : block_nodes) {
                    for (uint32_t w : adjacent[v]) {
                        if (block_stamp[w] != stamp) {
                            boundaries.push_back(v);
                            break;
                        }
                    }
                }
                std::sort(boundaries.begin(), boundaries.end());

                std::ostringstream line;
                line << "  - Block boundary=";
                if (boundaries.empty()) {
                    line << '-';
                } else {
                    for (size_t i = 0; i < boundaries.size(); ++i) {
                        if (i) line << ',';
                        line << nodes_[boundaries[i]].name;
                    }
                }
                const double average_degree = block_nodes.empty() ? 0.0 : 2.0 * static_cast<double>(block_edges) / block_nodes.size();
                line << " complex=" << (complex ? "yes" : "no")
                     << " complex_nodes=" << complex_nodes.size()
                     << " avg_degree=" << std::fixed << std::setprecision(2) << average_degree << '\n';
                debug_stream() << line.str();
            }

            if (complex) {
                tangle_blocks += complex_groups;
                for (uint32_t sid : complex_nodes) mark(sid, tangle_nodes);
            }
        };

        for (uint32_t root = 0; root < segment_count; ++root) {
            if (nodes_[root].deleted || discover[root] != UINT32_MAX) continue;
            discover[root] = low[root] = time++;
            dfs.push_back({root, 0});

            while (!dfs.empty()) {
                DfsFrame& frame = dfs.back();
                const uint32_t u = frame.node;
                if (frame.next < adjacent[u].size()) {
                    const uint32_t v = adjacent[u][frame.next++];
                    if (discover[v] == UINT32_MAX) {
                        parent[v] = u;
                        edge_stack.push_back((static_cast<uint64_t>(u) << 32) | v);
                        discover[v] = low[v] = time++;
                        dfs.push_back({v, 0});
                    } else if (v != parent[u] && discover[v] < discover[u]) {
                        low[u] = std::min(low[u], discover[v]);
                        edge_stack.push_back((static_cast<uint64_t>(u) << 32) | v);
                    }
                    continue;
                }

                dfs.pop_back();
                const uint32_t p = parent[u];
                if (p == UINT32_MAX) continue;
                low[p] = std::min(low[p], low[u]);
                if (low[u] >= discover[p]) emit_block(p, u);
            }
        }
    }

    uint64_t total = 0;
    for (const GfaNode& n : nodes_) {
        if (!n.deleted && n.is_complex) {
            ++total;
            log_stream() << "  - Complex node: " << n.name << " (length=" << n.length << ")\n";
        }
    }

    log_stream() << "  - Self-loop nodes: " << cycle_nodes << "\n";
    log_stream() << "  - Tangle blocks: " << tangle_blocks << "\n";
    log_stream() << "  - Tangle nodes: " << tangle_nodes << "\n";
    log_stream() << "  - Total complex nodes: " << total << "\n\n";
}

void GfaGraph::mark_segments_complex(const std::vector<uint32_t>& seg_ids)
{
    log_stream() << "Marking segments as complex ...\n";

    size_t marked = 0;

    for (uint32_t seg_id : seg_ids) {
        if (seg_id >= nodes_.size()) continue;
        if (nodes_[seg_id].deleted) continue;
        if (nodes_[seg_id].is_complex) continue;

        nodes_[seg_id].is_complex = true;
        ++marked;
    }

    log_stream() << "  - Segments marked as complex: " << marked << "\n\n";

    return;
}

namespace { inline uint64_t clamp0_i64(int64_t x){ return x>0 ? (uint64_t)x : 0ULL; } }
bool GfaGraph::shortest_distance_between_offsets(
    uint32_t v_from, uint32_t off_from,
    uint32_t v_to,   uint32_t off_to,
    uint64_t& out_dist,
    uint32_t step_cap
) const {
    const uint32_t V = (uint32_t)arc_indexs_.size();
    out_dist = std::numeric_limits<uint64_t>::max()/4;

    if (v_from >= V || v_to >= V) return false;
    const uint32_t sf = NodeHandle::get_segment_id(v_from);
    const uint32_t st = NodeHandle::get_segment_id(v_to);
    if (sf >= nodes_.size() || st >= nodes_.size()) return false;
    if (nodes_[sf].deleted || nodes_[st].deleted)   return false;

    // same segment
    if (v_from == v_to) {
        uint32_t len = nodes_[sf].length;
        uint32_t lo  = std::min(off_from, len);
        uint32_t hi  = std::min(off_to,   len);
        out_dist = (hi >= lo) ? (hi - lo) : 0ULL;
        return true;
    }

    const auto get_len = [&](uint32_t v)->uint32_t{
        uint32_t s = NodeHandle::get_segment_id(v);
        if (s >= nodes_.size() || nodes_[s].deleted) return 0u;
        return nodes_[s].length;
    };

    const uint64_t INF = std::numeric_limits<uint64_t>::max()/4;
    std::vector<uint64_t> dist(V, INF);

    // Record predecessors (only for reconstruction; shortest distance only)
    struct Prev{ uint32_t p=UINT32_MAX; int32_t ow=0; bool has=false; };
    std::vector<Prev> prev(V);

    // Initial cost: consume source tail first
    uint64_t init = clamp0_i64((int64_t)get_len(v_from) - (int64_t)off_from);
    dist[v_from] = init;

    using QN = std::pair<uint64_t,uint32_t>;
    std::priority_queue<QN, std::vector<QN>, std::greater<QN>> pq;
    pq.emplace(init, v_from);

    uint64_t best = INF;
    uint32_t pops = 0;

    while (!pq.empty()) {
        auto [d,u] = pq.top(); pq.pop();
        if (d != dist[u]) continue;
        if (++pops > step_cap) break;
        if (best != INF && d >= best) break; // No better path possible

        auto outs = getArcsFromVertex(u);
        for (const GfaArc* a : outs) {
            if (!a || a->get_del()) continue;
            uint32_t w = a->get_target_vertex_id();
            uint32_t wl = get_len(w);
            // If this edge reaches target, compute terminal cost
            if (w == v_to) {
                uint64_t add = clamp0_i64((int64_t)off_to - (int64_t)a->ow);
                uint64_t cand = d + add;
                if (cand < best) best = cand;
            }
            // Continue: append w's [ow, len)
            uint64_t addw = clamp0_i64((int64_t)wl - (int64_t)a->ow);
            uint64_t nd   = d + addw;
            if (nd < dist[w]) {
                dist[w] = nd;
                prev[w] = {u, a->ow, true};
                pq.emplace(nd, w);
            }
        }
    }
    if (best == INF) return false;
    out_dist = best;
    return true;
}


uint32_t GfaGraph::add_segment(
    const std::string& name,
    const std::string& sequence,
    bool name_check,
    const std::vector<uint32_t>& sample_ids
) {
    auto it = name_to_id_map_.find(name);
    if (it != name_to_id_map_.end()) {
        uint32_t existing_id = it->second;
        GfaNode& node = nodes_[existing_id];
        merge_sample_ids_into(node.sample_ids, sample_ids);

        const std::string& existing_seq = node.sequence;
        if (name_check) {
            if (existing_seq == sequence) return existing_id;

            if (node.deleted) {
                node.deleted = false;
                node.sequence = sequence;
                node.length = (sequence == "*") ? 0u : static_cast<uint32_t>(sequence.size());
                total_segment_length_ += node.length;
                return existing_id;
            }

            error_stream() << "Duplicated segment name with different sequence: " << name << "\n";
            std::exit(1);
        }

        return existing_id;
    }

    uint32_t id = static_cast<uint32_t>(total_segments_);
    name_to_id_map_[name] = id;
    nodes_.resize(id + 1);

    nodes_[id].name = name;
    nodes_[id].sequence = sequence;
    nodes_[id].length = (sequence == "*") ? 0u : static_cast<uint32_t>(sequence.size());
    nodes_[id].deleted = false;
    nodes_[id].sample_ids = sample_ids;

    total_segment_length_ += nodes_[id].length;
    ++total_segments_;

    arc_indexs_.resize(total_segments_ * 2, 0);
    return id;
}

bool GfaGraph::delete_segment(uint32_t seg_id) {
    if (seg_id >= nodes_.size()) return false;
    GfaNode &n = nodes_[seg_id];
    if (n.deleted) return true;

    n.deleted = true;
    total_segment_length_ -= n.length;
    n.length = 0;
    n.sequence.clear();

    uint32_t v_plus  = (seg_id << 1) | 0;
    uint32_t v_minus = (seg_id << 1) | 1;

    for (size_t ei : getArcsIdxFromVertex(v_plus))
        arcs_[ei].set_del(true);
    for (size_t ei : getArcsIdxToVertex(v_plus))
        arcs_[ei].set_del(true);
    for (size_t ei : getArcsIdxFromVertex(v_minus))
        arcs_[ei].set_del(true);
    for (size_t ei : getArcsIdxToVertex(v_minus))
        arcs_[ei].set_del(true);

    return true;
}

void GfaGraph::rebuild_after_edits() {
    fixArcLen();   // Refresh lv = slen-ov per arc and clamp ov
    cleanup();     // Remove deleted arcs
    arcSort();     // Sort by v
    arcIndex();    // Rebuild index
    fixSymmAdd();  // Add missing complements and align link_id
}

PathSequence GfaGraph::extend_left_from(
    uint32_t start_vertex,
    uint32_t min_bases,
    uint32_t max_steps
) const {
    PathSequence out;
    if (start_vertex >= arc_indexs_.size()) return out;

    static constexpr int kHardStepCap = 10000;

    // Avoid head inserts: collect chunks then reverse-append (latest on right)
    std::vector<std::string> chunks;
    chunks.reserve(8);

    // Append suffix of predecessor non-overlap region (near current boundary)
    auto append_suffix_of_prefix = [&](uint32_t v, uint32_t upto, uint32_t need) -> uint32_t {
        const uint32_t seg = NodeHandle::get_segment_id(v);
        if (seg >= nodes_.size()) return 0u;
        const auto& nd = nodes_[seg];
        if (nd.deleted || nd.length == 0) return 0u;

        const std::string oriented = get_oriented_sequence(Vertex(v));
        if (oriented.empty() || oriented == "*") return 0u;

        const uint32_t lim  = std::min<uint32_t>(upto, nd.length);
        if (lim == 0) return 0u;

        const uint32_t take = (need == 0) ? lim : std::min<uint32_t>(need, lim);
        const uint32_t beg  = lim - take; // Take suffix nearest boundary
        if (take) chunks.emplace_back(oriented.substr(beg, take));
        return take;
    };

    auto done = [&](uint32_t gathered_bases, uint32_t steps) -> bool {
        const bool base_ok = (min_bases == 0) ? true : (gathered_bases >= min_bases);
        const bool step_ok = (max_steps == 0) ? true : (steps          >= max_steps);
        return base_ok && step_ok;
    };

    // Instant appendable length: lv = |prev| - ov; unknown ov => full length
    auto instant_lv = [&](const GfaArc* a) -> uint32_t {
        if (!a) return 0u;
        const uint32_t prev_v   = a->get_source_vertex_id();
        const uint32_t prev_seg = NodeHandle::get_segment_id(prev_v);
        if (prev_seg >= nodes_.size()) return 0u;
        const auto& prev = nodes_[prev_seg];
        if (prev.deleted || prev.length == 0) return 0u;
        if (a->ov == INT32_MAX) return prev.length;
        const uint32_t ov_u = (a->ov < 0) ? 0u : static_cast<uint32_t>(a->ov);
        return (prev.length > ov_u) ? (prev.length - ov_u) : 0u;
    };

    uint32_t cur_v = start_vertex;
    uint32_t gathered = 0, steps = 0;

    for (int guard = 0; guard < kHardStepCap; ++guard) {
        if (done(gathered, steps)) break;

        auto in_arcs = getArcsToVertex(cur_v);
        if (in_arcs.empty()) break;

        // Roots of current tip (for intersection priority)
        std::vector<std::string> cur_roots;
        {
            const uint32_t cur_seg = NodeHandle::get_segment_id(cur_v);
            if (cur_seg < nodes_.size() && !nodes_[cur_seg].deleted) {
                gfaName::collect_roots_from_name(nodes_[cur_seg].name, cur_roots);
            }
        }

        // Prefer root intersection; otherwise max lv
        int chosen = -1; uint32_t best_inter = 0;
        int chosen_all = -1; uint32_t best_all = 0;

        for (int i = 0; i < (int)in_arcs.size(); ++i) {
            const GfaArc* a = in_arcs[i];
            if (!a || a->get_del()) continue;

            const uint32_t prev_v   = a->get_source_vertex_id();
            const uint32_t prev_seg = NodeHandle::get_segment_id(prev_v);
            if (prev_seg >= nodes_.size()) continue;
            const auto& prev = nodes_[prev_seg];
            if (prev.deleted || prev.length == 0) continue;

            const uint32_t lv = instant_lv(a);

            bool inter = false;
            if (!cur_roots.empty()) {
                std::vector<std::string> prev_roots; prev_roots.reserve(4);
                gfaName::collect_roots_from_name(prev.name, prev_roots);
                inter = gfaName::roots_intersect(cur_roots, prev_roots);
            }
            if (inter) {
                if (lv > best_inter) { best_inter = lv; chosen = i; }
            } else {
                if (lv > best_all)   { best_all   = lv; chosen_all = i; }
            }
        }
        if (chosen < 0) chosen = chosen_all;
        if (chosen < 0) break;

        const GfaArc* a = in_arcs[chosen];
        const uint32_t prev_v   = a->get_source_vertex_id();
        const uint32_t prev_seg = NodeHandle::get_segment_id(prev_v);
        if (prev_seg >= nodes_.size() || nodes_[prev_seg].deleted) break;

        const uint32_t lv = instant_lv(a);
        const uint32_t need_now = (min_bases == 0) ? 0u : (min_bases - std::min(min_bases, gathered));
        gathered += append_suffix_of_prefix(prev_v, lv, need_now);

        // Record path (exclude start; begin at first neighbor)
        out.vertexs.push_back(prev_v);

        cur_v = prev_v;
        ++steps;
    }

    // Reverse-append: latest on right
    for (auto it = chunks.rbegin(); it != chunks.rend(); ++it)
        out.sequence.append(*it);

    return out;
}

PathSequence GfaGraph::extend_right_from(
    uint32_t start_vertex,
    uint32_t min_bases,
    uint32_t max_steps
) const {
    PathSequence out;
    if (start_vertex >= arc_indexs_.size()) return out;

    static constexpr int kHardStepCap = 10000;

    // Append tail of target vertex (skip current edge ow)
    auto append_tail_from = [&](uint32_t v, int32_t ow, uint32_t need) -> uint32_t {
        const uint32_t seg = NodeHandle::get_segment_id(v);
        if (seg >= nodes_.size()) return 0u;
        const auto& nd = nodes_[seg];
        if (nd.deleted || nd.length == 0) return 0u;

        const std::string oriented = get_oriented_sequence(Vertex(v));
        if (oriented.empty() || oriented == "*") return 0u;

        uint32_t ow_u = (ow == INT32_MAX) ? 0u : (ow < 0 ? 0u : static_cast<uint32_t>(ow));
        if (ow_u > nd.length) ow_u = nd.length;

        const uint32_t avail = (nd.length > ow_u) ? (nd.length - ow_u) : 0u;
        const uint32_t take  = (need == 0) ? avail : std::min<uint32_t>(need, avail);
        if (take) out.sequence.append(oriented, ow_u, take);
        return take;
    };

    // Instant appendable tail: tail = |next| - ow; unknown ow => full length
    auto instant_tail = [&](const GfaArc* a) -> uint32_t {
        if (!a) return 0u;
        const uint32_t w   = a->get_target_vertex_id();
        const uint32_t seg = NodeHandle::get_segment_id(w);
        if (seg >= nodes_.size()) return 0u;
        const auto& nd = nodes_[seg];
        if (nd.deleted || nd.length == 0) return 0u;
        if (a->ow == INT32_MAX) return nd.length;
        const uint32_t ow_u = (a->ow < 0) ? 0u : static_cast<uint32_t>(a->ow);
        return (nd.length > ow_u) ? (nd.length - ow_u) : 0u;
    };

    auto done = [&](uint32_t gathered_bases, uint32_t steps) -> bool {
        const bool base_ok = (min_bases == 0) ? true : (gathered_bases >= min_bases);
        const bool step_ok = (max_steps == 0) ? true : (steps          >= max_steps);
        return base_ok && step_ok;
    };

    uint32_t cur_v = start_vertex;
    uint32_t gathered = 0, steps = 0;

    for (int guard = 0; guard < kHardStepCap; ++guard) {
        if (done(gathered, steps)) break;

        auto outs = getArcsFromVertex(cur_v);
        if (outs.empty()) break;

        // Roots of current tip (for intersection priority)
        std::vector<std::string> cur_roots;
        {
            const uint32_t seg = NodeHandle::get_segment_id(cur_v);
            if (seg < nodes_.size() && !nodes_[seg].deleted)
                gfaName::collect_roots_from_name(nodes_[seg].name, cur_roots);
        }

        // Prefer root intersection; otherwise max tail
        int chosen = -1; uint32_t best_inter = 0;
        int chosen_all = -1; uint32_t best_all = 0;

        for (int i = 0; i < (int)outs.size(); ++i) {
            const GfaArc* a = outs[i];
            if (!a || a->get_del()) continue;

            const uint32_t w_seg = NodeHandle::get_segment_id(a->get_target_vertex_id());
            if (w_seg >= nodes_.size()) continue;
            const auto& wnode = nodes_[w_seg];
            if (wnode.deleted || wnode.length == 0) continue;

            const uint32_t tail = instant_tail(a);

            bool inter = false;
            if (!cur_roots.empty()) {
                std::vector<std::string> nb; nb.reserve(4);
                gfaName::collect_roots_from_name(wnode.name, nb);
                inter = gfaName::roots_intersect(cur_roots, nb);
            }
            if (inter) {
                if (tail > best_inter) { best_inter = tail; chosen = i; }
            } else {
                if (tail > best_all)   { best_all   = tail; chosen_all = i; }
            }
        }
        if (chosen < 0) chosen = chosen_all;
        if (chosen < 0) break;

        const GfaArc* a = outs[chosen];
        const uint32_t w = a->get_target_vertex_id();

        const uint32_t need_now = (min_bases == 0) ? 0u : (min_bases - std::min(min_bases, gathered));
        gathered += append_tail_from(w, a->ow, need_now);

        // Record path (exclude start; begin at first neighbor)
        out.vertexs.push_back(w);

        cur_v = w;
        ++steps;
    }

    return out;
}

ExpandedSeqs GfaGraph::getSeqVec(uint32_t k) const {
    const uint32_t need0 = (k > 0 ? k - 1 : 0);
    const uint64_t n     = nodes_.size();

    ExpandedSeqs outs;
    outs.names.resize(n, "");
    outs.seqs.resize(n, "");
    outs.right_seqs.resize(n);

    // Trivial path: no extension
    if (need0 == 0 || n == 0) {
        for (uint32_t seg = 0; seg < n; ++seg) {
            outs.names[seg] = nodes_[seg].name;
            const auto& node = nodes_[seg];
            if (!node.deleted && node.sequence != "*") {
                const size_t real_len = (node.length > 0 && node.length <= node.sequence.size())
                                        ? size_t(node.length) : node.sequence.size();
                outs.seqs[seg] = std::string_view(node.sequence.data(), real_len);
            } else {
                outs.seqs[seg] = std::string_view();
            }
            outs.right_seqs[seg].clear();
        }
        return outs;
    }

    for (uint32_t seg = 0; seg < n; ++seg) {
        const auto& node = nodes_[seg];

        // name + own sequence (zero-copy view)
        outs.names[seg] = node.name;
        if (!node.deleted && node.sequence != "*") {
            const size_t real_len = (node.length > 0 && node.length <= node.sequence.size())
                                    ? size_t(node.length) : node.sequence.size();
            outs.seqs[seg] = std::string_view(node.sequence.data(), real_len);
        } else {
            outs.seqs[seg] = std::string_view();
        }

        if (node.deleted || node.length == 0) {
            outs.right_seqs[seg].clear();
            continue;
        }

        const uint32_t v = (seg << 1) | 0; // forward vertex

        // Build path and assembled sequence for (k-1) bp
        PathSequence p = extend_right_from(v, need0, 0);
        outs.right_seqs[seg].push_back(std::move(p.sequence)); // keep legacy output in sync
    }

    return outs;
}

std::string GfaGraph::format_overlap_field(int32_t ov, int32_t ow) const {
    if (ov == INT32_MAX || ow == INT32_MAX) return "*";
    if (ov == ow) return std::to_string(ov) + "M";
    return std::to_string(ov) + ":" + std::to_string(ow);
}

void GfaGraph::save_to_disk(
    const std::string& filename,
    bool write_paths,
    bool write_align,
    bool write_seq,
    const std::string& command_line
) const {
    SAVE saver(filename);

    std::ostringstream oss;

    // 0. Header
    oss << "H\tVN:Z:1.0\n";
    oss << "H\tTS:Z:liftasm\n";
    if (!command_line.empty()) {
        std::string cl = command_line;
        for (char& c : cl) {
            if (c == '\t' || c == '\n' || c == '\r') c = ' ';
        }
        oss << "H\tCL:Z:" << cl << "\n";
    }
    oss << "H\tSP:Z:" << (has_node_source_tag_ ? "present" : "absent") << "\n";
    oss << "H\tCO:Z:  LN:i   segment length\n";
    oss << "H\tCO:Z:  CX:i   complex node flag\n";
    oss << "H\tCO:Z:  SN:Z   node source/sample names\n";
    oss << "H\tCO:Z:  VT:Z   VCF-added variant category (--tag value)\n";
    oss << "H\tCO:Z:  VI:Z   VCF variant ID\n";
    saver.save(oss.str());
    oss.str("");
    oss.clear();

    // 1. S-lines
    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        const auto& n = nodes_[sid];
        if (n.deleted) continue;
        oss << "S\t" << n.name << '\t';
        if (n.sequence.empty() || n.sequence == "*" || !write_seq) oss << '*';
        else oss << n.sequence;
        oss << "\tLN:i:" << n.length;
        const std::string sn = sample_ids_to_csv_(n.sample_ids);
        if (!sn.empty()) oss << "\tSN:Z:" << sn;
        if (n.is_complex) oss << "\tCX:i:1";
        oss << GfaAuxParser::write_aux_as_text(n.aux.aux_data);
        oss << '\n';
        saver.save(oss.str()); oss.str(""); oss.clear();
    }

    // 2. L-lines
    for (const auto& a : arcs_) {
        if (a.get_del() || a.get_comp()) continue;

        uint32_t v = a.get_source_vertex_id();
        uint32_t w = a.get_target_vertex_id();
        uint32_t vseg = NodeHandle::get_segment_id(v);
        uint32_t wseg = NodeHandle::get_segment_id(w);
        bool vrev = NodeHandle::get_is_reverse(v);
        bool wrev = NodeHandle::get_is_reverse(w);

        const std::string& vname = nodes_[vseg].name;
        const std::string& wname = nodes_[wseg].name;

        oss << "L\t" << vname << '\t' << (vrev?'-':'+') << '\t'
            << wname << '\t' << (wrev?'-':'+') << '\t'
            << format_overlap_field(a.ov, a.ow);

        // aux
        uint64_t lid = a.get_link_id();
        if (lid < link_aux_.size()) {
            oss << GfaAuxParser::write_aux_as_text(link_aux_[lid].aux_data);
        }
        oss << '\n';
    }
    saver.save(oss.str()); oss.str(""); oss.clear();

    // 3. P-lines
    if (write_paths) {
        for (const auto& p : paths_) {
            oss << "P\t" << p.name << '\t';
            // segs: "name(+/-),name(+/-),..."
            for (size_t i = 0; i < p.segments.size(); ++i) {
                const auto& s = p.segments[i];
                const std::string& nm = nodes_[s.node_id].name;
                if (i) oss << ',';
                oss << nm << (s.is_reverse?'-':'+');
            }
            oss << "\t*\n";
            saver.save(oss.str()); oss.str(""); oss.clear();
        }
    }

    // 4. A-lines
    if (write_align) {
        for (const auto& a : alignments_) {
            oss << "A\t"
                << nodes_[a.unitig_node_id].name << '\t'
                << a.position_on_unitig << '\t'
                << (a.read_is_rev?'-':'+') << '\t'
                << a.read_name << '\t'
                << a.read_start_pos << '\t'
                << a.read_end_pos;
            if (!a.optional_tags.empty()) oss << '\t' << a.optional_tags;
            oss << '\n';
            saver.save(oss.str()); oss.str(""); oss.clear();
        }
    }

    return;
}

std::vector<std::string> GfaGraph::getAllSegmentNames() const {
    std::vector<std::string> names;
    names.reserve(nodes_.size());
    for (const auto& n : nodes_) names.push_back(n.name);
    return names;
}

// print summary statistics
void GfaGraph::print_graph_stats() const {
    const int label_width = 25, value_width = 12;
    log_stream() << std::left << std::setw(label_width) << "GFA Graph stats:" << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Segments (S):" << std::right << std::setw(value_width) << getNumNodes() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Links (unique):" << std::right << std::setw(value_width) << getNumUniqLinks() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Directed arcs:" << std::right << std::setw(value_width) << arcs_.size() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Max rank:" << std::right << std::setw(value_width) << getMaxRank() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Paths (P):" << std::right << std::setw(value_width) << getNumPaths() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Alignments (A):" << std::right << std::setw(value_width) << getNumAlignments() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Total segment length:" << std::right << std::setw(value_width) << total_segment_length_ << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Average segment len:" << std::right << std::setw(value_width)
                 << std::fixed << std::setprecision(3) << getAveSegmentsLength() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Total degree:" << std::right << std::setw(value_width) << getTotalDeg() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Max degree:" << std::right << std::setw(value_width) << getMaxDeg() << '\n';
    log_stream() << "  - " << std::left << std::setw(label_width) << "Average degree:" << std::right << std::setw(value_width)
                 << std::fixed << std::setprecision(3) << getAveDeg() << '\n' << "\n";
}

// debug
void GfaGraph::printArcList() const {
    log_stream() << "=== Arc list ===\n";
    for (size_t i = 0; i < arcs_.size(); ++i) {
        // if (arcs_[i].get_del() || arcs_[i].get_comp()) continue;
        if (arcs_[i].get_del()) continue;
        const GfaArc& a = arcs_[i];
        uint32_t v_id = a.get_source_vertex_id();
        uint32_t w_id = a.get_target_vertex_id();
        uint32_t seg_v = NodeHandle::get_segment_id(v_id);
        uint32_t seg_w = NodeHandle::get_segment_id(w_id);
        bool rev_v = NodeHandle::get_is_reverse(v_id);
        bool rev_w = NodeHandle::get_is_reverse(w_id);
        const std::string v_str = getNodeName(seg_v) + (rev_v ? "-" : "+");
        const std::string w_str = getNodeName(seg_w) + (rev_w ? "-" : "+");

        // Print arc with names and metadata
        log_stream() << "arc[" << i << "]: "
                     << v_str << "(v=" << v_id << ") -> " << w_str << "(w=" << w_id << ")"
                     << ", ov=" << a.ov << ", ow=" << a.ow
                     << ", link_id=" << a.get_link_id()
                     << ", rank=" << a.rank
                     << ", del=" << a.get_del()
                     << ", comp=" << a.get_comp()
                     << ", strong=" << a.get_strong()
                     << ", lv=" << a.get_source_segment_len()
                     << '\n';
    }
}
void GfaGraph::printArcIndex() const {
    log_stream() << "=== Index array ===\n";
    for (size_t i = 0; i < arc_indexs_.size(); ++i) {
        uint64_t val = arc_indexs_[i];
        uint32_t count = static_cast<uint32_t>(val);
        uint64_t start = val >> 32;
        log_stream() << "idx[" << i << "]: start=" << start << ", count=" << count << '\n';
    }
}


// Used to identify samples from S-line tags and -n/--name, and to manage sample ID mappings. (v0.1.3-r9 2026-06-03)
uint32_t GfaGraph::intern_sample_(const std::string& name) {
    if (name.empty()) return UINT32_MAX;

    auto it = sample_name_to_id_.find(name);
    if (it != sample_name_to_id_.end()) return it->second;

    uint32_t id = static_cast<uint32_t>(sample_id_to_name_.size());
    sample_id_to_name_.push_back(name);
    sample_name_to_id_.emplace(name, id);
    return id;
}

std::vector<uint32_t> GfaGraph::parse_sample_csv_(const std::string& csv) {
    std::vector<uint32_t> ids;
    std::stringstream ss(csv);
    std::string x;

    while (std::getline(ss, x, ',')) {
        if (x.empty()) continue;
        uint32_t id = intern_sample_(x);
        if (id != UINT32_MAX) ids.push_back(id);
    }

    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

std::string GfaGraph::sample_ids_to_csv_(const std::vector<uint32_t>& ids) const {
    std::string out;
    for (uint32_t id : ids) {
        if (id >= sample_id_to_name_.size()) continue;
        if (!out.empty()) out.push_back(',');
        out += sample_id_to_name_[id];
    }
    return out;
}

void GfaGraph::merge_sample_ids_into(std::vector<uint32_t>& dst, const std::vector<uint32_t>& src) const {
    if (src.empty()) return;
    if (dst.empty()) {
        dst = src;
        return;
    }

    std::vector<uint32_t> merged;
    merged.reserve(dst.size() + src.size());
    std::set_union(dst.begin(), dst.end(), src.begin(), src.end(), std::back_inserter(merged));
    dst.swap(merged);
}

void GfaGraph::merge_variant_tags_into(GfaAux& dst, const GfaAux& src) const {
    auto merge_one = [&](const char tag[2]) {
        std::vector<std::string> values;
        auto collect = [&](const std::vector<uint8_t>& aux) {
            const uint8_t* ptr = GfaAuxParser::get_aux_data(aux, tag);
            if (!ptr || *ptr != 'Z') return;
            std::string csv(reinterpret_cast<const char*>(ptr + 1));
            size_t beg = 0;
            while (beg <= csv.size()) {
                const size_t end = csv.find(',', beg);
                std::string value = csv.substr(beg, end == std::string::npos ? std::string::npos : end - beg);
                if (!value.empty()) values.push_back(std::move(value));
                if (end == std::string::npos) break;
                beg = end + 1;
            }
        };

        collect(dst.aux_data);
        collect(src.aux_data);
        if (values.empty()) return;
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());

        const uint8_t* old = GfaAuxParser::get_aux_data(dst.aux_data, tag);
        if (old) GfaAuxParser::delete_aux_tag(dst.aux_data, old);
        std::string csv = values.front();
        for (size_t i = 1; i < values.size(); ++i) csv += "," + values[i];
        GfaAuxParser::append_str_tag(dst.aux_data, tag, csv);
    };

    merge_one("VT");
    merge_one("VI");
}

bool GfaGraph::sample_ids_intersect(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) const {
    size_t i = 0;
    size_t j = 0;

    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) return true;
        if (a[i] < b[j]) ++i;
        else ++j;
    }

    return false;
}

bool GfaGraph::sample_ids_equal_or_empty(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) const {
    if (a.empty() || b.empty()) return true;
    return a == b;
}

bool GfaGraph::intersect_node_sample(uint32_t a_seg, uint32_t b_seg) const {
    if (a_seg >= nodes_.size() || b_seg >= nodes_.size()) return true;

    const auto& a = getNodeSampleIds(a_seg);
    const auto& b = getNodeSampleIds(b_seg);

    if (a.empty() || b.empty()) return true;
    return sample_ids_intersect(a, b);
}

bool GfaGraph::intersect_vertex_sample(uint32_t a_vtx, uint32_t b_vtx) const {
    return intersect_node_sample(
        NodeHandle::get_segment_id(a_vtx),
        NodeHandle::get_segment_id(b_vtx)
    );
}

bool GfaGraph::same_node_sample(uint32_t a_seg, uint32_t b_seg) const {
    if (a_seg >= nodes_.size() || b_seg >= nodes_.size()) return true;

    return sample_ids_equal_or_empty(
        getNodeSampleIds(a_seg),
        getNodeSampleIds(b_seg)
    );
}

bool GfaGraph::same_vertex_sample(uint32_t a_vtx, uint32_t b_vtx) const {
    return same_node_sample(
        NodeHandle::get_segment_id(a_vtx),
        NodeHandle::get_segment_id(b_vtx)
    );
}

const std::vector<std::string>& GfaGraph::getSampleNames() const { return sample_id_to_name_; }

const std::vector<uint32_t>& GfaGraph::getNodeSampleIds(uint32_t seg_id) const {
    static const std::vector<uint32_t> empty;
    return seg_id < nodes_.size() ? nodes_[seg_id].sample_ids : empty;
}
