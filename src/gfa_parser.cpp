#include "../include/gfa_parser.hpp"
#include "../include/gfa_walker.hpp"
#include "../include/gfa_name.hpp"

#include <iomanip>
#include <cstring>
#include <limits>
#include <iterator>
#include <stack>
#include <queue>
#include <unordered_set>
#include <algorithm>

using namespace wfa; // WFAlignerGapAffine

GfaGraph::GfaGraph()  = default;
GfaGraph::~GfaGraph() = default;

// Load GFA file into the graph structure
void GfaGraph::load_from_GFA(const std::string& filename) {
    log_stream() << "Loading genome graph from '" << filename << "' ...\n";

    /* pass-1: register S names */
    {
        GzChunkReader zr(filename);
        std::string line;
        while (zr.readLine(line)) {
            if (line.empty() || line[0] != 'S') continue;
            if (line.back() == '\n' || line.back() == '\r') line.pop_back();
            std::stringstream ss(line);
            char c; std::string name;
            ss >> c >> name;
            if (name.empty()) continue;
            if (name_to_id_map_.insert({name, total_segments_}).second)
                ++total_segments_;
        }
        zr.close();
    }
    nodes_.resize(total_segments_);

    /* pass-2: full parse */
    {
        GzChunkReader zr(filename);
        std::string line;
        while (zr.readLine(line)) {
            if (line.empty()) continue;
            if (line.back() == '\n' || line.back() == '\r') line.pop_back();
            std::stringstream ss(line);
            char type; ss >> type;
            ss.clear(); ss.str(line);
            switch (type) {
                case 'S': parseSLine(ss); break;
                case 'L': parseLLine(ss); break;
                case 'P': parsePLine(ss); break;
                case 'A': parseALine(ss); break;
                default:  break;
            }
        }
        zr.close();
    }

    if (forbid_overlap_ && has_overlap_) {
        error_stream() << "Overlap-style GFA is not supported. Please check the input file.\n";
        std::exit(1);
    }

    // Finalize the graph after loading
    finalize_();

    return;
}

uint32_t GfaGraph::get_or_add_segment(const std::string& name) {
    auto it = name_to_id_map_.find(name);
    if (it != name_to_id_map_.end()) return it->second;

    /* unseen before */
    uint32_t id = static_cast<uint32_t>(total_segments_);
    name_to_id_map_[name] = id;
    nodes_.resize(id + 1);
    nodes_[id].name = name;
    ++total_segments_;
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

std::string GfaGraph::get_oriented_sequence(uint32_t vertex) const {
    uint32_t seg_id = NodeHandle::get_segment_id(vertex);
    bool rev = NodeHandle::get_is_reverse(vertex);
    const std::string& seq = nodes_[seg_id].sequence;
    if (seq == "*") return std::string(); // unknown
    return rev ? seqUtils::revcomp(seq) : seq;
}

std::string GfaGraph::get_path_sequence(const std::vector<uint32_t>& path) const {
    if (path.empty()) return "*";

    size_t total_len = 0;
    for (uint32_t v : path) {
        uint32_t seg_id = NodeHandle::get_segment_id(v);
        total_len += getNodeLength(seg_id);
    }

    std::string out;
    out.reserve(total_len);

    for (uint32_t v : path) {
        std::string s = get_oriented_sequence(v);
        if (s.empty()) return "*";
        out += s;
    }
    return out;
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
    arcs_.push_back(arc);                                 // Now arcs_.size() is equal to link_aux_.size()

    return &arcs_.back();
}

/*------------------------------------------------------------*/
/*                      S-line parser                         */
/*------------------------------------------------------------*/
bool GfaGraph::parseSLine(std::stringstream& ss) {
    char c; std::string name, seq;
    ss >> c >> name >> seq;
    uint32_t len = 0; bool ln_tag_found = false;
    std::string rest; std::getline(ss, rest);

    auto aux_vec = GfaAuxParser::parse_aux_string(rest);

    // check for LN tag
    const uint8_t* pLN = GfaAuxParser::get_aux_data(aux_vec, "LN");
    if (pLN && *pLN == 'i') {
        // len = *reinterpret_cast<const int32_t*>(pLN + 1);
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

    // check for SR tag
    const uint8_t* pSR = GfaAuxParser::get_aux_data(aux_vec, "SR");
    if (pSR && *pSR == 'i') {
        int32_t rank = *reinterpret_cast<const int32_t*>(pSR + 1);
        if (static_cast<uint32_t>(rank) > max_rank_) {
            max_rank_ = static_cast<uint32_t>(rank);
        }
        GfaAuxParser::delete_aux_tag(aux_vec, pSR);
    }

    uint32_t id = get_or_add_segment(name);
    GfaNode& n = nodes_[id];
    n.name = name;
    n.sequence = seq;
    n.length = len;
    n.aux.aux_data = std::move(aux_vec);
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
            nodes_[v_seg].length = std::max<uint32_t>(nodes_[v_seg].length, static_cast<uint32_t>(ov + l1));
            prim->set_source_segment_len(nodes_[v_seg].length - ov);
        }
        GfaAuxParser::delete_aux_tag(laux.aux_data, pL1);
    }
    const uint8_t* pL2 = GfaAuxParser::get_aux_data(laux.aux_data, "L2");
    if (pL2 && *pL2 == 'i') {
        int32_t l2 = *reinterpret_cast<const int32_t*>(pL2 + 1);
        if (ow != INT32_MAX) {
            nodes_[w_seg].length = std::max<uint32_t>(nodes_[w_seg].length, static_cast<uint32_t>(ow + l2));
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
            warning_stream() << "Path references undefined segment '" << seg_name << "'\n";
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
        warning_stream() << "A-line refers undefined unitig '" << unitig << "'\n";
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
    for (auto& n : nodes_) {
        if (n.length == 0) {
            n.deleted = true;
            warning_stream() << "Delete zero-length segment '" << n.name << "'\n";
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
            warning_stream() << "[W] arc " << nodes_[NodeHandle::get_segment_id(a.get_source_vertex_id())].name
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
            a.set_del(true); continue;
        }
    
        uint32_t slen = nodes_[sseg].length;
        if (a.ov > static_cast<int32_t>(slen)) {
            warning_stream() << "ov longer than segment '" << nodes_[sseg].name << "': " << a.ov << " > " << slen << '\n';
            a.ov = slen;
        }
        a.set_source_segment_len(slen - a.ov);
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

// ---- GfaGraph::walk_dfs ----
std::vector<PathSequence>
GfaGraph::walk_dfs(
    uint32_t start_v,
    bool to_direction,
    uint32_t max_paths,
    uint32_t max_steps,
    uint32_t /*max_bases*/
) const {
    GfaWalker walker(*this);
    WalkDir dir = to_direction ? WalkDir::To : WalkDir::From;
    return walker.dfs(start_v, dir, max_paths, max_steps);
}

// ---- GfaGraph::walk_bfs ----
std::vector<PathSequence>
GfaGraph::walk_bfs(
    uint32_t start_v,
    bool to_direction,
    uint32_t max_paths,
    uint32_t max_steps,
    uint32_t /*max_bases*/
) const {
    GfaWalker walker(*this);
    WalkDir dir = to_direction ? WalkDir::To : WalkDir::From;
    return walker.bfs(start_v, dir, max_paths, max_steps);
}

 std::vector<std::vector<uint32_t>> GfaGraph::enumerate_paths_DFS(
    const uint32_t src, const uint32_t sink, const std::unordered_set<uint32_t>& region_set,
    const uint32_t max_depth, const uint32_t max_paths,
    const bool skip_comp, bool& hit_limits, const uint64_t DFS_guard
) const {
    GfaWalker walker(*this);
    return walker.enumerate_paths_DFS(src, sink, region_set, max_depth, max_paths, skip_comp, hit_limits, DFS_guard);
};


bool GfaGraph::is_connected_vertex(uint32_t v_from, uint32_t v_to, uint32_t step_cap) const {
    const uint32_t V = (uint32_t)arc_indexs_.size();
    if (v_from >= V || v_to >= V) return false;

    // Skip if source/target segments are deleted
    auto ok_seg = [&](uint32_t v)->bool{
        uint32_t s = NodeHandle::get_segment_id(v);
        return s < nodes_.size() && !nodes_[s].deleted;
    };
    if (!ok_seg(v_from) || !ok_seg(v_to)) return false;

    std::vector<uint8_t> vis(V, 0);
    std::queue<uint32_t> q;
    q.push(v_from); vis[v_from] = 1;

    uint32_t pops = 0;
    while (!q.empty()) {
        uint32_t u = q.front(); q.pop();
        if (u == v_to) return true;
        if (++pops > step_cap) break;

        auto outs = getArcsFromVertex(u);
        for (const GfaArc* a : outs) {
            if (!a || a->get_del()) continue;
            uint32_t w = a->get_target_vertex_id();
            if (!ok_seg(w)) continue;
            if (!vis[w]) { vis[w] = 1; q.push(w); }
        }
    }
    return false;
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


uint32_t GfaGraph::add_segment(const std::string& name, const std::string& sequence, bool name_check) {
    auto it = name_to_id_map_.find(name);
    if (it != name_to_id_map_.end()) {
        uint32_t existing_id = it->second;
        auto& node = nodes_[existing_id];
        const std::string& existing_seq = node.sequence;

        if (name_check) {
            if (existing_seq == sequence) {  // 2025-10-10: allow identical seq
                return existing_id;
            } else if (node.deleted) {  // 2025-10-10: previously deleted, now re-adding
                node.deleted = false;
                node.sequence = sequence;
                node.length = (sequence == "*") ? 0u : static_cast<uint32_t>(sequence.size());
                total_segment_length_ += node.length;
                return existing_id;
            } else {
                error_stream() << "Duplicated segment name with different sequence: " << name << "\n";
                error_stream() << "   - Existing: " << existing_seq << "\n";
                error_stream() << "   - New: " << sequence << "\n";
                std::exit(1);
            }
        } else {
            return existing_id;
        }
    }
    uint32_t id = static_cast<uint32_t>(total_segments_);
    name_to_id_map_[name] = id;
    nodes_.resize(id + 1);
    nodes_[id].name = name;
    nodes_[id].sequence = sequence;
    nodes_[id].length = (sequence=="*") ? 0u : (uint32_t)sequence.size();
    nodes_[id].deleted = false;
    total_segment_length_ += nodes_[id].length;
    ++total_segments_;

    // expand arc_indexs_ size: 2 vertexs per new segment
    arc_indexs_.resize(total_segments_*2, 0);
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

        const std::string oriented = get_oriented_sequence(v);
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

        const std::string oriented = get_oriented_sequence(v);
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
    bool write_paths, bool write_align, bool write_seq
) const {
    SAVE saver(filename);

    std::ostringstream oss;

    // 1. S-lines
    for (uint32_t sid = 0; sid < nodes_.size(); ++sid) {
        const auto& n = nodes_[sid];
        if (n.deleted) continue;
        oss << "S\t" << n.name << '\t';
        if (n.sequence.empty() || n.sequence == "*" || !write_seq) oss << '*';
        else oss << n.sequence;
        oss << "\tLN:i:" << n.length;
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
    log_stream() << "   - " << std::left << std::setw(label_width) << "Segments (S):" << std::right << std::setw(value_width) << getNumNodes() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Links (unique):" << std::right << std::setw(value_width) << getNumUniqLinks() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Directed arcs:" << std::right << std::setw(value_width) << arcs_.size() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Max rank:" << std::right << std::setw(value_width) << getMaxRank() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Paths (P):" << std::right << std::setw(value_width) << getNumPaths() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Alignments (A):" << std::right << std::setw(value_width) << getNumAlignments() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Total segment length:" << std::right << std::setw(value_width) << total_segment_length_ << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Average segment len:" << std::right << std::setw(value_width)
                 << std::fixed << std::setprecision(3) << getAveSegmentsLength() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Total degree:" << std::right << std::setw(value_width) << getTotalDeg() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Max degree:" << std::right << std::setw(value_width) << getMaxDeg() << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Average degree:" << std::right << std::setw(value_width)
                 << std::fixed << std::setprecision(3) << getAveDeg() << '\n';
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