#pragma once
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <cstring>

#include "kpaf.hpp"
#include "gfa_parser.hpp"

namespace mapconv {

struct MapRec {
    std::string unitig; uint32_t u_beg=0, u_end=0;
    std::string read;  uint32_t r_beg=0, r_end=0;
    char strand='+';
};

static inline bool maprec_less(const MapRec& a, const MapRec& b) {
    if (a.unitig != b.unitig) return a.unitig < b.unitig;
    if (a.u_beg  != b.u_beg ) return a.u_beg  < b.u_beg;
    if (a.u_end  != b.u_end ) return a.u_end  < b.u_end;
    if (a.read   != b.read  ) return a.read   < b.read;
    if (a.r_beg  != b.r_beg ) return a.r_beg  < b.r_beg;
    return a.r_end < b.r_end;
}

static inline bool has_suffix(const std::string& s, const char* suf) {
    const size_t n = std::strlen(suf);
    return s.size() >= n && s.compare(s.size()-n, n, suf) == 0;
}

static inline void write_one_map(std::ostream& out, const std::string& unitig, uint32_t ub, uint32_t ue, const std::string& read,   uint32_t rb, uint32_t re, char strand)
{
    out << unitig << ":" << ub << "-" << ue << "+\t" << read   << ":" << rb << "-" << re << strand << "\n";
}

// ---------- GFA2map ----------
static inline void gfa_to_map(const std::string& gfa_path, std::ostream& out) {
    log_stream() << "Converting GFA alignments to map format: " << gfa_path << "\n";

    GfaGraph g; g.load_from_GFA(gfa_path);
    const auto& aligns = g.getAlignments();

    std::vector<MapRec> rows;
    rows.reserve(aligns.size());

    for (const auto& a : aligns) {
        MapRec r;
        r.unitig = g.getNodeName(a.unitig_node_id);
        const uint32_t pos  = a.position_on_unitig;
        const uint32_t rlen = (a.read_end_pos >= a.read_start_pos) ? (a.read_end_pos - a.read_start_pos) : 0u;
        r.u_beg = pos; r.u_end = pos + rlen;
        r.read  = a.read_name;
        r.r_beg = a.read_start_pos; r.r_end = a.read_end_pos;
        r.strand = a.read_is_rev ? '-' : '+';
        rows.push_back(std::move(r));
    }

    std::sort(rows.begin(), rows.end(), maprec_less);
    for (const auto& r : rows)
        write_one_map(out, r.unitig, r.u_beg, r.u_end, r.read, r.r_beg, r.r_end, r.strand);
}

// ---------- PAF2map ----------
static inline void paf_emit_match_blocks(const paf::Record& p, std::ostream& out) {
    char cg_ty = 0; std::string_view cg_sv;
    const bool has_cg = paf::find_tag(p, "cg", cg_ty, cg_sv) && cg_ty=='Z' && !cg_sv.empty();

    if (!has_cg) {
        error_stream() << "Missing CG tag in PAF file. Please use minimap2 with the -c option.\n";
        std::exit(1);
    }

    uint32_t tpos= p.tstart, qofs=0, run_t0=0, run_q0=0, run_len=0;

    auto flush = [&]() {
        if (!run_len) return;
        const uint32_t ub = run_t0, ue = run_t0 + run_len;
        uint32_t rb, re;
        if (p.strand == '+') { rb = p.qstart + run_q0; re = rb + run_len; }
        else { re = p.qend - run_q0; rb = p.qend - (run_q0 + run_len); }
        write_one_map(out, p.tname, ub, ue, p.qname, rb, re, p.strand);
        run_len = 0;
    };
    auto start = [&]() { if (!run_len) { run_t0 = tpos; run_q0 = qofs; } };

    uint32_t num = 0;
    for (char c : cg_sv) {
        if (c >= '0' && c <= '9') { num = num*10u + uint32_t(c-'0'); continue; }
        const uint32_t len = num; num = 0;

        switch (c) {
            // case '=': case 'M': start(); run_len += len; tpos += len; qofs += len; break;
            // case 'X': flush(); tpos += len; qofs += len; break;
            case '=': case 'M': case 'X': start(); run_len += len; tpos += len; qofs += len; break;
            case 'I': case 'S': flush();  qofs += len; break;
            case 'D': case 'N': flush();  tpos += len; break;
            case 'H': case 'P': break;
            default:  flush(); break;
        }
    }
    flush();
}

static inline bool keep_paf_record(const paf::Record& r, bool primary_only, uint32_t min_len, int min_mapq) {
    if (!primary_only) return true;
    const uint32_t tlen = r.tend - r.tstart;
    const uint32_t qlen = r.qend - r.qstart;
    if (std::max(tlen, qlen) < min_len) return false;
    if (r.mapq < static_cast<uint32_t>(min_mapq)) return false;
    char ty=0; std::string_view v;
    if (!paf::find_tag(r, "tp", ty, v)) return true;
    if (ty!='A' || v.empty()) return true;
    return v[0]=='P';
}

static inline void paf_to_map(const std::string& paf_path, std::ostream& out, bool primary_only, uint32_t min_len, int min_mapq) {
    log_stream() << "Converting PAF to map format: " << paf_path << "\n";

    paf::Reader rd(paf_path);
    paf::Record r;
    while (rd.next(r)) {
        if (!keep_paf_record(r, primary_only, min_len, min_mapq)) continue;  // primary, alignment length and mapq filter
        paf_emit_match_blocks(r, out);
    }
}

static inline void file_to_map_auto(const std::string& in_path, const std::string& out_path, bool paf_primary_only, uint32_t min_len, int min_mapq) {
    std::ostream* out = &std::cout;
    std::ofstream fout;

    if (!out_path.empty() && out_path != "-") {
        fout.open(out_path);
        if (!fout) {
            error_stream() << out_path << ": No such file or directory\n";
            std::exit(1);
        }
        out = &fout;
    }

    if (has_suffix(in_path, ".gfa") || has_suffix(in_path, ".gfa.gz")) {
        gfa_to_map(in_path, *out);
    }
    else {
        paf_to_map(in_path, *out, paf_primary_only, min_len, min_mapq);
    }
}

} // namespace mapconv