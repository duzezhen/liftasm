#include "../include/seed_extend_scores.hpp"

#include <cmath>
#include <algorithm>

namespace seedExtend {

/* ---------------- identity ---------------- */
IdentityMetrics compute_cigar_stats(const std::vector<CIGAR::COp>& ops) {
    uint32_t matches = 0; // M/=
    uint32_t aligned_q = 0; // M+X+I
    uint32_t aligned_r = 0; // M+X+D
    uint32_t full_query = 0; // M+X+I+S+H
    uint32_t ins_bases = 0; // I
    uint32_t del_bases = 0; // D

    for (auto const &o : ops) {
        switch (o.op) {
            case '=':
            case 'M':
                matches     += o.len;
                aligned_q   += o.len;
                aligned_r   += o.len;
                full_query  += o.len;
                break;
            case 'X':
                aligned_q   += o.len;
                aligned_r   += o.len;
                full_query  += o.len;
                break;
            case 'I':
                aligned_q   += o.len;
                full_query  += o.len;
                ins_bases   += o.len;
                break;
            case 'D':
                aligned_r   += o.len;
                del_bases   += o.len;
                break;
            case 'S':
            case 'H':
                full_query  += o.len;
                break;
            default:
                break;
        }
    }

    IdentityMetrics idt;
    idt.matches    = matches;
    idt.mismatches = aligned_q - matches; // M+X - M
    idt.ins_bases  = ins_bases;
    idt.del_bases  = del_bases;
    idt.aln_len    = aligned_q + aligned_r - matches; // M+X+I + M+X+D - M
    idt.full_qry_len = full_query;
    idt.qry_local_identity  = aligned_q > 0 ? double(matches) / aligned_q : 0.0;
    idt.ref_local_identity  = aligned_r > 0 ? double(matches) / aligned_r : 0.0;
    idt.qry_global_identity = full_query > 0 ? double(matches) / full_query : 0.0;
    return idt;
}

/* ---------------- NM ---------------- */
void compute_nm(std::vector<FragAlign>& aligns) {
    if (aligns.empty()) return;

    for (auto& align : aligns) {
        align.tags.NM = align.identity_scores.mismatches + align.identity_scores.ins_bases + align.identity_scores.del_bases;
        if (align.tags.NM < 0) align.tags.NM = 0;
    }
}

/* ---------------- MAPQ ---------------- */
// mapq_raw = I × P × 40 × (1 – (s2 / s1)^2) × ln(s1)
inline int compute_mapq_core(int cm1, int cm2, int cs1, double qry_global_identity) {
    if (cm1 <= 0) return 0;
    if (cm2 < 0)  cm2 = 0;
    if (cm2 > cm1) cm2 = cm1;

    /* long-read penalty terms */
    float pen_s1 = (cm1 > 100 ? 1.0f : 0.01f * cm1);
    float pen_cm = (cs1 > 10 ? 1.0f : 0.1f * cs1);
    if (pen_s1 < pen_cm) pen_cm = pen_s1;

    float x = (float)cm2 / cm1;
    float raw = qry_global_identity * pen_cm * 40.0f * (1.0f - x * x) * std::log((float)cm1);

    int mapq = (raw < 0.0f ? 0 : (int)(raw + 0.5f));
    if (mapq > 60) mapq = 60;
    return mapq;
}

void compute_mapq(std::vector<FragAlign>& aligns) {
    if (aligns.empty()) return;

    for (size_t i = 0; i < aligns.size(); ++i) aligns[i].tags.tp = (i == 0 ? 'P' : 'S');

    int s1 = aligns[0].tags.cm;
    int s2 = (aligns.size() > 1 ? aligns[1].tags.cm : 0);
    int n_sub = (int)aligns.size() - 1;

    /* MAPQ (Primary) */
    double identity = aligns[0].identity_scores.qry_local_identity;
    int mapq = compute_mapq_core(s1, s2, aligns[0].tags.cs, identity);

    /* n_sub: 4.343 * ln(n_sub+1) （≈6.02*log10）*/
    mapq -= (int)(4.343f * std::log((float)n_sub + 1.0f) + 0.499f);
    if (mapq < 0) mapq = 0;

    aligns[0].MAPQ = mapq;
    /* MAPQ (Secondary) */
    for (size_t i = 1; i < aligns.size(); ++i)
        aligns[i].MAPQ = 0;
}

// --- mark supplementary alignments based on overlap with primary alignment ---
inline std::vector<bool> mark_supplementary(std::vector<FragAlign>& aligns, const double sec_pri_ratio) {
    std::vector<bool> is_supp(aligns.size(), false);
    if (aligns.empty()) return is_supp;

    // primary alignment
    aligns[0].tags.tp = 'P';
    int primary_score = aligns[0].tags.cm > 0 ? aligns[0].tags.cm : 1; 

    // secondary and supplementary alignments
    for (size_t i = 1; i < aligns.size(); ++i) {
        double score_ratio = double(aligns[i].tags.cm) / double(primary_score);
        if (score_ratio >= sec_pri_ratio) {
            is_supp[i] = false;
        } else {
            is_supp[i] = true;
        }
        aligns[i].tags.tp = 'S';
    }
    return is_supp;
}

// --- make SA:Z tag for a given index in the aligns vector ---
inline std::string make_SA_tag_for(
    size_t idx,
    const std::vector<FragAlign>& aligns
) {
    if (aligns.size() <= 1) return std::string();

    auto encode_one = [](const FragAlign& fa)->std::string {
        char strand = fa.align_rev ? '-' : '+';
        char buf[256];
        // rname,pos(1-based),strand,CIGAR,mapQ,NM
        std::snprintf(
            buf, sizeof(buf), "%s,%u,%c,%s,%u,%d",
            fa.ref_name.c_str(), fa.r_beg + 1, strand,
            fa.cigar.c_str(),
            (unsigned)(fa.MAPQ < 0 ? 0 : fa.MAPQ),
            fa.tags.NM
        );
        return std::string(buf);
    };

    std::string sa;
    sa.reserve(aligns.size() * 32);
    for (size_t i = 0; i < aligns.size(); ++i) {
        if (i == idx) continue;
        if (!sa.empty()) sa.push_back(';');
        sa += encode_one(aligns[i]);
    }
    return sa;
}

// --- assign flags and SA:Z tags to each FragAlign in the vector ---
void assign_flags_and_SA(
    std::vector<FragAlign>& aligns,
    const double sec_pri_ratio
) {
    if (aligns.empty()) return;

    // primary / secondary / supplementary
    auto is_supp = mark_supplementary(aligns, sec_pri_ratio);

    for (size_t i = 0; i < aligns.size(); ++i) {
        const auto& fa = aligns[i];

        // ---------- FLAG ----------
        uint16_t flag = 0;
        if (fa.align_rev) flag |= F_REVERSE;  // 0x0010
        if (i == 0) {
        } else {
            if (is_supp[i]) flag |= F_SUPPLEMENT; // 0x800
            else            flag |= F_SECONDARY;  // 0x100
        }
        aligns[i].flag = flag;

        // ---------- SA:Z ----------
        aligns[i].tags.SA = make_SA_tag_for(i, aligns);
    }
}

} // namespace seedExtend