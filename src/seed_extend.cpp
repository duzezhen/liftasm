#include "../include/seed_extend.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

using namespace wfa; // WFAlignerGapAffine

namespace seedExtend {

/* ================================================================== */
/* Fragment task (flanks + large gaps)                                */
/* ================================================================== */
static std::vector<FragTask>
build_frag_tasks_from_blocks(
    const mmidx::AnchorChain& ac,
    uint32_t read_len,
    uint32_t contig_len,
    const opt::ExtendOpts& ExtendOpts
) {
    std::vector<FragTask> T;

    bool chain_rev = ac.chain_rev;
    const auto& blocks = ac.blocks;

    if (blocks.empty()) return T;
    const auto& b0 = blocks.front();
    const auto& bN = blocks.back();

    // leading
    if (b0.q_beg > 0) {
        uint32_t q_beg = 0, q_end = b0.q_beg;
        uint32_t r_beg = 0, r_end = 0;
        r_end = b0.r_beg;
        r_beg = (r_end > ExtendOpts.flank_pad ? r_end - ExtendOpts.flank_pad : 0);
        q_beg = (q_end > ExtendOpts.flank_pad ? q_end - ExtendOpts.flank_pad : 0);
        if (q_beg < q_end && r_beg < r_end) T.push_back({r_beg, r_end, q_beg, q_end, F_LEAD, chain_rev});
    }

    // gaps between blocks
    for (size_t i=0; i+1<blocks.size(); ++i) {
        const auto& A = blocks[i];
        const auto& B = blocks[i+1];
        uint32_t rgap = (B.r_beg > A.r_end ? B.r_beg - A.r_end : 0);
        uint32_t qgap = (B.q_beg > A.q_end ? B.q_beg - A.q_end : 0);
        if (rgap > 0 || qgap > 0) {
            T.push_back({A.r_end, B.r_beg, A.q_end, B.q_beg, F_GAP, chain_rev});
        }
    }

    // trailing
    if (bN.q_end < read_len) {
        uint32_t q_beg = bN.q_end;
        uint32_t q_end = read_len;
        uint32_t r_beg, r_end;
        r_beg = bN.r_end;
        r_end = std::min(bN.r_end + ExtendOpts.flank_pad, contig_len);
        q_end = std::min(q_beg + ExtendOpts.flank_pad, read_len);
        if (q_beg < q_end && r_beg < r_end) T.push_back({r_beg, r_end, q_beg, q_end, F_TRAIL, chain_rev});
    }

    return T;
}

/* ================================================================== */
/* Run WFA fragment                                                     */
/* ================================================================== */
static void run_wfa_fragment(
    std::string_view rseq,
    std::string_view qseq,
    std::string& cigar_out,
    int& score_out, 
    const opt::ExtendOpts& ExtendOpts
) {
    WFAlignerGapAffine aligner(ExtendOpts.match, ExtendOpts.gap_open, ExtendOpts.gap_extend, ExtendOpts.mismatch, WFAligner::Alignment, WFAligner::MemoryHigh);
    aligner.alignEnd2End(rseq.data(), int(rseq.size()), qseq.data(), int(qseq.size()));
    int sc = aligner.getAlignmentScore();
    cigar_out = aligner.getCIGAR(true);
    score_out = sc;
}

/* ================================================================== */
/* Build full CIGAR: fragments + anchor blocks (length-correct)       */
/* ================================================================== */
static std::string build_full_cigar(
    const std::vector<mmidx::Anchor>& blocks,
    const std::vector<FragTask>& tasks,
    const std::vector<std::string>& tasks_cigs
) {
    std::vector<CIGAR::COp> ops, tmp;

    auto add_cig = [&](const std::string& cig) {
        if (cig=="*" || cig.empty()) return;
        tmp = CIGAR::parse(cig);
        for (auto const& o: tmp) CIGAR::append(ops, o.len, o.op);
    };
    auto add_block = [&](const mmidx::Anchor& b) {
        uint32_t len_r = b.r_end - b.r_beg;
        uint32_t len_q = b.q_end - b.q_beg;
        uint32_t m = std::min(len_r, len_q);
        uint32_t extra_r = len_r - m;
        uint32_t extra_q = len_q - m;
        if (m) CIGAR::append(ops, m, '=');
        if (extra_r) CIGAR::append(ops, extra_r, 'D');
        if (extra_q) CIGAR::append(ops, extra_q, 'I');
    };

    // gather leading
    size_t ti_lead = SIZE_MAX, ti_trail = SIZE_MAX;
    std::vector<size_t> ti_gap;
    for (size_t i = 0; i < tasks.size(); ++i) {
        if (tasks[i].kind == F_LEAD) ti_lead = i;
        else if (tasks[i].kind == F_TRAIL) ti_trail = i;
        else ti_gap.push_back(i);
    }

    if (ti_lead!=SIZE_MAX) add_cig(tasks_cigs[ti_lead]);

    for (size_t bi = 0; bi < blocks.size(); ++bi) {
        add_block(blocks[bi]); 
        if (bi < ti_gap.size()) add_cig(tasks_cigs[ti_gap[bi]]);
    }

    if (ti_trail!=SIZE_MAX) add_cig(tasks_cigs[ti_trail]);

    return CIGAR::pack(ops);
}

// ================================================================== */
// Compute alignment bounds from tasks and blocks
// ================================================================== */
static AlignmentBounds compute_bounds(
    const std::vector<FragTask>& tasks,
    const std::vector<mmidx::Anchor>& blocks
) {
    AlignmentBounds bounds{0, 0, 0, 0};

    if (!tasks.empty() && tasks.front().kind == F_LEAD) {
        bounds.r_beg = tasks.front().r_beg;
        bounds.q_beg = tasks.front().q_beg;
    } else {
        bounds.r_beg = blocks.front().r_beg;
        bounds.q_beg = blocks.front().q_beg;
    }

    if (!tasks.empty() && tasks.back().kind == F_TRAIL) {
        bounds.r_end = tasks.back().r_end;
        bounds.q_end = tasks.back().q_end;
    } else {
        bounds.r_end = blocks.back().r_end;
        bounds.q_end = blocks.back().q_end;
    }

    return bounds;
}

/* ================================================================== */
// This function trims leading/trailing deletes, clips query ends, and
// packs CIGAR into a string.
/* ================================================================== */
static FragAlign finalize_frag_align(
    AlignmentBounds bounds,
    std::string merged_cigar,
    int total_score,
    const std::string_view ref_name,
    const std::string_view read_name,
    const uint32_t& read_len,
    bool chain_rev,
    const opt::ExtendOpts& ExtendOpts
) {
    std::vector<CIGAR::COp> ops = CIGAR::parse(merged_cigar);
    // trim leading delete
    if (!ops.empty() && ops.front().op == 'D') {
        bounds.r_beg += ops.front().len;
        ops.erase(ops.begin());
    }
    // trim trailing delete
    if (!ops.empty() && ops.back().op == 'D') {
        bounds.r_end -= ops.back().len;
        ops.pop_back();
    }
    // clip query end: leading and trailing
    uint32_t clip_front = bounds.q_beg;
    uint32_t clip_back = read_len - bounds.q_end;
    if (clip_front > 0) CIGAR::prepend(ops, clip_front, ExtendOpts.hard_clip ? 'H' : 'S');
    if (clip_back > 0) CIGAR::append(ops, clip_back, ExtendOpts.hard_clip ? 'H' : 'S');

    std::string cigar = CIGAR::pack(ops);
    IdentityMetrics identity_scores = compute_cigar_stats(ops);

    if (identity_scores.full_qry_len != read_len) {
        error_stream() << read_name <<  ": CIGAR length mismatch: expected " << read_len
            << ", got " << identity_scores.full_qry_len << " in CIGAR: " << cigar << "\n";
        exit(1);
    }

    FragAlign out;
    out.ops       = std::move(ops);
    out.cigar     = std::move(cigar);
    out.tags.AS   = total_score;
    out.r_beg     = bounds.r_beg;
    out.r_end     = bounds.r_end;
    out.q_beg     = bounds.q_beg;
    out.q_end     = bounds.q_end;
    out.ref_name  = ref_name;
    out.qry_name  = read_name;
    out.align_rev = chain_rev;
    out.identity_scores  = identity_scores;

    return out;
}

static DynResResult dynamic_extend_end(
    EndDir dir,
    const std::string_view r_seq,
    const std::string_view q_seq,
    uint32_t anchor_q,               // LEAD: b0.q_beg；TRAIL: bN.q_end
    uint32_t anchor_r,               // LEAD: b0.r_beg；TRAIL: bN.r_end
    const opt::ExtendOpts& extendOpts
) {
    DynResResult out;
    if (!extendOpts.dyn_rescue_enable) return out;

    const uint32_t qlen = (uint32_t)q_seq.size();
    const uint32_t rlen = (uint32_t)r_seq.size();

    uint32_t q_beg=0, q_end=0, r_beg=0, r_end=0;
    if (dir == EndDir::LEAD) {
        q_end = anchor_q; q_beg = q_end;
        r_end = anchor_r; r_beg = r_end;
    } else { // TRAIL
        q_beg = anchor_q; q_end = q_beg;
        r_beg = anchor_r; r_end = r_beg;
    }

    const uint32_t q_min_lead  = (dir == EndDir::LEAD
                                  ? (q_end > extendOpts.dyn_max_query_bp ? q_end - extendOpts.dyn_max_query_bp : 0u)
                                  : q_beg);
    const uint32_t q_max_trail = (dir == EndDir::TRAIL
                                  ? std::min<uint32_t>(qlen, q_beg + extendOpts.dyn_max_query_bp)
                                  : q_end);
    const uint32_t r_min_lead  = (dir == EndDir::LEAD
                                  ? (r_end > extendOpts.dyn_max_ref_bp ? r_end - extendOpts.dyn_max_ref_bp : 0u)
                                  : r_beg);
    const uint32_t r_max_trail = (dir == EndDir::TRAIL
                                  ? std::min<uint32_t>(rlen, r_beg + extendOpts.dyn_max_ref_bp)
                                  : r_end);

    int prev_sc = 0;
    int steps = 0;

    while (steps < (int)extendOpts.dyn_max_steps) {
        uint32_t add = extendOpts.dyn_step_bp;

        if (dir == EndDir::LEAD) {
            const uint32_t q_margin = (q_beg > q_min_lead) ? (q_beg - q_min_lead) : 0u;
            const uint32_t r_margin = (r_beg > r_min_lead) ? (r_beg - r_min_lead) : 0u;
            if (q_margin == 0 || r_margin == 0) break;
            const uint32_t max_back = std::min(q_margin, r_margin);
            if (add > max_back) add = max_back;

            const uint32_t try_qb = q_beg - add;
            const uint32_t try_rb = r_beg - add;

            std::string_view qseg = q_seq.substr(try_qb, q_end - try_qb);
            std::string_view rseg = r_seq.substr(try_rb, r_end - try_rb);
            if (qseg.empty() || rseg.empty()) break;

            std::string cig; int sc = 0;
            run_wfa_fragment(rseg, qseg, cig, sc, extendOpts);

            if (sc - prev_sc < extendOpts.dyn_zdrop) break;

            q_beg = try_qb; r_beg = try_rb;
            out.used = true;
            out.q_beg = q_beg; out.q_end = q_end;
            out.r_beg = r_beg; out.r_end = r_end;
            out.cigar = std::move(cig); out.score = sc;
            prev_sc = sc;

        } else { // TRAIL
            const uint32_t q_margin = (q_max_trail > q_end) ? (q_max_trail - q_end) : 0u;
            const uint32_t r_margin = (r_max_trail > r_end) ? (r_max_trail - r_end) : 0u;
            if (q_margin == 0 || r_margin == 0) break;
            const uint32_t max_fwd = std::min(q_margin, r_margin);
            if (add > max_fwd) add = max_fwd;

            const uint32_t try_qe = q_end + add;
            const uint32_t try_re = r_end + add;

            std::string_view qseg = q_seq.substr(q_beg, try_qe - q_beg);
            std::string_view rseg = r_seq.substr(r_beg, try_re - r_beg);
            if (qseg.empty() || rseg.empty()) break;

            std::string cig; int sc = 0;
            run_wfa_fragment(rseg, qseg, cig, sc, extendOpts);

            if (sc - prev_sc < extendOpts.dyn_zdrop) break;

            q_end = try_qe; r_end = try_re;
            out.used = true;
            out.q_beg = q_beg; out.q_end = q_end;
            out.r_beg = r_beg; out.r_end = r_end;
            out.cigar = std::move(cig); out.score = sc;
            prev_sc = sc;
        }

        ++steps;
    }

    return out;
}

/* ================================================================== */
/* extend_chain_wfa()                                                 */
/* ================================================================== */
std::vector<FragAlign> extend_chain_wfa(
    const std::vector<mmidx::AnchorChain>& anchorChains, 
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs,
    const std::string_view read_name,
    const std::string_view read,
    const opt::ExtendOpts& ExtendOpts
) {
    std::vector<FragAlign> fragAligns;
    fragAligns.reserve(anchorChains.size());
    if (anchorChains.empty()) return fragAligns;

    const uint32_t read_len = read.size();

    uint32_t frag_id = 0;

    for (const auto& anchorChain : anchorChains) {
        if (anchorChain.blocks.empty()) return fragAligns; // empty anchor blocks
    
        // Reference & read info
        const uint32_t r_id = anchorChain.r_id;
        const bool chain_rev = anchorChain.chain_rev;
        const auto& blocks = anchorChain.blocks;
        const std::string ref_name = names[r_id];
        const std::string_view ref_seq = seqs[r_id];
        const uint32_t contig_len = ref_seq.size();

        debug_stream() << "Running WFA fragments for " << read_name << " and " << ref_name << "\n";

        // Tasks building
        std::vector<FragTask> tasks; tasks.reserve(blocks.size() + 2);  // +1 for leading/trailing
        std::vector<int> tasks_scores; tasks_scores.reserve(blocks.size() + 2);
        std::vector<std::string> tasks_cigs; tasks_cigs.reserve(blocks.size() + 2);

        /* ---- RC read (once) ---- */
        std::string read_use(read);
        if (chain_rev) read_use = seqUtils::revcomp(read);

        // dynamic extend leading and trailing
        std::string lead_cig,    trail_cig;
        int         lead_sc = 0, trail_sc = 0;
        uint32_t    lead_qb=0,   lead_rb=0;

        // Leading
        {
            const auto& b0 = blocks.front();
            auto res = dynamic_extend_end(EndDir::LEAD, ref_seq, read_use, b0.q_beg, b0.r_beg, ExtendOpts);
            if (res.used) {
                lead_cig = std::move(res.cigar);
                lead_sc  = res.score;
                lead_qb  = res.q_beg;
                lead_rb  = res.r_beg;
                tasks.push_back({res.r_beg, res.r_end, res.q_beg, res.q_end, F_LEAD, chain_rev});
                tasks_scores.push_back(lead_sc);
                tasks_cigs.push_back(lead_cig);
                debug_stream() << "  - ref=[" << res.r_beg << ", " << res.r_end << ")" << " qry=[" << res.q_beg << ", " << res.q_end << ") " << " CIGAR: " << lead_cig << " score: " << lead_sc << "\n";
            } else {
                lead_cig.clear();
            }
        }
        
        // build gap tasks
        std::vector<FragTask> gap_tasks; gap_tasks.reserve(blocks.size() - 1);
        for (size_t i = 0; i + 1 < blocks.size(); ++i) {
            const auto& A = blocks[i];
            const auto& B = blocks[i + 1];
            uint32_t rgap = (B.r_beg > A.r_end ? B.r_beg - A.r_end : 0);
            uint32_t qgap = (B.q_beg > A.q_end ? B.q_beg - A.q_end : 0);
            if (rgap > 0 || qgap > 0) {
                gap_tasks.push_back({A.r_end, B.r_beg, A.q_end, B.q_beg, F_GAP, chain_rev});
            }
        }
        // run WFA for gaps
        int total_score = 0;
        total_score += lead_sc;

        for (size_t i = 0; i < gap_tasks.size(); ++i) {
            auto& t = gap_tasks[i];
            tasks.push_back(t);

            // Check coordinates
            if (t.r_beg > t.r_end || t.q_beg > t.q_end) {
                error_stream() << "Invalid fragment coordinates for " << read_name << " and " << ref_name << ": "
                    << "ref=[" << t.r_beg << ", " << t.r_end << ") "
                    << "qry=[" << t.q_beg << ", " << t.q_end << ")\n";
                exit(1);
            } else if (t.r_beg == t.r_end && t.q_end > t.q_beg) {  // Insertion
                t.align_success = true;
                tasks_cigs.push_back(std::to_string(t.q_end - t.q_beg) + "I");
                tasks_scores.push_back(ExtendOpts.gap_open + ExtendOpts.gap_extend * (t.q_end - t.q_beg));
                continue;
            } else if (t.q_beg == t.q_end && t.r_end > t.r_beg) {  // Deletion
                t.align_success = true;
                tasks_cigs.push_back(std::to_string(t.r_end - t.r_beg) + "D");
                tasks_scores.push_back(ExtendOpts.gap_open + ExtendOpts.gap_extend * (t.r_end - t.r_beg));
                continue;
            } else if (t.r_beg == t.r_end && t.q_beg == t.q_end) {
                t.align_success = true;
                tasks_cigs.push_back("0=");
                tasks_scores.push_back(0);
                continue;
            }

            std::string_view rseg = ref_seq.substr(t.r_beg, t.r_end - t.r_beg);
            std::string_view qseg(read_use.data() + t.q_beg, t.q_end - t.q_beg);

            std::string cig; int sc = 0;
            run_wfa_fragment(rseg, qseg, cig, sc, ExtendOpts);
            tasks_cigs.push_back(cig);
            tasks_scores.push_back(sc);
            total_score += sc;
            debug_stream() << "  - ref=[" << t.r_beg << ", " << t.r_end << ")" << " qry=[" << t.q_beg << ", " << t.q_end << ") " << " CIGAR: " << cig << " score: " << sc << "\n";
        }

        // Trailing
        {
            const auto& bN = blocks.back();
            auto res = dynamic_extend_end(EndDir::TRAIL, ref_seq, read_use, bN.q_end, bN.r_end, ExtendOpts);
            if (res.used) {
                trail_cig = std::move(res.cigar);
                trail_sc  = res.score;
                tasks.push_back({res.r_beg, res.r_end, res.q_beg, res.q_end, F_TRAIL, chain_rev});
                tasks_scores.push_back(trail_sc);
                tasks_cigs.push_back(trail_cig);
                debug_stream() << "  - ref=[" << res.r_beg << ", " << res.r_end << ")" << " qry=[" << res.q_beg << ", " << res.q_end << ") " << " CIGAR: " << trail_cig << " score: " << trail_sc << "\n";
            } else {
                trail_cig.clear();
            }
        }
        total_score += trail_sc;

        /* ---- build final CIGAR ---- */
        std::string merged_cigar = build_full_cigar(blocks, tasks, tasks_cigs);

        /* ---- compute alignment bounds ---- */
        auto bounds = compute_bounds(tasks, blocks);

        /* ---- finalize fragment alignment ---- */
        FragAlign out = finalize_frag_align(
            bounds,
            merged_cigar,
            total_score,
            ref_name,
            read_name,
            read_len,
            chain_rev,
            ExtendOpts
        );

        out.r_id = r_id;
        // chain score and seed count
        out.tags.cm = anchorChain.score;
        out.tags.cs = anchorChain.seed_number;
        if (frag_id == 0) out.is_primary = true;

        fragAligns.push_back(std::move(out));

        ++frag_id;
    }
    return fragAligns;
}


// Compute alignment scores for each fragment
void cal_align_scores(
    std::vector<FragAlign>& fragAligns, 
    const double sec_pri_ratio
) {
    // sorted by cm and AS
    std::sort(fragAligns.begin(), fragAligns.end(),
        [](const FragAlign& a, const FragAlign& b){
            if (a.tags.cm != b.tags.cm) return a.tags.cm > b.tags.cm;
            return a.tags.AS > b.tags.AS;
    });

    compute_mapq(fragAligns);
    compute_nm(fragAligns);
    assign_flags_and_SA(fragAligns, sec_pri_ratio);  // flags and SA:Z tags
}

void filter_aligns(
    std::vector<FragAlign>& fragAligns,
    const opt::ExtendOpts& ExtendOpts, 
    const int sec_pri_num
) {
    if (fragAligns.empty()) return;

    std::vector<FragAlign> out;
    out.reserve(sec_pri_num);
    int kept_sec = 0;
    for (auto& fa : fragAligns) {
        uint16_t flag = fa.flag;
        if (is_primary(flag)) {
            out.push_back(std::move(fa));
        }
        else if (is_secondary(flag) && kept_sec < sec_pri_num) {
            out.push_back(std::move(fa));
            ++kept_sec;
        }
        if (kept_sec >= sec_pri_num) break;
    }
    fragAligns.swap(out);

    // std::erase_if(fragAligns, [&](auto& fragAlign) {
    //     return !(
    //         fragAlign.identity_scores.qry_local_identity >= ExtendOpts.min_qry_local_identity &&
    //         fragAlign.identity_scores.ref_local_identity >= ExtendOpts.min_ref_local_identity &&
    //         fragAlign.identity_scores.qry_global_identity >= ExtendOpts.min_qry_global_identity
    //     );
    // });
}

/* ------------------------------------------------------------------ */
/* SAM header                                                         */
/* ------------------------------------------------------------------ */
std::string format_sam_header(
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs,
    const std::string& cmd
) {
    std::ostringstream oss;
    oss << "@HD\tVN:1.6\tSO:unsorted\n";
    for (size_t r_id = 0; r_id < names.size(); ++r_id) {
        oss << "@SQ\tSN:" << names[r_id] << "\tLN:" << seqs[r_id].size() << "\n";
    }
    oss << "@PG\tID:" << program::name << "\tPN:" << program::name << "\tVN:" << program::version << "\tCL:" << cmd << "\n";
    return oss.str();
}

/* ------------------------------------------------------------------ */
// Write SAM record for each fragment alignment
/* ------------------------------------------------------------------ */
std::string format_sam_record(
    std::vector<seedExtend::FragAlign>& fragAligns,
    const std::string& read,
    const std::string& qual, 
    const double& sec_pri_ratio
) {
    if (fragAligns.empty()) return "";
    std::ostringstream oss;

    int s1 = fragAligns[0].tags.cm;
    int s2 = (fragAligns.size() > 1 ? fragAligns[1].tags.cm : 0);

    for (size_t i = 0; i < fragAligns.size(); ++i) {
        const auto& fa = fragAligns[i];
        uint16_t flag = fa.flag;

        std::string seq = fa.align_rev ? seqUtils::revcomp(read) : read;
        std::string qv  = fa.align_rev ? std::string(qual.rbegin(), qual.rend()) : qual;

        const auto& ops = fa.ops;
        uint32_t clip_front = (!ops.empty() && ops.front().op == 'H') ? ops.front().len : 0;
        uint32_t clip_back = (!ops.empty() && ops.back().op == 'H') ? ops.back().len : 0;
        if (clip_front) {
            seq = seq.substr(clip_front);
            if (clip_front < qv.size()) qv = qv.substr(clip_front);
        }
        if (clip_back) {
            seq.resize(seq.size() - clip_back);
            if (clip_back < qv.size()) qv.resize(qv.size() - clip_back);
        }
        if (seq.empty()) continue;

        // Signal-end sequencing
        const char* RNEXT = fa.RNEXT.empty() ? "*" : fa.RNEXT.c_str();
        const uint32_t PNEXT = fa.PNEXT;
        const int32_t  TLEN  = fa.TLEN;

        oss << fa.qry_name << '\t'  // QNAME
            << fa.flag     << '\t'  // FLAG
            << fa.ref_name << '\t'  // RNAME
            << (fa.r_beg + 1) << '\t'  // POS (1-based)
            << static_cast<unsigned>(fa.MAPQ) << '\t'  // MAPQ
            << fa.cigar   << '\t'  // CIGAR
            << (fa.RNEXT.empty() ? "*" : fa.RNEXT) << '\t'  // RNEXT
            << fa.PNEXT   << '\t'  // PNEXT
            << fa.TLEN    << '\t'  // TLEN
            << seq        << '\t'  // SEQ
            << qv;
        // tags
        oss << "\tAS:i:" << fa.tags.AS
            << "\tNM:i:" << fa.tags.NM
            << "\tcm:i:" << fa.tags.cm
            << "\tcs:i:" << fa.tags.cs
            << "\ttp:A:" << fa.tags.tp;
        if (i == 0) {
            oss << "\ts1:i:" << s1;
            if (s2 > 0) oss << "\ts2:i:" << s2;
        }
        if (!fa.tags.SA.empty()) {
            oss << "\tSA:Z:" << fa.tags.SA;
        }
        oss << '\n';
    }

    return oss.str();
}

/* ------------------------------------------------------------------ */
// Write PAF record for each fragment alignment
/* ------------------------------------------------------------------ */
static inline void paf_matches_and_blocklen_from_ops(
    const std::vector<CIGAR::COp>& ops,
    uint32_t& matches,
    uint32_t& blockLen
) {
    matches  = 0;
    blockLen = 0;
    bool has_eq_x = false;
    for (const auto& o : ops) {
        if (o.op == '=' || o.op == 'X') { has_eq_x = true; break; }
    }
    for (const auto& o : ops) {
        switch (o.op) {
            case '=': matches += o.len; blockLen += o.len; break;
            case 'X': blockLen += o.len; break;
            case 'M':
                // If no base-resolved ops, treat all 'M' as matches (approx)
                if (!has_eq_x) matches += o.len;
                blockLen += o.len;
                break;
            default: break; // I/D/S/H/N do not contribute to PAF mlen/blen
        }
    }
}

std::string format_paf_record(
    const std::vector<seedExtend::FragAlign>& fragAligns,
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs,
    const std::string& read
) {
    if (fragAligns.empty()) return "";

    std::ostringstream oss;
    const uint32_t qlen = static_cast<uint32_t>(read.size());
    // s1/s2: for primary alignment only (follow your SAM writer's convention)
    int s1 = fragAligns[0].tags.cm;
    int s2 = (fragAligns.size() > 1 ? fragAligns[1].tags.cm : 0);

    for (size_t i = 0; i < fragAligns.size(); ++i) {
        const auto& fa = fragAligns[i];

        // strand and query start/end on the ORIGINAL query strand (PAF rule)
        const char strand = fa.align_rev ? '-' : '+';
        uint32_t qs = fa.q_beg, q_end = fa.q_end;
        if (fa.align_rev) {
            qs = qlen - fa.q_end;  // convert RC coords back to fwd coords
            q_end = qlen - fa.q_beg;
        }

        // target length by name
        const uint32_t tlen = seqs[fa.r_id].size();

        // matches and block length from CIGAR ops
        const std::vector<CIGAR::COp>* pops = nullptr;
        std::vector<CIGAR::COp> tmp_ops;
        if (!fa.ops.empty()) {
            pops = &fa.ops;
        } else {
            tmp_ops = CIGAR::parse(fa.cigar);
            pops = &tmp_ops;
        }
        uint32_t mlen = 0, blen = 0;
        paf_matches_and_blocklen_from_ops(*pops, mlen, blen);

        // NM editing distance
        int nm_val = fa.tags.NM;

        // MAPQ in PAF is 0..255
        unsigned paf_mapq = (fa.MAPQ < 0 ? 0u : static_cast<unsigned>(fa.MAPQ));
        if (paf_mapq > 255) paf_mapq = 255;

        // ---- 12 mandatory PAF fields ----
        // 1: qname, 2: qlen, 3: qstart, 4: qend, 5: strand,
        // 6: tname, 7: tlen, 8: tstart, 9: tend,
        // 10: mlen, 11: blen, 12: mapq
        oss << fa.qry_name << '\t'
            << qlen << '\t'
            << qs   << '\t'
            << q_end   << '\t'
            << strand << '\t'
            << fa.ref_name << '\t'
            << tlen << '\t'
            << fa.r_beg << '\t'
            << fa.r_end << '\t'
            << mlen << '\t'
            << blen << '\t'
            << paf_mapq;

        // ---- Optional tags (minimap2-like) ----
        char tp = fa.tags.tp ? fa.tags.tp : (i == 0 ? 'P' : 'S');
        oss << "\ttp:A:" << tp;
        oss << "\tcm:i:" << fa.tags.cm;
        oss << "\tcs:i:" << fa.tags.cs;
        oss << "\tNM:i:" << fa.tags.NM;
        oss << "\tAS:i:" << fa.tags.AS;
        if (tp == 'P') {
            oss << "\ts1:i:" << s1;
            if (s2 > 0) oss << "\ts2:i:" << s2;
        }
        if (!fa.cigar.empty() && fa.cigar != "*") oss << "\tcg:Z:" << fa.cigar;
        oss << '\n';
    }
    return oss.str();
}

} // namespace seedExtend