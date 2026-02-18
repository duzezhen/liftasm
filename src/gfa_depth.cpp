#include <cassert>
#include <iostream>
#include <zlib.h>

#include "../include/progress_tracker.hpp"
#include "../include/logger.hpp"
#include "../include/gfa_depth.hpp"
#include "../include/kseq.h"
#include "../include/ThreadPool.hpp"
#include "../include/save.hpp"
#include "../include/options.hpp"

KSEQ_INIT(gzFile, gzread)


void GfaDepth::build_segment_offsets_(
    std::vector<uint32_t>& seg_offsets,
    uint64_t& total_len
) {
    const size_t nseg = nodes_.size();
    seg_offsets.assign(nseg + 1, 0);
    total_len = 0;
    for (size_t i = 0; i < nseg; ++i) {
        total_len += nodes_[i].length;
        seg_offsets[i + 1] = seg_offsets[i] + nodes_[i].length;
    }
}

void GfaDepth::fill_pos_depth_from_index_(
    const mmidx::MinimizerIndex& GIndex,
    const std::vector<uint32_t>& seg_offsets,
    std::vector<uint16_t>& pos_depth,
    std::vector<uint32_t>& pos_freq
) {
    const size_t nseg = nodes_.size();

    const auto& keys       = GIndex.keys();
    const auto& offs       = GIndex.offs();
    const auto& positions  = GIndex.positions();
    const auto& depths     = GIndex.depths();

    ProgressTracker tracker(keys.size()); // print every 10% internally
    for (uint64_t i = 0; i < keys.size(); ++i) {
        tracker.hit();  // update progress bar / print at intervals

        const uint16_t depth = (i < depths.size() ? depths[i] : 0u);
        const uint32_t beg   = offs[i];
        const uint32_t end   = offs[i + 1];
        if (beg >= end) continue;

        // All occurrences of this key; (end - beg) is the multiplicity of the key
        const uint32_t mult = end - beg;

        for (uint32_t j = beg; j < end; ++j) {
            const uint32_t seg_id = mmidx::seq_id(positions[j]);
            const uint32_t off    = mmidx::offset(positions[j]);
            if (seg_id >= nseg) continue;

            const uint32_t abs_pos = seg_offsets[seg_id] + off;
            if (abs_pos >= pos_depth.size()) continue;

            pos_depth[abs_pos] = depth;
            pos_freq[abs_pos]  = mult;
        }
    }
}

void GfaDepth::compute_segment_avg_depth_(
    const std::vector<uint32_t>& seg_offsets,
    const std::vector<uint16_t>& pos_depth,
    const std::vector<uint32_t>& pos_freq,
    std::vector<float>& seg_avg_depth
) {
    const size_t nseg = nodes_.size();
    seg_avg_depth.assign(nseg, 0.0f);

    for (size_t s = 0; s < nseg; ++s) {
        const auto& node = nodes_[s];
        const uint32_t len = node.length;
        if (len == 0) continue;

        uint64_t start = seg_offsets[s];
        uint64_t end   = seg_offsets[s + 1];
        if (end > pos_depth.size()) end = pos_depth.size();

        // Sum depth normalized by multiplicity, so multi-hit k-mers contribute less
        long double sum_depth = 0.0L;

        for (uint64_t p = start; p < end; ++p) {
            const uint32_t freq = pos_freq[p];
            if (freq == 0) continue; // no k-mer mapped here
            sum_depth += static_cast<long double>(pos_depth[p]) / static_cast<long double>(freq);
        }

        // Average per base of the segment (including positions w/o k-mers as 0)
        seg_avg_depth[s] = static_cast<float>(sum_depth / static_cast<long double>(len));
    }
}

void GfaDepth::write_kmer_depth_bed_(
    const std::string& out_file,
    const std::vector<uint32_t>& seg_offsets,
    const std::vector<uint16_t>& pos_depth,
    const std::vector<uint32_t>& pos_freq,
    const std::vector<float>&    seg_avg_depth
) {
    SAVE saver(out_file, 10 * 1024 * 1024);  // 10MB buffer
    std::string line; line.reserve(1 * 1024 * 1024);

    // Header: average depth section
    saver.save("## liftasm depth (from k-mer)\n");
    saver.save("## Segment\tLength\tAvgDepth\n");
    for (size_t s = 0; s < nodes_.size(); ++s) {
        line.clear();
        line += "# ";
        line += nodes_[s].name;                line += '\t';
        line += std::to_string(nodes_[s].length); line += '\t';

        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", seg_avg_depth[s]);
        line += buf; line += '\n';
        saver.save(line);
    }

    // Base-level section: chrom  pos  depth  freq
    saver.save("## chrom\tpos\tdepth\tfreq\n");
    for (size_t s = 0; s < nodes_.size(); ++s) {
        const auto& node = nodes_[s];
        if (node.length == 0) continue;

        const uint64_t start = seg_offsets[s];
        uint64_t end         = seg_offsets[s + 1];
        if (end > pos_depth.size()) end = pos_depth.size();

        for (uint64_t p = start; p < end; ++p) {
            line.clear();
            line += node.name;                         line += '\t';
            line += std::to_string(p - start);         line += '\t';
            line += std::to_string(pos_depth[p]);      line += '\t';
            line += std::to_string(pos_freq[p]);       line += '\n';
            saver.save(line);
        }
    }
}

bool GfaDepth::count_from_kmer(
    const mmidx::MinimizerIndex& GIndex,
    const std::string& out_file
) {
    log_stream() << "Computing node depth from k-mer index ..." << "\n";

    if (GIndex.size() == 0) {
        error_stream() << "k-mer index is empty. Please check the input GFA file for correctness.\n";
        std::exit(1);
    }

    // 1) Build segment offsets
    std::vector<uint32_t> seg_offsets;
    uint64_t total_len = 0;
    build_segment_offsets_(seg_offsets, total_len);

    // 2) Allocate base-level arrays and fill from index
    std::vector<uint16_t> pos_depth(total_len, 0);
    std::vector<uint32_t> pos_freq(total_len, 0);
    fill_pos_depth_from_index_(
        GIndex,
        seg_offsets,
        pos_depth,
        pos_freq
    );

    // 3) Compute per-segment average depth
    std::vector<float> seg_avg_depth;
    compute_segment_avg_depth_(
        seg_offsets,
        pos_depth,
        pos_freq,
        seg_avg_depth
    );

    // 4) Output BED file
    write_kmer_depth_bed_(out_file, seg_offsets, pos_depth, pos_freq, seg_avg_depth);

    return true;
}


bool GfaDepth::pass_gaf_filters_(const kgaf::Record& rec) {
    if (min_mapq_ > 0) {
        if (rec.mapq < min_mapq_) return false;
    }

    if (min_align_frac_ > 0.0) {
        const uint64_t qlen  = rec.qlen;
        const uint64_t qs    = rec.qstart;
        const uint64_t qe    = rec.qend;
        if (qlen == 0) return false;
        const uint64_t aln   = (qe > qs ? (qe - qs) : 0);
        const double   frac  = static_cast<double>(aln) / static_cast<double>(qlen);
        if (frac < min_align_frac_) return false;
    }

    return true;
}

bool GfaDepth::build_concatenated_path_(
    const std::vector<kgaf::PathStep>& steps,
    std::vector<Step>& P,
    uint64_t& path_len
) {
    P.clear(); path_len = 0;
    if (steps.empty()) return false;
    P.reserve(steps.size());
    for (const auto& st : steps) {
        auto it = name_to_id_map_.find(std::string(st.name));
        if (it == name_to_id_map_.end()) return false;
        const uint32_t sid  = (uint32_t)it->second;
        const uint32_t slen = nodes_[sid].length;
        P.push_back(Step{sid, slen, path_len, st.rev});
        path_len += slen;
    }
    return !P.empty();
}

bool GfaDepth::accumulate_cigar_(
    const kgaf::Record& rec,
    const std::vector<kgaf::CigarOp>& ops,
    const std::vector<Step>& P,
    const uint64_t path_len,
    const std::vector<uint32_t>& seg_offsets,
    std::vector<uint8_t>& coverage
) {
    if (P.empty() || ops.empty()) return false;

    uint64_t rpos = rec.tstart;
    if (rpos >= path_len) return false;

    int si = -1;
    for (int i = 0; i < (int)P.size(); ++i) {
        if (rpos < P[i].path_beg + P[i].seg_len) { si = i; break; }
    }
    if (si < 0) return false;

    uint32_t seg_id   = P[si].seg_id;
    uint32_t seg_len  = P[si].seg_len;
    uint64_t step_beg = P[si].path_beg;
    bool     rev      = P[si].rev;

    // rpos_in_step: reference position within current step
    uint32_t rpos_in_step = (uint32_t)(rpos - step_beg);

    // small helpers
    auto next_step = [&]() -> bool {
        ++si;
        if (si >= 0 && si < (int)P.size()) {
            seg_id   = P[si].seg_id;
            seg_len  = P[si].seg_len;
            step_beg = P[si].path_beg;
            rev      = P[si].rev;
            rpos_in_step = 0;
            return true;
        }
        return false; // end of path
    };

    // 5) walk the CIGAR on the concatenated path
    for (const auto& op : ops) {
        char c = op.op;
        uint32_t l = op.len;

        if (c=='M' || c=='=' || c=='X') {
            while (l > 0 && si >= 0 && si < (int)P.size()) {
                if (rpos_in_step >= seg_len) { if (!next_step()) break; continue; }
                uint32_t run = std::min<uint32_t>(l, seg_len - rpos_in_step);

                // map path-local positions to segment coordinates (consider reverse)
                if (!rev) {
                    // forward: seg_pos increases with rpos_in_step
                    uint32_t abs = seg_offsets[seg_id] + rpos_in_step;
                    uint8_t* p = coverage.data() + abs;
                    for (uint32_t k = 0; k < run; ++k) { uint8_t& v = p[k]; v = (v == 255 ? 255 : uint8_t(v + 1)); }
                } else {
                    // reverse: seg_pos decreases with rpos_in_step
                    // start index = (seg_len - 1) - rpos_in_step
                    int64_t idx = int64_t(seg_offsets[seg_id]) + int64_t(seg_len - 1 - rpos_in_step);
                    uint8_t* p = coverage.data() + idx;
                    for (uint32_t k = 0; k < run; ++k) { uint8_t& v = *p; v = (v == 255 ? 255 : uint8_t(v + 1)); --p; }
                }

                rpos += run; rpos_in_step += run; l -= run;
            }
        } else if (c=='D' || c=='N') {
            // reference-consuming only: just advance along path
            uint32_t remain = l;
            while (remain > 0 && si >= 0 && si < (int)P.size()) {
                if (rpos_in_step >= seg_len) { if (!next_step()) break; continue; }
                uint32_t adv = std::min<uint32_t>(remain, seg_len - rpos_in_step);
                rpos += adv; rpos_in_step += adv; remain -= adv;
            }
        } else {
            // I/S/H/P etc: no reference consume -> ignore
            continue;
        }

        if (si < 0 || si >= (int)P.size()) break; // ran off the path
    }
    return true;
}

void GfaDepth::write_gaf_depth_bed_(
    const std::string& out_bed,
    const std::vector<uint32_t>& seg_offsets,
    const std::vector<uint8_t>& coverage
) {
    const size_t n = nodes_.size();
    SAVE saver(out_bed, 10 * 1024 * 1024);
    std::string line; line.reserve(256);

    // header like k-mer version
    saver.save("## liftasm depth (from GAF)\n");
    saver.save("## Segment\tLength\tAvgDepth\n");

    // average depth per segment
    for (size_t s = 0; s < n; ++s) {
        uint32_t len = nodes_[s].length;
        uint64_t sum = 0;
        if (len) {
            const uint8_t* p = coverage.data() + seg_offsets[s];
            for (uint32_t i = 0; i < len; ++i) sum += p[i];
        }
        double avg = (len ? (double)sum / (double)len : 0.0);
        char buf[32]; std::snprintf(buf, sizeof(buf), "%.2f", avg);

        line.clear();
        line += "# ";
        line += nodes_[s].name;         line += '\t';
        line += std::to_string(len);    line += '\t';
        line += buf;                    line += '\n';
        saver.save(line);
    }

    // run-length coverage per segment: chrom  start  depth
    saver.save("## chrom\tpos\tdepth\n");
    for (size_t seg_id = 0; seg_id < n; ++seg_id) {
        const std::string& chrom = nodes_[seg_id].name;
        const uint8_t* base = coverage.data() + seg_offsets[seg_id];
        uint32_t seg_len = nodes_[seg_id].length;

        uint32_t start = 0;
        while (start < seg_len) {
            uint8_t d = base[start];

            line.clear();
            line += chrom;                  line += '\t';
            line += std::to_string(start);  line += '\t';
            line += std::to_string(int(d)); line += '\n';
            saver.save(line);

            start++;
        }
    }
}

bool GfaDepth::count_from_gaf(const std::string& gaf_path, const std::string& out_bed)
{
    log_stream() << "Computing node depth from GAF alignments ..." << "\n";

    // Build segment offsets
    const size_t n = nodes_.size();
    std::vector<uint32_t> seg_offsets;
    uint64_t total_len = 0;
    build_segment_offsets_(seg_offsets, total_len);
    std::vector<uint8_t> coverage(seg_offsets.back(), 0);
    

    // ---- open GAF ----
    kgaf::ReaderAuto reader(gaf_path);
    if (!reader.good()) {
        error_stream() << gaf_path << ": No such file or directory\n";
        std::exit(1);
    }

    auto tracker = ProgressTracker::Every(30000);  // print every 30,000 reads

    kgaf::Record                rec;
    std::vector<kgaf::PathStep> steps; steps.reserve(16);
    std::vector<kgaf::CigarOp>  ops;   ops.reserve(64);
    std::vector<Step>           P;     P.reserve(16);

    uint64_t n_rec = 0;
    while (reader.next(rec)) {
        tracker.hit();  // update progress
        ++n_rec;

        // filter out low-quality alignments
        if (!pass_gaf_filters_(rec)) {continue;}
        if (rec.path.empty() || rec.cg.empty()) continue;

        // 1) parse path: each step name is used as segment name AS-IS (no ow)
        steps.clear();
        kgaf::parse_path(rec.path, steps);

        uint64_t path_len = 0;
        if (!build_concatenated_path_(steps, P, path_len)) continue;
        if (rec.tstart >= path_len) continue;

        // 2) parse CIGAR
        ops.clear();
        kgaf::parse_cigar(rec.cg, ops);
        if (ops.empty()) continue;

        // 3) accumulate CIGAR ops onto segments in P
        accumulate_cigar_(rec, ops, P, path_len, seg_offsets, coverage);
    }
    tracker.finish();  // finish progress

    // ---- write BED (avg depth + run-length intervals) ----
    write_gaf_depth_bed_(out_bed, seg_offsets, coverage);

    return true;
}