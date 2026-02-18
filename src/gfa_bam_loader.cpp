#include "../include/gfa_bam_loader.hpp"
#include <htslib/hts.h>
#include <iostream>
#include <thread>
#include <atomic>
#include <algorithm>
#include <iomanip>

/*============================================================*/
/*             stringify ALL BAM auxiliary tags               */
/*============================================================*/
std::string stringify_all_aux_tags(const bam1_t* rec)
{
    std::ostringstream oss;
    const uint8_t* ptr = bam_aux_first(const_cast<bam1_t*>(rec));
    bool first = true;
    while (ptr) {
        if (!first) oss << '\t';
        first = false;
        oss << ptr[0] << ptr[1] << ':' << char(ptr[2]) << ':';
        switch (ptr[2]) {
            case 'A': oss << char(bam_aux2A(ptr)); break;
            case 'c': case 'C': case 's': case 'S': case 'i': case 'I': oss << bam_aux2i(ptr); break;
            case 'f': case 'd': oss << bam_aux2f(ptr); break;
            case 'Z': case 'H': oss << bam_aux2Z(ptr); break;
            default : oss << '?'; break;
        }
        ptr = bam_aux_next(const_cast<bam1_t*>(rec), ptr);
    }
    return oss.str();
}

/*============================================================*/
/*                    GfaBamLoader ctor/dtor                  */
/*============================================================*/
GfaBamLoader::GfaBamLoader(GfaGraph& g) : graph_(g)
{
    size_t n = graph_.nodes_.size();
    depth_total_bases_.assign(n, 0ull);
    depth_read_count_.assign(n, 0ull);

    // ---- compute offsets & allocate one contiguous coverage array ----
    offsets_.clear(); offsets_.resize(n + 1, 0);
    for (size_t i = 0; i < n; ++i) offsets_[i + 1] = offsets_[i] + graph_.nodes_[i].length;
    coverage_.resize(offsets_.back(), 0); // uint8_t default‑init to 0
}

GfaBamLoader::~GfaBamLoader() {}

/*============================================================*/
/*               internal helper: per‑record stats            */
/*============================================================*/
void GfaBamLoader::compute_record_stats(bam1_t* rec, uint32_t seg_id, uint8_t* cov_base)
{
    uint32_t* cigar   = bam_get_cigar(rec);
    uint32_t  n_cigar = rec->core.n_cigar;
    uint64_t  aligned_ref_bases = 0;

    uint32_t ref_pos = rec->core.pos;               // 0‑based
    uint8_t* cov_ptr = cov_base + offsets_[seg_id] + ref_pos;
    uint32_t seg_len = graph_.nodes_[seg_id].length;

    for (uint32_t k = 0; k < n_cigar; ++k) {
        int op  = bam_cigar_op(cigar[k]);
        int len = bam_cigar_oplen(cigar[k]);
        switch (op) {
        case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF: {
            aligned_ref_bases += len;
            // increment coverage, saturate at 255
            uint8_t* p = cov_ptr;
            uint8_t* q = p + std::min<uint32_t>(len, seg_len - ref_pos);
            while (p < q) { *p = (*p == 255 ? 255 : *p + 1); ++p; }
            ref_pos += len; cov_ptr += len;
            break; }
        case BAM_CDEL: case BAM_CREF_SKIP: // &1->read‑consuming, &2->ref‑consuming, &4->alignment-consuming
            ref_pos += len; cov_ptr += len; break;
        default: break; // I,S,H,P - ignore
        }
    }
    depth_total_bases_[seg_id] += aligned_ref_bases;
    depth_read_count_[seg_id]  += 1;
}

/*============================================================*/
/*              Public: single‑thread loader                   */
/*============================================================*/
bool GfaBamLoader::loadFromBam(
    const std::string& bam_path,
    bool keep_unmapped,
    int io_threads
) {
    log_stream() << "Loading BAM from '" << bam_path << "'\n";

    htsFile* fp = sam_open(bam_path.c_str(), "r");
    if (!fp) { error_stream() << "cannot open '" << bam_path << "'\n"; return false; }
    hts_set_threads(fp, io_threads); // multi‑thread BGZF decompression
    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) { sam_close(fp); error_stream() << "header read failed\n"; return false; }

    bam1_t* rec = bam_init1();
    uint64_t n_loaded = 0, n_skip_ref = 0, n_unmapped = 0;

    // direct pointer to coverage array for faster math
    uint8_t* cov_base = coverage_.data();

    while (sam_read1(fp, hdr, rec) >= 0) {
        uint16_t flag = rec->core.flag;
        if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue; // skip secondary/supplementary alignments
        if (rec->core.tid < 0) { ++n_unmapped; if (!keep_unmapped) continue; } // unmapped reads

        // ---- Find segment internal id ----
        const char* refname = (rec->core.tid >= 0 ? hdr->target_name[rec->core.tid] : "*");
        auto it = graph_.name_to_id_map_.find(refname);
        if (it == graph_.name_to_id_map_.end()) { ++n_skip_ref; continue; }
        uint32_t seg_id = static_cast<uint32_t>(it->second);

        compute_record_stats(rec, seg_id, cov_base);
        ++n_loaded;
        if (n_loaded % 100000 == 0) {
            log_stream() << "   - Loaded " << n_loaded << " alignments\n";
        }
    }

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(fp);

    log_stream() << "BAM summary: total loaded: " << n_loaded << "; unmapped: " << n_unmapped << "; missing ref: " << n_skip_ref << "\n";

    return true;
}

/*------------------------------------------------------------*/
/*                            Save                            */
/*------------------------------------------------------------*/
/**
 * @brief Saves the depth profile to a file in BED format.
 * 
 * @param out_put_file The output file path.
 * @return true if the save operation was successful, false otherwise.
 */
bool GfaBamLoader::save2BED(const std::string out_put_file) {
    log_stream() << "Wrote depth profile (bed format) to '" << out_put_file << "'\n";

    SAVE saver(out_put_file, 10 * 1024 * 1024);   // 10 MB cache

    std::string line;
    line.reserve(128);

    /* header */
    saver.save("## gfa_bam_loader v0.1\n");

    auto avgDepths = getAverageDepths();
    saver.save("## Segment\tReads\tDepth\n");
    for (size_t i = 0; i < graph_.nodes_.size(); ++i) {
        line.clear();
        line += "# ";
        line += graph_.nodes_[i].name;
        line += '\t';
        line += std::to_string(depth_read_count_[i]);
        line += '\t';

        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", avgDepths[i]);
        line += buf + std::string("\n");

        saver.save(line);
    }

    /* data lines */
    for (size_t seg_id = 0; seg_id < graph_.nodes_.size(); ++seg_id) {
        const std::string& chrom = graph_.nodes_[seg_id].name;
        const uint8_t* cov = coverage_.data() + offsets_[seg_id];
        uint32_t seg_len  = graph_.nodes_[seg_id].length;

        uint32_t start = 0;
        while (start < seg_len) {
            line.clear();
            uint8_t depth = cov[start];
            uint32_t end  = start;
            while (end + 1 < seg_len && cov[end + 1] == depth) ++end;

            line = chrom + '\t' +
                   std::to_string(start) + '\t' +
                   std::to_string(end)   + '\t' +
                   std::to_string(int(depth)) + '\n';
            saver.save(line);
            start = end + 1;
        }
    }
    return true;
}



/*============================================================*/
/*            pretty printing helper (header inline)          */
/*============================================================*/
std::vector<double> GfaBamLoader::getAverageDepths() const
{
    std::vector<double> out(graph_.nodes_.size(), 0.0);
    for (size_t i = 0; i < graph_.nodes_.size(); ++i) {
        uint32_t len = graph_.nodes_[i].length;
        if (len) out[i] = static_cast<double>(depth_total_bases_[i]) / len;
    }
    return out;
}

void GfaBamLoader::printDepthSummary(size_t max_rows) const
{
    std::cout << "\n==== Segment depth summary ====\n";
    auto depths = getAverageDepths();
    for (size_t i = 0; i < std::min(max_rows, depths.size()); ++i) {
        const auto& node = graph_.nodes_[i];
        std::cout << i << '\t' << node.name
                  << "\tlen=" << node.length
                  << "\treads=" << depth_read_count_[i]
                  << "\tdepth=" << std::fixed << std::setprecision(2) << depths[i] << '\n';
    }
    std::cout << "================================\n";
}

void GfaBamLoader::dumpNodeCoverage(uint32_t seg_id, uint32_t first, uint32_t last) const
{
    if (seg_id >= graph_.nodes_.size()) return;
    uint32_t len = graph_.nodes_[seg_id].length;
    first = std::max(first, 1u);
    last  = std::min(last, len);
    const uint8_t* p = coverage_.data() + offsets_[seg_id] + first - 1;
    for (uint32_t pos = first; pos <= last; ++pos, ++p)
        std::cout << pos << '\t' << int(*p) << '\n';
}