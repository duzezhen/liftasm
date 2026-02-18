#pragma once
/*
 * gfa_bam_loader.hpp  – BAM ➜ GFA alignment, depth & coverage loader
 * ----------------------------------------------------------------------------------------
 *  • Works with `GfaGraph` (friend‑class)
 *  • Provides:
 *       1. Per‑segment depth   (avg depth, aligned bases, read count)
 *       2. Per‑base  coverage  (uint8_t, flattened array)
 *       3. Optional multi‑thread load leveraging BAM index (BAI)
 *
 *  Usage (single thread):
 *     GfaBamLoader ld(graph);
 *     ld.loadFromBam("aln.bam");
 *     ld.save2File("depth.bed");
 *
 */

#include "gfa_parser.hpp"
#include "save.hpp"
#include <htslib/sam.h>
#include <htslib/hts.h>

#include <sstream>
#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <iomanip>
#include <limits>
#include <numeric>
#include <algorithm>

/*============================================================*/
/*            stringify ALL BAM auxiliary tags (debug)        */
/*============================================================*/
std::string stringify_all_aux_tags(const bam1_t* rec);

/*============================================================*/
/*                          GfaBamLoader                      */
/*============================================================*/
class GfaBamLoader {
public:
    explicit GfaBamLoader(GfaGraph& g);
    ~GfaBamLoader();

    /*------------------------------------------------------------*/
    /*              Single‑thread loader (fast path)              */
    /*------------------------------------------------------------*/
    bool loadFromBam(const std::string& bam_path, bool keep_unmapped = false, int io_threads = 4);

    /*------------------------------------------------------------*/
    /*                            Save                            */
    /*------------------------------------------------------------*/
    /**
     * Saves the depth profile to a BED format file.
     * @param out_put_file The output file path.
     * @return true if the save operation was successful, false otherwise.
     */
    bool save2BED(const std::string out_put_file);

    /*--------------------------------------------------------*/
    /*                  Accessors / helpers                   */
    /*--------------------------------------------------------*/
    std::vector<double> getAverageDepths() const;

    void printDepthSummary(size_t max_rows = 10) const;

    void dumpNodeCoverage(uint32_t seg_id, uint32_t first = 1, uint32_t last = 100) const;

private:
    /* helper – compute one BAM record, update shared arrays (single‑thread safe) */
    inline void compute_record_stats(bam1_t* rec, uint32_t seg_id, uint8_t* cov_base);

    GfaGraph& graph_;
    std::vector<uint64_t> depth_total_bases_; // Σ aligned ref bases
    std::vector<uint64_t> depth_read_count_; // primary alignments
    std::vector<uint32_t> offsets_; // prefix sum of segment lengths
    std::vector<uint8_t>  coverage_; // uint8 (0-255) per base
};  // GfaBamLoader