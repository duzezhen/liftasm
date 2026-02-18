#pragma once
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <numeric>

#include "gfa_parser.hpp"
#include "../include/kgaf.hpp"
#include "mmidx.hpp"


class GfaDepth : public GfaGraph {
public:
    GfaDepth(int min_mapq, double min_align_frac) {
        set_forbid_overlap(true);
        min_mapq_ = min_mapq;
        min_align_frac_ = min_align_frac;
    }

public:
    bool count_from_kmer(const mmidx::MinimizerIndex& GIndex, const std::string& out_file);

    bool count_from_gaf(const std::string& gaf_path, const std::string& out_bed);
    
protected:
    // filters for GAF alignments
    int    min_mapq_        = 0;     // minimum mapping quality (MAPQ) to consider an alignment
    double min_align_frac_  = 0.0;   // minimum aligned fraction (aln_len / qlen) to consider an alignment

protected:
    struct ReadRec { std::string name; std::string seq; std::string qual; };
    struct Step { uint32_t seg_id; uint32_t seg_len; uint64_t path_beg; bool rev; };

    void build_segment_offsets_(
        std::vector<uint32_t>& seg_offsets,
        uint64_t& total_len
    );


    void fill_pos_depth_from_index_(
        const mmidx::MinimizerIndex& GIndex,
        const std::vector<uint32_t>& seg_offsets,
        std::vector<uint16_t>& pos_depth,
        std::vector<uint32_t>& pos_freq
    );
    void compute_segment_avg_depth_(
        const std::vector<uint32_t>& seg_offsets,
        const std::vector<uint16_t>& pos_depth,
        const std::vector<uint32_t>& pos_freq,
        std::vector<float>& seg_avg_depth
    );
    void write_kmer_depth_bed_(
        const std::string& out_file,
        const std::vector<uint32_t>& seg_offsets,
        const std::vector<uint16_t>& pos_depth,
        const std::vector<uint32_t>& pos_freq,
        const std::vector<float>&    seg_avg_depth
    );

    bool pass_gaf_filters_(const kgaf::Record& rec);
    bool build_concatenated_path_(
        const std::vector<kgaf::PathStep>& steps,
        std::vector<Step>& P,
        uint64_t& path_len
    );
    bool accumulate_cigar_(
        const kgaf::Record& rec,
        const std::vector<kgaf::CigarOp>& ops,
        const std::vector<Step>& P,
        const uint64_t path_len,
        const std::vector<uint32_t>& seg_offsets,
        std::vector<uint8_t>& coverage
    );
    void write_gaf_depth_bed_(
        const std::string& out_bed,
        const std::vector<uint32_t>& seg_offsets,
        const std::vector<uint8_t>& coverage
    );
};
