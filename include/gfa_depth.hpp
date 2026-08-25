#pragma once
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <thread>
#include <utility>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <numeric>

#include "CommandOptions.hpp"
#include "gfa_parser.hpp"
#include "../include/kgaf.hpp"
#include "mmidx.hpp"


class GfaDepth : public GfaGraph {
public:
    explicit GfaDepth(DepthOpts params) : params_(std::move(params)) {
        set_forbid_overlap(params_.forbid_overlap);
    }

public:
    bool count_from_kmer(const mmidx::MinimizerIndex& GIndex, const std::string& out_file);

    bool count_from_gaf(const std::string& gaf_path, const std::string& out_file);

    bool count_from_A(const std::string& out_file);
    
protected:
    const DepthOpts params_;

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
        const std::string& out_file,
        const std::vector<uint32_t>& seg_offsets,
        const std::vector<uint8_t>& coverage
    );

    void write_A_depth_bed_(
        const std::string& out_file,
        const std::vector<uint32_t>& seg_offsets,
        const std::vector<uint32_t>& depth
    );
};
