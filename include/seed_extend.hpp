#pragma once
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include "mmidx.hpp"
#include "logger.hpp"

#include "seed_extend_types.hpp"
#include "seed_extend_scores.hpp"
#include "CIGAR.hpp"
#include "seq_utils.hpp"
#include "options.hpp"
#include "ProgramMetadata.hpp"
#include "bindings/cpp/WFAligner.hpp"


namespace seedExtend {

/* --------------------------------------------------------------- *
 * extend_chain_wfa()
 *
 * Parameters:
 *  • anchorChains        Vector of anchor chains to extend.
 *  • names               Vector of reference contig names.
 *  • seqs                Vector of reference contig sequences.
 *  • read_name           Read ID used in SAM.
 *  • read                Full read sequence (forward).
 *  • ExtendOpts        Extend parameters (gap, flank, zdrop).
 *
 *
 *  Return: FragAlign for the whole chain (CIGAR merged & unclipped)
 * --------------------------------------------------------------- */
std::vector<FragAlign> extend_chain_wfa(
    const std::vector<mmidx::AnchorChain>& anchorChains, 
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs,
    const std::string_view read_name,
    const std::string_view read,
    const opt::ExtendOpts& ExtendOpts
);

// Compute alignment scores for each fragment
void cal_align_scores(
    std::vector<FragAlign>& fragAligns, 
    const double sec_pri_ratio
);

// Filter aligns based on options
void filter_aligns(
    std::vector<FragAlign>& fragAligns,
    const opt::ExtendOpts& ExtendOpts, 
    const int sec_pri_num
);


/* ------------ SAM output helpers (text SAM, HTSlib friendly) -------- */
// write @HD/@SQ lines
std::string format_sam_header(
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs, 
    const std::string& cmd
);

// write one alignment line (FLAG set 0x10 if rev-strand)
std::string format_sam_record(
    std::vector<seedExtend::FragAlign>& fragAligns,
    const std::string& read, 
    const std::string& qual, 
    const double& sec_pri_ratio
);

// Write PAF record for each fragment alignment
std::string format_paf_record(
    const std::vector<seedExtend::FragAlign>& fragAligns,
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs,
    const std::string& read
);

} // namespace seedExtend