#pragma once
#include <vector>
#include <string>
#include <string_view>

#include "seed_extend_types.hpp"
#include "CIGAR.hpp"

namespace seedExtend {

IdentityMetrics compute_cigar_stats(const std::vector<CIGAR::COp>& ops);

void compute_nm(std::vector<FragAlign>& aligns);

void compute_mapq(std::vector<FragAlign>& fragAligns);

// --- assign flags and SA:Z tags to each FragAlign in the vector ---
void assign_flags_and_SA(
    std::vector<FragAlign>& aligns,
    const double sec_pri_ratio
);

} // namespace seedExtend