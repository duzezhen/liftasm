#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace gapfill_html {

struct Contig {
    uint32_t id{0}, component{0};
    uint8_t hap{0};
    uint64_t bp{0};
    int64_t start{0};
    bool reverse{false};
    std::string name, sample, first_vertex, last_vertex;
};

struct Connection {
    uint32_t component{0}, left{0}, right{0}, bridge{0};
    double phase_left{0.0}, phase_right{0.0};
    uint64_t phase_left_target_bp{0}, phase_left_bridge_bp{0};
    uint64_t phase_right_target_bp{0}, phase_right_bridge_bp{0};
    double boundary_left{0.0}, boundary_right{0.0};
    double identity_left{-1.0}, coverage_left{-1.0};
    double identity_right{-1.0}, coverage_right{-1.0};
    double sample_score{0.0}, phase_score{0.5};
    double alignment_score{1.0}, confidence{0.0};
    bool homolog_span{false};
    uint32_t sample_support{0}, spanning_samples{0};
    uint64_t left_cut_bp{0}, right_cut_bp{0};
    uint64_t bridge_left_bp{0}, bridge_right_bp{0};
    std::string status;
};

void save(
    const std::string& file,
    const std::vector<Contig>& contigs,
    const std::vector<Connection>& connections
);

}  // namespace gapfill_html
