#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace gapfill_html {

struct Contig {
    uint32_t id{0}, component{0};
    uint8_t hap{0};
    uint64_t bp{0}, backbone_bp{0};
    int64_t start{0};
    bool reverse{false};
    std::string name, sample, first_vertex, last_vertex;
};

struct Connection {
    // Spanning contig that proves the two target boundaries are connected.
    uint32_t component{0}, left{0}, right{0}, bridge{0};
    double phase_left{0.0}, phase_right{0.0};
    uint64_t phase_left_target_bp{0}, phase_left_bridge_bp{0};
    uint64_t phase_right_target_bp{0}, phase_right_bridge_bp{0};
    double phase_score{0.5};
    bool homolog_span{false};
    uint64_t left_cut_bp{0}, right_cut_bp{0}; // target contig cuts
    uint64_t bridge_left_bp{0}, bridge_right_bp{0};

    // Node-level path used for the actual gap sequence.
    uint64_t walk_length{0}, walk_left_cut_bp{0}, walk_right_cut_bp{0};
    std::string status;
};

void save(
    const std::string& file,
    const std::vector<Contig>& contigs,
    const std::vector<Connection>& connections
);

}  // namespace gapfill_html
