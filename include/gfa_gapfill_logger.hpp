#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class GfaGapfillDebugger {
public:
    struct BoundarySupport {
        double phase{-1.0};
        uint64_t first_bubble_bp{0}, second_bubble_bp{0};
    };

    // Contig preparation and candidate evidence.
    static void path(const std::string& name, const std::string& decision, const std::string& sample = "", uint8_t hap = 0, uint32_t component = UINT32_MAX, size_t nodes = 0, uint64_t bp = 0, const std::string& first = "", bool first_reverse = false, const std::string& last = "", bool last_reverse = false);
    static void component(const std::string& sample, uint32_t component, size_t paths);
    static void overlap(const std::string& a, const std::string& b, uint64_t shared_bp, uint64_t overlap_bp, double similarity, double a_fraction, double b_fraction, const BoundarySupport& left, const BoundarySupport& right);
    static void search(const std::string& sample, size_t paths, uint64_t tested, uint64_t low_support, uint64_t wrong_group, uint64_t anchor_failed, size_t candidates);
    static void candidate(const std::string& left, const std::string& right, const std::string& bridge, uint32_t component, double left_phase, double right_phase, uint64_t left_target_bp, uint64_t left_bridge_bp, uint64_t right_target_bp, uint64_t right_bridge_bp, bool homolog_span, double phase_score);

    // Candidate selection and node-level graph walk.
    static void selection(const std::string& left, const std::string& right, const std::string& bridge, const std::string& decision);
    static void node_walk(const std::string& left, const std::string& right, const std::string& bridge, uint64_t walk_length, uint64_t left_cut, uint64_t right_cut, const std::vector<std::string>& vertices);

    // Haplotype output.
    static void primary_candidate(uint32_t component, const std::string& sample, uint64_t component_bp, uint64_t hap1_ng50, uint64_t hap2_ng50);
    static void chain(const std::string& name, uint8_t hap, uint32_t component, size_t contigs, size_t gaps, uint64_t bp);
};
