#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

class GfaGapfillDebugger {
public:
    struct Source {
        std::string_view name;
        uint64_t length{0};
        bool reverse{false};
    };

    struct BoundarySupport {
        double minimizer{-1.0}, phase{-1.0};
        uint64_t first_bubble_bp{0}, second_bubble_bp{0};
    };

    struct Alignment {
        uint64_t target_begin{0}, target_end{0};
        uint64_t bridge_begin{0}, bridge_end{0};
        uint64_t query_cut{0}, reference_cut{0};
        uint32_t query_begin{0}, query_end{0};
        uint32_t reference_begin{0}, reference_end{0};
        uint32_t hits{0};
        uint8_t mapq{0};
        double identity{0.0}, coverage{0.0};
        bool accepted{false}, reverse{false};
    };

    static void path(const std::string& name, const std::string& decision, const std::string& sample = "", uint8_t hap = 0, uint32_t component = UINT32_MAX, size_t nodes = 0, uint64_t bp = 0, const std::string& first = "", bool first_reverse = false, const std::string& last = "", bool last_reverse = false);
    static void component(const std::string& sample, uint32_t component, size_t paths);
    static void overlap(const std::string& a, const std::string& b, uint64_t shared_bp, uint64_t overlap_bp, double similarity, double a_fraction, double b_fraction, const BoundarySupport& left, const BoundarySupport& right);
    static void search(const std::string& sample, size_t paths, uint64_t tested, uint64_t low_support, uint64_t wrong_group, uint64_t anchor_failed, size_t candidates);
    static void candidate(const std::string& left, const std::string& right, const std::string& bridge, uint32_t component, double left_phase, double right_phase, uint64_t left_target_bp, uint64_t left_bridge_bp, uint64_t right_target_bp, uint64_t right_bridge_bp, double left_boundary, double right_boundary, bool homolog_span, uint32_t sample_support, uint32_t spanning_samples, double sample_score, double phase_score, double alignment_score, double confidence);
    static void boundary_alignment(const Source& left, const Source& right, const Source& bridge, const Alignment& left_alignment, const Alignment& right_alignment);
    static void refined_boundary(const Source& left, uint64_t left_cut, const Source& right, uint64_t right_cut, const Source& bridge, uint64_t bridge_begin, uint64_t bridge_end, const BoundarySupport& left_support, const BoundarySupport& right_support, double left_alignment, double right_alignment);
    static void selection(const std::string& left, const std::string& right, const std::string& bridge, const std::string& decision, double confidence);
    static void assignment(const std::string& sample, uint32_t component, uint8_t hap, uint64_t hap1_bp, uint64_t hap2_bp, bool conflict);
    static void primary_candidate(uint32_t component, const std::string& sample, uint64_t component_bp, uint64_t hap1_ng50, uint64_t hap2_ng50);
    static void chain(const std::string& name, uint8_t hap, uint32_t component, size_t contigs, size_t gaps, uint64_t bp);
};
