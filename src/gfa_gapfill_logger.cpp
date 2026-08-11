#include "gfa_gapfill_logger.hpp"

#include "logger.hpp"

#include <iomanip>
#include <sstream>

namespace {

std::string metric(double value) {
    if (value < 0.0) return "NA";
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

std::string boundary_table_header() {
    std::ostringstream out;
    out << "  - " << std::left << std::setw(9) << "Boundary"
        << std::right << std::setw(8) << "Phase"
        << std::setw(18) << "Bubble bp"
        << std::setw(12) << "Minimizer"
        << std::setw(12) << "Alignment";
    return out.str();
}

std::string boundary_table_row(
    const char* label,
    const GfaGapfillDebugger::BoundarySupport& support,
    double alignment
) {
    std::ostringstream bubble;
    bubble << support.first_bubble_bp << '/' << support.second_bubble_bp;
    std::ostringstream out;
    out << "  - " << std::left << std::setw(9) << label
        << std::right << std::setw(8) << metric(support.phase)
        << std::setw(18) << bubble.str()
        << std::setw(12) << metric(support.minimizer)
        << std::setw(12) << metric(alignment);
    return out.str();
}

}  // namespace

void GfaGapfillDebugger::path(
    const std::string& name,
    const std::string& decision,
    const std::string& sample,
    uint8_t hap,
    uint32_t component,
    size_t nodes,
    uint64_t bp,
    const std::string& first,
    bool first_reverse,
    const std::string& last,
    bool last_reverse
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Path: " << name << '\n';
    debug_stream() << "  - Decision: " << decision << '\n';
    if (decision != "keep") return;
    debug_stream() << "  - Sample: " << sample << '\n';
    debug_stream() << "  - Component: " << component << " hap=" << static_cast<uint32_t>(hap) << '\n';
    debug_stream() << "  - Size: nodes=" << nodes << " bp=" << bp << '\n';
    debug_stream() << "  - Endpoints: " << first << (first_reverse ? '-' : '+')
                   << " -> " << last << (last_reverse ? '-' : '+') << '\n';
}

void GfaGapfillDebugger::component(const std::string& sample, uint32_t component, size_t paths) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Component: " << component << '\n';
    debug_stream() << "  - Sample: " << sample << '\n';
    debug_stream() << "  - Contig paths: " << paths << '\n';
}

void GfaGapfillDebugger::overlap(
    const std::string& a,
    const std::string& b,
    uint64_t shared_bp,
    uint64_t overlap_bp,
    double similarity,
    double a_fraction,
    double b_fraction,
    const BoundarySupport& left,
    const BoundarySupport& right
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Overlap support:\n";
    debug_stream() << "  - Paths: " << a << " <-> " << b << '\n';
    debug_stream() << "  - Whole-node bp: shared=" << shared_bp << " overlap=" << overlap_bp << " similarity=" << std::fixed << std::setprecision(3) << similarity << '\n';
    debug_stream() << "  - Contig overlap fractions: " << std::fixed << std::setprecision(3) << a_fraction << ", " << b_fraction << '\n';
    std::ostringstream left_line;
    left_line << std::fixed << std::setprecision(3) << "  - Left overlap boundary: minimizer=";
    if (left.minimizer < 0.0) left_line << "NA";
    else left_line << left.minimizer;
    left_line << " phase=";
    if (left.phase < 0.0) left_line << "NA";
    else left_line << left.phase;
    left_line << " bubble_bp=" << left.first_bubble_bp << '/' << left.second_bubble_bp;
    debug_stream() << left_line.str() << '\n';

    std::ostringstream right_line;
    right_line << std::fixed << std::setprecision(3) << "  - Right overlap boundary: minimizer=";
    if (right.minimizer < 0.0) right_line << "NA";
    else right_line << right.minimizer;
    right_line << " phase=";
    if (right.phase < 0.0) right_line << "NA";
    else right_line << right.phase;
    right_line << " bubble_bp=" << right.first_bubble_bp << '/' << right.second_bubble_bp;
    debug_stream() << right_line.str() << '\n';
}

void GfaGapfillDebugger::search(
    const std::string& sample,
    size_t paths,
    uint64_t tested,
    uint64_t low_support,
    uint64_t wrong_group,
    uint64_t anchor_failed,
    size_t candidates
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Candidate search: " << sample << '\n';
    debug_stream() << "  - Contig paths: " << paths << '\n';
    debug_stream() << "  - Pairs tested: " << tested << '\n';
    debug_stream() << "  - Rejected: low_support=" << low_support << " wrong_group=" << wrong_group << " anchor_failed=" << anchor_failed << '\n';
    debug_stream() << "  - Candidates: " << candidates << '\n';
}

void GfaGapfillDebugger::candidate(
    const std::string& left,
    const std::string& right,
    const std::string& bridge,
    uint32_t component,
    double left_phase,
    double right_phase,
    uint64_t left_target_bp,
    uint64_t left_bridge_bp,
    uint64_t right_target_bp,
    uint64_t right_bridge_bp,
    double left_boundary,
    double right_boundary,
    bool homolog_span,
    uint32_t sample_support,
    uint32_t spanning_samples,
    double sample_score,
    double phase_score,
    double alignment_score,
    double confidence
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Candidate:\n";
    debug_stream() << "  - Component: " << component << '\n';
    debug_stream() << "  - Left: " << left << '\n';
    debug_stream() << "  - Right: " << right << '\n';
    debug_stream() << "  - Bridge: " << bridge << '\n';
    const BoundarySupport left_support{
        left_boundary, left_phase, left_target_bp, left_bridge_bp
    };
    const BoundarySupport right_support{
        right_boundary, right_phase, right_target_bp, right_bridge_bp
    };
    debug_stream() << boundary_table_header() << '\n';
    debug_stream() << boundary_table_row("Left", left_support, -1.0) << '\n';
    debug_stream() << boundary_table_row("Right", right_support, -1.0) << '\n';
    debug_stream() << "  - Other hap spans gap: " << (homolog_span ? "yes" : "no") << '\n';
    debug_stream() << "  - Sample support: " << sample_support << '/' << spanning_samples << '\n';
    debug_stream() << "  - Confidence: sample=" << std::fixed << std::setprecision(3)
                   << sample_score << " phase=" << phase_score
                   << " alignment=" << alignment_score
                   << " combined=" << confidence << '\n';
}

void GfaGapfillDebugger::refined_boundary(
    const std::string& left, uint64_t left_cut,
    const std::string& right, uint64_t right_cut, uint64_t right_length,
    const std::string& bridge, uint64_t bridge_begin, uint64_t bridge_end,
    const BoundarySupport& left_support,
    const BoundarySupport& right_support,
    double left_alignment,
    double right_alignment
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Final gap boundary:\n";
    debug_stream() << "  - Left: " << left << ":0-" << left_cut << '\n';
    debug_stream() << "  - Right: " << right << ':' << right_cut << '-' << right_length << '\n';
    debug_stream() << "  - Bridge: " << bridge << ':' << bridge_begin << '-' << bridge_end << '\n';
    debug_stream() << boundary_table_header() << '\n';
    debug_stream() << boundary_table_row("Left", left_support, left_alignment) << '\n';
    debug_stream() << boundary_table_row("Right", right_support, right_alignment) << '\n';
}

void GfaGapfillDebugger::boundary_alignment(
    const std::string& left,
    const std::string& right,
    const std::string& bridge,
    const Alignment& left_alignment,
    const Alignment& right_alignment
) {
    if (!DEBUG_ENABLED) return;

    debug_stream() << "Boundary alignment: " << left << " -> " << right << " via " << bridge << '\n';
    const Alignment* alignments[2] = {&left_alignment, &right_alignment};
    const char* labels[2] = {"Left", "Right"};
    const std::string* names[2] = {&left, &right};
    for (uint32_t side = 0; side < 2; ++side) {
        const Alignment& alignment = *alignments[side];
        debug_stream() << "  - " << labels[side] << ": " << *names[side] << ':'
                       << alignment.target_begin << '-' << alignment.target_end
                       << " vs " << bridge << ':' << alignment.bridge_begin
                       << '-' << alignment.bridge_end
                       << " query=" << alignment.target_begin + alignment.query_begin
                       << '-' << alignment.target_begin + alignment.query_end
                       << " ref=" << alignment.bridge_begin + alignment.reference_begin
                       << '-' << alignment.bridge_begin + alignment.reference_end
                       << (alignment.reverse ? '-' : '+')
                       << " hits=" << alignment.hits
                       << " mapq=" << static_cast<uint32_t>(alignment.mapq)
                       << " identity=" << std::fixed << std::setprecision(3) << alignment.identity
                       << " coverage=" << alignment.coverage
                       << " " << (alignment.accepted ? "keep" : "drop")
                       << std::endl;
    }
}

void GfaGapfillDebugger::selection(
    const std::string& left,
    const std::string& right,
    const std::string& bridge,
    const std::string& decision,
    double confidence
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Selection: " << decision << '\n';
    debug_stream() << "  - Confidence: " << std::fixed << std::setprecision(2) << confidence << '\n';
    debug_stream() << "  - Left: " << left << '\n';
    debug_stream() << "  - Right: " << right << '\n';
    debug_stream() << "  - Bridge: " << bridge << '\n';
}

void GfaGapfillDebugger::assignment(
    const std::string& sample,
    uint32_t component,
    uint8_t hap,
    uint64_t hap1_bp,
    uint64_t hap2_bp,
    bool conflict
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Haplotype assignment:\n";
    debug_stream() << "  - Sample: " << sample << '\n';
    debug_stream() << "  - Component: " << component << '\n';
    debug_stream() << "  - Input bp: hap1=" << hap1_bp << " hap2=" << hap2_bp << '\n';
    debug_stream() << "  - Output haplotype: " << static_cast<uint32_t>(hap) << '\n';
    if (conflict) {
        debug_stream() << "  - Warning: more than two overlapping contigs\n";
    }
}

void GfaGapfillDebugger::primary_candidate(
    uint32_t component,
    const std::string& sample,
    uint64_t component_bp,
    uint64_t hap1_ng50,
    uint64_t hap2_ng50
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Primary candidate: component=" << component
                   << " sample=" << sample
                   << " component_bp=" << component_bp
                   << " hap1_NG50=" << hap1_ng50
                   << " hap2_NG50=" << hap2_ng50
                   << " mean_NG50=" << (hap1_ng50 + hap2_ng50) / 2
                   << '\n';
}

void GfaGapfillDebugger::chain(
    const std::string& name,
    uint8_t hap,
    uint32_t component,
    size_t contigs,
    size_t gaps,
    uint64_t bp
) {
    if (!DEBUG_ENABLED) return;
    debug_stream() << "Output contig: " << name << '\n';
    debug_stream() << "  - Component: " << component << " hap=" << static_cast<uint32_t>(hap) << '\n';
    debug_stream() << "  - Sources: contigs=" << contigs << " gaps=" << gaps << '\n';
    debug_stream() << "  - Length: " << bp << " bp\n";
}
