#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "gfa_parser.hpp"
#include "gfa_bubble_types.hpp"

namespace GfaBubble {

/* ================================================================================================================
 *                                       BUBBLE ANNOTATION START
 * ================================================================================================================ */
class BubbleAnnotation {
public:
    BubbleAnnotation(const std::vector<std::string>& sample_names, double max_mosaic_freq);

    void annotate(VcfRecord& record) const;

private:
    static std::string tissue_name_(std::string_view sample);
    static std::string format_frequency_(double value);

    double max_mosaic_freq_{0.3};
    std::vector<uint32_t> sample_tissues_;
    uint32_t tissue_count_{0};
};
/* ================================================================================================================
 *                                        BUBBLE ANNOTATION END
 * ================================================================================================================ */

} // namespace GfaBubble
