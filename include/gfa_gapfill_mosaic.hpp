#pragma once

#include "gfa_parser.hpp"
#include "gfa_bubble_types.hpp"

#include <cstdint>
#include <string>
#include <vector>


/* =================================================================================================================
 *                                         GAPFILL MOSAIC GRAPH START
 * ================================================================================================================= */
class GfaGapfillMosaic {
public:
    GfaGapfillMosaic(
        const GfaGraph& graph,
        const std::vector<GfaBubble::Bubble>& bubbles,
        const std::vector<GfaBubble::ReferencePath>& primary_paths,
        const std::vector<GfaBubble::MosaicSite>& mosaic_sites,
        uint64_t flank
    );

    void save(const std::string& prefix, const std::string& command_line) const;

private:
    const GfaGraph& graph_;
    const std::vector<GfaBubble::Bubble>& bubbles_;
    const std::vector<GfaBubble::ReferencePath>& primary_paths_;
    const std::vector<GfaBubble::MosaicSite>& mosaic_sites_;
    uint64_t flank_{100'000};
};
/* =================================================================================================================
 *                                          GAPFILL MOSAIC GRAPH END
 * ================================================================================================================= */
