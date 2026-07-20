#pragma once

#include <cstdint>
#include <string>

#include "gfa_parser.hpp"

class GfaSplitter : public GfaGraph {
public:
    uint32_t split_by_components(
        const std::string& prefix,
        bool write_paths = false,
        bool write_align = false,
        bool skip_comp = true
    );

private:
    bool save_component_to_disk_(
        const std::string& filename,
        uint32_t component_id,
        bool write_paths,
        bool write_align,
        bool write_seq
    ) const;
};