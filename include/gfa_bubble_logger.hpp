#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "gfa_parser.hpp"
#include "gfa_parser_types.hpp"


namespace GfaBubble {

class HomologousPathEnumeratorLogger {
public:
    explicit HomologousPathEnumeratorLogger(const GfaGraph& graph);

    std::string fmt(Vertex v) const;

    void header(const std::vector<Vertex>& srcs) const;
    void starts(const std::vector<Vertex>& starts) const;
    void compressed_starts(size_t before, const std::vector<Vertex>& starts) const;

    void similarity(
        double min_first,
        double max_first,
        double min_second,
        double max_second,
        double min_hi,
        double max_hi,
        bool any_pair_good
    ) const;

    void return_reason(const std::string& reason) const;
    void break_reason(const std::string& reason) const;
    void sequence(const std::string& seq) const;

    void kept(
        const std::vector<Vertex>& starts,
        const std::vector<Vertex>& ends,
        const std::vector<uint64_t>& lens,
        const std::vector<uint8_t>& keep
    ) const;

    void enumerate(
        Vertex start,
        Vertex end,
        const std::vector<std::vector<uint32_t>>& paths
    ) const;

    void stop_final_intersected_pair(
        Vertex a_start,
        Vertex a_end,
        Vertex b_start,
        Vertex b_end,
        uint32_t rank
    ) const;

    void stop_intersected_branch(
        Vertex start,
        Vertex end,
        Vertex by,
        uint32_t rank
    ) const;

private:
    const GfaGraph& graph_;
};

} // namespace GfaBubble
