#include "../include/gfa_bubble_logger.hpp"

#include "../include/logger.hpp"


namespace GfaBubble {

HomologousPathEnumeratorLogger::HomologousPathEnumeratorLogger(const GfaGraph& graph)
    : graph_(graph)
{}

std::string HomologousPathEnumeratorLogger::fmt(Vertex v) const
{
    return graph_.getNodeName(v.segment_id()) + (v.is_reverse() ? "-" : "+");
}

void HomologousPathEnumeratorLogger::header(const std::vector<Vertex>& srcs) const
{
    if (!DEBUG_ENABLED) return;

    debug_stream() << "Detect homologous paths from sources ...\n";
    debug_stream() << "  - srcs:\n";

    for (Vertex s : srcs) {
        debug_stream() << "    - " << fmt(s) << "\n";
    }
}

void HomologousPathEnumeratorLogger::starts(const std::vector<Vertex>& starts) const
{
    if (!DEBUG_ENABLED) return;

    debug_stream() << "  - starts:\n";

    for (Vertex s : starts) {
        debug_stream() << "    - " << fmt(s) << "\n";
    }
}

void HomologousPathEnumeratorLogger::compressed_starts(
    size_t before,
    const std::vector<Vertex>& starts
) const {
    if (!DEBUG_ENABLED) return;

    debug_stream() << "  - compressed starts by bubble groups: " << before << " -> " << starts.size() << "\n";

    for (Vertex s : starts) {
        debug_stream() << "    - " << fmt(s) << "\n";
    }
}

void HomologousPathEnumeratorLogger::similarity(
    double min_first,
    double max_first,
    double min_second,
    double max_second,
    double min_hi,
    double max_hi,
    bool any_pair_good
) const {
    if (!DEBUG_ENABLED) return;

    debug_stream() << "  - similarity:\n";
    debug_stream() << "    - first:  min=" << min_first  << " max=" << max_first  << "\n";
    debug_stream() << "    - second: min=" << min_second << " max=" << max_second << "\n";
    debug_stream() << "    - hi:     min=" << min_hi     << " max=" << max_hi     << "\n";
    debug_stream() << "    - any_pair_good: " << any_pair_good << "\n";
}

void HomologousPathEnumeratorLogger::return_reason(const std::string& reason) const
{
    if (!DEBUG_ENABLED) return;

    debug_stream() << "  - return: " << reason << "\n\n";
}

void HomologousPathEnumeratorLogger::break_reason(const std::string& reason) const
{
    if (!DEBUG_ENABLED) return;

    debug_stream() << "    - break: " << reason << "\n";
}

void HomologousPathEnumeratorLogger::sequence(const std::string& seq) const
{
    if (!DEBUG_ENABLED) return;

    debug_stream() << "  - seq: " << seq << "\n";
}

void HomologousPathEnumeratorLogger::kept(
    const std::vector<Vertex>& starts,
    const std::vector<Vertex>& ends,
    const std::vector<uint64_t>& lens,
    const std::vector<uint8_t>& keep
) const {
    if (!DEBUG_ENABLED) return;

    debug_stream() << "  - kept:\n";

    for (size_t i = 0; i < starts.size(); ++i) {
        if (!keep[i]) continue;

        debug_stream() << "    - " << fmt(starts[i]) << " -> " << fmt(ends[i]) << " len: " << lens[i] << "\n";
    }
}

void HomologousPathEnumeratorLogger::enumerate(
    Vertex start,
    Vertex end,
    const std::vector<std::vector<uint32_t>>& paths
) const {
    if (!DEBUG_ENABLED) return;

    debug_stream() << "  - enumerate:\n";
    debug_stream() << "    - " << fmt(start) << " -> " << fmt(end) << " paths: " << paths.size() << "\n";

    for (const std::vector<uint32_t>& path : paths) {
        std::string s;

        for (uint32_t vid : path) {
            s += fmt(Vertex(vid));
            s += " ";
        }

        debug_stream() << "      - " << s << "\n";
    }
}

} // namespace GfaBubble
