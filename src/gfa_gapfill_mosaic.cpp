#include "../include/gfa_gapfill_mosaic.hpp"

#include "../include/gfa_parser_AUX.hpp"
#include "../include/logger.hpp"
#include "../include/save.hpp"

#include <algorithm>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>


namespace {

struct Block {
    uint64_t begin{0};
    uint64_t end{0};
};

struct IndexedPath {
    const GfaBubble::ReferencePath* path{nullptr};
    std::vector<uint64_t> begins;
    uint64_t length{0};
    std::string sequence;
};

struct Occurrence {
    uint32_t path{0};
    uint32_t position{0};
};

struct Placement {
    uint32_t bubble{0};
    uint32_t path{0};
    uint32_t primary_allele{0};
    uint64_t left{0};
    uint64_t right{0};
    bool reverse{false};
};

struct Segment {
    std::string name;
    std::string sequence;
    std::string tags;
    uint64_t length{0};
};

struct Edge {
    uint32_t from{0};
    uint32_t to{0};
    bool from_reverse{false};
    bool to_reverse{false};

    auto key() const {
        return std::tie(from, from_reverse, to, to_reverse);
    }
};

struct PathRecord {
    std::string name;
    std::vector<std::pair<uint32_t, bool>> segments;
    std::string tags;
};

struct OutputGraph {
    std::vector<Segment> segments;
    std::vector<Edge> edges;
    std::vector<PathRecord> paths;
};

class PrimaryIndex {
public:
    PrimaryIndex(
        const GfaGraph& graph,
        const std::vector<GfaBubble::Bubble>& bubbles,
        const std::vector<GfaBubble::ReferencePath>& primary_paths
    ) : graph_(graph), bubbles_(bubbles)
    {
        std::unordered_set<uint32_t> anchors;
        anchors.reserve(bubbles.size() * 2);
        for (const auto& bubble : bubbles) {
            if (bubble.get_source() == UINT32_MAX || bubble.get_sink() == UINT32_MAX) continue;
            anchors.insert(bubble.get_source());
            anchors.insert(bubble.get_sink() ^ 1u);
        }
        occurrences_.reserve(anchors.size());

        paths_.reserve(primary_paths.size());
        for (const auto& path : primary_paths) {
            IndexedPath indexed;
            indexed.path = &path;
            indexed.begins.resize(path.vertices.size());

            uint64_t assembled = 0;
            bool valid = !path.vertices.empty();
            for (uint32_t i = 0; valid && i < path.vertices.size(); ++i) {
                const uint32_t vertex = path.vertices[i];
                const uint32_t sid = Vertex::get_segment_id(vertex);
                if (sid >= graph.getNumNodes() || graph.getNodeDeleted(sid)) {
                    valid = false;
                    break;
                }

                const uint64_t length = graph.getNodeLength(sid);
                const uint64_t overlap = i == 0 ? 0 : std::min<uint64_t>(
                    graph.get_edge_ow(path.vertices[i - 1], vertex), length
                );
                if (overlap > assembled) {
                    valid = false;
                    break;
                }
                indexed.begins[i] = assembled - overlap;
                assembled = indexed.begins[i] + length;
            }
            if (!valid || assembled == 0) continue;

            indexed.length = assembled;
            indexed.sequence = graph.get_path_sequence(path.vertices);
            if (indexed.sequence == "*" || indexed.sequence.size() != indexed.length) {
                indexed.sequence.clear();
            }

            const uint32_t path_id = static_cast<uint32_t>(paths_.size());
            paths_.push_back(std::move(indexed));
            for (uint32_t position = 0; position < path.vertices.size(); ++position) {
                const uint32_t vertex = path.vertices[position];
                if (anchors.contains(vertex)) occurrences_[vertex].push_back({path_id, position});
            }
        }
    }

    const std::vector<IndexedPath>& paths() const { return paths_; }

    std::vector<Placement> place_bubbles() const
    {
        std::vector<Placement> placements;
        placements.reserve(bubbles_.size());

        for (uint32_t bubble_id = 0; bubble_id < bubbles_.size(); ++bubble_id) {
            const auto& alleles = bubbles_[bubble_id].get_paths();
            Placement placement;
            uint32_t match_begin = 0, match_end = 0;
            bool found = false, ambiguous = false;

            for (uint32_t allele_id = 0; allele_id < alleles.size() && !ambiguous; ++allele_id) {
                const auto& allele = alleles[allele_id];
                if (allele.size() < 2) continue;

                for (bool reverse : {false, true}) {
                    const uint32_t first = reverse ? allele.back() ^ 1u : allele.front();
                    const auto occurrence = occurrences_.find(first);
                    if (occurrence == occurrences_.end()) continue;

                    for (const Occurrence& hit : occurrence->second) {
                        const IndexedPath& path = paths_[hit.path];
                        if (allele.size() > path.path->vertices.size() - hit.position) continue;

                        bool equal = true;
                        for (uint32_t i = 0; i < allele.size(); ++i) {
                            const uint32_t expected = reverse ? allele[allele.size() - i - 1] ^ 1u : allele[i];
                            if (path.path->vertices[hit.position + i] != expected) {
                                equal = false;
                                break;
                            }
                        }
                        if (!equal) continue;

                        const uint32_t end = hit.position + static_cast<uint32_t>(allele.size()) - 1;
                        const uint64_t first_pos = path.begins[hit.position]
                            + graph_.getNodeLength(Vertex::get_segment_id(path.path->vertices[hit.position])) - 1;
                        const uint64_t last_pos = path.begins[end];
                        const uint64_t source_pos = reverse ? last_pos : first_pos;
                        const uint64_t sink_pos = reverse ? first_pos : last_pos;

                        if (!found) {
                            placement = {
                                bubble_id, hit.path, allele_id,
                                std::min(source_pos, sink_pos) + 1,
                                std::max(source_pos, sink_pos), reverse
                            };
                            match_begin = hit.position;
                            match_end = end;
                            found = true;
                        } else if (placement.path != hit.path || match_begin != hit.position || match_end != end) {
                            ambiguous = true;
                            break;
                        }
                    }
                    if (ambiguous) break;
                }
            }

            if (found && !ambiguous && placement.left <= placement.right) {
                placements.push_back(placement);
            }
        }
        return placements;
    }

private:
    const GfaGraph& graph_;
    const std::vector<GfaBubble::Bubble>& bubbles_;
    std::vector<IndexedPath> paths_;
    std::unordered_map<uint32_t, std::vector<Occurrence>> occurrences_;
};

static std::vector<std::vector<Block>> build_blocks(
    const std::vector<IndexedPath>& paths,
    const std::vector<GfaBubble::MosaicSite>& sites,
    uint64_t flank
) {
    std::unordered_map<std::string, std::vector<Block>> raw;
    raw.reserve(paths.size());
    for (const auto& site : sites) {
        raw[site.chrom].push_back({site.begin, site.end});
    }

    std::vector<std::vector<Block>> blocks(paths.size());
    for (uint32_t path_id = 0; path_id < paths.size(); ++path_id) {
        const IndexedPath& path = paths[path_id];
        const auto found = raw.find(path.path->name);
        if (found == raw.end()) continue;

        for (const Block& site : found->second) {
            const uint64_t begin = site.begin > flank ? site.begin - flank : 0;
            const uint64_t end = site.end >= path.length || flank >= path.length - site.end
                ? path.length : site.end + flank;
            if (begin < end) blocks[path_id].push_back({begin, end});
        }

        auto& path_blocks = blocks[path_id];
        std::sort(path_blocks.begin(), path_blocks.end(), [](const Block& a, const Block& b) {
            return a.begin != b.begin ? a.begin < b.begin : a.end < b.end;
        });
        size_t kept = 0;
        for (const Block& block : path_blocks) {
            if (kept > 0 && block.begin <= path_blocks[kept - 1].end) {
                path_blocks[kept - 1].end = std::max(path_blocks[kept - 1].end, block.end);
            } else {
                path_blocks[kept++] = block;
            }
        }
        path_blocks.resize(kept);
    }
    return blocks;
}

static bool overlaps_block(const Placement& placement, const std::vector<Block>& blocks)
{
    const uint64_t end = std::max(placement.right, placement.left + 1);
    for (const Block& block : blocks) {
        if (placement.left < block.end && block.begin < end) return true;
        if (block.begin >= end) break;
    }
    return false;
}

static std::string sample_names(const GfaGraph& graph, uint32_t sid)
{
    std::string out;
    for (uint32_t sample : graph.getNodeSampleIds(sid)) {
        if (sample >= graph.getSampleNames().size()) continue;
        if (!out.empty()) out.push_back(',');
        out += graph.getSampleNames()[sample];
    }
    return out;
}

static std::string clean_command(std::string command)
{
    for (char& c : command) {
        if (c == '\t' || c == '\n' || c == '\r') c = ' ';
    }
    return command;
}

} // namespace


GfaGapfillMosaic::GfaGapfillMosaic(
    const GfaGraph& graph,
    const std::vector<GfaBubble::Bubble>& bubbles,
    const std::vector<GfaBubble::ReferencePath>& primary_paths,
    const std::vector<GfaBubble::MosaicSite>& mosaic_sites,
    uint64_t flank
) : graph_(graph), bubbles_(bubbles), primary_paths_(primary_paths),
    mosaic_sites_(mosaic_sites), flank_(flank)
{}

void GfaGapfillMosaic::save(const std::string& prefix, const std::string& command_line) const
{
    log_stream() << "Writing mosaic-focused graph to " << prefix << ".mosaic.gfa ...\n";

    const PrimaryIndex index(graph_, bubbles_, primary_paths_);
    const auto& paths = index.paths();
    const std::vector<Placement> placements = index.place_bubbles();
    const std::vector<std::vector<Block>> blocks = build_blocks(paths, mosaic_sites_, flank_);

    std::vector<std::vector<Placement>> retained(paths.size());
    for (const Placement& placement : placements) {
        if (placement.path < blocks.size() && overlaps_block(placement, blocks[placement.path])) {
            retained[placement.path].push_back(placement);
        }
    }

    OutputGraph output;
    std::unordered_map<uint32_t, uint32_t> original_nodes;
    original_nodes.reserve(bubbles_.size());

    auto add_original_node = [&](uint32_t sid) {
        const auto found = original_nodes.find(sid);
        if (found != original_nodes.end()) return found->second;

        const GfaNode* node = graph_.getNode(sid);
        if (!node || node->deleted) return UINT32_MAX;

        std::string tags;
        const std::string samples = sample_names(graph_, sid);
        if (!samples.empty()) tags += "\tSN:Z:" + samples;
        if (node->is_complex) tags += "\tCX:i:1";
        tags += GfaAuxParser::write_aux_as_text(node->aux.aux_data);

        const uint32_t id = static_cast<uint32_t>(output.segments.size());
        output.segments.push_back({node->name, node->sequence, std::move(tags), node->length});
        original_nodes.emplace(sid, id);
        return id;
    };

    size_t retained_bubbles = 0;
    size_t block_count = 0;
    for (uint32_t path_id = 0; path_id < paths.size(); ++path_id) {
        const IndexedPath& path = paths[path_id];
        block_count += blocks[path_id].size();

        std::vector<uint64_t> cuts{0, path.length};
        cuts.reserve(retained[path_id].size() * 2 + blocks[path_id].size() * 2 + 2);
        for (const Block& block : blocks[path_id]) {
            cuts.push_back(block.begin);
            cuts.push_back(block.end);
        }
        for (const Placement& placement : retained[path_id]) {
            cuts.push_back(placement.left);
            cuts.push_back(placement.right);
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

        std::vector<uint32_t> backbone;
        backbone.reserve(cuts.size() - 1);
        for (size_t i = 1; i < cuts.size(); ++i) {
            if (cuts[i - 1] >= cuts[i]) continue;
            const uint64_t begin = cuts[i - 1];
            const uint64_t end = cuts[i];
            const std::string name = path.path->name + ':' + std::to_string(begin) + '-' + std::to_string(end);
            const std::string sequence = path.sequence.empty() ? "*" : path.sequence.substr(begin, end - begin);
            backbone.push_back(static_cast<uint32_t>(output.segments.size()));
            output.segments.push_back({name, sequence, {}, end - begin});
        }

        for (size_t i = 1; i < backbone.size(); ++i) {
            output.edges.push_back({backbone[i - 1], backbone[i], false, false});
        }

        PathRecord primary_path;
        primary_path.name = path.path->name;
        for (uint32_t segment : backbone) primary_path.segments.push_back({segment, false});
        output.paths.push_back(std::move(primary_path));

        for (const Placement& placement : retained[path_id]) {
            const auto left_cut = std::lower_bound(cuts.begin(), cuts.end(), placement.left);
            const auto right_cut = std::lower_bound(cuts.begin(), cuts.end(), placement.right);
            if (left_cut == cuts.begin() || right_cut == cuts.end() || *left_cut != placement.left || *right_cut != placement.right) continue;

            const size_t left_index = static_cast<size_t>(left_cut - cuts.begin() - 1);
            const size_t right_index = static_cast<size_t>(right_cut - cuts.begin());
            if (left_index >= backbone.size() || right_index >= backbone.size()) continue;

            const auto& alleles = bubbles_[placement.bubble].get_paths();
            ++retained_bubbles;
            for (uint32_t allele_id = 0; allele_id < alleles.size(); ++allele_id) {
                if (allele_id == placement.primary_allele || alleles[allele_id].size() < 2) continue;

                std::vector<uint32_t> allele = alleles[allele_id];
                if (placement.reverse) {
                    std::reverse(allele.begin(), allele.end());
                    for (uint32_t& vertex : allele) vertex ^= 1u;
                }

                PathRecord route;
                route.name = "m" + std::to_string(placement.bubble) + ".p" + std::to_string(allele_id);
                const auto& bubble = bubbles_[placement.bubble];
                const uint32_t source = bubble.get_source();
                const uint32_t sink = bubble.get_sink();
                route.tags = "\tBI:i:" + std::to_string(placement.bubble)
                    + "\tPI:i:" + std::to_string(allele_id)
                    + "\tSR:Z:" + graph_.getNodeName(Vertex::get_segment_id(source))
                    + (Vertex::get_is_reverse(source) ? "-" : "+")
                    + "\tSK:Z:" + graph_.getNodeName(Vertex::get_segment_id(sink))
                    + (Vertex::get_is_reverse(sink) ? "-" : "+");
                route.segments.push_back({backbone[left_index], false});

                for (size_t i = 1; i + 1 < allele.size(); ++i) {
                    const uint32_t sid = Vertex::get_segment_id(allele[i]);
                    const uint32_t segment = add_original_node(sid);
                    if (segment == UINT32_MAX) continue;
                    route.segments.push_back({segment, Vertex::get_is_reverse(allele[i])});
                }
                route.segments.push_back({backbone[right_index], false});

                for (size_t i = 1; i < route.segments.size(); ++i) {
                    output.edges.push_back({
                        route.segments[i - 1].first, route.segments[i].first,
                        route.segments[i - 1].second, route.segments[i].second
                    });
                }
                output.paths.push_back(std::move(route));
            }
        }
    }

    std::sort(output.edges.begin(), output.edges.end(), [](const Edge& a, const Edge& b) {
        return a.key() < b.key();
    });
    output.edges.erase(std::unique(output.edges.begin(), output.edges.end(), [](const Edge& a, const Edge& b) {
        return a.key() == b.key();
    }), output.edges.end());

    SAVE with_seq(prefix + ".mosaic.gfa");
    SAVE no_seq(prefix + ".mosaic.noseq.gfa");
    std::string header = "H\tVN:Z:1.0\nH\tTS:Z:liftasm\n";
    if (!command_line.empty()) header += "H\tCL:Z:" + clean_command(command_line) + '\n';
    header += "H\tCO:Z:  BI:i   source bubble id\n";
    header += "H\tCO:Z:  SR:Z   original bubble source node\n";
    header += "H\tCO:Z:  SK:Z   original bubble sink node\n";
    with_seq.save(header);
    no_seq.save(header);

    for (const Segment& segment : output.segments) {
        const std::string tags = "\tLN:i:" + std::to_string(segment.length) + segment.tags + '\n';
        with_seq.save("S\t" + segment.name + '\t' + (segment.sequence.empty() ? "*" : segment.sequence) + tags);
        no_seq.save("S\t" + segment.name + "\t*" + tags);
    }
    for (const Edge& edge : output.edges) {
        const std::string line = "L\t" + output.segments[edge.from].name + '\t'
            + (edge.from_reverse ? '-' : '+') + '\t' + output.segments[edge.to].name + '\t'
            + (edge.to_reverse ? '-' : '+') + "\t0M\n";
        with_seq.save(line);
        no_seq.save(line);
    }
    for (const PathRecord& path : output.paths) {
        if (path.segments.empty()) continue;
        std::ostringstream line;
        line << "P\t" << path.name << '\t';
        for (size_t i = 0; i < path.segments.size(); ++i) {
            if (i) line << ',';
            line << output.segments[path.segments[i].first].name << (path.segments[i].second ? '-' : '+');
        }
        line << '\t';
        if (path.segments.size() == 1) line << '*';
        else {
            for (size_t i = 1; i < path.segments.size(); ++i) {
                if (i > 1) line << ',';
                line << "0M";
            }
        }
        line << path.tags << '\n';
        with_seq.save(line.str());
        no_seq.save(line.str());
    }

    log_stream() << "  - Mosaic blocks: " << block_count << "\n";
    log_stream() << "  - Bubbles retained: " << retained_bubbles << "\n";
    log_stream() << "  - Segments written: " << output.segments.size() << "\n\n";
}
