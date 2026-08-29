#include "gfa_gapfill.hpp"

#include "ThreadPool.hpp"
#include "logger.hpp"
#include "progress_tracker.hpp"

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

extern "C" {
#include "minimap.h"
}


uint32_t GfaGapfillMisassembly::source_id_(const GfaGapfill::Fragment& fragment) const {
    const std::string name = fragment.sample + ".hap" + std::to_string(fragment.hap);
    const auto found = graph_.source_ids_.find(name);
    return found == graph_.source_ids_.end() ? UINT32_MAX : found->second;
}

// Return the minimum node support that ends a terminal low-support run.
// Zero disables terminal misassembly detection for the component.
uint32_t GfaGapfillMisassembly::required_haplotype_support_(uint32_t component) const {
    if (graph_.params_.ms_haps > 0) return graph_.params_.ms_haps;
    const auto found = component_haplotypes_.find(component);
    const size_t total = found == component_haplotypes_.end() ? 0 : found->second.size();
    if (total <= 2) return 0;
    return total <= 4 ? 2 : 3;
}

uint32_t GfaGapfillMisassembly::node_support_(
    const GfaGapfill::Fragment& fragment,
    uint32_t segment
) const {
    const auto component = component_haplotypes_.find(fragment.component);
    if (component == component_haplotypes_.end() || segment >= graph_.nodes_.size()) return 0;

    const std::vector<uint32_t>& component_ids = component->second;
    const std::vector<uint32_t>& node_ids = graph_.nodes_[segment].sample_ids;
    size_t i = 0, j = 0;
    uint32_t support = 0;
    while (i < component_ids.size() && j < node_ids.size()) {
        if (component_ids[i] < node_ids[j]) ++i;
        else if (node_ids[j] < component_ids[i]) ++j;
        else {
            ++support;
            ++i;
            ++j;
        }
    }
    return support;
}

void GfaGapfillMisassembly::index_haplotypes_() {
    log_stream() << "Indexing haplotype support by component ...\n";
    component_haplotypes_.clear();

    // P paths are direct haplotype evidence. Add missing P support to the SN index.
    for (const GfaGapfill::Fragment& fragment : graph_.fragments_) {
        const uint32_t source = source_id_(fragment);
        if (source == UINT32_MAX || fragment.component == UINT32_MAX) continue;
        component_haplotypes_[fragment.component].push_back(source);
        for (uint32_t vertex : fragment.vertices) {
            const uint32_t segment = Vertex::get_segment_id(vertex);
            std::vector<uint32_t>& ids = graph_.nodes_[segment].sample_ids;
            const auto pos = std::lower_bound(ids.begin(), ids.end(), source);
            if (pos == ids.end() || *pos != source) ids.insert(pos, source);
        }
    }

    size_t memberships = 0, disabled = 0;
    for (auto& item : component_haplotypes_) {
        std::vector<uint32_t>& ids = item.second;
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        memberships += ids.size();
        disabled += required_haplotype_support_(item.first) == 0;
    }
    log_stream() << "  - Components with contig paths: " << component_haplotypes_.size() << '\n';
    log_stream() << "  - Component haplotype memberships: " << memberships << '\n';
    if (disabled > 0) {
        log_stream() << "  - Components skipped with at most two haplotypes: " << disabled << '\n';
    }
    std::cerr << '\n';
}

void GfaGapfillMisassembly::find_terminals_() {
    log_stream() << "Scanning contig ends for long low-support runs ...\n";
    terminals_.clear();
    fully_low_support_.clear();

    const auto original_range = [](const GfaGapfill::Fragment& fragment, uint32_t begin, uint32_t end) {
        uint64_t begin_bp = fragment.path_bp[begin];
        uint64_t end_bp = fragment.path_bp[end];
        if (fragment.reverse) {
            const uint64_t old_begin = begin_bp;
            begin_bp = fragment.length - end_bp;
            end_bp = fragment.length - old_begin;
        }
        return std::make_pair(begin_bp, end_bp);
    };

    const auto terminal_log = [&](const GfaGapfill::Fragment& fragment, uint32_t begin, uint32_t end, uint64_t length, bool normalized_left, uint32_t required_support) {
        const auto range = original_range(fragment, begin, end);
        const bool original_left = normalized_left != fragment.reverse;
        std::ostringstream line;
        line << "  - Low-support end: " << graph_.paths_[fragment.path_id].name
             << ':' << range.first << '-' << range.second
             << " | side=" << (original_left ? "left" : "right")
             << " | component=" << fragment.component
             << " | length=" << length
             << " | support<" << required_support << '\n';
        return line.str();
    };

    for (uint32_t fragment_id = 0; fragment_id < graph_.fragments_.size(); ++fragment_id) {
        const GfaGapfill::Fragment& fragment = graph_.fragments_[fragment_id];
        const uint32_t required_support = required_haplotype_support_(fragment.component);
        const uint32_t count = static_cast<uint32_t>(fragment.vertices.size());
        if (required_support == 0 || count == 0) continue;

        const auto has_low_support = [&](uint32_t pos) {
            return node_support_(fragment, Vertex::get_segment_id(fragment.vertices[pos])) < required_support;
        };

        // A left candidate is the maximal contiguous low-support prefix. It is
        // retained only when its sequence length reaches --ms_len.
        uint32_t left_anchor = 0;
        while (left_anchor < count && has_low_support(left_anchor)) ++left_anchor;
        if (left_anchor == count) {
            fully_low_support_.push_back(fragment_id);
            log_stream() << "  - Fully low-support contig: "
                         << graph_.paths_[fragment.path_id].name
                         << " | component=" << fragment.component
                         << " | length=" << fragment.length
                         << " | support<" << required_support << " throughout\n";
            continue;
        }
        const uint64_t left_bp = fragment.path_bp[left_anchor];
        if (left_anchor > 0 && left_bp >= graph_.params_.ms_len) {
            terminals_.push_back({
                fragment_id, 0, left_anchor, UINT32_MAX,
                left_bp, 0.0, true, false, false, {}
            });
            log_stream() << terminal_log(
                fragment, 0, left_anchor, left_bp, true, required_support
            );
        }

        // Apply the same definition to the maximal low-support suffix.
        uint32_t right_anchor = count - 1;
        while (right_anchor > 0 && has_low_support(right_anchor)) --right_anchor;
        if (right_anchor == 0 && has_low_support(0)) continue;
        const uint32_t right_begin = right_anchor + 1;
        const uint64_t right_bp = fragment.length - fragment.path_bp[right_begin];
        if (right_begin < count && right_bp >= graph_.params_.ms_len) {
            terminals_.push_back({
                fragment_id, right_begin, count, UINT32_MAX,
                right_bp, 0.0, false, false, false, {}
            });
            log_stream() << terminal_log(
                fragment, right_begin, count, right_bp, false, required_support
            );
        }
    }

    log_stream() << "  - Long low-support ends: " << terminals_.size() << '\n';
    if (!fully_low_support_.empty()) {
        log_stream() << "  - Fully low-support contigs: " << fully_low_support_.size() << '\n';
    }
    std::cerr << '\n';
}

std::vector<uint32_t> GfaGapfillMisassembly::select_references_() const {
    log_stream() << "Selecting trusted cross-component reference contigs ...\n";
    struct Group {
        uint32_t component{UINT32_MAX}, source{UINT32_MAX};
        uint64_t total_bp{0}, longest_bp{0};
        std::vector<uint32_t> fragments;
    };

    const auto key = [](uint32_t component, uint32_t source) {
        return (static_cast<uint64_t>(component) << 32) | source;
    };

    std::unordered_set<uint64_t> possible;
    possible.reserve(terminals_.size() + fully_low_support_.size());
    for (const Terminal& terminal : terminals_) {
        const GfaGapfill::Fragment& fragment = graph_.fragments_[terminal.fragment];
        possible.insert(key(fragment.component, source_id_(fragment)));
    }
    for (uint32_t fragment_id : fully_low_support_) {
        const GfaGapfill::Fragment& fragment = graph_.fragments_[fragment_id];
        possible.insert(key(fragment.component, source_id_(fragment)));
    }

    std::unordered_map<uint64_t, Group> groups;
    groups.reserve(graph_.fragments_.size());
    for (uint32_t fragment_id = 0; fragment_id < graph_.fragments_.size(); ++fragment_id) {
        const GfaGapfill::Fragment& fragment = graph_.fragments_[fragment_id];
        const uint32_t source = source_id_(fragment);
        if (source == UINT32_MAX || fragment.component == UINT32_MAX) continue;
        Group& group = groups[key(fragment.component, source)];
        group.component = fragment.component;
        group.source = source;
        group.total_bp += fragment.length;
        group.longest_bp = std::max(group.longest_bp, fragment.length);
        group.fragments.push_back(fragment_id);
    }

    std::unordered_map<uint32_t, const Group*> selected;
    for (const auto& item : groups) {
        const Group& group = item.second;
        if (possible.count(item.first)) continue;
        const auto found = selected.find(group.component);
        if (found == selected.end() ||
            group.total_bp > found->second->total_bp ||
            (group.total_bp == found->second->total_bp && group.longest_bp > found->second->longest_bp) ||
            (group.total_bp == found->second->total_bp && group.longest_bp == found->second->longest_bp && group.source < found->second->source)) {
            selected[group.component] = &group;
        }
    }

    std::vector<uint32_t> references;
    for (const auto& item : selected) {
        references.insert(
            references.end(), item.second->fragments.begin(), item.second->fragments.end()
        );
    }
    std::sort(references.begin(), references.end());
    log_stream() << "  - Trusted component haplotypes: " << selected.size() << '\n';
    log_stream() << "  - Trusted contigs: " << references.size() << "\n\n";
    return references;
}

void GfaGapfillMisassembly::confirm_terminals_(const std::vector<uint32_t>& references) {
    struct ReferenceChunk {
        uint32_t fragment{0};
        uint64_t begin{0}, end{0};
    };
    constexpr uint64_t CHUNK_BP = 20'000'000;
    constexpr uint64_t CHUNK_OVERLAP = 100'000;

    std::vector<std::string> target_names;
    std::vector<uint32_t> target_terminals;
    target_names.reserve(terminals_.size());
    target_terminals.reserve(terminals_.size());
    for (uint32_t i = 0; i < terminals_.size(); ++i) {
        Terminal& terminal = terminals_[i];
        terminal.sequence = graph_.fragment_sequence_(
            graph_.fragments_[terminal.fragment], terminal.begin, terminal.end
        );
        if (terminal.sequence.empty() || terminal.sequence == "*") continue;
        target_terminals.push_back(i);
        target_names.push_back("terminal" + std::to_string(i));
    }

    if (target_terminals.empty() || references.empty()) {
        warning_stream() << "  ! Cross-component checking skipped: no usable target or reference sequence\n\n";
        return;
    }

    mm_idxopt_t index_options;
    mm_mapopt_t map_options;
    mm_set_opt(nullptr, &index_options, &map_options);
    if (mm_set_opt(graph_.mm2_.preset.c_str(), &index_options, &map_options) < 0) {
        warning_stream() << "  ! Cross-component checking skipped: invalid minimap2 preset\n\n";
        return;
    }
    index_options.k = static_cast<short>(graph_.mm2_.k);
    index_options.w = static_cast<short>(graph_.mm2_.w);

    std::vector<const char*> sequences, names;
    sequences.reserve(target_terminals.size());
    names.reserve(target_names.size());
    for (size_t i = 0; i < target_terminals.size(); ++i) {
        sequences.push_back(terminals_[target_terminals[i]].sequence.c_str());
        names.push_back(target_names[i].c_str());
    }
    const int hpc = (index_options.flag & MM_I_HPC) ? 1 : 0;
    mm_idx_t* index = mm_idx_str(
        index_options.w, index_options.k, hpc, index_options.bucket_bits,
        static_cast<int>(sequences.size()), sequences.data(), names.data()
    );
    if (!index) {
        warning_stream() << "  ! Cross-component checking skipped: cannot build minimap2 index\n\n";
        return;
    }

    map_options.flag |= MM_F_CIGAR | MM_F_EQX;
    map_options.best_n = static_cast<short>(std::min<size_t>(
        std::numeric_limits<short>::max(), std::max<size_t>(graph_.mm2_.best_n, terminals_.size())
    ));
    map_options.zdrop = graph_.mm2_.zdrop;
    mm_mapopt_update(&map_options, index);

    // Keep long, repetitive contigs from becoming single-thread bottlenecks.
    std::vector<ReferenceChunk> chunks;
    for (uint32_t reference_id : references) {
        const uint64_t length = graph_.fragments_[reference_id].length;
        for (uint64_t begin = 0; begin < length; begin += CHUNK_BP) {
            chunks.push_back({
                reference_id, begin, std::min(length, begin + CHUNK_BP + CHUNK_OVERLAP)
            });
        }
    }

    using Intervals = std::vector<std::pair<uint32_t, uint32_t>>;
    std::vector<std::unordered_map<uint32_t, Intervals>> covered(terminals_.size());
    const auto collect = [&](std::vector<Alignment>&& alignments) {
        for (const Alignment& alignment : alignments) {
            covered[alignment.terminal][alignment.component].emplace_back(
                alignment.begin, alignment.end
            );
        }
    };

    log_stream() << "Checking low-support ends against other components ...\n";
    log_stream() << "  - Trusted sequence chunks: " << chunks.size() << '\n';
    ProgressTracker progress(chunks.size());
    ThreadPool pool(std::min<size_t>(graph_.params_.threads, chunks.size()));
    std::vector<std::vector<Alignment>> results(chunks.size());
    pool.parallel_for(chunks.size(), [&](size_t chunk_id) {
        const ReferenceChunk& chunk = chunks[chunk_id];
        std::vector<Alignment>& result = results[chunk_id];
        const GfaGapfill::Fragment& reference = graph_.fragments_[chunk.fragment];
        const bool has_other_component = std::any_of(
            target_terminals.begin(), target_terminals.end(), [&](uint32_t terminal_id) {
                return graph_.fragments_[terminals_[terminal_id].fragment].component != reference.component;
            }
        );
        if (!has_other_component) {
            progress.hit();
            return;
        }
        const std::string query = graph_.fragment_subsequence_(
            reference, chunk.begin, chunk.end
        );
        if (query.empty() || query == "*" || query.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
            progress.hit();
            return;
        }

        mm_tbuf_t* buffer = mm_tbuf_init();
        if (!buffer) {
            progress.hit();
            return;
        }
        int count = 0;
        const std::string query_name = graph_.paths_[reference.path_id].name + ':' +
            std::to_string(chunk.begin) + '-' + std::to_string(chunk.end);
        mm_reg1_t* hits = mm_map(
            index, static_cast<int>(query.size()), query.c_str(),
            &count, buffer, &map_options, query_name.c_str()
        );
        result.reserve(count);
        for (int i = 0; i < count; ++i) {
            const mm_reg1_t& hit = hits[i];
            if (!hit.p || hit.rid < 0 || static_cast<size_t>(hit.rid) >= target_terminals.size() ||
                hit.blen <= 0 || hit.re <= hit.rs || hit.mapq < graph_.mm2_.min_mapq) continue;
            const uint32_t terminal_id = target_terminals[hit.rid];
            const Terminal& terminal = terminals_[terminal_id];
            if (graph_.fragments_[terminal.fragment].component == reference.component) continue;
            const double identity = static_cast<double>(hit.mlen) / hit.blen;
            if (identity < graph_.mm2_.min_match) continue;
            const uint64_t aligned_span = std::max<int64_t>(
                hit.re - hit.rs, hit.qe - hit.qs
            );
            const uint64_t shorter = std::min<uint64_t>(terminal.sequence.size(), query.size());
            if (shorter == 0 || static_cast<double>(aligned_span) / shorter < graph_.mm2_.min_ali_ratio) continue;
            result.push_back({
                terminal_id, reference.component,
                static_cast<uint32_t>(hit.rs), static_cast<uint32_t>(hit.re)
            });
        }
        if (hits) {
            for (int i = 0; i < count; ++i) std::free(hits[i].p);
            std::free(hits);
        }
        mm_tbuf_destroy(buffer);
        progress.hit();
    });

    for (std::vector<Alignment>& result : results) {
        collect(std::move(result));
    }
    pool.stop();
    mm_idx_destroy(index);
    progress.finish();

    size_t confirmed = 0;
    for (uint32_t terminal_id = 0; terminal_id < terminals_.size(); ++terminal_id) {
        Terminal& terminal = terminals_[terminal_id];
        uint64_t best_bp = 0;
        uint32_t best_component = UINT32_MAX;
        for (auto& item : covered[terminal_id]) {
            Intervals& intervals = item.second;
            std::sort(intervals.begin(), intervals.end());
            uint64_t total = 0;
            uint32_t begin = 0, end = 0;
            bool have = false;
            for (const auto& interval : intervals) {
                if (!have) {
                    begin = interval.first;
                    end = interval.second;
                    have = true;
                } else if (interval.first <= end) {
                    end = std::max(end, interval.second);
                } else {
                    total += end - begin;
                    begin = interval.first;
                    end = interval.second;
                }
            }
            if (have) total += end - begin;
            if (total > best_bp || (total == best_bp && item.first < best_component)) {
                best_bp = total;
                best_component = item.first;
            }
        }
        terminal.coverage = terminal.sequence.empty() ? 0.0 :
            std::min(1.0, static_cast<double>(best_bp) / terminal.sequence.size());
        terminal.target_component = best_component;
        terminal.confirmed = terminal.coverage >= graph_.params_.ms_sim;
        confirmed += terminal.confirmed;

        if (!terminal.confirmed) {
            const GfaGapfill::Fragment& source = graph_.fragments_[terminal.fragment];
            uint64_t begin = source.path_bp[terminal.begin];
            uint64_t end = source.path_bp[terminal.end];
            if (source.reverse) {
                const uint64_t old_begin = begin;
                begin = source.length - end;
                end = source.length - old_begin;
            }
            std::ostringstream line;
            line << "  - Preserved low-support end: "
                 << graph_.paths_[source.path_id].name << ':' << begin << '-' << end;
            if (best_component == UINT32_MAX) line << " | target=none";
            else line << " | target=component" << best_component;
            line << " | coverage=" << std::fixed << std::setprecision(4)
                 << terminal.coverage << '\n';
            log_stream() << line.str();
        }
    }
    log_stream() << "  - Ends confirmed by one other component: " << confirmed << "\n\n";
}

std::vector<GfaGapfill::MisassemblyContig> GfaGapfillMisassembly::split_confirmed_() {
    log_stream() << "Splitting confirmed terminal misassemblies ...\n";
    const auto contig_name = [](const std::string& sample, uint8_t hap, size_t number) {
        std::ostringstream name;
        name << sample << ".h" << static_cast<uint32_t>(hap) << "ms" << std::setw(6) << std::setfill('0') << number << 'l';
        return name.str();
    };
    const auto original_range = [](const GfaGapfill::Fragment& fragment, uint32_t begin, uint32_t end) {
        uint64_t begin_bp = fragment.path_bp[begin];
        uint64_t end_bp = fragment.path_bp[end];
        if (fragment.reverse) {
            const uint64_t old_begin = begin_bp;
            begin_bp = fragment.length - end_bp;
            end_bp = fragment.length - old_begin;
        }
        return std::make_pair(begin_bp, end_bp);
    };

    // ------------------------------------------------ Collect confirmed ends ------------------------------------------------
    std::vector<uint32_t> splits;
    splits.reserve(terminals_.size() + fully_low_support_.size());
    for (uint32_t i = 0; i < terminals_.size(); ++i) {
        const Terminal& terminal = terminals_[i];
        if (!terminal.confirmed || terminal.sequence.empty()) continue;
        splits.push_back(i);
    }
    if (splits.empty()) {
        log_stream() << "  - Confirmed terminal misassemblies: 0\n\n";
        return {};
    }

    // ------------------------------------------------ Find contained whole contigs ------------------------------------------------
    // A fully low-support contig sharing a confirmed cut node belongs to the
    // same untrusted region, so the complete contig is split out.
    const uint32_t confirmed_ends = static_cast<uint32_t>(splits.size());
    std::vector<uint8_t> split_nodes(graph_.nodes_.size(), false);
    for (uint32_t split_id : splits) {
        const Terminal& split = terminals_[split_id];
        const GfaGapfill::Fragment& fragment = graph_.fragments_[split.fragment];
        for (uint32_t pos = split.begin; pos < split.end; ++pos) {
            split_nodes[Vertex::get_segment_id(fragment.vertices[pos])] = true;
        }
    }

    size_t contained = 0;
    for (uint32_t fragment_id : fully_low_support_) {
        const GfaGapfill::Fragment& fragment = graph_.fragments_[fragment_id];
        const bool shares_cut_node = std::any_of(
            fragment.vertices.begin(), fragment.vertices.end(), [&](uint32_t vertex) {
                return split_nodes[Vertex::get_segment_id(vertex)];
            }
        );
        if (!shares_cut_node) continue;

        std::string sequence = graph_.fragment_sequence_(
            fragment, 0, static_cast<uint32_t>(fragment.vertices.size())
        );
        if (sequence.empty() || sequence == "*") continue;
        terminals_.push_back({
            fragment_id, 0, static_cast<uint32_t>(fragment.vertices.size()),
            UINT32_MAX, fragment.length, 0.0,
            false, true, true, std::move(sequence)
        });
        splits.push_back(static_cast<uint32_t>(terminals_.size() - 1));
        ++contained;
    }

    log_stream() << "  - Confirmed terminal misassemblies: " << confirmed_ends << '\n';
    log_stream() << "  - Contained whole low-support contigs: " << contained << '\n';
    if (contained < fully_low_support_.size()) {
        log_stream() << "  - Whole low-support contigs preserved: " << fully_low_support_.size() - contained << '\n';
    }

    // ------------------------------------------------ Mark source paths ------------------------------------------------
    std::vector<uint32_t> retained_begin(graph_.fragments_.size(), 0);
    std::vector<uint32_t> retained_end;
    retained_end.reserve(graph_.fragments_.size());
    for (const GfaGapfill::Fragment& fragment : graph_.fragments_) {
        retained_end.push_back(static_cast<uint32_t>(fragment.vertices.size()));
    }
    for (uint32_t split_id : splits) {
        const Terminal& split = terminals_[split_id];
        if (split.begin == 0) retained_begin[split.fragment] = std::max(retained_begin[split.fragment], split.end);
        if (split.end == graph_.fragments_[split.fragment].vertices.size()) {
            retained_end[split.fragment] = std::min(retained_end[split.fragment], split.begin);
        }
        const GfaGapfill::Fragment& fragment = graph_.fragments_[split.fragment];
        for (uint32_t pos = split.begin; pos < split.end; ++pos) {
            split_nodes[Vertex::get_segment_id(fragment.vertices[pos])] = true;
        }
    }

    // ------------------------------------------------ Create independent misassembly contigs ------------------------------------------------
    std::unordered_set<std::string> path_names;
    path_names.reserve(graph_.paths_.size() + splits.size());
    for (const GfaPath& path : graph_.paths_) path_names.insert(path.name);

    std::vector<GfaGapfill::MisassemblyContig> records;
    records.reserve(splits.size());
    size_t number = 0;
    for (uint32_t split_id : splits) {
        const Terminal& split = terminals_[split_id];
        const GfaGapfill::Fragment& source = graph_.fragments_[split.fragment];
        std::string name;
        do {
            name = contig_name(source.sample, source.hap, ++number);
        } while (path_names.count(name) || graph_.name_to_id_map_.count(name));
        path_names.insert(name);

        std::vector<uint32_t> sample_ids;
        const uint32_t source_id = source_id_(source);
        if (source_id != UINT32_MAX) sample_ids.push_back(source_id);
        const uint32_t segment = graph_.add_segment(name, split.sequence, true, sample_ids);
        GfaPath path;
        path.name = name;
        path.segments.push_back({segment, false});
        graph_.paths_.push_back(std::move(path));

        const bool original_left = split.left != source.reverse;
        records.push_back({
            name, source.sample, graph_.paths_[source.path_id].name,
            split.whole ? "whole" : (original_left ? "left" : "right"),
            source.hap, source.component
        });

        const auto range = original_range(source, split.begin, split.end);
        if (split.whole) {
            log_stream() << "  - Whole low-support misassembly: "
                         << graph_.paths_[source.path_id].name
                         << ':' << range.first << '-' << range.second << '\n';
        } else {
            log_stream() << "  - Confirmed misassembly: " << graph_.paths_[source.path_id].name
                         << ':' << range.first << '-' << range.second
                         << " | target=component" << split.target_component
                         << " | coverage=" << std::fixed << std::setprecision(4)
                         << split.coverage << '\n';
        }
    }

    // ------------------------------------------------ Trim old P paths ------------------------------------------------
    // Trim every confirmed P path after all new contigs have been created.
    for (uint32_t i = 0; i < graph_.fragments_.size(); ++i) {
        if (retained_begin[i] == 0 && retained_end[i] == graph_.fragments_[i].vertices.size()) continue;
        const GfaGapfill::Fragment& source = graph_.fragments_[i];
        GfaPath& path = graph_.paths_[source.path_id];
        if (retained_begin[i] >= retained_end[i]) {
            path.segments.clear();
            continue;
        }

        std::vector<uint32_t> retained(
            source.vertices.begin() + retained_begin[i],
            source.vertices.begin() + retained_end[i]
        );
        if (source.reverse) {
            std::reverse(retained.begin(), retained.end());
            for (uint32_t& vertex : retained) vertex ^= 1u;
        }
        path.segments.clear();
        path.segments.reserve(retained.size());
        for (uint32_t vertex : retained) {
            path.segments.push_back({
                Vertex::get_segment_id(vertex), Vertex::get_is_reverse(vertex)
            });
        }
    }
    graph_.paths_.erase(
        std::remove_if(graph_.paths_.begin(), graph_.paths_.end(), [](const GfaPath& path) {
            return path.segments.empty();
        }),
        graph_.paths_.end()
    );

    // ------------------------------------------------ Remove unused source nodes ------------------------------------------------
    // Old nodes are removed only after no sample contig path uses them.
    std::vector<uint8_t> path_nodes(graph_.nodes_.size(), false);
    for (const GfaPath& path : graph_.paths_) {
        std::string sample;
        uint8_t hap = 0;
        if (!GfaGapfill::parse_path_name_(path.name, sample, hap)) continue;
        for (const PathSegment& segment : path.segments) {
            if (segment.node_id < path_nodes.size()) path_nodes[segment.node_id] = true;
        }
    }
    size_t removed = 0;
    for (uint32_t segment = 0; segment < split_nodes.size(); ++segment) {
        if (split_nodes[segment] && !path_nodes[segment] && graph_.delete_segment(segment)) ++removed;
    }
    if (removed > 0) {
        graph_.paths_.erase(
            std::remove_if(graph_.paths_.begin(), graph_.paths_.end(), [&](const GfaPath& path) {
                return std::any_of(path.segments.begin(), path.segments.end(), [&](const PathSegment& segment) {
                    return segment.node_id < graph_.nodes_.size() && graph_.nodes_[segment.node_id].deleted;
                });
            }),
            graph_.paths_.end()
        );
        graph_.rebuild_after_edits();
    }

    log_stream() << "  - Misassembly contigs created: " << records.size() << '\n';
    log_stream() << "  - Misassembly nodes removed: " << removed << "\n\n";
    return records;
}

std::vector<GfaGapfill::MisassemblyContig> GfaGapfillMisassembly::run() {
    index_haplotypes_();
    find_terminals_();
    if (terminals_.empty()) {
        return {};
    }

    const std::vector<uint32_t> references = select_references_();
    confirm_terminals_(references);
    return split_confirmed_();
}
