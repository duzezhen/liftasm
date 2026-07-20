#include "../include/gfa_augment.hpp"
#include "../include/kvcf.hpp"
#include "../include/gfa_name.hpp"
#include "../include/logger.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <unordered_set>

namespace {

std::string upper(std::string value) {
    for (char& c : value) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return value;
}

bool literal_allele(const std::string& allele) {
    if (allele.empty()) return false;
    for (unsigned char c : allele) {
        c = static_cast<unsigned char>(std::toupper(c));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') return false;
    }
    return true;
}

void append_op(std::vector<CIGAR::COp>& ops, uint32_t len, char op) {
    if (len == 0) return;
    if (!ops.empty() && ops.back().op == op) ops.back().len += len;
    else ops.emplace_back(len, op);
}

std::vector<CIGAR::COp> allele_ops(uint32_t ref_len, uint32_t alt_len) {
    const uint32_t paired = std::min(ref_len, alt_len);
    std::vector<CIGAR::COp> ops;
    ops.reserve(3);
    append_op(ops, paired, 'X');
    append_op(ops, ref_len - paired, 'D');
    append_op(ops, alt_len - paired, 'I');
    return ops;
}

} // namespace

GfaAugmenter::Stats GfaAugmenter::read_vcfs_(const std::vector<std::string>& files) {
    variants_.assign(nodes_.size(), {});
    Stats stats;

    for (const std::string& file : files) {
        kvcf::Reader reader(file);
        kvcf::Record rec;
        while (reader.next(rec)) {
            const auto name_it = name_to_id_map_.find(rec.chrom);
            if (name_it == name_to_id_map_.end()) {
                error_stream() << file + ":" + std::to_string(rec.line_no) + ": CHROM '" + rec.chrom + "' is not a GFA segment" << std::endl;
                std::exit(1);
            }

            const uint32_t sid = static_cast<uint32_t>(name_it->second);
            const GfaNode& node = nodes_[sid];
            if (node.deleted || node.sequence.empty() || node.sequence == "*") {
                error_stream() << file + ":" + std::to_string(rec.line_no) + ": segment '" + rec.chrom + "' has no sequence" << std::endl;
                std::exit(1);
            }
            if (!literal_allele(rec.ref)) {
                error_stream() << file + ":" + std::to_string(rec.line_no) + ": unsupported REF allele" << std::endl;
                std::exit(1);
            }
            for (const std::string& alt : rec.alts) {
                if (!literal_allele(alt)) {
                    error_stream() << file + ":" + std::to_string(rec.line_no) + ": unsupported symbolic/breakend ALT: " + alt << std::endl;
                    std::exit(1);
                }
            }

            const uint64_t beg64 = uint64_t(rec.pos) - 1;
            const uint64_t end64 = beg64 + rec.ref.size();
            if (end64 > node.length || end64 > node.sequence.size()) {
                error_stream() << file + ":" + std::to_string(rec.line_no) + ": REF interval exceeds segment '" + rec.chrom + "'" << std::endl;
                std::exit(1);
            }

            const uint32_t beg = static_cast<uint32_t>(beg64);
            const uint32_t end = static_cast<uint32_t>(end64);
            const std::string ref = upper(rec.ref);
            if (upper(node.sequence.substr(beg, ref.size())) != ref) {
                error_stream() << file + ":" + std::to_string(rec.line_no) + ": REF does not match segment '" + rec.chrom + "'" << std::endl;
                std::exit(1);
            }

            variants_[sid].push_back({sid, beg, end, rec.id, ref, rec.alts, file, rec.line_no});
            ++stats.records;
        }
    }
    return stats;
}

void GfaAugmenter::validate_variant_order_() {
    for (uint32_t sid = 0; sid < variants_.size(); ++sid) {
        auto& variants = variants_[sid];
        std::sort(variants.begin(), variants.end(), [](const VariantInput& a, const VariantInput& b) {
            return a.beg < b.beg || (a.beg == b.beg && a.end < b.end);
        });
        for (size_t i = 1; i < variants.size(); ++i) {
            if (variants[i].beg < variants[i - 1].end) {
                error_stream() << "overlapping VCF records on segment '" + nodes_[sid].name + "': " + variants[i - 1].source + ":" + std::to_string(variants[i - 1].line_no) + " and " + variants[i].source + ":" + std::to_string(variants[i].line_no) << std::endl;
                std::exit(1);
            }
        }
    }
}

void GfaAugmenter::create_alt_nodes_(const std::string& tag) {
    alt_nodes_.clear();
    uint64_t serial = 1;

    for (uint32_t sid = 0; sid < variants_.size(); ++sid) {
        for (const VariantInput& variant : variants_[sid]) {
            for (const std::string& raw_alt : variant.alts) {
                const std::string alt = upper(raw_alt);
                if (alt == variant.ref) continue;

                std::string temporary_name;
                do { temporary_name = "__fuguasm_vcf_" + std::to_string(serial++); }
                while (name_to_id_map_.find(temporary_name) != name_to_id_map_.end());

                const std::vector<uint32_t> sample_ids = nodes_[sid].sample_ids;
                const uint32_t alt_sid = add_segment(temporary_name, alt, true, sample_ids);
                GfaAuxParser::append_str_tag(nodes_[alt_sid].aux.aux_data, "VT", tag);
                if (!variant.id.empty() && variant.id != ".") {
                    GfaAuxParser::append_str_tag(nodes_[alt_sid].aux.aux_data, "VI", variant.id);
                }
                alt_nodes_.push_back({sid, alt_sid, variant.beg, variant.end, variant.id});
            }
        }
    }
}

void GfaAugmenter::inject_variant_alignments_() {
    for (const AltNode& alt : alt_nodes_) {
        cuts_[alt.ref_sid].v.push_back(alt.ref_beg);
        cuts_[alt.ref_sid].v.push_back(alt.ref_end);

        auto ops = allele_ops(alt.ref_end - alt.ref_beg, nodes_[alt.alt_sid].length);
        bubble_aligns_.emplace_back(
            nodes_[alt.ref_sid].name, nodes_[alt.alt_sid].name,
            alt.ref_beg, alt.ref_end, 0, nodes_[alt.alt_sid].length,
            Vertex::make_vertex(alt.ref_sid, false), Vertex::make_vertex(alt.alt_sid, false),
            60, std::move(ops)
        );
        bubble_aligns_.back().force_cuts = true;
    }
}

void GfaAugmenter::connect_alt_nodes_(const SegReplace::Expander& expander) {
    gfaName namer;

    auto materialized_vertex = [&](SegReplace::Seg interval) -> uint32_t {
        const uint32_t sid = static_cast<uint32_t>(SegReplace::Interval::seg_id(interval));
        if (sid >= nodes_.size()) {
            error_stream() << "invalid canonical VCF interval segment" << std::endl;
            std::exit(1);
        }
        const uint32_t beg = SegReplace::Interval::beg(interval);
        const uint32_t end = SegReplace::Interval::end(interval);
        const bool rev = SegReplace::Interval::is_reverse(interval);

        uint32_t output_sid = sid;
        if (beg != 0 || end != nodes_[sid].length) {
            const std::string name = namer.format_interval_name(nodes_[sid].name, beg, end, false);
            const auto it = name_to_id_map_.find(name);
            if (it == name_to_id_map_.end()) {
                error_stream() << "canonical VCF interval was not materialized: " + name << std::endl;
                std::exit(1);
            }
            output_sid = static_cast<uint32_t>(it->second);
        }
        return Vertex::make_vertex(output_sid, rev);
    };

    for (const AltNode& alt : alt_nodes_) {
        const SegReplace::Seg ref = SegReplace::Interval::pack(alt.ref_sid, alt.ref_beg, alt.ref_end, false);
        const SegReplace::Expansion canonical = expander.query(ref);
        if (canonical.empty()) {
            error_stream() << "empty canonical expansion for VCF interval" << std::endl;
            std::exit(1);
        }

        const uint32_t first = materialized_vertex(canonical.front());
        const uint32_t last = materialized_vertex(canonical.back());
        const uint32_t alt_vertex = Vertex::make_vertex(alt.alt_sid, false);

        std::vector<uint32_t> incoming;
        std::vector<uint32_t> outgoing;
        for (const GfaArc* arc : getArcsToVertex(first)) incoming.push_back(arc->get_source_vertex_id());
        for (const GfaArc* arc : getArcsFromVertex(last)) outgoing.push_back(arc->get_target_vertex_id());
        std::sort(incoming.begin(), incoming.end());
        incoming.erase(std::unique(incoming.begin(), incoming.end()), incoming.end());
        std::sort(outgoing.begin(), outgoing.end());
        outgoing.erase(std::unique(outgoing.begin(), outgoing.end()), outgoing.end());

        for (uint32_t source : incoming) add_arc(source, alt_vertex, 0, 0, -1, false);
        for (uint32_t target : outgoing) add_arc(alt_vertex, target, 0, 0, -1, false);
    }
    rebuild_after_edits();
}

void GfaAugmenter::rename_alt_nodes_() {
    uint64_t serial = 1;
    for (AltNode& alt : alt_nodes_) {
        GfaNode& node = nodes_[alt.alt_sid];
        if (node.deleted) continue;

        std::string name;
        do { name = "A" + std::to_string(serial++); }
        while (name_to_id_map_.find(name) != name_to_id_map_.end());
        name_to_id_map_.erase(node.name);
        node.name = name;
        name_to_id_map_[name] = alt.alt_sid;
    }
}

void GfaAugmenter::run_augment_pipeline_(const std::string& prefix) {
    bubble_aligns_.clear();
    cuts_.clear();
    rulemap_.clear();

    finalize_();
    prune_overlaps_();
    initialize_cuts_();
    overlaps_align_();
    dedup_aligns_(bubble_aligns_.size());
    inject_variant_alignments_();

    const auto groups = build_align_groups_();
    build_propagate_prune_cuts_(groups);
    if (DEBUG_ENABLED) print_cuts(cuts_);
    build_rulemap_(groups);

    SegReplace::Expander expander = build_SegReplace_();
    expander.save_map(prefix + ".augment.map");
    expand_and_rewire_edges_(expander);
    connect_alt_nodes_(expander);
    remove_unused_nodes_();
}

GfaAugmenter::Stats GfaAugmenter::augment(
    const std::vector<std::string>& files,
    const std::string& prefix,
    const std::string& tag
) {
    log_stream() << "Augmenting GFA graph with VCF ALT alleles ...\n\n";

    if (tag.empty() || tag.find_first_of("\t\r\n") != std::string::npos) {
        error_stream() << "variant tag must be non-empty and contain no tabs or newlines" << std::endl;
        std::exit(1);
    }

    Stats stats = read_vcfs_(files);
    if (stats.records == 0) return stats;
    validate_variant_order_();
    create_alt_nodes_(tag);
    stats.alleles = alt_nodes_.size();
    stats.split_segments = std::count_if(variants_.begin(), variants_.end(), [](const auto& v) { return !v.empty(); });
    if (alt_nodes_.empty()) return stats;

    run_augment_pipeline_(prefix);
    rename_alt_nodes_();

    log_stream() << "Augmented " << stats.alleles << " ALT allele(s) from " << stats.records << " VCF record(s) across " << stats.split_segments << " segment(s).\n\n";

    return stats;
}
