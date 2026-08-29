#pragma once

#include "gfa_deoverlapper.hpp"

#include <cstdint>
#include <string>
#include <vector>

class GfaAugmenter : public GfaDeoverlapper {
public:
    using GfaDeoverlapper::GfaDeoverlapper;

    struct Stats {
        uint64_t records{0};
        uint64_t alleles{0};
        uint64_t split_segments{0};
    };

    Stats augment(const std::vector<std::string>& vcf_files, const std::string& prefix, const std::string& tag);

protected:
    struct VariantInput {
        uint32_t beg{0};
        uint32_t end{0};
        std::string id;
        std::string ref;
        std::vector<std::string> alts;
    };

    struct AltNode {
        uint32_t ref_sid{0};
        uint32_t alt_sid{0};
        uint32_t ref_beg{0};
        uint32_t ref_end{0};
    };

    Stats read_vcfs_(const std::vector<std::string>& files);
    void sort_variants_();
    void create_alt_nodes_(const std::string& tag);
    void inject_variant_alignments_();
    void connect_alt_nodes_(const SegReplace::Expander& expander);
    void rename_alt_nodes_();
    void run_augment_pipeline_(const std::string& prefix);

    std::vector<std::vector<VariantInput>> variants_;
    std::vector<AltNode> alt_nodes_;
};
