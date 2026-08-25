#pragma once

#include "CommandOptions.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace rsite {

struct Enzyme {
    std::string name;         // e.g. EcoRI / Enzyme1
    std::string spec;         // e.g. G^AATTC / G^ANTC
    std::string pattern;      // spec without '^'
    uint32_t    cut_offset=0; // 0-based cut on forward motif
};

struct Site {
    uint32_t pos = 0;         // 0-based cut coordinate
    uint16_t enzyme_id = 0;   // index in enzymes_
    char     strand = '+';    // '+' direct match; '-' reverse-complement match
};

struct Nearest {
    bool     found = false;
    uint32_t pos   = 0;
    int64_t  dist  = 0;       // site_pos - query_pos
};

struct QueryPair {
    Nearest beg_site;
    Nearest end_site;
};

struct ChromResult {
    std::string chr;
    uint32_t chr_len = 0;
    std::vector<Site> sites;
    std::vector<uint32_t> cuts;
};

namespace detail {

struct CompiledEnzyme {
    std::string name;
    std::string spec;
    std::string pattern;
    uint32_t cut_offset = 0;
    uint32_t len = 0;
    uint32_t rev_cut_offset = 0;

    std::vector<uint8_t> fwd_mask;
    std::vector<uint8_t> rev_mask;
};

} // namespace detail

Enzyme parse_enzyme_spec(const std::string& spec, size_t order_1based);
std::vector<Enzyme> parse_enzyme_specs(const std::vector<std::string>& specs);

class Index {
public:
    Index() = default;

    void clear();

    void build(const ResSitOpts& params);

    size_t n_chrom() const;
    size_t n_enzyme() const;

    const std::vector<std::string>& chrom_order() const;

    const std::vector<uint32_t>* cuts(std::string_view chr) const;
    const std::vector<Site>* sites(std::string_view chr) const;

    Nearest nearest(std::string_view chr, uint32_t pos) const;
    QueryPair nearest(std::string_view chr, uint32_t beg, uint32_t end) const;
    QueryPair flanking(std::string_view chr, uint32_t beg, uint32_t end) const;

    void save_bed(const std::string& out_bed) const;

private:
    std::vector<detail::CompiledEnzyme> enzymes_;
    std::vector<std::string> chrom_order_;
    std::unordered_map<std::string, uint32_t> chrom_len_;
    std::unordered_map<std::string, std::vector<Site>> chrom_sites_;
    std::unordered_map<std::string, std::vector<uint32_t>> chrom_cuts_;

    ChromResult scan_one_chrom_(const std::string& chr, const std::string& seq, bool scan_revcomp) const;
};

} // namespace rsite
