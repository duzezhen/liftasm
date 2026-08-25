#pragma once
#include <cstdint>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>
#include "CommandOptions.hpp"
#include "gfa_parser.hpp"
#include "gfa_name.hpp"
#include "seq_utils.hpp"

// Directly inherit GfaGraph to use its protected API
class Gfa2fa : public GfaGraph {
public:
    explicit Gfa2fa(Gfa2FaOpts params) : params_(std::move(params)) {}

    // Export current graph (this) to FASTA; for segments with length < x,
    // extend y bp on both sides. The extension does NOT include the segment itself.
    void dump_to_stream(std::ostream& out) const;

    // Convenience: write to file
    bool dump_to_file(const std::string& fasta_path) const;

private:
    const Gfa2FaOpts params_;

private:
    // Helpers
    static uint32_t real_len_of_(const GfaNode& n);
    static void write_fasta_record_(
        std::ostream& out,
        const std::string& name,
        const std::string& seq,
        uint32_t wrap_width
    );
    bool node_good_(uint32_t seg) const;

    // Left/right extension (return extension only; never includes the core)
    std::string extend_left_seq_(uint32_t v_forward, uint32_t need_ybp) const;
    std::string extend_right_seq_(uint32_t v_forward, uint32_t need_ybp) const;
};
