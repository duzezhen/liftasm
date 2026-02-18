#pragma once
#include <string>
#include <vector>
#include <cstdint>
#include "gfa_parser.hpp"

// Simple extractor: read compact marker paths and emit concatenated sequences.
// Path grammar (now supports two forms):
//   A) Prefix markers: '<' = reverse, '>' = forward, e.g. "<s1>s2"
//   B) s1+,s2-,s3,s4+
//
class GfaSeq : public GfaGraph {
public:
    GfaSeq() = default;

    void extract_from_file(const std::string& path_file, std::vector<std::string>& paths, std::vector<std::string>& seqs) const;

    void save_to_file(const std::string& out_file, const std::vector<std::string>& paths, const std::vector<std::string>& seqs);

private:
    static inline bool is_marker_(char c) { return c=='<' || c=='>'; }
    static inline bool is_suffix_sign_(char c) { return c=='+' || c=='-'; }
    static void trim_(std::string& s);

    std::vector<std::string> extract_from_paths_(const std::vector<std::string>& paths) const;
};