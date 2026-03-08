#pragma once
#include <string>
#include <string_view>
#include <cstdint>
#include <vector>
#include <utility>

namespace CIGAR {

struct COp {
    uint32_t len;
    char     op;
    COp() = default;
    COp(uint32_t l, char c) : len(l), op(c) {}
};

/**
 * @brief Parse a CIGAR string into a vector of COp.
 */
std::vector<COp> parse(const std::string &cigar);
std::vector<COp> parse(std::string_view cigar);

/**
 * @brief Prepend one operation.
 */
void prepend(std::vector<COp>& ops, uint32_t len, char op);

/**
 * @brief Append one operation.
 */
void append(std::vector<COp>& ops, uint32_t len, char op);

/**
 * @brief Pack COp vector into CIGAR string.
 */
std::string pack(const std::vector<COp>& ops);


/**
 * @brief Convert BAM packed CIGAR array to string
 */
std::string to_string(const uint32_t* cigar, uint32_t n);

/**
 * @brief Reference span consumed by CIGAR string.
 * Count M/D/N/=/X
 */
uint32_t ref_span(std::string_view cigar);
uint32_t ref_span(const uint32_t* cigar, uint32_t n);

/**
 * @brief Query aligned span consumed by CIGAR string.
 * Count M/I/=/X, but not S/H.
 */
uint32_t query_span_aligned(std::string_view cigar);
uint32_t query_span_aligned(const uint32_t* cigar, uint32_t n);

/**
 * @brief Total read/query span represented by CIGAR string.
 * Count M/I/S/H/=/X
 */
uint32_t query_span(std::string_view cigar);
uint32_t query_span(const uint32_t* cigar, uint32_t n);

/**
 * @brief Match length from CIGAR string.
 * Count M/=/X
 */
uint32_t match_len(std::string_view cigar);
uint32_t match_len(const uint32_t* cigar, uint32_t n);

/**
 * @brief Split reference-consuming blocks by N from CIGAR string.
 * [beg1,end1), [beg2,end2), ...
 */
std::vector<std::pair<uint32_t,uint32_t>> ref_blocks(uint32_t ref_beg, std::string_view cigar);
std::vector<std::pair<uint32_t,uint32_t>> ref_blocks(uint32_t ref_beg, const uint32_t* cigar, uint32_t n);

/**
 * @brief Sum total bases in blocks.
 */
uint32_t total_block_bases(const std::vector<std::pair<uint32_t,uint32_t>>& blocks);

/**
 * @brief Compute aligned query interval [qb,qe) in forward-read coordinates.
 */
bool query_interval_fwd(std::string_view cigar, bool is_rev, uint32_t& qb, uint32_t& qe);
bool query_interval_fwd(const uint32_t* cigar, uint32_t n, bool is_rev, uint32_t& qb, uint32_t& qe);

} // namespace CIGAR