#pragma once
#include <string>
#include <cstdint>
#include <vector>

namespace CIGAR {

struct COp {
    uint32_t len;
    char     op;
    COp() = default;
    COp(uint32_t l, char c) : len(l), op(c) {}
};

/**
 * @brief Parse a CIGAR string into a vector of COp.
 * @date 2025-08-08
 * 
 * @param cig CIGAR string (e.g., "10M1I5D3M")
 * @param out Output vector of COp, each COp contains length and operation type.
 * 
 * @return Parsed vector of COp.
 */
std::vector<COp> parse(const std::string &cigar);

/**
 * @brief Prepend a COp to the beginning of the ops vector.
 * @date 2025-08-08
 * 
 * If the first operation is the same as the new one, it will increase its length.
 * Otherwise, it will insert a new COp at the front.
 * 
 * @param ops Vector of COp to modify.
 * @param len Length of the operation to prepend.
 * @param op Operation type (e.g., 'M', 'I', 'D').
 */
void prepend(std::vector<COp> &ops, uint32_t len, char op);

/**
 * @brief Append a COp to the end of the ops vector.
 * @date 2025-08-08
 * 
 * If the last operation is the same as the new one, it will increase its length.
 * Otherwise, it will add a new COp at the end.
 * 
 * @param ops Vector of COp to modify.
 * @param len Length of the operation to append.
 * @param op Operation type (e.g., 'M', 'I', 'D').
 */
void append(std::vector<COp> &ops, uint32_t len, char op);

/**
 * @brief Pack a vector of COp into a CIGAR string.
 * @date 2025-08-08
 * 
 * The output string will be in the format of "10M1I5D3M".
 * If the ops vector is empty, it returns "*".
 * 
 * @param ops Vector of COp to pack.
 * @return Packed CIGAR string.
 */
std::string pack(const std::vector<COp> &ops);

} // namespace CIGAR