#pragma once
#pragma once
#include <string>
#include <string_view>

namespace seqUtils {
    const unsigned char* rc_table();

    /**
     * @brief Return the reverse-complement of a DNA sequence.
     * @date 2025-08-05
     * 
     * @param seq Input DNA sequence (A,C,G,T,N,...)
     * 
     * @return Reverse-complemented sequence.
     */
    std::string revcomp(std::string_view seq);
    inline std::string revcomp(const std::string& seq) {
        return revcomp(std::string_view(seq));
    }
} // namespace seqUtils