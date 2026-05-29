#include "../include/seq_utils.hpp"
#include <array>

namespace seqUtils {

static std::array<unsigned char,256> make_rc_table() {
    std::array<unsigned char,256> t{};

    for (int i = 0; i < 256; ++i)
        t[i] = 'N';

    t['A'] = 'T'; t['T'] = 'A';
    t['C'] = 'G'; t['G'] = 'C';
    t['N'] = 'N';

    t['a'] = 'T'; t['t'] = 'A';
    t['c'] = 'G'; t['g'] = 'C';
    t['n'] = 'N';

    return t;
}

const unsigned char* rc_table() {
    static const auto table = make_rc_table();
    return table.data();
}

std::string revcomp(std::string_view seq) {
    auto table = rc_table();
    std::string out(seq.size(),'N');
    for (size_t i = 0, n = seq.size(); i < n; ++i)
        out[n-1-i] = static_cast<char>(table[static_cast<unsigned char>(seq[i])]);
    return out;
}

} // namespace seqUtils