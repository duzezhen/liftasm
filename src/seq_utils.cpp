#include "../include/seq_utils.hpp"
#include <array>

namespace seqUtils {

static std::array<unsigned char,256> make_rc_table() {
    std::array<unsigned char,256> t{};
    for (int i = 0; i < 256; ++i) t[i] = static_cast<unsigned char>(i);
    t['A']='T'; t['T']='A';
    t['C']='G'; t['G']='C';
    t['a']='t'; t['t']='a';
    t['c']='g'; t['g']='c';
    t['N']='N'; t['n']='n';
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