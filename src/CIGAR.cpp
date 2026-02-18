#include "../include/CIGAR.hpp"

namespace CIGAR {

std::vector<COp> parse(const std::string &cig) {
    std::vector<COp> cigar_ops;
    size_t L = cig.size();
    cigar_ops.clear();
    cigar_ops.reserve(L/2 + 1);
    uint32_t n = 0;
    for (char c : cig) {
        if (uint8_t(c - '0') <= 9) {
            n = n*10 + (c - '0');
        } else {
            cigar_ops.emplace_back(n, c);
            n = 0;
        }
    }
    return cigar_ops;
}

void prepend(std::vector<COp> &ops, uint32_t len, char op) {
    if (len == 0) return;
    if (!ops.empty() && ops.front().op == op) {
        ops.front().len += len;
    } else {
        ops.insert(ops.begin(), COp{len,op});
    }
}

void append(std::vector<COp> &ops, uint32_t len, char op) {
    if (len == 0) return;
    if (!ops.empty() && ops.back().op == op) {
        ops.back().len += len;
    } else {
        ops.push_back({len,op});
    }
}

std::string pack(const std::vector<COp> &ops) {
    std::string s;
    s.reserve(ops.size()*4);
    for (auto &o : ops) {
        s += std::to_string(o.len);
        s += o.op;
    }
    return s.empty() ? "*" : s;
}

} // namespace CIGAR