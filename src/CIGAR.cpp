#include "../include/CIGAR.hpp"

#include <cctype>

namespace CIGAR {

constexpr uint32_t OP_M = 0; // BAM_CMATCH
constexpr uint32_t OP_I = 1; // BAM_CINS
constexpr uint32_t OP_D = 2; // BAM_CDEL
constexpr uint32_t OP_N = 3; // BAM_CREF_SKIP
constexpr uint32_t OP_S = 4; // BAM_CSOFT_CLIP
constexpr uint32_t OP_H = 5; // BAM_CHARD_CLIP
constexpr uint32_t OP_P = 6; // BAM_CPAD
constexpr uint32_t OP_EQ = 7; // BAM_CEQUAL
constexpr uint32_t OP_X = 8; // BAM_CDIFF
constexpr uint32_t OP_B = 9; // BAM_CBACK

static inline uint32_t bam_cigar_op_u32(uint32_t c) {
    return c & 0xfu;
}

static inline uint32_t bam_cigar_oplen_u32(uint32_t c) {
    return c >> 4;
}

static inline char bam_cigar_opchr_u32(uint32_t op) {
    static const char table[] = "MIDNSHP=XB";
    return (op <= OP_B) ? table[op] : '?';
}

std::vector<COp> parse(const std::string& cig) {
    return parse(std::string_view(cig));
}

std::vector<COp> parse(std::string_view cig) {
    std::vector<COp> cigar_ops;
    const size_t L = cig.size();
    cigar_ops.reserve(L / 2 + 1);

    uint32_t n = 0;
    for (char c : cig) {
        if (uint8_t(c - '0') <= 9) {
            n = n * 10 + (uint32_t)(c - '0');
        } else {
            cigar_ops.emplace_back(n, c);
            n = 0;
        }
    }
    return cigar_ops;
}

void prepend(std::vector<COp>& ops, uint32_t len, char op) {
    if (len == 0) return;
    if (!ops.empty() && ops.front().op == op) {
        ops.front().len += len;
    } else {
        ops.insert(ops.begin(), COp{len, op});
    }
}

void append(std::vector<COp>& ops, uint32_t len, char op) {
    if (len == 0) return;
    if (!ops.empty() && ops.back().op == op) {
        ops.back().len += len;
    } else {
        ops.push_back(COp{len, op});
    }
}

std::string pack(const std::vector<COp>& ops) {
    std::string s;
    s.reserve(ops.size() * 4);
    for (const auto& o : ops) {
        s += std::to_string(o.len);
        s += o.op;
    }
    return s.empty() ? "*" : s;
}

std::string to_string(const uint32_t* cigar, uint32_t n) {
    if (!cigar || n == 0) return "*";
    std::string s;
    s.reserve(n * 5);
    for (uint32_t i = 0; i < n; ++i) {
        s += std::to_string(bam_cigar_oplen_u32(cigar[i]));
        s += bam_cigar_opchr_u32(bam_cigar_op_u32(cigar[i]));
    }
    return s;
}

uint32_t ref_span(std::string_view cigar) {
    uint32_t span = 0, num = 0;
    for (char c : cigar) {
        if (std::isdigit((unsigned char)c)) {
            num = num * 10 + (uint32_t)(c - '0');
        } else {
            switch (c) {
                case 'M':
                case 'D':
                case 'N':
                case '=':
                case 'X':
                    span += num;
                    break;
                default:
                    break;
            }
            num = 0;
        }
    }
    return span;
}

uint32_t ref_span(const uint32_t* cigar, uint32_t n) {
    if (!cigar) return 0;
    uint32_t span = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t op  = bam_cigar_op_u32(cigar[i]);
        const uint32_t len = bam_cigar_oplen_u32(cigar[i]);
        switch (op) {
            case OP_M:
            case OP_D:
            case OP_N:
            case OP_EQ:
            case OP_X:
                span += len;
                break;
            default:
                break;
        }
    }
    return span;
}

uint32_t query_span_aligned(std::string_view cigar) {
    uint32_t span = 0, num = 0;
    for (char c : cigar) {
        if (std::isdigit((unsigned char)c)) {
            num = num * 10 + (uint32_t)(c - '0');
        } else {
            switch (c) {
                case 'M':
                case 'I':
                case '=':
                case 'X':
                    span += num;
                    break;
                default:
                    break;
            }
            num = 0;
        }
    }
    return span;
}

uint32_t query_span_aligned(const uint32_t* cigar, uint32_t n) {
    if (!cigar) return 0;
    uint32_t span = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t op  = bam_cigar_op_u32(cigar[i]);
        const uint32_t len = bam_cigar_oplen_u32(cigar[i]);
        switch (op) {
            case OP_M:
            case OP_I:
            case OP_EQ:
            case OP_X:
                span += len;
                break;
            default:
                break;
        }
    }
    return span;
}

uint32_t query_span(std::string_view cigar) {
    uint32_t span = 0, num = 0;
    for (char c : cigar) {
        if (std::isdigit((unsigned char)c)) {
            num = num * 10 + (uint32_t)(c - '0');
        } else {
            switch (c) {
                case 'M':
                case 'I':
                case 'S':
                case 'H':
                case '=':
                case 'X':
                    span += num;
                    break;
                default:
                    break;
            }
            num = 0;
        }
    }
    return span;
}

uint32_t query_span(const uint32_t* cigar, uint32_t n) {
    if (!cigar) return 0;
    uint32_t span = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t op  = bam_cigar_op_u32(cigar[i]);
        const uint32_t len = bam_cigar_oplen_u32(cigar[i]);
        switch (op) {
            case OP_M:
            case OP_I:
            case OP_S:
            case OP_H:
            case OP_EQ:
            case OP_X:
                span += len;
                break;
            default:
                break;
        }
    }
    return span;
}

uint32_t match_len(std::string_view cigar) {
    uint32_t m = 0, num = 0;
    for (char c : cigar) {
        if (std::isdigit((unsigned char)c)) {
            num = num * 10 + (uint32_t)(c - '0');
        } else {
            switch (c) {
                case 'M':
                case '=':
                case 'X':
                    m += num;
                    break;
                default:
                    break;
            }
            num = 0;
        }
    }
    return m;
}

uint32_t match_len(const uint32_t* cigar, uint32_t n) {
    if (!cigar) return 0;
    uint32_t m = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t op  = bam_cigar_op_u32(cigar[i]);
        const uint32_t len = bam_cigar_oplen_u32(cigar[i]);
        switch (op) {
            case OP_M:
            case OP_EQ:
            case OP_X:
                m += len;
                break;
            default:
                break;
        }
    }
    return m;
}

std::vector<std::pair<uint32_t,uint32_t>> ref_blocks(uint32_t ref_beg, std::string_view cigar) {
    std::vector<std::pair<uint32_t,uint32_t>> blocks;
    if (cigar.empty() || cigar == "*") return blocks;

    blocks.reserve(4);

    uint32_t ref_pos = ref_beg;
    uint32_t blk_beg = ref_pos;
    uint32_t blk_len = 0;
    uint32_t num = 0;

    auto flush_block = [&]() {
        if (blk_len > 0) {
            blocks.push_back({blk_beg, blk_beg + blk_len});
            blk_len = 0;
        }
    };

    for (char c : cigar) {
        if (std::isdigit((unsigned char)c)) {
            num = num * 10 + (uint32_t)(c - '0');
            continue;
        }

        switch (c) {
            case 'M':
            case 'D':
            case '=':
            case 'X':
                if (blk_len == 0) blk_beg = ref_pos;
                blk_len += num;
                ref_pos += num;
                break;
            case 'N':
                flush_block();
                ref_pos += num;
                break;
            case 'I':
            case 'S':
            case 'H':
            case 'P':
                break;
            default:
                break;
        }
        num = 0;
    }

    flush_block();
    return blocks;
}

std::vector<std::pair<uint32_t,uint32_t>> ref_blocks(uint32_t ref_beg, const uint32_t* cigar, uint32_t n) {
    std::vector<std::pair<uint32_t,uint32_t>> blocks;
    if (!cigar || n == 0) return blocks;

    blocks.reserve(4);

    uint32_t ref_pos = ref_beg;
    uint32_t blk_beg = ref_pos;
    uint32_t blk_len = 0;

    auto flush_block = [&]() {
        if (blk_len > 0) {
            blocks.push_back({blk_beg, blk_beg + blk_len});
            blk_len = 0;
        }
    };

    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t op  = bam_cigar_op_u32(cigar[i]);
        const uint32_t len = bam_cigar_oplen_u32(cigar[i]);

        switch (op) {
            case OP_M:
            case OP_D:
            case OP_EQ:
            case OP_X:
                if (blk_len == 0) blk_beg = ref_pos;
                blk_len += len;
                ref_pos += len;
                break;
            case OP_N:
                flush_block();
                ref_pos += len;
                break;
            case OP_I:
            case OP_S:
            case OP_H:
            case OP_P:
                break;
            default:
                break;
        }
    }

    flush_block();
    return blocks;
}

uint32_t total_block_bases(const std::vector<std::pair<uint32_t,uint32_t>>& blocks) {
    uint32_t s = 0;
    for (const auto& iv : blocks) {
        if (iv.first < iv.second) s += (iv.second - iv.first);
    }
    return s;
}

bool query_interval_fwd(std::string_view cigar, bool is_rev, uint32_t& qb, uint32_t& qe) {
    qb = 0;
    qe = 0;
    if (cigar.empty() || cigar == "*") return false;

    auto aln_q = [](char op) {
        return op == 'M' || op == 'I' || op == '=' || op == 'X';
    };
    auto read_q = [](char op) {
        return op == 'M' || op == 'I' || op == 'S' || op == 'H' || op == '=' || op == 'X';
    };

    uint32_t qpos = 0;
    uint32_t L0   = 0;
    bool started  = false;
    uint32_t qb_st = 0, qe_st = 0;
    uint32_t num = 0;

    for (char c : cigar) {
        if (std::isdigit((unsigned char)c)) {
            num = num * 10 + (uint32_t)(c - '0');
            continue;
        }

        const char op = c;
        const uint32_t len = num;
        num = 0;

        if (read_q(op)) L0 += len;

        if (aln_q(op)) {
            if (!started) {
                qb_st = qpos;
                started = true;
            }
            qpos += len;
            qe_st = qpos;
        } else {
            if (read_q(op)) qpos += len;
        }
    }

    if (!started || qb_st >= qe_st) return false;

    if (is_rev) {
        qb = L0 - qe_st;
        qe = L0 - qb_st;
    } else {
        qb = qb_st;
        qe = qe_st;
    }
    return qb < qe;
}

bool query_interval_fwd(const uint32_t* cigar, uint32_t n, bool is_rev, uint32_t& qb, uint32_t& qe) {
    qb = 0;
    qe = 0;
    if (!cigar || n == 0) return false;

    auto aln_q = [](uint32_t op) {
        return op == OP_M || op == OP_I || op == OP_EQ || op == OP_X;
    };
    auto read_q = [](uint32_t op) {
        return op == OP_M || op == OP_I || op == OP_S || op == OP_H || op == OP_EQ || op == OP_X;
    };

    uint32_t qpos = 0;
    uint32_t L0   = 0;
    bool started  = false;
    uint32_t qb_st = 0, qe_st = 0;

    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t op  = bam_cigar_op_u32(cigar[i]);
        const uint32_t len = bam_cigar_oplen_u32(cigar[i]);

        if (read_q(op)) L0 += len;

        if (aln_q(op)) {
            if (!started) {
                qb_st = qpos;
                started = true;
            }
            qpos += len;
            qe_st = qpos;
        } else {
            if (read_q(op)) qpos += len;
        }
    }

    if (!started || qb_st >= qe_st) return false;

    if (is_rev) {
        qb = L0 - qe_st;
        qe = L0 - qb_st;
    } else {
        qb = qb_st;
        qe = qe_st;
    }
    return qb < qe;
}

} // namespace CIGAR