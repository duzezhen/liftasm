#pragma once
#include <cstdint>
#include <string>
#include <string_view>
#include <ostream>

#include "kio.hpp"

namespace paf {

struct Record {
    // 12 mandatory fields
    std::string qname;
    uint32_t    qlen   = 0;
    uint32_t    qstart = 0; // 0-based, inclusive
    uint32_t    qend   = 0; // 0-based, exclusive
    char        strand = '+'; // '+' or '-'
    std::string tname;
    uint32_t    tlen   = 0;
    uint32_t    tstart = 0;
    uint32_t    tend   = 0;
    uint32_t    nmatch = 0;
    uint32_t    alen   = 0;
    uint32_t    mapq   = 0;

    // Example: "cg:Z:10=1X5=\tNM:i:1\ttp:A:P"
    std::string opt;
};

static inline uint32_t parse_u32(std::string_view sv) {
    uint32_t v = 0;
    for (char c : sv) {
        if (c < '0' || c > '9') break;
        v = v * 10u + uint32_t(c - '0');
    }
    return v;
}

static inline bool split_12_fields(std::string_view v, std::string_view f[12], size_t& pos_after_12) {
    size_t pos = 0;
    for (int i = 0; i < 12; ++i) {
        size_t tab = v.find('\t', pos);
        if (tab == std::string_view::npos) {
            if (i != 11) return false;
            f[i] = v.substr(pos);
            pos = v.size();
        } else {
            f[i] = v.substr(pos, tab - pos);
            pos = tab + 1;
        }
    }
    pos_after_12 = pos;
    return true;
}

// Find a SAM-style tag in Record::opt (2-char tag).
// Return true if found; type_out gets 'Z','i','f','A',... value_out is the raw value view.
static inline bool find_tag(const Record& r, std::string_view tag2, char& type_out, std::string_view& value_out)
{
    if (tag2.size() != 2 || r.opt.empty()) return false;

    std::string_view s(r.opt);
    size_t i = 0;
    while (i < s.size()) {
        size_t j = s.find('\t', i);
        if (j == std::string_view::npos) j = s.size();
        std::string_view tok = s.substr(i, j - i);

        // tok: "cg:Z:...."
        if (tok.size() >= 5 &&
            tok[0] == tag2[0] && tok[1] == tag2[1] &&
            tok[2] == ':' && tok[4] == ':')
        {
            type_out  = tok[3];
            value_out = tok.substr(5);
            return true;
        }
        i = j + 1;
    }
    return false;
}

// Transform CIGAR for swapping reference/query.
// -  input: 10M5I10M2D3M
// - output: 3M2I10M5D10M
static inline std::string cigar_swap(std::string_view cg) {
    if (cg.empty()) return std::string();

    struct Op { uint32_t len; char op; };
    std::vector<Op> ops;
    ops.reserve(16);

    uint64_t num = 0;
    for (char c : cg) {
        if (c >= '0' && c <= '9') {
            num = num * 10 + uint64_t(c - '0');
            continue;
        }
        if (num == 0) {
            return std::string();
        }
        char op = c;
        if (op == 'I') op = 'D';
        else if (op == 'D') op = 'I';

        ops.push_back(Op{ (uint32_t)num, op });
        num = 0;
    }
    if (num != 0 || ops.empty()) return std::string();

    std::string out;
    out.reserve(cg.size());

    for (int i = (int)ops.size() - 1; i >= 0; --i) {
        out += std::to_string(ops[(size_t)i].len);
        out.push_back(ops[(size_t)i].op);
    }
    return out;
}

class Reader {
public:
    explicit Reader(const std::string& path) : lr_(path) {}

    bool next(Record& r) {
        std::string line;
        while (lr_.getline(line)) {
            if (line.empty() || line[0] == '#') continue;
            return parse_line_(line, r);
        }
        return false;
    }

private:
    kio::LineReader lr_;

    static bool parse_line_(const std::string& s, Record& r) {
        std::string_view v(s);
        std::string_view f[12];
        size_t pos = 0;

        if (!split_12_fields(v, f, pos)) {
            throw std::runtime_error(
                "PAF parse error: " + s + " (mandatory columns < 12)"
            );
        }

        r.qname  = std::string(f[0]);
        r.qlen   = parse_u32(f[1]);
        r.qstart = parse_u32(f[2]);
        r.qend   = parse_u32(f[3]);
        r.strand = f[4].empty() ? '+' : f[4][0];
        r.tname  = std::string(f[5]);
        r.tlen   = parse_u32(f[6]);
        r.tstart = parse_u32(f[7]);
        r.tend   = parse_u32(f[8]);
        r.nmatch = parse_u32(f[9]);
        r.alen   = parse_u32(f[10]);
        r.mapq   = parse_u32(f[11]);

        if (pos < v.size()) {
            // pos points to beginning of optional fields
            r.opt.assign(v.substr(pos));
        } else {
            r.opt.clear();
        }
        return true;
    }
};

static inline void write_paf(std::ostream& os, const Record& r) {
    os << r.qname  << '\t' << r.qlen   << '\t' << r.qstart << '\t' << r.qend << '\t'
       << r.strand << '\t'
       << r.tname  << '\t' << r.tlen   << '\t' << r.tstart << '\t' << r.tend << '\t'
       << r.nmatch << '\t' << r.alen   << '\t' << r.mapq;
    if (!r.opt.empty()) os << '\t' << r.opt;
    os << '\n';
}

} // namespace paf