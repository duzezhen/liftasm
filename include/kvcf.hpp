#pragma once

#include "kio.hpp"

#include <charconv>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace kvcf {

struct Record {
    std::string chrom;
    uint32_t pos{0};                 // 1-based VCF coordinate
    std::string id;
    std::string ref;
    std::vector<std::string> alts;
    uint64_t line_no{0};
};

class Reader {
public:
    explicit Reader(const std::string& path) : lines_(path) {}

    bool next(Record& record) {
        std::string line;
        while (lines_.getline(line)) {
            ++line_no_;
            if (line.empty() || line[0] == '#') continue;
            parse_(line, record);
            record.line_no = line_no_;
            return true;
        }
        return false;
    }

private:
    static std::vector<std::string_view> fields_(const std::string& line) {
        std::vector<std::string_view> out;
        size_t beg = 0;
        while (true) {
            const size_t end = line.find('\t', beg);
            out.emplace_back(line.data() + beg, (end == std::string::npos ? line.size() : end) - beg);
            if (end == std::string::npos) break;
            beg = end + 1;
        }
        return out;
    }

    void parse_(const std::string& line, Record& out) const {
        const auto f = fields_(line);
        if (f.size() < 5) fail_("expected at least 5 tab-separated columns");

        uint64_t pos = 0;
        const char* first = f[1].data();
        const char* last = first + f[1].size();
        const auto parsed = std::from_chars(first, last, pos);
        if (parsed.ec != std::errc{} || parsed.ptr != last || pos == 0 || pos > UINT32_MAX) {
            fail_("invalid POS");
        }

        out = {};
        out.chrom.assign(f[0]);
        out.pos = static_cast<uint32_t>(pos);
        out.id.assign(f[2]);
        out.ref.assign(f[3]);
        if (out.chrom.empty() || out.ref.empty() || out.ref == ".") fail_("missing CHROM or REF");

        size_t beg = 0;
        while (beg <= f[4].size()) {
            const size_t end = f[4].find(',', beg);
            const std::string_view alt = f[4].substr(beg, end == std::string_view::npos ? f[4].size() - beg : end - beg);
            if (!alt.empty() && alt != ".") out.alts.emplace_back(alt);
            if (end == std::string_view::npos) break;
            beg = end + 1;
        }
        if (out.alts.empty()) fail_("missing ALT");
    }

    [[noreturn]] void fail_(const std::string& message) const {
        throw std::runtime_error("VCF line " + std::to_string(line_no_) + ": " + message);
    }

    kio::LineReader lines_;
    uint64_t line_no_{0};
};

} // namespace kvcf
