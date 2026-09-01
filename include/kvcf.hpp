#pragma once

#include "kio.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace kvcf {

struct AlleleCall {
    uint32_t sample{0};
    uint32_t haplotype{0};
    uint32_t ploidy{0};
    uint32_t allele{0};
};

struct Record {
    std::string chrom;
    uint32_t pos{0};                 // 1-based VCF coordinate
    std::string id;
    std::string ref;
    std::vector<std::string> alts;
    std::vector<std::string> mutation_types;
    std::vector<AlleleCall> allele_calls;
    uint64_t line_no{0};
};

class Reader {
public:
    explicit Reader(const std::string& path) : path_(path), lines_(path) {}

    bool next(Record& record) {
        std::string line;
        while (lines_.getline(line)) {
            ++line_no_;
            if (line.empty()) continue;
            if (line[0] == '#') {
                if (line.starts_with("#CHROM")) parse_header_(line);
                continue;
            }
            parse_(line, record);
            record.line_no = line_no_;
            return true;
        }
        return false;
    }

    const std::vector<std::string>& samples() const noexcept { return samples_; }

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

    static size_t format_index_(std::string_view format, std::string_view key) {
        size_t index = 0;
        size_t begin = 0;
        while (true) {
            const size_t end = format.find(':', begin);
            if (format.substr(begin, end == std::string_view::npos ? format.size() - begin : end - begin) == key) {
                return index;
            }
            if (end == std::string_view::npos) return SIZE_MAX;
            begin = end + 1;
            ++index;
        }
    }

    static std::string_view format_value_(std::string_view sample, size_t index) {
        size_t begin = 0;
        while (index-- > 0) {
            const size_t end = sample.find(':', begin);
            if (end == std::string_view::npos) return {};
            begin = end + 1;
        }
        const size_t end = sample.find(':', begin);
        return sample.substr(begin, end == std::string_view::npos ? sample.size() - begin : end - begin);
    }

    void parse_header_(const std::string& line) {
        if (header_seen_) fail_("duplicate #CHROM header");
        const auto fields = fields_(line);
        if (fields.size() < 8 || fields[0] != "#CHROM") fail_("invalid #CHROM header");

        std::unordered_set<std::string> seen;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::string sample(fields[i]);
            if (sample.empty() || sample.find(',') != std::string::npos) fail_("invalid sample name");
            if (!seen.insert(sample).second) fail_("duplicate sample name '" + sample + "'");
            samples_.push_back(std::move(sample));
        }
        header_seen_ = true;
    }

    void parse_mt_(std::string_view info, Record& out) const {
        size_t begin = 0;
        while (begin <= info.size()) {
            const size_t end = info.find(';', begin);
            const std::string_view field = info.substr(
                begin, end == std::string_view::npos ? info.size() - begin : end - begin
            );
            if (field.starts_with("MT=") && out.mutation_types.empty()) {
                std::string_view values = field.substr(3);
                size_t value_begin = 0;
                while (value_begin <= values.size()) {
                    const size_t value_end = values.find(',', value_begin);
                    const std::string_view value = values.substr(
                        value_begin,
                        value_end == std::string_view::npos ? values.size() - value_begin : value_end - value_begin
                    );
                    out.mutation_types.emplace_back(value == "." ? "" : value);
                    if (value_end == std::string_view::npos) break;
                    value_begin = value_end + 1;
                }
            }
            if (end == std::string_view::npos) break;
            begin = end + 1;
        }
    }

    void parse_gt_(std::string_view gt, uint32_t sample, Record& out) const {
        if (gt.empty()) fail_("missing GT value");

        const uint32_t ploidy = static_cast<uint32_t>(
            1 + std::count_if(gt.begin(), gt.end(), [](char c) { return c == '/' || c == '|'; })
        );
        uint32_t haplotype = 0;
        size_t begin = 0;
        while (true) {
            const size_t end = gt.find_first_of("/|", begin);
            const std::string_view allele_text = gt.substr(
                begin, end == std::string_view::npos ? gt.size() - begin : end - begin
            );
            if (allele_text.empty()) fail_("invalid GT value");

            if (allele_text != ".") {
                uint64_t allele = 0;
                const char* first = allele_text.data();
                const char* last = first + allele_text.size();
                const auto parsed = std::from_chars(first, last, allele);
                if (parsed.ec != std::errc{} || parsed.ptr != last || allele > UINT32_MAX) {
                    fail_("invalid GT allele");
                }
                out.allele_calls.push_back({sample, haplotype, ploidy, static_cast<uint32_t>(allele)});
            }

            if (end == std::string_view::npos) break;
            begin = end + 1;
            ++haplotype;
        }
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

        out.chrom.clear();
        out.id.clear();
        out.ref.clear();
        out.alts.clear();
        out.mutation_types.clear();
        out.allele_calls.clear();
        out.chrom.assign(f[0]);
        out.pos = static_cast<uint32_t>(pos);
        out.id.assign(f[2]);
        out.ref.assign(f[3]);
        if (out.chrom.empty() || out.ref.empty() || out.ref == ".") fail_("missing CHROM or REF");

        size_t beg = 0;
        while (beg <= f[4].size()) {
            const size_t end = f[4].find(',', beg);
            const std::string_view alt = f[4].substr(beg, end == std::string_view::npos ? f[4].size() - beg : end - beg);
            if (alt.empty() || alt == ".") fail_("missing ALT");
            out.alts.emplace_back(alt);
            if (end == std::string_view::npos) break;
            beg = end + 1;
        }

        if (f.size() > 7 && f[7] != ".") parse_mt_(f[7], out);
        if (samples_.empty()) return;

        if (f.size() != samples_.size() + 9) fail_("sample column count does not match #CHROM header");
        const size_t gt_index = format_index_(f[8], "GT");
        if (gt_index == SIZE_MAX) return;
        for (uint32_t sample = 0; sample < samples_.size(); ++sample) {
            if (f[sample + 9] == ".") continue;
            parse_gt_(format_value_(f[sample + 9], gt_index), sample, out);
        }
    }

    [[noreturn]] void fail_(const std::string& message) const {
        throw std::runtime_error(path_ + ":" + std::to_string(line_no_) + ": " + message);
    }

    std::string path_;
    kio::LineReader lines_;
    bool header_seen_{false};
    std::vector<std::string> samples_;
    uint64_t line_no_{0};
};

} // namespace kvcf
