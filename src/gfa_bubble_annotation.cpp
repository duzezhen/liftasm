#include "../include/gfa_bubble_annotation.hpp"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <unordered_map>

namespace GfaBubble {

BubbleAnnotation::BubbleAnnotation(
    const std::vector<std::string>& sample_names,
    double max_mosaic_freq
) : max_mosaic_freq_(max_mosaic_freq)
{
    std::unordered_map<std::string, uint32_t> tissues;
    tissues.reserve(sample_names.size());
    sample_tissues_.reserve(sample_names.size());

    for (const std::string& sample : sample_names) {
        const auto [it, inserted] = tissues.emplace(
            tissue_name_(sample), static_cast<uint32_t>(tissues.size())
        );
        sample_tissues_.push_back(it->second);
        if (inserted) ++tissue_count_;
    }
}

void BubbleAnnotation::annotate(VcfRecord& record) const
{
    record.mosaic = false;
    if (record.alt.empty() || record.alt == "." || tissue_count_ == 0) {
        record.mt = ".";
        record.tf = ".";
        return;
    }

    const uint32_t allele_count = 2 + static_cast<uint32_t>(
        std::count(record.alt.begin(), record.alt.end(), ',')
    );
    std::vector<uint8_t> called(tissue_count_, 0);
    std::vector<uint8_t> support(static_cast<size_t>(tissue_count_) * allele_count, 0);

    const size_t sample_count = std::min(sample_tissues_.size(), record.sample_gt_tokens.size());
    for (size_t sample = 0; sample < sample_count; ++sample) {
        const uint32_t tissue = sample_tissues_[sample];
        const std::string& genotype = record.sample_gt_tokens[sample];

        for (size_t i = 0; i < genotype.size();) {
            if (!std::isdigit(static_cast<unsigned char>(genotype[i]))) {
                ++i;
                continue;
            }

            uint32_t allele = 0;
            do {
                allele = allele * 10 + static_cast<uint32_t>(genotype[i] - '0');
                ++i;
            } while (i < genotype.size() && std::isdigit(static_cast<unsigned char>(genotype[i])));

            called[tissue] = 1;
            if (allele < allele_count) {
                support[static_cast<size_t>(tissue) * allele_count + allele] = 1;
            }
        }
    }

    const uint32_t called_tissues = std::accumulate(called.begin(), called.end(), uint32_t(0));
    if (called_tissues == 0) {
        record.mt = ".";
        record.tf = ".";
        return;
    }

    std::vector<std::string> types;
    std::vector<std::string> frequencies;
    types.reserve(allele_count);
    frequencies.reserve(allele_count);

    for (uint32_t allele = 0; allele < allele_count; ++allele) {
        uint32_t carriers = 0;
        for (uint32_t tissue = 0; tissue < tissue_count_; ++tissue) {
            carriers += support[static_cast<size_t>(tissue) * allele_count + allele];
        }

        const double frequency = static_cast<double>(carriers) / called_tissues;
        frequencies.push_back(format_frequency_(frequency));

        constexpr double epsilon = 1e-12;
        if (carriers == 0) {
            types.emplace_back("unsupported");
        } else if (frequency <= max_mosaic_freq_ + epsilon) {
            types.emplace_back("mosaic_candidate");
            record.mosaic = true;
        } else if (frequency + epsilon >= 1.0 - max_mosaic_freq_) {
            types.emplace_back("germline_likely");
        } else {
            types.emplace_back("tissue_shared");
        }
    }

    auto join = [](const std::vector<std::string>& values) {
        std::string out;
        for (const std::string& value : values) {
            if (!out.empty()) out.push_back(',');
            out += value;
        }
        return out;
    };
    record.mt = join(types);
    record.tf = join(frequencies);
}

std::string BubbleAnnotation::tissue_name_(std::string_view sample)
{
    if (sample.ends_with(".hap1") || sample.ends_with(".hap2")) {
        sample.remove_suffix(5);
    }
    return std::string(sample);
}

std::string BubbleAnnotation::format_frequency_(double value)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

} // namespace GfaBubble
