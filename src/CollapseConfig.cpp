#include "../include/CollapseConfig.hpp"
#include "../include/logger.hpp"

#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <utility>

namespace {

[[noreturn]] void config_error(
    const char* section,
    const std::string& filename,
    size_t line_no,
    const std::string& reason,
    const std::string& sample = {}
) {
    error_stream() << "Invalid collapse input " << section << "\n";
    const std::string location = line_no > 0
        ? filename + ':' + std::to_string(line_no)
        : filename;
    error_stream() << "  - File: " << location << "\n";
    if (!sample.empty()) error_stream() << "  - Sample: " << sample << "\n";
    error_stream() << "  - Reason: " << reason << "\n";
    std::exit(1);
}

struct Columns {
    static constexpr size_t absent = std::numeric_limits<size_t>::max();

    size_t sample{absent};
    size_t gfa{absent};
    size_t hap1_gfa{absent};
    size_t hap2_gfa{absent};
    size_t vcf{absent};
    size_t count{0};

    bool has(size_t column) const noexcept { return column != absent; }
    bool ctg_mode() const noexcept { return has(hap1_gfa); }
};

Columns parse_header(
    const std::vector<std::string>& header,
    const std::string& filename,
    size_t line_no
) {
    Columns columns;
    columns.count = header.size();

    for (size_t i = 0; i < header.size(); ++i) {
        size_t* column = nullptr;
        if      (header[i] == "sample")   column = &columns.sample;
        else if (header[i] == "gfa")      column = &columns.gfa;
        else if (header[i] == "hap1_gfa") column = &columns.hap1_gfa;
        else if (header[i] == "hap2_gfa") column = &columns.hap2_gfa;
        else if (header[i] == "vcf")      column = &columns.vcf;
        else {
            config_error(
                "header", filename, line_no,
                "unknown column '" + header[i] +
                    "'; supported columns are sample, gfa, hap1_gfa, hap2_gfa, and vcf"
            );
        }

        if (*column != Columns::absent) {
            config_error(
                "header", filename, line_no,
                "duplicate column '" + header[i] + "'"
            );
        }
        *column = i;
    }

    if (!columns.has(columns.sample)) {
        config_error("header", filename, line_no, "missing required column 'sample'");
    }
    if (columns.has(columns.gfa) && (columns.has(columns.hap1_gfa) || columns.has(columns.hap2_gfa))) {
        config_error(
            "header", filename, line_no,
            "column 'gfa' cannot be combined with 'hap1_gfa' or 'hap2_gfa'"
        );
    }
    if (columns.has(columns.hap1_gfa) != columns.has(columns.hap2_gfa)) {
        config_error(
            "header", filename, line_no,
            "columns 'hap1_gfa' and 'hap2_gfa' must be provided together"
        );
    }
    if (!columns.has(columns.gfa) && !columns.ctg_mode()) {
        config_error(
            "header", filename, line_no,
            "expected column 'gfa' for UTG mode or both 'hap1_gfa' and 'hap2_gfa' for CTG mode"
        );
    }
    if (columns.has(columns.vcf) && !columns.ctg_mode()) {
        config_error(
            "header", filename, line_no,
            "column 'vcf' is only supported in CTG mode"
        );
    }

    return columns;
}

}  // namespace

CollapseConfig parse_collapse_config(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) config_error("file", filename, 0, "cannot open file");

    CollapseConfig config;
    std::unordered_map<std::string, size_t> sample_lines;
    Columns layout;
    bool header_read = false;
    size_t line_no = 0;
    std::string line;

    while (std::getline(input, line)) {
        ++line_no;

        std::istringstream fields(line);
        std::vector<std::string> values;
        for (std::string value; fields >> value;) values.push_back(std::move(value));

        if (values.empty() || values.front().front() == '#') continue;
        if (!header_read) {
            layout = parse_header(values, filename, line_no);
            header_read = true;
            continue;
        }

        const std::string sample = layout.sample < values.size()
            ? values[layout.sample]
            : std::string{};
        if (values.size() != layout.count) {
            config_error(
                "row", filename, line_no,
                "expected " + std::to_string(layout.count) +
                    " fields from the header, found " + std::to_string(values.size()),
                sample
            );
        }

        auto require_value = [&](size_t column, const char* name) -> const std::string& {
            if (values[column] == ".") {
                config_error(
                    "row", filename, line_no,
                    "column '" + std::string(name) + "' cannot be '.'",
                    sample
                );
            }
            return values[column];
        };

        require_value(layout.sample, "sample");
        if (layout.ctg_mode()) {
            require_value(layout.hap1_gfa, "hap1_gfa");
            require_value(layout.hap2_gfa, "hap2_gfa");
        } else {
            require_value(layout.gfa, "gfa");
        }

        const auto [it, inserted] = sample_lines.emplace(sample, line_no);
        if (!inserted) {
            config_error(
                "row", filename, line_no,
                "duplicate sample (first defined at line " + std::to_string(it->second) + ")",
                sample
            );
        }

        config.sampleNames.push_back(std::move(values[layout.sample]));
        if (layout.ctg_mode()) {
            config.hap1Files.push_back(std::move(values[layout.hap1_gfa]));
            config.hap2Files.push_back(std::move(values[layout.hap2_gfa]));
            if (layout.has(layout.vcf)) {
                std::string vcf = std::move(values[layout.vcf]);
                if (vcf == ".") vcf.clear();
                config.vcfFiles.push_back(std::move(vcf));
            }
        } else {
            config.gfaFiles.push_back(std::move(values[layout.gfa]));
        }
    }

    if (!header_read) config_error("header", filename, 0, "missing header");
    if (config.sampleNames.empty()) {
        config_error("file", filename, 0, "no sample rows after the header");
    }

    return config;
}
