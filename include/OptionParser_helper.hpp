#pragma once

#include "logger.hpp"

#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

static inline void parse_repeat_mask_(const char* s, CollapseOpts& opt) {
    unsigned a = 0, b = 0, c = 0;
    char extra = 0;

    if (!s || std::sscanf(s, "%u,%u,%u%c", &a, &b, &c, &extra) != 3) {
        error_stream() << "--repeat_mask expects MIN_LEN,MAX_PERIOD,MAX_MISMATCH, e.g. --repeat_mask 24,12,2\n";
        std::exit(1);
    }

    opt.repeat_mask_min_len = std::max(0u, a);
    opt.repeat_mask_max_period = std::max(0u, b);
    opt.repeat_mask_max_mismatch = std::max(0u, c);
}

static inline void parse_abnormal_seg_(const char* arg, CollapseOpts& opt) {
    std::string s(arg ? arg : "");
    const size_t p = s.find(',');

    if (p == std::string::npos) {
        error_stream() << "invalid --abnormal_seg, expected LEN,COUNT, got: " << s << "\n";
        std::exit(1);
    }

    opt.max_abnormal_cut_len = std::max<uint32_t>(0, static_cast<uint32_t>(std::stoul(s.substr(0, p))));
    opt.min_abnormal_cut_count = std::max<uint32_t>(0, static_cast<uint32_t>(std::stoul(s.substr(p + 1))));
}

static inline bool is_flag_(const char* s) { return s && s[0] == '-'; }

static inline std::string trim_copy_(const std::string& s) {
    const size_t first = s.find_first_not_of(" \t\n\r\f\v");
    if (first == std::string::npos) return "";
    const size_t last = s.find_last_not_of(" \t\n\r\f\v");
    return s.substr(first, last - first + 1);
}

static inline std::string lower_copy_(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

static inline bool is_size_unit_token_(const char* arg) {
    const std::string unit = lower_copy_(trim_copy_(arg ? arg : ""));
    return unit == "b" || unit == "k" || unit == "kb" || unit == "m" || unit == "mb" || unit == "g" || unit == "gb";
}

static inline std::string join_size_arg_(const char* arg, int argc, char** argv, int& optind) {
    std::string s(arg ? arg : "");
    if (optind < argc && !is_flag_(argv[optind]) && is_size_unit_token_(argv[optind])) {
        s += " ";
        s += argv[optind];
        ++optind;
    }
    return s;
}

static inline uint64_t parse_size_arg_(const char* arg, const char* opt_name) {
    std::string s(arg ? arg : "");
    auto die = [&]() {
        error_stream() << "invalid " << opt_name << ", expected non-negative number with optional unit " << "(k/kb/m/mb/g/gb), got: " << s << "\n";
        std::exit(1);
    };

    s = trim_copy_(s);
    if (s.empty()) die();

    errno = 0;
    char* end = nullptr;
    const long double value = std::strtold(s.c_str(), &end);
    if (end == s.c_str() || errno == ERANGE || !std::isfinite(value) || value < 0) die();

    while (*end && std::isspace(static_cast<unsigned char>(*end))) ++end;
    const std::string unit = lower_copy_(trim_copy_(std::string(end)));

    long double multiplier = 1.0L;
    if (unit.empty() || unit == "b") {
        multiplier = 1.0L;
    } else if (unit == "k" || unit == "kb") {
        multiplier = 1000.0L;
    } else if (unit == "m" || unit == "mb") {
        multiplier = 1000000.0L;
    } else if (unit == "g" || unit == "gb") {
        multiplier = 1000000000.0L;
    } else {
        die();
    }

    const long double scaled = value * multiplier;
    const long double rounded = std::floor(scaled + 0.5L);
    if (!std::isfinite(scaled) ||
        rounded > static_cast<long double>(std::numeric_limits<uint64_t>::max()) ||
        std::fabs(scaled - rounded) > 1e-6L) {
        die();
    }

    return static_cast<uint64_t>(rounded);
}

static inline uint64_t parse_size_arg_u64_(const char* arg, int argc, char** argv, int& optind, const char* opt_name) {
    const std::string s = join_size_arg_(arg, argc, argv, optind);
    return parse_size_arg_(s.c_str(), opt_name);
}

static inline uint32_t parse_size_arg_u32_(const char* arg, int argc, char** argv, int& optind, const char* opt_name) {
    const uint64_t value = parse_size_arg_u64_(arg, argc, argv, optind, opt_name);
    if (value > std::numeric_limits<uint32_t>::max()) {
        error_stream() << opt_name << " is too large for uint32_t: " << value << "\n";
        std::exit(1);
    }
    return static_cast<uint32_t>(value);
}

static inline uint16_t parse_size_arg_u16_(const char* arg, int argc, char** argv, int& optind, const char* opt_name) {
    const uint64_t value = parse_size_arg_u64_(arg, argc, argv, optind, opt_name);
    if (value > std::numeric_limits<uint16_t>::max()) {
        error_stream() << opt_name << " is too large for uint16_t: " << value << "\n";
        std::exit(1);
    }
    return static_cast<uint16_t>(value);
}

static inline std::vector<uint32_t> parse_size_list_u32_(const std::string& s, const char* opt_name) {
    std::vector<uint32_t> out;
    std::stringstream ss(s);
    std::string tok;

    while (std::getline(ss, tok, ',')) {
        tok = trim_copy_(tok);
        if (tok.empty()) continue;

        const uint64_t value = parse_size_arg_(tok.c_str(), opt_name);
        if (value > std::numeric_limits<uint32_t>::max()) {
            error_stream() << opt_name << " is too large for uint32_t: " << value << "\n";
            std::exit(1);
        }
        out.push_back(static_cast<uint32_t>(value));
    }

    return out;
}

static inline std::string format_size_arg_(uint64_t value) {
    if (value < 1000) return std::to_string(value);

    uint64_t divisor = 1000;
    const char* unit = "kb";
    if (value >= 1000000000ULL) {
        divisor = 1000000000ULL;
        unit = "Gb";
    } else if (value >= 1000000ULL) {
        divisor = 1000000ULL;
        unit = "Mb";
    }

    std::ostringstream oss;
    if (value % divisor == 0) {
        oss << (value / divisor);
    } else {
        oss << std::fixed << std::setprecision(2)
            << (static_cast<long double>(value) / static_cast<long double>(divisor));
        std::string s = oss.str();
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
        s += unit;
        return s;
    }
    oss << unit;
    return oss.str();
}

static inline std::string join_sizes_(const std::vector<uint32_t>& xs) {
    std::string out;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (i) out += ",";
        out += format_size_arg_(xs[i]);
    }
    return out;
}

static inline std::vector<int> parse_int_list_(const std::string& s) {
    std::vector<int> out;
    std::stringstream ss(s);
    std::string tok;

    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        out.push_back(std::stoi(tok));
    }

    return out;
}

static inline std::string join_ints_(const std::vector<int>& xs) {
    std::string out;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (i) out += ",";
        out += std::to_string(xs[i]);
    }
    return out;
}

static inline std::vector<uint32_t> parse_uint32_list_(const std::string& s, const char* opt_name) {
    std::vector<uint32_t> out;
    std::stringstream ss(s);
    std::string tok;

    while (std::getline(ss, tok, ',')) {
        tok = trim_copy_(tok);
        if (tok.empty()) continue;

        size_t used = 0;
        unsigned long value = 0;
        try {
            value = std::stoul(tok, &used);
        } catch (...) {
            error_stream() << "invalid " << opt_name << ", expected comma-separated non-negative integers, got: " << s << "\n";
            std::exit(1);
        }
        if (used != tok.size() || value > std::numeric_limits<uint32_t>::max()) {
            error_stream() << "invalid " << opt_name << ", expected comma-separated non-negative integers, got: " << s << "\n";
            std::exit(1);
        }
        out.push_back(static_cast<uint32_t>(value));
    }

    return out;
}

static inline std::string join_uints_(const std::vector<uint32_t>& xs) {
    std::string out;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (i) out += ",";
        out += std::to_string(xs[i]);
    }
    return out;
}

static inline std::vector<double> parse_double_list_(const std::string& s) {
    std::vector<double> out;
    std::stringstream ss(s);
    std::string tok;

    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        out.push_back(std::stod(tok));
    }

    return out;
}

static inline std::string join_doubles_(const std::vector<double>& xs) {
    std::ostringstream oss;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (i) oss << ",";
        oss << xs[i];
    }
    return oss.str();
}

static inline std::string format_double_(double x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

class HelpPrinter {
public:
    explicit HelpPrinter(std::ostream& out = std::cerr, int kOptWidth = 19, int kTypeWidth = 13) : out_(out), kOptWidth_(kOptWidth), kTypeWidth_(kTypeWidth) {}

    void section(const std::string& title) const {
        out_ << title << ":\n";
    }

    void line(const std::string& opt, const std::string& type, const std::string& desc) const {
        out_ << "  " << std::left << std::setw(kOptWidth_) << opt << std::setw(kTypeWidth_) << type << desc << "\n";
    }

    void note(const std::string& desc) const {
        line("", "", desc);
    }

    void blank() const {
        out_ << "\n";
    }

private:
    int kOptWidth_;
    int kTypeWidth_;
    std::ostream& out_;
};
