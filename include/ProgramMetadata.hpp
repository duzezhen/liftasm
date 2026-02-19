#pragma once
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>

namespace program {

inline constexpr const char* name        = "liftasm";
inline constexpr const char* description = "Alignment-based coordinate liftover for assemblies and GFA graphs.";
inline constexpr const char* version     = "0.1.0";
inline constexpr const char* build_date  = "2026/02/18";
inline constexpr const char* author      = "Zezhen Du";
inline constexpr const char* email       = "dzz0539@gmail.com or dzz0539@163.com";

inline std::string cmdline(int argc, char** argv) {
    std::ostringstream oss;
    std::copy(argv, argv + argc, std::ostream_iterator<const char*>(oss, " "));
    std::string s = oss.str();
    if (!s.empty() && s.back() == ' ')
        s.pop_back();
    return s;
}

} // namespace program