#pragma once
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>

namespace program {

inline constexpr const char* name        = "liftasm";
inline constexpr const char* description = "Merge, de-overlap, and collapse hifiasm GFA assembly graphs.";
inline constexpr const char* version     = "0.1.3-r22";
inline constexpr const char* build_date  = "2026/08/16";
inline constexpr const char* author      = "Zezhen Du";
inline constexpr const char* email       = "dzz0539@gmail.com or zezhen.du@yale.edu";

inline constexpr const char* version_note = "Adds homologous path exploration, enabling support for ultra-long graph (UL) or multiple GFA graphs.";

inline std::string cmdline(int argc, char** argv) {
    std::ostringstream oss;
    std::copy(argv, argv + argc, std::ostream_iterator<const char*>(oss, " "));
    std::string s = oss.str();
    if (!s.empty() && s.back() == ' ')
        s.pop_back();
    return s;
}

}