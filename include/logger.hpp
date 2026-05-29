// #pragma once
// #include <iostream>
// #include "get_time.hpp"

// inline constexpr std::size_t LOG_FUNC_COL_WIDTH = 18;
// inline constexpr char FILLER = '.';

// inline bool DEBUG_ENABLED;

// inline std::string pad_func_name(const char* func, std::size_t width = LOG_FUNC_COL_WIDTH, char filler = FILLER) {
//     std::string s(func ? func : "");
//     if (s.size() < width) {
//         s.append(width - s.size(), filler);
//     } else if (s.size() > width) {
//         s.resize(width);
//     }
//     return s;
// }

// // Format: [L::<func>::MM-DD HH:MM:SS] <message>
// #define log_stream() (std::cerr << "[I::" << getTime() << "::" << pad_func_name(__func__) << "] ")
// #define error_stream() (std::cerr << "[E::" << getTime() << "::" << pad_func_name(__func__) << "] ")
// #define warning_stream() (std::cerr << "[W::" << getTime() << "::" << pad_func_name(__func__) << "] ")


// struct NullBuf : std::streambuf {
//     int overflow(int c) override { return c; }
// };

// inline std::ostream& null_stream() {
//     static NullBuf buf;
//     static std::ostream ns(&buf);
//     return ns;
// }

// // // Format: [L::<FILE>:<LINE> <func>::] <message>
// // #define debug_stream() (DEBUG_ENABLED ? (std::cerr << "[D::" << __FILE__ << ":" << __LINE__ << " " << pad_func_name(__func__) << "] ") : null_stream())
// // Format: [L::<FILE>:<LINE> <func>::] <message>
// #define debug_stream() (DEBUG_ENABLED ? (std::cerr << "[D::" << __FILE__ << "::" << pad_func_name(__func__) << "] ") : null_stream())

// inline void set_debug(bool on) {
//     DEBUG_ENABLED = on;
// }



// 2026-04-07
#pragma once
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <sstream>
#include <iomanip>

inline auto START_WALL = std::chrono::steady_clock::now();
inline std::clock_t START_CPU = std::clock();

inline bool DEBUG_ENABLED = false;

inline std::string fmt_func(const char* func, std::size_t width = 10) {
    std::string s(func ? func : "");
    if (s.size() < width) s.append(width - s.size(), '.');
    else if (s.size() > width) s.resize(width);
    return s;
}

inline std::string fmt_func(double number, std::size_t width = 5) {
    std::string s(std::to_string(number));
    auto dot = s.find('.');
    size_t int_len = (dot == std::string::npos) ? s.size() : dot;
    width = std::max(width, int_len);
    if (s.size() < width) s.append(width - s.size(), '.');
    else if (s.size() > width) s.resize(width);
    return s;
}

inline double get_wall_sec() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now() - START_WALL).count();
}

inline double get_cpu_usage() {
    double wall = get_wall_sec();
    double cpu  = double(std::clock() - START_CPU) / CLOCKS_PER_SEC;
    if (wall <= 0.0) return 0.0;
    return cpu / wall;
}

// Format: [M::func::wall*usage]
#define log_stream() \
    (std::cerr << "[M::" \
               << fmt_func(__func__) << "::" \
               << fmt_func(get_wall_sec(), 6) << "*" \
               << fmt_func(get_cpu_usage(), 4) << "] ")

#define warning_stream() \
    (std::cerr << "[W::" \
               << fmt_func(__func__) << "::" \
               << fmt_func(get_wall_sec(), 6) << "*" \
               << fmt_func(get_cpu_usage(), 4) << "] ")

#define error_stream() \
    (std::cerr << "[E::" \
               << fmt_func(__func__) << "::" \
               << fmt_func(get_wall_sec(), 6) << "*" \
               << fmt_func(get_cpu_usage(), 4) << "] ")

struct NullBuf : std::streambuf {
    int overflow(int c) override { return c; }
};

inline std::ostream& null_stream() {
    static NullBuf buf;
    static std::ostream ns(&buf);
    return ns;
}

#define debug_stream() \
    (DEBUG_ENABLED ? \
        (std::cerr << "[D::" \
                   << fmt_func(__func__) << "::" \
                   << fmt_func(get_wall_sec(), 6) << "*" \
                   << fmt_func(get_cpu_usage(), 4) << "] ") \
        : null_stream())

inline void set_debug(bool on) {
    DEBUG_ENABLED = on;
}