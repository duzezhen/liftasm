#pragma once
#include <iostream>
#include "get_time.hpp"

inline constexpr std::size_t LOG_FUNC_COL_WIDTH = 18;
inline constexpr char FILLER = '.';

inline bool DEBUG_ENABLED;

inline std::string pad_func_name(const char* func, std::size_t width = LOG_FUNC_COL_WIDTH, char filler = FILLER) {
    std::string s(func ? func : "");
    if (s.size() < width) {
        s.append(width - s.size(), filler);
    } else if (s.size() > width) {
        s.resize(width);
    }
    return s;
}

// Format: [L::<func>::MM-DD HH:MM:SS] <message>
#define log_stream() (std::cerr << "[I::" << getTime() << "::" << pad_func_name(__func__) << "] ")
#define error_stream() (std::cerr << "[E::" << getTime() << "::" << pad_func_name(__func__) << "] ")
#define warning_stream() (std::cerr << "[W::" << getTime() << "::" << pad_func_name(__func__) << "] ")


struct NullBuf : std::streambuf {
    int overflow(int c) override { return c; }
};

inline std::ostream& null_stream() {
    static NullBuf buf;
    static std::ostream ns(&buf);
    return ns;
}

// Format: [L::<FILE>:<LINE> <func>::] <message>
#define debug_stream() (DEBUG_ENABLED ? (std::cerr << "[D::" << __FILE__ << ":" << __LINE__ << " " << pad_func_name(__func__) << "] ") : null_stream())

inline void set_debug(bool on) {
    DEBUG_ENABLED = on;
}