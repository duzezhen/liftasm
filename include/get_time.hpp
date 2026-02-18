#pragma once
#include <string>
#include <ctime>

// local time: "MM-DD HH:MM:SS"
inline std::string getTime() {
    std::time_t t = std::time(nullptr);
    struct tm tm_buf;
    localtime_r(&t, &tm_buf);
    char buf[20];
    std::strftime(buf, sizeof(buf), "%m-%d %H:%M:%S", &tm_buf);
    return buf;
}