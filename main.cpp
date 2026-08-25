#include <iomanip>
#include <iostream>
#include <string>

#include "include/CommandRunner.hpp"
#include "include/OptionParser.hpp"
#include "include/ProgramMetadata.hpp"
#include "include/get_time.hpp"
#include "include/logger.hpp"
#include "include/sys.hpp"

namespace {

bool is_flag(const char* value, const char* short_name, const char* long_name) {
    return value && (std::string(value) == short_name || std::string(value) == long_name);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        help(argv);
        return 1;
    }
    if (is_flag(argv[1], "-h", "--help")) {
        help(argv);
        return 0;
    }
    if (is_flag(argv[1], "-H", "--advanced")) {
        help(argv, true);
        return 0;
    }
    if (is_flag(argv[1], "-v", "--version")) {
        std::cerr << program::version << '\n';
        return 0;
    }

    const std::string command_line = program::cmdline(argc, argv);
    log_stream() << "Version: " << program::version << '\n';
    log_stream() << "Data: " << getTime() << '\n';
    log_stream() << "Command: " << command_line << "\n\n";

    const double start_time = realtime();
    const int status = run_command(argv[1], argc, argv, command_line);
    if (status != 0) return status;

    log_stream()
        << "Real time: " << std::fixed << std::setprecision(3)
        << (realtime() - start_time) << " sec; CPU: " << cputime()
        << " sec; Peak RSS: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB\n";
    return 0;
}
