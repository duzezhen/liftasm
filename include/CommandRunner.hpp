#pragma once

#include <string>

int run_command(
    const std::string& subcommand,
    int argc,
    char** argv,
    const std::string& command_line
);
