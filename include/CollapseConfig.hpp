#pragma once

#include <string>
#include <vector>

struct CollapseConfig {
    std::vector<std::string> sampleNames;
    std::vector<std::string> gfaFiles;
    std::vector<std::string> hap1Files;
    std::vector<std::string> hap2Files;
    std::vector<std::string> vcfFiles;
};

CollapseConfig parse_collapse_config(const std::string& filename);
