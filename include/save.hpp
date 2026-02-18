#ifndef SAVE_HPP
#define SAVE_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <memory>
#include <zlib.h>
#include <unordered_map>

#include "logger.hpp"
#include "get_time.hpp"

/**
 * @brief save result
 *
 * @param outputFileName
 *
 * @return
**/
class SAVE
{
private:
    // output stream
    std::string outputFileName_;

    // cached flag to avoid repeated string search
    bool is_gzip_{false};

    // file handles
    std::ofstream fpO;
    std::unique_ptr<gzFile_s, int(*)(gzFile)> gzfpO_{nullptr, gzclose};

    // internal write buffer
    std::string buffer_;
    size_t cache_size_{10 * 1024 * 1024};     // default 10 MB

    // disable copy
    SAVE(const SAVE&) = delete;
    SAVE& operator=(const SAVE&) = delete;

    /* flush internal buffer to file */
    void flush();
public:
    SAVE() = default;
    explicit SAVE(const std::string& outFileName, size_t cacheSize = 10 << 20); // default 10 MB
    ~SAVE();

    int save(const std::string& outTxt);
};

#endif
