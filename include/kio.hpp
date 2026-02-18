#pragma once
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <stdexcept>

#include <zlib.h>

namespace kio {

class File {
public:
    File() = default;
    ~File() { close(); }

    File(const File&) = delete;
    File& operator=(const File&) = delete;

    bool open(const std::string& path) {
        close();
        is_gz_ = endswith_(path, ".gz");
        if (!is_gz_) is_gz_ = endswith_(path, ".GZ");
        if (is_gz_) {
            gz_ = gzopen(path.c_str(), "rb");
            return gz_ != nullptr;
        } else {
            fp_ = std::fopen(path.c_str(), "rb");
            return fp_ != nullptr;
        }
    }

    void close() {
        if (gz_) { gzclose(gz_); gz_ = nullptr; }
        if (fp_) { std::fclose(fp_); fp_ = nullptr; }
        eof_ = false;
    }

    bool good() const { return (gz_ || fp_); }
    bool eof()  const { return eof_; }

    // Return number of bytes read; 0 => EOF
    int read(void* dst, int n) {
        if (!good() || n <= 0) return 0;
        int ret = 0;
        if (is_gz_) {
            ret = gzread(gz_, dst, n);
        } else {
            ret = (int)std::fread(dst, 1, (size_t)n, fp_);
        }
        if (ret <= 0) eof_ = true;
        return ret;
    }

private:
    static bool endswith_(const std::string& s, const char* suf) {
        const size_t n = std::strlen(suf);
        return s.size() >= n && s.compare(s.size() - n, n, suf) == 0;
    }

    bool  is_gz_ = false;
    bool  eof_   = false;
    gzFile gz_   = nullptr;
    FILE*  fp_   = nullptr;
};

// Fast buffered line reader (supports very long lines)
class LineReader {
public:
    explicit LineReader(const std::string& path, int bufsize = 1 << 16)
        : buf_((size_t)bufsize)
    {
        if (!f_.open(path)) {
            throw std::runtime_error("kio::LineReader: cannot open file: " + path);
        }
    }

    bool good() const { return f_.good(); }

    // Read next line; return false at EOF
    bool getline(std::string& out) {
        out.clear();
        for (;;) {
            // If buffer empty -> refill
            if (beg_ >= end_) {
                if (f_.eof()) return false;
                refill_();
                if (end_ == 0) return false;
            }

            // Search '\n' in [beg_, end_)
            const char* base = buf_.data();
            const char* p = (const char*)std::memchr(base + beg_, '\n', (size_t)(end_ - beg_));
            if (p) {
                const int i = (int)(p - base);
                // append [beg_, i)
                out.append(base + beg_, (size_t)(i - beg_));
                beg_ = i + 1;

                // trim trailing '\r'
                if (!out.empty() && out.back() == '\r') out.pop_back();
                return true;
            } else {
                // no newline in current buffer: append all and refill
                out.append(base + beg_, (size_t)(end_ - beg_));
                beg_ = end_;
            }
        }
    }

private:
    void refill_() {
        beg_ = 0;
        end_ = f_.read(buf_.data(), (int)buf_.size());
        if (end_ < 0) end_ = 0;
    }

    File f_;
    std::vector<char> buf_;
    int beg_ = 0;
    int end_ = 0;
};

} // namespace kio