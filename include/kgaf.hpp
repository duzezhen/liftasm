#pragma once

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>
#include <utility>
#include <istream>
#include <zlib.h>

/*
* kgaf::ReaderAuto r("aln.gaf.gz");  // auto-detect gzip vs plain
* kgaf::Record rec;
* while (r.next(rec)) {
*     // use rec.qname, rec.path, rec.cg, ...
* }
*/

namespace kgaf {

// ---------------- Basic types ----------------
struct Tag {
    char k1='?', k2='?';     // e.g. 'c','g' for cg
    char type='?';           // SAM type code: A/i/f/Z/H/...
    std::string_view val;    // value (without the "AA:T:" prefix)
};

struct Record {
    // Core 12 columns (GAF v1)
    std::string_view qname;
    uint32_t qlen=0, qstart=0, qend=0;
    char strand='?';
    std::string_view path;
    uint32_t tlen=0, tstart=0, tend=0;
    uint32_t n_match=0, aln_len=0;
    int32_t  mapq=-1;

    // Tags
    std::string_view cg;        // cg:Z:... (if present)
    std::vector<Tag> tags;      // other tags
};

struct PathStep {
    bool rev=false;             // '<' = reverse, '>' = forward
    std::string_view name;      // step name/id
};

struct CigarOp {
    char op='?';                // M,=,X,I,D,N,S,...
    uint32_t len=0;
};

// ---------------- Helpers ----------------

// Split a line into columns by '\t', mutate buf by writing NULs.
inline bool parse_line(char* buf, size_t len, Record& out, std::vector<Tag>& tmp_tags) {
    std::vector<char*> cols;
    cols.reserve(24);
    char* p = buf;
    char* e = buf + len;
    for (char* s = p; s <= e; ++p) {
        if (p == e || *p == '\t') {
            *p = '\0';
            cols.push_back(s);
            s = p + 1;
        }
    }
    if (cols.size() < 12) return false;

    auto to_u32 = [](const char* s)->uint32_t { return (uint32_t)std::strtoul(s, nullptr, 10); };
    auto to_i32 = [](const char* s)->int32_t  { return (int32_t) std::strtol(s, nullptr, 10); };

    out = Record{};
    out.qname   = std::string_view(cols[0]);
    out.qlen    = to_u32(cols[1]);
    out.qstart  = to_u32(cols[2]);
    out.qend    = to_u32(cols[3]);
    out.strand  = cols[4][0];
    out.path    = std::string_view(cols[5]);
    out.tlen    = to_u32(cols[6]);
    out.tstart  = to_u32(cols[7]);
    out.tend    = to_u32(cols[8]);
    out.n_match = to_u32(cols[9]);
    out.aln_len = to_u32(cols[10]);
    out.mapq    = to_i32(cols[11]);

    tmp_tags.clear();
    out.cg = std::string_view{};
    for (size_t i = 12; i < cols.size(); ++i) {
        const char* s = cols[i];
        size_t L = std::strlen(s);
        if (L < 5) continue;               // at least "AA:T:"
        if (s[2] != ':' || s[4] != ':') continue;

        Tag t;
        t.k1 = s[0]; t.k2 = s[1]; t.type = s[3];
        t.val = std::string_view(s + 5, L - 5);

        if (t.k1=='c' && t.k2=='g' && t.type=='Z') out.cg = t.val;
        else tmp_tags.push_back(t);
    }
    out.tags = tmp_tags; // copy views
    return true;
}

// Parse "path" like >52>53<88... into steps.
inline void parse_path(std::string_view path, std::vector<PathStep>& out) {
    out.clear();
    size_t i=0, n=path.size();
    while (i<n) {
        char c = path[i];
        if (c!='>' && c!='<') { ++i; continue; }
        bool rev = (c=='<');
        size_t j = i+1;
        while (j<n && path[j]!='>' && path[j]!='<' && path[j]!='\t' && path[j]!=' ') ++j;
        if (j > i+1) out.push_back({rev, path.substr(i+1, j-(i+1))});
        i = j;
    }
}

// Parse cg:Z:... into a vector of {len,op}.
inline void parse_cigar(std::string_view cg, std::vector<CigarOp>& out) {
    out.clear();
    uint32_t num = 0; bool have=false;
    for (char ch : cg) {
        if ((unsigned)(ch - '0') <= 9) { num = num*10 + (ch - '0'); have=true; }
        else {
            if (have) { out.push_back(CigarOp{ch,num}); num=0; have=false; }
        }
    }
}

// Zero-allocation iteration over CIGAR (useful for direct coverage accumulation)
template<class Fn>
inline void for_each_cigar(std::string_view cg, Fn&& fn) {
    uint32_t num = 0; bool have=false;
    for (char ch : cg) {
        if ((unsigned)(ch - '0') <= 9) { num = num*10 + (ch - '0'); have=true; }
        else {
            if (have) {
                if (!fn(num, ch)) return;
            }
            num=0; have=false;
        }
    }
}

// Trim trailing '\r' and/or '\n'
inline void rstrip_newline(std::string& s) {
    while (!s.empty() && (s.back()=='\n' || s.back()=='\r')) s.pop_back();
}

// Skip empty/comment lines; return false on EOF
template<class ReadLineFn, class ParseFn>
inline bool read_and_parse_loop(ReadLineFn&& readLine, ParseFn&& parse, Record& rec) {
    std::string* line = nullptr;
    while (true) {
        if (!readLine(line)) return false;
        if (!line || line->empty() || (*line)[0]=='#') continue;
        rstrip_newline(*line);
        if (line->empty()) continue;
        // parse in-place
        return parse(line->data(), line->size(), rec);
    }
}

// ---------------- Readers (explicit) ----------------

// 1) std::istream reader (no zlib)
class ReaderStream {
public:
    explicit ReaderStream(std::istream& is) : is_(is) {}
    bool next(Record& rec) {
        auto readLine = [&](std::string*& out)->bool {
            if (!std::getline(is_, buf_)) return false;
            out = &buf_; return true;
        };
        return read_and_parse_loop(
            readLine,
            [&](char* b, size_t l, Record& r){ tmp_tags_.clear(); return parse_line(b,l,r,tmp_tags_); },
            rec
        );
    }
private:
    std::istream& is_;
    std::string buf_;
    std::vector<Tag> tmp_tags_;
};

// 2) FILE* reader (fgets)
class ReaderFILE {
public:
    explicit ReaderFILE(FILE* fp, size_t chunk=65536) : fp_(fp), chunk_(chunk) { tmp_.resize(chunk_); }
    bool next(Record& rec) {
        auto readLine = [&](std::string*& out)->bool {
            buf_.clear();
            for (;;) {
                if (!std::fgets(tmp_.data(), (int)tmp_.size(), fp_)) {
                    if (buf_.empty()) return false; // EOF
                    break; // partial last line (no newline)
                }
                buf_.append(tmp_.data());
                if (!buf_.empty() && buf_.back()=='\n') break;
            }
            out = &buf_;
            return true;
        };
        return read_and_parse_loop(
            readLine,
            [&](char* b, size_t l, Record& r){ tmp_tags_.clear(); return parse_line(b,l,r,tmp_tags_); },
            rec
        );
    }
private:
    FILE* fp_;
    size_t chunk_;
    std::string buf_, tmp_;
    std::vector<Tag> tmp_tags_;
};

// 3) gzFile reader (zlib), using gzgets
class ReaderGz {
public:
    explicit ReaderGz(gzFile fp, size_t chunk=65536) : fp_(fp), chunk_(chunk) { tmp_.resize(chunk_); }
    bool next(Record& rec) {
        auto readLine = [&](std::string*& out)->bool {
            buf_.clear();
            for (;;) {
                if (!gzgets(fp_, tmp_.data(), (int)tmp_.size())) {
                    if (buf_.empty()) return false; // EOF
                    break; // partial last line
                }
                buf_.append(tmp_.data());
                if (!buf_.empty() && buf_.back()=='\n') break;
            }
            out = &buf_;
            return true;
        };
        return read_and_parse_loop(
            readLine,
            [&](char* b, size_t l, Record& r){ tmp_tags_.clear(); return parse_line(b,l,r,tmp_tags_); },
            rec
        );
    }
private:
    gzFile fp_;
    size_t chunk_;
    std::string buf_, tmp_;
    std::vector<Tag> tmp_tags_;
};

// ---------------- Auto-detecting reader ----------------

inline bool ends_with(const std::string& s, const char* suf) {
    size_t n = s.size(), m = std::strlen(suf);
    return (n>=m) && std::memcmp(s.data()+n-m, suf, m)==0;
}

// Content-based gzip detection (reads first 2 bytes)
inline bool file_is_gzip(const std::string& path) {
    FILE* fp = std::fopen(path.c_str(), "rb");
    if (!fp) return false;
    unsigned char mg[2] = {0,0};
    size_t n = std::fread(mg, 1, 2, fp);
    std::fclose(fp);
    return n == 2 && mg[0] == 0x1f && mg[1] == 0x8b;
}

// RAII reader that opens by path and auto-selects plain vs gzip.
class ReaderAuto {
public:
    // If path == "-", treat as plain stdin (most robust default).
    explicit ReaderAuto(const std::string& path, size_t chunk=65536) : chunk_(chunk) {
        open(path);
    }

    ~ReaderAuto() {
        cleanup();
    }

    // Return true if underlying handle is open
    bool good() const { return kind_ != K_FAIL; }

    // Read next record
    bool next(Record& rec) {
        switch (kind_) {
            case K_FILE: return r_file_->next(rec);
            case K_GZ  : return r_gz_->next(rec);
            default    : return false;
        }
    }

private:
    enum Kind { K_FAIL, K_FILE, K_GZ };

    void open(const std::string& path) {
        // stdin special case
        if (path == "-") {
            kind_ = K_FILE;
            fp_ = stdin;
            r_file_ = new ReaderFILE(fp_, chunk_);
            return;
        }

        // Try content-based detection first
        bool gz = file_is_gzip(path);

        // If content check failed to open/read, fall back to suffix hint
        if (!gz && ends_with(path, ".gz")) gz = true;

        if (gz) {
            gzfp_ = gzopen(path.c_str(), "rb");
            if (!gzfp_) { kind_ = K_FAIL; return; }
            r_gz_ = new ReaderGz(gzfp_, chunk_);
            kind_ = K_GZ;
        } else {
            fp_ = std::fopen(path.c_str(), "rb");
            if (!fp_) { kind_ = K_FAIL; return; }
            r_file_ = new ReaderFILE(fp_, chunk_);
            kind_ = K_FILE;
        }
    }

    void cleanup() {
        if (r_file_) { delete r_file_; r_file_ = nullptr; }
        if (r_gz_)   { delete r_gz_;   r_gz_   = nullptr; }
        if (fp_ && fp_ != stdin) { std::fclose(fp_); fp_ = nullptr; }
        if (gzfp_) { gzclose(gzfp_); gzfp_ = nullptr; }
        kind_ = K_FAIL;
    }

    Kind   kind_ = K_FAIL;
    size_t chunk_ = 65536;
    // underlying resources
    FILE*   fp_   = nullptr;
    gzFile  gzfp_ = nullptr;
    ReaderFILE* r_file_ = nullptr;
    ReaderGz*   r_gz_   = nullptr;
};

} // namespace kgaf