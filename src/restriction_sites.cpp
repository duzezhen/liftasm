#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <future>
#include <stdexcept>
#include <utility>

#include <zlib.h>

#include "../include/restriction_sites.hpp"
#include "../include/kseq.h"
#include "../include/logger.hpp"
#include "../include/ThreadPool.hpp"
#include "../include/progress_tracker.hpp"

#ifndef PROJECT_KSEQ_GZ_INIT_DONE
#define PROJECT_KSEQ_GZ_INIT_DONE
KSEQ_INIT(gzFile, gzread)
#endif

namespace rsite {
namespace detail {

// A=1, C=2, G=4, T=8
static inline uint8_t iupac_mask(char c) {
    switch ((char)std::toupper((unsigned char)c)) {
        case 'A': return 0x1;
        case 'C': return 0x2;
        case 'G': return 0x4;
        case 'T': return 0x8;
        case 'R': return 0x1 | 0x4;             // A/G
        case 'Y': return 0x2 | 0x8;             // C/T
        case 'S': return 0x2 | 0x4;             // C/G
        case 'W': return 0x1 | 0x8;             // A/T
        case 'K': return 0x4 | 0x8;             // G/T
        case 'M': return 0x1 | 0x2;             // A/C
        case 'B': return 0x2 | 0x4 | 0x8;       // C/G/T
        case 'D': return 0x1 | 0x4 | 0x8;       // A/G/T
        case 'H': return 0x1 | 0x2 | 0x8;       // A/C/T
        case 'V': return 0x1 | 0x2 | 0x4;       // A/C/G
        case 'N': return 0x1 | 0x2 | 0x4 | 0x8; // A/C/G/T
        default:  return 0;
    }
}

static inline uint8_t base_mask(char c) {
    switch ((char)std::toupper((unsigned char)c)) {
        case 'A': return 0x1;
        case 'C': return 0x2;
        case 'G': return 0x4;
        case 'T': return 0x8;
        default:  return 0;
    }
}

static inline char dna_comp(char c) {
    switch ((char)std::toupper((unsigned char)c)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'S': return 'S';
        case 'W': return 'W';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'B': return 'V';
        case 'D': return 'H';
        case 'H': return 'D';
        case 'V': return 'B';
        case 'N': return 'N';
        default:  return 'N';
    }
}

static inline std::string revcomp_iupac(std::string_view s) {
    std::string out;
    out.resize(s.size());
    for (size_t i = 0, n = s.size(); i < n; ++i) {
        out[i] = dna_comp(s[n - 1 - i]);
    }
    return out;
}

static inline std::string to_upper(std::string s) {
    for (char& c : s) c = (char)std::toupper((unsigned char)c);
    return s;
}

static inline bool match_mask(const char* seq, const std::vector<uint8_t>& motif_mask) {
    for (size_t i = 0; i < motif_mask.size(); ++i) {
        uint8_t b = base_mask(seq[i]);
        if ((b & motif_mask[i]) == 0) return false;
    }
    return true;
}

static inline std::string remove_caret(std::string_view s, uint32_t& cut_offset) {
    std::string out;
    out.reserve(s.size());

    bool seen = false;
    cut_offset = 0;

    for (char c : s) {
        if (c == '^') {
            if (seen) {
                error_stream() << "enzyme spec has multiple '^': " << std::string(s) << std::endl;
                std::exit(1);
            }
            seen = true;
            cut_offset = (uint32_t)out.size();
        } else {
            out.push_back((char)std::toupper((unsigned char)c));
        }
    }

    if (!seen) {
        error_stream() << "enzyme spec must contain '^': " << std::string(s) << std::endl;
        std::exit(1);
    }
    if (out.empty()) {
        error_stream() << "empty enzyme pattern: " << std::string(s) << std::endl;
        std::exit(1);
    }

    for (char c : out) {
        if (iupac_mask(c) == 0) {
            error_stream() << "invalid IUPAC base in enzyme spec: " << std::string(s) << std::endl;
            std::exit(1);
        }
    }
    return out;
}

static inline const std::unordered_map<std::string, std::string>& enzyme_name_db() {
    static const std::unordered_map<std::string, std::string> db = {
        {"A^AGCTT",   "HindIII"},
        {"^GATC",     "MboI/DpnII/Sau3AI"},
        {"A^GATCT",   "BglII"},
        {"CATG^",     "NlaIII"},
        {"G^ANTC",    "HinfI"},
        {"AG^CT",     "AluI"},
        {"C^CATGG",   "NcoI"},
        {"G^AATTC",   "EcoRI"},
        {"G^GATCC",   "BamHI"},
        {"CTGCA^G",   "PstI"},
        {"G^TCGAC",   "SalI"},
        {"A^CTAGT",   "SpeI"},
        {"T^CTAGA",   "XbaI"},
        {"CCC^GGG",   "SmaI"},
        {"GC^GGCCGC", "NotI"},
        {"GGTAC^C",   "KpnI"},
        {"CA^TATG",   "NdeI"},
        {"GAT^ATC",   "EcoRV"},
        {"C^TCGAG",   "XhoI"},
        {"T^CGA",     "TaqI"},
        {"C^CGG",     "MspI"},
        {"CCTGCA^GG", "SbfI"}
    };
    return db;
}

static inline CompiledEnzyme compile_enzyme(const Enzyme& ez) {
    CompiledEnzyme c;
    c.name = ez.name;
    c.spec = ez.spec;
    c.pattern = ez.pattern;
    c.cut_offset = ez.cut_offset;
    c.len = (uint32_t)ez.pattern.size();
    c.rev_cut_offset = c.len - c.cut_offset;

    c.fwd_mask.resize(c.len);
    for (size_t i = 0; i < ez.pattern.size(); ++i) {
        c.fwd_mask[i] = iupac_mask(ez.pattern[i]);
    }

    std::string rc = revcomp_iupac(ez.pattern);
    c.rev_mask.resize(c.len);
    for (size_t i = 0; i < rc.size(); ++i) {
        c.rev_mask[i] = iupac_mask(rc[i]);
    }
    return c;
}

}

namespace {

class FastaReader {
public:
    bool load(const std::string& path) {
        names.clear();
        seqs.clear();

        gzFile fp = gzopen(path.c_str(), "rb");
        if (!fp) {
            error_stream() << path << ": No such file or directory" << std::endl;
            std::exit(1);
        }

        kseq_t* ks = kseq_init(fp);
        if (!ks) {
            gzclose(fp);
            error_stream() << "kseq_init failed" << std::endl;
            std::exit(1);
        }

        while (kseq_read(ks) >= 0) {
            if (!ks->name.s || ks->name.l == 0) continue;

            std::string name(ks->name.s, ks->name.l);
            std::string s;
            s.assign(ks->seq.s ? ks->seq.s : "", (size_t)ks->seq.l);
            for (char& c : s) c = (char)std::toupper((unsigned char)c);

            names.push_back(name);
            seqs.emplace_back(std::move(name), std::move(s));
        }

        kseq_destroy(ks);
        gzclose(fp);
        return true;
    }

    std::vector<std::string> names;
    std::vector<std::pair<std::string, std::string>> seqs; // (name, seq)
};

}

Enzyme parse_enzyme_spec(const std::string& spec, size_t order_1based) {
    Enzyme ez;
    ez.spec = detail::to_upper(spec);
    ez.pattern = detail::remove_caret(ez.spec, ez.cut_offset);

    auto it = detail::enzyme_name_db().find(ez.spec);
    if (it != detail::enzyme_name_db().end()) {
        ez.name = it->second;
    } else {
        ez.name = "Enzyme" + std::to_string(order_1based);
    }
    return ez;
}

std::vector<Enzyme> parse_enzyme_specs(const std::vector<std::string>& specs) {
    std::vector<Enzyme> v;
    v.reserve(specs.size());
    for (size_t i = 0; i < specs.size(); ++i) {
        v.push_back(parse_enzyme_spec(specs[i], i + 1));
    }
    return v;
}

void Index::clear() {
    enzymes_.clear();
    chrom_order_.clear();
    chrom_len_.clear();
    chrom_sites_.clear();
    chrom_cuts_.clear();
}

void Index::build(const ResSitOpts& params)
{
    clear();

    log_stream() << "Building restriction site index from FASTA: " << params.genome_file << "\n";

    if (params.enzymes.empty()) {
        error_stream() << "no enzyme specs provided" << std::endl;
        std::exit(1);
    }

    auto parsed = parse_enzyme_specs(params.enzymes);
    enzymes_.reserve(parsed.size());
    for (const auto& ez : parsed) {
        enzymes_.push_back(detail::compile_enzyme(ez));
    }

    FastaReader fr;
    fr.load(params.genome_file);

    chrom_order_.reserve(fr.seqs.size());
    for (const auto& kv : fr.seqs) {
        chrom_order_.push_back(kv.first);
    }

    ProgressTracker prog(fr.seqs.size(), 10);

    ThreadPool pool(std::max(1, params.threads));
    std::vector<std::future<ChromResult>> futs;
    futs.reserve(fr.seqs.size());
    const bool scan_revcomp = params.scan_revcomp;

    for (const auto& kv : fr.seqs) {
        prog.hit();

        const std::string chr = kv.first;
        const std::string seq = kv.second;

        futs.emplace_back(
            pool.submit([this, chr, seq, scan_revcomp]() {
                return scan_one_chrom_(chr, seq, scan_revcomp);
            })
        );
    }

    for (auto& f : futs) {
        ChromResult r = f.get();
        chrom_len_[r.chr]   = r.chr_len;
        chrom_sites_[r.chr] = std::move(r.sites);
        chrom_cuts_[r.chr]  = std::move(r.cuts);
    }

    pool.stop();

    size_t total_sites = 0;
    size_t chrom_with_sites = 0;
    for (const auto& chr : chrom_order_) {
        auto it = chrom_cuts_.find(chr);
        if (it != chrom_cuts_.end() && !it->second.empty()) {
            ++chrom_with_sites;
            total_sites += it->second.size();
        }
    }

    log_stream() << "Restriction site index built\n";
    log_stream() << "  - Chromosomes: " << chrom_order_.size() << "\n";
    log_stream() << "  - Enzymes: " << enzymes_.size() << "\n";
    log_stream() << "  - Chromosomes with sites: " << chrom_with_sites << "\n";
    log_stream() << "  - Unique cut positions: " << total_sites << "\n";
}

size_t Index::n_chrom() const {
    return chrom_order_.size();
}

size_t Index::n_enzyme() const {
    return enzymes_.size();
}

const std::vector<std::string>& Index::chrom_order() const {
    return chrom_order_;
}

const std::vector<uint32_t>* Index::cuts(std::string_view chr) const {
    auto it = chrom_cuts_.find(std::string(chr));
    if (it == chrom_cuts_.end()) return nullptr;
    return &it->second;
}

const std::vector<Site>* Index::sites(std::string_view chr) const {
    auto it = chrom_sites_.find(std::string(chr));
    if (it == chrom_sites_.end()) return nullptr;
    return &it->second;
}

Nearest Index::nearest(std::string_view chr, uint32_t pos) const {
    Nearest ans;
    auto it = chrom_cuts_.find(std::string(chr));
    if (it == chrom_cuts_.end()) return ans;

    const auto& v = it->second;
    if (v.empty()) return ans;

    auto lb = std::lower_bound(v.begin(), v.end(), pos);

    if (lb == v.begin()) {
        ans.found = true;
        ans.pos   = *lb;
        ans.dist  = (int64_t)ans.pos - (int64_t)pos;
        return ans;
    }
    if (lb == v.end()) {
        ans.found = true;
        ans.pos   = v.back();
        ans.dist  = (int64_t)ans.pos - (int64_t)pos;
        return ans;
    }

    uint32_t r = *lb;
    uint32_t l = *(lb - 1);

    uint64_t dl = (l > pos) ? (uint64_t)(l - pos) : (uint64_t)(pos - l);
    uint64_t dr = (r > pos) ? (uint64_t)(r - pos) : (uint64_t)(pos - r);

    ans.found = true;
    if (dl <= dr) {
        ans.pos  = l;
        ans.dist = (int64_t)l - (int64_t)pos;
    } else {
        ans.pos  = r;
        ans.dist = (int64_t)r - (int64_t)pos;
    }
    return ans;
}

QueryPair Index::nearest(std::string_view chr, uint32_t beg, uint32_t end) const {
    QueryPair q;
    q.beg_site = nearest(chr, beg);
    q.end_site = nearest(chr, end);
    return q;
}

QueryPair Index::flanking(std::string_view chr, uint32_t beg, uint32_t end) const {
    QueryPair q;
    auto it = chrom_cuts_.find(std::string(chr));
    if (it == chrom_cuts_.end()) return q;

    const auto& v = it->second;
    if (v.empty()) return q;

    auto it1 = std::upper_bound(v.begin(), v.end(), beg);
    if (it1 != v.begin()) {
        uint32_t x = *(it1 - 1);
        q.beg_site.found = true;
        q.beg_site.pos   = x;
        q.beg_site.dist  = (int64_t)x - (int64_t)beg;
    }

    auto it2 = std::lower_bound(v.begin(), v.end(), end);
    if (it2 != v.end()) {
        uint32_t x = *it2;
        q.end_site.found = true;
        q.end_site.pos   = x;
        q.end_site.dist  = (int64_t)x - (int64_t)end;
    }

    return q;
}

void Index::save_bed(const std::string& out_bed) const {
    std::ofstream os(out_bed);
    if (!os) {
        error_stream() << out_bed << ": No such file or directory" << std::endl;
        std::exit(1);
    }

    for (const auto& chr : chrom_order_) {
        auto itc = chrom_cuts_.find(chr);
        auto its = chrom_sites_.find(chr);
        if (itc == chrom_cuts_.end() || its == chrom_sites_.end()) continue;

        for (const auto& s : its->second) {
            const auto& ez = enzymes_.at(s.enzyme_id);
            os << chr << '\t'
                << s.pos << '\t'
                << (s.pos + 1) << '\t'
                << (ez.name.empty() ? ez.spec : ez.name) << '\t'
                << s.strand << '\n';
        }
    }
}

ChromResult Index::scan_one_chrom_(const std::string& chr, const std::string& seq, bool scan_revcomp) const
{
    ChromResult out;
    out.chr = chr;
    out.chr_len = (uint32_t)seq.size();

    const size_t n = seq.size();

    for (uint16_t eid = 0; eid < enzymes_.size(); ++eid) {
        const auto& ez = enzymes_[eid];
        if (n < ez.len) continue;

        for (size_t i = 0; i + ez.len <= n; ++i) {
            const char* s = seq.data() + i;

            if (detail::match_mask(s, ez.fwd_mask)) {
                uint64_t cut = (uint64_t)i + (uint64_t)ez.cut_offset;
                if (cut <= n) {
                    Site x;
                    x.pos = (uint32_t)cut;
                    x.enzyme_id = eid;
                    x.strand = '+';
                    out.sites.push_back(x);
                    out.cuts.push_back(x.pos);
                }
            }

            if (scan_revcomp) {
                if (detail::match_mask(s, ez.rev_mask)) {
                    uint64_t cut = (uint64_t)i + (uint64_t)ez.rev_cut_offset;
                    if (cut <= n) {
                        Site x;
                        x.pos = (uint32_t)cut;
                        x.enzyme_id = eid;
                        x.strand = '-';
                        out.sites.push_back(x);
                        out.cuts.push_back(x.pos);
                    }
                }
            }
        }
    }

    std::sort(out.sites.begin(), out.sites.end(),
              [](const Site& a, const Site& b) {
                  if (a.pos != b.pos) return a.pos < b.pos;
                  if (a.enzyme_id != b.enzyme_id) return a.enzyme_id < b.enzyme_id;
                  return a.strand < b.strand;
              });

    std::sort(out.cuts.begin(), out.cuts.end());
    out.cuts.erase(std::unique(out.cuts.begin(), out.cuts.end()), out.cuts.end());

    return out;
}

} // namespace rsite
