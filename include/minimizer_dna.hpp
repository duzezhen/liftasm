#pragma once
/*
 * - Similarity estimators on:
 *     * Jaccard: |A∩B| / |A∪B|
 *     * Containment: |A∩B| / |A|
 *     * Max-containment: max(|A∩B|/|A|, |A∩B|/|B|)
 *
 * Requirements
 * - 1 <= k <= 32 (2-bit rolling encode fits 64-bit)
 * - w >= 1
 */

#include <algorithm>
#include <cstdint>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace minimizerdna {

// SplitMix64 mixer
static inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x =  x ^ (x >> 31);
    return x;
}

// Map nucleotide to 2-bit code; return 4 for invalid (e.g., 'N')
static inline uint8_t nt4(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4;
    }
}

struct Options {
    uint32_t k         = 17;     // k-mer length (1..32)
    uint32_t w         = 50;     // window length in bases (w >= k)
    bool     canonical = true;   // use canonical k-mer (min(fwd, rc))
    uint64_t seed      = 0x8a5cd789635d2dffULL;
};

struct Sketch {
    std::vector<uint64_t> vals;  // sorted ascending, unique

    size_t  size() const { return vals.size(); }
    bool   empty() const { return vals.empty(); }
    void   clear()       { vals.clear(); }

    static size_t intersect_count(const Sketch& A, const Sketch& B) {
        const auto& a = A.vals;
        const auto& b = B.vals;
        size_t i = 0, j = 0, inter = 0;
        while (i < a.size() && j < b.size()) {
            if (a[i] < b[j]) ++i;
            else if (b[j] < a[i]) ++j;
            else { ++inter; ++i; ++j; }
        }
        return inter;
    }

    static size_t union_count(const Sketch& A, const Sketch& B, size_t inter) {
        return A.size() + B.size() - inter;
    }

    // |A∩B| / |A∪B|
    double jaccard(const Sketch& other) const {
        if (empty() || other.empty()) return 0.0;
        const size_t inter = intersect_count(*this, other);
        const size_t uni   = size() + other.size() - inter;
        return uni ? double(inter) / double(uni) : 0.0;
    }

    // Containment A in B: |A∩B| / |A|
    double containment_in(const Sketch& superset) const {
        if (empty() || superset.empty()) return 0.0;
        const size_t inter = intersect_count(*this, superset);
        return size() ? double(inter) / double(size()) : 0.0;
    }

    // Both containments: |A∩B|/|A|, |A∩B|/|B|
    std::pair<double,double> containments(const Sketch& other) const {
        if (empty() || other.empty()) return {0.0, 0.0};
        const size_t inter = intersect_count(*this, other);
        const double p_B_in_A = double(inter) / double(size());
        const double p_A_in_B = double(inter) / double(other.size());
        return {p_B_in_A, p_A_in_B};
    }

    // Max containment: max(|A∩B|/|A|, |A∩B|/|B|)
    double max_containment(const Sketch& other) const {
        if (empty() || other.empty()) return 0.0;
        auto [p_B_in_A, p_A_in_B] = containments(other);
        return (p_B_in_A > p_A_in_B) ? p_B_in_A : p_A_in_B;
    }
};

class MinimizerBuilder {
public:
    explicit MinimizerBuilder(Options opt = {}) : opt_(opt) {
        if (opt_.k == 0 || opt_.k > 32) {
            throw std::invalid_argument("k must be in [1,32]");
        }
        if (opt_.w == 0) {
            throw std::invalid_argument("w must be >= 1");
        }
        mask_ = (opt_.k == 32) ? ~uint64_t(0) : ((uint64_t(1) << (2 * opt_.k)) - 1);
    }

    const Options& options() const { return opt_; }

    Sketch build(std::string_view s) const {
        std::vector<uint64_t> mins;
        mins.reserve(s.size() / std::max<size_t>(opt_.w, 1));

        uint64_t fwd = 0, rev = 0;
        uint32_t span = 0;
        uint64_t kmer_idx = 0;

        struct Item { uint64_t idx; uint64_t h; };
        std::deque<Item> dq;

        auto reset_run = [&](){
            span = 0; fwd = 0; rev = 0;
            kmer_idx = 0;
            dq.clear();
        };

        auto push_kmer = [&](uint64_t idx, uint64_t h){
            while (!dq.empty() && dq.back().h >= h) dq.pop_back();
            dq.push_back({idx, h});
            const uint64_t low = (idx >= (opt_.w - 1)) ? (idx - (opt_.w - 1)) : 0;
            while (!dq.empty() && dq.front().idx < low) dq.pop_front();
            if (idx + 1 >= opt_.w && !dq.empty()) {
                mins.push_back(dq.front().h);
            }
        };

        for (char c : s) {
            uint8_t v = nt4(c);
            if (v < 4) {
                fwd = ((fwd << 2) | v) & mask_;
                if (opt_.canonical) {
                    uint64_t comp = uint64_t(3 - v);
                    rev = (rev >> 2) | (comp << ((opt_.k - 1) * 2));
                }
                if (++span >= opt_.k) {
                    const uint64_t key = opt_.canonical ? (fwd < rev ? fwd : rev) : fwd;
                    const uint64_t h = mix64(key ^ opt_.seed);
                    push_kmer(kmer_idx++, h);
                }
            } else {
                reset_run();
            }
        }

        // Sort and deduplication
        std::sort(mins.begin(), mins.end());
        mins.erase(std::unique(mins.begin(), mins.end()), mins.end());
        mins.shrink_to_fit();

        Sketch sk;
        sk.vals = std::move(mins);
        return sk;
    }

    void build(std::string_view s, Sketch& sk) const {
        std::vector<uint64_t> mins;
        mins.reserve(s.size() / std::max<size_t>(opt_.w, 1));

        uint64_t fwd = 0, rev = 0;
        uint32_t span = 0;
        uint64_t kmer_idx = 0;

        struct Item { uint64_t idx; uint64_t h; };
        std::deque<Item> dq;

        auto reset_run = [&](){
            span = 0; fwd = 0; rev = 0;
            kmer_idx = 0;
            dq.clear();
        };

        auto push_kmer = [&](uint64_t idx, uint64_t h){
            while (!dq.empty() && dq.back().h >= h) dq.pop_back();
            dq.push_back({idx, h});
            const uint64_t low = (idx >= (opt_.w - 1)) ? (idx - (opt_.w - 1)) : 0;
            while (!dq.empty() && dq.front().idx < low) dq.pop_front();
            if (idx + 1 >= opt_.w && !dq.empty()) {
                mins.push_back(dq.front().h);
            }
        };

        for (char c : s) {
            uint8_t v = nt4(c);
            if (v < 4) {
                fwd = ((fwd << 2) | v) & mask_;
                if (opt_.canonical) {
                    uint64_t comp = uint64_t(3 - v);
                    rev = (rev >> 2) | (comp << ((opt_.k - 1) * 2));
                }
                if (++span >= opt_.k) {
                    const uint64_t key = opt_.canonical ? (fwd < rev ? fwd : rev) : fwd;
                    const uint64_t h = mix64(key ^ opt_.seed);
                    push_kmer(kmer_idx++, h);
                }
            } else {
                reset_run();
            }
        }

        // merge, sort and deduplication
        sk.vals.insert(sk.vals.end(), mins.begin(), mins.end());
        std::sort(sk.vals.begin(), sk.vals.end());
        sk.vals.erase(std::unique(sk.vals.begin(), sk.vals.end()), sk.vals.end());
        sk.vals.shrink_to_fit();
        return;
    }

    double jaccard(std::string_view a, std::string_view b) const {
        Sketch A = build(a);
        Sketch B = build(b);
        return A.jaccard(B);
    }
    double containment(std::string_view Aseq, std::string_view Bseq) const {
        Sketch A = build(Aseq);
        Sketch B = build(Bseq);
        return A.containment_in(B);
    }
    std::pair<double,double> containments(std::string_view a, std::string_view b) const {
        Sketch A = build(a);
        Sketch B = build(b);
        return A.containments(B);
    }
    std::pair<double,double> containments(const Sketch& A, const Sketch& B) const {
        return A.containments(B);
    }
    double max_containment(std::string_view a, std::string_view b) const {
        Sketch A = build(a);
        Sketch B = build(b);
        return A.max_containment(B);
    }
    double max_containment(const Sketch& A, const Sketch& B) const { return A.max_containment(B); }

private:
    Options  opt_{};
    uint64_t mask_{0};
};

inline Sketch build_sketch(std::string_view s, const Options& opt = {}) {
    return MinimizerBuilder(opt).build(s);
}
inline double jaccard_estimate(std::string_view a, std::string_view b, const Options& opt = {}) {
    MinimizerBuilder mb(opt);
    return mb.jaccard(a, b);
}
inline double containment_estimate(std::string_view A, std::string_view B, const Options& opt = {}) {
    MinimizerBuilder mb(opt);
    return mb.containment(A, B);
}
inline std::pair<double,double> containments_estimate(std::string_view a, std::string_view b, const Options& opt = {}) {
    MinimizerBuilder mb(opt);
    return mb.containments(a, b);
}
inline double max_containment_estimate(std::string_view a, std::string_view b, const Options& opt = {}) {
    MinimizerBuilder mb(opt);
    return mb.max_containment(a, b);
}

} // namespace minimizerdna