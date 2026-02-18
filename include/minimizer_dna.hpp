#pragma once
/*
 * minimizer_dna.hpp — Header-only minimizer-based similarity for DNA strings.
 *
 * What this does
 * --------------
 * - Compute all minimizers from a DNA sequence using (k, w) scheme:
 *     * roll k-mers, hash them, and take the minimum hash in each window of size w
 *     * use a monotonic deque to get O(n) time per sequence
 * - Return the **set of unique minimizer hashes** for the whole sequence
 * - Provide similarity estimators on **minimizer sets**:
 *     * Jaccard: |A∩B| / |A∪B|
 *     * Containment: |A∩B| / |A|
 *     * Max-containment: max(|A∩B|/|A|, |A∩B|/|B|)
 *
 * Why minimizers?
 * ---------------
 * - Far fewer items than all k-mers (~every window contributes 1 minimizer)
 * - Robust, widely used in genomics indexers (minimap2, mashmap, etc.)
 * - Great memory savings vs full k-mer sets, yet more faithful than random subsampling
 *
 * Requirements
 * ------------
 * - C++17 or later
 * - 1 <= k <= 32 (2-bit rolling encode fits 64-bit)
 * - w >= k
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

// ---------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------

/// SplitMix64 mixer — fast 64-bit mixing for keys (good enough here).
static inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x =  x ^ (x >> 31);
    return x;
}

/// Map nucleotide to 2-bit code; return 4 for invalid (e.g., 'N')
static inline uint8_t nt4(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4;
    }
}

// ---------------------------------------------------------------
// Options
// ---------------------------------------------------------------

struct Options {
    uint32_t k         = 17;     // k-mer length (1..32)
    uint32_t w         = 50;     // window length in bases (w >= k)
    bool     canonical = true;   // use canonical k-mer (min(fwd, rc))
    uint64_t seed      = 0x8a5cd789635d2dffULL; // hash seed
};

// ---------------------------------------------------------------
// Sketch (sorted unique minimizer hashes)
// ---------------------------------------------------------------

struct Sketch {
    std::vector<uint64_t> vals; // sorted ascending, unique

    size_t size() const { return vals.size(); }
    bool   empty() const { return vals.empty(); }
    void   clear()       { vals.clear(); }

    // Two-pointer intersection cardinality
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

    // Union cardinality (via sizes and intersection)
    static size_t union_count(const Sketch& A, const Sketch& B, size_t inter) {
        return A.size() + B.size() - inter;
    }

    // Jaccard on minimizer sets: |A∩B| / |A∪B|
    double jaccard(const Sketch& other) const {
        if (empty() || other.empty()) return 0.0;
        const size_t inter = intersect_count(*this, other);
        const size_t uni   = size() + other.size() - inter;
        return uni ? double(inter) / double(uni) : 0.0; // uni>0 here
    }

    // Containment A in B: |A∩B| / |A|
    double containment_in(const Sketch& superset) const {
        if (empty() || superset.empty()) return 0.0;
        const size_t inter = intersect_count(*this, superset);
        return size() ? double(inter) / double(size()) : 0.0; // size()>0 here
    }

    // Both containments: (A in B, B in A)
    std::pair<double,double> containments(const Sketch& other) const {
        if (empty() || other.empty()) return {0.0, 0.0};
        const size_t inter = intersect_count(*this, other);
        const double pA = double(inter) / double(size());
        const double pB = double(inter) / double(other.size());
        return {pA, pB};
    }

    // Max containment: max(|A∩B|/|A|, |A∩B|/|B|)
    double max_containment(const Sketch& other) const {
        if (empty() || other.empty()) return 0.0;
        auto [pA, pB] = containments(other); // already non-empty
        return (pA > pB) ? pA : pB;
    }
};

// ---------------------------------------------------------------
// Minimizer builder
//  - O(n) time using a monotonic deque of (kmer_index, hash)
//  - We produce one minimizer per valid window (ties kept as they appear);
//    final set is deduplicated (sort+unique).
// ---------------------------------------------------------------

class MinimizerBuilder {
public:
    explicit MinimizerBuilder(Options opt = {}) : opt_(opt) {
        if (opt_.k == 0 || opt_.k > 32)
            throw std::invalid_argument("k must be in [1,32]");
        if (opt_.w < opt_.k)
            throw std::invalid_argument("w must be >= k");
        mask_ = (opt_.k == 32) ? ~uint64_t(0) : ((uint64_t(1) << (2 * opt_.k)) - 1);
        m_win_k_ = opt_.w - opt_.k + 1; // number of k-mers per window
    }

    const Options& options() const { return opt_; }

    /// Build a minimizer set from a whole string_view
    Sketch build(std::string_view s) const {
        std::vector<uint64_t> mins;
        mins.reserve(s.size() / std::max<size_t>(opt_.w, 1)); // rough guess

        // Rolling k-mer state
        uint64_t fwd = 0, rev = 0;
        uint32_t span = 0;            // how many valid bases in current run
        uint64_t kmer_idx = 0;        // index of current k-mer

        // Monotonic deque: pairs of (kmer_idx, hash)
        struct Item { uint64_t idx; uint64_t h; };
        std::deque<Item> dq;

        auto reset_run = [&](){
            span = 0; fwd = 0; rev = 0;
            kmer_idx = 0;
            dq.clear();
        };

        auto push_kmer = [&](uint64_t idx, uint64_t h){
            // Maintain deque order by hash (increasing)
            while (!dq.empty() && dq.back().h >= h) dq.pop_back();
            dq.push_back({idx, h});
            // Remove elements outside the window (older than idx - (m_win_k_-1))
            const uint64_t low = (idx >= (m_win_k_ - 1)) ? (idx - (m_win_k_ - 1)) : 0;
            while (!dq.empty() && dq.front().idx < low) dq.pop_front();
            // If we have at least m_win_k_ k-mers, we can report a window minimizer
            if (idx + 1 >= m_win_k_ && !dq.empty()) {
                mins.push_back(dq.front().h);
            }
        };

        for (char c : s) {
            uint8_t v = nt4(c);
            if (v < 4) {
                // forward roll: append 2-bit
                fwd = ((fwd << 2) | v) & mask_;
                if (opt_.canonical) {
                    // reverse-complement rolling: shift right-by-2; insert complement at top bits
                    uint64_t comp = uint64_t(3 - v);
                    rev = (rev >> 2) | (comp << ((opt_.k - 1) * 2));
                }
                if (++span >= opt_.k) {
                    const uint64_t key = opt_.canonical ? (fwd < rev ? fwd : rev) : fwd;
                    const uint64_t h = mix64(key ^ opt_.seed);
                    push_kmer(kmer_idx++, h);
                }
            } else {
                // Break run at 'N' or invalid base
                reset_run();
            }
        }

        // Deduplicate minimizers to a set
        std::sort(mins.begin(), mins.end());
        mins.erase(std::unique(mins.begin(), mins.end()), mins.end());

        Sketch sk;
        sk.vals = std::move(mins);
        return sk;
    }

    // Convenience estimators on raw sequences
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
    double max_containment(std::string_view a, std::string_view b) const {
        Sketch A = build(a);
        Sketch B = build(b);
        return A.max_containment(B);
    }

private:
    Options  opt_{};
    uint64_t mask_{0};
    uint64_t m_win_k_{0}; // #k-mers per base-window (w - k + 1)
};

// ---------------------------------------------------------------
// Top-level helpers (pure functions), if you don't want to keep objects
// ---------------------------------------------------------------

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