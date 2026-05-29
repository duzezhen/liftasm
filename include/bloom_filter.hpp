#pragma once

/*
 *
 * Requirements for SketchT:
 *   - SketchT has member: std::vector<uint64_t> vals;
 *
 * Main API:
 *   fastbloom::BloomOptions opt;
 *   opt.bits = 1ULL << 27;   // 16 MiB
 *   opt.hashes = 4;
 *
 *   fastbloom::BloomSketch<minimizerdna::Sketch> x(seed_sketch, opt);
 *   x.add_sketch(next_node_sketch);
 *   auto c = x.bit_containments(y);
 *
 * Direction:
 *   c.first  = popcount(A & B) / popcount(A)
 *   c.second = popcount(A & B) / popcount(B)
 * 
 */

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <utility>
#include <vector>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

namespace fastbloom {

static inline uint64_t mix64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x =  x ^ (x >> 31);
    return x;
}

static inline uint64_t next_power_of_two_u64(uint64_t x)
{
    if (x <= 1) return 1;

    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;

    return x + 1;
}

struct BloomOptions {
    uint64_t bits = 1ULL << 27;   // 16 MiB default: 134,217,728 bits
    uint32_t hashes = 4;
    bool use_simd = true;
};

class Popcount {
public:
    static uint64_t words(const uint64_t* a, std::size_t n, bool use_simd)
    {
#if defined(__AVX2__)
        if (use_simd) return words_avx2_(a, n);
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
        if (use_simd) return words_neon_(a, n);
#endif
        return words_scalar_(a, n);
    }

    static uint64_t words_and(const uint64_t* a, const uint64_t* b, std::size_t n, bool use_simd)
    {
#if defined(__AVX2__)
        if (use_simd) return words_and_avx2_(a, b, n);
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
        if (use_simd) return words_and_neon_(a, b, n);
#endif
        return words_and_scalar_(a, b, n);
    }

private:
    static uint64_t words_scalar_(const uint64_t* a, std::size_t n)
    {
        uint64_t s = 0;

        for (std::size_t i = 0; i < n; ++i) {
#if defined(__GNUC__) || defined(__clang__)
            s += static_cast<uint64_t>(__builtin_popcountll(a[i]));
#else
            uint64_t x = a[i];
            while (x) {
                x &= x - 1;
                ++s;
            }
#endif
        }

        return s;
    }

    static uint64_t words_and_scalar_(const uint64_t* a, const uint64_t* b, std::size_t n)
    {
        uint64_t s = 0;

        for (std::size_t i = 0; i < n; ++i) {
            const uint64_t x = a[i] & b[i];

#if defined(__GNUC__) || defined(__clang__)
            s += static_cast<uint64_t>(__builtin_popcountll(x));
#else
            uint64_t y = x;
            while (y) {
                y &= y - 1;
                ++s;
            }
#endif
        }

        return s;
    }

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    static uint64_t sum_u8x16_(uint8x16_t v)
    {
#if defined(__aarch64__)
        return static_cast<uint64_t>(vaddvq_u8(v));
#else
        uint16x8_t s16 = vpaddlq_u8(v);
        uint32x4_t s32 = vpaddlq_u16(s16);
        uint64x2_t s64 = vpaddlq_u32(s32);
        return static_cast<uint64_t>(vgetq_lane_u64(s64, 0) + vgetq_lane_u64(s64, 1));
#endif
    }

    static uint64_t words_neon_(const uint64_t* a, std::size_t n)
    {
        const uint8_t* p = reinterpret_cast<const uint8_t*>(a);
        const std::size_t bytes = n * sizeof(uint64_t);

        uint64_t s = 0;
        std::size_t i = 0;

        for (; i + 16 <= bytes; i += 16) {
            const uint8x16_t x = vld1q_u8(p + i);
            s += sum_u8x16_(vcntq_u8(x));
        }

        for (; i < bytes; ++i) {
#if defined(__GNUC__) || defined(__clang__)
            s += static_cast<uint64_t>(__builtin_popcount(static_cast<unsigned>(p[i])));
#else
            uint8_t x = p[i];
            while (x) {
                x &= static_cast<uint8_t>(x - 1);
                ++s;
            }
#endif
        }

        return s;
    }

    static uint64_t words_and_neon_(const uint64_t* a, const uint64_t* b, std::size_t n)
    {
        const uint8_t* pa = reinterpret_cast<const uint8_t*>(a);
        const uint8_t* pb = reinterpret_cast<const uint8_t*>(b);
        const std::size_t bytes = n * sizeof(uint64_t);

        uint64_t s = 0;
        std::size_t i = 0;

        for (; i + 16 <= bytes; i += 16) {
            const uint8x16_t x = vandq_u8(vld1q_u8(pa + i), vld1q_u8(pb + i));
            s += sum_u8x16_(vcntq_u8(x));
        }

        for (; i < bytes; ++i) {
            const uint8_t x = static_cast<uint8_t>(pa[i] & pb[i]);

#if defined(__GNUC__) || defined(__clang__)
            s += static_cast<uint64_t>(__builtin_popcount(static_cast<unsigned>(x)));
#else
            uint8_t y = x;
            while (y) {
                y &= static_cast<uint8_t>(y - 1);
                ++s;
            }
#endif
        }

        return s;
    }
#endif

#if defined(__AVX2__)
    static uint64_t popcount_m256i_(__m256i v)
    {
        const __m256i lookup = _mm256_setr_epi8(
            0, 1, 1, 2, 1, 2, 2, 3,
            1, 2, 2, 3, 2, 3, 3, 4,
            0, 1, 1, 2, 1, 2, 2, 3,
            1, 2, 2, 3, 2, 3, 3, 4
        );

        const __m256i low_mask = _mm256_set1_epi8(0x0f);

        const __m256i lo = _mm256_and_si256(v, low_mask);
        const __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);

        const __m256i pc_lo = _mm256_shuffle_epi8(lookup, lo);
        const __m256i pc_hi = _mm256_shuffle_epi8(lookup, hi);
        const __m256i pc = _mm256_add_epi8(pc_lo, pc_hi);

        const __m256i zero = _mm256_setzero_si256();
        const __m256i sad = _mm256_sad_epu8(pc, zero);

        uint64_t tmp[4];
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(tmp), sad);

        return tmp[0] + tmp[1] + tmp[2] + tmp[3];
    }

    static uint64_t words_avx2_(const uint64_t* a, std::size_t n)
    {
        uint64_t s = 0;
        std::size_t i = 0;

        for (; i + 4 <= n; i += 4) {
            const __m256i x = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i));
            s += popcount_m256i_(x);
        }

        for (; i < n; ++i) {
            s += static_cast<uint64_t>(__builtin_popcountll(a[i]));
        }

        return s;
    }

    static uint64_t words_and_avx2_(const uint64_t* a, const uint64_t* b, std::size_t n)
    {
        uint64_t s = 0;
        std::size_t i = 0;

        for (; i + 4 <= n; i += 4) {
            const __m256i x = _mm256_and_si256(
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i)),
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + i))
            );

            s += popcount_m256i_(x);
        }

        for (; i < n; ++i) {
            s += static_cast<uint64_t>(__builtin_popcountll(a[i] & b[i]));
        }

        return s;
    }
#endif
};

class BitBloom {
public:
    BitBloom() = default;

    explicit BitBloom(const BloomOptions& opt)
    {
        reset(opt);
    }

    void reset(const BloomOptions& opt)
    {
        opt_ = opt;

        if (opt_.bits == 0) {
            opt_.bits = 1ULL << 27;
        }

        opt_.bits = next_power_of_two_u64(opt_.bits);

        if (opt_.hashes == 0) {
            opt_.hashes = 4;
        }

        opt_.hashes = std::max<uint32_t>(1, std::min<uint32_t>(opt_.hashes, 8));

        words_.assign(static_cast<std::size_t>(opt_.bits / 64), 0);

        mask_ = opt_.bits - 1;
        pop_valid_ = false;
        pop_cache_ = 0;
    }

    bool empty() const
    {
        return words_.empty();
    }

    uint64_t bit_size() const
    {
        return opt_.bits;
    }

    uint64_t byte_size() const
    {
        return opt_.bits / 8;
    }

    std::size_t word_size() const
    {
        return words_.size();
    }

    uint32_t hashes() const
    {
        return opt_.hashes;
    }

    void add(uint64_t h)
    {
        if (words_.empty()) return;

        uint64_t h1 = mix64(h ^ 0x9e3779b97f4a7c15ULL);
        uint64_t h2 = mix64(h ^ 0xbf58476d1ce4e5b9ULL);
        h2 |= 1ULL;

        for (uint32_t i = 0; i < opt_.hashes; ++i) {
            const uint64_t x = (h1 + uint64_t(i) * h2) & mask_;
            words_[static_cast<std::size_t>(x >> 6)] |= (1ULL << (x & 63));
        }

        pop_valid_ = false;
    }

    template<class VecLike>
    void add_values(const VecLike& vals)
    {
        for (uint64_t h : vals) {
            add(h);
        }
    }

    uint64_t popcount() const
    {
        if (!pop_valid_) {
            pop_cache_ = Popcount::words(words_.data(), words_.size(), opt_.use_simd);
            pop_valid_ = true;
        }

        return pop_cache_;
    }

    uint64_t intersection_popcount(const BitBloom& other) const
    {
        if (words_.empty() || other.words_.empty()) return 0;

        const std::size_t n = std::min(words_.size(), other.words_.size());

        return Popcount::words_and(
            words_.data(),
            other.words_.data(),
            n,
            opt_.use_simd && other.opt_.use_simd
        );
    }

private:
    BloomOptions opt_{};
    uint64_t mask_{0};

    std::vector<uint64_t> words_;

    mutable bool pop_valid_{false};
    mutable uint64_t pop_cache_{0};
};

template<class SketchT>
class BloomSketch {
public:
    BloomSketch() = default;

    BloomSketch(const SketchT& seed, const BloomOptions& opt)
        : bloom_(opt)
    {
        add_sketch(seed);
    }

    bool empty() const
    {
        return item_count_ == 0 || bloom_.empty();
    }

    uint64_t item_count() const
    {
        return item_count_;
    }

    uint64_t bit_popcount() const
    {
        return bloom_.popcount();
    }

    uint64_t bit_size() const
    {
        return bloom_.bit_size();
    }

    uint64_t byte_size() const
    {
        return bloom_.byte_size();
    }

    const BitBloom& bloom() const
    {
        return bloom_;
    }

    void add_sketch(const SketchT& sk)
    {
        if (sk.vals.empty()) return;

        bloom_.add_values(sk.vals);

        item_count_ += static_cast<uint64_t>(sk.vals.size());
    }

    std::pair<double, double> bit_containments(const BloomSketch& other) const
    {
        if (empty() || other.empty()) return {0.0, 0.0};

        const uint64_t inter = bloom_.intersection_popcount(other.bloom_);
        const uint64_t a_pop = bloom_.popcount();
        const uint64_t b_pop = other.bloom_.popcount();

        const double a_in_b = a_pop ? double(inter) / double(a_pop) : 0.0;
        const double b_in_a = b_pop ? double(inter) / double(b_pop) : 0.0;

        return {a_in_b, b_in_a};
    }

    double bit_jaccard(const BloomSketch& other) const
    {
        if (empty() || other.empty()) return 0.0;

        const uint64_t inter = bloom_.intersection_popcount(other.bloom_);
        const uint64_t a_pop = bloom_.popcount();
        const uint64_t b_pop = other.bloom_.popcount();

        const uint64_t uni = a_pop + b_pop - inter;
        return uni ? double(inter) / double(uni) : 0.0;
    }

private:
    BitBloom bloom_;
    uint64_t item_count_{0};
};

} // namespace fastbloom