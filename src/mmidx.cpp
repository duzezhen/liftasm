#include "../include/mmidx.hpp"
#include "./include/ThreadPool.hpp"
#include "../include/kseq.h"
#include "../include/progress_tracker.hpp"
#include "../include/gfa_parser.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <numeric>
#include <zlib.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

KSEQ_INIT(gzFile, gzread)

/* ---------- file-local utilities --------------------------------------- */
namespace {
inline uint8_t nt4(char c) noexcept {
    static const std::array<uint8_t,256> lut = []{
        std::array<uint8_t,256> t{}; t.fill(4);
        t['A']=t['a']=0; t['C']=t['c']=1; t['G']=t['g']=2; t['T']=t['t']=3;
        return t;
    }();
    return lut[static_cast<uint8_t>(c)];
}
inline uint64_t splitmix64(uint64_t x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
} // anon


namespace mmidx {

/* ------------------------------------ Constructor ------------------------------------ */
MinimizerIndex::MinimizerIndex(
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs,
    const std::vector<std::vector<std::string>>& right_seqs,
    opt::ChainOpts& chainOpts, 
    const opt::AnchorOpts& anchorOpts
) : names_(names), seqs_(seqs), right_seqs_(right_seqs), chainOpts_(chainOpts), anchorOpts_(anchorOpts) {
    if (chainOpts_.k == 0 || chainOpts_.k > 64)   throw std::invalid_argument("k must be 1..64");
    if (chainOpts_.w == 0 || chainOpts_.w >= 256) throw std::invalid_argument("w must be 1..255");
    mask_   = (__uint128_t(1) << (2 * chainOpts_.k)) - 1;
    shift1_ = 2 * (chainOpts_.k - 1);
}

/* ------------------------------------ Helpers ------------------------------------ */
uint64_t MinimizerIndex::hash128_(__uint128_t v) const
{ return splitmix64(uint64_t(v) ^ uint64_t(v >> 64)); }

uint32_t seed_q_fwd(
    const uint32_t q_off,
    bool chain_rev,
    const uint32_t& read_len,
    const uint32_t& k
) noexcept {
    if (!chain_rev) return q_off;
    return read_len - (q_off + k);
}

template<class Emit>
void MinimizerIndex::for_each_kmer_(
    std::string_view s,
    std::size_t k,
    __uint128_t mask,
    uint32_t shift1,
    Emit&& emit
) const {
    uint32_t span = 0;
    __uint128_t fwd = 0, rev = 0;

    for (uint32_t i = 0; i < s.size(); ++i) {
        uint8_t c = nt4(s[i]);
        if (c < 4) {
            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | (__uint128_t(3 - c) << shift1);

            if (++span >= k) {
                if (fwd == rev) continue;

                bool fwd_small = (fwd < rev);
                __uint128_t key = fwd_small ? fwd : rev;
                uint8_t dir     = fwd_small ? 0u  : 1u;

                uint32_t q_off = i + 1 - k;
                emit(key, q_off, dir);
            }
        } else {
            span = 0; fwd = 0; rev = 0;
        }
    }
}

template<class Emit>
void MinimizerIndex::for_each_minimizer_(
    std::string_view s,
    std::size_t k,
    std::size_t w,
    __uint128_t mask,
    uint32_t shift1,
    Emit&& emit
) const {
    std::array<Cand,256> buf{};
    uint8_t head = 0, tail = 0;
    uint32_t last_emit = ~0u, span = 0;
    __uint128_t fwd = 0, rev = 0;

    for (uint32_t i = 0; i < s.size(); ++i) {
        uint8_t c = nt4(s[i]);
        if (c < 4) {
            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | (__uint128_t(3 - c) << shift1);
            if (++span >= k) {
                if (fwd == rev) continue;

                bool fwd_small = (fwd < rev);
                __uint128_t key = fwd_small ? fwd : rev;
                uint8_t dir     = fwd_small ? 0u  : 1u;
                uint64_t h      = hash128_(key);

                while (head != tail && buf[uint8_t(tail - 1)].hash >= h)
                    tail = uint8_t(tail - 1);
                buf[tail] = {h, i, key, dir}; tail = uint8_t(tail + 1);

                while (head != tail && buf[head].off + w <= i)
                    head = uint8_t(head + 1);

                if (i + 1 >= k && i + 1 >= w) {
                // if (i + 1 >= k && i + 1 >= w && head != tail) {
                    const auto& m = buf[head];
                    if (m.off != last_emit) {
                        emit(m.key, m.off - k + 1, m.dir);
                        last_emit = m.off;
                    }
                }
            }
        } else { span = 0; fwd = rev = 0; head = tail = 0; }
    }
}


void MinimizerIndex::build_mm(bool expand_right) {
    if (names_.size() != seqs_.size()) {
        error_stream() << "Size mismatch between names and seqs.\n";
        std::exit(1);
    }
    if (expand_right && right_seqs_.size() != names_.size()) {
        error_stream() << "Size mismatch among names, seqs, and right_seqs.\n";
        std::exit(1);
    }

    if (chainOpts_.w < 1) {
        error_stream() << "Window size w must be at least 1.\n";
        std::exit(1);
    }
    if (chainOpts_.k < 1 || chainOpts_.k > 64) {
        error_stream() << "k must be in the range 1..64.\n";
        std::exit(1);
    }

    std::vector<MiniRec> minis;
    size_t est = 0;
    for (auto const &sv : seqs_) est += sv.size();
    minis.reserve(est / chainOpts_.w + 128);

    const size_t n_seqs    = names_.size();
    const size_t n_threads = chainOpts_.threads ? size_t(chainOpts_.threads)
                                                : std::max<size_t>(1, std::thread::hardware_concurrency());
    const size_t n_tasks   = n_threads;
    const size_t chunk     = (n_seqs + n_tasks - 1) / n_tasks;

    const uint32_t k    = static_cast<uint32_t>(chainOpts_.k);
    const uint32_t need = (k > 0 ? k - 1 : 0);

    ThreadPool pool(n_threads);
    std::vector<std::future<std::vector<MiniRec>>> futs;
    futs.reserve(n_tasks);

    for (size_t t = 0; t < n_tasks; ++t) {
        const uint32_t beg = uint32_t(t * chunk);
        const uint32_t end = uint32_t(std::min(n_seqs, size_t(beg + chunk)));
        if (beg >= end) break;

        futs.emplace_back(
            pool.submit([&, beg, end, need, expand_right]() -> std::vector<MiniRec> {
                std::vector<MiniRec> loc;
                size_t local_est = ((est / chainOpts_.w) / std::max<size_t>(1, n_tasks)) + 128;
                loc.reserve(local_est);

                for (uint32_t seq_id = beg; seq_id < end; ++seq_id) {
                    const auto &seqv = seqs_[seq_id];

                    // 1) collect minimizers from the main sequence
                    for_each_minimizer_(
                        seqv, chainOpts_.k, chainOpts_.w, mask_, shift1_,
                        [&]( __uint128_t key, uint32_t off, uint8_t dir ){
                            uint32_t packed = (off << 1) | dir;
                            loc.push_back({ key, ( (uint64_t)seq_id << 32 ) | packed });
                        }
                    );

                    // 2) collect right extensions
                    if (!expand_right) continue;
                    if (need == 0 || seqv.empty()) continue;

                    const auto &rex = right_seqs_[seq_id];
                    if (rex.empty()) continue;

                    const uint32_t core_len = (uint32_t)seqv.size();
                    const uint32_t tail_len = std::min<uint32_t>(need, core_len);
                    if (tail_len == 0) continue;

                    // k-1
                    const std::string tail(seqv.substr(core_len - tail_len));

                    for (const std::string_view &r : rex) {
                        if (r.empty()) continue;
                        std::string r_copy(r); 

                        // tail + r  (length ≤ 2k-2)
                        std::string tmp;
                        tmp.reserve(tail_len + r_copy.size());
                        tmp.append(tail);
                        tmp.append(r_copy);

                        // Only minimizers whose starting point falls at the tail of the main sequence (off < tail_len) are kept
                        const uint32_t base_off = core_len - tail_len;
                        for_each_kmer_(
                            tmp, chainOpts_.k, mask_, shift1_,
                            [&]( __uint128_t key, uint32_t off, uint8_t dir ){
                                if (off >= tail_len) return;
                                uint32_t g_off  = base_off + off;
                                uint32_t packed = (g_off << 1) | dir;
                                loc.push_back({ key, ( (uint64_t)seq_id << 32 ) | packed });
                            }
                        );
                    }
                }
                return loc;
            })
        );
    }

    // Collect results
    for (auto &f : futs) {
        std::vector<MiniRec> loc = f.get();
        minis.insert(minis.end(),
                     std::make_move_iterator(loc.begin()),
                     std::make_move_iterator(loc.end()));
    }
    pool.stop();

    // Sort + compress + build hash
    radix_sort128_(minis);

    keys_.clear(); offs_.clear(); positions_.clear();
    offs_.push_back(0);

    for (size_t i = 0; i < minis.size(); ) {
        size_t j = i + 1;
        while (j < minis.size() && minis[j].key == minis[i].key) ++j;

        keys_.push_back(minis[i].key);
        offs_.push_back(uint32_t(positions_.size() + (j - i)));
        for (size_t k = i; k < j; ++k) positions_.push_back(minis[k].pos);
        i = j;
    }

    build_hash_();

    // Set mid_occ
    chainOpts_.mid_occ = calc_max_occ(chainOpts_.mid_occ_frac);
}


/* ------------------------------------ Sort ------------------------------------ */
void MinimizerIndex::radix_sort128_(std::vector<MiniRec>& v)
{
    constexpr int  BITS = 16;
    constexpr size_t B  = 1u << BITS;
    static thread_local std::vector<MiniRec> tmp;
    tmp.resize(v.size());
    static thread_local std::array<uint32_t,B> cnt;

    auto pass = [&](int shift, bool hi) {
        cnt.fill(0);
        for (const auto& r : v) {
            uint64_t part = hi ? uint64_t(r.key >> 64) : uint64_t(r.key);
            ++cnt[(part >> shift) & (B - 1)];
        }
        uint32_t sum = 0;
        for (auto& c : cnt) { uint32_t t = c; c = sum; sum += t; }
        for (const auto& r : v) {
            uint64_t part = hi ? uint64_t(r.key >> 64) : uint64_t(r.key);
            tmp[cnt[(part >> shift) & (B - 1)]++] = r;
        }
        v.swap(tmp);
    };
    for (int s = 0; s < 64; s += BITS) pass(s, false);
    for (int s = 0; s < 64; s += BITS) pass(s, true);
}

/* ------------------------------------ build open-address hash ------------------------------------ */
void MinimizerIndex::build_hash_() {
    if (keys_.empty()) { ht_.clear(); ht_mask_ = 0; return; }
    uint64_t sz = 1ULL; while (sz < keys_.size() * 2) sz <<= 1;
    ht_.assign(sz, 0u); ht_mask_ = sz - 1;

    for (uint32_t i = 0; i < keys_.size(); ++i) {
        uint64_t p = hash128_(keys_[i]) & ht_mask_;
        while (ht_[p]) p = (p + 1) & ht_mask_;
        ht_[p] = i + 1;
    }
}

// ----------------------- High-frequency minimizer filtering -----------------------
uint32_t MinimizerIndex::calc_max_occ(float f) const {
    // If no filtering requested or no keys, return
    if (f <= 0.0f || keys_.empty()) return std::numeric_limits<uint32_t>::max();

    // Build the per-key occurrence array: occ_i = offs_[i+1] - offs_[i]
    std::vector<uint32_t> occ;
    occ.reserve(keys_.size());
    for (size_t i = 0; i < keys_.size(); ++i) {
        uint32_t c = offs_[i + 1] - offs_[i];
        occ.push_back(c);
    }
    if (occ.empty()) return std::numeric_limits<uint32_t>::max();

    // (1 - f) quantile, clamped into [0, n-1]
    const size_t n = occ.size();
    size_t kth = static_cast<size_t>((1.0 - double(f)) * double(n));
    if (kth >= n) kth = n - 1;

    // nth_element puts the kth element into its final (sorted) position
    std::nth_element(occ.begin(), occ.begin() + kth, occ.end());
    uint32_t thres = occ[kth] + 1;

    return thres;
}

/* ------------------------------------ Query ------------------------------------ */
span_pos MinimizerIndex::query(std::string_view kmer) const {
    if (kmer.size() != chainOpts_.k) {
        error_stream() << "k-mer size mismatch: expected " << chainOpts_.k << ", got " << kmer.size() << '\n';
        return {};
    }
    return lookup_hash(canonical_key(kmer));
}

span_pos MinimizerIndex::lookup_hash(__uint128_t key) const {
    if (ht_.empty()) return {};
    uint64_t p = hash128_(key) & ht_mask_;
    for (;; p = (p + 1) & ht_mask_) {
        uint32_t t = ht_[p];
        if (t == 0) return {};
        uint32_t idx = t - 1;
        if (keys_[idx] == key)
            return { positions_.data() + offs_[idx], positions_.data() + offs_[idx + 1] };
    }
}

/* ------------------------------------ canonical_key ------------------------------------ */
__uint128_t MinimizerIndex::canonical_key(std::string_view s) const {
    __uint128_t f = 0, r = 0;
    for (char ch : s) {
        uint8_t x = nt4(ch);
        if (x >= 4) throw std::invalid_argument("invalid base");
        f = ((f << 2) | x) & mask_;
        r = (r >> 2) | (__uint128_t(3 - x) << shift1_);
    }
    return f < r ? f : r;
}

/* ------------------------------------ collect_seeds ------------------------------------ */
std::vector<MM128> MinimizerIndex::collect_seeds(
    std::string_view read, bool keep_same_strand_only
) const {
    std::vector<MM128> mm128s;

    uint32_t read_len = static_cast<uint32_t>(read.size());

    std::vector<MiniOnQuery> qmini;
    qmini.reserve(read_len / (chainOpts_.w ? chainOpts_.w : 1) + 8);

    for_each_minimizer_(read, chainOpts_.k, chainOpts_.w, mask_, shift1_,
        [&]( __uint128_t key, uint32_t q_off, uint8_t dir ){
            qmini.push_back({key, q_off, dir});
        }
    );

    for (auto& q : qmini) {
        span_pos hit = lookup_hash(q.key);
        if (!hit) continue;
        if (hit.size() > chainOpts_.mid_occ) continue;  // Skip high-frequency minimizers

        for (pos_t p : hit) {
            const uint8_t ref_dir = static_cast<uint8_t>(is_rev(p));  // 0:FWD, 1:REV

            if (keep_same_strand_only && q.dir != ref_dir) continue;

            const uint8_t  rev   = static_cast<uint8_t>(q.dir ^ ref_dir);
            const uint32_t r_id  = seq_id(p);
            const uint32_t r_off = offset(p);
            const uint16_t fhi   = 0;
            const uint8_t  segid = 0;
            const uint8_t  qspan = static_cast<uint8_t>(chainOpts_.k);
            const uint32_t q_off = seed_q_fwd(q.q_off, rev, read_len, chainOpts_.k);

            // x: rev<<63 | tid<<32 | tpos
            // y: flags_hi<<48 | segid<<40 | qspan<<32 | q_off
            mm128s.push_back(MM128::make(rev, r_id, r_off, fhi, segid, qspan, q_off));
        }
    }

    return mm128s;
}

void MinimizerIndex::sort_seeds_by_ref(std::vector<MM128>& seeds) const {
    std::sort(seeds.begin(), seeds.end(), [](const MM128& A, const MM128& B){
        if (A.x != B.x) return A.x < B.x;
        return A.y < B.y;
    });
}
void MinimizerIndex::sort_seeds_by_qry(std::vector<MM128>& seeds) const {
    std::sort(seeds.begin(), seeds.end(), [](const MM128& A, const MM128& B){
        if (A.y != B.y) return A.y < B.y;
        return A.x < B.x;
    });
}


void MinimizerIndex::count_depth(const std::vector<std::string>& reads)
{
    log_stream() << "Counting k-mer depth from reads ...\n";

    if (ht_.empty() || keys_.empty()) { depths_.clear(); return; }

    depths_.assign(keys_.size(), 0);

    const size_t n_threads = chainOpts_.threads ? size_t(chainOpts_.threads) : 1;
    const size_t QUEUE_CAP   = 64;
    const size_t BATCH_READS = 2048;

    static constexpr size_t N_SHARDS = 1024;
    std::array<std::mutex, N_SHARDS> shard_mtx;

    struct ReadBatch {
        std::vector<std::string> seqs;
    };

    std::queue<std::unique_ptr<ReadBatch>> q;
    std::mutex mtx;
    std::condition_variable cv_prod, cv_cons;
    std::atomic<bool> producer_done{false};

    ThreadPool pool(n_threads);
    std::vector<std::future<void>> futs;
    futs.reserve(n_threads);

    auto worker = [&]() {
        for (;;) {
            std::unique_ptr<ReadBatch> batch;
            {
                std::unique_lock<std::mutex> lk(mtx);
                cv_cons.wait(lk, [&] { return producer_done.load(std::memory_order_acquire) || !q.empty(); });
                if (q.empty()) {
                    if (producer_done.load(std::memory_order_acquire)) break;
                    else continue;
                }
                batch = std::move(q.front());
                q.pop();
                lk.unlock();
                cv_prod.notify_one();
            }

            for (const auto& seq : batch->seqs) {
                for_each_kmer_(seq, chainOpts_.k, mask_, shift1_,
                    [&](__uint128_t key, uint32_t off, uint8_t dir){
                        uint64_t p = hash128_(key) & ht_mask_;
                        for (;; p = (p + 1) & ht_mask_) {
                            uint32_t t = ht_[p];
                            if (t == 0) break;
                            uint32_t idx = t - 1;
                            if (keys_[idx] == key) {
                                std::lock_guard<std::mutex> lk(shard_mtx[idx & (N_SHARDS - 1)]);
                                uint16_t& v = depths_[idx];
                                if (v != 0xFFFF) {
                                    uint32_t nv = uint32_t(v) + 1;
                                    v = (nv > 0xFFFF) ? 0xFFFF : uint16_t(nv);
                                }
                                break;
                            }
                        }
                    }
                );
            }
        }
    };

    for (size_t t = 0; t < n_threads; ++t) {
        futs.emplace_back(pool.submit(worker));
    }

    auto flush_batch = [&](std::unique_ptr<ReadBatch>& batch) {
        if (!batch || batch->seqs.empty()) return;
        std::unique_lock<std::mutex> lk(mtx);
        cv_prod.wait(lk, [&]{ return q.size() < QUEUE_CAP; });
        q.push(std::move(batch));
        lk.unlock();
        cv_cons.notify_one();
    };

    auto tracker = ProgressTracker::Every(30000);  // print every 30,000 reads
    for (const auto& path : reads) {
        gzFile fp = gzopen(path.c_str(), "r");
        if (!fp) {
            error_stream() << path << ": No such file or directory\n";
            std::exit(1);
        }
        kseq_t* ks = kseq_init(fp);

        std::unique_ptr<ReadBatch> batch = std::make_unique<ReadBatch>();
        batch->seqs.reserve(BATCH_READS);

        while (kseq_read(ks) >= 0) {
            tracker.hit();  // update progress
            batch->seqs.emplace_back(ks->seq.s, ks->seq.l);
            if (batch->seqs.size() >= BATCH_READS) {
                flush_batch(batch);
                batch = std::make_unique<ReadBatch>();
                batch->seqs.reserve(BATCH_READS);
            }
        }
        flush_batch(batch);

        kseq_destroy(ks);
        gzclose(fp);
    }
    tracker.finish();  // finish progress

    producer_done.store(true, std::memory_order_release);
    cv_cons.notify_all();

    for (auto& f : futs) f.get();
    pool.stop();
}

/* ========================================================================
 *                               CHAINING
 * ===================================================================== */
int32_t MinimizerIndex::mm_comput_sc_(
    const MM128* ai, const MM128* aj,
    int32_t max_dist_x, int32_t max_dist_y, int32_t bw,
    float chn_pen_gap, float chn_pen_skip,
    int is_cdna, int n_seg
) {
    int32_t dq = int32_t(ai->y) - int32_t(aj->y);
    int32_t sidi = ai->q_span();
    int32_t sidj = aj->q_span();
    if (dq <= 0 || dq > max_dist_x) return INT32_MIN;

    int32_t dr = int32_t(ai->x) - int32_t(aj->x);
    if (sidi == sidj && (dr == 0 || dq > max_dist_y)) return INT32_MIN;

    int32_t dd = dr > dq ? (dr - dq) : (dq - dr);
    if (sidi == sidj && dd > bw) return INT32_MIN;
    if (n_seg > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) return INT32_MIN;

    int32_t dg = (dr < dq) ? dr : dq;
    int32_t q_span = aj->q_span();
    int32_t sc = (q_span < dg ? q_span : dg);

    if (dd || dg > q_span) {
        float lin_pen = chn_pen_gap * float(dd) + chn_pen_skip * float(dg);
        float log_pen = (dd >= 1) ? float(mm_log2_i32_(dd + 1)) : 0.0f;
        if (is_cdna || sidi != sidj) {
            if (sidi != sidj && dr == 0) ++sc;
            else if (dr > dq || sidi != sidj) sc -= int(std::min(lin_pen, log_pen));
            else sc -= int(lin_pen + .5f * log_pen);
        } else {
            sc -= int(lin_pen + .5f * log_pen);
        }
    }
    return sc;
}

int64_t MinimizerIndex::mm_chain_bk_end_(
    int32_t max_drop,
    const std::vector<int32_t>& f,
    const std::vector<int64_t>& p,
    std::vector<int32_t>& t,
    const std::vector<std::pair<int32_t,int64_t>>& z, // z[k] = {f[i], i}
    int64_t k
) {
    int64_t i = z[k].second, end_i = -1, max_i = i;
    int32_t max_s = 0;
    if (i < 0 || t[i] != 0) return i;
    do {
        t[i] = 2;
        end_i = i = p[i];
        int32_t s = (i < 0) ? z[k].first : (int32_t)z[k].first - f[i];
        if (s > max_s) { max_s = s; max_i = i; }
        else if (max_s - s > max_drop) break;
    } while (i >= 0 && t[i] == 0);
    for (i = z[k].second; i >= 0 && i != end_i; i = p[i]) t[i] = 0;
    return max_i;
}

std::vector<Sc_Len> MinimizerIndex::mm_chain_backtrack_ (
    int64_t n,
    const std::vector<int32_t>& f, const std::vector<int64_t>& p,
    std::vector<int32_t>& v, std::vector<int32_t>& t,
    int32_t min_cnt, int32_t min_sc, int32_t max_drop,
    int32_t& n_u, int32_t& n_v
) {
    std::vector<std::pair<int32_t,int64_t>> z; z.reserve(n);
    for (int64_t i=0;i<n;++i) if (f[i] >= min_sc) z.push_back({f[i], i});
    if (z.empty()) { n_u = n_v = 0; return {}; }

    std::sort(z.begin(), z.end(), [](auto& A, auto& B){ return (uint64_t)A.first < (uint64_t)B.first; });

    std::fill(t.begin(), t.end(), 0);
    int64_t k; n_u = n_v = 0;

    for (k = (int64_t)z.size() - 1; k >= 0; --k) {
        if (t[z[k].second] == 0) {
            int64_t n_v0 = n_v, end_i;
            int32_t sc;
            end_i = mm_chain_bk_end_(max_drop, f, p, t, z, k);
            for (int64_t i = z[k].second; i != end_i; i = p[i]) ++n_v, t[i] = 1;
            sc = (end_i < 0) ? z[k].first : (int32_t)z[k].first - f[end_i];
            if (sc >= min_sc && n_v > n_v0 && (n_v - n_v0) >= min_cnt) ++n_u;
            else n_v = n_v0;
        }
    }

    std::vector<Sc_Len> u(n_u);
    std::fill(t.begin(), t.end(), 0);
    int32_t ucnt = 0; n_v = 0;
    for (k = (int64_t)z.size() - 1; k >= 0; --k) {
        if (t[z[k].second] == 0) {
            int64_t n_v0 = n_v, end_i;
            int32_t sc;
            end_i = mm_chain_bk_end_(max_drop, f, p, t, z, k);
            for (int64_t i = z[k].second; i != end_i; i = p[i])
                v[n_v++] = (int32_t)i, t[i] = 1;
            sc = (end_i < 0) ? z[k].first : (int32_t)z[k].first - f[end_i];
            if (sc >= min_sc && n_v > n_v0 && (n_v - n_v0) >= min_cnt)
                u[ucnt++].sc_len = (uint64_t(uint32_t(sc)) << 32) | uint32_t(n_v - n_v0);
            else n_v = n_v0;
        }
    }
    n_u = ucnt;
    return u;
}

void MinimizerIndex::mm_compact_sort_ (
    const std::vector<MM128>& a,
    const std::vector<int32_t>& v,
    const std::vector<Sc_Len>& u,
    std::vector<MM128>& b_out,
    std::vector<Sc_Len>& u_sorted
) {
    b_out.clear(); u_sorted.clear();

    std::vector<MM128> b_tmp; b_tmp.reserve(v.size());
    {
        int64_t k = 0;
        for (size_t i = 0; i < u.size(); ++i) {
            int32_t len = u[i].len();
            int64_t k0 = k;
            for (int32_t j = 0; j < len; ++j)
                b_tmp.push_back(a[ v[k0 + (len - j - 1)] ]);
            k += len;
        }
    }

    struct Key { uint64_t first_x; int32_t chain_idx; int64_t offset; int32_t len; };
    std::vector<Key> keys; keys.reserve(u.size());
    {
        int64_t off = 0;
        for (size_t i = 0; i < u.size(); ++i) {
            int32_t len = u[i].len();
            keys.push_back({ b_tmp[off].x, (int32_t)i, off, len });
            off += len;
        }
    }
    std::sort(keys.begin(), keys.end(), [](const Key& A, const Key& B) { return A.first_x < B.first_x; });

    b_out.reserve(b_tmp.size());
    u_sorted.resize(u.size());
    for (size_t i = 0; i < keys.size(); ++i) {
        const auto& K = keys[i];
        u_sorted[i] = u[K.chain_idx];
        for (int32_t j = 0; j < K.len; ++j)
            b_out.push_back(b_tmp[K.offset + j]);
    }
}

void MinimizerIndex::chain_dp(
    const std::vector<MM128>& seeds,        // input anchors (sorted by x,y)
    std::vector<MM128>& b_out,              // output anchors (compacted & sorted per chain)
    std::vector<Sc_Len>& u_out              // output chains (score<<32 | len), sorted
) const {
    const int64_t n = (int64_t)seeds.size();
    b_out.clear(); u_out.clear();
    if (n == 0) return;

    if (chainOpts_.max_dist_x < chainOpts_.bw) chainOpts_.max_dist_x = chainOpts_.bw;
    if (chainOpts_.max_dist_y < chainOpts_.bw && !chainOpts_.is_cdna) chainOpts_.max_dist_y = chainOpts_.bw;
    int32_t max_drop = chainOpts_.is_cdna ? INT32_MAX : chainOpts_.bw;

    // DP arrays
    std::vector<int32_t> f(n), v(n), t(n, 0);
    std::vector<int64_t> p(n);

    // DP loop
    int32_t mmax_f = 0;
    int64_t st = 0, max_ii = -1;
    for (int64_t i = 0; i < n; ++i) {
        int64_t max_j = -1;
        int32_t q_span = seeds[i].q_span();
        int32_t max_f  = q_span;
        int32_t n_skip = 0;

        while (st < i && ( seeds[i].r_id() != seeds[st].r_id() || (int64_t)seeds[i].x > (int64_t)seeds[st].x + chainOpts_.max_dist_x )) ++st;
        if (i - st > chainOpts_.max_iter) st = i - chainOpts_.max_iter;

        int64_t j;
        for (j = i - 1; j >= st; --j) {
            int32_t sc = mm_comput_sc_(
                &seeds[i], &seeds[j],
                chainOpts_.max_dist_x, chainOpts_.max_dist_y, chainOpts_.bw,
                chainOpts_.chn_pen_gap, chainOpts_.chn_pen_skip,
                chainOpts_.is_cdna, chainOpts_.n_seg
            );
            if (sc == INT32_MIN) continue;
            sc += f[j];
            if (sc > max_f) {
                max_f = sc; max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == (int32_t)i) {
                if (++n_skip > chainOpts_.max_skip) break;
            }
            if (p[j] >= 0) t[p[j]] = (int32_t)i;
        }
        int64_t end_j = j;

        if (max_ii < 0 || (int64_t)seeds[i].x - (int64_t)seeds[max_ii].x > (int64_t)chainOpts_.max_dist_x) {
            int32_t best = INT32_MIN; max_ii = -1;
            for (j = i - 1; j >= st; --j)
                if (best < f[j]) { best = f[j]; max_ii = j; }
        }
        if (max_ii >= 0 && max_ii < end_j) {
            int32_t tmp = mm_comput_sc_(
                &seeds[i], &seeds[max_ii],
                chainOpts_.max_dist_x, chainOpts_.max_dist_y, chainOpts_.bw,
                chainOpts_.chn_pen_gap, chainOpts_.chn_pen_skip,
                chainOpts_.is_cdna, chainOpts_.n_seg
            );
            if (tmp != INT32_MIN && max_f < tmp + f[max_ii]) {
                max_f = tmp + f[max_ii]; max_j = max_ii;
            }
        }
        f[i] = max_f; p[i] = max_j;
        v[i] = (max_j >= 0 && v[max_j] > max_f) ? v[max_j] : max_f;
        if (max_ii < 0 || ( (int64_t)seeds[i].x - (int64_t)seeds[max_ii].x <= (int64_t)chainOpts_.max_dist_x && f[max_ii] < f[i]) ) max_ii = i;
        if (mmax_f < max_f) mmax_f = max_f;
    }

    // Backtrack to collect chains (u: meta, v: indices)
    int32_t n_u = 0, n_v = 0;
    std::vector<Sc_Len> u = mm_chain_backtrack_(n, f, p, v, t, chainOpts_.min_cnt, chainOpts_.min_sc, max_drop, n_u, n_v);
    if (n_u == 0) return; // no chains

    // Compact and sort per-chain
    mm_compact_sort_(seeds, v, u, b_out, u_out);
}

/* --------------------------------------------------------------------
 * build anchor chains from seed chains
* ------------------------------------------------------------------ */
std::vector<AnchorChain> MinimizerIndex::build_anchor_from_chains(
    const std::vector<MM128>& b_out,              // output anchors (compacted & sorted per chain)
    const std::vector<Sc_Len>& u_out,             // output chains (score<<32 | len), sorted
    const uint32_t& read_len
) const {
    std::vector<AnchorChain> anchors;
    anchors.resize(u_out.size());
    uint32_t ci = 0;
    uint64_t offset = 0;
    for (const auto& sc_len : u_out) {
        const int32_t score = sc_len.score();
        const uint32_t len = sc_len.len();

        auto it_b = b_out.begin() + offset;
        auto it_e = it_b + len;

        auto& ac = anchors[ci];
        ac.r_id = it_b->r_id();
        ac.chain_rev = it_b->dir() == 1;
        ac.seed_number = len;
        ac.score = score;
        auto& blocks = ac.blocks;
        blocks.reserve(len);

        for (auto it = it_b; it != it_e; ++it) {
            uint32_t r_beg = it->r_off(), r_end = r_beg + chainOpts_.k;
            uint32_t q_beg = it->q_off();
            uint32_t q_end = q_beg + chainOpts_.k;
            if (blocks.empty()) {
                blocks.push_back({r_beg, r_end, q_beg, q_end});
            } else {
                auto& cur = blocks.back();
                int64_t rgap = int64_t(r_beg) - int64_t(cur.r_end);
                int64_t qgap = int64_t(q_beg) - int64_t(cur.q_end);
                if (rgap < 0 || qgap < 0) {  // negative gap, skip
                    continue;
                } else if (rgap == 0 && qgap == 0) {
                    cur.r_end = std::max(cur.r_end, r_end);
                    cur.q_end = std::max(cur.q_end, q_end);
                } else {
                    blocks.push_back({r_beg, r_end, q_beg, q_end});
                }
            }
        }
        ci++;
        offset += len;
    }
    return anchors;
}

/* --------------------------------------------------------------------
 * Cross-chain merging:
 *   Sort chains by (r_id, chain_rev, first_block.r_beg) and merge
 *   adjacent if they belong to same contig/strand and gaps satisfy
 *   the same small_gap/small_slop criteria used inside one chain.
* ------------------------------------------------------------------ */
void MinimizerIndex::anchors_coordinate_merge_(
    std::vector<Anchor>& anchors
) const {
    if (anchors.size() <= 1) return;

    // ---- sort by (r_beg, q_beg) --------------------
    std::sort(anchors.begin(), anchors.end(), [&](const Anchor& A, const Anchor& B) {
        if (A.r_beg != B.r_beg) return A.r_beg < B.r_beg;
        return A.q_beg < B.q_beg;
    });

    // ---- merge anchor chains ---------------------------------------
    size_t i = 0;
    while (i < anchors.size() - 1) {
        auto& last = anchors[i];
        auto& cur = anchors[i + 1];

        if (last.r_end >= cur.r_beg && last.q_end >= cur.q_beg) {  // overlapping, merge
            // last_r: [------->    ]         last_q: [------->    ]
            // cur_r : [     ------>]         cur_q : [     ------>]
            last.r_beg = std::min(last.r_beg, cur.r_beg);
            last.r_end = std::max(last.r_end, cur.r_end);
            last.q_beg = std::min(last.q_beg, cur.q_beg);
            last.q_end = std::max(last.q_end, cur.q_end);
            anchors.erase(anchors.begin() + i + 1);
        } else if (last.r_end < cur.r_beg && last.q_end < cur.q_beg) {  // non-overlapping, move to next anchor
            // last_r: [------->          ]         last_q: [------->          ]
            // cur_r : [           ------>]         cur_q : [           ------>]
            ++i; // move to next anchor
        } else {  // one of the anchor is overlapping, and the other is not, juest delete the next anchor (07/27/2025)
            // last_r: [------->     ]         last_q: [------->          ]
            // cur_r : [      ------>]         cur_q : [           ------>]
            anchors.erase(anchors.begin() + i + 1);
        }
    }
}

std::vector<AnchorChain> MinimizerIndex::merge_anchor_chains(
    const std::vector<AnchorChain>& anchors
) const {
    if (anchors.size() <= 1) return anchors;

    // ---- helpers -------------------------------------------------
    auto gap = [](uint32_t last_end, uint32_t cur_beg) -> int64_t { // positive: forward gap; negative: overlap
        return (int64_t)cur_beg - (int64_t)last_end;
    };
    auto smallGap2D = [](const int64_t rgap, const int64_t qgap, const opt::AnchorOpts& params) -> bool {
        return (rgap >= 0 && qgap >= 0 && rgap <= params.small_gap && qgap <= params.small_gap);
    };
    auto smallSlop2D = [](const int64_t rgap, const int64_t qgap, const opt::AnchorOpts& params) -> bool {
        return (rgap >= 0 && qgap >= 0 && std::abs(rgap - qgap) <= params.small_slop);
    };
    auto increase2D = [](const Anchor& last_blk, const Anchor& cur_blk) -> bool {
        return last_blk.r_end < cur_blk.r_beg  && last_blk.q_end < cur_blk.q_beg;
    };
    auto contains2D = [](const Anchor& last_blk, const Anchor& cur_blk, const int64_t& rgap, const int64_t& qgap) -> bool {
        return rgap == qgap && cur_blk.r_beg <= last_blk.r_end && cur_blk.r_beg >= last_blk.r_beg && cur_blk.q_beg <= last_blk.q_end && cur_blk.q_beg >= last_blk.q_beg;
    };
    auto equal2D = [](const Anchor& last_blk, const Anchor& cur_blk) -> bool {
        return last_blk.r_beg == cur_blk.r_beg && last_blk.r_end == cur_blk.r_end && last_blk.q_beg == cur_blk.q_beg && last_blk.q_end == cur_blk.q_end;
    };

    // ---- sort by (r_id, strand, first r_beg) --------------------
    std::vector<size_t> order(anchors.size());
    std::iota(order.begin(), order.end(), size_t(0));
    std::sort(order.begin(), order.end(), [&](size_t A, size_t B) {
        const auto& a = anchors[A]; const auto& b = anchors[B];
        if (a.r_id != b.r_id) return a.r_id < b.r_id;
        if (a.chain_rev != b.chain_rev) return a.chain_rev < b.chain_rev;
        return a.blocks.front().r_beg < b.blocks.front().r_beg;
    });

    std::vector<AnchorChain> merged;
    merged.reserve(anchors.size());

    // ---- merge anchor chains ------------------------------------
    for (size_t oi = 0; oi < order.size(); ++oi) {
        const AnchorChain& cur = anchors[order[oi]];
        if (merged.empty()) {
            merged.push_back(cur);
            continue;
        }
        AnchorChain& last = merged.back();
        
        if (last.r_id == cur.r_id && last.chain_rev == cur.chain_rev) {
            // try to append each block from cur into last
            auto& last_blk = last.blocks.back();
            auto& cur_blk = cur.blocks.front();
            int64_t rgap = gap(last_blk.r_end, cur_blk.r_beg);
            int64_t qgap = gap(last_blk.q_end, cur_blk.q_beg);

            bool small_gap_bool = smallGap2D(rgap, qgap, anchorOpts_);
            bool small_slop_bool = smallSlop2D(rgap, qgap, anchorOpts_);
            bool increasing_bool = increase2D(last_blk, cur_blk);

            if (small_gap_bool && small_slop_bool && increasing_bool) {
                // Case 1: Small-gap append
                // last: [--->     ]
                // cur : [     --->]
                // rgap == qgap > 0 and <= small_gap
                last.blocks.insert(last.blocks.end(), cur.blocks.begin(), cur.blocks.end());
                last.seed_number += cur.seed_number;
                last.score       += cur.score;
            } else if (contains2D(last_blk, cur_blk, rgap, qgap)) {
                // Case 2: Coordinate merge (overlap or touch)
                // last: [------->    ]
                // cur : [     ------>]
                last.blocks.insert(last.blocks.end(), cur.blocks.begin(), cur.blocks.end());
                anchors_coordinate_merge_(last.blocks);
                last.seed_number += cur.seed_number;
                last.score       += cur.score;
            } else {
                // Case 3: Too far apart — start a new chain
                // last: [--->        ]
                // cur : [        --->]
                merged.push_back(cur);
            }
        } else {
            // Case 4: Different contig or orientation — always split
            // last: rev=0 [--->]
            // cur : rev=1 [<---]
            merged.push_back(cur);
        }
    }

    return merged;
}

std::vector<std::vector<uint32_t>> MinimizerIndex::group_anchor_chains_by_q_overlap(
    const std::vector<AnchorChain>& anchors
) const {
    std::vector<std::vector<uint32_t>> groups;  // output groups of anchor indexes, such as: [[1,4,5], [2,3], ...].
    if (anchors.empty()) return groups;

    const float MIN_OVLP = anchorOpts_.min_group_overlap_frac;

    struct Span {
        uint32_t idx;    // index in anchors
        uint32_t q_beg;
        uint32_t q_end;
        uint32_t len;
        int      score = 0;
    };

    std::vector<Span> spans;
    spans.reserve(anchors.size());

    for (uint32_t i = 0; i < anchors.size(); ++i) {
        const auto& ac = anchors[i];
        if (ac.blocks.empty()) continue;

        uint32_t q_min = UINT32_MAX; uint32_t q_max = 0;

        for (const auto& b : ac.blocks) {
            q_min = std::min({q_min, b.q_beg, b.q_end});
            q_max = std::max({q_max, b.q_beg, b.q_end});
        }

        if (q_max > q_min) {
            spans.push_back(Span{ i, q_min, q_max, q_max - q_min, ac.score });
        }
    }
    if (spans.empty()) return groups;

    std::sort(spans.begin(), spans.end(),
        [](const Span& A, const Span& B) {
            if (A.score != B.score) return A.score > B.score;
            if (A.len != B.len)     return A.len > B.len;
            if (A.q_beg != B.q_beg) return A.q_beg < B.q_beg;
            return A.q_end < B.q_end;
        }
    );

    std::vector<uint8_t> used(anchors.size(), 0);

    for (size_t si = 0; si < spans.size(); ++si) {
        const Span& seed = spans[si];
        uint32_t seed_idx = seed.idx;
        if (used[seed_idx]) continue;

        std::vector<uint32_t> group;
        group.reserve(8);
        group.push_back(seed_idx);
        used[seed_idx] = 1;

        uint32_t sb = seed.q_beg;
        uint32_t se = seed.q_end;
        uint32_t sl = seed.len;

        for (size_t sj = si + 1; sj < spans.size(); ++sj) {
            const Span& sp = spans[sj];
            uint32_t idx = sp.idx;
            if (used[idx]) continue;

            uint32_t ovlp_beg = std::max(sb, sp.q_beg);
            uint32_t ovlp_end = std::min(se, sp.q_end);
            if (ovlp_end <= ovlp_beg) continue;

            uint32_t ovlp_len = ovlp_end - ovlp_beg;
            uint32_t short_len = std::min(sl, sp.len);
            if (short_len == 0) continue;

            float frac = static_cast<float>(ovlp_len) / static_cast<float>(short_len);

            bool same_group;
            if (MIN_OVLP > 0.0f)
                same_group = (frac >= MIN_OVLP);
            else
                same_group = (ovlp_len > 0);

            if (!same_group) continue;

            group.push_back(idx);
            used[idx] = 1;
        }

        std::sort(group.begin(), group.end(),
            [&](uint32_t ia, uint32_t ib) {
                const auto& A = anchors[ia];
                const auto& B = anchors[ib];
                if (A.score != B.score) return A.score > B.score;
                return A.seed_number > B.seed_number;
            }
        );

        groups.push_back(std::move(group));
    }

    return groups;
}

std::vector<std::vector<uint32_t>> MinimizerIndex::filter_anchor_chains(
    std::vector<std::vector<uint32_t>> groups, 
    int sec_pri_num
) const {
    std::vector<std::vector<uint32_t>> filtered;  // output groups of anchor indexes, such as: [[1,4,5], [2,3], ...].

    if (groups.empty()) return filtered;

    filtered.reserve(groups.size());

    for (auto& group : groups) {
        if (group.size() <= static_cast<size_t>(sec_pri_num)) {
            filtered.push_back(std::move(group));
        } else {
            std::vector<uint32_t> sub_group;
            sub_group.reserve(sec_pri_num);
            for (int i = 0; i < sec_pri_num; ++i) {
                sub_group.push_back(group[i]);
            }
            filtered.push_back(std::move(sub_group));
        }
    }

    uint16_t max_kept = anchorOpts_.max_kept;

    if (max_kept > 0 && filtered.size() > max_kept) {
        filtered.resize(max_kept);
    }

    return filtered;
}

/* ------------------------------------ Debug print ------------------------------------ */
void MinimizerIndex::print_seeds(
    const std::string_view read_name,
    uint32_t read_len, 
    const std::vector<MM128>& seeds
) const {
    debug_stream() << "-- Seeds (" << seeds.size() << ") --\n";
    debug_stream() << "Contig\tRead\tr_off\tq_off\tdir\n";
    debug_stream() << std::string(50, '-') << "\n";

    for (auto const &seed : seeds) {
        uint8_t dir   = uint8_t(seed.dir());
        uint32_t r_id = uint32_t(seed.r_id());
        uint32_t r_off= uint32_t(seed.r_off());
        uint32_t q_off= uint32_t(seed.q_off());
        const auto &contig_name = names_[r_id];
        char strand = dir ? '-' : '+';
        debug_stream() << contig_name << "\t" << read_name << "\t" << r_off << "\t" << q_off << "\t" << strand << "\n";
    }
    debug_stream() << "\n\n";
}

void MinimizerIndex::print_chains(
    const std::string_view read_name,
    uint32_t read_len,
    const std::vector<MM128>& b_out,
    const std::vector<Sc_Len>& u_out
) const {
    using std::setw;
    using std::left;

    const int posW   = 12;

    debug_stream() << "-- Chain seeds (" << u_out.size() << ") --\n";

    uint64_t offset = 0;
    for (size_t ci = 0; ci < u_out.size(); ++ci) {
        const int32_t  score = u_out[ci].score();
        const uint32_t len   = u_out[ci].len();

        if (len == 0) continue;

        auto it_b = b_out.begin() + offset;
        auto it_e = it_b + len;
        offset   += len;

        const auto  contig_name = names_[it_b->r_id()];
        const char  strand      = it_b->dir() ? '-' : '+';

        debug_stream() << "Chain " << ci << " (" << contig_name << " vs " << read_name << ")\n";
        debug_stream() << "   - score: "  << score  << '\n';
        debug_stream() << "   - strand: " << strand << '\n';
        debug_stream() << "   - seeds: "  << len    << '\n';

        debug_stream()
            << "      "
            << setw(posW)  << left << "r_off"
            << setw(posW)  << left << "q_off"
            << setw(posW)  << left << "(r_gap, q_gap, gap_diff)"
            << '\n';

        uint32_t last_r_off = 0;
        uint32_t last_q_off = 0;
        for (auto it = it_b; it != it_e; ++it) {
            const uint32_t r_off  = it->r_off();
            const uint32_t q_off  = it->q_off();
            int64_t        r_gap  = (it != it_b) ? static_cast<int64_t>(r_off) - static_cast<int64_t>(last_r_off) : 0;
            int64_t        q_gap  = (it != it_b) ? static_cast<int64_t>(q_off) - static_cast<int64_t>(last_q_off) : 0;
            int64_t        gap_diff = r_gap - q_gap;
            last_r_off = r_off;
            last_q_off = q_off;

            debug_stream()
                << "      "
                << setw(posW)  << left << r_off
                << setw(posW)  << left << q_off
                << "(" << r_gap << "," << q_gap << "," << gap_diff << ")"
                << '\n';
        }

        debug_stream() << '\n';
    }

    debug_stream() << "\n\n";
}

void MinimizerIndex::print_anchors(
    const std::string_view read_name,
    uint32_t read_len, 
    const std::vector<AnchorChain>& anchorChains, 
    const std::vector<std::vector<uint32_t>>& groups,
    std::string type
) const {
    if (!type.empty()) debug_stream() << "-- Anchor Blocks (" << type << ") --\n";
    else debug_stream() << "-- Anchor Blocks (" << anchorChains.size() << ") --\n";

    for (uint32_t gi = 0; gi < groups.size(); ++gi) {
        const auto& group = groups[gi];
        debug_stream() << "Anchor Group " << gi << " (size: " << group.size() << "):\n";
        for (size_t i = 0; i < group.size(); ++i) {
            uint32_t ai = group[i];
            const auto &ac = anchorChains[ai];
            const auto &contig_name = names_[ac.r_id];
            debug_stream() << "  Anchor " << ai << " (" << read_name << " and " << contig_name << ")\n";
            uint32_t r_beg = UINT32_MAX; uint32_t r_end = 0;
            uint32_t q_beg = UINT32_MAX; uint32_t q_end = 0;
            for (auto const &b : ac.blocks) {
                r_beg = std::min({r_beg, b.r_beg, b.r_end});
                r_end = std::max({r_end, b.r_beg, b.r_end});
                q_beg = std::min({q_beg, b.q_beg, b.q_end});
                q_end = std::max({q_end, b.q_beg, b.q_end});
            }

            std::string score_info  = "score: " + std::to_string(ac.score);
            std::string strand_info = "strand: " + std::string(1, ac.chain_rev ? '-' : '+');
            std::string seed_info   = "seeds: " + std::to_string(ac.seed_number);
            std::string length_info = "length: " + std::to_string(r_end > r_beg ? r_end - r_beg : r_beg - r_end);
            std::string ref_range   = "ref: [" + std::to_string(r_beg) + ", " + std::to_string(r_end) + ")";
            std::string qry_range   = "qry: [" + std::to_string(q_beg) + ", " + std::to_string(q_end) + ")";

            debug_stream() << "   - " << score_info << std::endl;
            debug_stream() << "   - " << strand_info << std::endl;
            debug_stream() << "   - " << seed_info << std::endl;
            debug_stream() << "   - " << length_info << std::endl;
            debug_stream() << "   - " << ref_range << std::endl;
            debug_stream() << "   - " << qry_range << std::endl << std::endl;
        }
        debug_stream() << "\n";
    }

    debug_stream() << "\n\n";
}

void MinimizerIndex::print_index_stats() const {
    const uint64_t n_seqs = static_cast<uint64_t>(names_.size());
    uint64_t total_len = 0;
    for (auto const& seq : seqs_) total_len += static_cast<uint64_t>(seq.size());

    const uint64_t distinct = static_cast<uint64_t>(keys_.size());
    const uint64_t total_minimizers = static_cast<uint64_t>(positions_.size());

    // Count singletons and average occurrences
    uint64_t singletons = 0;
    for (size_t i = 0; i < keys_.size(); ++i) {
        uint32_t c = offs_[i + 1] - offs_[i];
        if (c == 1) ++singletons;
    }

    const double pct_singleton = (distinct ? (100.0 * double(singletons) / double(distinct)) : 0.0);
    const double avg_occ       = (distinct ? (double(total_minimizers)   / double(distinct)) : 0.0);
    const double avg_spacing   = (total_minimizers ? (double(total_len)  / double(total_minimizers)) : 0.0);

    const int label_width = 25, value_width = 12;

    // Minimizer Index stats
    log_stream() << std::left << std::setw(label_width) << "Minimizer Index stats:" << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "k-mer size:" << std::right << std::setw(value_width) << chainOpts_.k << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Window size:" << std::right << std::setw(value_width) << chainOpts_.w << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Sequences number:" << std::right << std::setw(value_width) << n_seqs << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Distinct minimizers:" << std::right << std::setw(value_width) << distinct << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Singletons (%):" << std::right << std::setw(value_width) << std::fixed << std::setprecision(2) << pct_singleton << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Avg occurrences:" << std::right << std::setw(value_width) << std::fixed << std::setprecision(3) << avg_occ << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Avg spacing:" << std::right << std::setw(value_width) << avg_spacing << '\n';
    log_stream() << "   - " << std::left << std::setw(label_width) << "Total length:" << std::right << std::setw(value_width) << total_len << '\n';

    if (chainOpts_.mid_occ_frac > 0.0f) {
        log_stream() << "   - " << std::left << std::setw(label_width) << "Drop fraction:" << std::right << std::setw(value_width) << std::setprecision(4) << chainOpts_.mid_occ_frac << '\n';
        log_stream() << "   - " << std::left << std::setw(label_width) << "Max occ threshold:" << std::right << std::setw(value_width) << chainOpts_.mid_occ << '\n';
    }
}






inline bool MinimizerIndex::is_connected_vertex_(uint32_t v_from, uint32_t v_to) const {
    if (!graph_) return false;
    return graph_->is_connected_vertex(v_from, v_to, /*step_cap=*/200000);
}

// ============= Cross-segment transition score A->B =============
// Notes: dq uses read gap + previous segment length (A.r_end-A.r_beg);
//        dr uses shortest added bases on the graph (ow-aware) from (A.r_id, A.r_end) to (B.r_id, B.r_beg);
//        scoring follows mm2 style, treating (sidi != sidj) as true (cross-segment).
int32_t MinimizerIndex::score_transition_seg_(const SegNode& A, const SegNode& B) const {
    std::cerr << "[D::score_transition_seg_] A: (r_id=" << A.r_id << ", rev=" << A.rev << ", r_beg=" << A.r_beg << ", r_end=" << A.r_end << ", q_beg=" << A.q_beg << ", q_end=" << A.q_end << ")\n";
    std::cerr << "[D::score_transition_seg_] B: (r_id=" << B.r_id << ", rev=" << B.rev << ", r_beg=" << B.r_beg << ", r_end=" << B.r_end << ", q_beg=" << B.q_beg << ", q_end=" << B.q_end << ")\n";

    // --- Fetch segment lengths (needed to convert forward coords to oriented offsets) ---
    auto get_len = [&](uint32_t seg)->uint32_t{
        const GfaNode* n = graph_->getNode(seg);
        return (n && !n->deleted) ? n->length : 0u;
    };
    const uint32_t lenA = get_len(A.r_id);
    const uint32_t lenB = get_len(B.r_id);

    // Base read gap
    int64_t dq = (int64_t)B.q_beg - (int64_t)A.q_end;
    if (dq <= 0 || dq > chainOpts_.max_dist_x) return INT32_MIN;
    int32_t dq32 = (int32_t)dq;

    // Compute dr: shortest added bases on the graph, from A end to B start
    // Assume forward vertex (GfaGraph adds complements); for strictness use A/B direction.
    if (!graph_) return INT32_MIN;
    const uint32_t ori_A   = A.rev ? 1u : 0u;
    const uint32_t v_from  = (A.r_id << 1) | ori_A;
    const uint32_t ori_B   = B.rev ? 1u : 0u;
    const uint32_t v_to    = (B.r_id << 1) | ori_B;

    
    // --- Map forward coordinates to *oriented* offsets ---
    // Forward strand (rev==false):
    //   off_from = A.r_end                      (continue after the covered tail)
    //   off_to   = B.r_beg
    // Reverse strand (rev==true), interval [r_beg, r_end) -> oriented [len-r_end, len-r_beg):
    //   off_from = lenA - A.r_beg               (continue after covered tail in reverse orientation)
    //   off_to   = lenB - B.r_end               (stop right at the first used base on B)
    uint32_t off_from = A.rev ? (lenA > A.r_beg ? (lenA - A.r_beg) : 0u) : A.r_end;
    uint32_t off_to   = A.rev ? (lenB > B.r_end ? (lenB - B.r_end) : 0u) : B.r_beg;

    // --- Optional fast pruning using oriented connectivity ---
    if (!graph_->is_connected_vertex(v_from, v_to, /*step_cap=*/200000)) return INT32_MIN;

    // --- Shortest added-bases distance on the graph (ow-rule aware, oriented) ---
    uint64_t dist = 0;
    bool ok = graph_->shortest_distance_between_offsets(
        v_from, off_from,
        v_to,   off_to,
        dist, /*step_cap=*/200000
    );
    std::cerr << "[D::score_transition_seg_] from (v=" << v_from << ", off=" << off_from << ") to (v=" << v_to << ", off=" << off_to << "): dist=" << dist << "\n";
    if (!ok) return INT32_MIN;
    if (dist > (uint64_t)INT32_MAX) dist = (uint64_t)INT32_MAX;
    int32_t dr = (int32_t)dist;

    // pair scoring (treat A/B as different segments: sidi != sidj)
    // Reuse the formula idea, without bandwidth/same-seg limits
    int32_t dd = (dr > dq32) ? (dr - dq32) : (dq32 - dr);
    int32_t dg = (dr < dq32) ? dr : dq32;

    // Representative q_span (use k for cross-segment)
    int32_t q_span = (int32_t)chainOpts_.k;
    int32_t sc = (q_span < dg ? q_span : dg);

    if (dd || dg > q_span) {
        float lin_pen = chainOpts_.chn_pen_gap * (float)dd + chainOpts_.chn_pen_skip * (float)dg;
        float log_pen = (dd >= 1) ? (float)mm_log2_i32_(dd + 1) : 0.0f;
        // Cross-segment: sidi != sidj holds, use the looser branch
        if (dr == 0) ++sc;
        else sc -= (int)std::min(lin_pen, log_pen);
    }
    return sc;
}

// ============= Second-stage cross-segment chaining: DP in read order =============
std::vector<AnchorChain>
MinimizerIndex::cross_segment_chaining(const std::vector<AnchorChain>& per_seg_chains) const {
    std::vector<AnchorChain> empty_ret;
    if (per_seg_chains.empty() || !graph_) return empty_ret;

    // 1) Compress AnchorChain to SegNode, sort by read q_beg
    std::vector<SegNode> nodes; nodes.reserve(per_seg_chains.size());
    for (uint32_t i = 0; i < (uint32_t)per_seg_chains.size(); ++i) {
        const auto& ac = per_seg_chains[i];
        if (ac.blocks.empty()) continue;

        const auto& first = ac.blocks.front();
        const auto& last  = ac.blocks.back();

        SegNode s;
        s.r_id     = ac.r_id;
        s.rev      = ac.chain_rev;
        s.r_beg    = first.r_beg;
        s.r_end    = last.r_end;
        s.q_beg    = first.q_beg;
        s.q_end    = last.q_end;
        s.intra_sc = ac.score;
        s.seeds    = ac.seed_number;
        s.src_idx  = i;
        nodes.push_back(s);
    }
    if (nodes.empty()) return empty_ret;

    std::sort(nodes.begin(), nodes.end(), [](const SegNode& A, const SegNode& B){
        if (A.rev != B.rev) return A.rev < B.rev;       // group by direction
        if (A.q_beg != B.q_beg) return A.q_beg < B.q_beg; // then by read start
        return A.q_end < B.q_end;
    });

    // 2) Run DP per direction (no mixing)
    auto solve_one_dir = [&](bool rev_flag) -> AnchorChain {
        // Collect node indices with same direction
        std::vector<uint32_t> idx;
        for (uint32_t i=0;i<nodes.size();++i) if (nodes[i].rev == rev_flag) idx.push_back(i);
        if (idx.empty()) return AnchorChain{}; // empty

        const int64_t m = (int64_t)idx.size();
        std::vector<int32_t> F(m);      // Best score
        std::vector<int64_t> P(m, -1);  // Best predecessor index in idx; -1 for no predecessor
        int32_t bestF = INT32_MIN; int64_t bestI = -1;

        // DP main loop (window pruning; simple q window)
        for (int64_t ii = 0; ii < m; ++ii) {
            const SegNode& cur = nodes[idx[ii]];
            int32_t best_here = cur.intra_sc; // at least intra-seg score
            int64_t  prev     = -1;

            // Simple window: limit lookback count/read distance (avoid O(N^2))
            int64_t jj = ii - 1;
            int64_t looked = 0;
            while (jj >= 0 && looked < chainOpts_.max_iter) {
                const SegNode& pre = nodes[idx[jj]];
                // Read is monotonic (reverse already grouped by rev)
                // Limit read gap to avoid too far
                if ((int64_t)cur.q_beg - (int64_t)pre.q_end > (int64_t)chainOpts_.max_dist_x * 2) {
                    // Too far, break (q-sorted)
                    break;
                }

                int32_t trans = score_transition_seg_(pre, cur);
                if (trans != INT32_MIN) {
                    int32_t cand = F[jj] + trans + cur.intra_sc;
                    if (cand > best_here) { best_here = cand; prev = jj; }
                }
                --jj; ++looked;
            }
            F[ii] = best_here; P[ii] = prev;
            if (best_here > bestF) { bestF = best_here; bestI = ii; }
        }

        // 3) Backtrack best path (empty => return empty AnchorChain)
        if (bestI < 0) return AnchorChain{};

        std::vector<uint32_t> path_ids; path_ids.reserve(64);
        for (int64_t k = bestI; k >= 0; k = P[k]) path_ids.push_back(idx[k]);
        std::reverse(path_ids.begin(), path_ids.end());

        // 4) Assemble cross-seg AnchorChain: append blocks and merge coords
        AnchorChain out;
        out.r_id       = nodes[path_ids.front()].r_id;   // Note: r_id is not meaningful here; keep first seg id
        out.chain_rev  = rev_flag;
        out.seed_number= 0;
        out.score      = 0;

        for (uint32_t id : path_ids) {
            const auto& src = per_seg_chains[nodes[id].src_idx];
            out.seed_number += src.seed_number;
            out.score       += src.score;
            out.blocks.insert(out.blocks.end(), src.blocks.begin(), src.blocks.end());
        }
        anchors_coordinate_merge_(out.blocks);

        return out;
    };

    AnchorChain best_pos = solve_one_dir(false);
    AnchorChain best_neg = solve_one_dir(true);

    std::vector<AnchorChain> ret;
    // Return only the higher-score chain; push both if desired
    if (!best_pos.blocks.empty() && !best_neg.blocks.empty()) {
        ret.push_back( (best_pos.score >= best_neg.score) ? best_pos : best_neg );
    } else if (!best_pos.blocks.empty()) ret.push_back(best_pos);
    else if (!best_neg.blocks.empty())   ret.push_back(best_neg);
    return ret;
}



// Small helper: clamp to [0, hi]
static inline uint32_t clampu_(uint32_t x, uint32_t hi) {
    return x > hi ? hi : x;
}

// Seed-level graph-aware scoring:
// - Same segment -> reuse mm_comput_sc_()
// - Different segment -> graph shortest path + dq += len(prev_seg)
int32_t MinimizerIndex::mm_comput_sc_graph_(const MM128* ai, const MM128* aj) const {
    // Base read gap
    int32_t dq_base = (int32_t)ai->q_off() - (int32_t)aj->q_off();
    if (dq_base <= 0) return INT32_MIN;

    // Same segment: use original mm2 scoring (bandwidth/distance limits)
    if (ai->r_id() == aj->r_id()) {
        return mm_comput_sc_(
            ai, aj,
            chainOpts_.max_dist_x, chainOpts_.max_dist_y, chainOpts_.bw,
            chainOpts_.chn_pen_gap, chainOpts_.chn_pen_skip,
            chainOpts_.is_cdna, chainOpts_.n_seg
        );
    }

    // === Cross-segment case ===
    if (!graph_) return INT32_MIN;

    const bool rev = (ai->dir() != 0);   // 0:'+' 1:'-'
    const uint32_t k = (uint32_t)chainOpts_.k;

    // Segment length: use index sequence length
    auto seg_len = [&](uint32_t rid)->uint32_t {
        return (rid < seqs_.size()) ? (uint32_t)seqs_[rid].size() : 0u;
    };
    const uint32_t lenA = seg_len(aj->r_id());
    const uint32_t lenB = seg_len(ai->r_id());
    if (lenA == 0 || lenB == 0) return INT32_MIN;

    // Requested rule: dq += previous segment length
    int32_t dq = int32_t(ai->q_off()) - int32_t(aj->q_off());
    if (dq <= 0 || dq > (int64_t)chainOpts_.max_dist_x * 2) return INT32_MIN;

    // Oriented vertices
    const uint32_t ori    = rev ? 1u : 0u;
    const uint32_t v_from = (aj->r_id() << 1) | ori;
    const uint32_t v_to   = (ai->r_id() << 1) | ori;

    // Fast connectivity pruning
    if (!is_connected_vertex_(v_from, v_to)) return INT32_MIN;

    // forward -> oriented offsets
    // Start point (leave A covered tail):
    //   '+': off_from = aj->r_off() + k
    //   '-': off_from = lenA - aj->r_off()
    uint32_t off_from = rev
        ? clampu_(lenA >= aj->r_off() ? (lenA - aj->r_off()) : 0u, lenA)
        : clampu_(aj->r_off() + k, lenA);

    // Target point (enter B covered start):
    //   '+': off_to = ai->r_off()
    //   '-': off_to = lenB - (ai->r_off() + k)
    uint32_t off_to = rev
        ? ( (ai->r_off() + k <= lenB) ? (lenB - (ai->r_off() + k)) : 0u )
        : clampu_(ai->r_off(), lenB);

    // Graph shortest added-bases path (ow-aware)
    uint64_t dist = 0;
    if (!graph_->shortest_distance_between_offsets(v_from, off_from, v_to, off_to, dist, /*cap=*/200000))
        return INT32_MIN;
    if (dist > (uint64_t)INT32_MAX) dist = (uint64_t)INT32_MAX;
    const int32_t dr = (int32_t)dist;

    // mm2-style scoring (cross-segment as sidi!=sidj)
    const int32_t dd = (dr > dq) ? (dr - dq) : (dq - dr);
    const int32_t dg = (dr < dq) ? dr : dq;
    const int32_t q_span = (int32_t)k;

    int32_t sc = (q_span < dg ? q_span : dg);
    if (dd || dg > q_span) {
        float lin_pen = chainOpts_.chn_pen_gap * (float)dd + chainOpts_.chn_pen_skip * (float)dg;
        float log_pen = (dd >= 1) ? (float)mm_log2_i32_(dd + 1) : 0.0f;
        if (dr == 0) ++sc;  // Small reward for perfect junction
        else sc -= (int)std::min(lin_pen, log_pen);
    }
    return sc;
}

// ==================== New: segment -> 64-bit bitmap ====================
static inline uint64_t seg_bit_(uint32_t seg) noexcept {
    // Uses splitmix64 in this file; 64-bit bitmap, O(1) lookup
    return 1ull << (uint32_t)(splitmix64((uint64_t)seg) & 63);
}


// Same pattern as chain_dp(), but:
// - candidate window clipped by read coord (q_off), no same r_id restriction
// - scoring uses mm_comput_sc_graph_(), allowing cross-segment transitions
void MinimizerIndex::chain_dp_graph(
    const std::vector<MM128>& seeds,   // must be nondecreasing in read coords
    std::vector<MM128>& b_out,
    std::vector<Sc_Len>& u_out
) const {
    b_out.clear(); u_out.clear();
    const int64_t n = (int64_t)seeds.size();
    if (n == 0) return;

    if (chainOpts_.max_dist_x < chainOpts_.bw) chainOpts_.max_dist_x = chainOpts_.bw;
    if (chainOpts_.max_dist_y < chainOpts_.bw && !chainOpts_.is_cdna) chainOpts_.max_dist_y = chainOpts_.bw;
    const int32_t max_drop = chainOpts_.is_cdna ? INT32_MAX : chainOpts_.bw;

    std::vector<int32_t> f(n), v(n), t(n, 0);
    std::vector<int64_t> p(n);
    std::vector<uint64_t> used(n, 0); // bitmap: closed segments (already left)

    int64_t st = 0, max_ii = -1;
    int32_t mmax_f = 0;

    for (int64_t i = 0; i < n; ++i) {
        int64_t max_j = -1;
        int32_t q_span = seeds[i].q_span();
        int32_t max_f  = q_span;
        int32_t n_skip = 0;
        uint64_t best_used = 0; // best-path closed-seg bitmap for i

        // Simple window (by read distance & count)
        while (st < i && ( (int64_t)seeds[i].q_off() - (int64_t)seeds[st].q_off() > (int64_t)chainOpts_.max_dist_x * 2 ))
            ++st;
        if (i - st > chainOpts_.max_iter) st = i - chainOpts_.max_iter;

        int64_t j;
        for (j = i - 1; j >= st; --j) {
            const MM128* ai = &seeds[i];
            const MM128* aj = &seeds[j];

            uint64_t new_used;
            int32_t sc;

            if (ai->r_id() == aj->r_id()) {
                // Same segment: must be same direction; segment not closed
                if (ai->dir() != aj->dir()) continue;
                if (used[j] & seg_bit_(ai->r_id())) continue;

                const bool plus = (ai->dir() == 0);
                const int32_t dri = (int32_t)ai->r_off() - (int32_t)aj->r_off();
                if (plus ? (dri <= 0) : (dri >= 0)) continue;

                new_used = used[j];
                sc = mm_comput_sc_(
                    ai, aj,
                    chainOpts_.max_dist_x, chainOpts_.max_dist_y, chainOpts_.bw,
                    chainOpts_.chn_pen_gap, chainOpts_.chn_pen_skip,
                    chainOpts_.is_cdna, chainOpts_.n_seg
                );
            } else {
                // Cross-segment: do not return to closed seg; leaving aj closes aj
                if (used[j] & seg_bit_(ai->r_id())) continue;
                new_used = used[j] | seg_bit_(aj->r_id());
                sc = mm_comput_sc_graph_(ai, aj);
            }

            if (sc == INT32_MIN) continue;
            sc += f[j];

            if (sc > max_f) {
                max_f = sc; max_j = j;
                best_used = new_used;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == (int32_t)i) {
                if (++n_skip > chainOpts_.max_skip) break;
            }
            if (p[j] >= 0) t[p[j]] = (int32_t)i;
        }
        int64_t end_j = j;

        // Heuristic max_ii
        if (max_ii < 0 || (int64_t)seeds[i].q_off() - (int64_t)seeds[max_ii].q_off() > (int64_t)chainOpts_.max_dist_x) {
            int32_t best = INT32_MIN; max_ii = -1;
            for (j = i - 1; j >= st; --j) if (best < f[j]) { best = f[j]; max_ii = j; }
        }
        if (max_ii >= 0 && max_ii < end_j) {
            int32_t tmp;
            uint64_t tmp_used = 0;

            if (seeds[i].r_id() == seeds[max_ii].r_id()) {
                if (seeds[i].dir() != seeds[max_ii].dir()) tmp = INT32_MIN;
                else if (used[max_ii] & seg_bit_(seeds[i].r_id())) tmp = INT32_MIN;
                else {
                    // monotonic constraint within same segment
                    const bool plus = (seeds[i].dir() == 0);
                    const int32_t dri = (int32_t)seeds[i].r_off() - (int32_t)seeds[max_ii].r_off();
                    if (plus ? (dri <= 0) : (dri >= 0)) {
                        tmp = INT32_MIN;
                    } else {
                        tmp = mm_comput_sc_(
                            &seeds[i], &seeds[max_ii],
                            chainOpts_.max_dist_x, chainOpts_.max_dist_y, chainOpts_.bw,
                            chainOpts_.chn_pen_gap, chainOpts_.chn_pen_skip,
                            chainOpts_.is_cdna, chainOpts_.n_seg
                        );
                        tmp_used = used[max_ii];
                    }
                }
            } else {
                if (used[max_ii] & seg_bit_(seeds[i].r_id())) tmp = INT32_MIN;
                else {
                    tmp = mm_comput_sc_graph_(&seeds[i], &seeds[max_ii]);
                    tmp_used = used[max_ii] | seg_bit_(seeds[max_ii].r_id());
                }
            }

            if (tmp != INT32_MIN && max_f < tmp + f[max_ii]) {
                max_f = tmp + f[max_ii];
                max_j = max_ii;
                best_used = tmp_used;
            }
        }

        f[i] = max_f; p[i] = max_j;
        used[i] = best_used;                   // store closed-seg bitmap
        v[i] = (max_j >= 0 && v[max_j] > max_f) ? v[max_j] : max_f;

        if (max_ii < 0 || ( (int64_t)seeds[i].q_off() - (int64_t)seeds[max_ii].q_off() <= (int64_t)chainOpts_.max_dist_x && f[max_ii] < f[i]) )
            max_ii = i;

        if (mmax_f < max_f) mmax_f = max_f;
    }

    // Backtrack + compact sort
    int32_t n_u = 0, n_v = 0;
    std::vector<Sc_Len> u = mm_chain_backtrack_(n, f, p, v, t, chainOpts_.min_cnt, chainOpts_.min_sc, max_drop, n_u, n_v);
    if (n_u == 0) return;

    mm_compact_sort_(seeds, v, u, b_out, u_out);
}


} // namespace mmidx