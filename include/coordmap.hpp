#pragma once
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <optional>
#include <limits>
#include <cassert>
#include <cctype>
#include <cmath>
#include <unordered_set>

#include "logger.hpp"

namespace coordmap {

struct Interval { uint32_t beg=0, end=0; };  // half-open [beg,end)

struct Record {
    uint32_t a_ctg=0, b_ctg=0;
    Interval a{}, b{};      // Genomic coords
    bool     rev=false;     // Reverse strand
    int8_t   a_dir=+1;      // A-side consumption: +1 increasing, -1 decreasing
    int8_t   b_dir=+1;      // Same for B side
};

struct Hit {
    uint32_t ctg;        // Target contig id
    uint32_t beg, end;   // [beg,end) on target contig
    bool     rev;        // This hit is on the reverse strand
    uint32_t qbeg, qend; // [qbeg,qend) on query sequence
};

inline std::string to_string(const Hit& h, std::string_view ctg_name) {
    std::ostringstream os;
    os << ctg_name << ":" << h.beg << "-" << h.end << (h.rev?'-':'+');
    return os.str();
}

// ------------------------- Binner -------------------------
struct Binner {
    // 16K,128K,1M,8M,64M,512M (covers up to 2^29)
    static constexpr uint32_t BIN_FIRST_SHIFT = 14; // 16K
    static constexpr uint32_t BIN_NEXT_SHIFT  = 3;

    static inline void reg2bins(uint32_t beg, uint32_t end, std::vector<uint32_t>& vec) {
        if (end == 0) return;
        if (end > (1u<<29)) end = (1u<<29);
        --end; // [beg,end) -> [beg,end] for binning

        static constexpr uint32_t kLevelOffsets[] = {0, 1, 9, 73, 585, 4681};

        uint32_t s = beg >> BIN_FIRST_SHIFT;
        uint32_t e = end >> BIN_FIRST_SHIFT;

        for (uint32_t i = kLevelOffsets[0] + s; i <= kLevelOffsets[0] + e; ++i) vec.push_back(i);
        for (int lvl = 1; lvl < 6; ++lvl) {
            s >>= BIN_NEXT_SHIFT;
            e >>= BIN_NEXT_SHIFT;
            for (uint32_t i = kLevelOffsets[lvl] + s; i <= kLevelOffsets[lvl] + e; ++i) vec.push_back(i);
        }
    }
};

inline uint32_t pack_ref(uint32_t rec_idx, bool side_is_B) { return (rec_idx<<1) | (side_is_B?1u:0u); }
inline uint32_t unpack_rec(uint32_t ref) { return ref >> 1; }
inline bool     unpack_side_is_B(uint32_t ref) { return (ref & 1u)!=0; }

// ------------------------- parsing helpers -------------------------
struct Piece {               // Normalized piece on a single contig (e.g., utg000032l:1744739-1745955:0-3+ -> utg000032l:1744739-1744742+)
    uint32_t ctg = 0;
    uint32_t beg = 0, end = 0;
    int8_t   dir = +1;       // +1 increasing coordinates, -1 decreasing coordinates
    uint32_t len() const { return end - beg; }
};

struct Slice { uint32_t beg = 0, end = 0; bool has = false; };

struct SidePlan {              // e.g., utg000032l:1745955-1744739:0-3+
    std::vector<Piece> pieces; // Concatenated pieces in written order (beg, end, len, dir=+1/-1) (e.g., utg000032l, 1744739, 1745955, -1)
    Slice clip;                // Sub-interval on virtual sequence (e.g., 0, 3)
    bool  sign_plus = true;    // Trailing +/- (walk direction) (e.g., +)
    uint64_t total_len = 0;    // Total concatenated length
};

class CoordMap {
public:
    bool load(const std::string& path) {
        clear();
        std::ifstream in(path);
        if (!in) {
            error_stream() << path << ": No such file or directory" << std::endl;
            std::exit(1);
        }
        std::string line; line.reserve(512);
        uint32_t n_line = 0;

        std::vector<Record> tmp; tmp.reserve(1<<20);

        while (std::getline(in, line)) {
            ++n_line; trim(line);
            if (line.empty() || line[0]=='#') continue;

            std::string_view L, R;
            if (!split_two_tokens(line, L, R)) {
                error_stream() << path << ": Bad line (need two columns) at line " << n_line << std::endl;
                std::exit(1);
            }

            // parse sides
            SidePlan A, B;
            if (!parse_side(L, A, n_line)) return false;
            if (!parse_side(R, B, n_line)) return false;

            std::vector<Piece> As = slice_on_virtual(A);
            std::vector<Piece> Bs = slice_on_virtual(B);

            apply_sign_to_consumption(As, A.sign_plus);
            apply_sign_to_consumption(Bs, B.sign_plus);

            // Length check
            uint64_t lenA = total_len(As), lenB = total_len(Bs);
            if (lenA != lenB) {
                error_stream() << path << ": Length mismatch at line " << n_line << std::endl;
                std::exit(1);
            }
            if (lenA == 0) continue;

            fuse_to_records(As, Bs, tmp);
        }
        in.close();

        recs_ = std::move(tmp);
        normalize_records(recs_);
        build_bin_index();
        return true;
    }

    /**
     * Find hits for query region [qbeg, qend) on contig
     * 
     * @param ctg             Query contig name
     * @param qbeg            Query start position (0-based, inclusive)
     * @param qend            Query end position (0-based, exclusive)
     * @param max_hops        Maximum number of hops (0=direct hits only, 1=second layer, etc.)
     * @param max_fanout      Maximum number of hits to expand per node (limits branching)
     * @param min_len         Minimum absolute hit length in bp (shorter hits are filtered)
     * @param min_frac        Minimum relative hit length as fraction of query length
     * @param max_total_hits  Global hits limit to prevent explosion
     * @param merge_indels    Whether to merge indels to output a single hit
     * @return                Vector of Hit structs with target positions and orientations
     */
    std::vector<Hit> map_range(
        std::string_view ctg, uint32_t qbeg, uint32_t qend, 
        int max_hops = 15, size_t max_fanout = 512, uint32_t min_len = 15,
        double min_frac = 0.1, size_t max_total_hits = 2000, bool merge_indels = false
    ) const {
        std::vector<Hit> out;
        if (qbeg >= qend) return out;

        auto it = name2id_.find(std::string(ctg));
        if (it == name2id_.end()) return out;
        const uint32_t src_ctg = it->second;

        // visited: (ctg,beg,end,rev) to avoid loops
        struct Key { uint32_t ctg, beg, end; bool rev; };
        struct KeyHash {
            size_t operator()(const Key& k) const noexcept {
                uint64_t x = (uint64_t(k.ctg) << 32) ^ uint64_t(k.beg);
                x ^= (uint64_t(k.end) + 0x9e3779b97f4a7c15ULL + (x<<6) + (x>>2));
                x ^= (k.rev ? 0x9ae16a3b2f90404fULL : 0ULL);
                return (size_t)x;
            }
        };
        struct KeyEq {
            bool operator()(const Key& a, const Key& b) const noexcept {
                return a.ctg==b.ctg && a.beg==b.beg && a.end==b.end && a.rev==b.rev;
            }
        };
        std::unordered_map<Key, uint8_t, KeyHash, KeyEq> visited;
        visited.reserve(128);
        visited.emplace(Key{src_ctg, qbeg, qend, false}, 1);  // mark source

        // BFS queue
        struct Node {
            uint32_t ctg, beg, end;     // current node interval (on its contig)
            uint32_t src_beg, src_end;  // mapped back to "original source contig" sub-interval
            bool     rev_to_src;        // relative to source, whether reversed
            int      depth;             // number of hops
        };
        std::vector<Node> q; q.reserve(128);
        q.push_back(Node{src_ctg, qbeg, qend, qbeg, qend, false, 0});

        auto to_src_interval = [](const Node& cur, uint32_t sub_qbeg, uint32_t sub_qend) -> std::pair<uint32_t,uint32_t> {
            if (!cur.rev_to_src) {
                uint32_t s0 = cur.src_beg + (sub_qbeg - cur.beg);
                uint32_t s1 = cur.src_beg + (sub_qend - cur.beg);
                return {s0, s1};
            } else {
                uint32_t s0 = cur.src_end - (sub_qend - cur.beg);
                uint32_t s1 = cur.src_end - (sub_qbeg - cur.beg);
                return {s0, s1};
            }
        };

        auto min_keep_len = [&](uint32_t cur_len)->uint32_t {
            if (cur_len == 0) return 0;
            uint32_t m = min_len;
            if (min_frac > 0.0) {
                uint32_t f = (uint32_t)std::ceil(cur_len * min_frac);
                if (f > m) m = f;
            }
            if (m > cur_len) m = cur_len;  // For short queries
            if (m == 0) m = 1;
            return m;
        };

        size_t total_nodes = 0;
        for (size_t qi = 0; qi < q.size(); ++qi) {
            if (++total_nodes > max_total_hits) break;

            const Node cur = q[qi];
            const uint32_t cur_len = cur.end - cur.beg;

            // Hop: map current interval to other contigs
            auto step_hits = map_range_internal(cur.ctg, cur.beg, cur.end);
            if (step_hits.empty()) continue;

            // 1. First hop (depth==0): don't filter
            // 2. Next hops (depth>0): filter using min_keep_len(cur_len)
            if (cur.depth > 0) {
                const uint32_t min_keep = min_keep_len(cur_len);
                size_t out_i = 0;
                for (auto& h : step_hits) {
                    uint32_t tlen = h.end  - h.beg;
                    uint32_t qlen = h.qend - h.qbeg;
                    if (tlen >= min_keep && qlen >= min_keep)
                        step_hits[out_i++] = h;
                }
                step_hits.resize(out_i);
                if (step_hits.empty()) continue;
            }

            // Merge overlapping/adjacent hits to reduce fan-out
            step_hits = merge_adjacent_hits(step_hits);

            // Keep top K longest hits for expansion
            std::vector<size_t> order; order.reserve(step_hits.size());
            for (size_t i = 0; i < step_hits.size(); ++i) order.push_back(i);
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b){
                return (step_hits[a].end-step_hits[a].beg) > (step_hits[b].end-step_hits[b].beg);
            });

            size_t expand_limit = std::min(order.size(), max_fanout);

            for (size_t oi = 0; oi < order.size(); ++oi) {
                const Hit& h1 = step_hits[order[oi]];

                if (h1.ctg == src_ctg) continue;
                if (h1.ctg == cur.ctg) continue;

                const bool cum_rev = (cur.rev_to_src ^ h1.rev);

                Key k{h1.ctg, h1.beg, h1.end, cum_rev};
                if (visited.find(k) != visited.end()) continue;
                visited.emplace(k, 1);

                // Remap sub-interval back to source coordinates
                auto [src_sub_b, src_sub_e] = to_src_interval(cur, h1.qbeg, h1.qend);

                // output
                Hit out_h;
                out_h.ctg  = h1.ctg;
                out_h.beg  = h1.beg;
                out_h.end  = h1.end;
                out_h.rev  = cum_rev;
                out_h.qbeg = src_sub_b;
                out_h.qend = src_sub_e;
                out.push_back(out_h);

                // Depth limit
                if (cur.depth >= max_hops) continue;

                // Only allow the top K longest hits to enter the next level, limiting fan-out
                if (oi < expand_limit) {
                    q.push_back(Node{h1.ctg, h1.beg, h1.end, src_sub_b, src_sub_e, cum_rev, cur.depth+1});
                }
            }
        }

        // Compress the output again to reduce fragmentation
        out = merge_adjacent_hits(out);
        if (merge_indels) {
            std::vector<Hit> merged_hits = hits_chaining(out);
            if (merged_hits.empty()) return out;
            double fracTmp = double(merged_hits[0].qend - merged_hits[0].qbeg) / double(qend - qbeg);
            if (fracTmp >= min_frac) { return merged_hits; }
        }

        return out;
    }

    // Merge adjacent/overlapping hits on the same contig and strand
    std::vector<Hit> merge_adjacent_hits(std::vector<Hit> hits) const {
        if (hits.size() <= 1) return hits;

        // Sort by (rev, ctg, qbeg)
        std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
            if (a.rev  != b.rev)  return a.rev  < b.rev;
            if (a.ctg  != b.ctg)  return a.ctg  < b.ctg;
            if (a.qbeg != b.qbeg) return a.qbeg < b.qbeg;
            if (a.qend != b.qend) return a.qend < b.qend;
            if (a.beg  != b.beg)  return a.beg  < b.beg;
            return a.end < b.end;
        });

        std::vector<Hit> out;
        out.reserve(hits.size());

        Hit cur = hits[0];
        for (size_t i = 1; i < hits.size(); ++i) {
            const Hit& nx = hits[i];

            if (cur.ctg != nx.ctg || cur.rev != nx.rev) {  // Different contig or strand
                out.push_back(cur); cur = nx; continue;
            }

            uint32_t q_gap, t_gap;
            q_gap = cur.qend >= nx.qbeg ? cur.qend - nx.qbeg : UINT32_MAX;
            if (!cur.rev) { t_gap = cur.end >= nx.beg ? cur.end - nx.beg : UINT32_MAX; }  // forward
            else          { t_gap = nx.end >= cur.beg ? nx.end - cur.beg : UINT32_MAX; }  // reverse

            if (q_gap == UINT32_MAX || t_gap == UINT32_MAX || q_gap != t_gap) { out.push_back(cur); cur = nx; continue; }

            cur.qbeg = std::min(cur.qbeg, nx.qbeg);
            cur.qend = std::max(cur.qend, nx.qend);
            cur.beg  = std::min(cur.beg , nx.beg );
            cur.end  = std::max(cur.end , nx.end );
        }
        out.push_back(cur);
        return out;
    }

    // DP and chaining
    static inline std::vector<Hit> hits_chaining(
        const std::vector<Hit>& hits,
        uint32_t max_lookback = 10,
        uint32_t max_gap = 20000,
        uint32_t max_off_diff = 500,
        double gap_pen =0.1
    ) {
        auto len = [](const Hit& h)->uint32_t {
            uint32_t qlen = (h.qend > h.qbeg) ? (h.qend - h.qbeg) : 0u;
            qlen = std::max(qlen, 30u);
            return qlen;
        };
        auto offset = [](const Hit& h)->int64_t {
            return (!h.rev) ? (int64_t)h.beg - (int64_t)h.qbeg : (int64_t)h.end + (int64_t)h.qbeg;
        };
        auto iabs64 = [](int64_t x)->uint64_t { return (uint64_t)(x < 0 ? -x : x); };

        // copy and sort
        std::vector<Hit> v = hits;
        v.erase(std::remove_if(v.begin(), v.end(), [&](const Hit& h){ return len(h) == 0; }), v.end());
        std::sort(v.begin(), v.end(), [&](const Hit& a, const Hit& b){
            if (a.ctg  != b.ctg)  return a.ctg  < b.ctg;
            if (a.rev  != b.rev)  return a.rev  < b.rev;
            if (a.qbeg != b.qbeg) return a.qbeg < b.qbeg;
            if (a.qend != b.qend) return a.qend < b.qend;
            return (!a.rev) ? (a.beg < b.beg) : (a.end < b.end);
        });

        std::vector<Hit> out;
        if (v.empty()) return out;

        const int NEG = std::numeric_limits<int>::min() / 4;
        std::vector<uint8_t> alive(v.size(), 1);
        size_t alive_cnt = v.size();

        // Chaining loop
        while (alive_cnt) {
            int best_score = NEG;
            int best_i = -1;

            const size_t n = v.size();
            std::vector<int> dp(n, NEG), pre(n, -1);

            // DP
            for (size_t i = 0; i < n; ++i) {
                if (!alive[i]) continue;
                const Hit& hi = v[i];
                const uint32_t li = len(hi);
                dp[i] = (int)li;

                const size_t j0 = (i > max_lookback ? i - max_lookback : 0);
                for (size_t jj = i; jj-- > j0; ) {
                    if (!alive[jj]) continue;
                    const Hit& hj = v[jj];
                    if (hj.ctg != hi.ctg || hj.rev != hi.rev) continue;

                    if (hi.qbeg < hj.qend) continue;
                    const uint32_t qgap = hi.qbeg - hj.qend;
                    if (qgap > max_gap) continue;

                    uint32_t tgap = 0;
                    if (!hi.rev) {
                        if (hi.beg < hj.end) continue;
                        tgap = hi.beg - hj.end;
                    } else {
                        if (hj.beg < hi.end) continue;
                        tgap = hj.beg - hi.end;
                    }
                    if (tgap > max_gap) continue;

                    if (qgap != 0 && tgap != 0 && iabs64(offset(hi) - offset(hj)) > max_off_diff) continue;

                    double pen_d = 0.0;
                    if (qgap == 0 || tgap == 0) {
                        pen_d = 0.1 * std::sqrt((double)(qgap + tgap));
                    } else {
                        pen_d = gap_pen * (double)(qgap + tgap);
                    }
                    const int trans = (int)li - (int)std::lround(pen_d);

                    if (dp[jj] != NEG && dp[jj] + trans > dp[i]) {
                        dp[i] = dp[jj] + trans;
                        pre[i] = (int)jj;
                    }
                }

                if (dp[i] > best_score) { best_score = dp[i]; best_i = (int)i; }
            }

            if (best_i < 0) break;

            // Backtrack
            uint32_t qmin = v[(size_t)best_i].qbeg, qmax = v[(size_t)best_i].qend;
            uint32_t tmin = std::min(v[(size_t)best_i].beg, v[(size_t)best_i].end);
            uint32_t tmax = std::max(v[(size_t)best_i].beg, v[(size_t)best_i].end);

            Hit merged = v[(size_t)best_i];

            for (int k = best_i; k >= 0; k = pre[(size_t)k]) {
                if (!alive[(size_t)k]) break;
                const Hit& h = v[(size_t)k];

                qmin = std::min(qmin, h.qbeg);
                qmax = std::max(qmax, h.qend);
                tmin = std::min(tmin, std::min(h.beg, h.end));
                tmax = std::max(tmax, std::max(h.beg, h.end));

                alive[(size_t)k] = 0;
                --alive_cnt;

                if (pre[(size_t)k] < 0) break;
            }

            merged.qbeg = qmin; merged.qend = qmax;
            merged.beg  = tmin; merged.end  = tmax;
            out.push_back(merged);
        }

        // Sort by length descending
        std::sort(out.begin(), out.end(), [](const Hit& a, const Hit& b){
            uint32_t alen = a.qend - a.qbeg;
            uint32_t blen = b.qend - b.qbeg;
            if (alen != blen) return alen > blen;
            if (a.ctg != b.ctg) return a.ctg < b.ctg;
            if (a.rev != b.rev) return a.rev < b.rev;
            return a.qbeg < b.qbeg;
        });

        return out;
    }

    // Internal data
    const std::string& contig_name(uint32_t id) const { return id < id2name_.size()? id2name_[id] : kEmpty_; }
    const std::uint32_t contig_id(std::string_view name) const {
        auto it = name2id_.find(std::string(name));
        return it != name2id_.end() ? it->second : std::numeric_limits<uint32_t>::max();
    }
    size_t num_records() const { return recs_.size(); }

    void clear() {
        recs_.clear();
        name2id_.clear(); id2name_.clear();
        idx_.clear();
    }

private:
    // Trim leading and trailing whitespace
    static inline void trim(std::string& s) {
        auto ns = [](unsigned char c){ return !std::isspace(c); };
        auto b = std::find_if(s.begin(), s.end(), ns);
        auto e = std::find_if(s.rbegin(), s.rend(), ns).base();
        if (b >= e) { s.clear(); return; }
        s.assign(b, e);
    }
    // Split line into two whitespace-delimited tokens (e.g., foo bar baz -> A=foo, B=bar)
    static inline bool split_two_tokens(const std::string& line, std::string_view& A, std::string_view& B) {
        size_t i = 0, n = line.size();
        while (i < n && std::isspace((unsigned char)line[i])) ++i;
        size_t j = i; while (j < n && !std::isspace((unsigned char)line[j])) ++j;
        if (i == j) return false; A = std::string_view(line).substr(i, j - i);
        while (j < n && std::isspace((unsigned char)line[j])) ++j;
        size_t k = j; while (k < n && !std::isspace((unsigned char)line[k])) ++k;
        if (j == k) return false; B = std::string_view(line).substr(j, k - j);
        return true;
    }
    // Trim leading and trailing whitespace from a std::string_view in place
    static inline void trim_view(std::string_view& v) {
        while (!v.empty() && std::isspace((unsigned char)v.front())) v.remove_prefix(1);
        while (!v.empty() && std::isspace((unsigned char)v.back()))  v.remove_suffix(1);
    }
    // Convert a decimal string_view to uint32_t, with whitespace trimming
    static bool sv_to_u32(std::string_view sv, uint32_t& out) {
        uint64_t x = 0;
        trim_view(sv);
        if (sv.empty()) return false;
        for (char c : sv) {
            if (c < '0' || c > '9') return false;
            x = x*10 + (c - '0');
            if (x > std::numeric_limits<uint32_t>::max()) return false;
        }
        out = (uint32_t)x; return true;
    }
    // Parse "num-num" format (e.g., "12-34") into two uint32_t values.
    static bool parse_numdash_num(std::string_view sv, uint32_t& a, uint32_t& b) {
        trim_view(sv);
        size_t m = sv.find('-'); if (m == std::string_view::npos) return false;
        std::string_view s1 = sv.substr(0, m), s2 = sv.substr(m + 1);
        trim_view(s1); trim_view(s2);
        return sv_to_u32(s1, a) && sv_to_u32(s2, b);
    }

    // Intern a contig name to a stable numeric id (create if missing)
    uint32_t intern_ctg(const std::string& name) {
        auto it = name2id_.find(name);
        if (it != name2id_.end()) return it->second;
        uint32_t id = (uint32_t)id2name_.size();
        id2name_.push_back(name);
        name2id_.emplace(name, id);
        if (id >= idx_.size()) idx_.resize(id+1);
        return id;
    }

    // Parse one side expression into a SidePlan.
    //  e.g., h1tg000001l:347957-348031;h1tg000002l:12400-12345:0-50+
    //         - [h1tg000001l, 347957,348031,+1], [h1tg000002l,12345,12400,-1], clip=[0,50], dir=+
    bool parse_side(std::string_view sv, SidePlan& out, uint32_t line_no) {
        out = SidePlan{};
        trim_view(sv);
        if (sv.empty()) {
            error_stream() << "Empty side at line " << line_no << std::endl;
            std::exit(1);
        }

        // directional
        char tail = sv.back();
        if (tail == '+' || tail == '-') { out.sign_plus = (tail == '+'); sv.remove_suffix(1); trim_view(sv); }

        // Split concatenation pieces by ';'.
        std::vector<std::string_view> parts;
        size_t start = 0;
        while (start < sv.size()) {
            size_t pos = sv.find(';', start);
            if (pos == std::string_view::npos) { parts.push_back(sv.substr(start)); break; }
            parts.push_back(sv.substr(start, pos - start));
            start = pos + 1;
        }
        if (parts.empty()) {
            error_stream() << "Bad side at line " << line_no << std::endl;
            std::exit(1);
        }

        // Detect optional global slice in the last piece (e.g., h1tg000002l:12400-12345:0-50+)
        std::string_view last = parts.back(); trim_view(last);
        size_t first_colon = last.find(':');
        if (first_colon == std::string_view::npos) {
            error_stream() << "Missing ':' in side at line " << line_no << std::endl;
            std::exit(1);
        }
        size_t last_colon = last.rfind(':');
        std::optional<Slice> gslice;
        if (last_colon != first_colon) {
            std::string_view sub = last.substr(last_colon + 1);
            uint32_t sb = 0, se = 0;
            if (!parse_numdash_num(sub, sb, se) || sb > se) {
                error_stream() << "Bad slice at line " << line_no << std::endl;
                std::exit(1);
            }
            gslice = Slice{sb, se, true};
            parts.back() = last.substr(0, last_colon);
        }

        // Parse each piece as name:range
        for (auto p : parts) {
            trim_view(p);
            size_t c = p.find(':');
            if (c == std::string_view::npos) {
                error_stream() << "Bad piece (no ':') at line " << line_no << std::endl;
                std::exit(1);
            }
            std::string name(p.substr(0,c));
            trim(name);
            std::string_view rg = p.substr(c + 1); trim_view(rg);

            uint32_t x = 0, y = 0;
            if (!parse_numdash_num(rg, x, y) || x == y) {
                error_stream() << "Bad or zero-length range at line " << line_no << std::endl;
                std::exit(1);
            }
            uint32_t lo = std::min(x, y), hi = std::max(x, y);
            int8_t dir  = (y > x ? +1 : -1);
            Piece pc;
            pc.ctg = intern_ctg(name);
            pc.beg = lo; pc.end = hi; pc.dir = dir;
            out.pieces.push_back(pc);
            out.total_len += (hi-lo);
        }
        if (gslice.has_value()) out.clip = *gslice;
        else out.clip = Slice{0, (uint32_t)out.total_len, false};
        return true;
    }

    // Find the real pieces corresponding to the virtual slice [clip.beg, clip.end) on the concatenated sequence.
    // Example:
    //   pieces (in virtual L->R order):
    //     - [A, 100,110, +1]   // len=10, virtual span [0,10)
    //     - [B, 200,208, -1]   // len=8,  virtual span [10,18)  (reverse on real coords)
    //     - [C,  50, 70, +1]   // len=20, virtual span [18,38)
    //   clip = [6,26)          // slice on the virtual concatenated sequence
    //   output sub-pieces (still virtual L->R order):
    //     - [A, 106,110, +1]   // takes virtual [6,10)  from piece1  => real [100+6, 100+10)
    //     - [B, 200,208, -1]   // takes virtual [10,18) from piece2  => real [end-8, end-0) = [208-8,208)
    //     - [C,  50, 58, +1]   // takes virtual [18,26) from piece3  => real [50+0, 50+8)
    static std::vector<Piece> slice_on_virtual(const SidePlan& sp) {
        uint64_t total = sp.total_len;
        uint32_t qs = sp.clip.has ? sp.clip.beg : 0;
        uint32_t qe = sp.clip.has ? sp.clip.end : (uint32_t)total;
        if (qs > qe) std::swap(qs, qe);
        if (qe > total) qe = (uint32_t)total;

        std::vector<Piece> out;
        if (qs == qe) return out;

        uint64_t acc = 0;
        for (const auto& pc : sp.pieces) {
            uint32_t L = pc.len();
            uint64_t seg_s = acc, seg_e = acc + L;
            acc = seg_e;

            // Intersect with [qs, qe)
            if (!(qs < seg_e && seg_s < qe)) continue;
            uint32_t off0 = (uint32_t)(std::max<uint64_t>(qs, seg_s) - seg_s);
            uint32_t off1 = (uint32_t)(std::min<uint64_t>(qe, seg_e) - seg_s);
            // Map [off0,off1) (virtual offset) to "real coordinates"
            Piece sub; sub.ctg = pc.ctg; sub.dir = pc.dir;
            if (pc.dir > 0) {
                sub.beg = pc.beg + off0;
                sub.end = pc.beg + off1;
            } else {
                sub.beg = pc.end - off1;
                sub.end = pc.end - off0;
            }
            out.push_back(sub);
        }
        return out;
    }

    // Example:
    //   v = [ [A,100,110,+1], [B,200,208,-1] ]
    //   sign_plus=false => v becomes [ [B,200,208,+1], [A,100,110,-1] ].
    static void apply_sign_to_consumption(std::vector<Piece>& v, bool sign_plus) {
        if (sign_plus) return;
        std::reverse(v.begin(), v.end());
        for (auto& p : v) p.dir = -p.dir;
    }
    // Sum of lengths of pieces
    static uint64_t total_len(const std::vector<Piece>& v) {
        uint64_t s = 0; for (auto& p : v) s += (p.end - p.beg); return s;
    }
    // Example:
    //   A: [ [A1,100,110,+1], [A2,200,206,-1] ]   (total 10+6=16)
    //   B: [ [B1,500,508,+1], [B2,900,908,-1] ]   (total  8+8=16)
    //   Output Records (delta chunks): 8, 2, 6
    //     r1: A=[100,108) dir+   B=[500,508) dir+
    //     r2: A=[108,110) dir+   B=[906,908) dir-
    //     r3: A=[200,206) dir-   B=[900,906) dir-
    void fuse_to_records(const std::vector<Piece>& A, const std::vector<Piece>& B, std::vector<Record>& outRecs) {
        size_t ia = 0, ib = 0;
        uint32_t a_pos = (A.empty() ? 0 : (A[0].dir > 0 ? A[0].beg : A[0].end));
        uint32_t b_pos = (B.empty() ? 0 : (B[0].dir > 0 ? B[0].beg : B[0].end));

        while (ia < A.size() && ib < B.size()) {
            const Piece& pa = A[ia];
            const Piece& pb = B[ib];
            assert(pa.dir == +1 || pa.dir == -1);
            assert(pb.dir == +1 || pb.dir == -1);
            if (pa.ctg == pb.ctg) {
                if (warned_same_ctg_.insert(pa.ctg).second){
                    warning_stream() << "Reference and query contig names must be different: " << contig_name(pa.ctg) << "\n";
                }
            }

            uint32_t a_rem = (pa.dir > 0 ? (pa.end - a_pos) : (a_pos - pa.beg));
            uint32_t b_rem = (pb.dir > 0 ? (pb.end - b_pos) : (b_pos - pb.beg));
            uint32_t delta = std::min(a_rem, b_rem);
            assert(delta > 0);

            Record r; r.a_ctg = pa.ctg; r.b_ctg = pb.ctg;
            r.a_dir = pa.dir; r.b_dir = pb.dir;
            r.rev   = (r.a_dir != r.b_dir);

            // Take delta from A
            if (pa.dir > 0) { r.a.beg = a_pos;         r.a.end = a_pos + delta; a_pos += delta; }
            else            { r.a.beg = a_pos - delta; r.a.end = a_pos;         a_pos -= delta; }

            // Take delta from B
            if (pb.dir > 0) { r.b.beg = b_pos;         r.b.end = b_pos + delta; b_pos += delta; }
            else            { r.b.beg = b_pos - delta; r.b.end = b_pos;         b_pos -= delta; }

            outRecs.push_back(r);

            // Move to next piece when the current one is fully consumed
            if ((pa.dir > 0 && a_pos == pa.end) || (pa.dir < 0 && a_pos == pa.beg)) {
                ++ia;
                if (ia < A.size()) a_pos = (A[ia].dir > 0 ? A[ia].beg : A[ia].end);
            }
            if ((pb.dir > 0 && b_pos == pb.end) || (pb.dir < 0 && b_pos == pb.beg)) {
                ++ib;
                if (ib < B.size()) b_pos = (B[ib].dir > 0 ? B[ib].beg : B[ib].end);
            }
        }
    }

    // ==================== Index ====================
    struct ContigIndex { std::unordered_map<uint32_t, std::vector<uint32_t>> bin2refs; };  // key=bin, value=packed ref ids

    static void normalize_records(std::vector<Record>& v) {
        if (v.empty()) return;

        struct Item { Record r; long long inv; };
        std::vector<Item> w; w.reserve(v.size());

        auto inv_of = [](const Record& r) -> long long {
            if (r.a_dir == r.b_dir) {
                // forward: pB = pA - d，d = (a_dir>0 ? A.beg-B.beg : A.end-B.end)
                return (r.a_dir > 0)
                    ? (long long)r.a.beg - (long long)r.b.beg
                    : (long long)r.a.end - (long long)r.b.end;
            } else {
                // reverse: pB = -pA + s，s = (a_dir>0 ? A.beg+B.end : A.end+B.beg)
                return (r.a_dir > 0)
                    ? (long long)r.a.beg + (long long)r.b.end
                    : (long long)r.a.end + (long long)r.b.beg;
            }
        };

        for (auto& r : v) w.push_back({r, inv_of(r)});

        std::sort(w.begin(), w.end(), [](const Item& a, const Item& b){
            if (a.r.a_ctg != b.r.a_ctg) return a.r.a_ctg < b.r.a_ctg;
            if (a.r.b_ctg != b.r.b_ctg) return a.r.b_ctg < b.r.b_ctg;
            if (a.r.a_dir != b.r.a_dir) return a.r.a_dir < b.r.a_dir;
            if (a.r.b_dir != b.r.b_dir) return a.r.b_dir < b.r.b_dir;
            if (a.inv != b.inv)         return a.inv < b.inv;
            if (a.r.a.beg != b.r.a.beg) return a.r.a.beg < b.r.a.beg;
            return a.r.a.end < b.r.a.end;
        });

        auto same_group = [](const Item& x, const Item& y){
            return x.r.a_ctg == y.r.a_ctg && x.r.b_ctg == y.r.b_ctg && x.r.a_dir == y.r.a_dir && x.r.b_dir == y.r.b_dir && x.inv == y.inv;
        };

        auto map_A_to_B = [](const Record& proto, long long inv, uint32_t Abeg, uint32_t Aend) -> std::pair<uint32_t,uint32_t> {
            if (proto.a_dir == proto.b_dir) {
                long long Bb = (long long)Abeg - inv;
                long long Be = (long long)Aend - inv;
                return { (uint32_t)Bb, (uint32_t)Be };
            } else {
                // reverse: B = [-Aend + s, -Abeg + s]
                long long Bb = -(long long)Aend + inv;
                long long Be = -(long long)Abeg + inv;
                if (Bb < Be) return { (uint32_t)Bb, (uint32_t)Be };
                else         return { (uint32_t)Be, (uint32_t)Bb };
            }
        };

        std::vector<Record> out; out.reserve(w.size());
        Item cur = w[0];

        for (size_t i = 1; i < w.size(); ++i) {
            const Item& nx = w[i];
            if (same_group(cur, nx) && nx.r.a.beg <= cur.r.a.end) {
                uint32_t Ab = std::min(cur.r.a.beg, nx.r.a.beg);
                uint32_t Ae = std::max(cur.r.a.end, nx.r.a.end);
                cur.r.a.beg = Ab; cur.r.a.end = Ae;
                auto [Bb, Be] = map_A_to_B(cur.r, cur.inv, Ab, Ae);
                cur.r.b.beg = Bb; cur.r.b.end = Be;
            } else {
                out.push_back(cur.r);
                cur = nx;
            }
        }
        out.push_back(cur.r);
        v.swap(out);
    }

    void build_bin_index() {
        size_t approx_bins = recs_.size() * 10u;
        for (auto& ci : idx_) ci.bin2refs.reserve(approx_bins / std::max<size_t>(1, idx_.size()));

        std::vector<uint32_t> bins; bins.reserve(64);
        for (uint32_t i = 0; i < recs_.size(); ++i) {
            const Record& r = recs_[i];
            // A -> B
            bins.clear(); Binner::reg2bins(r.a.beg, r.a.end, bins);
            auto& Amap = idx_[r.a_ctg].bin2refs;
            for (uint32_t b : bins) Amap[b].push_back(pack_ref(i, false));
            // B -> A
            bins.clear(); Binner::reg2bins(r.b.beg, r.b.end, bins);
            auto& Bmap = idx_[r.b_ctg].bin2refs;
            for (uint32_t b : bins) Bmap[b].push_back(pack_ref(i, true));
        }
    }

    static inline bool overlaps(uint32_t a0, uint32_t a1, uint32_t b0, uint32_t b1) {
        return (a0 < b1) && (b0 < a1);
    }

    // query
    std::vector<Hit> map_range_internal(uint32_t q_ctg, uint32_t qbeg, uint32_t qend) const {
        std::vector<Hit> out;
        if (qbeg >= qend) return out;

        std::vector<uint32_t> qb; qb.reserve(16);
        Binner::reg2bins(qbeg, qend, qb);

        std::vector<uint32_t> cands; cands.reserve(qb.size()*4);
        const auto& cmap = idx_[q_ctg].bin2refs;
        for (uint32_t b : qb) {
            auto it = cmap.find(b);
            if (it == cmap.end()) continue;
            const auto& v = it->second;
            cands.insert(cands.end(), v.begin(), v.end());
        }
        if (cands.empty()) return out;
        std::sort(cands.begin(), cands.end());
        cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

        out.reserve(cands.size());
        for (uint32_t ref : cands) {
            uint32_t rid = unpack_rec(ref);
            bool side_is_B = unpack_side_is_B(ref);
            const Record& r = recs_[rid];

            bool query_on_A = (!side_is_B && r.a_ctg == q_ctg);
            bool query_on_B = ( side_is_B && r.b_ctg == q_ctg);
            if (!query_on_A && !query_on_B) continue;

            const Interval& Iq = query_on_A ? r.a : r.b;
            if (!overlaps(qbeg, qend, Iq.beg, Iq.end)) continue;

            uint32_t ibeg = std::max(qbeg, Iq.beg);
            uint32_t iend = std::min(qend, Iq.end);

            int8_t qdir = query_on_A ? r.a_dir : r.b_dir;
            int8_t odir = query_on_A ? r.b_dir : r.a_dir;
            const Interval& Io = query_on_A ? r.b : r.a;

            uint32_t off0, off1;
            if (qdir > 0) {
                off0 = ibeg - Iq.beg;
                off1 = iend - Iq.beg;
            } else {
                off0 = Iq.end - iend;
                off1 = Iq.end - ibeg;
            }

            Hit h;
            h.ctg = query_on_A ? r.b_ctg : r.a_ctg;
            h.rev = (qdir != odir);

            h.qbeg = ibeg;
            h.qend = iend;

            if (odir > 0) {
                h.beg = Io.beg + off0;
                h.end = Io.beg + off1;
            } else {
                h.beg = Io.end - off1;
                h.end = Io.end - off0;
            }
            out.push_back(h);
        }
        return out;
    }

private:
    std::vector<Record> recs_;
    std::unordered_map<std::string, uint32_t> name2id_;
    std::vector<std::string> id2name_;
    std::vector<ContigIndex> idx_;
    static inline const std::string kEmpty_{};
    inline static std::unordered_set<uint32_t> warned_same_ctg_;
};

} // namespace coordmap