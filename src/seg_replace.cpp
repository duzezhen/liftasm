#include "../include/seg_replace.hpp"

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <stdexcept>
#include <tuple>
#include <fstream>

namespace SegReplace {

// ───────────────────────────── Interval ─────────────────────────────
Interval::u128 Interval::pack(
    uint64_t chrId, uint32_t beg, uint32_t end, bool strand
) {
    if (chrId >= (1ULL << 63)){
        error_stream() << "chr id out of range (< 2^63 required)\n";
        std::exit(1);
    }
    if (beg > 0xFFFFFFFFu || end > 0xFFFFFFFFu) {
        error_stream() << "pos out of range (<= 2^32-1 required)\n";
        std::exit(1);
    }
    if (beg >= end) {
        error_stream() << "Invalid range: chrId=" << chrId << " beg=" << beg << " end=" << end << "\n";
        std::exit(1);
    }

    u128 v = 0;
    v |= static_cast<u128>(chrId) << 65;
    v |= static_cast<u128>(beg)   << 33;
    v |= static_cast<u128>(end)   << 1;
    v |= static_cast<u128>(strand ? 1 : 0);
    return v;
}

uint64_t        Interval::seg_id(u128 v)        { return static_cast<uint64_t>(v >> 65); }
uint32_t        Interval::beg(u128 v)           { return static_cast<uint32_t>((v >> 33) & 0xFFFFFFFFu); }
uint32_t        Interval::end(u128 v)           { return static_cast<uint32_t>((v >> 1)  & 0xFFFFFFFFu); }
bool            Interval::is_reverse(u128 v)    { return static_cast<bool>(v & 1); }
uint32_t        Interval::len(u128 v)           { return end(v) - beg(v); }
Interval::u128  Interval::toggle_strand(u128 v) { return v ^ static_cast<u128>(1); }
Interval::u128  Interval::canonical(u128 v)     { return is_reverse(v) ? toggle_strand(v) : v; }

std::vector<Interval::u128> Interval::toggle_strand_chain(std::vector<Interval::u128> vs) {  // sample to reverse_and_toggle();
    std::reverse(vs.begin(), vs.end());
    for (auto& v : vs) v = Interval::toggle_strand(v);
    return vs;
}

std::string Interval::to_string(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x) {
        uint32_t d = static_cast<uint32_t>(x % 10);
        s.push_back(static_cast<char>('0' + d));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::string Interval::format(u128 v, const std::vector<std::string>& names) {
    std::string s;
    uint64_t chrId = seg_id(v);
    if (chrId < names.size()) s += names[chrId] + ":";
    else s += std::to_string(seg_id(v)) + ":";
    s += std::to_string(beg(v)) + "-";
    s += std::to_string(end(v));
    s += is_reverse(v) ? "-" : "+";
    return s;
}
std::string Interval::format(u128 v) {
    std::string s;
    s += std::to_string(seg_id(v)) + ":";
    s += std::to_string(beg(v)) + "-";
    s += std::to_string(end(v));
    s += is_reverse(v) ? "-" : "+";
    return s;
}

std::string Interval::name(u128 v) {
    std::string s;
    s += std::to_string(seg_id(v)) + ":";
    s += std::to_string(beg(v)) + "-";
    s += std::to_string(end(v));
    return s;
}

// ───────────────────────────── Expander ─────────────────────────────
void Expander::reverse_and_toggle(Expansion& v) {
    std::reverse(v.begin(), v.end());
    for (auto& e : v) e = Interval::toggle_strand(e);
}

void Expander::build_index() {
    std::unordered_set<Seg, U128Hash, U128Eq> all;
    for (const auto& kv : *rules_) {
        all.insert(kv.first);
        for (const auto& x : kv.second) all.insert(x);
    }
    for (const auto& s : all) (void)expand_(s);
}

void Expander::remove_from_index_by_key(const Seg& s) {
    memo_.erase(s);
    Seg sr = Interval::toggle_strand(s);
    memo_.erase(sr);
}

Expansion Expander::query(Seg s) const {
    if (auto it = memo_.find(s); it != memo_.end())
        return it->second;

    Seg sr = Interval::toggle_strand(s);
    if (auto it = memo_.find(sr); it != memo_.end()) {
        Expansion out = it->second;
        reverse_and_toggle(out);
        return out;
    }
    return Expansion{ s };  // Fallback: identity
}

void Expander::print_expansion_direct(const Seg& seg, const Expansion& v) const {
    std::string log = Interval::format(seg, names_) + " =>";
    for (const auto& x : v) log += " " + Interval::format(x, names_);
    log += "\n";
    log_stream() << log;
}

void Expander::print_expansion(const Seg& seg) {
    Expansion v = query(seg);
    print_expansion_direct(seg, v);
}

    // sort index and print
void Expander::print_index() const {
    std::vector<std::pair<Seg, Expansion>> items(memo_.begin(), memo_.end());

    auto seg_less = [](const auto& a, const auto& b) {
        return std::tuple<uint64_t, uint32_t, uint32_t, bool>(
                    Interval::seg_id(a.first),
                    Interval::beg(a.first),
                    Interval::end(a.first),
                    Interval::is_reverse(a.first)) <
                std::tuple<uint64_t, uint32_t, uint32_t, bool>(
                    Interval::seg_id(b.first),
                    Interval::beg(b.first),
                    Interval::end(b.first),
                    Interval::is_reverse(b.first));
    };
    std::sort(items.begin(), items.end(), seg_less);

    for (const auto& kv : items) print_expansion_direct(kv.first, kv.second);
}

// Expand a segment
Expansion Expander::expand_(Seg s) {
    const bool is_rev = Interval::is_reverse(s);
    if (is_rev) s = Interval::toggle_strand(s);

    // Memo lookup
    if (auto it = memo_.find(s); it != memo_.end()) {
        Expansion out = it->second;
        if (is_rev) reverse_and_toggle(out);
        return out;
    }

    // Cycle detection
    if (in_stack_.count(s)) {
        error_stream() << "Cycle detected at: " << Interval::format(s) << "\n";
        std::exit(1);
    }
    in_stack_.insert(s);

    Expansion out;
    auto rit = rules_->find(s);
    if (rit == rules_->end() || rit->second.empty()) {
        // Leaf
        out.push_back(s);
    } else {
        // Expand children & concatenate
        for (const auto& child : rit->second) {
            Expansion sub = expand_(child);
            out.insert(out.end(), sub.begin(), sub.end());
        }
        // Remove consecutive duplicates
        out.erase(std::unique(out.begin(), out.end()), out.end());
    }

    in_stack_.erase(s);
    // Cache forward-strand results (non-trivial only)
    if (!is_rev && !out.empty() && out.front() != s)
        memo_.emplace(s, out);

    if (is_rev) reverse_and_toggle(out);
    return out;
}


bool Expander::save_map(const std::string& path) {
    if (memo_.empty()) return false;

    std::ofstream ofs(path, std::ios::out | std::ios::binary);
    if (!ofs) return false;

    std::vector<std::pair<Seg, Expansion>> items(memo_.begin(), memo_.end());

    auto seg_less = [](const auto& a, const auto& b) {
        return std::tuple<uint64_t, uint32_t, uint32_t, bool>(
                    Interval::seg_id(a.first),
                    Interval::beg(a.first),
                    Interval::end(a.first),
                    Interval::is_reverse(a.first)) <
               std::tuple<uint64_t, uint32_t, uint32_t, bool>(
                    Interval::seg_id(b.first),
                    Interval::beg(b.first),
                    Interval::end(b.first),
                    Interval::is_reverse(b.first));
    };
    std::sort(items.begin(), items.end(), seg_less);

    for (const auto& kv : items) {
        const Seg        s   = kv.first;
        const Expansion& exp = kv.second;

        const uint64_t sid = Interval::seg_id(s);
        const uint32_t sb  = Interval::beg(s);
        const uint32_t se  = Interval::end(s);
        const bool     sr  = Interval::is_reverse(s);
        const uint32_t sl  = se - sb;
        if (sl == 0 || exp.empty()) continue;

        uint32_t used = 0;
        for (Seg t : exp) {
            if (used >= sl) break;

            const uint32_t tl = Interval::len(t);
            if (tl == 0) continue;

            const uint32_t take = (sl - used < tl) ? (sl - used) : tl;

            // Source sub-interval
            Seg s_sub = Interval::pack(sid, sb + used, sb + used + take, sr);

            // Target sub-interval
            Seg t_sub;
            if (take == tl) {
                t_sub = t;
            } else {
                const uint64_t tid = Interval::seg_id(t);
                const uint32_t tb  = Interval::beg(t);
                const uint32_t te  = Interval::end(t);
                const bool     tr  = Interval::is_reverse(t);
                if (!tr)
                    t_sub = Interval::pack(tid, tb, tb + take, false);
                else
                    t_sub = Interval::pack(tid, te - take, te, true);
            }

            ofs << Interval::format(s_sub, names_) << '\t'
                << Interval::format(t_sub, names_) << '\n';

            used += take;
        }
    }

    ofs.flush();
    return static_cast<bool>(ofs);
}

} // namespace SegReplace