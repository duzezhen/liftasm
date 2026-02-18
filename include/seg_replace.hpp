#pragma once
// seg_replace.hpp — DAG-based replacement expansion with memoization.
//
// Layout (current):
//   [127..65] chrId  (63 bits)
//   [64 ..33] beg    (32 bits)
//   [32 .. 1] end    (32 bits)
//   [0]       strand (1 bit)    0:'+' 1:'-'
//
// Features
// ────────────────────────────────────────────────────────────────
// • Interval: pack / unpack / toggle / format / parse helpers
// • Expander: memoized DFS with cycle detection
// • build_index(): pre-compute full expansion for all segments
// • query():  prefer index, auto-fill on miss (heavy-query friendly)
// • Printing helpers for quick inspection
//

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "logger.hpp"

namespace SegReplace {

// ────────────────────────────────────────────────────────────────
// Interval — 128-bit encoded genomic segment helpers
// ────────────────────────────────────────────────────────────────
struct Interval {
    using u128 = unsigned __int128;

    // Encode and decode
    static u128                 pack(uint64_t chrId, uint32_t beg, uint32_t end, bool strand);
    static uint64_t             seg_id(u128 v);
    static uint32_t             beg(u128 v);
    static uint32_t             end(u128 v);
    static bool                 is_reverse(u128 v);
    static uint32_t             len(u128 v);
    static u128                 toggle_strand(u128 v);                      // Toggle strand
    static std::vector<u128>    toggle_strand_chain(std::vector<u128> vs);  // Toggle strand for a chain (sample to reverse_and_toggle)
    static u128                 canonical(u128 v);                          // Return forward-strand
    static std::string          to_string(u128 x);                          // Convert to decimal string (ostream lacks __int128 support)
    static std::string          format(u128 v, const std::vector<std::string>& names);  // format_interval
    static std::string          format(u128 v);  // format_interval
    static std::string          name(u128 v);
};

// ────────────────────────────────────────────────────────────────
// Hash / equality helpers for u128
// ────────────────────────────────────────────────────────────────
struct U128Hash {
    std::size_t operator()(Interval::u128 x) const noexcept {
        std::uint64_t hi = static_cast<std::uint64_t>(x >> 64);
        std::uint64_t lo = static_cast<std::uint64_t>(x);
        // Simple 64-bit mix
        return std::hash<std::uint64_t>()(
            hi ^ (lo + 0x9e3779b97f4a7c15ULL + (hi << 6) + (hi >> 2)));
    }
};
struct U128Eq {
    bool operator()(Interval::u128 a, Interval::u128 b) const noexcept { return a == b; }
};

// Type aliases for readability
using Seg       = Interval::u128;
using Expansion = std::vector<Seg>;
using RuleMap   = std::unordered_map<Seg, Expansion, U128Hash, U128Eq>;

// ────────────────────────────────────────────────────────────────
// Expander — memoized DFS with cycle detection
// ────────────────────────────────────────────────────────────────
class Expander {
public:
    Expander(const RuleMap& r, std::vector<std::string> names = {}) : rules_(&r), names_(std::move(names)) {}

    static void reverse_and_toggle(Expansion& v);  // Reverse order & flip strand for reverse-requested expansions.
    void build_index();  // Pre-compute expansion for all segments in rules (forward-strand)
    void remove_from_index_by_key(const Seg& s);  // Remove s and its reverse from memo

    const RuleMap& index_view() const { return memo_; }  // Read-only view of current memo
    void invalidate_all() { memo_.clear(); in_stack_.clear(); }  // Invalidate entire cache (call after rules mutate).

    Expansion query(Seg s) const;  // query

    void print_expansion_direct(const Seg& seg, const Expansion& v) const;  // Print “seg => seg1 seg2 …” (already expanded)
    void print_expansion(const Seg& seg);  // print expansion
    void print_index() const;  // sort index and print

    bool save_map(const std::string& path);

private:
    const std::vector<std::string> names_;
    const RuleMap* rules_;
    RuleMap memo_;
    std::unordered_set<Seg, U128Hash, U128Eq>            in_stack_;

    // Expand a segment
    Expansion expand_(Seg s);
};

} // namespace SegReplace