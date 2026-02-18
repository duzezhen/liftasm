#pragma once
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <cctype>

#include "logger.hpp"

class gfaName {

public:
    struct Piece {
        std::string root;
        uint64_t    lo{0}, hi{0}; // canonical ends (unordered input normalized)
        bool        rev{false};   // input written as hi-lo → reverse along path
        uint64_t    len{0};       // |hi - lo|
    };
    struct Span { std::string root; uint64_t start{0}, end{0}; };

public:
    explicit gfaName(char out_delim = ';')
    : out_delim_(out_delim == '+' ? '+' : ';') {}

    // Format a half-open slice [beg, end) along the composite path.
    // Each token keeps its own direction (lo->hi or hi->lo).
    // If is_rev is true, traverse the whole slice in reverse path order:
    //   reverse token order and flip each sub-span.
    // utg000657l:1715746-1748768;utg001413l:68019-69628:0-34630 -> utg000657l:1715746-1748768;utg001413l:68019-69627
    std::string format_interval_name(
        const std::string& base,
        uint64_t beg,
        uint64_t end,
        bool is_rev
    ) const {
        if (!(beg < end)) {
            error_stream() << "require beg < end, got [" << beg << "," << end << ")\n";
            std::exit(1);
        }

        // Plain root → no inherent per-token direction.
        if (base.find(':') == std::string::npos) {
            std::ostringstream oss;
            if (!is_rev) oss << base << ":" << beg << "-" << end;   // ascending
            else         oss << base << ":" << end << "-" << beg;   // descending
            return oss.str();
        }

        // Parse pieces and preserve per-token direction.
        auto pieces = parse_composite_with_dir_(base);
        if (pieces.empty()) {
            // Fallback if parse failed
            std::ostringstream oss;
            if (!is_rev) oss << base << ":" << beg << "-" << end;
            else         oss << base << ":" << end << "-" << beg;
            return oss.str();
        }

        // Total path length is the sum of piece lengths (|hi-lo|).
        uint64_t total_len = 0;
        for (const auto& p : pieces) total_len += p.len;
        if (end > total_len) {
            error_stream() << "slice exceeds base length: " + base + " (length=" + std::to_string(total_len) + ")" + " slice=[" + std::to_string(beg) + "," + std::to_string(end) + ")";
            std::exit(1);
        }

        // Build forward-path spans honoring each piece direction.
        std::vector<Span> spans; // path order = base order; span may ascend or descend
        {
            uint64_t acc = 0;
            for (const auto& p : pieces) {
                uint64_t s = std::max<uint64_t>(beg, acc);
                uint64_t e = std::min<uint64_t>(end, acc + p.len);
                if (e > s) {
                    uint64_t off_s = s - acc;  // offsets along piece path
                    uint64_t off_e = e - acc;
                    if (!p.rev) {
                        // piece forward: lo -> hi
                        uint64_t a = p.lo + off_s;
                        uint64_t b = p.lo + off_e;         // a < b
                        spans.push_back(Span{p.root, a, b});
                    } else {
                        // piece reverse: hi -> lo
                        uint64_t a = p.hi - off_s;
                        uint64_t b = p.hi - off_e;         // a > b
                        spans.push_back(Span{p.root, a, b});
                    }
                }
                acc += p.len;
                if (acc >= end) break;
            }
        }

        // Merge adjacent spans if contiguous on the same root in path order.
        std::vector<Span> merged;
        for (const auto& sp : spans) {
            if (!merged.empty() && can_merge_(merged.back(), sp)) {
                auto& m = merged.back();

                // Compute union on normalized coordinates
                uint64_t mL = std::min(m.start, m.end);
                uint64_t mH = std::max(m.start, m.end);
                uint64_t sL = std::min(sp.start, sp.end);
                uint64_t sH = std::max(sp.start, sp.end);
                uint64_t uL = std::min(mL, sL);
                uint64_t uH = std::max(mH, sH);

                // Preserve display direction of the accumulated span
                const bool m_asc = (m.start < m.end);
                if (m_asc) { m.start = uL; m.end = uH; }
                else       { m.start = uH; m.end = uL; }
            } else {
                merged.push_back(sp);
            }
        }

        return render_many_(merged);
    }

    // Merge a chain of names into one composite name.
    // Each name can be composite with mixed per-token directions.
    std::string merge_chain_names(
        const std::vector<std::string>& names
    ) const {
        // Flatten all names into forward-path spans (each span keeps its own direction).
        std::vector<Span> spans;
        for (size_t i = 0; i < names.size(); ++i) {
            auto v = flatten_full_name_with_dir_(names[i]);
            spans.insert(spans.end(), v.begin(), v.end());
        }

        // Merge consecutive spans if they overlap or touch on the same root.
        // The merged span keeps the orientation (start/end ordering) of the first span in the group.
        std::vector<Span> merged;
        for (const auto& sp : spans) {
            if (!merged.empty() && can_merge_(merged.back(), sp)) {
                auto& m = merged.back();

                // Compute union on normalized coordinates
                uint64_t mL = std::min(m.start, m.end);
                uint64_t mH = std::max(m.start, m.end);
                uint64_t sL = std::min(sp.start, sp.end);
                uint64_t sH = std::max(sp.start, sp.end);
                uint64_t uL = std::min(mL, sL);
                uint64_t uH = std::max(mH, sH);

                // Preserve display direction of the accumulated span
                const bool m_asc = (m.start < m.end);
                if (m_asc) { m.start = uL; m.end = uH; }
                else       { m.start = uH; m.end = uL; }
            } else {
                merged.push_back(sp);
            }
        }

        return render_many_(merged);
    }

    std::string force_name_dir(const std::string& name, bool is_rev) const {
        if (name.empty()) return name;

        char delim = detect_delim(name);
        auto toks = split_(name, delim);
        if (toks.empty()) return name;

        std::ostringstream oss;

        auto emit_one = [&](const std::string& tk) {
            std::string root; uint64_t lo = 0, hi = 0; bool rev_token = false;
            if (!parse_token_dir_(tk, root, lo, hi, rev_token)) {
                oss << tk;
                return;
            }
            if (!is_rev) {
                if (!rev_token) oss << root << ":" << lo << "-" << hi;
                else            oss << root << ":" << hi << "-" << lo;
            } else {
                if (!rev_token) oss << root << ":" << hi << "-" << lo;
                else            oss << root << ":" << lo << "-" << hi;
            }
        };

        if (!is_rev) {
            for (size_t i = 0; i < toks.size(); ++i) {
                if (i) oss << out_delim_;
                emit_one(toks[i]);
            }
        } else {
            for (size_t k = toks.size(); k-- > 0; ) {
                if (k + 1 != toks.size()) oss << out_delim_;
                emit_one(toks[k]);
            }
        }

        return oss.str();
    }

    // Detect delimiter between composite name tokens (';' or '+').
    static char detect_delim(const std::string& s) {
        if (s.find(';') != std::string::npos) return ';';
        if (s.find('+') != std::string::npos) return '+';
        return ';';
    }

    // Collect root tokens from a (possibly composite) segment name.
    // For tokens like "root:lo-hi", keep only "root".
    static void collect_roots_from_name(const std::string& name, std::vector<std::string>& out) {
        out.clear();
        if (name.empty()) return;

        const bool has_delim = (name.find(';') != std::string::npos) || (name.find('+') != std::string::npos);
        if (!has_delim) {
            size_t p = name.find(':');
            out.emplace_back(p == std::string::npos ? name : name.substr(0, p));
            return;
        }
        const char d = detect_delim(name);
        size_t i = 0;
        while (i < name.size()) {
            size_t j = name.find(d, i);
            if (j == std::string::npos) j = name.size();
            if (j > i) {
                const std::string tk = name.substr(i, j - i);
                size_t p = tk.find(':');
                out.emplace_back(p == std::string::npos ? tk : tk.substr(0, p));
            }
            i = j + 1;
        }
    }

    // Return true if A ∩ B is non-empty.
    static bool roots_intersect(
        const std::vector<std::string>& A,
        const std::vector<std::string>& B
    ) {
        for (const auto& x : A)
            for (const auto& y : B)
                if (x == y) return true;
        return false;
    }

    // utg000657l:1715746-1748768;utg001413l:68019-69627 -> [utg000657l:1715746-1748768, utg001413l:68019-69627]
    std::vector<Piece> parse_composite_with_dir(const std::string& base) const {
        return parse_composite_with_dir_(base);
    }

private:
    char out_delim_{';'};

    // --- string helpers ---
    static inline void trim_(std::string& x) {
        size_t i = 0, j = x.size();
        while (i < j && std::isspace(static_cast<unsigned char>(x[i]))) ++i;
        while (j > i && std::isspace(static_cast<unsigned char>(x[j-1]))) --j;
        x.assign(x, i, j - i);
    }
    static inline std::vector<std::string> split_(const std::string& base, char delim) {
        std::vector<std::string> toks;
        size_t i = 0;
        while (i < base.size()) {
            size_t j = base.find(delim, i);
            if (j == std::string::npos) j = base.size();
            std::string tk = base.substr(i, j - i);
            trim_(tk);
            if (!tk.empty()) toks.emplace_back(std::move(tk));
            i = j + 1;
        }
        return toks;
    }

    // Parse "root:lo-hi" or "root:hi-lo"; return lo,hi plus rev flag.
    static inline bool parse_token_dir_(
        const std::string& s,
        std::string& root, uint64_t& lo, uint64_t& hi, bool& rev
    ) {
        size_t c = s.find(':'), d = (c==std::string::npos) ? std::string::npos : s.find('-', c+1);
        if (c == std::string::npos || d == std::string::npos || d <= c+1) return false;
        root = s.substr(0, c);
        std::string a = s.substr(c+1, d-(c+1));
        std::string b = s.substr(d+1);
        trim_(root); trim_(a); trim_(b);
        if (root.empty() || a.empty() || b.empty()) return false;
        uint64_t x=0,y=0;
        try { x = std::stoull(a); y = std::stoull(b); } catch (...) { return false; }
        rev = (x > y);                // written as hi-lo → reverse along path
        lo  = std::min(x,y);
        hi  = std::max(x,y);
        return true;
    }

    // Parse composite with per-token direction preserved.
    // utg000657l:1715746-1748768;utg001413l:68019-69627 -> [utg000657l:1715746-1748768, utg001413l:68019-69627]
    std::vector<Piece> parse_composite_with_dir_(const std::string& base) const {
        if (base.find(':') == std::string::npos) return {};
        // Single token
        if (base.find(';') == std::string::npos && base.find('+') == std::string::npos) {
            std::string r; uint64_t lo=0, hi=0; bool rev=false;
            if (!parse_token_dir_(base, r, lo, hi, rev)) return {};
            return { Piece{r, lo, hi, rev, hi - lo} };
        }
        // Composite
        char delim = detect_delim(base);
        std::vector<Piece> out;
        for (auto& tk : split_(base, delim)) {
            std::string r; uint64_t lo=0, hi=0; bool rev=false;
            if (!parse_token_dir_(tk, r, lo, hi, rev)) {
                error_stream() << "invalid token: " + tk + " in composite name: " + base + "\n";
                std::exit(1);
            }
            out.push_back(Piece{std::move(r), lo, hi, rev, hi - lo});
        }
        return out;
    }

    static inline std::string render_one_(const std::string& root, uint64_t a, uint64_t b) {
        std::ostringstream oss; oss << root << ":" << a << "-" << b; return oss.str();
    }
    std::string render_many_(const std::vector<Span>& v) const {
        std::ostringstream oss;
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) oss << out_delim_;
            oss << v[i].root << ":" << v[i].start << "-" << v[i].end;
        }
        return oss.str();
    }
    static inline bool can_merge_(const Span& a, const Span& b) {
        if (a.root != b.root) return false;

        uint64_t aL = std::min(a.start, a.end);
        uint64_t aH = std::max(a.start, a.end);
        uint64_t bL = std::min(b.start, b.end);
        uint64_t bH = std::max(b.start, b.end);

        const bool overlap = (std::max(aL, bL) < std::min(aH, bH));
        if (overlap) {
            error_stream()
                << "overlapping spans on root '" << a.root
                << "': [" << a.start << "," << a.end << ") vs ["
                << b.start << "," << b.end << ")\n";
            std::exit(1);
        }

        return a.end == b.start || a.start == b.end;
    }

    // Expand a full name into forward-path spans (respect per-token direction).
    std::vector<Span> flatten_full_name_with_dir_(const std::string& name) const {
        if (name.find(':') == std::string::npos)
            return { Span{name, 0, 0} }; // unknown length, kept as degenerate
        auto pcs = parse_composite_with_dir_(name);
        if (pcs.empty()) return { Span{name, 0, 0} };
        std::vector<Span> v;
        for (const auto& p : pcs) {
            if (!p.rev) v.push_back(Span{p.root, p.lo, p.hi});
            else        v.push_back(Span{p.root, p.hi, p.lo});
        }
        return v;
    }
};