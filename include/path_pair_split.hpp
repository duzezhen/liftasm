#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace path_pair_split {

struct COp {
    uint32_t len = 0;
    char op = 0;
    COp() = default;
    COp(uint32_t l, char c) : len(l), op(c) {}
};

inline bool consumes_a(char op) { return op == 'M' || op == '=' || op == 'X' || op == 'D'; }
inline bool consumes_b(char op) { return op == 'M' || op == '=' || op == 'X' || op == 'I'; }

struct Pos { size_t node = 0; uint32_t off = 0; };
struct Span { Pos beg; Pos end; };

inline Span whole_span(size_t n) { return {{0, 0}, {n, 0}}; }

enum class AnchorKind : uint8_t { Name, Compared };

inline uint8_t anchor_priority(AnchorKind k) {
    return k == AnchorKind::Compared ? 0 : 1;
}

struct IntervalEdge {
    std::string from_name;
    uint32_t from_beg = 0, from_end = 0;
    std::string to_name;
    uint32_t to_beg = 0, to_end = 0;
};

inline std::string edge_key(const IntervalEdge& e) {
    return e.from_name + '\t' + std::to_string(e.from_beg) + '\t' + std::to_string(e.from_end) + '\t' + e.to_name + '\t' + std::to_string(e.to_beg) + '\t' + std::to_string(e.to_end);
}

class ComparedIndex {
public:
    using TargetMap = std::unordered_map<std::string, std::vector<IntervalEdge>>;

    void add(const std::string& a, uint32_t ab, uint32_t ae, const std::string& b, uint32_t bb, uint32_t be)
    {
        add_direct_({a, ab, ae, b, bb, be});
        add_direct_({b, bb, be, a, ab, ae});
    }

    const TargetMap* targets(const std::string& from_name) const {
        auto it = by_pair_.find(from_name);
        if (it == by_pair_.end()) return nullptr;
        return &it->second;
    }

private:
    std::unordered_set<std::string> seen_;
    std::unordered_map<std::string, TargetMap> by_pair_;

    bool add_direct_(const IntervalEdge& e) {
        if (e.from_beg >= e.from_end || e.to_beg >= e.to_end) return false;

        std::string k = edge_key(e);
        if (!seen_.insert(k).second) return false;

        by_pair_[e.from_name][e.to_name].push_back(e);
        return true;
    }
};

template<class VertexT, class LenOf, class RevOf, class OverlapOf>
class PathLayout {
public:
    /**
     * @brief Construct a layout that maps between path coordinates and absolute coordinates.
     *
     * Absolute coordinates are on the concatenated path sequence after removing
     * each node's prefix overlap with its previous node.
     */
    PathLayout(const std::vector<VertexT>& p, LenOf len_of, RevOf rev_of, OverlapOf overlap_of)
        : path_(&p), len_of_(len_of), rev_of_(rev_of), overlap_of_(overlap_of)
    {
        offsets_.reserve(p.size() + 1);
        clips_.reserve(p.size());

        uint64_t cur = 0;
        offsets_.push_back(0);

        for (size_t i = 0; i < p.size(); ++i) {
            const uint32_t L = len_of_(p[i]);
            uint32_t clip = 0;

            if (i > 0) {
                clip = overlap_of_(p[i - 1], p[i]);
                if (clip > L) clip = L;
            }

            clips_.push_back(clip);

            cur += uint64_t(L - clip);
            if (cur > UINT32_MAX) throw std::runtime_error("path length exceeds uint32_t");

            offsets_.push_back(uint32_t(cur));
        }
    }

    uint32_t total_len() const { return offsets_.empty() ? 0 : offsets_.back(); }

    uint32_t clip(size_t i) const {
        return i < clips_.size() ? clips_[i] : 0u;
    }

    uint32_t abs(Pos p) const {
        if (!path_ || path_->empty()) return 0;
        if (p.node >= path_->size()) return total_len();

        const uint32_t c = clip(p.node);
        const uint32_t L = len(p.node);
        const uint32_t off = std::min<uint32_t>(std::max<uint32_t>(p.off, c), L);

        return offsets_[p.node] + (off - c);
    }

    Pos pos(uint32_t x) const {
        if (!path_ || path_->empty()) return {};
        if (x >= total_len()) return {path_->size(), 0};

        auto it = std::upper_bound(offsets_.begin(), offsets_.end(), x);
        size_t i = size_t((it - offsets_.begin()) - 1);

        return {i, clip(i) + (x - offsets_[i])};
    }

    uint32_t len(size_t i) const { return len_of_((*path_)[i]); }

    std::pair<uint32_t,uint32_t> to_forward(size_t i, uint32_t ob, uint32_t oe) const {
        uint32_t L = len(i);
        if (ob > oe || oe > L) return {0, 0};
        if (!rev_of_((*path_)[i])) return {ob, oe};
        return {L - oe, L - ob};
    }

    std::pair<uint32_t,uint32_t> to_oriented(size_t i, uint32_t fb, uint32_t fe) const {
        uint32_t L = len(i);
        if (fb > fe || fe > L) return {0, 0};
        if (!rev_of_((*path_)[i])) return {fb, fe};
        return {L - fe, L - fb};
    }

private:
    const std::vector<VertexT>* path_;
    LenOf len_of_;
    RevOf rev_of_;
    OverlapOf overlap_of_;
    std::vector<uint32_t> offsets_;
    std::vector<uint32_t> clips_;
};

template<class VertexT>
struct SplitTask {
    enum class Type : uint8_t { Gap, EqualAnchor, ComparedAnchor };

    Type type = Type::Gap;
    AnchorKind kind = AnchorKind::Name;
    Span a_span, b_span;
    std::vector<VertexT> sub_a, sub_b;
    bool whole_node = false;

    bool usable_by_node_aligner() const {
        return type == Type::Gap && whole_node && !sub_a.empty() && !sub_b.empty();
    }
};

/**
 * @brief Return whether a span consists only of whole nodes.
 *
 * A span is considered whole-node if:
 *   - beg.off == 0
 *   - end.off == 0
 *   - beg.node <= end.node
 *
 * Such spans can be sliced directly as path sub-vectors.
 *
 * @param s Span in path coordinates.
 * @return true if the span covers complete nodes only.
 */
inline bool span_is_whole_nodes(const Span& s) {
    return s.beg.off == 0 && s.end.off == 0 && s.beg.node <= s.end.node;
}

/**
 * @brief Extract a whole-node subpath from a path.
 *
 * Coordinate system:
 *   - Input span is expressed in path coordinates.
 *   - The span must satisfy span_is_whole_nodes().
 *
 * Example:
 *   path = [A,B,C,D]
 *   span = { {1,0}, {3,0} }
 *   result = [B,C]
 *
 * @param path Input path.
 * @param s    Whole-node span.
 * @return Subpath, or empty vector if the span is not whole-node.
 */
template<class VertexT>
std::vector<VertexT> slice_whole_nodes(const std::vector<VertexT>& path, Span s) {
    if (!span_is_whole_nodes(s)) return {};
    if (s.beg.node > path.size() || s.end.node > path.size()) return {};
    return {path.begin() + s.beg.node, path.begin() + s.end.node};
}

template<class VertexT, class NameOf, class LenOf, class RevOf, class OverlapOf>
void add_cigar(
    ComparedIndex* cmp,
    const std::vector<VertexT>& a,
    const std::vector<VertexT>& b,
    NameOf name_of,
    LenOf len_of,
    RevOf rev_of,
    OverlapOf overlap_of,
    const std::vector<COp>& cigar,
    uint32_t beg_a,
    uint32_t beg_b, 
    uint32_t min_compared_len = 1
) {
    PathLayout<VertexT, LenOf, RevOf, OverlapOf> la(a, len_of, rev_of, overlap_of);
    PathLayout<VertexT, LenOf, RevOf, OverlapOf> lb(b, len_of, rev_of, overlap_of);

    auto add_compared_rect_at = [&](uint32_t sa0, uint32_t sa1, uint32_t sb0, uint32_t sb1) {
        if (!cmp || sa0 >= sa1 || sb0 >= sb1) return;

        if (std::min(sa1 - sa0, sb1 - sb0) < min_compared_len) return;

        Pos xa0 = la.pos(sa0);
        Pos xa1 = la.pos(sa1);
        Pos xb0 = lb.pos(sb0);
        Pos xb1 = lb.pos(sb1);

        if (xa0.node >= a.size() || xb0.node >= b.size()) return;

        if (xa1.node == a.size()) xa1 = {a.size() - 1, la.len(a.size() - 1)};
        if (xb1.node == b.size()) xb1 = {b.size() - 1, lb.len(b.size() - 1)};

        if (xa0.node != xa1.node || xb0.node != xb1.node) return;

        auto fa = la.to_forward(xa0.node, xa0.off, xa1.off);
        auto fb = lb.to_forward(xb0.node, xb0.off, xb1.off);

        auto name_a = name_of(a[xa0.node]);
        auto name_b = name_of(b[xb0.node]);

        if (name_a != name_b) {
            cmp->add(name_a, fa.first, fa.second, name_b, fb.first, fb.second);
        }
    };

    uint32_t pa = beg_a;
    uint32_t pb = beg_b;
    const uint32_t aln_beg_a = pa;
    const uint32_t aln_beg_b = pb;

    for (const auto& op : cigar) {
        if (!op.len) continue;

        if (consumes_a(op.op)) pa += op.len;
        if (consumes_b(op.op)) pb += op.len;
    }

    add_compared_rect_at(aln_beg_a, pa, aln_beg_b, pb);
}


template<class VertexT, class NameOf, class LenOf, class RevOf, class OverlapOf>
std::vector<SplitTask<VertexT>> split_path_pair(
    const std::vector<VertexT>& a,
    const std::vector<VertexT>& b,
    NameOf name_of,
    LenOf len_of,
    RevOf rev_of,
    OverlapOf overlap_of,
    ComparedIndex* cmp = nullptr
) {
    if (a.empty() || b.empty()) return {};

    struct Anchor {
        Span a, b;
        AnchorKind kind;
        uint32_t ab = 0, ae = 0;
        uint32_t bb = 0, be = 0;
    };

    std::vector<Anchor> anchors;
    anchors.reserve(a.size() + 8);

    PathLayout<VertexT, LenOf, RevOf, OverlapOf> la(a, len_of, rev_of, overlap_of);
    PathLayout<VertexT, LenOf, RevOf, OverlapOf> lb(b, len_of, rev_of, overlap_of);

    const Span a_span = whole_span(a.size());
    const Span b_span = whole_span(b.size());

    const uint32_t a0 = 0, a1 = la.total_len();
    const uint32_t b0 = 0, b1 = lb.total_len();

    std::vector<std::string> a_names(a.size());
    std::vector<std::string> b_names(b.size());

    std::unordered_map<std::string, std::vector<size_t>> a_by_name;
    std::unordered_map<std::string, std::vector<size_t>> b_by_name;

    a_by_name.reserve(a.size() * 2 + 1);
    b_by_name.reserve(b.size() * 2 + 1);

    for (size_t i = 0; i < a.size(); ++i) {
        a_names[i] = name_of(a[i]);
        a_by_name[a_names[i]].push_back(i);
    }

    for (size_t j = 0; j < b.size(); ++j) {
        b_names[j] = name_of(b[j]);
        b_by_name[b_names[j]].push_back(j);
    }

    auto push_anchor = [&](Span as, Span bs, AnchorKind kind) {
        const uint32_t ab = la.abs(as.beg);
        const uint32_t ae = la.abs(as.end);
        const uint32_t bb = lb.abs(bs.beg);
        const uint32_t be = lb.abs(bs.end);

        if (ae <= ab || be <= bb) return;
        if (kind != AnchorKind::Compared && ae - ab != be - bb) return;

        anchors.push_back({as, bs, kind, ab, ae, bb, be});
    };

    // Add anchors based on node name matches.
    {
        size_t min_j = 0;

        for (size_t i = 0; i < a.size(); ++i) {
            auto it = b_by_name.find(a_names[i]);
            if (it == b_by_name.end()) continue;

            const auto& js = it->second;
            auto jt = std::lower_bound(js.begin(), js.end(), min_j);
            if (jt == js.end()) continue;

            const size_t j = *jt;

            push_anchor(
                {{i, la.clip(i)}, {i + 1, 0}},
                {{j, lb.clip(j)}, {j + 1, 0}},
                AnchorKind::Name
            );

            min_j = j + 1;
        }
    }

    // Add anchors based on the ComparedIndex.
    if (cmp) {
        for (const auto& akv : a_by_name) {
            const auto* targets = cmp->targets(akv.first);
            if (!targets) continue;

            const auto& a_nodes = akv.second;

            for (const auto& tkv : *targets) {
                auto bit = b_by_name.find(tkv.first);
                if (bit == b_by_name.end()) continue;

                const auto& b_nodes = bit->second;

                for (const IntervalEdge& e : tkv.second) {
                    if (e.from_beg >= e.from_end || e.to_beg >= e.to_end) continue;

                    for (size_t i : a_nodes) {
                        if (e.from_end > len_of(a[i])) continue;

                        auto ao = la.to_oriented(i, e.from_beg, e.from_end);
                        ao.first = std::max(ao.first, la.clip(i));
                        if (ao.first >= ao.second) continue;

                        for (size_t j : b_nodes) {
                            if (e.to_end > len_of(b[j])) continue;

                            auto bo = lb.to_oriented(j, e.to_beg, e.to_end);
                            bo.first = std::max(bo.first, lb.clip(j));
                            if (bo.first >= bo.second) continue;

                            push_anchor(
                                {{i, ao.first}, {i, ao.second}},
                                {{j, bo.first}, {j, bo.second}},
                                AnchorKind::Compared
                            );
                        }
                    }
                }
            }
        }
    }

    auto make_task = [&](typename SplitTask<VertexT>::Type type, Span as, Span bs, AnchorKind kind) {
        SplitTask<VertexT> t;
        t.type = type;
        t.kind = kind;
        t.a_span = as;
        t.b_span = bs;
        t.whole_node = span_is_whole_nodes(as) && span_is_whole_nodes(bs);
        t.sub_a = slice_whole_nodes(a, as);
        t.sub_b = slice_whole_nodes(b, bs);
        return t;
    };

    if (anchors.empty()) {
        return {make_task(SplitTask<VertexT>::Type::Gap, a_span, b_span, AnchorKind::Name)};
    }

    std::sort(anchors.begin(), anchors.end(), [](const Anchor& x, const Anchor& y) {
        if (x.ab != y.ab) return x.ab < y.ab;
        if (x.bb != y.bb) return x.bb < y.bb;

        const uint32_t xlen = x.ae - x.ab;
        const uint32_t ylen = y.ae - y.ab;
        if (xlen != ylen) return xlen > ylen;

        return anchor_priority(x.kind) < anchor_priority(y.kind);
    });

    std::vector<SplitTask<VertexT>> out;
    out.reserve(anchors.size() * 2 + 1);

    uint32_t ca = a0, cb = b0;

    for (const auto& an : anchors) {
        if (an.ab < ca || an.bb < cb || an.ae > a1 || an.be > b1) continue;

        if (an.ab > ca && an.bb > cb) {
            out.emplace_back(make_task(
                SplitTask<VertexT>::Type::Gap,
                {la.pos(ca), la.pos(an.ab)},
                {lb.pos(cb), lb.pos(an.bb)},
                AnchorKind::Name
            ));
        }

        auto type = (an.kind == AnchorKind::Compared)
            ? SplitTask<VertexT>::Type::ComparedAnchor
            : SplitTask<VertexT>::Type::EqualAnchor;

        out.emplace_back(make_task(type, an.a, an.b, an.kind));

        ca = an.ae;
        cb = an.be;
    }

    if (ca < a1 && cb < b1) {
        out.emplace_back(make_task(
            SplitTask<VertexT>::Type::Gap,
            {la.pos(ca), la.pos(a1)},
            {lb.pos(cb), lb.pos(b1)},
            AnchorKind::Name
        ));
    }

    return out;
}

template<class VertexT, class NameOf>
static inline std::string join_path_names(
    const std::vector<VertexT>& path,
    NameOf name_of,
    const char* sep = ";"
) {
    std::string s;
    for (const auto& v : path) {
        if (!s.empty()) s += sep;
        s += name_of(v);
    }
    return s;
}

template<class VertexT, class NameOf, class LenOf, class RevOf, class OverlapOf>
static inline std::string join_span_names(
    const std::vector<VertexT>& path,
    const Span& span,
    NameOf name_of,
    LenOf len_of,
    RevOf rev_of,
    OverlapOf overlap_of,
    const char* sep = ";"
) {
    std::string s;
    if (path.empty() || span.beg.node >= path.size()) return s;

    const size_t last = (span.end.off == 0) ? span.end.node : span.end.node + 1;

    for (size_t i = span.beg.node; i < last && i < path.size(); ++i) {
        const uint32_t L = len_of(path[i]);
        const uint32_t clip = (i > 0) ? std::min<uint32_t>(overlap_of(path[i - 1], path[i]), L) : 0u;

        uint32_t beg = (i == span.beg.node) ? span.beg.off : clip;
        uint32_t end = L;

        if (i == span.end.node) {
            end = span.end.off;
        } else if (i + 1 == span.end.node && span.end.off == 0) {
            end = L;
        }

        beg = std::max(beg, clip);

        if (beg >= end || end > L) continue;

        uint32_t fbeg = beg;
        uint32_t fend = end;

        if (rev_of(path[i])) {
            fbeg = L - beg;
            fend = L - end;
        }

        if (!s.empty()) s += sep;
        s += name_of(path[i]);
        s += ":";
        s += std::to_string(fbeg);
        s += "-";
        s += std::to_string(fend);
    }

    return s;
}

template<class VertexT>
static inline const char* task_type_name(const SplitTask<VertexT>& t) {
    if (t.type == SplitTask<VertexT>::Type::Gap)              return "Gap           ";
    if (t.type == SplitTask<VertexT>::Type::EqualAnchor)      return "EqualAnchor   ";
                                                               return "ComparedAnchor";
}

} // namespace path_pair_split
