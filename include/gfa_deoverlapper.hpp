#pragma once
#include <unordered_set>
#include <iterator>
#include <utility>
#include <optional>

#include "gfa_parser.hpp"
#include "seg_replace.hpp"
#include "logger.hpp"
#include "seq_utils.hpp"
#include "options.hpp"
#include "gfa_name.hpp"
#include "aligner.hpp"

// GfaDeoverlapper
// ------------------------------------------------------------------
// Purpose:
//   Transform an "overlap-based" GFA into a "non-overlap" GFA by
//   cutting segments at consistent boundaries and rebuilding edges.
//
// Workflow:
//   1. Prune inconsistent overlaps (forward vs reverse mismatches).
//   2. Initialize cut points for all segments (0, L, and overlap windows).
//   3. Iteratively propagate cut points between connected segments
//      until convergence (or max iterations).
//   4. Build canonical replacement rules (SegReplace).
//   5. Expand and rewire edges based on rules to eliminate overlaps.
//   6. Reconstruct the GFA with non-overlapping segments.
//
// Usage:
//   GfaDeoverlapper deoverlapper;
//   deoverlapper.load_from_GFA(FILE);
//   deoverlapper.deoverlap();
//   deoverlapper.write_to_GFA(FILE);

class GfaDeoverlapper : public GfaGraph {
public:

    struct Cuts {
        std::vector<uint32_t> v;
        std::string           name;
    };

public:

    GfaDeoverlapper(int min_eq, int max_prop_iters) {
        set_forbid_overlap(false);
        MIN_EQ_FOR_CUT_ = min_eq;
        MAX_PROPAGATION_ITERS_ = max_prop_iters;
    }

    void set_opts(
        const opt::ChainOpts&  chain, const opt::AnchorOpts& anchor,
        const opt::ExtendOpts& extend, const opt::AlignOpts&  align, 
        bool use_wfa
    ) {
        chainOpts_ = chain; anchorOpts_ = anchor; extendOpts_ = extend; alignOpts_ = align; use_wfa_ = use_wfa;
    }

    void deoverlap(const std::string& prefix);

protected:  // data will be used in rulemap building
    /* -------------------------------------------- collapse bubbles -------------------------------------------- */
    struct MmWfaHit {
        uint32_t r_beg = 0, r_end = 0;     // alignment start/end on reference (0-based, end-exclusive)
        uint32_t q_beg = 0, q_end = 0;     // alignment start/end on query (0-based, end-exclusive, excluded soft-clip)
        std::string cigar;
    };
    struct BubbleAlignment {
        std::string             name_a;   // segment name of A
        std::string             name_b;   // segment name of B
        // Coordinates are on the forward strand of the segments
        std::uint32_t           beg_a;    // alignment start on A (0-based, inclusive)
        std::uint32_t           end_a;    // alignment end on A (0-based, exclusive)
        std::uint32_t           beg_b;    // alignment start on B (0-based, inclusive), excluded soft-clip
        std::uint32_t           end_b;    // alignment end on B (0-based, exclusive), excluded soft-clip
        uint32_t                v_a;      // [31:1]=segment id, [0]=rev, A vertex id
        uint32_t                v_b;      // [31:1]=segment id, [0]=rev, B vertex id
        std::vector<CIGAR::COp> ops;      // alignment cigar (sega vs segb)

        BubbleAlignment() = default;

        BubbleAlignment(std::string na, std::string nb, std::uint32_t ba, std::uint32_t ea, std::uint32_t bb, std::uint32_t eb, uint32_t va, uint32_t vb, std::vector<CIGAR::COp>&& o)
            : name_a(std::move(na)), name_b(std::move(nb)), beg_a(ba), end_a(ea), beg_b(bb), end_b(eb), v_a(va), v_b(vb), ops(std::move(o)) {}
    };

    std::vector<BubbleAlignment> bubble_aligns_;  // Pairwise alignments between bubble branches

    std::vector<Cuts> cuts_;  // Cutting points for each segment
    SegReplace::RuleMap rulemap_;  // Replacement rules for segments

    int   MIN_EQ_FOR_CUT_ ;          // Threshold: only when the length of '=' in CIGAR >= 100, a cut point is made / rulemap is established
    int   MAX_PROPAGATION_ITERS_;    // Maximum number of propagation iterations
    opt::ChainOpts  chainOpts_;
    opt::AnchorOpts anchorOpts_;
    opt::ExtendOpts extendOpts_;
    opt::AlignOpts  alignOpts_;
    bool            use_wfa_ = false;  // whether to use WFA for alignment (default: use minimap2)

protected:
    // offset -> index
    static inline uint32_t find_index_by_offset(int32_t offset, const std::vector<int32_t>& offsets) {
        auto it = std::lower_bound(offsets.begin(), offsets.end(), offset);
        if (it != offsets.end() && *it == offset) {
            return static_cast<uint32_t>(it - offsets.begin());
        }
        return UINT32_MAX;
    }

    inline uint64_t encode_edge_u64_(uint32_t v, uint32_t w) {
        return (uint64_t(v) << 32) | w;
    };
    inline std::pair<uint32_t, uint32_t> decode_edge_u64_(uint64_t e) {
        uint32_t v = (uint32_t)(e >> 32);
        uint32_t w = (uint32_t)(e & 0xFFFFFFFF);
        return {v, w};
    };

    static inline uint32_t clamp32_(int64_t x, uint32_t L) {
        if (x < 0) return 0u;
        if ((uint64_t)x > (uint64_t)L) return L;
        return (uint32_t)x;
    }

    // Check if an expansion is an identity mapping
    static inline bool is_identity_(const SegReplace::Expansion& e, SegReplace::Seg k) {
        return e.size() == 1 && e[0] == k;
    }

    // touchable: same chr & same direction, and end-to-end
    static inline bool touchable_(SegReplace::Seg a, SegReplace::Seg b) {
        if (SegReplace::Interval::seg_id(a) != SegReplace::Interval::seg_id(b)) return false;
        if (SegReplace::Interval::is_reverse(a) != SegReplace::Interval::is_reverse(b)) return false;
        uint32_t ab = SegReplace::Interval::beg(a), ae = SegReplace::Interval::end(a);
        uint32_t bb = SegReplace::Interval::beg(b), be = SegReplace::Interval::end(b);
        return (ae == bb) || (be == ab);
    }

    // On segment v, take the overlap window [vb, ve) based on the direction
    static inline std::pair<uint32_t,uint32_t> v_overlap_pos_(uint32_t Lseg, uint32_t OV, bool rev) {
        if (!rev) return {Lseg - OV, Lseg};
        else      return {0u, OV};
    }

    // On segment w, take the overlap window [wb, we) based on the direction
    static inline std::pair<uint32_t,uint32_t> w_overlap_pos_(uint32_t Lseg, uint32_t OW, bool rev) {
        if (!rev) return {0u, OW};
        else      return {Lseg - OW, Lseg};
    }

    bool prefer_v_as_leaf_(uint32_t v_seg_id, uint32_t w_seg_id, bool& v_is_leaf, bool by_name=false) const {
        if (v_seg_id == w_seg_id) return false;

        if (!by_name) {
            v_is_leaf = v_seg_id < w_seg_id;
            return true;
        }

        // by name
        if (v_seg_id >= nodes_.size() || w_seg_id >= nodes_.size()) {
            error_stream() << "Invalid segment IDs provided: " 
                   << v_seg_id << ", " << w_seg_id 
                   << ". Segment ID exceeds current node list size (" 
                   << nodes_.size() << ").\n";
            std::exit(1);
        }
        const std::string& vn = nodes_[v_seg_id].name;
        const std::string& wn = nodes_[w_seg_id].name;

        v_is_leaf = vn < wn;
    
        return true;
    }

    bool merge_and_dedup_cuts_(std::vector<uint32_t>& cuts, const std::vector<uint32_t>& new_cuts) {
        if (new_cuts.empty()) return false;

        // Remember original information
        const bool was_rev = cuts.size() < 2 ? false : cuts.front() > cuts.back();
        const size_t old_sz = cuts.size();

        // merge and dedup
        cuts.insert(cuts.end(), new_cuts.begin(), new_cuts.end());
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

        // Restore direction
        if (was_rev) { std::reverse(cuts.begin(), cuts.end()); }

        return cuts.size() != old_sz;
    }

    void emit_node_expansion_edges_(
        const SegReplace::Expansion& node_expansion,
        std::unordered_set<uint64_t>& seen_edges,
        gfaName& namer,
        const std::string& v_name,
        const std::string& w_name,
        bool v_is_rev,
        bool w_is_rev,
        std::optional<std::pair<uint32_t,uint32_t>> v_span = std::nullopt,
        std::optional<std::pair<uint32_t,uint32_t>> w_span = std::nullopt
    );


protected:
    // Generate and filter the "forward and reverse inconsistent overlap" mapping, and mark the edges that need to be removed in arcs_
    void prune_overlaps_();

    // Initialize all segment cut points (0, L, and "window endpoints"), and sort and deduplicate
    void initialize_cuts_();

    // Align v-slice (reference) against w-slice (query) and pack as BubbleAlignment.
    std::vector<MmWfaHit> filter_aligns_(std::vector<MmWfaHit> a);
    std::vector<MmWfaHit> align_wfa_(
        const std::string& v_name,
        const std::string& w_name,
        const std::string& v_seq_slice,
        const std::string& w_seq_slice
    );
    std::vector<MmWfaHit> align_mm2_(
        const std::string& v_name,
        const std::string& w_name,
        const std::string& v_seq_slice,
        const std::string& w_seq_slice
    );
    std::vector<BubbleAlignment> align_and_pack_(
        uint32_t            v_vertex,           // packed vertex id for A
        uint32_t            w_vertex,           // packed vertex id for B
        const std::string&  v_name,             // segment name of A
        const std::string&  w_name,             // segment name of B
        const std::string&  v_seq_slice,        // oriented slice sequence of A (already rev-comp if needed)
        const std::string&  w_seq_slice,        // oriented slice sequence of B (already rev-comp if needed)
        uint32_t            vb, uint32_t ve,    // slice coordinates on A in plus strand: [vb, ve)
        uint32_t            wb, uint32_t we     // slice coordinates on B in plus strand: [wb, we)
    );
    // Align overlaps between edges
    void overlaps_align_();

    // Deduplicate alignments (keep the longest among identical segment pairs)
    void dedup_aligns_();

    // Build cuts_ from alignments
    void build_cuts_from_cigar_();

    // Propagate cut points (align mutual projections within the window by offset)
    void propagate_cuts_();

    // Print/verify expansion index (retain original log behavior)
    static void rulemap_verify(const std::vector<GfaNode>& nodes, const SegReplace::RuleMap& idx);

    // Build replacement rules rulemap (retain original logic and logs)
    void build_rulemap_();
    SegReplace::Expander build_SegReplace_(bool filter_abnormal = true);

    // Expand and replace edges based on rule index (retain original logs and structure updates), and finally rebuild
    void expand_and_rewire_edges_(const SegReplace::Expander& ex);

    // Remove unused nodes
    void remove_unused_nodes_();

    // Print cut points (use after initialize_cuts_ and propagate_cuts_)
    static void print_cuts(const std::vector<Cuts>& cuts);

    void build_rules_from_pair_windows_(
        SegReplace::RuleMap& rulemap,
        // trunk
        uint32_t trunk_seg_id, bool trunk_is_rev,
        const std::vector<uint32_t>& cut_in_overlap_trunk,
        // leaf
        uint32_t leaf_seg_id, bool leaf_is_rev,
        const std::vector<uint32_t>& new_cuts_leaf,
        const std::vector<int32_t>&  new_offsets_leaf,
        // leaf older cuts (in overlap region only)
        const std::vector<uint32_t>& cut_in_overlap_leaf
    ) const;

    // Convert an overlap window v->w and the cuts at both ends into rules and write them into the rulemap.
    // This is essentially a functionalized version of your original large code block.
    // Convention: the overlap of v/w is [vb,ve), [wb,we) (half-open interval), and cuts are already sorted and deduplicated (the function will still sort/unique as needed internally).
    //
    // Parameter description:
    //   v_seg_id, w_seg_id : segment id
    //   v_is_rev, w_is_rev : orientation
    //   v_name, w_name     : segment name (for logging only)
    //   vb,ve, wb,we       : overlap window on both ends of the edge in their respective sequences [beg,end)
    //   v_cuts, w_cuts     : cut point vectors at both ends (will be read and new cut points may be appended during the process, but will not be sorted/unique at the end of this function)
    //   rulemap            : rule table (output)
    void record_rulemap(
        uint32_t v_seg_id, uint32_t w_seg_id,
        bool v_is_rev, bool w_is_rev,
        const std::string& v_name, const std::string& w_name,
        uint32_t vb, uint32_t ve, uint32_t wb, uint32_t we,
        SegReplace::RuleMap& rulemap
    );

    // Refine the "cross-boundary merging → re-expansion" within the chain
    void refine_chain_(SegReplace::Expansion& chain, const SegReplace::Expander& ex) const;
};