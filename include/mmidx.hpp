#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>
#include <utility>
#include <cmath>

#include "mmidx_types.hpp"
#include "seq_utils.hpp"
#include "options.hpp"
#include "logger.hpp"

class GfaGraph;

/*
 * 
 * Usage:
 *   std::vector<std::pair<std::string,std::string_view>> name_seqs = vector<sequence_name, sequence>
 *   mmidx::MinimizerIndex idx(name_seqs, k, w, threads);
 *   idx.build();
 *   std::vector<mmidx::MinimizerIndex::Seed> seeds;
 *   idx.collect_seeds(read, seeds);
 *   auto chains = idx.chain_seeds(seeds); // multi-chain per contig
 *   auto anchorChains = idx.build_anchor_blocks(seeds, chains, read.size(), k);
 */

namespace mmidx {

/* ----------------------------- class ----------------------------- */
class MinimizerIndex {
public:
    /* build the minimizer index */
    MinimizerIndex(
        const std::vector<std::string>& names, 
        const std::vector<std::string_view>& seqs, 
        const std::vector<std::vector<std::string>>& right_seqs,
        opt::ChainOpts& chainOpts, 
        const opt::AnchorOpts& anchorOpts
    );

    /* build minimizer table */
    // void build_mm();
    void build_mm(bool expand_right = false);

    /*  High-frequency minimizer filtering */
    uint32_t calc_max_occ(float f) const;

    // ----------------------- seeding -----------------------
    // false = all seeds; true = only seeds on the same strand as the read
    std::vector<MM128> collect_seeds(std::string_view read, bool keep_same_strand_only = false) const;

    // Sort anchors by target coordinate then query coordinate (required by DP)
    void sort_seeds_by_ref(std::vector<MM128>& seeds) const;
    void sort_seeds_by_qry(std::vector<MM128>& seeds) const;

    /* k-mer exact query */
    span_pos query(std::string_view kmer) const;
    span_pos lookup_hash(__uint128_t key) const;
    __uint128_t canonical_key(std::string_view s) const;

    // count k-mer depth from reads
    void count_depth(const std::vector<std::string>& reads);

    // ----------------------- seeding-chaining -----------------------
    void chain_dp(
        const std::vector<MM128>& seeds,
        std::vector<MM128>& b_out,
        std::vector<Sc_Len>& u_out
    ) const;

    /* -------- merge seeds with overlapping coordinates in the chain into a large segment -------- */
    std::vector<AnchorChain> build_anchor_from_chains(
        const std::vector<MM128>& b_out,
        const std::vector<Sc_Len>& u_out,
        const uint32_t& read_len
    ) const;
    std::vector<AnchorChain> merge_anchor_chains(
        const std::vector<AnchorChain>& ac
    ) const;

    std::vector<std::vector<uint32_t>> group_anchor_chains_by_q_overlap(
        const std::vector<AnchorChain>& anchors
    ) const;

    std::vector<std::vector<uint32_t>> filter_anchor_chains(
        std::vector<std::vector<uint32_t>> groups, 
        int sec_pri_num
    ) const;

    void print_seeds(
        const std::string_view read_name,
        uint32_t read_len, 
        const std::vector<MM128>& seeds
    ) const;

    void print_chains(
        const std::string_view read_name,
        uint32_t read_len, 
        const std::vector<MM128>& b_out,
        const std::vector<Sc_Len>& u_out
    ) const;

    void print_anchors(
        const std::string_view read_name,
        uint32_t read_len, 
        const std::vector<AnchorChain>& anchorChains, 
        const std::vector<std::vector<uint32_t>>& groups,
        std::string type=""
    ) const;

    void print_index_stats() const;

    void attach_graph(const GfaGraph* g) { graph_ = g; }
    std::vector<AnchorChain> cross_segment_chaining(const std::vector<AnchorChain>& per_seg_chains) const;
    void chain_dp_graph(const std::vector<MM128>& seeds, std::vector<MM128>& b_out, std::vector<Sc_Len>& u_out) const;

public:
    /* getters */
    std::size_t k() const { return chainOpts_.k; }
    std::size_t w() const { return chainOpts_.w; }
    std::size_t size() const { return keys_.size(); }
    const std::vector<__uint128_t>& keys()      const { return keys_; }
    const std::vector<uint32_t>&    offs()      const { return offs_; }
    const std::vector<pos_t>&       positions() const { return positions_; }
    const std::vector<uint16_t>&    depths()    const { return depths_; }

    const GfaGraph* graph_{nullptr};

private:
    struct Cand    { uint64_t hash; uint32_t off; __uint128_t key; uint8_t dir; };
    struct MiniOnQuery { __uint128_t key; uint32_t q_off; uint8_t dir; };

    /* seq vector (name, seq-view) */
    const std::vector<std::string> names_;
    const std::vector<std::string_view> seqs_;
    const std::vector<std::vector<std::string>> right_seqs_; // right extensions for each sequence
    opt::ChainOpts&  chainOpts_;
    const opt::AnchorOpts& anchorOpts_;

    /* bit masks for rolling k-mer */
    __uint128_t  mask_;
    uint32_t     shift1_;

    /* minimizer tables */
    std::vector<__uint128_t> keys_;       // minimizer keys
    std::vector<uint32_t>    offs_;       // offsets in positions_
    std::vector<pos_t>       positions_;  // sequence ID, offset and direction in the sequence
    std::vector<uint16_t>    depths_;     // k-mer depths

    /* hash table */
    std::vector<uint32_t> ht_;
    uint64_t              ht_mask_ = 0;

private:
    /* internal build helpers */
    uint64_t hash128_(__uint128_t v) const;
    void radix_sort128_(std::vector<MiniRec>& v);
    void build_hash_();

    template<class Emit> void for_each_kmer_(
        std::string_view s,
        std::size_t k,
        __uint128_t mask,
        uint32_t shift1,
        Emit&& emit
    ) const;

    template<class Emit> void for_each_minimizer_(
        std::string_view s,
        std::size_t k,
        std::size_t w,
        __uint128_t mask,
        uint32_t shift1,
        Emit&& emit
    ) const;

    // ----------------------- seeding-chaining -----------------------
    static inline int32_t mm_log2_i32_(int32_t x) { if (x <= 1) return 0; return int32_t(std::log2(double(x))); };
    static inline int32_t mm_comput_sc_(
        const MM128* ai, const MM128* aj,
        int32_t max_dist_x, int32_t max_dist_y, int32_t bw,
        float chn_pen_gap, float chn_pen_skip,
        int is_cdna, int n_seg
    );
    static inline int64_t mm_chain_bk_end_(
        int32_t max_drop,
        const std::vector<int32_t>& f,
        const std::vector<int64_t>& p,
        std::vector<int32_t>& t,
        const std::vector<std::pair<int32_t,int64_t>>& z, // z[k] = {f[i], i}
        int64_t k
    );
    static inline std::vector<Sc_Len> mm_chain_backtrack_(
        int64_t n,
        const std::vector<int32_t>& f, const std::vector<int64_t>& p,
        std::vector<int32_t>& v, std::vector<int32_t>& t,
        int32_t min_cnt, int32_t min_sc, int32_t max_drop,
        int32_t& n_u, int32_t& n_v
    );
    static inline void mm_compact_sort_(
        const std::vector<MM128>& a,
        const std::vector<int32_t>& v,
        const std::vector<Sc_Len>& u,
        std::vector<MM128>& b_out,
        std::vector<Sc_Len>& u_sorted
    );

    // ----------------------- Merge anchor with overlapping coordinates -----------------------
    void anchors_coordinate_merge_(std::vector<Anchor>& anchors) const;


    int32_t score_transition_seg_(const SegNode& A, const SegNode& B) const;
    inline bool is_connected_vertex_(uint32_t v_from, uint32_t v_to) const;
    // Graph-aware seed scoring: same segment uses mm_comput_sc_, cross-segment uses graph API
    int32_t mm_comput_sc_graph_(const MM128* ai, const MM128* aj) const;
};

} // namespace mmidx