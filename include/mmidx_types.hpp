#pragma once

namespace mmidx {

// ----------------------- Build -------------------
// [63:32] = seq_id, [31:1] = offset, [0] = direction (0=fwd,1=rev complement)
using pos_t = uint64_t;
inline uint32_t seq_id(pos_t p)  { return uint32_t(p>>32); }
inline uint32_t offset(pos_t p)  { return uint32_t((p>>1)&0x7fffffffU); }
inline uint8_t  is_rev(pos_t p)  { return uint8_t(p & 1); }

struct MiniRec { __uint128_t key; pos_t pos; };  // pos: [63:32] = seq_id, [31:1] = offset, [0] = direction
/* used by query() */
struct span_pos {
    const pos_t* b;
    const pos_t* e;
    const pos_t* begin() const { return b; }
    const pos_t* end()   const { return e; }
    size_t       size()  const { return size_t(e-b); }
    const pos_t& operator[](size_t i) const { return b[i]; }
    operator bool() const { return b!=e; }
};


// ----------------------- Seeding -------------------
struct MM128 {
    // The dir here refers to whether ref and query are in the same direction
    // x: [63]=dir (0=fwd,1=rev), [62..32]=r_id, [31..0]=r_off
    // y: [63..48]=flags_hi, [47..40]=seg_id, [39..32]=q_span, [31..0]=q_off
    // q_off is the query offset in the forward or reverse strand according to dir
    uint64_t x, y;

    // --- get helpers ---
    inline uint8_t  dir()   const { return uint8_t(x >> 63); }
    inline uint32_t r_id()  const { return uint32_t((x >> 32) & 0x7FFFFFFF); }
    inline uint32_t r_off() const { return uint32_t(x & 0xFFFFFFFF); }

    inline uint16_t flags_hi() const { return uint16_t(y >> 48); }
    inline uint8_t  seg_id()   const { return uint8_t((y >> 40) & 0xFF); }
    inline uint8_t  q_span()   const { return uint8_t((y >> 32) & 0xFF); }
    inline uint32_t q_off()    const { return uint32_t(y & 0xFFFFFFFF); }

    // --- pack helpers ---
    static inline uint64_t pack_x(uint8_t dir, uint32_t r_id, uint32_t r_off) {
        return (uint64_t(dir & 1) << 63) | (uint64_t(r_id & 0x7FFFFFFF) << 32) |  uint64_t(r_off);
    }
    static inline uint64_t pack_y(uint16_t flags_hi, uint8_t seg_id, uint8_t q_span, uint32_t q_off) {
        return (uint64_t(flags_hi) << 48) | (uint64_t(seg_id) << 40) | (uint64_t(q_span) << 32) | uint64_t(q_off);
    }
    static inline MM128 make(uint8_t dir, uint32_t r_id, uint32_t r_off, uint16_t flags_hi, uint8_t seg_id, uint8_t q_span, uint32_t q_off) {
        return MM128{ pack_x(dir, r_id, r_off), pack_y(flags_hi, seg_id, q_span, q_off) };
    }
};  // 128-bit minimizer key


// ----------------------- Chaining -------------------
struct Sc_Len {
    uint64_t sc_len; // [63:32] = score, [31:0] = length
    inline int32_t score() const { return int32_t(sc_len >> 32); }
    inline uint32_t len() const { return uint32_t(sc_len); }
};


// ----------------------- Anchors -------------------
/* merged non-overlapping anchor block */
struct Anchor {
    // NOTE:
    // - 0-base, [inclusive, exclusive)
    // - query coordinates are relative to the forward or reverse strand according to chain_rev
    uint32_t r_beg, r_end;
    uint32_t q_beg, q_end;
};
struct AnchorChain {
    std::vector<Anchor> blocks;
    uint32_t            r_id;
    uint32_t            seed_number = 0;   // number of seeds in this anchor chain
    int                 score = 0;         // score of the anchor chain
    bool                chain_rev;
};

struct SegNode {
    uint32_t r_id;     // segment id
    bool     rev;      // Chain direction vs read: true=rev, false=fwd
    uint32_t r_beg, r_end;
    uint32_t q_beg, q_end;
    int32_t  intra_sc; // Chain score within this segment (reuses AnchorChain.score)
    uint32_t seeds;    // Seed count within this segment (for tie-breaks)
    // Index for backtracking (position of this SegNode's AnchorChain in input)
    uint32_t src_idx;
};

}  // namespace mmidx