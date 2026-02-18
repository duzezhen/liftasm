#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <stdexcept>

#include "gfa_parser_AUX.hpp"

// NodeHandle: (segment_id, is_reverse)
struct NodeHandle {
    uint64_t id{0};
    bool is_reverse{false};
    uint64_t get_vertex_id() const            { return (id << 1) | (is_reverse ? 1 : 0); }
    static uint64_t get_segment_id(uint64_t v){ return v >> 1; }
    static bool get_is_reverse(uint64_t v)    { return v & 1ULL; }
    bool operator==(const NodeHandle& o) const{ return id==o.id && is_reverse==o.is_reverse; }
};
/* hash */
namespace std {
template<> struct hash<NodeHandle>{
    size_t operator()(const NodeHandle& h) const noexcept {
        return hash<uint64_t>()(h.id) ^ (hash<bool>()(h.is_reverse) << 1);
    }
};
}

// segment structure (S-line)
struct GfaNode {
    std::string name;
    std::string sequence{"*"};
    uint32_t    length{0};
    bool        deleted{false};
    GfaAux      aux;

    std::vector<uint8_t> coverage;  // coverage per base, only used in GfaBamLoader
};

// path segment structure (P-line)
struct PathSegment { uint64_t node_id{0}; bool is_reverse{false}; };
struct GfaPath     { std::string name; std::vector<PathSegment> segments; };

// alignment structure (A-line)
struct GfaAlignment {
    uint64_t unitig_node_id{0};
    uint32_t position_on_unitig{0};
    bool read_is_rev{false};
    std::string read_name;
    uint32_t read_start_pos{0}, read_end_pos{0};
    std::string optional_tags;
};

/* Arc structure (L-line)
       |<--- lv --->|<-- ov -->|
    v: ------------------------>
                    ||overlap|||
                 w: -------------------------->
                    |<-- ow -->|<---- lw ---->|
*/
struct GfaArc {
    uint64_t v_lv{0};      // [63:33]=segment id, [32]=rev, [31:0]=lv
    uint32_t w{0};         // [31:1]=segment id, [0]=rev, target vertex id
    int32_t  rank{-1};
    int32_t  ov{INT32_MAX}, ow{INT32_MAX};

    /* flags & link id */
    static constexpr uint64_t FLAG_COMP   = 1ULL << 0;
    static constexpr uint64_t FLAG_DEL    = 1ULL << 1;
    static constexpr uint64_t FLAG_STRONG = 1ULL << 2;
    static constexpr uint64_t FLAG_MASK   = 0x7ULL;      // lowest 3 bits
    uint64_t link_id_and_flags{0}; // [63:3]=link id, [2]=strong, [1]=del, [0]=comp

    void set_comp(bool x)   { x ? link_id_and_flags |= FLAG_COMP   : link_id_and_flags &= ~FLAG_COMP; }
    void set_del(bool x)    { x ? link_id_and_flags |= FLAG_DEL    : link_id_and_flags &= ~FLAG_DEL; }
    void set_strong(bool x) { x ? link_id_and_flags |= FLAG_STRONG : link_id_and_flags &= ~FLAG_STRONG; }
    bool get_comp() const   { return link_id_and_flags & FLAG_COMP;   }
    bool get_del() const    { return link_id_and_flags & FLAG_DEL;    }
    bool get_strong() const { return link_id_and_flags & FLAG_STRONG; }

    void     set_link_id(uint64_t id) { link_id_and_flags = (id << 3) | (link_id_and_flags & FLAG_MASK); }
    uint64_t get_link_id() const      { return link_id_and_flags >> 3; }

    uint32_t get_source_vertex_id() const          { return static_cast<uint32_t>(v_lv >> 32); }
    uint32_t get_source_segment_len() const        { return static_cast<uint32_t>(v_lv); }
    void     set_source_segment_len(uint32_t len)  { v_lv = (v_lv & 0xffffffff00000000ULL) | len; }
    uint32_t get_target_vertex_id() const          { return w; }

    uint32_t get_source_segment_id() const         { return get_source_vertex_id() >> 1; }
    uint32_t get_target_segment_id() const         { return get_target_vertex_id() >> 1; }
    bool get_source_is_reverse() const             { return get_source_vertex_id() & 1; }
    bool get_target_is_reverse() const             { return get_target_vertex_id() & 1; }

    bool operator<(const GfaArc& o) const          { return v_lv < o.v_lv; }
};

// For path assembly without P-lines
struct PathSequence {
    std::vector<uint32_t> vertexs;  // vertex-id chain
    std::string sequence;           // assembled sequence using ow rule
};

// For expanded sequences with k-1 bases from outgoing edges
struct ExpandedSeqs {
    std::vector<std::string> names = {};
    std::vector<std::string_view> seqs = {};
    std::vector<std::vector<std::string>> right_seqs = {};    // k-1 bases expanded from the right
};