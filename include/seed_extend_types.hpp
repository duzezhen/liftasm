#pragma once
#include <cstdint>
#include <string>

#include "CIGAR.hpp"

namespace seedExtend {

// ------------------------------------------------------------------
// Alignment bounds: 0-based [beg, end) on ref and query
// ------------------------------------------------------------------
struct AlignmentBounds {
    uint32_t r_beg = 0;
    uint32_t r_end = 0;
    uint32_t q_beg = 0;
    uint32_t q_end = 0;
};

// ------------------------------------------------------------------
// Fragments to align around anchor blocks
//   • r_beg, r_end: 0-based [beg, end) on ref
//   • q_beg, q_end: 0-based [beg, end) on query, Note: if chain_rev is true, q_beg and q_end are relative to the reverse-complemented query
//   • kind: F_LEAD, F_GAP, F_TRAIL
//   • rev: true if the chain is reverse-complemented
//   • align_success: true if WFA succeeded for this fragment
// ------------------------------------------------------------------
enum FragKind : uint8_t { F_LEAD, F_GAP, F_TRAIL };
struct FragTask {
    uint32_t r_beg, r_end;   // ref segment
    uint32_t q_beg, q_end;   // query segment (already forward or RC)
    FragKind kind;
    bool     rev           = false;  // chain orientation
    bool     align_success = false;  // did WFA succeed?
};

enum class EndDir { LEAD, TRAIL };
struct DynResResult {
    bool        used = false;
    uint32_t    r_beg = 0, r_end = 0;
    uint32_t    q_beg = 0, q_end = 0;
    std::string cigar;
    int         score = 0;
};

// ------------------------------------------------------------------
// Identity metrics
//   • qry_local_identity = matches / (M+I)
//   • ref_local_identity = matches / (M+D)
//   • qry_global_identity = matches / full_query_length (M+I+S+H)
//   • full_qry_len = M + I + S + H
// ------------------------------------------------------------------
struct IdentityMetrics {
    uint32_t       matches = 0;      // number of matches
    uint32_t       mismatches = 0;   // number of mismatches (approx if no X/=)
    uint32_t       ins_bases = 0;    // total inserted bases
    uint32_t       del_bases = 0;    // total deleted bases
    uint32_t       aln_len = 0;      // M+I+D
    uint32_t       full_qry_len = 0;
    double         qry_local_identity  = 0.0;
    double         ref_local_identity  = 0.0;
    double         qry_global_identity = 0.0;
};

// ------------------------------------------------------------------
// Final alignment of one chain (possibly split into fragments)
// ------------------------------------------------------------------
enum SamFlag : uint16_t {
    F_PAIRED       = 0x0001, // template is paired
    F_PROPER_PAIR  = 0x0002, // read mapped in proper pair
    F_UNMAP        = 0x0004, // read unmapped
    F_MUNMAP       = 0x0008, // mate unmapped
    F_REVERSE      = 0x0010, // read reverse complemented
    F_MREVERSE     = 0x0020, // mate reverse complemented
    F_READ1        = 0x0040, // first in pair
    F_READ2        = 0x0080, // second in pair
    F_SECONDARY    = 0x0100, // secondary alignment
    F_QCFAIL       = 0x0200, // QC fail
    F_DUP          = 0x0400, // optical/PCR duplicate
    F_SUPPLEMENT   = 0x0800  // supplementary alignment (chimeric/split)
};
static inline bool is_primary(uint16_t flag) { return (flag & (F_SECONDARY | F_SUPPLEMENT)) == 0; }
static inline bool is_secondary(uint16_t flag) { return (flag & F_SECONDARY) != 0; }
static inline bool is_supplementary(uint16_t flag) { return (flag & F_SUPPLEMENT) != 0; }
struct samTAGs {
    int NM = 0;          // edit distance
    int AS = 0;          // alignment score
    int cm = 0;          // chain score
    int cs = 0;          // chain seed count
    char tp = 'P';       // primary or secondary alignment
    std::string SA;      // SA:Z tag for supplementary alignments
};
struct FragAlign {
    std::string      qry_name;          // read name
    uint16_t         flag = 0;
    std::string      ref_name;          // contig name
    uint32_t         r_beg=0, r_end=0;  // 0-based [beg, end) on ref
    double           MAPQ = 0.0;        // mapping quality
    std::string      cigar;             // final CIGAR string
    std::string      RNEXT;             // next read name (for paired-end)
    uint32_t         PNEXT = 0;         // next position (1-based)
    uint32_t         TLEN = 0;          // template length (for paired-end)
    std::string      qry_seq;           // query sequence
    std::string      qry_qual;          // query quality
    // --- optional stats for SAM tags ---
    samTAGs          tags;               // SAM tags
    IdentityMetrics  identity_scores;    // per-chain identity metrics

    uint32_t         r_id = 0;           // reference id
    uint32_t         q_beg=0, q_end=0;   // Don't include soft clips
    std::vector<CIGAR::COp> ops;         // CIGAR operations
    bool             align_rev  = false;
    bool             is_primary = false;
};

} // namespace seedExtend