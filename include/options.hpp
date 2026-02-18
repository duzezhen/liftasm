#pragma once
#include <cstdint>
#include <vector>

namespace opt {

// Chain parameters
struct ChainOpts {
    uint16_t k;
    uint16_t w;
    uint16_t threads;

    int bw;             // diagonal bandwidth
    int max_skip;       // MAX_SKIP
    int max_iter;       // hard cap for backscan per anchor
    int max_dist_x;     // break segment if ref gap too large
    int max_dist_y;     // break segment if query gap too large
    
    int gap_ext;        // gap extend slope (used in log space)
    int anchor_score;   // base score per seed
    int min_n_seeds;    // report chain minimum seeds
    int min_score;      // report chain minimum score
    bool use_log_gap;   // use log2 for gap length scaling

    int n_seg;
    int32_t min_cnt;    // minimum number of seeds in a chain
    int32_t min_sc;     // minimum score of a chain
    float chn_pen_gap;  // gap open (rough)
    float chn_pen_skip; // skip penalty

    int is_cdna;        // is cDNA (0: no, 1: yes)

    float mid_occ_frac;
    float q_occ_frac;
	int32_t min_mid_occ, max_mid_occ;
	int32_t mid_occ;     // ignore seeds with occurrences above this threshold
	int32_t max_occ, max_max_occ, occ_dist;
	int64_t mini_batch_size; // size of a batch of query bases to process in parallel
	int64_t max_sw_mat;
};

// Anchor parameters
struct AnchorOpts {
    uint32_t small_slop;  // slop = std::abs(rgap - qgap)
    uint32_t small_gap;   // gap = cur_beg - last_end
    uint16_t max_kept;    // maximum number of anchor groups to keep for alignment
    float min_group_overlap_frac;  // minimum query overlap fraction to consider anchors in the same group
};

/* ------------------------- Extend parameters ------------------------- */
struct ExtendOpts {
    uint16_t k;
    uint16_t w;
    uint16_t threads;

    uint32_t flank_pad;               // Extra sequence (bp) to include beyond chain ends for flank SW
    int      match;                   // Match score for WFA
    int      mismatch;                // Mismatch penalty for WFA
    int      gap_open;                // Gap opening penalty for WFA
    int      gap_extend;              // Gap extension penalty for WFA
    int      min_align_score;         // drop an alignment if the score of the max scoring segment is below this threshold
    bool     hard_clip;               // Hard-clip the read to the alignment bounds

    double min_qry_local_identity;    // minimum identity of the local query segment to be considered for alignment
    double min_ref_local_identity;    // minimum identity of the local reference segment to be considered for alignment
    double min_qry_global_identity;   // minimum identity of the global query segment to be considered for alignment

    bool     dyn_rescue_enable;       // enable dynamic rescue
    uint32_t dyn_step_bp;             // extend step size in base pairs
    int      dyn_zdrop;               // z-drop: stop if current score falls behind best historical score by this threshold
    uint32_t dyn_max_query_bp;        // maximum query extension length
    uint32_t dyn_max_ref_bp;          // maximum reference extension length
    int      dyn_min_gain;            // minimum score gain per step
    int      dyn_max_steps;           // maximum number of steps to try, to prevent excessive runtime
};

struct AlignOpts {
    uint16_t k;
    uint16_t w;
    uint16_t threads;
    uint32_t buffer_size;
    uint32_t queue_cap;

    double sec_pri_ratio;   // secondary alignment to primary alignment ratio
    int    sec_pri_num;     // number of secondary alignments to report
    bool   outPAF;          // output PAF format instead of SAM
};

// Options initialization presets
void init_opts(uint16_t k, uint16_t w, double sec_pri_ratio, int sec_pri_num, bool outPAF, uint16_t t, ChainOpts& chainOpts, AnchorOpts& anchorOpts, ExtendOpts& extendOpts, AlignOpts& alignOpts);

struct Preset {
    static void map_ont(ChainOpts&, AnchorOpts&, ExtendOpts&, AlignOpts&);
    static void map_hifi(ChainOpts&, AnchorOpts&, ExtendOpts&, AlignOpts&);
    static void map_illumina(ChainOpts&, AnchorOpts&, ExtendOpts&, AlignOpts&);
    static void map_asm_5(ChainOpts&, AnchorOpts&, ExtendOpts&, AlignOpts&);
    static void map_other(ChainOpts&, AnchorOpts&, ExtendOpts&, AlignOpts&);
};

} // namespace opt