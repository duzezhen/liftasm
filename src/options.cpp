#include "../include/options.hpp"
#include <algorithm>

namespace opt {

void init_opts(uint16_t k, uint16_t w, double sec_pri_ratio, int sec_pri_num, bool outPAF, uint16_t t, ChainOpts& c, AnchorOpts& a, ExtendOpts& e, AlignOpts& al) {
    al.k = k;
    al.w = w;
    al.threads = t;
    al.buffer_size = 10;
    al.queue_cap = 1e3;
    al.sec_pri_ratio = sec_pri_ratio;
    al.sec_pri_num = sec_pri_num;
    al.outPAF = outPAF;

    c.k = k;
    c.w = w;
    c.threads = t;
    c.mid_occ_frac = 2e-4f;
    c.q_occ_frac = 0.01f;
    c.bw = 500;
    c.max_skip = 25;
    c.max_iter = 500;
    c.max_dist_x = 5e4;
    c.max_dist_y = 5e4;
    c.gap_ext = 2;
    c.anchor_score = 1;
    c.min_n_seeds = 3;
    c.min_score = 3;
    c.use_log_gap = true;
    c.n_seg = 1;
    c.min_cnt = 3;
    c.min_sc = 40;
    c.chn_pen_gap = 0.8f;
    c.chn_pen_skip = 0.0f;
    c.is_cdna = 0; // not cDNA by default

    a.small_slop = 4;
    a.small_gap = 50;
    a.max_kept = 10;
    a.min_group_overlap_frac = 0.01f;

    e.k = k;
    e.w = w;
    e.threads = t;
    e.flank_pad = 30;
    e.match = 0;
    e.mismatch = 4;
    e.gap_open = 4;
    e.gap_extend = 2;
    e.min_align_score = c.min_score * e.match;
    e.min_qry_local_identity = 0.90f;
    e.min_ref_local_identity = 0.90f;
    e.min_qry_global_identity = 0.90f;
    e.hard_clip = false;
    e.dyn_rescue_enable = true;
    e.dyn_step_bp = 1e3;
    e.dyn_zdrop = -1e3;
    e.dyn_max_query_bp = 2e4;
    e.dyn_max_ref_bp = 2e4;
    e.dyn_min_gain = 15;
    e.dyn_max_steps = 20;
}

void Preset::map_ont(ChainOpts& c, AnchorOpts& a, ExtendOpts& e, AlignOpts& al) {
    c.bw = /*500*/10; c.chn_pen_gap = 5; c.gap_ext = 2;

    e.min_qry_local_identity = 0.85f;
    e.min_ref_local_identity = 0.85f;
    e.min_qry_global_identity = 0.85f;
    e.dyn_zdrop = -1000; e.gap_open = 16; e.gap_extend = 1;
    return;
}

void Preset::map_hifi(ChainOpts& c, AnchorOpts& a, ExtendOpts& e, AlignOpts& al) {
    c.bw = /*1000*/10; c.chn_pen_gap = 5; c.gap_ext = 1; c.min_score = 10;
    a.small_slop = 100; a.small_gap = 5e3;

    e.mismatch = 4;
    e.gap_open = 6;
    e.gap_extend = 2;
    e.min_align_score = 200;
    e.dyn_zdrop = -1000;
    return;
}

void Preset::map_illumina(ChainOpts& c, AnchorOpts& a, ExtendOpts& e, AlignOpts& al) {
    c.bw = /*200*/10; c.chn_pen_gap = 4; c.gap_ext = 2;
    e.dyn_zdrop = -100; e.gap_open = 4; e.gap_extend = 2;
    return;
}

void Preset::map_asm_5(ChainOpts& c, AnchorOpts& a, ExtendOpts& e, AlignOpts& al) {
    c.bw          = 2e4;
    c.max_dist_x  = 3e4;
    c.max_dist_y  = 3e4;
    c.max_iter    = 1e5;
    c.max_skip    = 500;

    c.min_cnt     = 3;
    c.min_sc      = 10;
    c.chn_pen_gap = 3.0f;
    c.chn_pen_skip= 0.0f;
    c.mid_occ_frac = 1e-3f;

    c.is_cdna = 0; c.n_seg = 1;

    a.small_gap  = 1.5e4;
    a.small_slop = 1.5e4;
    a.max_kept   = 30;
    a.min_group_overlap_frac = 0.0f;

    e.match      = 0;
    e.mismatch   = 4;
    e.gap_open   = 6;
    e.gap_extend = 1;

    e.flank_pad  = 200;
    e.dyn_rescue_enable   = true;
    e.dyn_step_bp         = 2e3;
    e.dyn_max_query_bp    = 1.5e4;
    e.dyn_max_ref_bp      = 1.5e4;
    e.dyn_max_steps       = 250;
    e.dyn_zdrop           = -1e3;
    e.min_qry_local_identity  = 0.80f;
    e.min_ref_local_identity  = 0.80f;
    e.min_qry_global_identity = 0.80f;

    // e.max_wfa_bp  = 200'000;   // rseg+qseg > 200 kb 不跑 WFA，直接用 =/I/D 退化
    // e.max_gap_wfa = 150'000;   // 单向 gap > 150 kb 直接 I/D，不跑 WFA

    al.sec_pri_num = 1;
}

void Preset::map_other(ChainOpts& c, AnchorOpts& a, ExtendOpts& e, AlignOpts& al) {
    c.bw = 500; c.chn_pen_gap = 5; c.gap_ext = 2;
    e.dyn_zdrop = -1000; e.gap_open = 8; e.gap_extend = 1;
    return;
}

} // namespace opt