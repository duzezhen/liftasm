#pragma once

#include "options.hpp"

#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>

struct GlobalOpts {
    /**
     * @brief To reduce memort usage, change k and w from 15 and 10 to 25 and 30
     * @date 2026-04-20
     * @version 0.1.3
     */
    int kmerLen    = 25;
    int minimizerW = 30;

    int threads   = 4;
    int IOthreads = 4;
    bool debug    = false;
};

struct StatOpts {
    std::vector<std::string> gfaFiles;  // input GFA
};

struct SeqOpts {
    std::vector<std::string> gfaFiles;  // input GFA
    std::string pathFile;       // input path file
    std::string outFile = "";   // output extracted sequences
};

struct Gfa2FaOpts {
    std::vector<std::string> gfaFiles;  // input
    std::string outFile = "";  // output FASTA file for gfa2fa
    int   min_len_xbp = 0;     // Extend only when segment length < x; x=0 disables
    int   extend_ybp  = 50;    // Extend y bp on both ends
    int   wrap_width  = 60;    // FASTA line width; 0 means no wrap
    bool  skip_unknown = true; // Skip segments that are "*" or length=0
};

struct DepthOpts {
    std::vector<std::string> gfaFiles;  // input GFA
    std::vector<std::string> reads;  // input read files (fastq/fasta, gz)
    std::string gafFile;             // input GAF alignment file
    std::string outFile = "";
    bool   base_depth = false; // whether to compute base-level depth instead of segment-level depth
    int    min_mapq = 20;      // minimum mapping quality to keep
    double min_frac = 0.7;     // minimum aligned-length / read-length ratio to keep
};

struct BubbleOpts {
    std::vector<std::string> gfaFiles;  // input
    std::vector<std::string> gfaNames;  // will be added to S-line as "SN:Z:" tag in output stats file
    std::string out_prefix = "out";     // output prefix

    // VCF output
    bool write_vcf = false;
    std::string ref_file;               // reference FASTA for REF sequence
    std::string paf_file;               // PAF used to liftover contig coordinates to reference
    uint32_t ali_min_len  = 50000;      // minimum alignment length to keep in PAF (liftover)
    uint32_t ali_min_mapq = 5;          // minimum mapping quality to keep in PAF (liftover)

    uint32_t max_depth = 1e5; // maximum depth for bubble detection
    uint16_t max_paths = 20;  // maximum number of paths to enumerate per bubble
    uint32_t min_len   = 0;   // minimum sequence length for bubble detection
    uint32_t min_num   = 0;   // minimum number of nodes in a bubble

    uint64_t DFS_guard = 1e6;  // DFS state guard

    // complex-source detection (is_complex_source_)
    uint16_t cx_branch_degree = 3;    // minimum degree counted as a branch
    uint16_t cx_hub_degree    = 100;  // minimum degree counted as a hub
    uint32_t cx_min_nodes     = 200;  // minimum nodes in a large complex block
    uint32_t cx_min_branches  = 24;   // minimum irreducible branches in a complex block

    double path_diff = 0.05;        // sequence difference threshold for path clustering

    uint32_t stall_round_limit = 2; // stop DFS path search after this many rounds without a new unique path

    bool keep_nested = false;       // whether to keep bubbles contained within larger bubbles
    bool check_complex = true;      // detect and exclude complex graph regions

    // homologous path options
    double   same_sim        = 0.8;   // minimum similarity for homologous paths from the same source
    uint64_t same_min_len    = 1e6;   // minimum total length of each same-source homologous path (default: 1Mb)

    uint32_t diff_min_src    = 5e5;   // minimum source length for searching homologous paths between different sources (default: 500 kb)
    double   diff_sim        = 0.9;   // minimum similarity for homologous paths between different sources
    uint64_t diff_min_len    = 3e6;   // minimum total sequence length in "different source" homologous paths (default: 3Mb)

    uint32_t homo_num        = 1;     // max paths per source branch
    uint32_t homo_k          = 17;    // k-mer size for minimizer used in homologous path search
    uint32_t homo_w          = 100;   // window size for minimizer used in homologous path search

    uint32_t homo_extend_bp  = 1e6;         // bp extended per homologous-path extension round (default: 1Mb)
    uint64_t homo_bloom_bits = 1ULL << 27;  // Bloom filter size in bits for homologous-path sketches
    uint32_t homo_bloom_hash = 4;           // number of Bloom filter hash functions
};

struct CollapseOpts {
    std::vector<std::string> gfaFiles;  // input GFA
    std::vector<std::string> gfaNames;  // will be added to S-line as "SN:Z:" tag in output stats file
    std::string prefix = "out";         // output file prefix
    // uint32_t iterations = 5;            // number of collapse rounds
    uint32_t iterations = 3;            // number of collapse rounds
    uint32_t long_node_len = 1e7;       // split nodes longer than this; 0 disables (default: 10Mb)

    std::vector<double>   min_jaccards     = {0.8};  // collapse
    std::vector<double>   same_sims        = {0.9};
    std::vector<uint32_t> same_min_lens    = {1000000};
    std::vector<uint32_t> diff_min_srcs    = {500000};
    std::vector<double>   diff_sims        = {0.9};
    std::vector<uint32_t> diff_min_lens    = {3000000};  // (default: 3Mb)
    // std::vector<int>      min_eqs          = {1000, 100, 50, 3, 3};  // collapse defaults
    std::vector<int>      min_eqs          = {1000, 3, 3};  // collapse defaults
    int                   min_eq           = 3;  // deoverlap default

    uint16_t max_iters = 100;  // max iterations for cut propagation

    uint32_t repeat_mask_min_len      = 24;  // minimum tandem-repeat span to mask
    uint32_t repeat_mask_max_period   = 6;   // maximum tandem-repeat period to test
    uint32_t repeat_mask_max_mismatch = 2;   // maximum mismatches allowed in a repeat window
    uint32_t repeat_norm_len          = 1000;

    // Used to prune abnormal cut points
    uint32_t max_abnormal_cut_len   = 4;
    uint32_t min_abnormal_cut_count = 5;

    uint32_t min_trans_len = 3;  // Minimum interval length to allow transitive replacement expansion (i.e. if A->B and B->C, then A->C is allowed only if the replacement interval is >= this length)

    double min_match_ratio = 0.8;   // minimum CIGAR match ratio, alignment score large than this value will be saved
    double min_ali_ratio = 0.001;   // minimum aligned length ratio to keep (deoverlap default)
    // std::vector<double> min_ali_ratios = {0.05, 0.01, 0.005, 0.001};   // minimum aligned length ratio to keep (collapse default)
    std::vector<double> min_ali_ratios = {0.05, 0.01, 0.001};   // minimum aligned length ratio to keep (collapse default)
    uint8_t min_mapq = 20;          // minimum mapping quality to keep
    uint32_t all_pair_len = 1000;   // align all pairs when every path is shorter than this

    std::string mm2_preset = "asm5";
};

struct File2mapOpts {
    std::vector<std::string> inputFiles;
    std::string outFile          = "";
    bool        paf_primary_only = true;
    uint32_t    min_len          = 50000;  // minimum alignment length to keep
    int         min_mapq         = 5;
};

struct CoordMapOpts {
    int      max_hops       = 8;
    uint32_t max_fanout     = 512;
    uint32_t min_len        = 15;
    double   min_frac       = 0.10;
    uint32_t max_total_hits = 2000;
};

struct LiftoverOpts {
    std::vector<std::string> mapFiles;  // coordinate mapping file
    std::string pafFile;  // input PAF alignment file
    std::string bedFile;  // input BED file (required if not checking)
    std::string referenceFile;  // reference FASTA file (required if checking)
    std::string outFile = "";

    std::string regex;

    // normal liftover
    double min_frac = 0.9;  // for BED: require (h.qend-h.qbeg)/(qend-qbeg) >= min_frac

    // extend liftover
    uint32_t flank_win   = 5000;   // window size for each extension (bp)
    uint32_t max_flank   = 20000;  // max extension on each side (bp)
    uint32_t max_gap     = 20000;  // max allowed difference between target and query gap
    uint16_t max_hit     = 5;      // max extended liftover results to keep

    // PAF
    uint32_t min_mapq  = 5;
    uint32_t min_len   = 10000;

    // check-mode
    bool do_check = false;
    uint32_t win = 500;
    uint32_t step = 500;
    uint32_t max_examples = 20;
};

struct MapqBoostOpts {
    std::vector<std::string> mapFiles;  // coordinate mapping files
    std::string in_bam, out_bam;
    uint32_t batch_size = 20000;
    uint8_t mapq_low = 3, mapq_cap = 60;
    bool name_check = false;

    double sub_ovlp_frac = 0.95;
    double K_mapq = 60.0;
    double K_rs   = 50.0;
    double K_as   = 2000.0;
    double K_ml   = 2000.0;
    double K_nm   = 100.0;
    double W_mapq = 1000.0;
    double W_rs   = 100.0;
    double W_as   = 100.0;
    double W_ml   = 1.0;
    double W_nm   = 1.0;
    double close_as_eps = 0.01;
};

struct MapOpts {
    std::vector<std::string> gfaFiles;  // input
    std::vector<std::string> reads;  // input read files (fastq/fasta, gz)
    std::string outFile = "";        // output (SAM/PAF)
    // alignment scoring / filtering / format
    std::string preset   = "hifi";  // mapping preset: hifi/ont/illumina/other
    double sec_pri_ratio = 0.8;     // minimum allowed ratio of secondary-to-primary score
    int    sec_pri_num   = 5;       // number of secondary alignments to keep
    bool   outPAF        = false;   // alignment format: PAF if true, else SAM

    bool use_wfa   = false; // whether to use WFA for alignment, if not use minimap2

    int  zdrop     = 5e3;
    bool zdrop_set = false;

    opt::ChainOpts  chainOpts;
    opt::AnchorOpts anchorOpts;
    opt::ExtendOpts extendOpts;
    opt::AlignOpts  alignOpts;
};

struct ResSitOpts {
    std::string genome_file;
    std::vector<std::string> enzymes;  // e.g. G^AATTC A^AGCTT G^ANTC
    std::string out_bed;
    bool scan_revcomp = false;
    int threads = 4;
};

struct SplitOpts {
    std::vector<std::string> gfaFiles;
    std::vector<std::string> gfaNames;
    std::string prefix = "out";
};

struct AugmentOpts {
    std::vector<std::string> gfa_files;
    std::vector<std::string> gfa_names;
    std::vector<std::string> vcf_files;
    std::string prefix = "out";
    std::string tag = "VCF";
};

struct CleanOpts {
    std::string utg_file;
    std::vector<std::string> ctg_files;
    std::string prefix = "out";
    uint32_t min_reads = 5;
    double min_purity = 0.8;
    double min_comp_overlap = 0.05;
};

enum class ToolMode {
    stat,       // print graph statistics
    bubble,     // identify bubble in graph
    seq,        // Extract sequences along paths from a GFA
    depth,      // compute depth information from sequencing data for a GFA
    deoverlap,  // de-overlap graph
    collapse,   // collapse unitigs and bubbles (requires no-overlap graph, output by deoverlap)
    align,      // align reads to graph (requires no-overlap graph, output by collapse)
    gfa2fa,     // extract sequences along paths from a GFA
    file2map,   // convert GFA/PAF to map format
    liftover,   // Liftover coordinates using a mapping file
    mapq_boost, // adjust mapping quality
    ressit,     // restrict cuts in genome
    split,      // split graph by connected components
    augment,    // augment a GFA with VCF alleles
    clean,      // remove weak UTG links using contig read components
};

struct AppConfig {
    ToolMode         mode = ToolMode::stat;
    GlobalOpts       global;
    StatOpts         stat;
    SeqOpts          seq;
    Gfa2FaOpts       gfa2fa;
    DepthOpts        depth;
    BubbleOpts       bubble;
    CollapseOpts     collapse;
    CoordMapOpts     coordmap;
    File2mapOpts     file2map;
    LiftoverOpts     liftover;
    MapqBoostOpts    homq;
    MapOpts          map;
    ResSitOpts       ressit;
    SplitOpts        split;
    AugmentOpts      augment;
    CleanOpts        clean;
};

void help(char** argv, bool advanced=false);

// statistics about a GFA file
AppConfig main_stat(int argc, char** argv);
void help_stat(char** argv);

// Extract sequences along paths from a GFA
AppConfig main_seq(int argc, char** argv);
void help_seq(char** argv);

// extract sequences along paths from a GFA
AppConfig main_gfa2fa(int argc, char** argv);
void help_gfa2fa(char** argv);

// compute depth information from sequencing data for a GFA
AppConfig main_depth(int argc, char** argv);
void help_depth(char** argv);

// identify bubble in graph
AppConfig main_bubble(int argc, char** argv);
void help_bubble(char** argv, bool advanced=false);

// de-overlap graph
AppConfig main_deoverlap(int argc, char** argv);
void help_deoverlap(char** argv, bool advanced = false);

// collapse unitigs and bubbles (requires no-overlap graph, output by deoverlap)
AppConfig main_collapse(int argc, char** argv);
void help_collapse(char** argv, bool advanced=false);

// convert GFA/PAF to map format
AppConfig main_file2map(int argc, char** argv);
void help_file2map(char** argv);

// Liftover coordinates using a mapping file
AppConfig main_liftover(int argc, char** argv);
void help_liftover(char** argv, bool advanced=false);

// adjust mapping quality
AppConfig main_mapq_boost(int argc, char** argv);
void help_mapq_boost(char** argv, bool advanced=false);

// align reads to graph (requires no-overlap graph, output by collapse)
AppConfig main_align(int argc, char** argv);
void help_align(char** argv);

// restrict cuts
AppConfig main_res_cut(int argc, char** argv);
void help_res_cut(char** argv);

// split graph by connected components
AppConfig main_split(int argc, char** argv);
void help_split(char** argv);

// augment a GFA with VCF alleles
AppConfig main_augment(int argc, char** argv);
void help_augment(char** argv, bool advanced=false);

// remove weak UTG links using contig read components
AppConfig main_clean(int argc, char** argv);
void help_clean(char** argv);

void finalize_opt_cfg(AppConfig& cfg);
