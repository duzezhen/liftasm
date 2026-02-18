#pragma once

#include <string>
#include <vector>
#include <cstdint>

/* ===================== I/O option ===================== */
struct IOInOpts {
    std::string gfaFile;             // input GFA (original graph)
    std::string pathFile;            // input path file (one path per line)
    std::vector<std::string> reads;  // input read files (fastq/fasta, gz)
    std::string gafFile;             // input GAF alignment file
    std::string inputFile;           // input file (generic)
    std::string pafFile;             // input PAF alignment file (liftover)
    std::string mapFile;             // coordinate mapping file (liftover)
    std::string bedFile;             // input BED file (liftover, required if not checking)
    std::string referenceFile;       // reference FASTA file (liftover, required if checking)
};

struct IOOutOpts {
    std::string gfaDepth = "";       // output depth information to file
    std::string gfaSeq = "";         // output extracted sequences
    std::string gfaDeoverlap = "";   // output GFA after de-overlap
    std::string gfaCollapse = "";    // output GFA after unitig collapse
    std::string gfaBubbles = "";     // output GFA after bubble collapse
    std::string gfafinal = "";       // final graph if pipeline combines steps
    std::string alignOut = "";       // alignments output (SAM/PAF)
    std::string gfa2faOut = "";      // output FASTA file for gfa2fa
    std::string prefix = "out";      // output file prefix
    std::string outFile = "";        // output file name (generic)
};

/* ===================== Global ===================== */
struct GlobalOpts {
    int kmerLen    = 15;
    int minimizerW = 10;

    bool use_wfa   = false; // whether to use WFA for alignment, if not use minimap2

    int threads   = 4;
    int IOthreads = 4;
    bool debug    = false;
};

/* ===================== Bubble options ===================== */
struct BubbleOpts {
    int max_depth = 2e3;    // maximum depth for bubble detection
    int max_paths = 20;     // maximum number of paths to enumerate per bubble
    uint32_t min_len  = 0;  // minimum sequence length for bubble detection
    uint32_t min_num  = 0;  // minimum number of nodes in a bubble

    uint64_t DFS_guard = 1e4;  // DFS state guard

    // complex-source detection (is_complex_source_)
    uint16_t cx_depth       = 6;    // local BFS depth
    uint32_t cx_nodes       = 60;   // visited nodes limit
    uint32_t cx_branches    = 12;   // branching nodes limit
    int      cx_deg_branch  = 3;    // degree >= this counts as branching
    int      cx_deg_hub     = 8;    // degree >= this => hub => complex
    int      cx_deg_cap     = 12;   // stop counting degree after this

    // path clustering (pick_representative_paths)
    double   path_diff = 0.05;  // small than this is same cluster; larger than this is different haplotype
};

/* ===================== Depth options ===================== */
struct DepthOpts {
    int    min_mapq = 20;      // minimum mapping quality to keep
    double min_frac = 0.7;     // minimum aligned-length / read-length ratio to keep
};

/* ===================== Collapse options ===================== */
struct CollapseOpts {
    double min_jaccard    = 0.8;  // Jaccard threshold to align branches
    int    min_eq         = 3;    // minimal '=' length to make cuts / rules
    int    max_iters      = 10;   // max iterations for cut propagation

    double min_new_frac = 0.10;   // align submission dedup in bubbles_align_()
};

/* ===================== Index / chaining / alignment ===================== */
struct MapOpts {
    // alignment scoring / filtering / format
    std::string preset   = "hifi";  // mapping preset: hifi/ont/illumina/other
    double sec_pri_ratio = 0.8;     // minimum allowed ratio of secondary-to-primary score
    int    sec_pri_num   = 5;       // number of secondary alignments to keep
    bool   outPAF        = false;   // alignment format: PAF if true, else SAM
};

struct Gfa2FaOpts {
    int   min_len_xbp = 0;     // Extend only when segment length < x; x=0 disables
    int   extend_ybp  = 50;    // Extend y bp on both ends
    int   wrap_width  = 60;    // FASTA line width; 0 means no wrap
    bool  skip_unknown = true; // Skip segments that are "*" or length=0
};

struct File2mapOpts {
    bool paf_primary_only = true;
    uint32_t min_len = 50000;  // minimum alignment length to keep
    int min_mapq = 5;
};

struct CoordMapOpts {
    int      max_hops       = 15;
    uint32_t max_fanout     = 512;
    uint32_t min_len        = 15;
    double   min_frac       = 0.10;
    uint32_t max_total_hits = 2000;
};

struct LiftoverOpts {
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
    std::string in_bam, out_bam;
    uint32_t batch_size = 20000;
    uint8_t mapq_low = 3, mapq_new = 60;
    double min_frac = 0.95;
    int min_equiv = 1;
    bool name_check = false;
};

/* ===================== Tool mode ===================== */
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
    mapq_boost  // adjust mapping quality
};

/* ===================== Whole configuration ===================== */
struct AppConfig {
    ToolMode       mode = ToolMode::stat;
    IOInOpts       in;
    IOOutOpts      out;
    GlobalOpts     global;
    BubbleOpts     bubble;
    DepthOpts      depth;
    CollapseOpts   collapse;
    MapOpts        map;
    Gfa2FaOpts     gfa2fa;
    CoordMapOpts   coordmap;
    File2mapOpts   file2map;
    LiftoverOpts   liftover;
    MapqBoostOpts  homq;
};

void help(char** argv);

// statistics about a GFA file
AppConfig main_stat(int argc, char** argv);
void help_stat(char** argv);

// identify bubble in graph
AppConfig main_bubble(int argc, char** argv);
void help_bubble(char** argv);

// Extract sequences along paths from a GFA
AppConfig main_seq(int argc, char** argv);
void help_seq(char** argv);

// compute depth information from sequencing data for a GFA
AppConfig main_depth(int argc, char** argv);
void help_depth(char** argv);

// de-overlap graph
AppConfig main_deoverlap(int argc, char** argv);
void help_deoverlap(char** argv);

// collapse unitigs and bubbles (requires no-overlap graph, output by deoverlap)
AppConfig main_collapse(int argc, char** argv);
void help_collapse(char** argv);

// align reads to graph (requires no-overlap graph, output by collapse)
AppConfig main_align(int argc, char** argv);
void help_align(char** argv);

// extract sequences along paths from a GFA
AppConfig main_gfa2fa(int argc, char** argv);
void help_gfa2fa(char** argv);

// convert GFA/PAF to map format
AppConfig main_file2map(int argc, char** argv);
void help_file2map(char** argv);

// Liftover coordinates using a mapping file
AppConfig main_liftover(int argc, char** argv);
void help_liftover(char** argv, const AppConfig& cfg, bool print_all=false);

// adjust mapping quality
AppConfig main_mapq_boost(int argc, char** argv);
void help_mapq_boost(char** argv);