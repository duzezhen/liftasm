#include "../include/OptionParser.hpp"
#include "../include/ProgramMetadata.hpp"
#include "../include/logger.hpp"

#include <getopt.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iomanip>

static inline bool is_flag_(const char* s) { return s && s[0] == '-'; }

static void validate_and_print(AppConfig& cfg) {
    using std::left;
    using std::setw;

    auto die    = [&](const char* msg){ error_stream() << msg << "\n"; std::exit(1); };
    auto ensure = [&](bool ok, const char* msg){ if (!ok) die(msg); };
    auto onoff  = [](bool b){ return b ? "ON" : "OFF"; };

    const int KEYW = 25;

    auto kv = [&](const char* key, const auto& val) {
        std::ostringstream oss;
        oss << left << setw(KEYW) << key << ": " << val;
        log_stream() << oss.str() << "\n";
    };

    if (cfg.global.debug) { cfg.global.threads = 1; set_debug(true); }
    ensure(cfg.global.threads >= 1, "threads must be >= 1");
    if (cfg.global.use_wfa) {
        ensure(cfg.global.kmerLen >= 1 && cfg.global.kmerLen <= 61, "k-mer length must be in [1,61] when using WFA");
    } else {
        ensure(cfg.global.kmerLen > 0 && cfg.global.kmerLen <= 28, "k-mer length must be in [1,28] when using mm2");
    }
    ensure(cfg.global.minimizerW >= 1 && cfg.global.minimizerW < 256, "minimizer window must be in (0,256)");

    const char* mode =
        cfg.mode == ToolMode::stat       ? "stat" :
        cfg.mode == ToolMode::bubble     ? "bubble" :
        cfg.mode == ToolMode::seq        ? "seq" :
        cfg.mode == ToolMode::depth      ? "depth" :
        cfg.mode == ToolMode::deoverlap  ? "deoverlap" :
        cfg.mode == ToolMode::collapse   ? "collapse" :
        cfg.mode == ToolMode::align      ? "align" :
        cfg.mode == ToolMode::gfa2fa     ? "gfa2fa" :
        cfg.mode == ToolMode::file2map   ? "file2map" :
        cfg.mode == ToolMode::liftover   ? "liftover" :
        cfg.mode == ToolMode::mapq_boost ? "mapq_boost" : "unknown";

    kv("Version", program::version);
    kv("Mode", mode);
    kv("Threads", cfg.global.threads);
    kv("Debug", onoff(cfg.global.debug));

    switch (cfg.mode) {
        case ToolMode::stat:
            ensure(!cfg.in.gfaFile.empty(), "--gfa is required");
            break;

        case ToolMode::bubble:
            ensure(!cfg.in.gfaFile.empty(),        "--gfa is required");
            ensure(cfg.bubble.max_depth >= 1,      "--depth must be >= 1");
            ensure(cfg.bubble.max_paths >= 1,      "--paths must be >= 1");
            ensure(cfg.bubble.DFS_guard >= 1,      "--DFS_guard must be >= 1");
            ensure(cfg.bubble.cx_depth >= 1,       "--cx_depth must be >= 1");
            ensure(cfg.bubble.cx_nodes >= 1,       "--cx_nodes must be >= 1");
            ensure(cfg.bubble.cx_branches >= 1,    "--cx_branches must be >= 1");
            ensure(cfg.bubble.cx_deg_branch >= 1,  "--cx_deg_branch must be >= 1");
            ensure(cfg.bubble.cx_deg_hub >= 1,     "--cx_deg_hub must be >= 1");
            ensure(cfg.bubble.cx_deg_cap >= 1,     "--cx_deg_cap must be >= 1");
            ensure(cfg.bubble.path_diff >= 0.0 && cfg.bubble.path_diff <= 1.0, "--path_diff must be in [0,1]");
            ensure(cfg.bubble.min_len >= 0,        "--len must be >= 0");
            ensure(cfg.bubble.min_num >= 0,        "--num must be >= 0");
            kv("Maximum depth",            cfg.bubble.max_depth);
            kv("Maximum paths",            cfg.bubble.max_paths);
            kv("Max DFS states",           cfg.bubble.DFS_guard);
            kv("Complex BFS depth",        cfg.bubble.cx_depth);
            kv("Complex node cap",         cfg.bubble.cx_nodes);
            kv("Complex branch cap",       cfg.bubble.cx_branches);
            kv("Branch degree threshold",  cfg.bubble.cx_deg_branch);
            kv("Hub degree threshold",     cfg.bubble.cx_deg_hub);
            kv("Degree count cap",         cfg.bubble.cx_deg_cap);
            kv("Path diff threshold",      cfg.bubble.path_diff);
            kv("Minimum length",           cfg.bubble.min_len);
            kv("Minimum number of nodes",  cfg.bubble.min_num);
            break;

        case ToolMode::seq:
            ensure(!cfg.in.gfaFile.empty(),  "--gfa is required");
            ensure(!cfg.in.pathFile.empty(), "--path is required");
            kv("GFA input", cfg.in.gfaFile);
            kv("Path file", cfg.in.pathFile);
            if (!cfg.out.gfaSeq.empty()) kv("Output FASTA ", cfg.out.gfaSeq);
            break;

        case ToolMode::depth:
            ensure(!cfg.in.gfaFile.empty(), "--gfa is required");
            ensure(!( !cfg.in.gafFile.empty() && !cfg.in.reads.empty() ), "Options --gaf and --read are mutually exclusive");
            ensure(!( cfg.in.gafFile.empty() && cfg.in.reads.empty() ), "You must provide either --gaf or --read");
            ensure(cfg.depth.min_frac >= 0.0 && cfg.depth.min_frac <= 1.0, "--min_frac must be in [0,1]");
            ensure(cfg.depth.min_mapq >= 0   && cfg.depth.min_mapq <= 60,  "--min_mapq must be in [0,60]");
            kv("K-mer size",        cfg.global.kmerLen);
            kv("Minimizer window",  cfg.global.minimizerW);
            kv("Minimum MAPQ",      cfg.depth.min_mapq);
            kv("Minimum fraction",  cfg.depth.min_frac);
            break;

        case ToolMode::deoverlap:
            ensure(!cfg.in.gfaFile.empty(),     "--gfa is required");
            ensure(!cfg.out.prefix.empty(),     "--prefix is required");
            ensure(cfg.collapse.min_eq    >= 1, "--min_eq must be >= 1");
            ensure(cfg.collapse.max_iters >= 1, "--max_iters must be >= 1");
            ensure(cfg.collapse.min_eq    >= 1, "--min_eq must be >= 1");
            ensure(cfg.collapse.max_iters >= 1, "--max_iters must be >= 1");
            kv("K-mer size",          cfg.global.kmerLen);
            kv("Minimizer window",    cfg.global.minimizerW);
            kv("Alignment method",   (cfg.global.use_wfa ? "WFA" : "mm2"));
            kv("Minimum '=' length",  cfg.collapse.min_eq);
            kv("Max iterations",      cfg.collapse.max_iters);
            kv("Output prefix",       cfg.out.prefix);
            break;

        case ToolMode::collapse:
            ensure(!cfg.in.gfaFile.empty(), "--gfa is required");
            ensure(!cfg.out.prefix.empty(), "--prefix is required");
            ensure(cfg.collapse.min_jaccard >= 0.0 && cfg.collapse.min_jaccard <= 1.0, "--min_jaccard must be in [0,1]");
            ensure(cfg.collapse.min_eq      >= 1, "--min_eq must be >= 1");
            ensure(cfg.collapse.max_iters   >= 1, "--max_iters must be >= 1");
            ensure(cfg.collapse.min_new_frac >= 0.0 && cfg.collapse.min_new_frac <= 1.0, "--min_new_frac must be in [0,1]");
            ensure(cfg.bubble.max_depth     >= 1, "--depth must be >= 1");
            ensure(cfg.bubble.max_paths     >= 1, "--paths must be >= 1");
            ensure(cfg.bubble.DFS_guard     >= 1, "--DFS_guard must be >= 1");
            ensure(cfg.bubble.cx_depth      >= 1, "--cx_depth must be >= 1");
            ensure(cfg.bubble.cx_nodes      >= 1, "--cx_nodes must be >= 1");
            ensure(cfg.bubble.cx_branches   >= 1, "--cx_branches must be >= 1");
            ensure(cfg.bubble.cx_deg_branch >= 1, "--cx_deg_branch must be >= 1");
            ensure(cfg.bubble.cx_deg_hub    >= 1, "--cx_deg_hub must be >= 1");
            ensure(cfg.bubble.cx_deg_cap    >= 1, "--cx_deg_cap must be >= 1");
            ensure(cfg.bubble.path_diff     >= 0.0 && cfg.bubble.path_diff <= 1.0, "--path_diff must be in [0,1]");
            kv("K-mer size",               cfg.global.kmerLen);
            kv("Minimizer window",         cfg.global.minimizerW);
            kv("Alignment method",        (cfg.global.use_wfa ? "WFA" : "mm2"));
            kv("Minimum Jaccard",          cfg.collapse.min_jaccard);
            kv("Minimum '=' length",       cfg.collapse.min_eq);
            kv("Max iterations",           cfg.collapse.max_iters);
            kv("Minimum new fraction",     cfg.collapse.min_new_frac);
            kv("Maximum depth",            cfg.bubble.max_depth);
            kv("Maximum paths",            cfg.bubble.max_paths);
            kv("Max DFS states",           cfg.bubble.DFS_guard);
            kv("Complex BFS depth",        cfg.bubble.cx_depth);
            kv("Complex node cap",         cfg.bubble.cx_nodes);
            kv("Complex branch cap",       cfg.bubble.cx_branches);
            kv("Branch degree threshold",  cfg.bubble.cx_deg_branch);
            kv("Hub degree threshold",     cfg.bubble.cx_deg_hub);
            kv("Degree count cap",         cfg.bubble.cx_deg_cap);
            kv("Path diff threshold",      cfg.bubble.path_diff);
            kv("Output prefix",            cfg.out.prefix);
            break;

        case ToolMode::align:
            ensure(!cfg.in.gfaFile.empty(), "--gfa is required");
            ensure(!cfg.in.reads.empty(),   "--read is required");
            ensure(cfg.map.sec_pri_ratio > 0.0 && cfg.map.sec_pri_ratio <= 1.0, "--secondary must be in (0,1]");
            ensure(cfg.map.sec_pri_num   >= 0, "-N must be >= 0");
            kv("K-mer size",        cfg.global.kmerLen);
            kv("Minimizer window",  cfg.global.minimizerW);
            kv("Preset",            cfg.map.preset);
            kv("Sec/pri ratio",     cfg.map.sec_pri_ratio);
            kv("Sec kept",          cfg.map.sec_pri_num);
            kv("Output",            (cfg.map.outPAF ? "PAF" : "SAM"));
            break;

        case ToolMode::gfa2fa:
            ensure(!cfg.in.gfaFile.empty(),   "--gfa is required");
            ensure(!cfg.out.gfa2faOut.empty(),"--output is required");
            ensure(cfg.gfa2fa.min_len_xbp >= 0, "--min_len must be >= 0");
            ensure(cfg.gfa2fa.extend_ybp  >= 0, "--extend must be >= 0");
            ensure(cfg.gfa2fa.wrap_width  >= 0, "--wrap must be >= 0");
            kv("Min length (X)",   cfg.gfa2fa.min_len_xbp);
            kv("Extend (Y)",       cfg.gfa2fa.extend_ybp);
            kv("Wrap",             cfg.gfa2fa.wrap_width);
            kv("Skip unknown",     (cfg.gfa2fa.skip_unknown ? "true" : "false"));
            break;

        case ToolMode::file2map:
            ensure(!cfg.in.inputFile.empty(), "--input is required");
            ensure(!cfg.out.prefix.empty(),   "--prefix is required");
            ensure(cfg.file2map.min_len >= 0, "--min_len must be >= 0");
            ensure(cfg.file2map.min_mapq >= 0 && cfg.file2map.min_mapq <= 60, "--min_mapq must be in [0,60]");
            kv("Input file", cfg.in.inputFile);
            kv("PAF primary only", (cfg.file2map.paf_primary_only ? "true" : "false"));
            kv("Minimum alignment length", cfg.file2map.min_len);
            kv("Minimum mapping quality", cfg.file2map.min_mapq);
            kv("Output file", cfg.out.prefix + ".map");
            break;

        case ToolMode::liftover:
            ensure(!cfg.in.mapFile.empty(), "-m/--map is required");
            if (cfg.liftover.do_check) { ensure(!cfg.in.referenceFile.empty(), "-r/--ref is required when using --check"); }
            else { ensure(!cfg.in.bedFile.empty(), "-b/--bed is required when not using --check"); }
            ensure(cfg.liftover.win            >= 1, "--win must be >= 1");
            ensure(cfg.liftover.step           >= 1, "--step must be >= 1");
            ensure(cfg.liftover.min_frac       >= 0.0 && cfg.liftover.min_frac <= 1.0, "--min_frac must be in [0,1]");
            ensure(cfg.liftover.max_examples   >= 1, "--max_examples must be >= 1");
            ensure(cfg.liftover.flank_win      >= 1, "--flank_win must be >= 1");
            ensure(cfg.liftover.max_flank      >= 1, "--max_flank must be >= 1");
            ensure(cfg.liftover.max_gap        >= 0, "--max_gap must be >= 0");
            ensure(cfg.liftover.max_hit        >= 1, "--max_hit must be >= 1");
            ensure(cfg.coordmap.max_hops       >= 0, "--cm_max_hops must be >= 0");
            ensure(cfg.coordmap.max_fanout     >= 1, "--cm_max_fanout must be >= 1");
            ensure(cfg.coordmap.min_len        >= 1, "--cm_min_len must be >= 1");
            ensure(cfg.coordmap.min_frac       >= 0.0 && cfg.coordmap.min_frac <= 1.0, "--cm_min_frac must be in [0,1]");
            ensure(cfg.coordmap.max_total_hits >= 1, "--cm_max_hits must be >= 1");
            ensure (cfg.liftover.min_len >= 0, "--min_len must be >= 0");
            ensure (cfg.liftover.min_mapq >= 0 && cfg.liftover.min_mapq <= 60, "--min_mapq must be in [0,60]");
            kv("Coordinate map", cfg.in.mapFile);
            if (cfg.liftover.do_check) { kv("Reference file", cfg.in.referenceFile); kv("Check mode", "ON"); }
            else { kv("BED file", cfg.in.bedFile); kv("Check mode", "OFF"); }
            if (!cfg.in.pafFile.empty()) { kv("PAF file", cfg.in.pafFile); }
            kv("Output file",       cfg.out.outFile);
            kv("Regex filter",     (cfg.liftover.regex.empty() ? "NONE" : cfg.liftover.regex));
            kv("Min fraction",      cfg.liftover.min_frac);
            kv("Flank window",      cfg.liftover.flank_win);
            kv("Flank max",         cfg.liftover.max_flank);
            kv("Gap max",           cfg.liftover.max_gap);
            kv("Max number",        cfg.liftover.max_hit);
            // kv("Window size",       cfg.liftover.win);
            // kv("Step size",         cfg.liftover.step);
            kv("Max hops",          cfg.coordmap.max_hops);
            kv("Max fanout",        cfg.coordmap.max_fanout);
            kv("Min keep length",   cfg.coordmap.min_len);
            kv("Min keep fraction", cfg.coordmap.min_frac);
            kv("Max total hits",    cfg.coordmap.max_total_hits);
            break;

        case ToolMode::mapq_boost:
            ensure(!cfg.in.mapFile.empty(),    "-m/--map is required");
            ensure(!cfg.homq.in_bam.empty(),   "-i/--in is required");
            ensure(!cfg.homq.out_bam.empty(),  "-o/--out is required");
            ensure(cfg.global.IOthreads >= 1,  "--io_threads must be >= 1");
            ensure(cfg.global.threads  >= 1,   "-t/--threads must be >= 1");
            ensure(cfg.homq.batch_size >= 1,   "--batch must be >= 1");
            ensure(cfg.homq.min_frac > 0.0 && cfg.homq.min_frac <= 1.0, "--min_frac must be in (0,1]");
            ensure(cfg.homq.mapq_low >= 0 && cfg.homq.mapq_low <= 255,   "--mapq_low must be in [0,255]");
            ensure(cfg.homq.mapq_new >= 0 && cfg.homq.mapq_new <= 255,   "--mapq_new must be in [0,255]");
            ensure(cfg.homq.min_equiv >= 1,    "--min_equivalents must be >= 1");
            ensure(cfg.coordmap.max_hops       >= 1, "--cm_max_hops must be >= 1");
            ensure(cfg.coordmap.max_fanout     >= 1, "--cm_max_fanout must be >= 1");
            ensure(cfg.coordmap.min_len        >= 1, "--cm_min_len must be >= 1");
            ensure(cfg.coordmap.min_frac       >= 0.0 && cfg.coordmap.min_frac <= 1.0, "--cm_min_frac must be in [0,1]");
            ensure(cfg.coordmap.max_total_hits >= 1, "--cm_max_hits must be >= 1");
            kv("Coordinate map",          cfg.in.mapFile);
            kv("Input BAM",               cfg.homq.in_bam);
            kv("Output BAM",              cfg.homq.out_bam);
            kv("Worker threads",          cfg.global.threads);
            kv("HTS I/O threads",         cfg.global.IOthreads);
            kv("Batch size",              cfg.homq.batch_size);
            kv("Only raise if MAPQ<=Q",   (int)cfg.homq.mapq_low);
            kv("New MAPQ value",          (int)cfg.homq.mapq_new);
            kv("Min covered fraction",    cfg.homq.min_frac);
            kv("Min equivalent contigs",  cfg.homq.min_equiv);
            kv("Name check",              onoff(cfg.homq.name_check));
            kv("Max hops",                cfg.coordmap.max_hops);
            kv("Max fanout",              cfg.coordmap.max_fanout);
            kv("Min keep length",         cfg.coordmap.min_len);
            kv("Min keep fraction",       cfg.coordmap.min_frac);
            kv("Max total hits",          cfg.coordmap.max_total_hits);
            break;
    }
}

void help(char** argv) {
    std::cerr 
        << "Usage: " << argv[0] << " <subcommand> [options]\n\n"
        << program::description << "\n"
        << "Version: " << program::version << "\n"
        << "Date:    " << program::build_date << "\n"
        // << "Author:  " << program::author << "\n"
        << "\n"
        << "Subcommands:\n"
        << "  stat        collect statistics about a GFA file\n"
        << "  bubble      detect bubbles in a GFA\n"
        << "  seq         extract sequences along paths from a GFA\n"
        << "  depth       compute depth information from sequencing data for a GFA\n"
        << "  deoverlap   convert an overlap-based GFA to a non-overlap GFA\n"
        << "  collapse    collapse homologous sequences in bubbles into a single node\n"
        << "  align       align sequencing data to a GFA\n"
        << "  gfa2fa      export all segments to FASTA with two-sided extension for short nodes\n"
        << "  file2map    convert GFA/PAF to map format for liftover\n"
        << "  liftover    liftover coordinates using a coordinate map\n"
        << "  mapq_boost  raise MAPQ for alignments in homologous regions using a coordinate map\n\n";
}

void help_stat(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE\n\n"
        << "Collect statistics about a GFA file\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa    FILE    input GFA file\n\n"
        << "General Options:\n"
        << "  -h, --help           show this help\n\n";
}

AppConfig main_stat(int argc, char** argv) {
    if (argc < 3) { help_stat(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::stat;

    const struct option long_opts[] = {
        {"gfa",     required_argument, nullptr, 'g'},
        {"help",    no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:h";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': cfg.in.gfaFile = optarg; break;
            case 'h': help_stat(argv); std::exit(0);
            default:  help_stat(argv); std::exit(1);
        }
    }
    if (cfg.in.gfaFile.empty()) { error_stream() << "--gfa is required\n"; help_stat(argv); std::exit(1); }
    return cfg;
}

void help_bubble(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE -o FILE [options]\n\n"
        << "Detect bubbles in a GFA graph and write them to a GFA file\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa           FILE    input GFA file\n"
        << "  -o, --output        FILE    output GFA containing detected bubbles [stdout]\n\n"
        << "Detection Options:\n"
        << "      --depth         INT     maximum depth for bubble detection [" << BubbleOpts().max_depth << "]\n"
        << "      --paths         INT     maximum number of paths to enumerate per bubble [" << BubbleOpts().max_paths << "]\n"
        << "      --DFS_guard     INT     max DFS states [" << BubbleOpts().DFS_guard << "]\n"
        << "      --cx_depth      INT     local BFS depth for complexity check [" << BubbleOpts().cx_depth << "]\n"
        << "      --cx_nodes      INT     local visited node cap [" << BubbleOpts().cx_nodes << "]\n"
        << "      --cx_branches   INT     local branching node cap [" << BubbleOpts().cx_branches << "]\n"
        << "      --cx_deg_branch INT     degree>=this counts as branching [" << BubbleOpts().cx_deg_branch << "]\n"
        << "      --cx_deg_hub    INT     degree>=this counts as hub (complex) [" << BubbleOpts().cx_deg_hub << "]\n"
        << "      --cx_deg_cap    INT     stop counting degree after this [" << BubbleOpts().cx_deg_cap << "]\n"
        << "      --path_diff     FLOAT   node difference ratio threshold for haplotype separation [" << BubbleOpts().path_diff << "]\n"
        << "      --min_len       INT     minimum total sequence length of a bubble for output (0 = no filter) [" << BubbleOpts().min_len << "]\n"
        << "      --min_num       INT     minimum number of nodes inside a bubble required for output (0 = no filter) [" << BubbleOpts().min_num << "]\n\n"
        << "General Options:\n"
        << "  -t, --threads       INT     threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug                 debug mode (forces threads=1)\n"
        << "  -h, --help                  show this help\n\n";
}

AppConfig main_bubble(int argc, char** argv) {
    if (argc < 3) { help_bubble(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::bubble;

    const struct option long_opts[] = {
        {"gfa",          required_argument, nullptr, 'g'},
        {"output",       required_argument, nullptr, 'o'},
        {"depth",        required_argument, nullptr, 1001},
        {"paths",        required_argument, nullptr, 1002},
        {"DFS_guard",    required_argument, nullptr, 1003},
        {"cx_depth",     required_argument, nullptr, 1004},
        {"cx_nodes",     required_argument, nullptr, 1005},
        {"cx_branches",  required_argument, nullptr, 1006},
        {"cx_deg_branch",required_argument, nullptr, 1007},
        {"cx_deg_hub",   required_argument, nullptr, 1008},
        {"cx_deg_cap",   required_argument, nullptr, 1009},
        {"path_diff",    required_argument, nullptr, 1010},
        {"min_len",      required_argument, nullptr, 1011},
        {"min_num",      required_argument, nullptr, 1012},
        {"threads",      required_argument, nullptr, 't'},
        {"debug",        no_argument,       nullptr, 'd'},
        {"help",         no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:o:t:dh";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g':  cfg.in.gfaFile          = optarg;            break;
            case 'o':  cfg.out.gfaBubbles      = optarg;            break;
            case 1001: cfg.bubble.max_depth    = std::stoi(optarg); break;
            case 1002: cfg.bubble.max_paths    = std::stoi(optarg); break;
            case 1003: cfg.bubble.DFS_guard    = (uint64_t)std::stoull(optarg); break;
            case 1004: cfg.bubble.cx_depth     = (uint16_t)std::stoul(optarg);  break;
            case 1005: cfg.bubble.cx_nodes     = (uint32_t)std::stoul(optarg);  break;
            case 1006: cfg.bubble.cx_branches  = (uint32_t)std::stoul(optarg);  break;
            case 1007: cfg.bubble.cx_deg_branch= std::stoi(optarg); break;
            case 1008: cfg.bubble.cx_deg_hub   = std::stoi(optarg); break;
            case 1009: cfg.bubble.cx_deg_cap   = std::stoi(optarg); break;
            case 1010: cfg.bubble.path_diff    = std::stod(optarg); break;
            case 1011: cfg.bubble.min_len      = std::stoi(optarg); break;
            case 1012: cfg.bubble.min_num      = std::stoi(optarg); break;
            case 't':  cfg.global.threads      = std::stoi(optarg); break;
            case 'd':  cfg.global.debug        = true;              break;
            case 'h': help_bubble(argv); std::exit(0);
            default:  help_bubble(argv); std::exit(1);
        }
    }

    validate_and_print(cfg);

    return cfg;
}

void help_seq(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE -p FILE [-o FILE]\n\n"
        << "Extract sequences for compact marker paths (one path per line)\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa       FILE    input GFA file\n"
        << "  -p, --path      FILE    text file of paths (one per line)\n"
        << "                           (path format: '>S1<S2>S3' or 'S1+,S2-,S3+')\n"
        << "  -o, --output    FILE    output file, fasta format [stdout]\n\n"
        << "General Options:\n"
        << "  -h, --help              show this help\n\n";
}

AppConfig main_seq(int argc, char** argv) {
    if (argc < 3) { help_seq(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::seq;

    const struct option long_opts[] = {
        {"gfa",       required_argument, nullptr, 'g'},
        {"path",      required_argument, nullptr, 'p'},
        {"output",    required_argument, nullptr, 'o'},
        {"help",      no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:p:o:h";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': cfg.in.gfaFile       = optarg; break;
            case 'p': cfg.in.pathFile      = optarg; break;
            case 'o': cfg.out.gfaSeq       = optarg; break;
            case 'h': help_seq(argv); std::exit(0);
            default:  help_seq(argv); std::exit(1);
        }
    }
    
    validate_and_print(cfg);
    return cfg;
}

void help_depth(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE [--gaf FILE | -r FILE ...] [options]\n\n"
        << "Compute node depth information from sequencing data for a GFA graph\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa         FILE    input GFA file\n"
        << "  -r, --read        FILE    input sequencing data file(s) for **k-mer mode**\n"
        << "                             (depth calculated from k-mer coverage)\n"
        << "      --gaf         FILE    input GAF alignment file for **alignment mode**\n"
        << "                             (depth calculated directly from alignment)\n"
        << "  -o, --output      FILE    output depth information to file [stdout]\n\n"
        << "Index options:\n"
        << "  -k, --kmer        INT     k-mer length [" << GlobalOpts().kmerLen << "]\n\n"
        << "GFA loading options:\n"
        << "  --min_mapq        INT     minimum mapping quality to keep [" << DepthOpts().min_mapq << "]\n"
        << "  --min_frac        FLOAT   minimum aligned-length / read-length ratio to keep [" << DepthOpts().min_frac << "]\n\n"
        << "General Options:\n"
        << "  -t, --threads     INT     threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug               debug mode (forces threads=1)\n"
        << "  -h, --help                show this help\n\n";
}

AppConfig main_depth(int argc, char** argv) {
    if (argc < 3) { help_depth(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::depth;

    const struct option long_opts[] = {
        {"gfa",       required_argument, nullptr, 'g'},
        {"read",      required_argument, nullptr, 'r'},
        {"gaf",       required_argument, nullptr, 1000},
        {"output",    required_argument, nullptr, 'o'},
        {"kmer",      required_argument, nullptr, 'k'},
        {"min_mapq",  required_argument, nullptr, 1001},
        {"min_frac",  required_argument, nullptr, 1002},
        {"threads",   required_argument, nullptr, 't'},
        {"debug",     no_argument,       nullptr, 'd'},
        {"help",      no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:r:o:k:t:dh";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': cfg.in.gfaFile = optarg; break;
            case 'r': {
                cfg.in.reads.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.in.reads.emplace_back(argv[optind]); ++optind;
                }
                break;
            }
            case 1000: cfg.in.gafFile        = optarg;                           break;
            case 'o':  cfg.out.gfaDepth      = optarg;                           break;
            case 'k':  cfg.global.kmerLen    = std::stoi(optarg);                break;
            case 1001: cfg.depth.min_mapq    = std::max(0, std::atoi(optarg));   break;
            case 1002: cfg.depth.min_frac    = std::max(0.0, std::atof(optarg)); break;
            case 't':  cfg.global.threads    = std::max(1, std::stoi(optarg));   break;
            case 'd':  cfg.global.debug      = true;                             break;
            case 'h':  help_depth(argv); std::exit(0);
            default:   help_depth(argv); std::exit(1);
        }
    }
    
    validate_and_print(cfg);

    return cfg;
}

void help_deoverlap(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE [options]\n\n"
        << "Convert an overlap-based GFA into a non-overlap GFA\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa              FILE      input GFA file\n"
        << "  -p, --prefix           STR       output prefix [" << IOOutOpts().prefix << "]\n"
        << "                                    - prefix.deoverlap.gfa and prefix.deoverlap.noseq.gfa\n\n"
        << "Index options:\n"
        << "  -k, --kmer             INT       k-mer length [" << GlobalOpts().kmerLen << "]\n"
        << "  -w, --window           INT       minimizer window [" << GlobalOpts().minimizerW << "]\n\n"
        << "Align options:\n"
        << "      --use_wfa                    whether to use WFA for alignment, if not use mm2\n\n"
        << "Collapse options:\n"
        << "      --min_eq           INT       minimum match length to add cut points at segment ends [" << CollapseOpts().min_eq << "]\n"
        << "      --max_iters        INT       maximum iterations for cut point propagation [" << CollapseOpts().max_iters << "]\n\n"
        << "General Options:\n"
        << "  -t, --threads          INT       threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug                      debug mode\n"
        << "  -h, --help                       show this help\n\n";
}

AppConfig main_deoverlap(int argc, char** argv) {
    if (argc < 3) { help_deoverlap(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::deoverlap;

    const struct option long_opts[] = {
        {"gfa",        required_argument, nullptr, 'g'},
        {"prefix",     required_argument, nullptr, 'p'},
        {"kmer",       required_argument, nullptr, 'k'},
        {"window",     required_argument, nullptr, 'w'},
        {"use_wfa",    no_argument,       nullptr, 2001},
        {"min_eq",     required_argument, nullptr, 3002},
        {"max_iters",  required_argument, nullptr, 3003},
        {"threads",    required_argument, nullptr, 't'},
        {"debug",      no_argument,       nullptr, 'd'},
        {"help",       no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:p:k:w:t:dh";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g':   cfg.in.gfaFile              = optarg;                         break;
            case 'p':   cfg.out.prefix              = optarg;                         break;
            case 'k':   cfg.global.kmerLen          = std::max(1, std::stoi(optarg)); break;
            case 'w':   cfg.global.minimizerW       = std::max(1, std::stoi(optarg)); break;
            case 2001:  cfg.global.use_wfa          = true;                           break;
            case 3002:  cfg.collapse.min_eq         = std::stoi(optarg);              break;
            case 3003:  cfg.collapse.max_iters      = std::stoi(optarg);              break;
            case 't':   cfg.global.threads          = std::max(1, std::stoi(optarg)); break;
            case 'd':   cfg.global.debug            = true;                           break;
            case 'h':   help_deoverlap(argv); std::exit(0);
            default:    help_deoverlap(argv); std::exit(1);
        }
    }
    
    validate_and_print(cfg);

    return cfg;
}

void help_collapse(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE [options]\n\n"
        << "Collapse homologous sequences in bubbles into a single node\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa              FILE      input GFA file (typically de-overlapped)\n"
        << "  -p, --prefix           STR       output prefix [" << IOOutOpts().prefix << "]\n"
        << "                                    - prefix.collapse.gfa and prefix.collapse.noseq.gfa\n\n"
        << "Index options:\n"
        << "  -k, --kmer             INT       k-mer length [" << GlobalOpts().kmerLen << "]\n"
        << "  -w, --window           INT       minimizer window [" << GlobalOpts().minimizerW << "]\n\n"
        << "Align options:\n"
        << "      --use_wfa                    whether to use WFA for alignment, if not use mm2\n\n"
        << "Collapse options:\n"
        << "      --min_jaccard      FLOAT     nodes with similarity above this value will be collapsed [" << CollapseOpts().min_jaccard << "]\n"
        << "      --min_eq           INT       minimum match length to add cut points at segment ends [" << CollapseOpts().min_eq << "]\n"
        << "      --max_iters        INT       maximum iterations for cut point propagation [" << CollapseOpts().max_iters << "]\n"
        << "      --min_new_frac     FLOAT     submit align only if the unaligned fraction ≥ threshold [" << CollapseOpts().min_new_frac << "]\n\n"
        << "Bubble Detection Options:\n"
        << "      --depth            INT       maximum depth for bubble detection [" << BubbleOpts().max_depth << "]\n"
        << "      --paths            INT       maximum number of paths to enumerate per bubble [" << BubbleOpts().max_paths << "]\n"
        << "      --DFS_guard        INT       max DFS states [" << BubbleOpts().DFS_guard << "]\n"
        << "      --cx_depth         INT       local BFS depth [" << BubbleOpts().cx_depth << "]\n"
        << "      --cx_nodes         INT       local visited node cap [" << BubbleOpts().cx_nodes << "]\n"
        << "      --cx_branches      INT       branching node cap [" << BubbleOpts().cx_branches << "]\n"
        << "      --cx_deg_branch    INT       degree>=this counts as branching [" << BubbleOpts().cx_deg_branch << "]\n"
        << "      --cx_deg_hub       INT       degree>=this counts as hub (complex) [" << BubbleOpts().cx_deg_hub << "]\n"
        << "      --cx_deg_cap       INT       stop counting degree after this [" << BubbleOpts().cx_deg_cap << "]\n"
        << "      --path_diff        FLOAT     node difference ratio threshold for haplotype separation [" << BubbleOpts().path_diff << "]\n\n"
        << "General Options:\n"
        << "  -t, --threads          INT       threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug                      debug mode (forces threads=1)\n"
        << "  -h, --help                       show this help\n\n";
}

AppConfig main_collapse(int argc, char** argv) {
    if (argc < 3) { help_collapse(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::collapse;

    const struct option long_opts[] = {
        {"gfa",             required_argument, nullptr, 'g'},
        {"prefix",          required_argument, nullptr, 'p'},
        {"kmer",            required_argument, nullptr, 'k'},
        {"window",          required_argument, nullptr, 'w'},
        {"use_wfa",         no_argument,       nullptr, 2001},
        {"min_jaccard",     required_argument, nullptr, 3001},
        {"min_eq",          required_argument, nullptr, 3002},
        {"max_iters",       required_argument, nullptr, 3003},
        {"min_new_frac",    required_argument, nullptr, 3004},
        {"depth",           required_argument, nullptr, 4001},
        {"paths",           required_argument, nullptr, 4002},
        {"DFS_guard",       required_argument, nullptr, 4003},
        {"cx_depth",        required_argument, nullptr, 4004},
        {"cx_nodes",        required_argument, nullptr, 4005},
        {"cx_branches",     required_argument, nullptr, 4006},
        {"cx_deg_branch",   required_argument, nullptr, 4007},
        {"cx_deg_hub",      required_argument, nullptr, 4008},
        {"cx_deg_cap",      required_argument, nullptr, 4009},
        {"path_diff",       required_argument, nullptr, 4010},
        {"min_len",         required_argument, nullptr, 4011},
        {"min_num",         required_argument, nullptr, 4012},
        {"threads",         required_argument, nullptr, 't'},
        {"debug",           no_argument,       nullptr, 'd'},
        {"help",            no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:p:k:w:t:dh";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g':   cfg.in.gfaFile              = optarg;                         break;
            case 'p':   cfg.out.prefix              = optarg;                         break;
            case 'k':   cfg.global.kmerLen          = std::max(1, std::stoi(optarg)); break;
            case 'w':   cfg.global.minimizerW       = std::max(1, std::stoi(optarg)); break;
            case 2001:  cfg.global.use_wfa          = true;                           break;
            case 3001:  cfg.collapse.min_jaccard    = std::stod(optarg);              break;
            case 3002:  cfg.collapse.min_eq         = std::stoi(optarg);              break;
            case 3003:  cfg.collapse.max_iters      = std::stoi(optarg);              break;
            case 3004:  cfg.collapse.min_new_frac   = std::stod(optarg);              break;
            case 4001:  cfg.bubble.max_depth        = std::stoi(optarg);              break;
            case 4002:  cfg.bubble.max_paths        = std::stoi(optarg);              break;
            case 4003:  cfg.bubble.DFS_guard        = (uint64_t)std::stoull(optarg);  break;
            case 4004:  cfg.bubble.cx_depth         = (uint16_t)std::stoul(optarg);   break;
            case 4005:  cfg.bubble.cx_nodes         = (uint32_t)std::stoul(optarg);   break;
            case 4006:  cfg.bubble.cx_branches      = (uint32_t)std::stoul(optarg);   break;
            case 4007:  cfg.bubble.cx_deg_branch    = std::stoi(optarg);              break;
            case 4008:  cfg.bubble.cx_deg_hub       = std::stoi(optarg);              break;
            case 4009:  cfg.bubble.cx_deg_cap       = std::stoi(optarg);              break;
            case 4010:  cfg.bubble.path_diff        = std::stod(optarg);              break;
            case 4011:  cfg.bubble.min_len          = (uint32_t)std::stoul(optarg);   break;
            case 4012:  cfg.bubble.min_num          = (uint32_t)std::stoul(optarg);   break;
            case 't':   cfg.global.threads          = std::max(1, std::stoi(optarg)); break;
            case 'd':   cfg.global.debug            = true;                           break;
            case 'h':   help_collapse(argv); std::exit(0);
            default:    help_collapse(argv); std::exit(1);
        }
    }
    
    validate_and_print(cfg);

    return cfg;
}

void help_align(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE -r FILE ... [options]\n\n"
        << "Align sequencing data to a GFA\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa         FILE    input GFA file\n"
        << "  -r, --read        FILE    input sequencing data file(s), one argument can take multiple files\n"
        << "  -p, --paf                 output PAF (default SAM)\n"
        << "  -o, --output      FILE    output alignments to file [stdout]\n\n"
        << "Index options:\n"
        << "  -k, --kmer        INT     k-mer length [" << GlobalOpts().kmerLen << "]\n"
        << "  -w, --window      INT     minimizer window [" << GlobalOpts().minimizerW << "]\n\n"
        << "Alignment options:\n"
        << "  -S, --secondary   FLOAT   min ratio secondary/primary [" << MapOpts().sec_pri_ratio << "]\n"
        << "  -N                INT     number of secondary to keep [" << MapOpts().sec_pri_num << "]\n"
        << "  -P, --preset      STR     hifi/ont/illumina/other [" << MapOpts().preset << "]\n\n"
        << "General Options:\n"
        << "  -t, --threads     INT     threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug               debug mode (forces threads=1)\n"
        << "  -h, --help                show this help\n\n";
}

AppConfig main_align(int argc, char** argv) {
    if (argc < 3) { help_align(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::align;

    const struct option long_opts[] = {
        {"gfa",       required_argument, nullptr, 'g'},
        {"read",      required_argument, nullptr, 'r'},
        {"paf",       no_argument,       nullptr, 'p'},
        {"output",    required_argument, nullptr, 'o'},
        {"secondary", required_argument, nullptr, 'S'},
        {"N",         required_argument, nullptr, 'N'},
        {"preset",    required_argument, nullptr, 'P'},
        {"kmer",      required_argument, nullptr, 'k'},
        {"window",    required_argument, nullptr, 'w'},
        {"threads",   required_argument, nullptr, 't'},
        {"debug",     no_argument,       nullptr, 'd'},
        {"help",      no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:r:po:S:N:P:k:w:t:dh";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': cfg.in.gfaFile = optarg; break;
            case 'r': {
                cfg.in.reads.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.in.reads.emplace_back(argv[optind]); ++optind;
                }
                break;
            }
            case 'p': cfg.map.outPAF        = true;                           break;
            case 'o': cfg.out.alignOut      = optarg;                         break;
            case 'S': cfg.map.sec_pri_ratio = std::stod(optarg);              break;
            case 'N': cfg.map.sec_pri_num   = std::stoi(optarg);              break;
            case 'P': cfg.map.preset        = optarg;                         break;
            case 'k': cfg.global.kmerLen    = std::stoi(optarg);              break;
            case 'w': cfg.global.minimizerW = std::stoi(optarg);              break;
            case 't': cfg.global.threads    = std::max(1, std::stoi(optarg)); break;
            case 'd': cfg.global.debug      = true;                           break;
            case 'h': help_align(argv); std::exit(0);
            default:  help_align(argv); std::exit(1);
        }
    }
    
    validate_and_print(cfg);

    return cfg;
}

void help_gfa2fa(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE [options]\n\n"
        << "Export all segments of a GFA to FASTA, extending both ends by Y bp\n"
        << "for segments whose effective length is shorter than X bp.\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa          FILE     input GFA file\n"
        << "  -o, --output       FILE     output FASTA file [stdout]\n\n"
        << "Extension options:\n"
        << "      --min_len      UINT     only extend nodes with len < X [" << Gfa2FaOpts().min_len_xbp << "]\n"
        << "      --extend       UINT     extend Y bp on both sides [" << Gfa2FaOpts().extend_ybp << "]\n"
        << "      --wrap         UINT     FASTA wrap width; 0 = no wrap [" << Gfa2FaOpts().wrap_width << "]\n"
        << "      --keep_unknown          include nodes with '*' or zero length\n"
        << "                               (by default such nodes are skipped)\n\n"
        << "General Options:\n"
        << "  -h, --help                   show this help\n\n";
}

// ---------- gfa2fa main ----------
AppConfig main_gfa2fa(int argc, char** argv) {
    if (argc < 3) { help_gfa2fa(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::gfa2fa;

    // Long options for this subcommand
    const struct option long_opts[] = {
        {"gfa",         required_argument, nullptr, 'g'},
        {"output",      required_argument, nullptr, 'o'},
        {"min_len",     required_argument, nullptr, 3001},
        {"extend",      required_argument, nullptr, 3002},
        {"wrap",        required_argument, nullptr, 3003},
        {"keep_unknown",no_argument,       nullptr, 3004},
        {"help",        no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    // Short options only for the common flags in this subcommand
    const char* short_opts = "g:o:h";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': cfg.in.gfaFile           = optarg;            break;
            case 'o': cfg.out.gfa2faOut        = optarg;            break; // reuse gfa2faOut as FASTA-out
            case 3001: cfg.gfa2fa.min_len_xbp  = std::stoi(optarg); break;
            case 3002: cfg.gfa2fa.extend_ybp   = std::stoi(optarg); break;
            case 3003: cfg.gfa2fa.wrap_width   = std::stoi(optarg); break;
            case 3004: cfg.gfa2fa.skip_unknown = false;             break; // include '*' or zero-len nodes
            case 'h':  help_gfa2fa(argv); std::exit(0);
            default:   help_gfa2fa(argv); std::exit(1);
        }
    }

    validate_and_print(cfg);
    return cfg;
}

void help_file2map(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -i FILE [options]\n\n"
        << "convert GFA/PAF to map format for liftover\n\n"
        << "PAF input notes:\n"
        << "  * For PAF input, generate PAF with minimap2 in 1-to-1 mode:\n"
        << "     minimap2 -t8 -cx asm5 --secondary=no reference.fa query.fa > aln.paf\n\n"
        << "Input/Output:\n"
        << "  -i, --input        FILE     input PAF/GFA file\n"
        << "  -o, --output       FILE     output file name (default: stdout)\n\n"
        << "PAF options:\n"
        << "      --all                   include secondary/supplementary, (default: only primary tp:A:P)\n"
        << "      --min_len      INT      minimum alignment length to keep [" << File2mapOpts().min_len << "]\n"
        << "      --min_mapq     INT      minimum mapping quality to keep [" << File2mapOpts().min_mapq << "]\n\n"
        << "General Options:\n"
        << "  -h, --help                  show this help\n\n";
}

AppConfig main_file2map(int argc, char** argv) {
    if (argc < 3) { help_file2map(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::file2map;

    const struct option long_opts[] = {
        {"input",    required_argument, nullptr, 'i'},
        {"output",   required_argument, nullptr, 'o'},
        {"all",      no_argument,       nullptr, 1001},
        {"min_len",  required_argument, nullptr, 1002},
        {"min_mapq", required_argument, nullptr, 1003},
        {"help",     no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "i:o:h";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'i': cfg.in.inputFile = optarg; break;
            case 'o': cfg.out.outFile   = optarg; break;
            case 1001: cfg.file2map.paf_primary_only = false; break;
            case 1002: cfg.file2map.min_len = std::max(0, std::stoi(optarg)); break;
            case 1003: cfg.file2map.min_mapq = std::atoi(optarg); break;
            case 'h': help_file2map(argv); std::exit(0);
            default:  help_file2map(argv); std::exit(1);
        }
    }

    if (cfg.out.outFile.empty())  cfg.out.outFile = "-";

    validate_and_print(cfg);
    return cfg;
}

static AppConfig make_liftover_defaults() {
    AppConfig cfg;
    cfg.mode = ToolMode::liftover;
    cfg.coordmap.max_hops = 0;
    return cfg;
}

void help_liftover(char** argv, const AppConfig& cfg, bool print_all) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -m FILE [options]\n\n"
        << "LiftOver BED coordinates using a .map file, which produced by deoverlap/collapse/file2map\n\n"
        << "Input/Output:\n"
        << "  -m, --map           FILE     input mapping file (.map)\n"
        << "      --paf           FILE     input PAF file to assist liftover (optional)\n"
        << "  -b, --bed           FILE     input BED file, required if not checking\n"
        << "  -o, --output        FILE     output file name (default: stdout)\n"
        << "  -r, --ref           FILE     FASTA/FASTQ containing sequences for src+dst contigs, required if checking\n\n"
        << "Map Options:\n"
        << "      --regex         STR      keep hits whose label matches regex\n"
        << "                                (\"h\\d+tg\" to match h1tg, h2tg, ...)\n"
        << "      --min_frac      FLOAT    minimum fraction of the BED interval that must be lifted over [" << cfg.liftover.min_frac << "]\n"
        << "                                (0.9 means at least 90% of the BED region is mapped)\n"
        << "      --flank_win     INT      extension window size around BED ends to identify insertions [" << cfg.liftover.flank_win << "]\n"
        << "      --max_flank     INT      max extension around BED ends to identify insertions [" << cfg.liftover.max_flank << "]\n"
        << "      --max_gap       INT      max allowed |ref_gap - query_gap| when calling INS [" << cfg.liftover.max_gap << "]\n"
        << "      --max_hit       INT      max extended liftover results to keep [" << cfg.liftover.max_hit << "]\n\n"
        << "PAF Options:\n"
        << "      --min_len       INT      minimum alignment length for liftover [" << cfg.liftover.min_len << "]\n"
        << "      --min_mapq      INT      minimum mapping quality for liftover [" << cfg.liftover.min_mapq << "]\n\n"
        << "Check Options:\n"
        << "      --check                  run sequence consistency check, requires -r\n";

    if (print_all) {
        std::cerr
            << "  -w, --win           INT      window size [" << LiftoverOpts().win << "]\n"
            << "  -s, --step          INT      step size [" << LiftoverOpts().step << "]\n"
            << "      --max_examples  INT      max reported examples [" << LiftoverOpts().max_examples << "]\n";
    }
        
    std::cerr 
        << std::endl
        << "CoordMap Options:\n"
        << "      --cm_max_hops   INT      max hop depth when chaining names [" << cfg.coordmap.max_hops << "]\n"
        << "                                (0 allows source -> A; 1 allows source -> A -> B)\n";

    if (print_all) {
        std::cerr
            << "      --cm_max_fanout INT      max fanout per search [" << cfg.coordmap.max_fanout << "]\n"
            << "      --cm_min_len    INT      min length (bp) to keep a hit for hops [" << cfg.coordmap.min_len << "]\n"
            << "      --cm_min_frac   FLOAT    min fraction of current interval length to keep a hit [" << cfg.coordmap.min_frac << "]\n"
            << "                                (threshold = max(min_len, cur_len * min_frac))\n"
            << "      --cm_max_hits   INT      global BFS node cap [" << cfg.coordmap.max_total_hits << "]\n";
    }

    std::cerr
        << std::endl
        << "General Options:\n"
        << "  -t, --threads       INT      worker threads [" << GlobalOpts().threads << "]\n"
        << "  -h, --help                   show this help\n"
        << "      --help_all               show all options\n\n";
}

AppConfig main_liftover(int argc, char** argv) {
    AppConfig cfg = make_liftover_defaults();

    if (argc < 3) { help_liftover(argv, cfg); std::exit(1); }

    const struct option long_opts[] = {
        {"map",           required_argument, nullptr, 'm'},
        {"paf",           required_argument, nullptr, 0001},
        {"bed",           required_argument, nullptr, 'b'},
        {"output",        required_argument, nullptr, 'o'},
        {"ref",           required_argument, nullptr, 'r'},
        {"regex",         required_argument, nullptr, 1001},
        {"min_frac",      required_argument, nullptr, 1002},
        {"flank_win",     required_argument, nullptr, 1003},
        {"max_flank",     required_argument, nullptr, 1004},
        {"max_gap",       required_argument, nullptr, 1005},
        {"max_hit",       required_argument, nullptr, 1006},
        {"min_len",       required_argument, nullptr, 2001},
        {"min_mapq",      required_argument, nullptr, 2002},
        {"check",         no_argument,       nullptr, 3001},
        {"win",           required_argument, nullptr, 'w'},
        {"step",          required_argument, nullptr, 's'},
        {"max_examples",  required_argument, nullptr, 3002},
        {"cm_max_hops",   required_argument, nullptr, 4001},
        {"cm_max_fanout", required_argument, nullptr, 4002},
        {"cm_min_len",    required_argument, nullptr, 4003},
        {"cm_min_frac",   required_argument, nullptr, 4004},
        {"cm_max_hits",   required_argument, nullptr, 4005},
        {"threads",       required_argument, nullptr, 't'},
        {"help",          no_argument,       nullptr, 'h'},
        {"help_all",      no_argument,       nullptr, 5001},
        {0,0,0,0}
    };
    const char* short_opts = "m:b:o:r:w:s:t:h";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'm':  cfg.in.mapFile              = optarg;                                          break;
            case 0001: cfg.in.pafFile              = optarg;                                          break;
            case 'b':  cfg.in.bedFile              = optarg;                                          break;
            case 'o':  cfg.out.outFile             = optarg;                                          break;
            case 'r':  cfg.in.referenceFile        = optarg;                                          break;
            case 1001: cfg.liftover.regex          = optarg;                                          break;
            case 1002: cfg.liftover.min_frac       = std::max(0.0, std::min(1.0, std::atof(optarg))); break;
            case 1003: cfg.liftover.flank_win      = (uint32_t)std::stoul(optarg);                    break;
            case 1004: cfg.liftover.max_flank      = (uint32_t)std::stoul(optarg);                    break;
            case 1005: cfg.liftover.max_gap        = (uint32_t)std::stoul(optarg);                    break;
            case 1006: cfg.liftover.max_hit        = (uint16_t)std::stoul(optarg);                    break;
            case 2001: cfg.liftover.min_len        = std::max(0, std::atoi(optarg));                  break;
            case 2002: cfg.liftover.min_mapq       = std::atoi(optarg);                               break;
            case 3001: cfg.liftover.do_check       = true;                                            break;
            case 'w':  cfg.liftover.win            = (uint32_t)std::stoul(optarg);                    break;
            case 's':  cfg.liftover.step           = (uint32_t)std::stoul(optarg);                    break;
            case 3002: cfg.liftover.max_examples   = (uint32_t)std::stoul(optarg);                    break;
            case 4001: cfg.coordmap.max_hops       = std::max(0, std::atoi(optarg));                  break;
            case 4002: cfg.coordmap.max_fanout     = (uint32_t)std::max(1, std::atoi(optarg));        break;
            case 4003: cfg.coordmap.min_len        = (uint32_t)std::max(1, std::atoi(optarg));        break;
            case 4004: cfg.coordmap.min_frac       = std::max(0.0, std::atof(optarg));                break;
            case 4005: cfg.coordmap.max_total_hits = (uint32_t)std::max(1, std::atoi(optarg));        break;
            case 't':  cfg.global.threads          = std::max(1, std::stoi(optarg));                  break;
            case 'h':  help_liftover(argv, cfg);       std::exit(0);
            case 5001: help_liftover(argv, cfg, true); std::exit(0);
            default:   help_liftover(argv, cfg);       std::exit(1);
        }
    }

    if (cfg.out.outFile.empty())  cfg.out.outFile = "-";

    validate_and_print(cfg);
    return cfg;
}

void help_mapq_boost(char** argv) {
    static const char* kTail =
        "    | liftasm mapq_boost -m merge.map \\\n"
        "    | samtools view -b -q 10 -o out.boost.bam\n\n";

    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -m FILE [options]\n\n"
        << "Raise MAPQ for alignments whose genome regions are homologous according to a coordinate map (generated by deoverlap/collapse/file2map)\n\n"
        << "Input requirements (IMPORTANT):\n"
        << "  * All alignments must be present (primary + secondary).\n"
        << "     This includes secondary alignments recorded as separate records\n"
        << "      OR alternative hits recorded in tags (e.g. BWA's XA tag in the primary record).\n"
        << "     In short, any competing alignment that causes the primary MAPQ to be low\n"
        << "      (or zero) must be visible to this program, either as separate records or encoded in tags.\n"
        << "  * All alignments of the same read must appear consecutively.\n"
        << "     If your BAM/SAM comes directly from the aligner (without coordinate sorting), \n"
        << "      records for the same read are usually already grouped together.\n"
        << "     If unsure, run `samtools sort -n` (slower, more I/O). Alternatively, add --name_check to make the program verify the order.\n\n"
        << "Example pipeline:\n"
        << "  # BWA\n"
        << "  bwa mem -a HG002.hap1_hap2.fa read.1.fq.gz read.2.fq.gz \\\n"
        << kTail
        << "  # HISAT2\n"
        << "  hisat2 -x INDEX -1 read.1.fq.gz -2 read.2.fq.gz \\\n"
        << kTail
        << "  # STAR\n"
        << "  STAR --genomeDir ./ --readFilesIn read.1.fq.gz read.2.fq.gz --readFilesCommand zcat \\\n"
        << "       --outSAMtype BAM Unsorted --outStd BAM_Unsorted \\\n"
        << kTail
        << "Input/Output:\n"
        << "  -m, --map             FILE    coordinate mapping file (generated by deoverlap/collapse/file2map)\n"
        << "  -i, --in              FILE    input SAM/BAM/CRAM (default: stdin)\n"
        << "  -o, --out             FILE    output SAM/BAM (default: stdout)\n\n"
        << "Boost Options:\n"
        << "      --batch           INT     batch size [" << MapqBoostOpts().batch_size << "]\n"
        << "      --mapq_low        INT     only raise if MAPQ<=Q (no larger than 255) [" << (int)MapqBoostOpts().mapq_low << "]\n"
        << "      --mapq_new        INT     new MAPQ value (no larger than 255) [" << (int)MapqBoostOpts().mapq_new << "]\n"
        << "      --min_frac        FLOAT   min coverage fraction (no larger than 1) [" << MapqBoostOpts().min_frac << "]\n"
        << "      --min_equivalents INT     min #equivalent contigs [" << MapqBoostOpts().min_equiv << "]\n"
        << "      --name_check              check whether QNAMEs are name-sorted (off by default)\n\n"
        << "CoordMap Options:\n"
        << "      --cm_max_hops     INT     BFS depth limit [" << CoordMapOpts().max_hops << "]\n"
        << "                                 (e.g. 1 allows: source -> A; 2 allows: source -> A -> B)\n"
        << "      --cm_max_fanout   INT     max fanout per search [" << CoordMapOpts().max_fanout << "]\n"
        << "      --cm_min_len      INT     min length (bp) to keep a hit for hops [" << CoordMapOpts().min_len << "]\n"
        << "      --cm_min_frac     FLOAT   min fraction of current interval length to keep a hit [" << CoordMapOpts().min_frac << "]\n"
        << "                                 (threshold = max(min_len, cur_len * min_frac))\n"
        << "      --cm_max_hits     INT     global BFS node cap [" << CoordMapOpts().max_total_hits << "]\n\n"
        << "General Options:\n"
        << "  -t, --threads         INT     worker threads [" << GlobalOpts().threads << "]\n"
        << "      --io_threads      INT     HTS I/O threads [" << GlobalOpts().IOthreads << "]\n"
        << "  -h, --help                    show this help\n\n";
}

AppConfig main_mapq_boost(int argc, char** argv) {
    if (argc < 3) { help_mapq_boost(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::mapq_boost;

    const struct option long_opts[] = {
        {"map",             required_argument, nullptr, 'm'},
        {"in",              required_argument, nullptr, 'i'},
        {"out",             required_argument, nullptr, 'o'},
        {"batch",           required_argument, nullptr, 1001},
        {"mapq_low",        required_argument, nullptr, 1002},
        {"mapq_new",        required_argument, nullptr, 1003},
        {"min_frac",        required_argument, nullptr, 1004},
        {"min_equivalents", required_argument, nullptr, 1005},
        {"name_check",      no_argument,       nullptr, 1006},
        {"cm_max_hops",     required_argument, nullptr, 2001},
        {"cm_max_fanout",   required_argument, nullptr, 2002},
        {"cm_min_len",      required_argument, nullptr, 2003},
        {"cm_min_frac",     required_argument, nullptr, 2004},
        {"cm_max_hits",     required_argument, nullptr, 2005},
        {"threads",         required_argument, nullptr, 't'},
        {"io_threads",      required_argument, nullptr, 3001},
        {"help",            no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "m:i:o:t:h";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'm':  cfg.in.mapFile              = optarg;                                         break;
            case 'i':  cfg.homq.in_bam             = optarg;                                         break;
            case 'o':  cfg.homq.out_bam            = optarg;                                         break;
            case 1001: cfg.homq.batch_size         = std::max(1, std::atoi(optarg));                 break;
            case 1002: cfg.homq.mapq_low           = (uint8_t)std::clamp(std::atoi(optarg), 0, 255);  break;
            case 1003: cfg.homq.mapq_new           = (uint8_t)std::clamp(std::atoi(optarg), 0, 255); break;
            case 1004: cfg.homq.min_frac           = std::max(0.0, std::atof(optarg));               break;
            case 1005: cfg.homq.min_equiv          = std::max(1, std::atoi(optarg));                 break;
            case 1006: cfg.homq.name_check         = true;                                           break;
            case 2001: cfg.coordmap.max_hops       = std::max(1, std::atoi(optarg));                 break;
            case 2002: cfg.coordmap.max_fanout     = (uint32_t)std::max(1, std::atoi(optarg));       break;
            case 2003: cfg.coordmap.min_len        = (uint32_t)std::max(1, std::atoi(optarg));       break;
            case 2004: cfg.coordmap.min_frac       = std::max(0.0, std::atof(optarg));               break;
            case 2005: cfg.coordmap.max_total_hits = (uint32_t)std::max(1, std::atoi(optarg));       break;
            case 3001: cfg.global.IOthreads        = std::max(1, std::atoi(optarg));                 break;
            case 't':  cfg.global.threads          = std::max(1, std::atoi(optarg));                 break;
            case 'h': help_mapq_boost(argv); std::exit(0);
            default: help_mapq_boost(argv);  std::exit(1);
        }
    }

    if (cfg.homq.in_bam.empty())  cfg.homq.in_bam  = "-";
    if (cfg.homq.out_bam.empty()) cfg.homq.out_bam = "-";

    validate_and_print(cfg);
    return cfg;
}