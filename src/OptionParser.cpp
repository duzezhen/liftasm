#include "../include/OptionParser.hpp"
#include "../include/ProgramMetadata.hpp"
#include "../include/logger.hpp"
#include "../include/get_time.hpp"

#include <getopt.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iomanip>


static inline bool is_flag_(const char* s) { return s && s[0] == '-'; }

static inline void parse_repeat_mask_(const char* s, CollapseOpts& opt) {
    unsigned a = 0, b = 0, c = 0;
    char extra = 0;

    if (!s || std::sscanf(s, "%u,%u,%u%c", &a, &b, &c, &extra) != 3) {
        error_stream() << "--repeat_mask expects MIN_LEN,MAX_PERIOD,MAX_MISMATCH, e.g. --repeat_mask 24,12,2\n";
        std::exit(1);
    }

    opt.repeat_mask_min_len = std::max(0u, a);
    opt.repeat_mask_max_period = std::max(0u, b);
    opt.repeat_mask_max_mismatch = std::max(0u, c);
}

static void parse_abnormal_seg_(const char* arg, CollapseOpts& opt) {
    std::string s(arg ? arg : "");
    const size_t p = s.find(',');

    if (p == std::string::npos) {
        error_stream() << "invalid --abnormal_seg, expected LEN,COUNT, got: " << s << "\n";
        std::exit(1);
    }

    opt.max_abnormal_cut_len = std::max<uint32_t>(0, static_cast<uint32_t>(std::stoul(s.substr(0, p))));
    opt.min_abnormal_cut_count = std::max<uint32_t>(0, static_cast<uint32_t>(std::stoul(s.substr(p + 1))));
}

static void validate_and_print(int argc, char** argv, AppConfig& cfg) {
    using std::left;
    using std::setw;

    auto die    = [&](const std::string& msg){ error_stream() << msg << "\n"; std::exit(1); };
    auto ensure = [&](bool ok, const std::string& msg){ if (!ok) die(msg); };
    auto onoff  = [](bool b){ return b ? "ON" : "OFF"; };

    const int KEYW = 25;

    auto kv = [&](const std::string& key, const auto& val) {
        std::ostringstream oss;
        oss << left << setw(KEYW) << key << ": " << val;
        log_stream() << oss.str() << "\n";
    };

    if (cfg.global.debug) { cfg.global.threads = 1; set_debug(true); }
    ensure(cfg.global.threads >= 1, "threads must be >= 1");
    if (cfg.map.use_wfa) {
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
        cfg.mode == ToolMode::mapq_boost ? "mapq_boost" :
        cfg.mode == ToolMode::ressit     ? "ressit" :
        "unknown";

    switch (cfg.mode) {
        case ToolMode::stat:
            ensure(
                !cfg.stat.gfaFiles.empty(),
                "-g/--gfa is required"
            );
            break;

        case ToolMode::seq:
            ensure(
                !cfg.seq.gfaFiles.empty(),
                "-g/--gfa is required"
            );
            ensure(
                !cfg.seq.pathFile.empty(), 
                "-p/--path is required"
            );
            break;

        case ToolMode::gfa2fa:
            ensure(
                !cfg.gfa2fa.gfaFiles.empty(), 
                "-g/--gfa is required"
            );
            ensure(
                cfg.gfa2fa.min_len_xbp >= 0, 
                "--min_len must be >= 0"
            );
            ensure(
                cfg.gfa2fa.extend_ybp >= 0, 
                "--extend must be >= 0"
            );
            ensure(
                cfg.gfa2fa.wrap_width >= 0, 
                "--wrap must be >= 0"
            );
            break;

        case ToolMode::depth:
            ensure(
                !cfg.depth.gfaFiles.empty(), 
                "-g/--gfa is required"
            );
            ensure(
                !(!cfg.depth.gafFile.empty() && !cfg.depth.reads.empty()), 
                "Options --gaf and -r/--read are mutually exclusive"
            );
            ensure(
                cfg.depth.min_frac >= 0.0 && cfg.depth.min_frac <= 1.0,
                "--min_frac must be in [0,1]"
            );
            ensure(
                cfg.depth.min_mapq >= 0 && cfg.depth.min_mapq <= 60,
                "--min_mapq must be in [0,60]"
            );
            break;

        case ToolMode::bubble:
            ensure(
                !cfg.bubble.gfaFiles.empty(),
                "-g/--gfa is required"
            );
            ensure(
                !cfg.bubble.out_prefix.empty(), 
                "-p/--prefix is required"
            );
            if (cfg.bubble.write_vcf) {
                ensure(
                    cfg.bubble.ali_min_len >= 0,
                    "--ali_min_len must be >= 0"
                );
                ensure(
                    cfg.bubble.ali_min_mapq >= 0,
                    "--ali_min_mapq must be >= 0"
                );
            }
            ensure(
                cfg.bubble.max_depth >= 1,
                "--depth must be >= 1"
            );
            ensure(
                cfg.bubble.max_paths >= 1,
                "--paths must be >= 1"
            );
            ensure(
                cfg.bubble.DFS_guard >= 1,
                "--DFS_guard must be >= 1"
            );
            ensure(
                cfg.bubble.cx_depth >= 0,
                "--cx_depth must be >= 0"
            );
            ensure(
                cfg.bubble.cx_nodes >= 1,
                "--cx_nodes must be >= 1"
            );
            ensure(
                cfg.bubble.cx_branches >= 1,
                "--cx_branches must be >= 1"
            );
            ensure(
                cfg.bubble.cx_deg_branch >= 1,
                "--cx_deg_branch must be >= 1"
            );
            ensure(
                cfg.bubble.cx_deg_hub >= 1,
                "--cx_deg_hub must be >= 1"
            );
            ensure(
                cfg.bubble.path_diff >= 0.0 && cfg.bubble.path_diff <= 1.0, 
                "--path_diff must be in [0,1]"
            );
            ensure(
                cfg.bubble.stall_round_limit >= 0, 
                "--stall_rounds must be >= 0"
            );
            ensure(
                cfg.bubble.min_len >= 0,
                "--min_len must be >= 0"
            );
            ensure(
                cfg.bubble.min_num >= 0,
                "--min_num must be >= 0"
            );
            ensure(
                cfg.bubble.same_sim >= 0.0 && cfg.bubble.same_sim <= 1.0, 
                "--same_sim must be in [0,1]"
            );
            ensure(
                cfg.bubble.same_min_num >= 0, 
                "--same_min_num must be >= 0"
            );
            ensure(
                cfg.bubble.same_min_len >= 0,
                "--same_min_len must be >= 0"
            );
            ensure(
                cfg.bubble.diff_min_src >= 0,
                "--diff_min_src must be >= 0"
            );
            ensure(
                cfg.bubble.homo_num >= 1,
                "--hnum must be >= 1"
            );
            ensure(
                cfg.bubble.diff_sim >= 0.0 && cfg.bubble.diff_sim <= 1.0, 
                "--diff_sim must be in [0,1]"
            );
            ensure(
                cfg.bubble.diff_min_num >= 0, 
                "--diff_min_num must be >= 0"
            );
            ensure(
                cfg.bubble.diff_min_len >= 0,
                "--diff_min_len must be >= 0"
            );
            ensure(
                cfg.bubble.homo_k >= 1 && cfg.bubble.homo_k <= 28, 
                "--hk must be in [1,28]"
            );
            ensure(
                cfg.bubble.homo_w > 0, 
                "--hw must be > 0"
            );
            ensure(
                cfg.bubble.homo_boot_bp > 0, 
                "--hboot must be > 0"
            );
            ensure(
                cfg.bubble.homo_cnt_len > 0, 
                "--hcnt must be > 0"
            );
            ensure(
                cfg.bubble.homo_bloom_bits > 0, 
                "--hbits must be > 0"
            );
            ensure(
                cfg.bubble.homo_bloom_hash > 0, 
                "--hhash must be > 0"
            );

            break;

        case ToolMode::deoverlap:
            ensure(
                !cfg.collapse.gfaFiles.empty(), 
                "-g/--gfa is required"
            );
            ensure(
                !cfg.collapse.prefix.empty(),
                "-p/--prefix is required"
            );
            ensure(
                cfg.collapse.min_eq >= 1,
                "--min_eq must be >= 1"
            );
            ensure(
                cfg.collapse.max_iters >= 1,
                "--max_iters must be >= 1"
            );
            ensure(
                cfg.collapse.min_abnormal_cut_count >= 0,
                "--abnormal_cut MAX_LEN should be >= 0"
            );
            ensure(
                cfg.collapse.max_abnormal_cut_len >= 0,
                "--abnormal_cut MIN_COUNT should be >= 0"
            );
            ensure(
                cfg.map.zdrop >= 0,
                "--zdrop must be >= 0"
            );
            ensure (
                cfg.collapse.min_match_ratio >= 0.0 && cfg.collapse.min_match_ratio <= 1.0,
                "--min_match must be in [0,1]"
            );
            break;

        case ToolMode::collapse:
            ensure(
                !cfg.collapse.gfaFiles.empty(), 
                "-g/--gfa is required"
            );
            ensure(
                !cfg.collapse.prefix.empty(), 
                "-p/--prefix is required"
            );
            ensure(
                cfg.collapse.min_jaccard >= 0.0 && cfg.collapse.min_jaccard <= 1.0,
                "--min_jaccard must be in [0,1]"
            );
            ensure(
                cfg.collapse.min_eq >= 1, 
                "--min_eq must be >= 1"
            );
            ensure(
                cfg.collapse.max_iters >= 1, 
                "--max_iters must be >= 1"
            );
            ensure(
                cfg.collapse.repeat_mask_min_len >= 0,
                "--repeat_mask MIN_LEN must be >= 0"
            );
            ensure(
                cfg.collapse.repeat_mask_max_period >= 0,
                "--repeat_mask MAX_PERIOD must be >= 0"
            );
            ensure(
                cfg.collapse.repeat_mask_max_period <= cfg.collapse.repeat_mask_min_len,
                "--repeat_mask MAX_PERIOD should be <= MIN_LEN"
            );
            ensure(
                cfg.collapse.min_abnormal_cut_count >= 0,
                "--abnormal_cut MAX_LEN should be >= 0"
            );
            ensure(
                cfg.collapse.max_abnormal_cut_len >= 0,
                "--abnormal_cut MIN_COUNT should be >= 0"
            );
            ensure(
                cfg.bubble.max_depth >= 1, 
                "--depth must be >= 1"
            );
            ensure(
                cfg.bubble.max_paths >= 1, 
                "--paths must be >= 1"
            );
            ensure(
                cfg.bubble.DFS_guard >= 1, 
                "--DFS_guard must be >= 1"
            );
            ensure(
                cfg.bubble.cx_depth >= 0, 
                "--cx_depth must be >= 0"
            );
            ensure(
                cfg.bubble.cx_nodes >= 1, 
                "--cx_nodes must be >= 1"
            );
            ensure(
                cfg.bubble.cx_branches >= 1, 
                "--cx_branches must be >= 1"
            );
            ensure(
                cfg.bubble.cx_deg_branch >= 1, 
                "--cx_deg_branch must be >= 1"
            );
            ensure(
                cfg.bubble.cx_deg_hub >= 1, 
                "--cx_deg_hub must be >= 1"
            );
            ensure(
                cfg.bubble.path_diff >= 0.0 && cfg.bubble.path_diff <= 1.0, 
                "--path_diff must be in [0,1]"
            );
            ensure(
                cfg.bubble.stall_round_limit >= 0, 
                "--stall_rounds must be >= 0"
            );
            ensure(
                cfg.bubble.same_sim >= 0.0 && cfg.bubble.same_sim <= 1.0, 
                "--same_sim must be in [0,1]"
            );
            ensure(
                cfg.bubble.same_min_num >= 0, 
                "--same_min_num must be >= 0"
            );
            ensure(
                cfg.bubble.same_min_len >= 0, 
                "--same_min_len must be >= 0"
            );
            ensure(
                cfg.bubble.diff_min_src >= 0, 
                "--diff_min_src must be >= 0"
            );
            ensure(
                cfg.bubble.homo_num >= 1,
                "--hnum must be >= 1"
            );
            ensure(
                cfg.bubble.diff_sim >= 0.0 && cfg.bubble.diff_sim <= 1.0, 
                "--diff_sim must be in [0,1]"
            );
            ensure(
                cfg.bubble.diff_min_num >= 0, 
                "--diff_min_num must be >= 0"
            );
            ensure(
                cfg.bubble.diff_min_len >= 0, 
                "--diff_min_len must be >= 0"
            );
            ensure(
                cfg.bubble.homo_k >= 1 && cfg.bubble.homo_k <= 28, 
                "--hk must be in [1,28]"
            );
            ensure(
                cfg.bubble.homo_w > 0, 
                "--hw must be > 0"
            );
            ensure(
                cfg.bubble.homo_boot_bp > 0, "--hboot must be > 0"
            );
            ensure(
                cfg.bubble.homo_cnt_len > 0, "--hcnt must be > 0"
            );
            ensure(
                cfg.bubble.homo_bloom_bits > 0, "--hbits must be > 0"
            );
            ensure(
                cfg.bubble.homo_bloom_hash > 0, "--hhash must be > 0"
            );
            ensure(
                cfg.map.zdrop >= 0, 
                "--zdrop must be >= 0"
            );
            ensure (
                cfg.collapse.min_match_ratio >= 0.0 && cfg.collapse.min_match_ratio <= 1.0,
                "--min_match must be in [0,1]"
            );
            break;

        case ToolMode::file2map:
            ensure(
                !cfg.file2map.inputFiles.empty(), 
                "-i/--input is required"
            );
            ensure(
                cfg.file2map.min_len >= 0, 
                "--min_len must be >= 0"
            );
            ensure(
                cfg.file2map.min_mapq >= 0 && cfg.file2map.min_mapq <= 60, 
                "--min_mapq must be in [0,60]"
            );
            break;

        case ToolMode::liftover:
            ensure(
                !cfg.liftover.mapFiles.empty(), 
                "-m/--map is required"
            );
            if (cfg.liftover.do_check) {
                ensure(
                    !cfg.liftover.referenceFile.empty(), 
                    "-r/--ref is required when using --check"
                );
            } else {
                ensure(
                    !cfg.liftover.bedFile.empty(), 
                    "-b/--bed is required when not using --check"
                );
            }
            ensure(
                cfg.liftover.min_frac >= 0.0 && cfg.liftover.min_frac <= 1.0, 
                "--min_frac must be in [0,1]"
            );
            ensure(
                cfg.liftover.flank_win >= 1, 
                "--flank_win must be >= 1"
            );
            ensure(
                cfg.liftover.max_flank >= 1, 
                "--max_flank must be >= 1"
            );
            ensure(
                cfg.liftover.max_gap >= 0, 
                "--max_gap must be >= 0"
            );
            ensure(
                cfg.liftover.max_hit >= 1, 
                "--max_hit must be >= 1"
            );
            ensure(
                cfg.liftover.min_len >= 0, 
                "--min_len must be >= 0"
            );
            ensure(
                cfg.liftover.min_mapq >= 0 && cfg.liftover.min_mapq <= 60, 
                "--min_mapq must be in [0,60]"
            );
            ensure(
                cfg.liftover.win >= 1, 
                "-w/--win must be >= 1"
            );
            ensure(
                cfg.liftover.step >= 1, 
                "-s/--step must be >= 1"
            );
            ensure(
                cfg.liftover.max_examples >= 1, 
                "--max_examples must be >= 1"
            );
            ensure(
                cfg.coordmap.max_hops >= 1, 
                "--cm_max_hops must be >= 1"
            );
            ensure(
                cfg.coordmap.max_fanout >= 1, 
                "--cm_max_fanout must be >= 1"
            );
            ensure(
                cfg.coordmap.min_len >= 1, 
                "--cm_min_len must be >= 1"
            );
            ensure(
                cfg.coordmap.min_frac >= 0.0 && cfg.coordmap.min_frac <= 1.0, 
                "--cm_min_frac must be in [0,1]"
            );
            ensure(
                cfg.coordmap.max_total_hits >= 1, 
                "--cm_max_hits must be >= 1"
            );
            break;

        case ToolMode::mapq_boost:
            ensure(
                !cfg.homq.mapFiles.empty(),
                "-m/--map is required"
            );
            ensure(
                !cfg.homq.in_bam.empty(),
                "-i/--in is required"
            );
            ensure(
                !cfg.homq.out_bam.empty(),
                "-o/--out is required"
            );
            ensure(
                cfg.homq.batch_size >= 1,
                "--batch must be >= 1"
            );
            ensure(
                cfg.homq.mapq_low >= 0 && cfg.homq.mapq_low <= 255,
                "--mapq_low must be in [0,255]"
            );
            ensure(
                cfg.homq.mapq_cap >= 0 && cfg.homq.mapq_cap <= 255,
                "--mapq_cap must be in [0,255]"
            );
            ensure(
                cfg.coordmap.max_hops >= 1, 
                "--cm_max_hops must be >= 1"
            );
            ensure(
                cfg.coordmap.max_fanout >= 1, 
                "--cm_max_fanout must be >= 1"
            );
            ensure(
                cfg.coordmap.min_len >= 1, 
                "--cm_min_len must be >= 1"
            );
            ensure(
                cfg.coordmap.min_frac >= 0.0 && cfg.coordmap.min_frac <= 1.0, 
                "--cm_min_frac must be in [0,1]"
            );
            ensure(
                cfg.coordmap.max_total_hits >= 1, 
                "--cm_max_hits must be >= 1"
            );
            ensure(
                cfg.homq.sub_ovlp_frac >= 0.0 && cfg.homq.sub_ovlp_frac <= 1.0, 
                "--sub_ovlp_frac must be in [0,1]"
            );
            ensure(
                cfg.homq.close_as_eps > 0.0 && cfg.homq.close_as_eps <= 1.0, 
                "--close_as_eps must be in (0,1]"
            );
            ensure(
                cfg.global.threads >= 1,
                "-t/--threads must be >= 1"
            );
            ensure(
                cfg.global.IOthreads >= 1,
                "--io_threads must be >= 1"
            );
            if ((!cfg.ressit.genome_file.empty() || cfg.ressit.enzymes.size() > 0)) {
                ensure(!cfg.ressit.genome_file.empty(), "-g/--genome is required");
                ensure(!cfg.ressit.enzymes.empty(), "-e/--enzyme is required");
                for (size_t i = 0; i < cfg.ressit.enzymes.size(); ++i) {
                    const auto& ez = cfg.ressit.enzymes[i];
                    ensure(
                        !ez.empty(), 
                        "empty enzyme spec found"
                    );
                    ensure(
                        (std::count(ez.begin(), ez.end(), '^') == 1), 
                        "enzyme spec must contain exactly one '^': " + ez
                    );
                }
            }
            break;

        case ToolMode::align:
            ensure(
                !cfg.map.gfaFiles.empty(), 
                "-g/--gfa is required"
            );
            ensure(
                !cfg.map.reads.empty(),
                "-r/--read is required"
            );
            ensure(
                cfg.map.sec_pri_ratio > 0.0 && cfg.map.sec_pri_ratio <= 1.0, 
                "--secondary must be in (0,1]"
            );
            ensure(
                cfg.map.sec_pri_num >= 0, 
                "-N must be >= 0"
            );
            break;

        case ToolMode::ressit:
            ensure(
                !cfg.ressit.genome_file.empty(), 
                "-g/--genome is required"
            );
            ensure(
                !cfg.ressit.enzymes.empty(), 
                "-e/--enzyme is required"
            );
            for (size_t i = 0; i < cfg.ressit.enzymes.size(); ++i) {
                const auto& ez = cfg.ressit.enzymes[i];
                ensure(
                    !ez.empty(), 
                    "empty enzyme spec found"
                );
                ensure(
                    (std::count(ez.begin(), ez.end(), '^') == 1), 
                    "enzyme spec must contain exactly one '^': " + ez
                );
            }
            break;
    }
}

void help(char** argv, bool advanced) {
    std::cerr 
        << "Usage: " << argv[0] << " <subcommand> [options]\n\n"
        << program::description << "\n"
        << "Version: " << program::version << "\n"
        << "Date:    " << program::build_date << "\n";

    if (advanced) {
        std::cerr << "Note:    " << program::version_note << "\n";
    }

    std::cerr << "\n"
        << "Subcommands:\n"
        << "  stat        collect statistics about a GFA file\n"
        << "  seq         extract sequences along paths from a GFA\n"
        << "  gfa2fa      export all segments to FASTA with two-sided extension for short nodes\n"
        << "  depth       compute depth information from sequencing data for a GFA\n"
        << "  bubble      detect bubbles in a GFA\n"
        << "  collapse    de-overlap and collapse homologous sequences into single nodes\n"
        << "  file2map    convert GFA/PAF to map format for liftover\n"
        << "  liftover    liftover coordinates using a coordinate map\n"
        << "  mapq_boost  raise MAPQ for alignments in homologous regions using a coordinate map\n";

    if (advanced) {
        std::cerr
            << "  deoverlap   convert an overlap-based GFA to a non-overlap GFA\n"
            << "  align       align sequencing data to a GFA (beta)\n"
            << "  ressit      find restriction sites in a genome\n";
    }

    std::cerr << std::endl
        << "General Options:\n"
        << "  -h, --help          show basic options\n"
        << "  -H, --advanced      show advanced options\n\n";
}

void help_stat(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ...\n\n"
        << "Collect statistics about GFA file\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa    FILE ...    input GFA file(s)\n\n"
        << "General Options:\n"
        << "  -h, --help               show basic options\n\n";
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
            case 'g': {
                cfg.stat.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.stat.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'h': {
                help_stat(argv);
                std::exit(0);
            }
            default: {
                help_stat(argv);
                std::exit(1);
            }
        }
    }
    
    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_seq(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... -p FILE [options]\n\n"
        << "Extract sequences for compact marker paths (one path per line)\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa       FILE ...    input GFA file(s)\n"
        << "  -p, --path      FILE        text file of paths (one per line)\n"
        << "                               * path format: '>S1<S2>S3' or 'S1+,S2-,S3+'\n"
        << "  -o, --output    FILE        output file, fasta format [stdout]\n\n"
        << "General Options:\n"
        << "  -h, --help                  show basic options\n\n";
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
            case 'g': {
                cfg.seq.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.seq.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'p': {
                cfg.seq.pathFile = optarg;
                break;
            }
            case 'o': {
                cfg.seq.outFile = optarg;
                break;
            }
            case 'h': {
                help_seq(argv);
                std::exit(0);
            }
            default: {
                help_seq(argv);
                std::exit(1);
            }
        }
    }
    
    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_gfa2fa(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Export all segments of a GFA to FASTA, extending both ends by Y bp\n"
        << "for segments whose effective length is shorter than X bp.\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa          FILE ...    input GFA file(s)\n"
        << "  -o, --output       FILE        output FASTA file [stdout]\n\n"
        << "Extension options:\n"
        << "      --min_len      INT         only extend nodes with len < X [" << Gfa2FaOpts().min_len_xbp << "]\n"
        << "      --extend       INT         extend Y bp on both sides [" << Gfa2FaOpts().extend_ybp << "]\n"
        << "      --wrap         INT         FASTA wrap width; 0 = no wrap [" << Gfa2FaOpts().wrap_width << "]\n"
        << "      --keep_unknown             include nodes with '*' or zero length\n"
        << "                                  * by default such nodes are skipped\n\n"
        << "General Options:\n"
        << "  -h, --help                     show basic options\n\n";
}

AppConfig main_gfa2fa(int argc, char** argv) {
    if (argc < 3) { help_gfa2fa(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::gfa2fa;

    const struct option long_opts[] = {
        {"gfa",         required_argument, nullptr, 'g'},
        {"output",      required_argument, nullptr, 'o'},
        {"min_len",     required_argument, nullptr, 1001},
        {"extend",      required_argument, nullptr, 1002},
        {"wrap",        required_argument, nullptr, 1003},
        {"keep_unknown",no_argument,       nullptr, 1004},
        {"help",        no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:o:h";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.gfa2fa.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.gfa2fa.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'o': {
                cfg.gfa2fa.outFile = optarg;
                break;
            }
            case 1001: {
                cfg.gfa2fa.min_len_xbp = std::stoi(optarg);
                break;
            }
            case 1002: {
                cfg.gfa2fa.extend_ybp = std::stoi(optarg);
                break;
            }
            case 1003: {
                cfg.gfa2fa.wrap_width = std::stoi(optarg);
                break;
            }
            case 1004: {  // include '*' or zero-len nodes
                cfg.gfa2fa.skip_unknown = false;
                break;
            }
            case 'h': {
                help_gfa2fa(argv);
                std::exit(0);
            }
            default: {
                help_gfa2fa(argv);
                std::exit(1);
            }
        }
    }

    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_depth(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Compute node depth information from GFA (A-lines) / reads / GAF for a GFA graph\n\n"
        << "Three ways to compute depth:\n"
        << "  1. k-mer mode (provide -r):\n"
        << "       " << argv[0] << " " << argv[1] << " -g graph.gfa -r reads.fq -o graph.depth\n"
        << "  2. alignment mode (provide --gaf):\n"
        << "       " << argv[0] << " " << argv[1] << " -g graph.gfa --gaf aln.gaf -o graph.depth\n"
        << "  3. GFA-only mode (depth is calculated from 'A' lines in the GFA file):\n"
        << "       " << argv[0] << " " << argv[1] << " -g graph.gfa -o graph.depth\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa         FILE ...    input GFA file(s)\n"
        << "  -r, --read        FILE        input sequencing data file(s) for k-mer mode\n"
        << "                                 * depth calculated from k-mer coverage\n"
        << "      --gaf         FILE        input GAF alignment file for alignment mode\n"
        << "                                 * depth calculated directly from alignment\n"
        << "  -o, --output      FILE        output depth information to file [stdout]\n"
        << "      --base_depth              output per-base depth\n\n"
        << "Index options:\n"
        << "  -k, --kmer        INT         k-mer length [" << GlobalOpts().kmerLen << "]\n\n"
        << "GFA loading options:\n"
        << "  --min_mapq        INT         minimum mapping quality to keep [" << DepthOpts().min_mapq << "]\n"
        << "  --min_frac        FLOAT       minimum aligned-length / read-length ratio to keep [" << DepthOpts().min_frac << "]\n\n"
        << "General Options:\n"
        << "  -t, --threads     INT         number of threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug                   debug mode (forces threads=1)\n"
        << "  -h, --help                    show basic options\n\n";
}

AppConfig main_depth(int argc, char** argv) {
    if (argc < 3) { help_depth(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::depth;

    const struct option long_opts[] = {
        {"gfa",        required_argument, nullptr, 'g'},
        {"read",       required_argument, nullptr, 'r'},
        {"gaf",        required_argument, nullptr, 0001},
        {"output",     required_argument, nullptr, 'o'},
        {"base_depth", no_argument,       nullptr, 0002},

        {"kmer",       required_argument, nullptr, 'k'},

        {"min_mapq",   required_argument, nullptr, 2001},
        {"min_frac",   required_argument, nullptr, 2002},

        {"threads",    required_argument, nullptr, 't'},
        {"debug",      no_argument,       nullptr, 'd'},
        {"help",       no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:r:o:k:t:dh";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.depth.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.depth.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'r': {
                cfg.depth.reads.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.depth.reads.emplace_back(argv[optind]); ++optind;
                }
                break;
            }
            case 0001: {
                cfg.depth.gafFile = optarg;
                break;
            }
            case 'o': {
                cfg.depth.outFile = optarg;
                break;
            }
            case 0002: {
                cfg.depth.base_depth = true;
                break;
            }

            case 'k': {
                cfg.global.kmerLen = std::stoi(optarg);
                break;
            }

            case 2001: {
                cfg.depth.min_mapq = std::max(0, std::atoi(optarg));
                break;
            }
            case 2002: {
                cfg.depth.min_frac = std::max(0.0, std::atof(optarg));
                break;
            }

            case 't': {
                cfg.global.threads = std::max(1, std::stoi(optarg));
                break;
            }
            case 'd': {
                cfg.global.debug = true;
                break;
            }
            case 'h': {
                help_depth(argv);
                std::exit(0);
            }
            default: {
                help_depth(argv);
                std::exit(1);
            }
        }
    }
    
    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_bubble(char** argv, bool advanced) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Detect bubbles in GFA graph and write them to GFA/VCF file\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa             FILE ...    input GFA file(s)\n"
        << "  -p, --prefix          STR         output prefix [" << BubbleOpts().out_prefix << "]\n"
        << "                                     * prefix.bubbles.gfa, prefix.bubbles.noseq.gfa (and prefix.bubbles.vcf if --write_vcf is set)\n"
        << "      --write_vcf                   also output bubbles in VCF format\n"
        << "      --ref             FILE        reference FASTA used to extract REF sequence for VCF\n"
        << "      --paf             FILE        PAF file for liftover of contig coordinates to reference\n"
        << "                                     * must contain cg tag; e.g. minimap2 -cx asm5 ref.fa ctg.fa > aln.paf\n"
        << "      --ali_min_mapq    INT         minimum mapping quality for PAF liftover [" << BubbleOpts().ali_min_mapq << "]\n"
        << "      --ali_min_len     INT         minimum alignment length for PAF liftover [" << BubbleOpts().ali_min_len << "]\n\n"
        << "Detection Options:\n"
        << "      --depth           INT         maximum DFS depth for path exploration inside a bubble [" << BubbleOpts().max_depth << "]\n"
        << "      --paths           INT         maximum number of DFS paths to explore per bubble [" << BubbleOpts().max_paths << "]\n";
    if (advanced) {
        std::cerr
            << "      --DFS_guard       INT         max DFS states [" << BubbleOpts().DFS_guard << "]\n"
            << "      --cx_depth        INT         BFS depth limit for source complexity check; 0 disables this check [" << BubbleOpts().cx_depth << "]\n"
            << "      --cx_nodes        INT         max number of nodes visited in local BFS; exceeding marks source as complex [" << BubbleOpts().cx_nodes << "]\n"
            << "      --cx_branches     INT         max number of branching nodes allowed; exceeding marks source as complex [" << BubbleOpts().cx_branches << "]\n"
            << "      --cx_deg_branch   INT         node degree threshold to be considered branching [" << BubbleOpts().cx_deg_branch << "]\n"
            << "      --cx_deg_hub      INT         node degree threshold to be considered a hub (immediately complex) [" << BubbleOpts().cx_deg_hub << "]\n";
    }
    std::cerr
        << "      --path_diff       FLOAT       node difference ratio threshold for haplotype separation [" << BubbleOpts().path_diff << "]\n"
        << "      --stall_rounds    INT         stop DFS path search after this many rounds without a new unique path [" << BubbleOpts().stall_round_limit << "]\n"
        << "      --min_len         INT         minimum total sequence length of a bubble for output (0 = no filter) [" << BubbleOpts().min_len << "]\n"
        << "      --min_num         INT         minimum number of nodes inside a bubble required for output (0 = no filter) [" << BubbleOpts().min_num << "]\n"
        << "      --keep_nested                 keep bubbles contained within larger bubbles\n\n";
    if (advanced) {
        std::cerr
            << "Homologous Path Options:\n"
            << "      --same_sim        FLOAT       minimum similarity for homologous paths from the same source [" << BubbleOpts().same_sim << "]\n"
            << "      --same_num        INT         minimum number of nodes in each same-source homologous path [" << BubbleOpts().same_min_num << "]\n"
            << "      --same_len        INT         minimum total length of each same-source homologous path [" << BubbleOpts().same_min_len << "]\n"
            << "      --diff_src_len    INT         minimum source length for searching homologous paths between different sources [" << BubbleOpts().diff_min_src << "]\n"
            << "      --diff_sim        FLOAT       minimum similarity for homologous paths between different sources [" << BubbleOpts().diff_sim << "]\n"
            << "      --diff_num        INT         minimum number of nodes in each different-source homologous path [" << BubbleOpts().diff_min_num << "]\n"
            << "      --diff_len        INT         minimum total length of each different-source homologous path [" << BubbleOpts().diff_min_len << "]\n"
            << "      --hnum            INT         max paths per source branch [" << BubbleOpts().homo_num << "]\n"
            << "      --hk              INT         k-mer size for minimizer used in homologous path search [" << BubbleOpts().homo_k << "]\n"
            << "      --hw              INT         window size for minimizer used in homologous path search [" << BubbleOpts().homo_w << "]\n"
            << "      --hboot           INT         bootstrap bp before first homologous-path similarity check [" << BubbleOpts().homo_boot_bp << "]\n"
            << "      --hcnt            INT         minimum node length counted as one homologous-path extension step [" << BubbleOpts().homo_cnt_len << "]\n"
            << "      --hbits           INT         Bloom filter size in bits for homologous-path sketches [" << BubbleOpts().homo_bloom_bits << "]\n"
            << "      --hhash           INT         number of Bloom filter hash functions for homologous-path sketches [" << BubbleOpts().homo_bloom_hash << "]\n\n";

    }
    std::cerr
        << "General Options:\n"
        << "  -t, --threads         INT         number of threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug                       debug mode (forces threads=1)\n"
        << "  -h, --help                        show basic options\n"
        << "  -H, --advanced                    show advanced options\n\n";
}

AppConfig main_bubble(int argc, char** argv) {
    if (argc < 3) { help_bubble(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::bubble;

    const struct option long_opts[] = {
        {"gfa",             required_argument, nullptr, 'g'},
        {"prefix",          required_argument, nullptr, 'p'},
        {"write_vcf",       no_argument,       nullptr, 0001},
        {"ref",             required_argument, nullptr, 0002},
        {"paf",             required_argument, nullptr, 0003},
        {"ali_min_mapq",    required_argument, nullptr, 0004},
        {"ali_min_len",     required_argument, nullptr, 0005},

        {"depth",           required_argument, nullptr, 1001},
        {"paths",           required_argument, nullptr, 1002},
        {"DFS_guard",       required_argument, nullptr, 1003},
        {"cx_depth",        required_argument, nullptr, 1004},
        {"cx_nodes",        required_argument, nullptr, 1005},
        {"cx_branches",     required_argument, nullptr, 1006},
        {"cx_deg_branch",   required_argument, nullptr, 1007},
        {"cx_deg_hub",      required_argument, nullptr, 1008},
        {"path_diff",       required_argument, nullptr, 1009},
        {"stall_rounds",    required_argument, nullptr, 1010},
        {"min_len",         required_argument, nullptr, 1011},
        {"min_num",         required_argument, nullptr, 1012},
        {"keep_nested",     no_argument,       nullptr, 1013},

        {"same_sim",        required_argument, nullptr, 2001},
        {"same_num",        required_argument, nullptr, 2002},
        {"same_len",        required_argument, nullptr, 2003},
        {"diff_src_len",    required_argument, nullptr, 2004},
        {"diff_sim",        required_argument, nullptr, 2005},
        {"diff_num",        required_argument, nullptr, 2006},
        {"diff_len",        required_argument, nullptr, 2007},
        {"hnum",            required_argument, nullptr, 2008},
        {"hk",              required_argument, nullptr, 2009},
        {"hw",              required_argument, nullptr, 2010},
        {"hboot",           required_argument, nullptr, 2011},
        {"hcnt",            required_argument, nullptr, 2012},
        {"hbits",           required_argument, nullptr, 2013},
        {"hhash",           required_argument, nullptr, 2014},


        {"threads",         required_argument, nullptr, 't'},
        {"debug",           no_argument,       nullptr, 'd'},
        {"help",            no_argument,       nullptr, 'h'},
        {"advanced",        no_argument,       nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "g:p:t:dhH";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.bubble.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.bubble.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'p': {
                cfg.bubble.out_prefix = optarg;
                break;
            } 
            case 0001: {
                cfg.bubble.write_vcf = true;
                break;
            }
            case 0002: {
                cfg.bubble.ref_file = optarg;
                break;
            }
            case 0003: {
                cfg.bubble.paf_file = optarg;
                break;
            }
            case 0004: {
                cfg.bubble.ali_min_mapq = (uint32_t)std::stoul(optarg);
                break;
            }
            case 0005: {
                cfg.bubble.ali_min_len = (uint32_t)std::stoul(optarg);
                break;
            }

            case 1001: {
                cfg.bubble.max_depth = std::stoi(optarg);
                break;
            }
            case 1002: {
                cfg.bubble.max_paths = std::stoi(optarg);
                break;
            }
            case 1003: {
                cfg.bubble.DFS_guard = (uint64_t)std::stoull(optarg);
                break;
            }
            case 1004: {
                cfg.bubble.cx_depth = (uint16_t)std::stoul(optarg);
                break;
            }
            case 1005: {
                cfg.bubble.cx_nodes = (uint32_t)std::stoul(optarg);
                break;
            }
            case 1006: {
                cfg.bubble.cx_branches = (uint32_t)std::stoul(optarg);
                break;
            }
            case 1007: {
                cfg.bubble.cx_deg_branch = std::stoi(optarg);
                break;
            }
            case 1008: {
                cfg.bubble.cx_deg_hub = std::stoi(optarg);
                break;
            }
            case 1009: {
                cfg.bubble.path_diff = std::stod(optarg);
                break;
            }
            case 1010: {
                cfg.bubble.stall_round_limit = (uint32_t)std::max(0, std::stoi(optarg));
                break;
            }
            case 1011: {
                cfg.bubble.min_len = std::stoi(optarg);
                break;
            }
            case 1012: {
                cfg.bubble.min_num = std::stoi(optarg);
                break;
            }
            case 1013: {
                cfg.bubble.keep_nested = true;
                break;
            }

            case 2001: {
                cfg.bubble.same_sim = std::stod(optarg);
                break;
            }
            case 2002: {
                cfg.bubble.same_min_num = (uint32_t)std::stoul(optarg);
                break;
            }
            case 2003: {
                cfg.bubble.same_min_len = (uint64_t)std::stoull(optarg);
                break;
            }
            case 2004: {
                cfg.bubble.diff_min_src = (uint32_t)std::stoul(optarg);
                break;
            }
            case 2005: {
                cfg.bubble.diff_sim = std::stod(optarg);
                break;
            }
            case 2006: {
                cfg.bubble.diff_min_num = (uint32_t)std::stoul(optarg);
                break;
            }
            case 2007: {
                cfg.bubble.diff_min_len = (uint64_t)std::stoull(optarg);
                break;
            }
            case 2008: {
                cfg.bubble.homo_num = std::max((uint32_t)1, (uint32_t)std::stoul(optarg));
                break;
            }
            case 2009: {
                cfg.bubble.homo_k = (uint32_t)std::max(1, std::stoi(optarg));
                break;
            }
            case 2010: {
                cfg.bubble.homo_w = (uint32_t)std::max(1, std::stoi(optarg));
                break;
            }
            case 2011: {
                cfg.bubble.homo_boot_bp = (uint32_t)std::stoul(optarg);
                break;
            }
            case 2012: {
                cfg.bubble.homo_cnt_len = (uint32_t)std::stoul(optarg);
                break;
            }
            case 2013: {
                cfg.bubble.homo_bloom_bits = (uint64_t)std::stoull(optarg);
                break;
            }
            case 2014: {
                cfg.bubble.homo_bloom_hash = (uint32_t)std::stoul(optarg);
                break;
            }


            case 't': {
                cfg.global.threads = std::max(1, std::stoi(optarg));
                break;
            }
            case 'd': {
                cfg.global.debug = true;
                break;
            }
            case 'h': {
                help_bubble(argv);
                std::exit(0);
            }
            case 'H': {
                help_bubble(argv, true);
                std::exit(0);
            }
            default: {
                help_bubble(argv);
                std::exit(1);
            }
        }
    }

    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_deoverlap(char** argv, bool advanced) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Convert an overlap-based GFA into a non-overlap GFA\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa              FILE ...    input GFA file(s)\n"
        << "  -p, --prefix           STR         output prefix [" << CollapseOpts().prefix << "]\n"
        << "                                      * prefix.deoverlap.gfa, prefix.deoverlap.noseq.gfa and prefix.deoverlap.map\n\n"
        << "Index options:\n"
        << "  -k, --kmer             INT         k-mer length [" << GlobalOpts().kmerLen << "]\n"
        << "  -w, --window           INT         minimizer window [" << GlobalOpts().minimizerW << "]\n\n"
        << "Align options:\n";
        if (advanced) {
            std::cerr
                << "      --use_wfa                      whether to use WFA for alignment, if not use mm2 (beta)\n";
        }
    std::cerr
        << "  -z, --zdrop            INT         Z-drop score [" << MapOpts().extendOpts.dyn_zdrop << "]\n"
        << "      --min_match        FLOAT       minimum CIGAR match ratio (keep alignment only if match fraction ≥ threshold) [" << CollapseOpts().min_match_ratio << "]\n\n"
        << "Collapse options:\n"
        << "      --min_eq           INT         minimum match length to add cut points at segment ends [" << CollapseOpts().min_eq << "]\n"
        << "      --max_iters        INT         maximum iterations for cut point propagation [" << CollapseOpts().max_iters << "]\n";
    if (advanced) {
        std::cerr
            << "      --abnormal_cut     INT,INT     prune abnormal cut points after propagation as MAX_LEN,MIN_COUNT [" << CollapseOpts().max_abnormal_cut_len << "," << CollapseOpts().min_abnormal_cut_count << "]\n"
            << "                                      * e.g. 0-3-6-9-100 -> 0-9-100 when 4,2 due to consecutive short segments\n";
    }
    std::cerr << "\n"
        << "General Options:\n"
        << "  -t, --threads          INT         number of threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug                        debug mode\n"
        << "  -h, --help                         show basic options\n\n";
}

AppConfig main_deoverlap(int argc, char** argv) {
    if (argc < 3) { help_deoverlap(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::deoverlap;

    const struct option long_opts[] = {
        {"gfa",          required_argument, nullptr, 'g'},
        {"prefix",       required_argument, nullptr, 'p'},

        {"kmer",         required_argument, nullptr, 'k'},
        {"window",       required_argument, nullptr, 'w'},

        {"use_wfa",      no_argument,       nullptr, 2001},
        {"zdrop",        required_argument, nullptr, 'z'},
        {"min_match",    required_argument, nullptr, 2002},

        {"min_eq",       required_argument, nullptr, 3002},
        {"max_iters",    required_argument, nullptr, 3003},
        {"abnormal_seg", required_argument, nullptr, 3004},

        {"threads",      required_argument, nullptr, 't'},
        {"debug",        no_argument,       nullptr, 'd'},
        {"help",         no_argument,       nullptr, 'h'},
        {"advanced",     no_argument,       nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "g:p:k:w:z:t:dhH";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.collapse.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.collapse.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'p': {
                cfg.collapse.prefix = optarg;
                break;
            }

            case 'k': {
                cfg.global.kmerLen = std::max(1, std::stoi(optarg));
                break;
            }
            case 'w': {
                cfg.global.minimizerW = std::max(1, std::stoi(optarg));
                break;
            }

            case 2001: {
                cfg.map.use_wfa = true;
                break;
            }
            case 'z': {
                cfg.map.zdrop = std::max(1, std::stoi(optarg));
                cfg.map.zdrop_set = true;
                break;
            }
            case 2002: {
                cfg.collapse.min_match_ratio = std::max(0.0, std::stod(optarg));
                break;
            }

            case 3002: {
                cfg.collapse.min_eq = std::stoi(optarg);
                break;
            }
            case 3003: {
                cfg.collapse.max_iters = std::stoi(optarg);
                break;
            }
            case 3004: {
                parse_abnormal_seg_(optarg, cfg.collapse);
                break;
            }

            case 't': {
                cfg.global.threads = std::max(1, std::stoi(optarg));
                break;
            }
            case 'd': {
                cfg.global.debug = true; break;
            }
            case 'h': {
                help_deoverlap(argv);
                std::exit(0);
            }
            case 'H': {
                help_deoverlap(argv, true);
                std::exit(0);
            }
            default : {
                help_deoverlap(argv);
                std::exit(1);
            }
        }
    }
    
    validate_and_print(argc, argv, cfg);
    finalize_opt_cfg(cfg);

    return cfg;
}

void help_collapse(char** argv, bool advanced) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "De-overlap (if needed) and collapse homologous sequences in GFA into single nodes\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa              FILE ...     input GFA file(s) (overlap or no-overlap)\n"
        << "  -p, --prefix           STR          output prefix [" << CollapseOpts().prefix << "]\n"
        << "                                       * prefix.collapse.gfa, prefix.collapse.noseq.gfa and prefix.collapse.map\n\n"
        << "Index options:\n"
        << "  -k, --kmer             INT          k-mer length [" << GlobalOpts().kmerLen << "]\n"
        << "  -w, --window           INT          minimizer window [" << GlobalOpts().minimizerW << "]\n\n"
        << "Align options:\n";
        if (advanced) {
            std::cerr
                << "      --use_wfa                       whether to use WFA for alignment, if not use mm2 (beta)\n";
        }
    std::cerr
        << "  -z, --zdrop            INT          Z-drop score [" << MapOpts().extendOpts.dyn_zdrop << "]\n"
        << "      --min_match        FLOAT        minimum CIGAR match ratio (keep alignment only if match fraction ≥ threshold) [" << CollapseOpts().min_match_ratio << "]\n\n"
        << "Collapse options:\n"
        << "      --min_jaccard      FLOAT        nodes with similarity above this value will be collapsed [" << CollapseOpts().min_jaccard << "]\n"
        << "      --min_eq           INT          minimum match length to add cut points at segment ends [" << CollapseOpts().min_eq << "]\n"
        << "      --max_iters        INT          maximum iterations for cut point propagation [" << CollapseOpts().max_iters << "]\n";
        if (advanced) {
            std::cerr
                << "      --repeat_mask      INT,INT,INT  tandem-repeat masking as MIN_LEN,MAX_PERIOD,MAX_MISMATCH [" << CollapseOpts().repeat_mask_min_len << "," << CollapseOpts().repeat_mask_max_period << "," << CollapseOpts().repeat_mask_max_mismatch << "]\n"
                << "                                       * e.g. ATATATGATCG contains a 9 bp repeat with period 2 (AT) and 1 mismatch (G), detectable by 8,12,1\n"
                << "      --abnormal_cut     INT,INT      prune abnormal cut points after propagation as MAX_LEN,MIN_COUNT [" << CollapseOpts().max_abnormal_cut_len << "," << CollapseOpts().min_abnormal_cut_count << "]\n"
                << "                                       * e.g. 0-3-6-9-100 -> 0-9-100 when 4,2 due to consecutive short segments\n";
        }
    std::cerr << "\n"
        << "Bubble Detection Options:\n"
        << "      --depth            INT          maximum DFS depth for path exploration inside a bubble [" << BubbleOpts().max_depth << "]\n"
        << "      --paths            INT          maximum number of DFS paths to explore per bubble [" << BubbleOpts().max_paths << "]\n";
    if (advanced) {
        std::cerr
            << "      --DFS_guard        INT          max DFS states [" << BubbleOpts().DFS_guard << "]\n"
            << "      --cx_depth         INT          BFS depth limit for source complexity check; 0 disables this check [" << BubbleOpts().cx_depth << "]\n"
            << "      --cx_nodes         INT          max number of nodes visited in local BFS; exceeding marks source as complex [" << BubbleOpts().cx_nodes << "]\n"
            << "      --cx_branches      INT          max number of branching nodes allowed; exceeding marks source as complex [" << BubbleOpts().cx_branches << "]\n"
            << "      --cx_deg_branch    INT          node degree threshold to be considered branching [" << BubbleOpts().cx_deg_branch << "]\n"
            << "      --cx_deg_hub       INT          node degree threshold to be considered a hub (immediately complex) [" << BubbleOpts().cx_deg_hub << "]\n";
    }
    std::cerr
        << "      --path_diff        FLOAT        node difference ratio threshold for haplotype separation [" << BubbleOpts().path_diff << "]\n"
        << "      --stall_rounds     INT          stop DFS path search after this many rounds without a new unique path [" << BubbleOpts().stall_round_limit << "]\n\n";
    if (advanced) {
        std::cerr
            << "Homologous Path Detection Options:\n"
            << "      --same_sim         FLOAT        minimum similarity for homologous paths from the same source [" << BubbleOpts().same_sim << "]\n"
            << "      --same_num         INT          minimum number of nodes in each same-source homologous path [" << BubbleOpts().same_min_num << "]\n"
            << "      --same_len         INT          minimum total length of each same-source homologous path [" << BubbleOpts().same_min_len << "]\n"
            << "      --diff_src_len     INT          minimum source length for searching homologous paths between different sources [" << BubbleOpts().diff_min_src << "]\n"
            << "      --diff_sim         FLOAT        minimum similarity for homologous paths between different sources [" << BubbleOpts().diff_sim << "]\n"
            << "      --diff_num         INT          minimum number of nodes in each different-source homologous path [" << BubbleOpts().diff_min_num << "]\n"
            << "      --diff_len         INT          minimum total length of each different-source homologous path [" << BubbleOpts().diff_min_len << "]\n"
            << "      --hnum             INT          max paths per source branch [" << BubbleOpts().homo_num << "]\n"
            << "      --hk               INT          k-mer size for minimizer used in homologous path search [" << BubbleOpts().homo_k << "]\n"
            << "      --hw               INT          window size for minimizer used in homologous path search [" << BubbleOpts().homo_w << "]\n"
            << "      --hboot            INT          bootstrap bp before first homologous-path similarity check [" << BubbleOpts().homo_boot_bp << "]\n"
            << "      --hcnt             INT          minimum node length counted as one homologous-path extension step [" << BubbleOpts().homo_cnt_len << "]\n"
            << "      --hbits            INT          Bloom filter size in bits for homologous-path sketches [" << BubbleOpts().homo_bloom_bits << "]\n"
            << "      --hhash            INT          number of Bloom filter hash functions for homologous-path sketches [" << BubbleOpts().homo_bloom_hash << "]\n\n";

    }
    std::cerr
        << "General Options:\n"
        << "  -t, --threads          INT          number of threads (~5GB RAM per thread) [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug                         debug mode (forces threads=1)\n"
        << "  -h, --help                          show basic options\n"
        << "  -H, --advanced                      show advanced options\n\n";
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
        {"zdrop",           required_argument, nullptr, 'z'},
        {"min_match",       required_argument, nullptr, 2002},

        {"min_jaccard",     required_argument, nullptr, 3001},
        {"min_eq",          required_argument, nullptr, 3002},
        {"max_iters",       required_argument, nullptr, 3003},
        {"repeat_mask",     required_argument, nullptr, 3004},
        {"abnormal_seg",    required_argument, nullptr, 3005},

        {"depth",           required_argument, nullptr, 4001},
        {"paths",           required_argument, nullptr, 4002},
        {"DFS_guard",       required_argument, nullptr, 4003},
        {"cx_depth",        required_argument, nullptr, 4004},
        {"cx_nodes",        required_argument, nullptr, 4005},
        {"cx_branches",     required_argument, nullptr, 4006},
        {"cx_deg_branch",   required_argument, nullptr, 4007},
        {"cx_deg_hub",      required_argument, nullptr, 4008},
        {"path_diff",       required_argument, nullptr, 4009},
        {"stall_rounds",    required_argument, nullptr, 4010},

        {"same_sim",        required_argument, nullptr, 5001},
        {"same_num",        required_argument, nullptr, 5002},
        {"same_len",        required_argument, nullptr, 5003},
        {"diff_src_len",    required_argument, nullptr, 5004},
        {"diff_sim",        required_argument, nullptr, 5005},
        {"diff_num",        required_argument, nullptr, 5006},
        {"diff_len",        required_argument, nullptr, 5007},
        {"hnum",            required_argument, nullptr, 5108},
        {"hk",              required_argument, nullptr, 5109},
        {"hw",              required_argument, nullptr, 5110},
        {"hboot",           required_argument, nullptr, 5011},
        {"hcnt",            required_argument, nullptr, 5012},
        {"hbits",           required_argument, nullptr, 5013},
        {"hhash",           required_argument, nullptr, 5014},


        {"threads",         required_argument, nullptr, 't'},
        {"debug",           no_argument,       nullptr, 'd'},
        {"help",            no_argument,       nullptr, 'h'},
        {"advanced",        no_argument,       nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "g:p:k:w:z:t:dhH";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.collapse.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.collapse.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'p': {
                cfg.collapse.prefix = optarg;
                break;
            }

            case 'k': {
                cfg.global.kmerLen = std::max(1, std::stoi(optarg));
                break;
            }
            case 'w': {
                cfg.global.minimizerW = std::max(1, std::stoi(optarg));
                break;
            }

            case 2001: {
                cfg.map.use_wfa = true;
                break;
            }
            case 'z': {
                cfg.map.zdrop = std::max(1, std::stoi(optarg));
                cfg.map.zdrop_set = true;
                break;
            }
            case 2002: {
                cfg.collapse.min_match_ratio = std::max(0.0, std::stod(optarg));
                break;
            }

            case 3001: {
                cfg.collapse.min_jaccard = std::stod(optarg);
                break;
            }
            case 3002: {
                cfg.collapse.min_eq = std::stoi(optarg);
                break;
            }
            case 3003: {
                cfg.collapse.max_iters = std::stoi(optarg);
                break;
            }
            case 3004: {
                parse_repeat_mask_(optarg, cfg.collapse);
                break;
            }
            case 3005: {
                parse_abnormal_seg_(optarg, cfg.collapse);
                break;
            }

            case 4001: {
                cfg.bubble.max_depth = std::stoi(optarg);
                break;
            }
            case 4002: {
                cfg.bubble.max_paths = std::stoi(optarg);
                break;
            }
            case 4003: {
                cfg.bubble.DFS_guard = (uint64_t)std::stoull(optarg);
                break;
            }
            case 4004: {
                cfg.bubble.cx_depth = (uint16_t)std::stoul(optarg);
                break;
            }
            case 4005: {
                cfg.bubble.cx_nodes = (uint32_t)std::stoul(optarg);
                break;
            }
            case 4006: {
                cfg.bubble.cx_branches = (uint32_t)std::stoul(optarg);
                break;
            }
            case 4007: {
                cfg.bubble.cx_deg_branch = std::stoi(optarg);
                break;
            }
            case 4008: {
                cfg.bubble.cx_deg_hub = std::stoi(optarg);
                break;
            }
            case 4009: {
                cfg.bubble.path_diff = std::stod(optarg);
                break;
            }
            case 4010: {
                cfg.bubble.stall_round_limit = (uint32_t)std::max(0, std::stoi(optarg));
                break;
            }

            case 5001: {
                cfg.bubble.same_sim = std::stod(optarg);
                break;
            }
            case 5002: {
                cfg.bubble.same_min_num = (uint32_t)std::stoul(optarg);
                break;
            }
            case 5003: {
                cfg.bubble.same_min_len = (uint64_t)std::stoull(optarg);
                break;
            }
            case 5004: {
                cfg.bubble.diff_min_src = (uint32_t)std::stoul(optarg);
                break;
            }
            case 5005: {
                cfg.bubble.diff_sim = std::stod(optarg);
                break;
            }
            case 5006: {
                cfg.bubble.diff_min_num = (uint32_t)std::stoul(optarg);
                break;
            }
            case 5007: {
                cfg.bubble.diff_min_len = (uint64_t)std::stoull(optarg);
                break;
            }
            case 5108: {
                cfg.bubble.homo_num = std::max((uint32_t)1, (uint32_t)std::stoul(optarg));
                break;
            }
            case 5109: {
                cfg.bubble.homo_k = (uint32_t)std::max(1, std::stoi(optarg));
                break;
            }
            case 5110: {
                cfg.bubble.homo_w = (uint32_t)std::max(1, std::stoi(optarg));
                break;
            }
            case 5011: {
                cfg.bubble.homo_boot_bp = (uint32_t)std::stoul(optarg);
                break;
            }
            case 5012: {
                cfg.bubble.homo_cnt_len = (uint32_t)std::stoul(optarg);
                break;
            }
            case 5013: {
                cfg.bubble.homo_bloom_bits = (uint64_t)std::stoull(optarg);
                break;
            }
            case 5014: {
                cfg.bubble.homo_bloom_hash = (uint32_t)std::stoul(optarg);
                break;
            }

            case 't': {
                cfg.global.threads = std::max(1, std::stoi(optarg));
                break;
            }
            case 'd': {
                cfg.global.debug = true;
                break;
            }
            case 'h': {
                help_collapse(argv);
                std::exit(0);
            }
            case 'H': {
                help_collapse(argv, true);
                std::exit(0);
            }
            default: {
                help_collapse(argv);
                std::exit(1);
            }
        }
    }
    
    validate_and_print(argc, argv, cfg);
    finalize_opt_cfg(cfg);

    return cfg;
}

void help_file2map(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -i FILE ... [options]\n\n"
        << "convert GFA/PAF to map format for liftover\n\n"
        << "PAF input notes:\n"
        << "  * For PAF input, generate PAF with minimap2 in 1-to-1 mode:\n"
        << "     - minimap2 -t8 --eqx -cx asm5 --secondary=no reference.fa query.fa > aln.paf\n"
        << "  * IMPORTANT: reference and query sequence names must be different.\n"
        << "     - For example, avoid using the same name like 'chr1' in both files.\n"
        << "       Rename them to something like 'chr1_ref' and 'chr1_qry'.\n"
        << "     - Sequence names must NOT contain ':', '-' or '+' characters,\n"
        << "       to avoid confusion with map coordinate formats such as 'chr:beg-end+'.\n\n"
        << "Input/Output:\n"
        << "  -i, --input        FILE ...    input PAF/GFA file(s)\n"
        << "  -o, --output       FILE        output file name [stdout]\n\n"
        << "PAF options:\n"
        << "      --all                      include secondary/supplementary, (default: only primary tp:A:P)\n"
        << "      --min_len      INT         minimum alignment length to keep [" << File2mapOpts().min_len << "]\n"
        << "      --min_mapq     INT         minimum mapping quality to keep [" << File2mapOpts().min_mapq << "]\n\n"
        << "General Options:\n"
        << "  -h, --help                     show basic options\n\n";
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
            case 'i': {
                cfg.file2map.inputFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.file2map.inputFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'o': {
                cfg.file2map.outFile = optarg;
                break;
            }
            case 1001: {
                cfg.file2map.paf_primary_only = false;
                break;
            }
            case 1002: {
                cfg.file2map.min_len = std::max(0, std::stoi(optarg));
                break;
            }
            case 1003: {
                cfg.file2map.min_mapq = std::atoi(optarg);
                break;
            }
            case 'h': {
                help_file2map(argv);
                std::exit(0);
            }
            default: {
                help_file2map(argv);
                std::exit(1);
            }
        }
    }

    if (cfg.file2map.outFile.empty()) cfg.file2map.outFile = "-";

    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_liftover(char** argv, bool advanced) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -m FILE ... -b FILE [options]\n\n"
        << "LiftOver BED coordinates using a .map file, which produced by deoverlap/collapse/file2map\n\n"
        << "Input/Output:\n"
        << "  -m, --map           FILE ...    coordinate mapping file(s) (generated by deoverlap/collapse/file2map)\n"
        << "      --paf           FILE        input PAF file to assist liftover (optional)\n"
        << "  -b, --bed           FILE        input BED file, required if not checking\n"
        << "  -o, --output        FILE        output file name [stdout]\n"
        << "  -r, --ref           FILE        FASTA/FASTQ containing sequences for src+dst contigs, required if checking\n\n"
        << "Map Options:\n"
        << "      --regex         STR         keep hits whose label matches regex\n"
        << "                                   * \"h\\d+tg\" to match h1tg, h2tg, ...\n"
        << "      --min_frac      FLOAT       minimum fraction of the BED interval that must be lifted over [" << LiftoverOpts().min_frac << "]\n"
        << "                                   * 0.9 means at least 90% of the BED region is mapped\n";
    if (advanced) {
        std::cerr
            << "      --flank_win     INT         extension window size around BED ends to identify insertions [" << LiftoverOpts().flank_win << "]\n"
            << "      --max_flank     INT         max extension around BED ends to identify insertions [" << LiftoverOpts().max_flank << "]\n"
            << "      --max_gap       INT         max allowed |ref_gap - query_gap| when calling INS [" << LiftoverOpts().max_gap << "]\n"
            << "      --max_hit       INT         max extended liftover results to keep [" << LiftoverOpts().max_hit << "]\n";
    }
    std::cerr
        << std::endl
        << "PAF Options:\n"
        << "      --min_len       INT         minimum alignment length for liftover [" << LiftoverOpts().min_len << "]\n"
        << "      --min_mapq      INT         minimum mapping quality for liftover [" << LiftoverOpts().min_mapq << "]\n\n"
        << "Check Options:\n"
        << "      --check                     run sequence consistency check, requires -r\n";
    if (advanced) {
        std::cerr
            << "  -w, --win           INT         window size [" << LiftoverOpts().win << "]\n"
            << "  -s, --step          INT         step size [" << LiftoverOpts().step << "]\n"
            << "      --max_examples  INT         max reported examples [" << LiftoverOpts().max_examples << "]\n";
    }  
    std::cerr 
        << std::endl
        << "CoordMap Options:\n"
        << "      --cm_max_hops   INT         max hop depth when chaining names [" << CoordMapOpts().max_hops << "]\n"
        << "                                   * e.g. 1 allows: source -> A; 2 allows: source -> A -> B\n";

    if (advanced) {
        std::cerr
            << "      --cm_max_fanout INT         max fanout per search [" << CoordMapOpts().max_fanout << "]\n"
            << "      --cm_min_len    INT         min length (bp) to keep a hit for hops [" << CoordMapOpts().min_len << "]\n"
            << "      --cm_min_frac   FLOAT       min fraction of current interval length to keep a hit [" << CoordMapOpts().min_frac << "]\n"
            << "                                   * threshold = max(min_len, cur_len * min_frac)\n"
            << "      --cm_max_hits   INT         global BFS node cap [" << CoordMapOpts().max_total_hits << "]\n";
    }

    std::cerr
        << std::endl
        << "General Options:\n"
        << "  -t, --threads       INT         number of threads [" << GlobalOpts().threads << "]\n"
        << "  -h, --help                      show basic options\n"
        << "  -H, --advanced                  show advanced options\n\n";
}

AppConfig main_liftover(int argc, char** argv) {
    if (argc < 3) { help_liftover(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::liftover;

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
        {"advanced",      no_argument,       nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "m:b:o:r:w:s:t:hH";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'm': {
                cfg.liftover.mapFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.liftover.mapFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 0001: {
                cfg.liftover.pafFile = optarg;
                break;
            }
            case 'b': {
                cfg.liftover.bedFile = optarg;
                break;
            }
            case 'o': {
                cfg.liftover.outFile = optarg;
                break;
            }
            case 'r': {
                cfg.liftover.referenceFile = optarg;
                break;
            }

            case 1001: {
                cfg.liftover.regex = optarg;
                break;
            }
            case 1002: {
                cfg.liftover.min_frac = std::max(0.0, std::min(1.0, std::atof(optarg)));
                break;
            }
            case 1003: {
                cfg.liftover.flank_win = (uint32_t)std::stoul(optarg);
                break;
            }
            case 1004: {
                cfg.liftover.max_flank = (uint32_t)std::stoul(optarg);
                break;
            }
            case 1005: {
                cfg.liftover.max_gap = (uint32_t)std::stoul(optarg);
                break;
            }
            case 1006: {
                cfg.liftover.max_hit = (uint16_t)std::stoul(optarg);
                break;
            }

            case 2001: {
                cfg.liftover.min_len = std::max(0, std::atoi(optarg));
                break;
            }
            case 2002: {
                cfg.liftover.min_mapq = std::atoi(optarg);
                break;
            }

            case 3001: {
                cfg.liftover.do_check = true;
                break;
            }
            case 'w': {
                cfg.liftover.win = (uint32_t)std::stoul(optarg);
                break;
            }
            case 's': {
                cfg.liftover.step = (uint32_t)std::stoul(optarg);
                break;
            }
            case 3002: {
                cfg.liftover.max_examples = (uint32_t)std::stoul(optarg);
                break;
            }

            case 4001: {
                cfg.coordmap.max_hops = std::max(1, std::atoi(optarg));
                break;
            }
            case 4002: {
                cfg.coordmap.max_fanout = (uint32_t)std::max(1, std::atoi(optarg));
                break;
            }
            case 4003: {
                cfg.coordmap.min_len = (uint32_t)std::max(1, std::atoi(optarg));
                break;
            }
            case 4004: {
                cfg.coordmap.min_frac = std::max(0.0, std::atof(optarg));
                break;
            }
            case 4005: {
                cfg.coordmap.max_total_hits = (uint32_t)std::max(1, std::atoi(optarg));
                break;
            }

            case 't': {
                cfg.global.threads = std::max(1, std::stoi(optarg));
                break;
            }
            case 'h': {
                help_liftover(argv);
                std::exit(0);
            }
            case 'H': {
                help_liftover(argv, true);
                std::exit(0);
            }
            default: {
                help_liftover(argv);
                std::exit(1);
            }
        }
    }

    if (cfg.liftover.outFile.empty())  cfg.liftover.outFile = "-";

    validate_and_print(argc, argv, cfg);
    return cfg;
}

void help_mapq_boost(char** argv, bool advanced) {
    static const char* kTail =
        "    | liftasm mapq_boost -m merge.map \\\n"
        "    | samtools view -b -q 10 -o out.boost.bam\n";

    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -m FILE ... [options]\n\n"
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
        << "  bwa mem HG002.hap1_hap2.fa read.1.fq.gz read.2.fq.gz \\\n"
        << kTail
        << "  # HISAT2\n"
        << "  hisat2 -x INDEX -1 read.1.fq.gz -2 read.2.fq.gz \\\n"
        << kTail
        << "  # STAR\n"
        << "  STAR --genomeDir ./ --readFilesIn read.1.fq.gz read.2.fq.gz --readFilesCommand zcat \\\n"
        << "       --outSAMtype BAM Unsorted --outStd BAM_Unsorted \\\n"
        << kTail
        << "  # minimap2 (Pore-C; large -N required due to many secondary alignments)\n"
        << "  minimap2 -ax map-ont -N 1000 HG002.hap1_hap2.fa pore-c.fq.gz \\\n"
        << kTail
        << std::endl
        << "Input/Output:\n"
        << "  -m, --map             FILE ...    coordinate mapping file(s) (generated by deoverlap/collapse/file2map)\n"
        << "  -i, --in              FILE        input SAM/BAM [stdin]\n"
        << "  -o, --out             FILE        output SAM/BAM [stdout]\n"
        << "  -g, --genome          FILE        input genome FASTA for restriction-site search (optional)\n"
        << "                                     * filter alignments for enzyme-cut data, e.g. Hi-C, Pore-C, CiFi\n\n"
        << "Boost Options:\n";
    if (advanced) {
        std::cerr << "      --batch           INT         batch size [" << MapqBoostOpts().batch_size << "]\n";
    }
    std::cerr
        << "      --mapq_low        INT         only raise if MAPQ<=Q (no larger than 255) [" << (int)MapqBoostOpts().mapq_low << "]\n"
        << "      --mapq_cap        INT         upper bound of boosted MAPQ (no larger than 255) [" << (int)MapqBoostOpts().mapq_cap << "]\n"
        << "      --name_check                  check whether QNAMEs are name-sorted (off by default)\n\n"
        << "CoordMap Options:\n"
        << "      --cm_max_hops     INT         BFS depth limit [" << CoordMapOpts().max_hops << "]\n"
        << "                                     * e.g. 1 allows: source -> A; 2 allows: source -> A -> B\n";
    if (advanced) {
        std::cerr
            << "      --cm_max_fanout   INT         max fanout per search [" << CoordMapOpts().max_fanout << "]\n"
            << "      --cm_min_len      INT         min length (bp) to keep a hit for hops [" << CoordMapOpts().min_len << "]\n"
            << "      --cm_min_frac     FLOAT       min fraction of current interval length to keep a hit [" << CoordMapOpts().min_frac << "]\n"
            << "                                     * threshold = max(min_len, cur_len * min_frac)\n"
            << "      --cm_max_hits     INT         global BFS node cap [" << CoordMapOpts().max_total_hits << "]\n";
    }
    std::cerr
        << std::endl
        << "Subgrouping Options:\n"
        << "      --sub_ovlp_frac   FLOAT       minimum overlap fraction on read to cluster alignments into the same subgroup [" << MapqBoostOpts().sub_ovlp_frac << "]\n"
        << "                                     * alignments overlapping the same read region are evaluated together\n";
    if (advanced) {
        std::cerr
            << "Scoring Options, used to rank alignments within a subgroup and select the best candidate for boosting:\n"
            << "      --K_mapq          FLOAT       MAPQ score scaling factor [" << MapqBoostOpts().K_mapq << "]\n"
            << "      --K_rs            FLOAT       restriction-site distance scaling factor [" << MapqBoostOpts().K_rs << "]\n"
            << "      --K_as            FLOAT       AS score scaling factor [" << MapqBoostOpts().K_as << "]\n"
            << "      --K_ml            FLOAT       match-length score scaling factor [" << MapqBoostOpts().K_ml << "]\n"
            << "      --K_nm            FLOAT       NM score scaling factor [" << MapqBoostOpts().K_nm << "]\n"
            << "      --W_mapq          FLOAT       MAPQ term weight [" << MapqBoostOpts().W_mapq << "]\n"
            << "      --W_rs            FLOAT       restriction-site term weight [" << MapqBoostOpts().W_rs << "]\n"
            << "      --W_as            FLOAT       AS term weight [" << MapqBoostOpts().W_as << "]\n"
            << "      --W_ml            FLOAT       match-length term weight [" << MapqBoostOpts().W_ml << "]\n"
            << "      --W_nm            FLOAT       NM penalty weight [" << MapqBoostOpts().W_nm << "]\n"
            << "      --close_as_eps    FLOAT       max relative AS difference to consider two alignments equally scored for MAPQ recalculation [" << MapqBoostOpts().close_as_eps << "]\n";
    }
    std::cerr
        << std::endl
        << "Restriction-site Options (optional):\n"
        << "  -e, --enzyme          STR         enzyme recognition site with cut mark '^', can take multiple arguments after one -e\n"
        << "                                     * e.g. G^AATTC  A^AGCTT  G^ANTC\n"
        << "      --rc                          scan reverse-complement motifs\n\n"
        << "General Options:\n"
        << "  -t, --threads         INT         number of threads [" << GlobalOpts().threads << "]\n"
        << "      --io_threads      INT         HTS I/O threads [" << GlobalOpts().IOthreads << "]\n"
        << "  -d, --debug                       debug mode (forces threads=1)\n"
        << "  -h, --help                        show basic options\n"
        << "  -H, --advanced                    show advanced options\n\n";
}

AppConfig main_mapq_boost(int argc, char** argv) {
    if (argc < 3) { help_mapq_boost(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::mapq_boost;

    const struct option long_opts[] = {
        {"map",             required_argument, nullptr, 'm'},
        {"in",              required_argument, nullptr, 'i'},
        {"out",             required_argument, nullptr, 'o'},
        {"genome",          required_argument, nullptr, 'g'},

        {"batch",           required_argument, nullptr, 1001},
        {"mapq_low",        required_argument, nullptr, 1002},
        {"mapq_cap",        required_argument, nullptr, 1003},
        {"name_check",      no_argument,       nullptr, 1004},

        {"cm_max_hops",     required_argument, nullptr, 2001},
        {"cm_max_fanout",   required_argument, nullptr, 2002},
        {"cm_min_len",      required_argument, nullptr, 2003},
        {"cm_min_frac",     required_argument, nullptr, 2004},
        {"cm_max_hits",     required_argument, nullptr, 2005},

        {"sub_ovlp_frac",   required_argument, nullptr, 3001},
        {"K_mapq",          required_argument, nullptr, 3002},
        {"K_rs",            required_argument, nullptr, 3003},
        {"K_as",            required_argument, nullptr, 3004},
        {"K_ml",            required_argument, nullptr, 3005},
        {"K_nm",            required_argument, nullptr, 3006},
        {"W_mapq",          required_argument, nullptr, 3007},
        {"W_rs",            required_argument, nullptr, 3008},
        {"W_as",            required_argument, nullptr, 3009},
        {"W_ml",            required_argument, nullptr, 3010},
        {"W_nm",            required_argument, nullptr, 3011},
        {"close_as_eps",    required_argument, nullptr, 3012},

        {"enzyme",          required_argument, nullptr, 'e'},
        {"rc",              no_argument,       nullptr, 4001},

        {"threads",         required_argument, nullptr, 't'},
        {"io_threads",      required_argument, nullptr, 5001},
        {"debug",           no_argument,       nullptr, 'd'},
        {"help",            no_argument,       nullptr, 'h'},
        {"advanced",        no_argument,       nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "m:i:o:g:e:t:dhH";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'm': {
                cfg.homq.mapFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.homq.mapFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'i': {
                cfg.homq.in_bam = optarg;
                break;
            }
            case 'o': {
                cfg.homq.out_bam = optarg;
                break;
            }
            case 'g': {
                cfg.ressit.genome_file = optarg;
                break;
            }

            case 1001: {
                cfg.homq.batch_size = std::max(1, std::atoi(optarg));
                break;
            }
            case 1002: {
                cfg.homq.mapq_low = (uint8_t)std::clamp(std::atoi(optarg), 0, 255);
                break;
            }
            case 1003: {
                cfg.homq.mapq_cap = (uint8_t)std::clamp(std::atoi(optarg), 0, 255); 
                break;
            }
            case 1004: {
                cfg.homq.name_check = true;
                break;
            }

            case 2001: {
                cfg.coordmap.max_hops = std::max(1, std::atoi(optarg));
                break;
            }
            case 2002: {
                cfg.coordmap.max_fanout = (uint32_t)std::max(1, std::atoi(optarg));
                break;
            }
            case 2003: {
                cfg.coordmap.min_len = (uint32_t)std::max(1, std::atoi(optarg));
                break;
            }
            case 2004: {
                cfg.coordmap.min_frac = std::max(0.0, std::atof(optarg));
                break;
            }
            case 2005: {
                cfg.coordmap.max_total_hits = (uint32_t)std::max(1, std::atoi(optarg));
                break;
            }

            case 3001: {
                cfg.homq.sub_ovlp_frac = std::max(0.0, std::atof(optarg));
                break;
            }
            case 3002: {
                cfg.homq.K_mapq = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3003: {
                cfg.homq.K_rs = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3004: {
                cfg.homq.K_as = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3005: {
                cfg.homq.K_ml = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3006: {
                cfg.homq.K_nm = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3007: {
                cfg.homq.W_mapq = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3008: {
                cfg.homq.W_rs = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3009: {
                cfg.homq.W_as = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3010: {
                cfg.homq.W_ml = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3011: {
                cfg.homq.W_nm = std::max(1.0, std::atof(optarg));
                break;
            }
            case 3012: {
                cfg.homq.close_as_eps = std::max(1e-6, std::atof(optarg));
                break;
            }

            case  'e': {
                cfg.ressit.enzymes.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.ressit.enzymes.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 4001: {
                cfg.ressit.scan_revcomp = true;
                break;
            }

            case 't': {
                cfg.global.threads = std::max(1, std::atoi(optarg));
                break;
            }
            case 5001: {
                cfg.global.IOthreads = std::max(1, std::atoi(optarg));
                break;
            }
            case 'd': {
                cfg.global.debug = true;
                break;
            }
            case 'h': {
                help_mapq_boost(argv);
                std::exit(0);
            }
            case 'H': {
                help_mapq_boost(argv, true);
                std::exit(0);
            }
            default: {
                help_mapq_boost(argv);
                std::exit(1);
            }
        }
    }

    if (cfg.homq.in_bam.empty())  cfg.homq.in_bam  = "-";
    if (cfg.homq.out_bam.empty()) cfg.homq.out_bam = "-";

    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_align(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... -r FILE ... [options]\n\n"
        << "Align sequencing data to a GFA\n\n"
        << "Input/Output:\n"
        << "  -g, --gfa         FILE ...    input GFA file(s)\n"
        << "  -r, --read        FILE ...    input sequencing data file(s)\n"
        << "  -p, --paf                     output PAF (default SAM)\n"
        << "  -o, --output      FILE        output alignments to file [stdout]\n\n"
        << "Index options:\n"
        << "  -k, --kmer        INT         k-mer length [" << GlobalOpts().kmerLen << "]\n"
        << "  -w, --window      INT         minimizer window [" << GlobalOpts().minimizerW << "]\n\n"
        << "Alignment options:\n"
        << "  -S, --secondary   FLOAT       min ratio secondary/primary [" << MapOpts().sec_pri_ratio << "]\n"
        << "  -N                INT         number of secondary to keep [" << MapOpts().sec_pri_num << "]\n"
        << "  -P, --preset      STR         hifi/ont/illumina/other [" << MapOpts().preset << "]\n\n"
        << "General Options:\n"
        << "  -t, --threads     INT         number of threads [" << GlobalOpts().threads << "]\n"
        << "  -d, --debug                   debug mode (forces threads=1)\n"
        << "  -h, --help                    show basic options\n\n";
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
            case 'g': {
                cfg.map.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.map.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'r': {
                cfg.map.reads.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.map.reads.emplace_back(argv[optind]); ++optind;
                }
                break;
            }
            case 'p': {
                cfg.map.outPAF = true;
                break;
            }
            case 'o': {
                cfg.map.outFile = optarg;
                break;
            }
            case 'S': {
                cfg.map.sec_pri_ratio = std::stod(optarg);
                break;
            }
            case 'N': {
                cfg.map.sec_pri_num = std::stoi(optarg);
                break;
            }
            case 'P': {
                cfg.map.preset = optarg;
                break;
            }
            case 'k': {
                cfg.global.kmerLen = std::stoi(optarg);
                break;
            }
            case 'w': {
                cfg.global.minimizerW = std::stoi(optarg);
                break;
            }
            case 't': {
                cfg.global.threads = std::max(1, std::stoi(optarg));
                break;
            }
            case 'd': {
                cfg.global.debug = true;
                break;
            }
            case 'h': {
                help_align(argv);
                std::exit(0);
            }
            default: {
                help_align(argv);
                std::exit(1);
            }
        }
    }
    
    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_res_cut(char** argv) {
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... -e STR ... [options]\n\n"
        << "Build restriction cut-site index from a genome\n\n"
        << "Input/Output:\n"
        << "  -g, --genome      FILE    input genome FASTA\n"
        << "  -o, --output      FILE    output BED file (optional)\n\n"
        << "Scan options:\n"
        << "  -e, --enzyme      STR     enzyme recognition site with cut mark '^', can take multiple arguments after one -e\n"
        << "                             * e.g. G^AATTC  A^AGCTT  G^ANTC\n"
        << "      --rc                  scan reverse-complement motifs\n\n"
        << "General options:\n"
        << "  -t, --threads     INT     number of threads [" << GlobalOpts().threads << "]\n"
        << "  -h, --help                show basic options\n\n";
}

AppConfig main_res_cut(int argc, char** argv) {
    if (argc < 3) { help_res_cut(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::ressit;

    const struct option long_opts[] = {
        {"genome",    required_argument, nullptr, 'g'},
        {"output",    required_argument, nullptr, 'o'},
        {"enzyme",    required_argument, nullptr, 'e'},
        {"rc",        no_argument,       nullptr, 1001},
        {"threads",   required_argument, nullptr, 't'},
        {"help",      no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:o:e:t:h";

    auto is_flag_ = [](const char* s) -> bool {
        return s && s[0] == '-';
    };

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.ressit.genome_file = optarg;
                break;
            }
            case 'o': {
                cfg.ressit.out_bed = optarg;
                break;
            }
            case 'e': {
                cfg.ressit.enzymes.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.ressit.enzymes.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 1001: {
                cfg.ressit.scan_revcomp = true;
                break;
            }
            case 't': {
                cfg.global.threads = std::max(1, std::stoi(optarg));
                break;
            }
            case 'h': {
                help_res_cut(argv);
                std::exit(0);
            }
            default: {
                help_res_cut(argv);
                std::exit(1);
            }
        }
    }

    validate_and_print(argc, argv, cfg);
    return cfg;
}

void finalize_opt_cfg(AppConfig& cfg) {
    init_opts(
        cfg.global.kmerLen, cfg.global.minimizerW,
        cfg.map.sec_pri_ratio, cfg.map.sec_pri_num,
        cfg.map.outPAF, cfg.global.threads,
        cfg.map.chainOpts, cfg.map.anchorOpts, cfg.map.extendOpts, cfg.map.alignOpts
    );

    if (cfg.mode == ToolMode::collapse || cfg.mode == ToolMode::deoverlap) {
        opt::Preset::map_asm_5(cfg.map.chainOpts, cfg.map.anchorOpts, cfg.map.extendOpts, cfg.map.alignOpts);
    }

    if (cfg.map.zdrop_set) {
        cfg.map.extendOpts.dyn_zdrop = cfg.map.zdrop;
    }
}
