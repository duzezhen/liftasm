#include "../include/OptionParser.hpp"
#include "../include/OptionParser_helper.hpp"
#include "../include/CollapseConfig.hpp"
#include "../include/ProgramMetadata.hpp"
#include "../include/logger.hpp"
#include "../include/get_time.hpp"

#include <getopt.h>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

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
        cfg.mode == ToolMode::split      ? "split" :
        cfg.mode == ToolMode::augment    ? "augment" :
        cfg.mode == ToolMode::clean      ? "clean" :
        cfg.mode == ToolMode::gapfill    ? "gapfill" :
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
                cfg.bubble.cx_branch_degree >= 3,
                "--cx_branch_deg must be >= 3"
            );
            ensure(
                cfg.bubble.cx_hub_degree >= cfg.bubble.cx_branch_degree,
                "--cx_hub_deg must be >= --cx_branch_deg"
            );
            ensure(
                cfg.bubble.cx_min_nodes >= 1,
                "--cx_min_nodes must be >= 1"
            );
            ensure(
                cfg.bubble.cx_min_branches >= 1,
                "--cx_min_branches must be >= 1"
            );
            ensure(
                cfg.bubble.path_sim >= 0.0 && cfg.bubble.path_sim <= 1.0,
                "--path_sim must be in [0,1]"
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
                cfg.bubble.homo_extend_bp > 0, 
                "--hextend must be > 0"
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
                cfg.collapse.min_trans_len >= 0,
                "--min_trans_len must be >= 0"
            );
            ensure(
                cfg.map.zdrop >= 0,
                "--zdrop must be >= 0"
            );
            ensure(
                cfg.collapse.min_match_ratio >= 0.0 && cfg.collapse.min_match_ratio <= 1.0,
                "--min_match must be in [0,1]"
            );
            ensure(
                cfg.collapse.min_ali_ratio >= 0.0 && cfg.collapse.min_ali_ratio <= 1.0,
                "--min_ali_ratio must be in [0,1]"
            );
            ensure(
                cfg.collapse.min_mapq <= 60,
                "--min_mapq must be <= 60"
            );
            ensure(
                cfg.collapse.trim_min_len > 0,
                "--trim_len must be > 0"
            );
            ensure(
                cfg.collapse.trim_max_overlap >= 0.0 && cfg.collapse.trim_max_overlap <= 1.0,
                "--trim_ovlp must be in [0,1]"
            );
            ensure(
                cfg.collapse.mm2_preset == "asm5" || cfg.collapse.mm2_preset == "asm10" || cfg.collapse.mm2_preset == "asm20" || cfg.collapse.mm2_preset == "sr" || cfg.collapse.mm2_preset == "lr:hq",
                "-x/--preset must be one of: asm5, asm10, asm20, sr, lr:hq"
            );
            break;

        case ToolMode::collapse:
            ensure(
                !cfg.collapse.configFile.empty(),
                "-i/--input is required"
            );
            ensure(
                cfg.collapse.input.hap1Files.size() == cfg.collapse.input.hap2Files.size(),
                "invalid CTG input: haplotype file counts differ"
            );
            ensure(
                cfg.collapse.input.vcfFiles.empty() || !cfg.collapse.input.hap1Files.empty(),
                "VCF input is only supported in CTG mode"
            );
            ensure(
                cfg.collapse.input.vcfFiles.empty() || cfg.collapse.input.vcfFiles.size() == cfg.collapse.input.hap1Files.size(),
                "invalid CTG input: VCF and sample counts differ"
            );
            ensure(
                cfg.collapse.ctg_anchor_coverage >= 0.0 && cfg.collapse.ctg_anchor_coverage <= 1.0,
                "--ctg_anchor_cov must be in [0,1]"
            );
            ensure(
                cfg.collapse.ctg_min_coverage >= 0.0 && cfg.collapse.ctg_min_coverage <= 1.0,
                "--ctg_cov must be in [0,1]"
            );
            ensure(
                cfg.collapse.ctg_end_fraction >= 0.0 && cfg.collapse.ctg_end_fraction <= 1.0,
                "--ctg_end must be in [0,1]"
            );
            ensure(
                !cfg.collapse.ctg_anchor_only || !cfg.collapse.input.hap1Files.empty(),
                "--anchor_only is only supported in CTG mode"
            );
            ensure(
                cfg.collapse.input.hap1Files.empty() || cfg.collapse.input.sampleNames.size() == cfg.collapse.input.hap1Files.size(),
                "invalid CTG input: sample and GFA counts differ"
            );
            ensure(
                !cfg.collapse.input.hap1Files.empty() || cfg.collapse.input.sampleNames.size() == cfg.collapse.input.gfaFiles.size(),
                "invalid UTG input: sample and GFA counts differ"
            );
            ensure(
                !cfg.collapse.prefix.empty(), 
                "-p/--prefix is required"
            );
            ensure(
                cfg.collapse.iterations >= 1,
                "--iterations must be >= 1"
            );
            ensure(
                !cfg.collapse.min_jaccards.empty(),
                "--min_jaccard must contain at least one value"
            );
            for (double x : cfg.collapse.min_jaccards) {
                ensure(
                    x >= 0.0 && x <= 1.0,
                    "--min_jaccard values must be in [0,1]"
                );
            }
            ensure(
                !cfg.collapse.min_eqs.empty(), 
                "--min_eq must contain at least one value"
            );
            for (int x : cfg.collapse.min_eqs) {
                ensure(
                    x > 0, 
                    "--min_eq values must be positive"
                );
            }
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
                cfg.collapse.repeat_mask_max_mismatch >= 0,
                "--repeat_mask MAX_MISMATCH must be >= 0"
            );
            ensure(
                cfg.collapse.repeat_norm_len >= 0,
                "--repeat_norm_len must be >= 0"
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
                cfg.collapse.min_trans_len >= 0,
                "--min_trans_len must be >= 0"
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
                cfg.bubble.cx_branch_degree >= 3,
                "--cx_branch_deg must be >= 3"
            );
            ensure(
                cfg.bubble.cx_hub_degree >= cfg.bubble.cx_branch_degree,
                "--cx_hub_deg must be >= --cx_branch_deg"
            );
            ensure(
                cfg.bubble.cx_min_nodes >= 1,
                "--cx_min_nodes must be >= 1"
            );
            ensure(
                cfg.bubble.cx_min_branches >= 1,
                "--cx_min_branches must be >= 1"
            );
            ensure(
                cfg.bubble.path_sim >= 0.0 && cfg.bubble.path_sim <= 1.0,
                "--path_sim must be in [0,1]"
            );
            ensure(
                cfg.bubble.stall_round_limit >= 0, 
                "--stall_rounds must be >= 0"
            );
            ensure(
                !cfg.collapse.same_sims.empty(),
                "--same_sim must contain at least one value"
            );
            for (double x : cfg.collapse.same_sims) {
                ensure(
                    x >= 0.0 && x <= 1.0,
                    "--same_sim values must be in [0,1]"
                );
            }
            ensure(
                !cfg.collapse.same_min_lens.empty(),
                "--same_len must contain at least one value"
            );
            ensure(
                !cfg.collapse.diff_min_srcs.empty(),
                "--diff_src_len must contain at least one value"
            );
            ensure(
                cfg.bubble.homo_num >= 1,
                "--hnum must be >= 1"
            );
            ensure(
                !cfg.collapse.diff_sims.empty(),
                "--diff_sim must contain at least one value"
            );
            for (double x : cfg.collapse.diff_sims) {
                ensure(
                    x >= 0.0 && x <= 1.0,
                    "--diff_sim values must be in [0,1]"
                );
            }
            ensure(
                !cfg.collapse.diff_min_lens.empty(),
                "--diff_len must contain at least one value"
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
                cfg.bubble.homo_extend_bp > 0,
                "--hextend must be > 0"
            );
            ensure(
                cfg.bubble.homo_bloom_bits > 0,
                "--hbits must be > 0"
            );
            ensure(
                cfg.bubble.homo_bloom_hash > 0,
                "--hhash must be > 0"
            );
            ensure(
                cfg.map.zdrop >= 0,
                "--zdrop must be >= 0"
            );
            ensure(
                !cfg.collapse.min_match_ratios.empty(),
                "--min_match must contain at least one value"
            );
            for (double x : cfg.collapse.min_match_ratios) {
                ensure(
                    x >= 0.0 && x <= 1.0,
                    "--min_match values must be in [0,1]"
                );
            }
            ensure(
                !cfg.collapse.min_ali_ratios.empty(),
                "--min_ali_ratio must contain at least one value"
            );
            for (double x : cfg.collapse.min_ali_ratios) {
                ensure(
                    x >= 0.0 && x <= 1.0,
                    "--min_ali_ratio values must be in [0,1]"
                );
            }
            ensure(
                cfg.collapse.min_mapq <= 60,
                "--min_mapq must be <= 60"
            );
            ensure(
                cfg.collapse.trim_min_len > 0,
                "--trim_len must be > 0"
            );
            ensure(
                cfg.collapse.trim_max_overlap >= 0.0 && cfg.collapse.trim_max_overlap <= 1.0,
                "--trim_ovlp must be in [0,1]"
            );
            ensure(
                cfg.collapse.mm2_preset == "asm5" || cfg.collapse.mm2_preset == "asm10" || cfg.collapse.mm2_preset == "asm20" || cfg.collapse.mm2_preset == "sr" || cfg.collapse.mm2_preset == "lr:hq",
                "-x/--preset must be one of: asm5, asm10, asm20, sr, lr:hq"
            );
            break;

        case ToolMode::file2map:
            ensure(
                !cfg.file2map.input_files.empty(),
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

        case ToolMode::split:
            ensure(
                !cfg.split.gfaFiles.empty(),
                "-g/--gfa is required"
            );
            ensure(
                cfg.split.gfaNames.empty() || cfg.split.gfaNames.size() == cfg.split.gfaFiles.size(),
                "-n/--name must have the same number of values as -g/--gfa"
            );
            break;

        case ToolMode::augment:
            ensure(
                !cfg.augment.gfa_files.empty(), 
                "-g/--gfa is required"
            );
            ensure(
                !cfg.augment.vcf_files.empty(), 
                "-v/--vcf is required"
            );
            ensure(
                cfg.augment.gfa_names.empty() || cfg.augment.gfa_names.size() == cfg.augment.gfa_files.size(),
                   "-n/--name must have the same number of values as -g/--gfa"
                );
            ensure(
                !cfg.augment.prefix.empty(), 
                "-p/--prefix is required"
            );
            ensure(
                !cfg.augment.tag.empty(), 
                "--tag must not be empty"
            );
            ensure(
                cfg.collapse.min_match_ratio >= 0.0 && cfg.collapse.min_match_ratio <= 1.0,
                "--min_match must be in [0,1]"
            );
            ensure(
                cfg.collapse.min_ali_ratio >= 0.0 && cfg.collapse.min_ali_ratio <= 1.0,
                "--min_ali_ratio must be in [0,1]"
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
                cfg.collapse.mm2_preset == "asm5" || cfg.collapse.mm2_preset == "asm10" || cfg.collapse.mm2_preset == "asm20" || cfg.collapse.mm2_preset == "sr" || cfg.collapse.mm2_preset == "lr:hq",
                "-x/--preset must be one of: asm5, asm10, asm20, sr, lr:hq"
            );
            break;

        case ToolMode::clean:
            ensure(
                !cfg.clean.utg_file.empty(), 
                "-g/--gfa is required"
            );
            ensure(
                !cfg.clean.ctg_files.empty(), 
                "-c/--ctg is required"
            );
            ensure(
                !cfg.clean.prefix.empty(), 
                "-p/--prefix is required"
            );
            ensure(
                cfg.clean.min_reads >= 1, 
                "--min_reads must be >= 1"
            );
            ensure(
                cfg.clean.min_purity >= 0.0 && cfg.clean.min_purity <= 1.0,
                "--min_purity must be in [0,1]"
            );
            ensure(
                cfg.clean.min_comp_overlap >= 0.0 && cfg.clean.min_comp_overlap <= 1.0,
                "--min_comp_overlap must be in [0,1]"
            );
            break;

        case ToolMode::gapfill:
            ensure(
                !cfg.gapfill.gfa_file.empty(), 
                "-g/--gfa is required"
            );
            ensure(
                !cfg.gapfill.prefix.empty(), 
                "-p/--prefix is required"
            );
            ensure(
                cfg.gapfill.max_depth >= 1, 
                "--depth must be >= 1"
            );
            ensure(
                cfg.gapfill.max_paths >= 1, 
                "--paths must be >= 1"
            );
            ensure(
                cfg.gapfill.DFS_guard >= 1, 
                "--DFS_guard must be >= 1"
            );
            ensure(
                cfg.gapfill.path_sim >= 0.0 && cfg.gapfill.path_sim <= 1.0,
                "--path_sim must be in [0,1]"
            );
            ensure(
                cfg.gapfill.phase_path_len > 0,
                "--phase_len must be > 0"
            );
            ensure(
                cfg.gapfill.phase_win > 0,
                "--phase_win must be > 0"
            );
            ensure(
                cfg.gapfill.min_overlap > 0, 
                "--min_overlap must be > 0"
            );
            ensure(
                cfg.gapfill.min_contig > 0, 
                "--min_contig must be > 0"
            );
            ensure(
                cfg.gapfill.max_gap > 0, 
                "--max_gap must be > 0"
            );
            ensure(
                cfg.gapfill.min_similarity >= 0.0 && cfg.gapfill.min_similarity <= 1.0,
                "--min_similarity must be in [0,1]"
            );
            ensure(
                cfg.gapfill.max_overlap >= 0.0 && cfg.gapfill.max_overlap <= 1.0,
                "--max_overlap must be in [0,1]"
            );
            ensure(
                cfg.gapfill.ms_sim >= 0.0 && cfg.gapfill.ms_sim <= 1.0,
                "--ms_sim must be in [0,1]"
            );
            ensure(
                cfg.gapfill.dedup_similarity >= 0.0 && cfg.gapfill.dedup_similarity <= 1.0,
                "--dedup_sim must be in [0,1]"
            );
            ensure(
                cfg.gapfill.dedup_component > 0,
                "--dedup_component must be > 0"
            );
            ensure(
                cfg.gapfill.min_match >= 0.0 && cfg.gapfill.min_match <= 1.0,
                "--min_match must be in [0,1]"
            );
            ensure(
                cfg.gapfill.min_ali_ratio >= 0.0 && cfg.gapfill.min_ali_ratio <= 1.0,
                "--min_ali_ratio must be in [0,1]"
            );
            ensure(
                cfg.gapfill.mm2_preset == "asm5" || cfg.gapfill.mm2_preset == "asm10" || cfg.gapfill.mm2_preset == "asm20" || cfg.gapfill.mm2_preset == "sr" || cfg.gapfill.mm2_preset == "lr:hq",
                "--preset must be asm5, asm10, asm20, sr, or lr:hq"
            );
            break;
    }
}

void help(char** argv, bool advanced) {
    HelpPrinter hp(std::cerr, 18, 0);
    std::cerr 
        << "Usage: " << argv[0] << " <subcommand> [options]\n\n"
        << program::description << "\n\n"
        << "Version: " << program::version << "\n"
        << "Date:    " << program::build_date << "\n";

    if (advanced) {
        std::cerr << "Note:    " << program::version_note << "\n";
    }

    hp.blank();
    
    hp.section("Subcommands");
    hp.line("stat", "", "collect statistics about a GFA file");
    hp.line("seq", "", "extract sequences along paths from a GFA");
    hp.line("bubble", "", "detect bubbles in a GFA");
    hp.line("augment", "", "augment a GFA with VCF alleles as variant bubbles");
    hp.line("clean", "", "remove weak UTG links using contig read components (beta)");
    hp.line("collapse", "", "de-overlap and collapse homologous sequences into single nodes");
    hp.line("gapfill", "", "fill sample contig gaps with source-aware graph walks");
    hp.line("file2map", "", "convert GFA/PAF to map format for liftover");
    hp.line("liftover", "", "liftover coordinates using a coordinate map");
    hp.line("mapq_boost", "", "raise MAPQ for alignments in homologous regions using a coordinate map");

    if (advanced) {
        hp.line("gfa2fa", "", "export all segments to FASTA with two-sided extension for short nodes");
        hp.line("depth", "", "compute depth information from sequencing data for a GFA");
        hp.line("deoverlap", "", "convert an overlap-based GFA to a non-overlap GFA");
        hp.line("align", "", "align sequencing data to a GFA (beta)");
        hp.line("ressit", "", "find restriction sites in a genome");
        hp.line("split", "", "split graph into multiple GFAs based on components");
    }

    hp.blank();
    
    hp.section("General Options");
    hp.line("-h, --help", "", "show basic options");
    hp.line("-H, --advanced", "", "show advanced options");
    
    hp.blank();
}

void help_stat(char** argv) {
    HelpPrinter hp(std::cerr, 13, 11);
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ...\n\n"
        << "Collect statistics about GFA file\n";
    
    hp.blank();

    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    
    hp.blank();
    
    hp.section("General Options");
    hp.line("-h, --help", "", "show basic options");
    
    hp.blank();
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
    HelpPrinter hp(std::cerr, 16, 13);
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... -p FILE [options]\n\n"
        << "Extract sequences for compact marker paths or open walks (one record per line)\n";
    
    hp.blank();

    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    hp.line("-p, --path", "FILE", "text file of paths/open walks (one per line)");
    hp.note(" * path format: '>S1<S2>S3' or 'S1+,S2-,S3+'");
    hp.note(" * open-walk format: 'S1<TAB>+/-<TAB>LEN'");
    hp.line("-o, --output", "FILE", "output file, fasta format [stdout]");
    
    hp.blank();
    
    hp.section("General Options");
    hp.line("-h, --help", "", "show basic options");
    
    hp.blank();
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
    HelpPrinter hp;
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Export all segments of a GFA to FASTA, extending both ends by Y bp\n"
        << "for segments whose effective length is shorter than X bp.\n";
    
    hp.blank();
    
    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    hp.line("-o, --output", "FILE", "output FASTA file [stdout]");
    
    hp.blank();
    
    hp.section("Extension options");
    hp.line("--min_len", "INT", "only extend nodes with len < X [" + std::to_string(Gfa2FaOpts().min_len_xbp) + "]");
    hp.line("--extend", "INT", "extend Y bp on both sides [" + std::to_string(Gfa2FaOpts().extend_ybp) + "]");
    hp.line("--wrap", "INT", "FASTA wrap width; 0 = no wrap [" + std::to_string(Gfa2FaOpts().wrap_width) + "]");
    hp.line("--keep_unknown", "", "include nodes with '*' or zero length");
    hp.note(" * by default such nodes are skipped");
    
    hp.blank();
    
    hp.section("General Options");
    hp.line("-h, --help", "", "show basic options");
    
    hp.blank();
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
    HelpPrinter hp;
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Compute node depth information from GFA (A-lines) / reads / GAF for a GFA graph\n\n"
        << "Three ways to compute depth:\n"
        << "  1. k-mer mode (provide -r):\n"
        << "       " << argv[0] << " " << argv[1] << " -g graph.gfa -r reads.fq -o graph.depth\n"
        << "  2. alignment mode (provide --gaf):\n"
        << "       " << argv[0] << " " << argv[1] << " -g graph.gfa --gaf aln.gaf -o graph.depth\n"
        << "  3. GFA-only mode (depth is calculated from 'A' lines in the GFA file):\n"
        << "       " << argv[0] << " " << argv[1] << " -g graph.gfa -o graph.depth\n";
    
    hp.blank();
    
    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    hp.line("-r, --read", "FILE", "input sequencing data file(s) for k-mer mode");
    hp.note(" * depth calculated from k-mer coverage");
    hp.line("--gaf", "FILE", "input GAF alignment file for alignment mode");
    hp.note(" * depth calculated directly from alignment");
    hp.line("-o, --output", "FILE", "output depth information to file [stdout]");
    hp.line("--base_depth", "", "output per-base depth");
    
    hp.blank();
    
    hp.section("Index options");
    hp.line("-k, --kmer", "INT", "k-mer length [" + std::to_string(GlobalOpts().kmerLen) + "]");
    
    hp.blank();
    
    hp.section("GFA loading options");
    hp.line("--min_mapq", "INT", "minimum mapping quality (MAPQ) to keep [" + std::to_string(DepthOpts().min_mapq) + "]");
    hp.line("--min_frac", "FLOAT", "minimum aligned-length / read-length ratio to keep [" + format_double_(DepthOpts().min_frac) + "]");
    
    hp.blank();
    
    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-d, --debug", "", "debug mode (forces threads=1)");
    hp.line("-h, --help", "", "show basic options");
    
    hp.blank();
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
    HelpPrinter hp;
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Detect bubbles in GFA graph and write them to GFA/VCF file\n";
    
    hp.blank();

    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    hp.line("-n, --name", "STR ...", "sample names; otherwise reuse liftasm tags or use filename prefixes");
    hp.line("-p, --prefix", "STR", "output prefix [" + BubbleOpts().out_prefix + "]");
    hp.note(" * <prefix>.bubbles.{gfa,noseq.gfa,vcf}");
    hp.line("--write_vcf", "", "also output bubbles in VCF format");
    hp.line("--ref", "FILE", "reference FASTA used to extract REF sequence for VCF");
    hp.line("--paf", "FILE", "PAF file for liftover of contig coordinates to reference");
    hp.note(" * must contain cg tag; e.g. minimap2 -cx asm5 ref.fa ctg.fa > aln.paf");
    hp.line("--ali_min_mapq", "INT", "minimum mapping quality for PAF liftover [" + format_size_arg_(BubbleOpts().ali_min_mapq) + "]");
    hp.line("--ali_min_len", "INT", "minimum alignment length for PAF liftover [" + format_size_arg_(BubbleOpts().ali_min_len) + "]");
    
    hp.blank();
    
    hp.section("Detection Options");
    hp.line("--depth", "INT", "maximum DFS depth for path exploration inside a bubble [" + format_size_arg_(BubbleOpts().max_depth) + "]");
    hp.line("--paths", "INT", "maximum number of DFS paths to explore per bubble [" + format_size_arg_(BubbleOpts().max_paths) + "]");
    if (advanced) {
        hp.line("--DFS_guard", "INT", "max DFS states [" + format_size_arg_(BubbleOpts().DFS_guard) + "]");
    }
    hp.line("--path_sim", "FLOAT", "minimum minimizer Jaccard for path clustering [" + format_double_(BubbleOpts().path_sim) + "]");
    hp.line("--stall_rounds", "INT", "stop DFS path search after this many rounds without a new unique path [" + format_size_arg_(BubbleOpts().stall_round_limit) + "]");
    hp.line("--min_len", "INT", "minimum total sequence length of a bubble for output (0 = no filter) [" + format_size_arg_(BubbleOpts().min_len) + "]");
    hp.line("--min_num", "INT", "minimum number of nodes inside a bubble required for output (0 = no filter) [" + format_size_arg_(BubbleOpts().min_num) + "]");
    hp.line("--keep_nested", "", "keep bubbles contained within larger bubbles");
    hp.line("--no_cx", "", "disable complex graph region detection");
    if (advanced) {
        hp.line("--cx_branch_deg", "INT", "minimum node degree counted as a branch [" + format_size_arg_(BubbleOpts().cx_branch_degree) + "]");
        hp.line("--cx_hub_deg", "INT", "minimum node degree counted as a hub [" + format_size_arg_(BubbleOpts().cx_hub_degree) + "]");
        hp.line("--cx_min_nodes", "INT", "minimum nodes for a branch-rich block to be complex [" + format_size_arg_(BubbleOpts().cx_min_nodes) + "]");
        hp.line("--cx_min_branches", "INT", "minimum branches in a complex block [" + format_size_arg_(BubbleOpts().cx_min_branches) + "]");
        
        hp.blank();
        
        hp.section("Homologous Path Options");
        hp.line("--same_sim", "FLOAT", "minimum similarity for homologous paths from the same source [" + format_double_(BubbleOpts().same_sim) + "]");
        hp.line("--same_len", "INT", "minimum total length of each same-source homologous path [" + format_size_arg_(BubbleOpts().same_min_len) + "]");
        hp.line("--diff_src_len", "INT", "minimum source length for searching homologous paths between different sources [" + format_size_arg_(BubbleOpts().diff_min_src) + "]");
        hp.line("--diff_sim", "FLOAT", "minimum similarity for homologous paths between different sources [" + format_double_(BubbleOpts().diff_sim) + "]");
        hp.line("--diff_len", "INT", "minimum total length of each different-source homologous path [" + format_size_arg_(BubbleOpts().diff_min_len) + "]");
        hp.line("--hnum", "INT", "max paths per source branch [" + format_size_arg_(BubbleOpts().homo_num) + "]");
        hp.line("--hk", "INT", "k-mer size for minimizer used in homologous path search [" + format_size_arg_(BubbleOpts().homo_k) + "]");
        hp.line("--hw", "INT", "window size for minimizer used in homologous path search [" + format_size_arg_(BubbleOpts().homo_w) + "]");
        hp.line("--hextend", "INT", "bp extended per homologous-path extension round [" + format_size_arg_(BubbleOpts().homo_extend_bp) + "]");
        hp.line("--hbits", "INT", "Bloom filter size in bits for homologous-path sketches [" + format_size_arg_(BubbleOpts().homo_bloom_bits) + "]");
        hp.line("--hhash", "INT", "number of Bloom filter hash functions for homologous-path sketches [" + format_size_arg_(BubbleOpts().homo_bloom_hash) + "]");
    }
    
    hp.blank();
    
    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-d, --debug", "", "debug mode (forces threads=1)");
    hp.line("-h, --help", "", "show basic options");
    hp.line("-H, --advanced", "", "show advanced options");
    
    hp.blank();
}

AppConfig main_bubble(int argc, char** argv) {
    if (argc < 3) { help_bubble(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::bubble;

    const struct option long_opts[] = {
        {"gfa",             required_argument, nullptr, 'g'},
        {"name",            required_argument, nullptr, 'n'},
        {"prefix",          required_argument, nullptr, 'p'},
        {"write_vcf",       no_argument,       nullptr, 0001},
        {"ref",             required_argument, nullptr, 0002},
        {"paf",             required_argument, nullptr, 0003},
        {"ali_min_mapq",    required_argument, nullptr, 0004},
        {"ali_min_len",     required_argument, nullptr, 0005},

        {"depth",           required_argument, nullptr, 1001},
        {"paths",           required_argument, nullptr, 1002},
        {"DFS_guard",       required_argument, nullptr, 1003},
        {"path_sim",        required_argument, nullptr, 1004},
        {"stall_rounds",    required_argument, nullptr, 1005},
        {"min_len",         required_argument, nullptr, 1006},
        {"min_num",         required_argument, nullptr, 1007},
        {"keep_nested",     no_argument,       nullptr, 1008},
        {"no_cx",           no_argument,       nullptr, 1009},
        {"cx_branch_deg",   required_argument, nullptr, 1010},
        {"cx_hub_deg",      required_argument, nullptr, 1011},
        {"cx_min_nodes",    required_argument, nullptr, 1012},
        {"cx_min_branches", required_argument, nullptr, 1013},

        {"same_sim",        required_argument, nullptr, 2001},
        {"same_len",        required_argument, nullptr, 2002},
        {"diff_src_len",    required_argument, nullptr, 2003},
        {"diff_sim",        required_argument, nullptr, 2004},
        {"diff_len",        required_argument, nullptr, 2005},
        {"hnum",            required_argument, nullptr, 2006},
        {"hk",              required_argument, nullptr, 2007},
        {"hw",              required_argument, nullptr, 2008},
        {"hextend",         required_argument, nullptr, 2009},
        {"hbits",           required_argument, nullptr, 2010},
        {"hhash",           required_argument, nullptr, 2011},

        {"threads",         required_argument, nullptr, 't'},
        {"debug",           no_argument,       nullptr, 'd'},
        {"help",            no_argument,       nullptr, 'h'},
        {"advanced",        no_argument,       nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "g:n:p:t:dhH";

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
            case 'n': {
                cfg.bubble.gfaNames.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.bubble.gfaNames.emplace_back(argv[optind]);
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
                cfg.bubble.ali_min_mapq = parse_size_arg_u32_(optarg, argc, argv, optind, "--ali_min_mapq");
                break;
            }
            case 0005: {
                cfg.bubble.ali_min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--ali_min_len");
                break;
            }

            case 1001: {
                cfg.bubble.max_depth = parse_size_arg_u32_(optarg, argc, argv, optind, "--depth");
                break;
            }
            case 1002: {
                cfg.bubble.max_paths = parse_size_arg_u16_(optarg, argc, argv, optind, "--paths");
                break;
            }
            case 1003: {
                cfg.bubble.DFS_guard = parse_size_arg_u64_(optarg, argc, argv, optind, "--DFS_guard");
                break;
            }
            case 1004: {
                cfg.bubble.path_sim = std::stod(optarg);
                break;
            }
            case 1005: {
                cfg.bubble.stall_round_limit = parse_size_arg_u32_(optarg, argc, argv, optind, "--stall_rounds");
                break;
            }
            case 1006: {
                cfg.bubble.min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--min_len");
                break;
            }
            case 1007: {
                cfg.bubble.min_num = parse_size_arg_u32_(optarg, argc, argv, optind, "--min_num");
                break;
            }
            case 1008: {
                cfg.bubble.keep_nested = true;
                break;
            }
            case 1009: {
                cfg.bubble.check_complex = false;
                break;
            }
            case 1010: {
                cfg.bubble.cx_branch_degree = parse_size_arg_u16_(optarg, argc, argv, optind, "--cx_branch_deg");
                break;
            }
            case 1011: {
                cfg.bubble.cx_hub_degree = parse_size_arg_u16_(optarg, argc, argv, optind, "--cx_hub_deg");
                break;
            }
            case 1012: {
                cfg.bubble.cx_min_nodes = parse_size_arg_u32_(optarg, argc, argv, optind, "--cx_min_nodes");
                break;
            }
            case 1013: {
                cfg.bubble.cx_min_branches = parse_size_arg_u32_(optarg, argc, argv, optind, "--cx_min_branches");
                break;
            }

            case 2001: {
                cfg.bubble.same_sim = std::stod(optarg);
                break;
            }
            case 2002: {
                cfg.bubble.same_min_len = parse_size_arg_u64_(optarg, argc, argv, optind, "--same_len");
                break;
            }
            case 2003: {
                cfg.bubble.diff_min_src = parse_size_arg_u32_(optarg, argc, argv, optind, "--diff_src_len");
                break;
            }
            case 2004: {
                cfg.bubble.diff_sim = std::stod(optarg);
                break;
            }
            case 2005: {
                cfg.bubble.diff_min_len = parse_size_arg_u64_(optarg, argc, argv, optind, "--diff_len");
                break;
            }
            case 2006: {
                cfg.bubble.homo_num = parse_size_arg_u32_(optarg, argc, argv, optind, "--hnum");
                break;
            }
            case 2007: {
                cfg.bubble.homo_k = parse_size_arg_u32_(optarg, argc, argv, optind, "--hk");
                break;
            }
            case 2008: {
                cfg.bubble.homo_w = parse_size_arg_u32_(optarg, argc, argv, optind, "--hw");
                break;
            }
            case 2009: {
                cfg.bubble.homo_extend_bp = parse_size_arg_u32_(optarg, argc, argv, optind, "--hextend");
                break;
            }
            case 2010: {
                cfg.bubble.homo_bloom_bits = parse_size_arg_u64_(optarg, argc, argv, optind, "--hbits");
                break;
            }
            case 2011: {
                cfg.bubble.homo_bloom_hash = parse_size_arg_u32_(optarg, argc, argv, optind, "--hhash");
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
    HelpPrinter hp(std::cerr, 20, 13);
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Convert an overlap-based GFA into a non-overlap GFA\n";

    hp.blank();
    
    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    hp.line("-n, --name", "STR ...", "sample names; otherwise reuse liftasm tags or use filename prefixes");
    hp.line("-p, --prefix", "STR", "output prefix [" + CollapseOpts().prefix + "]");
    hp.note(" * <prefix>.deoverlap.{gfa,noseq.gfa,map}");

    hp.blank();
    
    hp.section("Index options");
    hp.line("-k, --kmer", "INT", "k-mer length [" + std::to_string(GlobalOpts().kmerLen) + "]");
    hp.line("-w, --window", "INT", "minimizer window [" + std::to_string(GlobalOpts().minimizerW) + "]");

    hp.blank();
    
    hp.section("Align options");
    hp.line("-x, --preset", "STR", "mapping preset: asm5/asm10/asm20/sr/lr:hq [" + CollapseOpts().mm2_preset + "]");
    if (advanced) {
        hp.line("--use_wfa", "", "whether to use WFA for alignment, if not use mm2 (beta)");
    }
    hp.line("-z, --zdrop", "INT", "Z-drop score [" + format_size_arg_(MapOpts().extendOpts.dyn_zdrop) + "]");
    hp.line("--min_match", "FLOAT", "minimum match fraction in the alignment CIGAR [" + format_double_(CollapseOpts().min_match_ratio) + "]");
    hp.line("--min_ali_ratio", "FLOAT", "minimum aligned-length fraction relative to the sequence length [" + format_double_(CollapseOpts().min_ali_ratio) + "]");
    hp.line("--min_mapq", "INT", "minimum mapping quality (MAPQ) to keep (no larger than 60) [" + std::to_string(int(CollapseOpts().min_mapq)) + "]");
    if (advanced) {
        hp.line("--trim_len", "INT", "minimum length of both alignments to trim a small overlap [" + format_size_arg_(CollapseOpts().trim_min_len) + "]");
        hp.line("--trim_ovlp", "FLOAT", "maximum overlap fraction of the shorter alignment to trim [" + format_double_(CollapseOpts().trim_max_overlap) + "]");
    }
    
    hp.blank();
    
    hp.section("Collapse options");
    hp.line("--min_eq", "INT", "minimum match length to add cut points at segment [" + std::to_string(CollapseOpts().min_eq) + "]");
    hp.line("--max_iters", "INT", "maximum iterations for cut point propagation [" + format_size_arg_(CollapseOpts().max_iters) + "]");
    if (advanced) {
        hp.line("--abnormal_cut", "INT,INT", "prune abnormal cut points after propagation as MAX_LEN,MIN_COUNT [" + std::to_string(CollapseOpts().max_abnormal_cut_len) + "," + std::to_string(CollapseOpts().min_abnormal_cut_count) + "]");
        hp.note(" * e.g. 0-3-6-9-100 -> 0-9-100 when 4,2 due to consecutive short segments");
        hp.line("--min_trans_len", "INT", "minimum interval length to allow transitive replacement expansion [" + format_size_arg_(CollapseOpts().min_trans_len) + "]");
    }

    hp.blank();
    
    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-d, --debug", "", "debug mode");
    hp.line("-h, --help", "", "show basic options");
    hp.line("-H, --advanced", "", "show advanced options");
    
    hp.blank();
}

AppConfig main_deoverlap(int argc, char** argv) {
    if (argc < 3) { help_deoverlap(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::deoverlap;

    const struct option long_opts[] = {
        {"gfa",          required_argument, nullptr, 'g'},
        {"name",         required_argument, nullptr, 'n'},
        {"prefix",       required_argument, nullptr, 'p'},

        {"kmer",         required_argument, nullptr, 'k'},
        {"window",       required_argument, nullptr, 'w'},

        {"preset",       required_argument, nullptr, 'x'},
        {"use_wfa",      no_argument,       nullptr, 2001},
        {"zdrop",        required_argument, nullptr, 'z'},
        {"min_match",    required_argument, nullptr, 2002},
        {"min_ali_ratio",required_argument, nullptr, 2003},
        {"min_mapq",     required_argument, nullptr, 2004},
        {"trim_len",     required_argument, nullptr, 2005},
        {"trim_ovlp",    required_argument, nullptr, 2006},

        {"min_eq",       required_argument, nullptr, 3001},
        {"max_iters",    required_argument, nullptr, 3002},
        {"abnormal_cut", required_argument, nullptr, 3003},
        {"min_trans_len",required_argument, nullptr, 3004},

        {"threads",      required_argument, nullptr, 't'},
        {"debug",        no_argument,       nullptr, 'd'},
        {"help",         no_argument,       nullptr, 'h'},
        {"advanced",     no_argument,       nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "g:n:p:k:w:x:z:t:dhH";

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
            case 'n': {
                cfg.collapse.gfaNames.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.collapse.gfaNames.emplace_back(argv[optind]);
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

            case 'x': {
                cfg.collapse.mm2_preset = optarg;
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
            case 2003: {
                cfg.collapse.min_ali_ratio = std::max(0.0, std::stod(optarg));
                break;
            }
            case 2004: {
                const unsigned long v = std::stoul(optarg);
                if (v > 60) {
                    error_stream() << "--min_mapq must be <= 60\n";
                    std::exit(1);
                }
                cfg.collapse.min_mapq = static_cast<uint8_t>(v);
                break;
            }
            case 2005: {
                cfg.collapse.trim_min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--trim_len");
                break;
            }
            case 2006: {
                cfg.collapse.trim_max_overlap = std::stod(optarg);
                break;
            }

            case 3001: {
                cfg.collapse.min_eq = std::stoi(optarg);
                break;
            }
            case 3002: {
                cfg.collapse.max_iters = std::stoi(optarg);
                break;
            }
            case 3003: {
                parse_abnormal_seg_(optarg, cfg.collapse);
                break;
            }
            case 3004: {
                cfg.collapse.min_trans_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--min_trans_len");
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
    HelpPrinter hp(std::cerr, 20, 20);

    std::cerr
        << "Usage: " << argv[0] << " collapse -i FILE [options]\n\n"
        << "De-overlap (if needed) and collapse homologous sequences in GFA into single nodes\n";

    hp.blank();

    hp.section("Input/Output");
    hp.line("-i, --input", "FILE", "header-driven, whitespace-delimited sample/GFA table");
    hp.note(" * UTG config: sample  gfa");
    hp.note(" * CTG config: sample  hap1_gfa  hap2_gfa  [vcf]");
    hp.line("-p, --prefix", "STR", "output prefix [" + CollapseOpts().prefix + "]");
    hp.note(" * UTG: <prefix>.iterN.collapse.{gfa,noseq.gfa,map}");
    hp.note(" * CTG: <prefix>.<sample>.collapse.{gfa,noseq.gfa,map} and <prefix>.iterN.collapse.{gfa,noseq.gfa,map}");
    hp.note(" * contig mode always writes filtered and raw component alignments as *.paf and *.all.paf");

    hp.blank();
    
    hp.section("Index options");
    hp.line("-k, --kmer", "INT", "k-mer length [" + std::to_string(GlobalOpts().kmerLen) + "]");
    hp.line("-w, --window", "INT", "minimizer window [" + std::to_string(GlobalOpts().minimizerW) + "]");

    hp.blank();
    
    hp.section("Align options");
    hp.line("-x, --preset", "STR", "mapping preset: asm5/asm10/asm20/sr/lr:hq [" + CollapseOpts().mm2_preset + "]");
    if (advanced) {
        hp.line("--use_wfa", "", "whether to use WFA for alignment, if not use mm2 (beta)");
    }
    hp.line("-z, --zdrop", "INT", "Z-drop score [" + format_size_arg_(MapOpts().extendOpts.dyn_zdrop) + "]");
    hp.line("--min_match", "FLOAT[,FLOAT...]", "minimum CIGAR match fractions by iteration [" + join_doubles_(CollapseOpts().min_match_ratios) + "]; contig mode uses the first value");
    hp.line("--min_ali_ratio", "FLOAT[,FLOAT...]", "minimum aligned-length fractions by iteration [" + join_doubles_(CollapseOpts().min_ali_ratios) + "]; contig mode uses the first value");
    hp.line("--min_mapq", "INT", "minimum mapping quality (MAPQ) to keep (no larger than 60) [" + std::to_string(int(CollapseOpts().min_mapq)) + "]");
    if (advanced) {
        hp.line("--all_pair_len", "INT", "align all path pairs when every path is shorter than this [" + format_size_arg_(CollapseOpts().all_pair_len) + "]");
        hp.line("--trim_len", "INT", "minimum length of both alignments to trim a small overlap [" + format_size_arg_(CollapseOpts().trim_min_len) + "]");
        hp.line("--trim_ovlp", "FLOAT", "maximum overlap fraction of the shorter alignment to trim [" + format_double_(CollapseOpts().trim_max_overlap) + "]");
    }
    hp.line("--anchor_only", "", "align only shared-read anchor regions in contig collapse");
    
    hp.blank();
    
    hp.section("Collapse options");
    hp.line("--iterations", "INT", "number of collapse iterations [" + std::to_string(CollapseOpts().iterations) + "]");
    hp.note(" * lists shorter than --iterations reuse their last value");
    hp.line("--min_jaccard", "FLOAT[,FLOAT...]", "Jaccard thresholds by collapse iteration [" + join_doubles_(CollapseOpts().min_jaccards) + "]; contig mode uses the first value");
    hp.line("--min_eq", "INT[,INT...]", "minimum match lengths for cut points [" + join_ints_(CollapseOpts().min_eqs) + "]; contig mode uses the first value");
    hp.line("--max_iters", "INT", "maximum iterations for cut point propagation [" + format_size_arg_(CollapseOpts().max_iters) + "]");
    hp.line("--ctg_min_len", "INT", "ignore shorter input contigs in contig mode; 0 disables [" + format_size_arg_(CollapseOpts().ctg_min_len) + "]");
    if (advanced) {
        hp.line("--repeat_mask", "INT,INT,INT", "tandem-repeat detection as MIN_LEN,MAX_PERIOD,MAX_MISMATCH [" + std::to_string(CollapseOpts().repeat_mask_min_len) + "," + std::to_string(CollapseOpts().repeat_mask_max_period) + "," + std::to_string(CollapseOpts().repeat_mask_max_mismatch) + "]");
        hp.note(" * UTG masks repeat matches during alignment; CTG normalizes repeat-only P-line bubbles after the final collapse");
        hp.note(" * e.g. ATATATGATCG contains a 9 bp repeat with period 2 (AT) and 1 mismatch (G), detectable by 8,12,1");
        hp.line("--repeat_norm_len", "INT", "max single-node bubble path length for UTG repeat normalization; 0 disables [" + format_size_arg_(CollapseOpts().repeat_norm_len) + "]");
        hp.line("--abnormal_cut", "INT,INT", "prune abnormal cut points after propagation as MAX_LEN,MIN_COUNT [" + std::to_string(CollapseOpts().max_abnormal_cut_len) + "," + std::to_string(CollapseOpts().min_abnormal_cut_count) + "]");
        hp.note(" * e.g. 0-3-6-9-100 -> 0-9-100 when 4,2 due to consecutive short segments");
        hp.line("--min_trans_len", "INT", "minimum interval length to allow transitive replacement expansion [" + format_size_arg_(CollapseOpts().min_trans_len) + "]");
    }

    hp.blank();
    
    hp.section("Bubble Detection Options");
    hp.line("--depth", "INT", "maximum DFS depth for path exploration inside a bubble [" + format_size_arg_(BubbleOpts().max_depth) + "]");
    hp.line("--paths", "INT", "maximum DFS paths per bubble; 0 = max(20, twice the input sample count) [0]");
    if (advanced) {
        hp.line("--DFS_guard", "INT", "max DFS states [" + format_size_arg_(BubbleOpts().DFS_guard) + "]");
    }
    hp.line("--path_sim", "FLOAT", "minimum minimizer Jaccard for path clustering [" + format_double_(BubbleOpts().path_sim) + "]");
    hp.line("--stall_rounds", "INT", "stop DFS path search after this many rounds without a new unique path [" + format_size_arg_(BubbleOpts().stall_round_limit) + "]");
    hp.line("--no_cx", "", "disable complex graph region detection");
    if (advanced) {
        hp.line("--cx_branch_deg", "INT", "minimum node degree counted as a branch [" + format_size_arg_(BubbleOpts().cx_branch_degree) + "]");
        hp.line("--cx_hub_deg", "INT", "minimum node degree counted as a hub [" + format_size_arg_(BubbleOpts().cx_hub_degree) + "]");
        hp.line("--cx_min_nodes", "INT", "minimum nodes for a branch-rich block to be complex [" + format_size_arg_(BubbleOpts().cx_min_nodes) + "]");
        hp.line("--cx_min_branches", "INT", "minimum branches in a complex block [" + format_size_arg_(BubbleOpts().cx_min_branches) + "]");
    }

    if (advanced) {
        hp.blank();
        hp.section("Homologous Path Detection Options");
        hp.line("--same_sim", "FLOAT[,FLOAT...]", "same-source homologous path similarity thresholds by collapse iteration [" + join_doubles_(CollapseOpts().same_sims) + "]");
        hp.line("--same_len", "INT[,INT...]", "minimum total lengths of each same-source homologous path by collapse iteration [" + join_sizes_(CollapseOpts().same_min_lens) + "]");
        hp.line("--diff_src_len", "INT[,INT...]", "minimum source lengths for different-source search by collapse iteration [" + join_sizes_(CollapseOpts().diff_min_srcs) + "]");
        hp.line("--diff_sim", "FLOAT[,FLOAT...]", "different-source homologous path similarity thresholds by collapse iteration [" + join_doubles_(CollapseOpts().diff_sims) + "]");
        hp.line("--diff_len", "INT[,INT...]", "minimum total lengths of each different-source homologous path by collapse iteration [" + join_sizes_(CollapseOpts().diff_min_lens) + "]");
        hp.line("--hnum", "INT", "max paths per source branch [" + format_size_arg_(BubbleOpts().homo_num) + "]");
        hp.line("--hk", "INT", "k-mer size for minimizer used in homologous path search [" + format_size_arg_(BubbleOpts().homo_k) + "]");
        hp.line("--hw", "INT", "window size for minimizer used in homologous path search [" + format_size_arg_(BubbleOpts().homo_w) + "]");
        hp.line("--hextend", "INT", "bp extended per homologous-path extension round [" + format_size_arg_(BubbleOpts().homo_extend_bp) + "]");
        hp.line("--hbits", "INT", "Bloom filter size in bits for homologous-path sketches [" + format_size_arg_(BubbleOpts().homo_bloom_bits) + "]");
        hp.line("--hhash", "INT", "number of Bloom filter hash functions for homologous-path sketches [" + format_size_arg_(BubbleOpts().homo_bloom_hash) + "]");
        hp.line("--ctg_anchor_cov", "FLOAT", "shared-read coverage required to consider two contigs homologous [" + format_double_(CollapseOpts().ctg_anchor_coverage) + "]");
        hp.line("--ctg_short", "INT", "short contig length for one-partner filtering; 0 disables [" + format_size_arg_(CollapseOpts().ctg_short_len) + "]");
        hp.line("--ctg_cov", "FLOAT", "minimum aligned span coverage for component merging [" + format_double_(CollapseOpts().ctg_min_coverage) + "]");
        hp.line("--ctg_end", "FLOAT", "maximum terminal overhang fraction for end-to-end support; 0 disables [" + format_double_(CollapseOpts().ctg_end_fraction) + "]");
    }

    hp.blank();
    
    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads (~5GB RAM per thread) [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-d, --debug", "", "debug mode (forces threads=1)");
    hp.line("-h, --help", "", "show basic options");
    hp.line("-H, --advanced", "", "show advanced options");
    
    hp.blank();
}

AppConfig main_collapse(int argc, char** argv) {
    if (argc < 3) { help_collapse(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::collapse;
    cfg.bubble.max_paths = 0;

    const struct option long_opts[] = {
        {"input",           required_argument, nullptr, 'i'},
        {"prefix",          required_argument, nullptr, 'p'},

        {"kmer",            required_argument, nullptr, 'k'},
        {"window",          required_argument, nullptr, 'w'},

        {"preset",          required_argument, nullptr, 'x'},
        {"use_wfa",         no_argument,       nullptr, 2001},
        {"zdrop",           required_argument, nullptr, 'z'},
        {"min_match",       required_argument, nullptr, 2002},
        {"min_ali_ratio",   required_argument, nullptr, 2003},
        {"min_mapq",        required_argument, nullptr, 2004},
        {"all_pair_len",    required_argument, nullptr, 2005},
        {"trim_len",        required_argument, nullptr, 2006},
        {"trim_ovlp",       required_argument, nullptr, 2007},
        {"anchor_only",     no_argument,       nullptr, 2008},

        {"iterations",      required_argument, nullptr, 3001},
        {"min_jaccard",     required_argument, nullptr, 3002},
        {"min_eq",          required_argument, nullptr, 3003},
        {"max_iters",       required_argument, nullptr, 3004},
        {"ctg_min_len",     required_argument, nullptr, 3005},
        {"repeat_mask",     required_argument, nullptr, 3006},
        {"repeat_norm_len", required_argument, nullptr, 3007},
        {"abnormal_cut",    required_argument, nullptr, 3008},
        {"min_trans_len",   required_argument, nullptr, 3009},

        {"depth",           required_argument, nullptr, 4001},
        {"paths",           required_argument, nullptr, 4002},
        {"DFS_guard",       required_argument, nullptr, 4003},
        {"path_sim",        required_argument, nullptr, 4004},
        {"stall_rounds",    required_argument, nullptr, 4005},
        {"no_cx",           no_argument,       nullptr, 4006},
        {"cx_branch_deg",   required_argument, nullptr, 4007},
        {"cx_hub_deg",      required_argument, nullptr, 4008},
        {"cx_min_nodes",    required_argument, nullptr, 4009},
        {"cx_min_branches", required_argument, nullptr, 4010},

        {"same_sim",        required_argument, nullptr, 5001},
        {"same_len",        required_argument, nullptr, 5002},
        {"diff_src_len",    required_argument, nullptr, 5003},
        {"diff_sim",        required_argument, nullptr, 5004},
        {"diff_len",        required_argument, nullptr, 5005},
        {"hnum",            required_argument, nullptr, 5006},
        {"hk",              required_argument, nullptr, 5007},
        {"hw",              required_argument, nullptr, 5008},
        {"hextend",         required_argument, nullptr, 5009},
        {"hbits",           required_argument, nullptr, 5010},
        {"hhash",           required_argument, nullptr, 5011},
        {"ctg_anchor_cov",  required_argument, nullptr, 5012},
        {"ctg_short",       required_argument, nullptr, 5013},
        {"ctg_cov",         required_argument, nullptr, 5014},
        {"ctg_end",         required_argument, nullptr, 5015},

        {"threads",         required_argument, nullptr, 't'},
        {"debug",           no_argument,       nullptr, 'd'},
        {"help",            no_argument,       nullptr, 'h'},
        {"advanced",        no_argument,       nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "i:p:k:w:x:z:t:dhH";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'i': {
                if (!cfg.collapse.configFile.empty()) {
                    error_stream() << "-i/--input may be specified only once\n";
                    std::exit(1);
                }
                cfg.collapse.configFile = optarg;
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

            case 'x': {
                cfg.collapse.mm2_preset = optarg;
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
                cfg.collapse.min_match_ratios = parse_double_list_(optarg);
                if (!cfg.collapse.min_match_ratios.empty()) {
                    cfg.collapse.min_match_ratio = cfg.collapse.min_match_ratios.front();
                }
                break;
            }
            case 2003: {
                cfg.collapse.min_ali_ratios = parse_double_list_(optarg);
                if (!cfg.collapse.min_ali_ratios.empty()) cfg.collapse.min_ali_ratio = cfg.collapse.min_ali_ratios.front();
                break;
            }
            case 2004: {
                const unsigned long v = std::stoul(optarg);
                if (v > 60) {
                    error_stream() << "--min_mapq must be <= 60\n";
                    std::exit(1);
                }
                cfg.collapse.min_mapq = static_cast<uint8_t>(v);
                break;
            }
            case 2005: {
                cfg.collapse.all_pair_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--all_pair_len");
                break;
            }
            case 2006: {
                cfg.collapse.trim_min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--trim_len");
                break;
            }
            case 2007: {
                cfg.collapse.trim_max_overlap = std::stod(optarg);
                break;
            }
            case 2008: {
                cfg.collapse.ctg_anchor_only = true;
                break;
            }

            case 3001: {
                cfg.collapse.iterations = parse_size_arg_u32_(optarg, argc, argv, optind, "--iterations");
                break;
            }
            case 3002: {
                cfg.collapse.min_jaccards = parse_double_list_(optarg);
                break;
            }
            case 3003: {
                cfg.collapse.min_eqs = parse_int_list_(optarg);

                while (optind < argc && !is_flag_(argv[optind])) {
                    std::vector<int> more = parse_int_list_(argv[optind]);
                    cfg.collapse.min_eqs.insert(cfg.collapse.min_eqs.end(), more.begin(), more.end());
                    ++optind;
                }

                if (!cfg.collapse.min_eqs.empty()) cfg.collapse.min_eq = cfg.collapse.min_eqs.front();

                break;
            }
            case 3004: {
                cfg.collapse.max_iters = parse_size_arg_u16_(optarg, argc, argv, optind, "--max_iters");
                break;
            }
            case 3005: {
                cfg.collapse.ctg_min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--ctg_min_len");
                break;
            }
            case 3006: {
                parse_repeat_mask_(optarg, cfg.collapse);
                break;
            }
            case 3007: {
                cfg.collapse.repeat_norm_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--repeat_norm_len");
                break;
            }
            case 3008: {
                parse_abnormal_seg_(optarg, cfg.collapse);
                break;
            }
            case 3009: {
                cfg.collapse.min_trans_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--min_trans_len");
                break;
            }

            case 4001: {
                cfg.bubble.max_depth = parse_size_arg_u32_(optarg, argc, argv, optind, "--depth");
                break;
            }
            case 4002: {
                cfg.bubble.max_paths = parse_size_arg_u16_(optarg, argc, argv, optind, "--paths");
                break;
            }
            case 4003: {
                cfg.bubble.DFS_guard = parse_size_arg_u64_(optarg, argc, argv, optind, "--DFS_guard");
                break;
            }
            case 4004: {
                cfg.bubble.path_sim = std::stod(optarg);
                break;
            }
            case 4005: {
                cfg.bubble.stall_round_limit = parse_size_arg_u32_(optarg, argc, argv, optind, "--stall_rounds");
                break;
            }
            case 4006: {
                cfg.bubble.check_complex = false;
                break;
            }
            case 4007: {
                cfg.bubble.cx_branch_degree = parse_size_arg_u16_(optarg, argc, argv, optind, "--cx_branch_deg");
                break;
            }
            case 4008: {
                cfg.bubble.cx_hub_degree = parse_size_arg_u16_(optarg, argc, argv, optind, "--cx_hub_deg");
                break;
            }
            case 4009: {
                cfg.bubble.cx_min_nodes = parse_size_arg_u32_(optarg, argc, argv, optind, "--cx_min_nodes");
                break;
            }
            case 4010: {
                cfg.bubble.cx_min_branches = parse_size_arg_u32_(optarg, argc, argv, optind, "--cx_min_branches");
                break;
            }

            case 5001: {
                cfg.collapse.same_sims = parse_double_list_(optarg);
                break;
            }
            case 5002: {
                cfg.collapse.same_min_lens = parse_size_list_u32_(optarg, "--same_len");
                break;
            }
            case 5003: {
                cfg.collapse.diff_min_srcs = parse_size_list_u32_(optarg, "--diff_src_len");
                break;
            }
            case 5004: {
                cfg.collapse.diff_sims = parse_double_list_(optarg);
                break;
            }
            case 5005: {
                cfg.collapse.diff_min_lens = parse_size_list_u32_(optarg, "--diff_len");
                break;
            }
            case 5006: {
                cfg.bubble.homo_num = parse_size_arg_u32_(optarg, argc, argv, optind, "--homo_num");
                break;
            }
            case 5007: {
                cfg.bubble.homo_k = parse_size_arg_u32_(optarg, argc, argv, optind, "--hk");
                break;
            }
            case 5008: {
                cfg.bubble.homo_w = parse_size_arg_u32_(optarg, argc, argv, optind, "--hw");
                break;
            }
            case 5009: {
                cfg.bubble.homo_extend_bp = parse_size_arg_u32_(optarg, argc, argv, optind, "--hextend");
                break;
            }
            case 5010: {
                cfg.bubble.homo_bloom_bits = parse_size_arg_u64_(optarg, argc, argv, optind, "--hbits");
                break;
            }
            case 5011: {
                cfg.bubble.homo_bloom_hash = parse_size_arg_u32_(optarg, argc, argv, optind, "--hhash");
                break;
            }
            case 5012: {
                cfg.collapse.ctg_anchor_coverage = std::stod(optarg);
                break;
            }
            case 5013: {
                cfg.collapse.ctg_short_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--ctg_short");
                break;
            }
            case 5014: {
                cfg.collapse.ctg_min_coverage = std::stod(optarg);
                break;
            }
            case 5015: {
                cfg.collapse.ctg_end_fraction = std::stod(optarg);
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

    if (!cfg.collapse.configFile.empty()) {
        cfg.collapse.input = parse_collapse_config(cfg.collapse.configFile);
    }

    if (cfg.bubble.max_paths == 0) {
        const size_t input_count = cfg.collapse.input.hap1Files.empty()
            ? cfg.collapse.input.gfaFiles.size()
            : cfg.collapse.input.hap1Files.size();
        const size_t automatic = std::max<size_t>(20, input_count * 2);
        cfg.bubble.max_paths = static_cast<uint16_t>(
            std::min<size_t>(automatic, std::numeric_limits<uint16_t>::max())
        );
    }
    
    validate_and_print(argc, argv, cfg);

    finalize_opt_cfg(cfg);

    return cfg;
}

void help_file2map(char** argv) {
    HelpPrinter hp;
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
        << "       to avoid confusion with map coordinate formats such as 'chr:beg-end+'.\n";
    
    hp.blank();

    hp.section("Input/Output");
    hp.line("-i, --input", "FILE ...", "input PAF/GFA file(s)");
    hp.line("-o, --output", "FILE", "output file name [stdout]");
    
    hp.blank();
    
    hp.section("PAF options");
    hp.line("--all", "", "include secondary/supplementary, (default: only primary tp:A:P)");
    hp.line("--min_len", "INT", "minimum alignment length to keep [" + format_size_arg_(File2mapOpts().min_len) + "]");
    hp.line("--min_mapq", "INT", "minimum mapping quality (MAPQ) to keep [" + std::to_string(File2mapOpts().min_mapq) + "]");
    
    hp.blank();
    
    hp.section("General Options");
    hp.line("-h, --help", "", "show basic options");
    
    hp.blank();
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
                cfg.file2map.input_files.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.file2map.input_files.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'o': {
                cfg.file2map.output_file = optarg;
                break;
            }
            case 1001: {
                cfg.file2map.paf_primary_only = false;
                break;
            }
            case 1002: {
                cfg.file2map.min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--min_len");
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

    if (cfg.file2map.output_file.empty()) cfg.file2map.output_file = "-";

    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_liftover(char** argv, bool advanced) {
    HelpPrinter hp;
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -m FILE ... -b FILE [options]\n\n"
        << "LiftOver BED coordinates using a .map file, which produced by deoverlap/collapse/file2map\n";

    hp.blank();

    hp.section("Input/Output");
    hp.line("-m, --map", "FILE ...", "coordinate mapping file(s) (generated by deoverlap/collapse/file2map)");
    hp.line("--paf", "FILE", "input PAF file to assist liftover (optional)");
    hp.line("-b, --bed", "FILE", "input BED file, required if not checking");
    hp.line("-o, --output", "FILE", "output file name [stdout]");
    hp.line("-r, --ref", "FILE", "FASTA/FASTQ containing sequences for src+dst contigs, required if checking");
    
    hp.blank();
    
    hp.section("Map Options");
    hp.line("--regex", "STR", "keep hits whose label matches regex");
    hp.note(" * \"h\\d+tg\" to match h1tg, h2tg, ...");
    hp.line("--min_frac", "FLOAT", "minimum fraction of the BED interval that must be lifted over [" + format_double_(LiftoverOpts().min_frac) + "]");
    hp.note(" * 0.9 means at least 90% of the BED region is mapped");
    if (advanced) {
        hp.line("--flank_win", "INT", "extension window size around BED ends to identify insertions [" + format_size_arg_(LiftoverOpts().flank_win) + "]");
        hp.line("--max_flank", "INT", "max extension around BED ends to identify insertions [" + format_size_arg_(LiftoverOpts().max_flank) + "]");
        hp.line("--max_gap", "INT", "max allowed |ref_gap - query_gap| when calling INS [" + format_size_arg_(LiftoverOpts().max_gap) + "]");
        hp.line("--max_hit", "INT", "max extended liftover results to keep [" + format_size_arg_(LiftoverOpts().max_hit) + "]");
    }
    
    hp.blank();
    
    hp.section("PAF Options");
    hp.line("--min_len", "INT", "minimum alignment length for liftover [" + format_size_arg_(LiftoverOpts().min_len) + "]");
    hp.line("--min_mapq", "INT", "minimum mapping quality for liftover [" + format_size_arg_(LiftoverOpts().min_mapq) + "]");
    
    hp.blank();
    
    hp.section("Check Options");
    hp.line("--check", "", "run sequence consistency check, requires -r");
    if (advanced) {
        hp.line("-w, --win", "INT", "window size [" + format_size_arg_(LiftoverOpts().win) + "]");
        hp.line("-s, --step", "INT", "step size [" + format_size_arg_(LiftoverOpts().step) + "]");
        hp.line("--max_examples", "INT", "max reported examples [" + format_size_arg_(LiftoverOpts().max_examples) + "]");
    }  
    hp.blank();
    
    hp.section("CoordMap Options");
    hp.line("--cm_max_hops", "INT", "max hop depth when chaining names [" + std::to_string(CoordMapOpts().max_hops) + "]");
    hp.note(" * e.g. 1 allows: source -> A; 2 allows: source -> A -> B");

    if (advanced) {
        hp.line("--cm_max_fanout", "INT", "max fanout per search [" + format_size_arg_(CoordMapOpts().max_fanout) + "]");
        hp.line("--cm_min_len", "INT", "min length (bp) to keep a hit for hops [" + format_size_arg_(CoordMapOpts().min_len) + "]");
        hp.line("--cm_min_frac", "FLOAT", "min fraction of current interval length to keep a hit [" + format_double_(CoordMapOpts().min_frac) + "]");
        hp.note(" * threshold = max(min_len, cur_len * min_frac)");
        hp.line("--cm_max_hits", "INT", "global BFS node cap [" + format_size_arg_(CoordMapOpts().max_total_hits) + "]");
    }

    hp.blank();
    
    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-h, --help", "", "show basic options");
    hp.line("-H, --advanced", "", "show advanced options");
    
    hp.blank();
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
                cfg.liftover.flank_win = parse_size_arg_u32_(optarg, argc, argv, optind, "--flank_win");
                break;
            }
            case 1004: {
                cfg.liftover.max_flank = parse_size_arg_u32_(optarg, argc, argv, optind, "--max_flank");
                break;
            }
            case 1005: {
                cfg.liftover.max_gap = parse_size_arg_u32_(optarg, argc, argv, optind, "--max_gap");
                break;
            }
            case 1006: {
                cfg.liftover.max_hit = parse_size_arg_u16_(optarg, argc, argv, optind, "--max_hit");
                break;
            }

            case 2001: {
                cfg.liftover.min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--min_len");
                break;
            }
            case 2002: {
                cfg.liftover.min_mapq = parse_size_arg_u32_(optarg, argc, argv, optind, "--min_mapq");
                break;
            }

            case 3001: {
                cfg.liftover.do_check = true;
                break;
            }
            case 'w': {
                cfg.liftover.win = parse_size_arg_u32_(optarg, argc, argv, optind, "-w");
                break;
            }
            case 's': {
                cfg.liftover.step = parse_size_arg_u32_(optarg, argc, argv, optind, "-s");
                break;
            }
            case 3002: {
                cfg.liftover.max_examples = parse_size_arg_u32_(optarg, argc, argv, optind, "--max_examples");
                break;
            }

            case 4001: {
                cfg.coordmap.max_hops = std::max(1, std::atoi(optarg));
                break;
            }
            case 4002: {
                cfg.coordmap.max_fanout = parse_size_arg_u32_(optarg, argc, argv, optind, "--cm_max_fanout");
                break;
            }
            case 4003: {
                cfg.coordmap.min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--cm_min_len");
                break;
            }
            case 4004: {
                cfg.coordmap.min_frac = std::max(0.0, std::atof(optarg));
                break;
            }
            case 4005: {
                cfg.coordmap.max_total_hits = parse_size_arg_u32_(optarg, argc, argv, optind, "--cm_max_hits");
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
    HelpPrinter hp;
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
        << kTail;
    
    hp.blank();

    hp.section("Input/Output");
    hp.line("-m, --map", "FILE ...", "coordinate mapping file(s) (generated by deoverlap/collapse/file2map)");
    hp.line("-i, --in", "FILE", "input SAM/BAM [stdin]");
    hp.line("-o, --out", "FILE", "output SAM/BAM [stdout]");
    hp.line("-g, --genome", "FILE", "input genome FASTA for restriction-site search (optional)");
    hp.note(" * filter alignments for enzyme-cut data, e.g. Hi-C, Pore-C, CiFi");
    
    hp.blank();
    
    hp.section("Boost Options");
    if (advanced) {
        hp.line("--batch", "INT", "batch size [" + format_size_arg_(MapqBoostOpts().batch_size) + "]");
    }
    hp.line("--mapq_low", "INT", "only raise if MAPQ<=Q (no larger than 255) [" + std::to_string(int(MapqBoostOpts().mapq_low)) + "]");
    hp.line("--mapq_cap", "INT", "upper bound of boosted MAPQ (no larger than 255) [" + std::to_string(int(MapqBoostOpts().mapq_cap)) + "]");
    hp.line("--name_check", "", "check whether QNAMEs are name-sorted (off by default)");
    
    hp.blank();
    
    hp.section("CoordMap Options");
    hp.line("--cm_max_hops", "INT", "BFS depth limit [" + std::to_string(CoordMapOpts().max_hops) + "]");
    hp.note(" * e.g. 1 allows: source -> A; 2 allows: source -> A -> B");
    if (advanced) {
        hp.line("--cm_max_fanout", "INT", "max fanout per search [" + format_size_arg_(CoordMapOpts().max_fanout) + "]");
        hp.line("--cm_min_len", "INT", "min length (bp) to keep a hit for hops [" + format_size_arg_(CoordMapOpts().min_len) + "]");
        hp.line("--cm_min_frac", "FLOAT", "min fraction of current interval length to keep a hit [" + format_double_(CoordMapOpts().min_frac) + "]");
        hp.note(" * threshold = max(min_len, cur_len * min_frac)");
        hp.line("--cm_max_hits", "INT", "global BFS node cap [" + format_size_arg_(CoordMapOpts().max_total_hits) + "]");
    }
    
    hp.blank();
    
    hp.section("Subgrouping Options");
    hp.line("--sub_ovlp_frac", "FLOAT", "minimum overlap fraction on read to cluster alignments into the same subgroup [" + format_double_(MapqBoostOpts().sub_ovlp_frac) + "]");
    hp.note(" * alignments overlapping the same read region are evaluated together");
    
    if (advanced) {
        hp.blank();
        hp.section("Scoring Options, used to rank alignments within a subgroup and select the best candidate for boosting");
        hp.line("--K_mapq", "FLOAT", "MAPQ score scaling factor [" + format_double_(MapqBoostOpts().K_mapq) + "]");
        hp.line("--K_rs", "FLOAT", "restriction-site distance scaling factor [" + format_double_(MapqBoostOpts().K_rs) + "]");
        hp.line("--K_as", "FLOAT", "AS score scaling factor [" + format_double_(MapqBoostOpts().K_as) + "]");
        hp.line("--K_ml", "FLOAT", "match-length score scaling factor [" + format_double_(MapqBoostOpts().K_ml) + "]");
        hp.line("--K_nm", "FLOAT", "NM score scaling factor [" + format_double_(MapqBoostOpts().K_nm) + "]");
        hp.line("--W_mapq", "FLOAT", "MAPQ term weight [" + format_double_(MapqBoostOpts().W_mapq) + "]");
        hp.line("--W_rs", "FLOAT", "restriction-site term weight [" + format_double_(MapqBoostOpts().W_rs) + "]");
        hp.line("--W_as", "FLOAT", "AS term weight [" + format_double_(MapqBoostOpts().W_as) + "]");
        hp.line("--W_ml", "FLOAT", "match-length term weight [" + format_double_(MapqBoostOpts().W_ml) + "]");
        hp.line("--W_nm", "FLOAT", "NM penalty weight [" + format_double_(MapqBoostOpts().W_nm) + "]");
        hp.line("--close_as_eps", "FLOAT", "max relative AS difference to consider two alignments equally scored for MAPQ recalculation [" + format_double_(MapqBoostOpts().close_as_eps) + "]");
    }
    
    hp.blank();
    
    hp.section("Restriction-site Options (optional)");
    hp.line("-e, --enzyme", "STR", "enzyme recognition site with cut mark '^', can take multiple arguments after one -e");
    hp.note(" * e.g. G^AATTC  A^AGCTT  G^ANTC");
    hp.line("--rc", "", "scan reverse-complement motifs");
    
    hp.blank();
    
    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("--io_threads", "INT", "HTS I/O threads [" + std::to_string(GlobalOpts().IOthreads) + "]");
    hp.line("-d, --debug", "", "debug mode (forces threads=1)");
    hp.line("-h, --help", "", "show basic options");
    hp.line("-H, --advanced", "", "show advanced options");
    
    hp.blank();
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
                cfg.homq.batch_size = parse_size_arg_u32_(optarg, argc, argv, optind, "--batch");
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
                cfg.coordmap.max_fanout = parse_size_arg_u32_(optarg, argc, argv, optind, "--cm_max_fanout");
                break;
            }
            case 2003: {
                cfg.coordmap.min_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--cm_min_len");
                break;
            }
            case 2004: {
                cfg.coordmap.min_frac = std::max(0.0, std::atof(optarg));
                break;
            }
            case 2005: {
                cfg.coordmap.max_total_hits = parse_size_arg_u32_(optarg, argc, argv, optind, "--cm_max_hits");
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
    HelpPrinter hp;
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... -r FILE ... [options]\n\n"
        << "Align sequencing data to a GFA\n";

    hp.blank();

    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    hp.line("-r, --read", "FILE ...", "input sequencing data file(s)");
    hp.line("-p, --paf", "", "output PAF (default SAM)");
    hp.line("-o, --output", "FILE", "output alignments to file [stdout]");
    
    hp.blank();

    hp.section("Index options");
    hp.line("-k, --kmer", "INT", "k-mer length [" + std::to_string(GlobalOpts().kmerLen) + "]");
    hp.line("-w, --window", "INT", "minimizer window [" + std::to_string(GlobalOpts().minimizerW) + "]");
    
    hp.blank();

    hp.section("Alignment options");
    hp.line("-x, --preset", "STR", "hifi/ont/illumina/other [" + MapOpts().preset + "]");
    hp.line("-s, --secondary", "FLOAT", "min ratio secondary/primary [" + format_double_(MapOpts().sec_pri_ratio) + "]");
    hp.line("-N", "INT", "number of secondary to keep [" + std::to_string(MapOpts().sec_pri_num) + "]");
    
    hp.blank();

    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-d, --debug", "", "debug mode (forces threads=1)");
    hp.line("-h, --help", "", "show basic options");
    
    hp.blank();
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

        {"kmer",      required_argument, nullptr, 'k'},
        {"window",    required_argument, nullptr, 'w'},

        {"preset",    required_argument, nullptr, 'x'},
        {"secondary", required_argument, nullptr, 's'},
        {"N",         required_argument, nullptr, 'N'},

        {"threads",   required_argument, nullptr, 't'},
        {"debug",     no_argument,       nullptr, 'd'},
        {"help",      no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:r:po:k:w:x:s:N:t:dh";

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

            case 'k': {
                cfg.global.kmerLen = std::stoi(optarg);
                break;
            }
            case 'w': {
                cfg.global.minimizerW = std::stoi(optarg);
                break;
            }

            case 'x': {
                cfg.map.preset = optarg;
                break;
            }
            case 's': {
                cfg.map.sec_pri_ratio = std::stod(optarg);
                break;
            }
            case 'N': {
                cfg.map.sec_pri_num = std::stoi(optarg);
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
    HelpPrinter hp;
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... -e STR ... [options]\n\n"
        << "Build restriction cut-site index from a genome\n";
    
    hp.blank();
    
    hp.section("Input/Output");
    hp.line("-g, --genome", "FILE", "input genome FASTA");
    hp.line("-o, --output", "FILE", "output BED file (optional)");
    
    hp.blank();
    
    hp.section("Scan options");
    hp.line("-e, --enzyme", "STR", "enzyme recognition site with cut mark '^', can take multiple arguments after one -e");
    hp.note(" * e.g. G^AATTC  A^AGCTT  G^ANTC");
    hp.line("--rc", "", "scan reverse-complement motifs");
    
    hp.blank();
    
    hp.section("General options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-h, --help", "", "show basic options");
    
    hp.blank();
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

void help_split(char** argv) {
    HelpPrinter hp;
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... [options]\n\n"
        << "Split graph by connected components\n";

    hp.blank();

    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    hp.line("-n, --name", "STR ...", "sample names; otherwise reuse liftasm tags or use filename prefixes");
    hp.line("-p, --prefix", "STR", "output prefix [" + SplitOpts().prefix + "]");
    hp.note(" * <prefix>.compN.{gfa,noseq.gfa}");
    
    hp.blank();
    
    hp.section("General options");
    hp.line("-h, --help", "", "show basic options");
    
    hp.blank();
}

AppConfig main_split(int argc, char** argv) {
    if (argc < 3) { help_split(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::split;

    const struct option long_opts[] = {
        {"gfa",    required_argument, nullptr, 'g'},
        {"name",   required_argument, nullptr, 'n'},
        {"prefix", required_argument, nullptr, 'p'},
        {"help",   no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:n:p:h";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.split.gfaFiles.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.split.gfaFiles.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'n': {
                cfg.split.gfaNames.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.split.gfaNames.emplace_back(argv[optind]);
                    ++optind;
                }
                break;
            }
            case 'p': {
                cfg.split.prefix = optarg;
                break;
            }
            case 'h': {
                help_split(argv);
                std::exit(0);
            }
            default: {
                help_split(argv);
                std::exit(1);
            }
        }
    }

    validate_and_print(argc, argv, cfg);

    return cfg;
}

void help_augment(char** argv, bool advanced) {
    HelpPrinter hp(std::cerr, 20, 13);
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE ... -v FILE ... [options]\n\n"
        << "Augment GFA segments with VCF alleles as variant bubbles\n";

    hp.blank();

    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE ...", "input GFA file(s)");
    hp.line("-n, --name", "STR ...", "sample names; otherwise reuse liftasm tags or use filename prefixes");
    hp.line("-v, --vcf", "FILE ...", "input VCF or VCF.GZ file(s)");
    hp.line("-p, --prefix", "STR", "output prefix [" + AugmentOpts().prefix + "]");
    hp.note(" * <prefix>.augment.{gfa,noseq.gfa,map}");
    hp.line("--tag", "STR", "VT:Z tag value written on augmented ALT nodes [" + AugmentOpts().tag + "]");
    
    hp.blank();
    
    hp.section("Index options");
    hp.line("-k, --kmer", "INT", "k-mer length [" + std::to_string(GlobalOpts().kmerLen) + "]");
    hp.line("-w, --window", "INT", "minimizer window [" + std::to_string(GlobalOpts().minimizerW) + "]");
    
    hp.blank();
    
    hp.section("Align options");
    hp.line("-x, --preset", "STR", "mapping preset: asm5/asm10/asm20/sr/lr:hq [" + CollapseOpts().mm2_preset + "]");
    if (advanced) hp.line("--use_wfa", "", "whether to use WFA for alignment, if not use mm2 (beta)");
    hp.line("-z, --zdrop", "INT", "Z-drop score [" + format_size_arg_(MapOpts().extendOpts.dyn_zdrop) + "]");
    hp.line("--min_match", "FLOAT", "minimum match fraction in the alignment CIGAR [" + format_double_(CollapseOpts().min_match_ratio) + "]");
    hp.line("--min_ali_ratio", "FLOAT", "minimum aligned-length fraction relative to the sequence length [" + format_double_(CollapseOpts().min_ali_ratio) + "]");
    hp.line("--min_mapq", "INT", "minimum mapping quality (MAPQ) to keep (no larger than 60) [" + std::to_string(int(CollapseOpts().min_mapq)) + "]");
    
    hp.blank();
    
    hp.section("Collapse options");
    hp.line("--min_eq", "INT", "minimum match length to add cut points at segment [" + std::to_string(CollapseOpts().min_eq) + "]");
    hp.line("--max_iters", "INT", "maximum iterations for cut point propagation [" + format_size_arg_(CollapseOpts().max_iters) + "]");
    if (advanced) {
        hp.line("--abnormal_cut", "INT,INT", "prune abnormal cut points after propagation as MAX_LEN,MIN_COUNT [" + std::to_string(CollapseOpts().max_abnormal_cut_len) + "," + std::to_string(CollapseOpts().min_abnormal_cut_count) + "]");
        hp.note(" * e.g. 0-3-6-9-100 -> 0-9-100 when 4,2 due to consecutive short segments");
        hp.line("--min_trans_len", "INT", "minimum interval length to allow transitive replacement expansion [" + format_size_arg_(CollapseOpts().min_trans_len) + "]");
    }
    
    hp.blank();
    
    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-d, --debug", "", "debug mode");
    hp.line("-h, --help", "", "show basic options");
    hp.line("-H, --advanced", "", "show advanced options");
    
    hp.blank();
}

AppConfig main_augment(int argc, char** argv) {
    if (argc < 3) { help_augment(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::augment;
    const struct option long_opts[] = {
        {"gfa",           required_argument, nullptr, 'g'},
        {"name",          required_argument, nullptr, 'n'},
        {"vcf",           required_argument, nullptr, 'v'},
        {"prefix",        required_argument, nullptr, 'p'},
        {"tag",           required_argument, nullptr, 0001},

        {"kmer",          required_argument, nullptr, 'k'},
        {"window",        required_argument, nullptr, 'w'},

        {"preset",        required_argument, nullptr, 'x'},
        {"use_wfa",       no_argument,       nullptr, 2001},
        {"zdrop",         required_argument, nullptr, 'z'},
        {"min_match",     required_argument, nullptr, 2002},
        {"min_ali_ratio", required_argument, nullptr, 2003},
        {"min_mapq",      required_argument, nullptr, 2004},

        {"min_eq",        required_argument, nullptr, 3001},
        {"max_iters",     required_argument, nullptr, 3002},
        {"abnormal_cut",  required_argument, nullptr, 3003},
        {"min_trans_len", required_argument, nullptr, 3004},

        {"threads",       required_argument, nullptr, 't'},
        {"debug",         no_argument,       nullptr, 'd'},
        {"help",          no_argument,       nullptr, 'h'},
        {"advanced",      no_argument,      nullptr, 'H'},
        {0,0,0,0}
    };
    const char* short_opts = "g:n:v:p:k:w:x:z:t:dhH";

    int idx = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.augment.gfa_files.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) cfg.augment.gfa_files.emplace_back(argv[optind++]);
                break;
            }
            case 'n': {
                cfg.augment.gfa_names.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) cfg.augment.gfa_names.emplace_back(argv[optind++]);
                break;
            }
            case 'v': {
                cfg.augment.vcf_files.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) cfg.augment.vcf_files.emplace_back(argv[optind++]);
                break;
            }
            case 'p': {
                cfg.augment.prefix = optarg; 
                break;
            }
            case 0001: {
                cfg.augment.tag = optarg; 
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

            case 'x': {
                cfg.collapse.mm2_preset = optarg; 
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
            case 2003: {
                cfg.collapse.min_ali_ratio = std::max(0.0, std::stod(optarg)); 
                break;
            }
            case 2004: {
                const unsigned long value = std::stoul(optarg);
                if (value > 60) { error_stream() << "--min_mapq must be <= 60\n"; std::exit(1); }
                cfg.collapse.min_mapq = static_cast<uint8_t>(value);
                break;
            }

            case 3001: {
                cfg.collapse.min_eq = std::stoi(optarg); 
                break;
            }
            case 3002: {
                cfg.collapse.max_iters = std::stoi(optarg); 
                break;
            }
            case 3003: {
                parse_abnormal_seg_(optarg, cfg.collapse); 
                break;
            }
            case 3004: {
                cfg.collapse.min_trans_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--min_trans_len"); 
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
                help_augment(argv, false); 
                std::exit(0);
            }
            case 'H': {
                help_augment(argv, true); 
                std::exit(0);
            }
            default: {
                help_augment(argv); 
                std::exit(1);
            }
        }
    }
    validate_and_print(argc, argv, cfg);
    finalize_opt_cfg(cfg);
    return cfg;
}

void help_clean(char** argv) {
    HelpPrinter hp(std::cerr, 20, 13);
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g UTG.gfa -c CTG1.gfa CTG2.gfa [options]\n\n"
        << "Remove weak UTG links using read components from contig GFAs (beta)\n";

    hp.blank();

    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE", "input UTG GFA containing read A-lines");
    hp.line("-c, --ctg", "FILE ...", "contig GFA file(s) used together to build read components");
    hp.line("-p, --prefix", "STR", "output prefix [" + CleanOpts().prefix + "]");
    hp.note(" * <prefix>.clean.{gfa,noseq.gfa,reads.tsv}");

    hp.blank();

    hp.section("Cleaning options");
    hp.line("--min_reads", "INT", "minimum labeled reads required at each link endpoint [" + std::to_string(CleanOpts().min_reads) + "]");
    hp.line("--min_purity", "FLOAT", "minimum dominant-component fraction at each endpoint [" + format_double_(CleanOpts().min_purity) + "]");
    hp.line("--min_comp_overlap", "FLOAT", "minimum UTG-node containment between connected contig components [" + format_double_(CleanOpts().min_comp_overlap) + "]");
    
    hp.blank();

    hp.section("General Options");
    hp.line("-d, --debug", "", "debug mode");
    hp.line("-h, --help", "", "show options");
    
    hp.blank();
}

AppConfig main_clean(int argc, char** argv) {
    if (argc < 3) { help_clean(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::clean;

    const struct option long_opts[] = {
        {"gfa",              required_argument, nullptr, 'g'},
        {"ctg",              required_argument, nullptr, 'c'},
        {"prefix",           required_argument, nullptr, 'p'},

        {"min_reads",        required_argument, nullptr, 1001},
        {"min_purity",       required_argument, nullptr, 1002},
        {"min_comp_overlap", required_argument, nullptr, 1003},

        {"debug",            no_argument,  nullptr, 'd'},
        {"help",             no_argument,   nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:c:p:dh";

    int idx = 0;
    int c = 0;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.clean.utg_file = optarg; 
                break;
            }
            case 'c':{
                cfg.clean.ctg_files.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.clean.ctg_files.emplace_back(argv[optind++]);
                }
                break;
            }
            case 'p': {
                cfg.clean.prefix = optarg; 
                break;
            }

            case 1001: {
                cfg.clean.min_reads = static_cast<uint32_t>(std::stoul(optarg)); 
                break;
            }
            case 1002: {
                cfg.clean.min_purity = std::stod(optarg); 
                break;
            }
            case 1003: {
                cfg.clean.min_comp_overlap = std::stod(optarg); 
                break;
            }

            case 'd': {
                cfg.global.debug = true; 
                break;
            }
            case 'h': {
                help_clean(argv); 
                std::exit(0);
            }
            default: {
                help_clean(argv); 
                std::exit(1);
            }
        }
    }

    validate_and_print(argc, argv, cfg);
    return cfg;
}

void help_gapfill(char** argv) {
    HelpPrinter hp(std::cerr, 20, 13);
    std::cerr
        << "Usage: " << argv[0] << " " << argv[1] << " -g FILE [options]\n\n"
        << "Fill supported contig gaps with source-aware node walks between trusted graph anchors.\n";

    hp.blank();

    hp.section("Input/Output");
    hp.line("-g, --gfa", "FILE", "final contig-mode collapse GFA");
    hp.line("-p, --prefix", "STR", "output prefix [" + GapfillOpts().prefix + "]");
    hp.note(" * <prefix>.<sample>.hap{1,2}.gapfill.{gfa,noseq.gfa}");
    hp.note(" * <prefix>.primary.hap{1,2}.gapfill.{gfa,noseq.gfa}");
    hp.note(" * <prefix>.primary.hap{1,2}.gapfill.bubbles.vcf");
    hp.note(" * <prefix>.primary.hap{1,2}.lowQ.gapfill.{gfa,noseq.gfa}");
    hp.note(" * <prefix>.gapfill.{tsv,relocations.tsv,html}");

    hp.blank();

    hp.section("Index options");
    hp.line("-k, --kmer", "INT", "minimizer k-mer size [" + std::to_string(GlobalOpts().kmerLen) + "]");
    hp.line("-w, --window", "INT", "minimizer window size [" + std::to_string(GlobalOpts().minimizerW) + "]");

    hp.blank();

    hp.section("Align options");
    hp.line("-x, --preset", "STR", "mapping preset: asm5/asm10/asm20/sr/lr:hq [" + CollapseOpts().mm2_preset + "]");
    hp.line("-z, --zdrop", "INT", "minimap2 dynamic z-drop [" + std::to_string(MapOpts().zdrop) + "]");
    hp.line("--min_match", "FLOAT", "minimum matching-base fraction in an alignment [" + format_double_(GapfillOpts().min_match) + "]");
    hp.line("--min_ali_ratio", "FLOAT", "minimum aligned fraction of the shorter sequence [" + format_double_(GapfillOpts().min_ali_ratio) + "]");
    hp.line("--min_mapq", "INT", "minimum mapping quality [" + std::to_string(GapfillOpts().min_mapq) + "]");

    hp.blank();

    hp.section("Bubble options");
    hp.line("--depth", "INT", "maximum graph distance searched for one bubble or gap walk [" + format_size_arg_(GapfillOpts().max_depth) + "]");
    hp.line("--paths", "INT", "maximum alternative paths kept for one bubble [" + format_size_arg_(GapfillOpts().max_paths) + "]");
    hp.line("--DFS_guard", "INT", "stop one bubble or gap walk after this many visited states [" + format_size_arg_(GapfillOpts().DFS_guard) + "]");
    hp.line("--path_sim", "FLOAT", "minimum minimizer Jaccard for path clustering [" + format_double_(GapfillOpts().path_sim) + "]");
    hp.line("--stall_rounds", "INT", "stop after this many rounds find no new path [" + format_size_arg_(GapfillOpts().stall_round_limit) + "]");

    hp.blank();

    hp.section("Phase options");
    hp.line("--full_phase", "STR ...", "samples assembled with chromosome-scale phase information, e.g. trio-binned or Hi-C phased");
    hp.note(" * Full-phase targets use only full-phase samples as spanning evidence and may exchange haplotypes");
    hp.line("--phase_len", "INT", "use a bubble only when all internal paths are at most this many bp [" + std::to_string(GapfillOpts().phase_path_len) + "]");
    hp.line("--phase_win", "INT", "local window used once for bubble phase evidence [" + format_size_arg_(GapfillOpts().phase_win) + "]");
    
    hp.blank();

    hp.section("Gap options");
    hp.line("--min_contig", "INT", "do not use contigs shorter than this [" + format_size_arg_(GapfillOpts().min_contig) + "]");
    hp.line("--max_gap", "INT", "largest graph gap that may be filled [" + format_size_arg_(GapfillOpts().max_gap) + "]");
    hp.line("--max_overlap", "FLOAT", "largest allowed overlap between the two target contigs [" + format_double_(GapfillOpts().max_overlap) + "]");
    hp.line("--min_overlap", "INT", "bridge must share at least this many bp with each target contig [" + format_size_arg_(GapfillOpts().min_overlap) + "]");
    hp.line("--min_similarity", "FLOAT", "minimum shared-node similarity for homologous spanning evidence [" + format_double_(GapfillOpts().min_similarity) + "]");

    hp.blank();

    hp.section("Misassembly options");
    hp.line("--ms_len", "INT", "minimum low-support sequence at a contig end [" + format_size_arg_(GapfillOpts().ms_len) + "]");
    hp.line("--ms_sim", "FLOAT", "minimum fraction covered by one other component [" + format_double_(GapfillOpts().ms_sim) + "]");
    hp.line("--ms_haps", "INT|auto", "minimum supporting haplotypes in each component [auto]");
    hp.note(" * auto: disabled for <=2 haplotypes; 1 for 3-4; 2 otherwise");

    hp.blank();

    hp.section("Deduplication options");
    hp.line("--dedup_sim", "FLOAT", "minimum aligned sequence fraction for removing redundant contigs or components [" + format_double_(GapfillOpts().dedup_similarity) + "]");
    hp.line("--dedup_component", "INT", "largest primary component checked against larger contigs of the same haplotype [" + format_size_arg_(GapfillOpts().dedup_component) + "]");
    
    hp.blank();

    hp.section("General Options");
    hp.line("-t, --threads", "INT", "number of threads [" + std::to_string(GlobalOpts().threads) + "]");
    hp.line("-d, --debug", "", "debug mode");
    hp.line("-h, --help", "", "show options");
    
    hp.blank();
}

AppConfig main_gapfill(int argc, char** argv) {
    if (argc < 3) { help_gapfill(argv); std::exit(1); }

    AppConfig cfg;
    cfg.mode = ToolMode::gapfill;

    const struct option long_opts[] = {
        {"gfa",                    required_argument, nullptr, 'g'},
        {"prefix",                 required_argument, nullptr, 'p'},

        {"kmer",                   required_argument, nullptr, 'k'},
        {"window",                 required_argument, nullptr, 'w'},

        {"preset",                 required_argument, nullptr, 'x'},
        {"zdrop",                  required_argument, nullptr, 'z'},
        {"min_match",              required_argument, nullptr, 2001},
        {"min_ali_ratio",          required_argument, nullptr, 2002},
        {"min_mapq",               required_argument, nullptr, 2003},

        {"depth",                  required_argument, nullptr, 3001},
        {"paths",                  required_argument, nullptr, 3002},
        {"DFS_guard",              required_argument, nullptr, 3003},
        {"path_sim",               required_argument, nullptr, 3004},
        {"stall_rounds",           required_argument, nullptr, 3005},

        {"full_phase",             required_argument, nullptr, 4001},
        {"phase_len",              required_argument, nullptr, 4002},
        {"phase_win",              required_argument, nullptr, 4003},

        {"min_contig",             required_argument, nullptr, 5002},
        {"max_gap",                required_argument, nullptr, 5003},
        {"max_overlap",            required_argument, nullptr, 5004},
        {"min_overlap",            required_argument, nullptr, 5005},
        {"min_similarity",         required_argument, nullptr, 5006},

        {"ms_len",                 required_argument, nullptr, 6001},
        {"ms_sim",                 required_argument, nullptr, 6002},
        {"ms_haps",                required_argument, nullptr, 6003},

        {"dedup_sim",              required_argument, nullptr, 7001},
        {"dedup_component",        required_argument, nullptr, 7002},

        {"threads",                required_argument, nullptr, 't'},
        {"debug",                  no_argument,       nullptr, 'd'},
        {"help",                   no_argument,       nullptr, 'h'},
        {0,0,0,0}
    };
    const char* short_opts = "g:p:k:w:x:z:t:dh";

    int idx = 0;
    int c = 0;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &idx)) != -1) {
        switch (c) {
            case 'g': {
                cfg.gapfill.gfa_file = optarg; 
                break;
            }
            case 'p': {
                cfg.gapfill.prefix = optarg; 
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
            case 'x': {
                cfg.gapfill.mm2_preset = optarg;
                break;
            }
            case 'z': {
                cfg.map.zdrop = std::max(1, std::stoi(optarg)); 
                cfg.map.zdrop_set = true; 
                break;
            }
            case 2001: {
                cfg.gapfill.min_match = std::stod(optarg); 
                break;
            }
            case 2002: {
                cfg.gapfill.min_ali_ratio = std::stod(optarg); 
                break;
            }
            case 2003: {
                const unsigned long value = std::stoul(optarg);
                if (value > 60) {
                    error_stream() << "--min_mapq must be <= 60\n";
                    std::exit(1);
                }
                cfg.gapfill.min_mapq = static_cast<uint8_t>(value);
                break;
            }
            case 3001: {
                cfg.gapfill.max_depth = parse_size_arg_u32_(optarg, argc, argv, optind, "--depth");
                break;
            }
            case 3002: {
                cfg.gapfill.max_paths = parse_size_arg_u16_(optarg, argc, argv, optind, "--paths");
                break;
            }
            case 3003: {
                cfg.gapfill.DFS_guard = parse_size_arg_u64_(optarg, argc, argv, optind, "--DFS_guard");
                break;
            }
            case 3004: {
                cfg.gapfill.path_sim = std::stod(optarg);
                break;
            }
            case 3005: {
                cfg.gapfill.stall_round_limit = parse_size_arg_u32_(optarg, argc, argv, optind, "--stall_rounds");
                break;
            }

            case 4001: {
                cfg.gapfill.full_phase_samples.emplace_back(optarg);
                while (optind < argc && !is_flag_(argv[optind])) {
                    cfg.gapfill.full_phase_samples.emplace_back(argv[optind++]);
                }
                break;
            }
            case 4002: {
                cfg.gapfill.phase_path_len = parse_size_arg_u32_(optarg, argc, argv, optind, "--phase_len");
                break;
            }
            case 4003: {
                cfg.gapfill.phase_win = parse_size_arg_u64_(optarg, argc, argv, optind, "--phase_win");
                break;
            }
            case 5002: {
                cfg.gapfill.min_contig = parse_size_arg_u64_(optarg, argc, argv, optind, "--min_contig");
                break;
            }
            case 5003: {
                cfg.gapfill.max_gap = parse_size_arg_u64_(optarg, argc, argv, optind, "--max_gap");
                break;
            }
            case 5004: {
                cfg.gapfill.max_overlap = std::stod(optarg);
                break;
            }
            case 5005: {
                cfg.gapfill.min_overlap = parse_size_arg_u64_(optarg, argc, argv, optind, "--min_overlap");
                break;
            }
            case 5006: {
                cfg.gapfill.min_similarity = std::stod(optarg);
                break;
            }
            case 6001: {
                cfg.gapfill.ms_len = parse_size_arg_u64_(optarg, argc, argv, optind, "--ms_len");
                break;
            }
            case 6002: {
                cfg.gapfill.ms_sim = std::stod(optarg);
                break;
            }
            case 6003: {
                cfg.gapfill.ms_haps = std::string(optarg) == "auto" ? 0 :
                    parse_size_arg_u32_(optarg, argc, argv, optind, "--ms_haps");
                break;
            }

            case 7001: {
                cfg.gapfill.dedup_similarity = std::stod(optarg);
                break;
            }
            case 7002: {
                cfg.gapfill.dedup_component = parse_size_arg_u64_(optarg, argc, argv, optind, "--dedup_component");
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
                help_gapfill(argv); 
                std::exit(0);
            }
            default: {
                help_gapfill(argv); 
                std::exit(1);
            }
        }
    }

    std::sort(cfg.gapfill.full_phase_samples.begin(), cfg.gapfill.full_phase_samples.end());
    cfg.gapfill.full_phase_samples.erase(
        std::unique(cfg.gapfill.full_phase_samples.begin(), cfg.gapfill.full_phase_samples.end()),
        cfg.gapfill.full_phase_samples.end()
    );

    validate_and_print(argc, argv, cfg);
    finalize_opt_cfg(cfg);
    return cfg;
}

void finalize_opt_cfg(AppConfig& cfg) {
    opt::InitParams params;
    params.k = cfg.global.kmerLen;
    params.w = cfg.global.minimizerW;
    params.sec_pri_ratio = cfg.map.sec_pri_ratio;
    params.sec_pri_num = cfg.map.sec_pri_num;
    params.out_paf = cfg.map.outPAF;
    params.threads = cfg.global.threads;

    opt::AlignmentOptions options = opt::init_opts(params);
    cfg.map.chainOpts = std::move(options.chain);
    cfg.map.anchorOpts = std::move(options.anchor);
    cfg.map.extendOpts = std::move(options.extend);
    cfg.map.alignOpts = std::move(options.align);

    if (cfg.mode == ToolMode::collapse || cfg.mode == ToolMode::deoverlap || cfg.mode == ToolMode::augment || cfg.mode == ToolMode::gapfill) {
        opt::Preset::map_asm_5(cfg.map.chainOpts, cfg.map.anchorOpts, cfg.map.extendOpts, cfg.map.alignOpts);
    }

    if (cfg.map.zdrop_set) {
        cfg.map.extendOpts.dyn_zdrop = cfg.map.zdrop;
    }
}
