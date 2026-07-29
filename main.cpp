// g++ main.cpp src/*.cpp -o liftasm -std=c++20 -march=native -O3 -fopenmp -I /home/zd233/zd233/00-software/00-project/liftasm/minimap2 /home/zd233/zd233/00-software/00-project/liftasm/minimap2/libminimap2.a /home/zd233/zd233/00-software/WFA2-lib/build/libwfa2cpp.a /home/zd233/zd233/00-software/WFA2-lib/build/libwfa2.a /home/zd233/zd233/00-software/htslib-1.22.1/build/libhts.a -lz -lbz2 -lcrypto -lssl -lpthread -ldl -lstdc++fs

#include <iostream>
#include <string>
#include <iomanip>

#include "include/OptionParser.hpp"
#include "include/OptionParser_helper.hpp"
#include "include/options.hpp"
#include "include/gfa_parser.hpp"
#include "include/gfa_depth.hpp"
#include "include/gfa_seq.hpp"
#include "include/gfa_deoverlapper.hpp"
#include "include/gfa_collapser.hpp"
#include "include/gfa_ctg_collapser.hpp"
#include "include/gfa_bubble.hpp"
#include "include/gfa_gfa2fa.hpp"
#include "include/file2map.hpp"
#include "include/liftover.hpp"
#include "include/coordmap.hpp"
#include "include/MapqBoost.hpp"
#include "include/restriction_sites.hpp"
#include "include/gfa_split.hpp"
#include "include/gfa_augment.hpp"
#include "include/gfa_clean.hpp"
#include "include/aligner.hpp"
#include "include/sys.hpp"

using namespace wfa;

static inline bool is_top_help_flag(const char* s) { return (std::string(s) == "-h" || std::string(s) == "--help"); }
static inline bool is_top_help_advanced_flag(const char* s) { return (std::string(s) == "-H" || std::string(s) == "--advanced"); }
static inline bool is_top_ver_flag(const char* s) { return (std::string(s) == "-v" || std::string(s) == "--version"); }

int main(int argc, char** argv) {
    if (argc < 2) { help(argv); return 1; }
    if (is_top_help_flag(argv[1])) { help(argv); return 0; }
    if (is_top_help_advanced_flag(argv[1])) { help(argv, true); return 0; }
    if (is_top_ver_flag(argv[1]))  { std::cerr << program::version << "\n"; return 0; }

    // Informations
    const std::string command_line = program::cmdline(argc, argv);
    log_stream() << "Version: " << program::version << "\n";
    log_stream() << "Data: " << getTime() << "\n";
    log_stream() << "Command: " << command_line << "\n" << "\n";

    // Dispatch by subcommand
    const std::string sub = argv[1];

    // timing
    double realtime0 = realtime();

    if (sub == "stat") {
        AppConfig cfg = main_stat(argc, argv);
        GfaGraph G;
        G.load_from_GFA(cfg.stat.gfaFiles);
    } else if (sub == "seq") {
        AppConfig cfg = main_seq(argc, argv);
        GfaSeq G;
        G.load_from_GFA(cfg.seq.gfaFiles);
        std::vector<std::string> paths, seqs;
        G.extract_from_file(cfg.seq.pathFile, paths, seqs);
        G.save_to_file(cfg.seq.outFile, paths, seqs);
    } else if (sub == "gfa2fa") {
        AppConfig cfg = main_gfa2fa(argc, argv);
        Gfa2fa G(cfg.gfa2fa.min_len_xbp, cfg.gfa2fa.extend_ybp, cfg.gfa2fa.wrap_width, cfg.gfa2fa.skip_unknown);
        G.load_from_GFA(cfg.gfa2fa.gfaFiles);
        G.dump_to_file(cfg.gfa2fa.outFile);
    } else if (sub == "depth") {
        AppConfig cfg = main_depth(argc, argv);
        opt::ChainOpts  chainOpts;
        opt::AnchorOpts anchorOpts;
        opt::ExtendOpts extendOpts;
        opt::AlignOpts  alignOpts;
        init_opts(
            cfg.global.kmerLen, 
            cfg.global.minimizerW,
            cfg.map.sec_pri_ratio, 
            cfg.map.sec_pri_num,
            cfg.map.outPAF, 
            cfg.global.threads,
            chainOpts, 
            anchorOpts, 
            extendOpts, 
            alignOpts
        );

        // Load graph
        bool forbid_overlap = (!cfg.depth.reads.empty() || !cfg.depth.gafFile.empty());
        GfaDepth G(cfg.depth.min_mapq, cfg.depth.min_frac, cfg.depth.base_depth, forbid_overlap);
        G.load_from_GFA(cfg.depth.gfaFiles);

        if (cfg.depth.reads.size() > 0) {
            auto name_seqs = G.getSeqVec(cfg.global.kmerLen);
            mmidx::MinimizerIndex GIndex(name_seqs.names, name_seqs.seqs, name_seqs.right_seqs, chainOpts, anchorOpts);
            GIndex.build_mm(true);
            GIndex.print_index_stats();
            GIndex.count_depth(cfg.depth.reads);
            G.count_from_kmer(GIndex, cfg.depth.outFile);
        } else if (cfg.depth.gafFile.size() > 0) {
            G.count_from_gaf(cfg.depth.gafFile, cfg.depth.outFile);
        } else {
            G.count_from_A(cfg.depth.outFile);
        }
    } else if (sub == "bubble") {
        AppConfig cfg = main_bubble(argc, argv);
        GfaGraph G;
        G.load_from_GFA(cfg.bubble.gfaFiles, cfg.bubble.gfaNames);
        G.build_nodes_connectivity_index();
        G.build_vertex_topological_index();
        if (cfg.bubble.check_complex) {
            G.detect_complex_regions(
                cfg.bubble.cx_branch_degree, cfg.bubble.cx_hub_degree,
                cfg.bubble.cx_min_nodes, cfg.bubble.cx_min_branches, /*skip_comp=*/false
            );
        }

        GfaBubble::GfaBubbleFinder finder(
            G, 
            cfg.bubble.max_depth, 
            cfg.bubble.max_paths, 
            cfg.bubble.DFS_guard, 
            cfg.bubble.path_diff, 
            cfg.bubble.stall_round_limit, 
            /*skip_comp=*/false, 
            cfg.bubble.keep_nested, 
            cfg.bubble.same_sim, 
            cfg.bubble.same_min_len,
            cfg.bubble.diff_min_src, 
            cfg.bubble.diff_sim, 
            cfg.bubble.diff_min_len, 
            cfg.bubble.homo_num, 
            cfg.bubble.homo_k, 
            cfg.bubble.homo_w, 
            cfg.bubble.homo_extend_bp, 
            cfg.bubble.homo_bloom_bits, 
            cfg.bubble.homo_bloom_hash,
            cfg.global.threads
        );

        finder.find_bubbles();
        finder.find_homologous_paths();

        if (DEBUG_ENABLED) finder.print_bubbles();
        if (DEBUG_ENABLED) finder.print_homologous_paths();

        finder.save_bubble_as_gfa(cfg.bubble.out_prefix + ".bubbles.gfa", cfg.bubble.min_len, cfg.bubble.min_num, /*write_seq=*/true, command_line);
        finder.save_bubble_as_gfa(cfg.bubble.out_prefix + ".bubbles.noseq.gfa", cfg.bubble.min_len, cfg.bubble.min_num, /*write_seq=*/false, command_line);
        if (cfg.bubble.write_vcf) {
            finder.save_bubble_as_vcf(cfg.bubble.out_prefix, cfg.bubble.paf_file, cfg.bubble.ref_file, cfg.bubble.ali_min_mapq, cfg.bubble.ali_min_len);
        }
    } else if (sub == "deoverlap") {
        AppConfig cfg = main_deoverlap(argc, argv);

        GfaDeoverlapper G(
            cfg.collapse.min_eq, 
            cfg.collapse.min_match_ratio, 
            cfg.collapse.min_ali_ratio,
            cfg.collapse.min_mapq,
            cfg.collapse.trim_min_len,
            cfg.collapse.trim_max_overlap,
            cfg.collapse.max_iters, 
            cfg.collapse.max_abnormal_cut_len, 
            cfg.collapse.min_abnormal_cut_count, 
            cfg.collapse.min_trans_len,
            cfg.collapse.mm2_preset
        );

        G.set_opts(
            cfg.map.chainOpts, 
            cfg.map.anchorOpts, 
            cfg.map.extendOpts, 
            cfg.map.alignOpts, 
            cfg.map.use_wfa
        );
        G.load_from_GFA(cfg.collapse.gfaFiles, cfg.collapse.gfaNames);
        G.print_graph_stats();

        G.deoverlap(cfg.collapse.prefix);

        G.save_to_disk(cfg.collapse.prefix + ".deoverlap.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/true, command_line);
        G.save_to_disk(cfg.collapse.prefix + ".deoverlap.noseq.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/false, command_line);
    } else if (sub == "collapse") {
        AppConfig cfg = main_collapse(argc, argv);
        const bool ctg_mode = !cfg.collapse.hap1Files.empty();
        std::vector<std::string> cur_gfa_files;
        std::vector<std::string> cur_gfa_names;

        if (ctg_mode) {
            GfaCtgCollapser G(
                cfg.collapse.min_jaccards.front(),
                cfg.collapse.min_eqs.front(),
                cfg.collapse.min_match_ratios.front(),
                cfg.collapse.min_ali_ratios.front(),
                cfg.collapse.min_mapq,
                cfg.collapse.trim_min_len,
                cfg.collapse.trim_max_overlap,
                cfg.collapse.max_iters,
                cfg.bubble.homo_k,
                cfg.bubble.homo_w,
                cfg.collapse.all_pair_len,
                cfg.collapse.repeat_mask_min_len,
                cfg.collapse.repeat_mask_max_period,
                cfg.collapse.repeat_mask_max_mismatch,
                cfg.collapse.max_abnormal_cut_len,
                cfg.collapse.min_abnormal_cut_count,
                cfg.collapse.min_trans_len,
                cfg.collapse.mm2_preset
            );
            G.set_opts(
                cfg.map.chainOpts,
                cfg.map.anchorOpts,
                cfg.map.extendOpts,
                cfg.map.alignOpts,
                cfg.map.use_wfa
            );
            G.set_component_filter(cfg.collapse.ctg_min_coverage, cfg.collapse.ctg_end_fraction);
            cur_gfa_files = G.collapse_samples(
                cfg.collapse.hap1Files,
                cfg.collapse.hap2Files,
                cfg.collapse.gfaNames,
                cfg.collapse.vcfFiles,
                cfg.collapse.prefix,
                cfg.collapse.ctg_min_reads,
                cfg.collapse.ctg_short_len,
                cfg.collapse.ctg_min_len,
                cfg.collapse.ctg_anchor_only,
                cfg.collapse.paf,
                command_line
            );
        } else {
            cur_gfa_files = cfg.collapse.gfaFiles;
            cur_gfa_names = cfg.collapse.gfaNames;
        }

        auto iter_value = [](const auto& xs, size_t iter) {
            return xs[std::min(iter, xs.size() - 1)];
        };

        for (size_t iter = 0; iter < cfg.collapse.iterations; ++iter) {
            const int min_eq = iter_value(cfg.collapse.min_eqs, iter);
            const double min_jaccard = iter_value(cfg.collapse.min_jaccards, iter);
            const double min_match_ratio = iter_value(cfg.collapse.min_match_ratios, iter);
            const double min_ali_ratio = iter_value(cfg.collapse.min_ali_ratios, iter);
            const double same_sim = iter_value(cfg.collapse.same_sims, iter);
            const uint32_t same_min_len = iter_value(cfg.collapse.same_min_lens, iter);
            const uint32_t diff_min_src = iter_value(cfg.collapse.diff_min_srcs, iter);
            const double diff_sim = iter_value(cfg.collapse.diff_sims, iter);
            const uint32_t diff_min_len = iter_value(cfg.collapse.diff_min_lens, iter);

            log_stream() << "Collapse iteration " << (iter + 1) << "/" << cfg.collapse.iterations << " ...\n";
            log_stream() << "  - min_eq       : " << format_size_arg_(min_eq) << "\n";
            log_stream() << "  - min_jaccard  : " << format_double_(min_jaccard) << "\n";
            log_stream() << "  - min_match    : " << format_double_(min_match_ratio) << "\n";
            log_stream() << "  - min_ali_ratio: " << format_double_(min_ali_ratio) << "\n";
            log_stream() << "  - same_sim     : " << format_double_(same_sim) << "\n";
            log_stream() << "  - same_len     : " << format_size_arg_(same_min_len) << "\n";
            log_stream() << "  - diff_src_len : " << format_size_arg_(diff_min_src) << "\n";
            log_stream() << "  - diff_sim     : " << format_double_(diff_sim) << "\n";
            log_stream() << "  - diff_len     : " << format_size_arg_(diff_min_len) << "\n\n";

            GfaCtgCollapser G(
                min_jaccard,
                min_eq,
                min_match_ratio,
                min_ali_ratio,
                cfg.collapse.min_mapq,
                cfg.collapse.trim_min_len,
                cfg.collapse.trim_max_overlap,
                cfg.collapse.max_iters,
                cfg.bubble.homo_k,
                cfg.bubble.homo_w,
                cfg.collapse.all_pair_len,
                cfg.collapse.repeat_mask_min_len,
                cfg.collapse.repeat_mask_max_period,
                cfg.collapse.repeat_mask_max_mismatch,
                cfg.collapse.max_abnormal_cut_len,
                cfg.collapse.min_abnormal_cut_count,
                cfg.collapse.min_trans_len,
                cfg.collapse.mm2_preset
            );

            G.set_opts(
                cfg.map.chainOpts,
                cfg.map.anchorOpts,
                cfg.map.extendOpts,
                cfg.map.alignOpts,
                cfg.map.use_wfa
            );
            G.load_from_GFA(cur_gfa_files, cur_gfa_names);
            G.build_nodes_connectivity_index();
            G.build_vertex_topological_index();
            if (iter == 0 && cfg.bubble.check_complex) {
                G.detect_complex_regions(
                    cfg.bubble.cx_branch_degree,
                    cfg.bubble.cx_hub_degree,
                    cfg.bubble.cx_min_nodes,
                    cfg.bubble.cx_min_branches,
                    /*skip_comp=*/false
                );
            }

            GfaBubble::GfaBubbleFinder finder(
                G,
                cfg.bubble.max_depth,
                cfg.bubble.max_paths,
                cfg.bubble.DFS_guard,
                cfg.bubble.path_diff,
                cfg.bubble.stall_round_limit,
                /*skip_comp=*/false,
                /*keep_nested=*/false,
                same_sim,
                same_min_len,
                diff_min_src,
                diff_sim,
                diff_min_len,
                cfg.bubble.homo_num,
                cfg.bubble.homo_k,
                cfg.bubble.homo_w,
                cfg.bubble.homo_extend_bp,
                cfg.bubble.homo_bloom_bits,
                cfg.bubble.homo_bloom_hash,
                cfg.global.threads
            );

            finder.find_bubbles();
            const auto& bubbles = finder.get_bubbles();
            if (DEBUG_ENABLED) finder.print_bubbles();

            finder.find_homologous_paths();
            const auto& homologous_paths = finder.get_homologous_paths();
            if (DEBUG_ENABLED) finder.print_homologous_paths();

            // Mark complex
            std::vector<uint32_t> complex_segs = finder.collect_orientation_conflict_segments();
            G.mark_segments_complex(complex_segs);

            // Repeat normalization
            const bool is_last_iter = (iter + 1 == cfg.collapse.iterations);
            if (is_last_iter && cfg.collapse.repeat_norm_len > 0) {
                G.normalize_homopolymer_bubbles(
                    bubbles,
                    cfg.collapse.repeat_norm_len
                );

                G.disable_repeat_mask();
            }

            const std::string iter_prefix = cfg.collapse.prefix + ".iter" + std::to_string(iter + 1);

            G.collapse_homologous_seq(bubbles, homologous_paths, iter_prefix);
            if (G.getNumPaths() > 0) G.rebuild_component_paths();

            const std::string iter_gfa = iter_prefix + ".collapse.gfa";
            const std::string iter_noseq = iter_prefix + ".collapse.noseq.gfa";
            const std::string iter_bubble_noseq = iter_prefix + ".bubbles.noseq.gfa";

            const bool write_paths = G.getNumPaths() > 0;
            G.save_to_disk(iter_gfa, /*write_paths=*/write_paths, /*write_align=*/false, /*write_seq=*/true, command_line);
            G.save_to_disk(iter_noseq, /*write_paths=*/write_paths, /*write_align=*/false, /*write_seq=*/false, command_line);
            finder.save_bubble_as_gfa(iter_bubble_noseq, cfg.bubble.min_len, cfg.bubble.min_num, /*write_seq=*/false, command_line);

            cur_gfa_files.assign(1, iter_gfa);

            cur_gfa_names.assign(1, "");
        }
    } else if (sub == "file2map") {
        AppConfig cfg = main_file2map(argc, argv);
        mapconv::file_to_map_auto(cfg.file2map.inputFiles, cfg.file2map.outFile, cfg.file2map.paf_primary_only, cfg.file2map.min_len, cfg.file2map.min_mapq);
    } else if (sub == "liftover") {
        AppConfig cfg = main_liftover(argc, argv);
        liftover::RunOpts opt = liftover::set_opts(
            cfg.liftover.mapFiles, 
            cfg.liftover.bedFile, 
            cfg.liftover.referenceFile, 
            cfg.liftover.outFile,
            cfg.liftover.regex, 
            cfg.liftover.min_frac, 
            cfg.liftover.flank_win, 
            cfg.liftover.max_flank, 
            cfg.liftover.max_gap, 
            cfg.liftover.max_hit,
            cfg.liftover.pafFile, 
            cfg.liftover.min_mapq, 
            cfg.liftover.min_len,
            cfg.liftover.do_check, 
            cfg.liftover.win, 
            cfg.liftover.step, 
            cfg.liftover.max_examples, 
            cfg.coordmap.max_hops, 
            cfg.coordmap.max_fanout, 
            cfg.coordmap.min_len, 
            cfg.coordmap.min_frac, 
            cfg.coordmap.max_total_hits, 
            cfg.global.threads
        );
        std::vector<std::string> blocks = liftover::liftover(opt);
        liftover::save_liftover_results(opt.out_file, blocks);
    } else if (sub == "mapq_boost") {
        AppConfig cfg = main_mapq_boost(argc, argv);
        coordmap::CoordMap coormap_idx;
        coormap_idx.load(cfg.homq.mapFiles);

        // Restriction-site index (optional)
        rsite::Index rsite_idx;
        if (!cfg.ressit.genome_file.empty() || cfg.ressit.enzymes.size() > 0) {
            rsite_idx.build(cfg.ressit.genome_file, cfg.ressit.enzymes, cfg.ressit.scan_revcomp, cfg.global.threads);
        }

        mapqboost::MapqBooster booster(
            coormap_idx, 
            rsite_idx, 
            cfg.homq.batch_size, 
            cfg.homq.mapq_low, 
            cfg.homq.mapq_cap, 
            cfg.homq.name_check,
            cfg.coordmap.max_hops, 
            cfg.coordmap.max_fanout, 
            cfg.coordmap.min_len, 
            cfg.coordmap.min_frac, 
            cfg.coordmap.max_total_hits,
            cfg.homq.sub_ovlp_frac, 
            cfg.homq.K_mapq, 
            cfg.homq.K_rs, 
            cfg.homq.K_as, 
            cfg.homq.K_ml, 
            cfg.homq.K_nm, 
            cfg.homq.W_mapq, 
            cfg.homq.W_rs, 
            cfg.homq.W_as, 
            cfg.homq.W_ml, 
            cfg.homq.W_nm, 
            cfg.homq.close_as_eps,
            cfg.global.threads, 
            cfg.global.IOthreads
        );

        booster.run(cfg.homq.in_bam, cfg.homq.out_bam);
    } else if (sub == "align") {
        AppConfig cfg = main_align(argc, argv);
        opt::ChainOpts  chainOpts;
        opt::AnchorOpts anchorOpts;
        opt::ExtendOpts extendOpts;
        opt::AlignOpts  alignOpts;
        init_opts(
            cfg.global.kmerLen, 
            cfg.global.minimizerW,
            cfg.map.sec_pri_ratio, 
            cfg.map.sec_pri_num,
            cfg.map.outPAF, 
            cfg.global.threads,
            chainOpts, 
            anchorOpts, 
            extendOpts, 
            alignOpts
        );

        // Preset
        if      (cfg.map.preset == "ont")      opt::Preset::map_ont(chainOpts, anchorOpts, extendOpts, alignOpts);
        else if (cfg.map.preset == "hifi")     opt::Preset::map_hifi(chainOpts, anchorOpts, extendOpts, alignOpts);
        else if (cfg.map.preset == "illumina") opt::Preset::map_illumina(chainOpts, anchorOpts, extendOpts, alignOpts);
        else if (cfg.map.preset == "asm5")     opt::Preset::map_asm_5(chainOpts, anchorOpts, extendOpts, alignOpts);
        else                                   opt::Preset::map_other(chainOpts, anchorOpts, extendOpts, alignOpts);

        // Build index
        GfaGraph G;
        G.load_from_GFA(cfg.map.gfaFiles);
        G.build_nodes_connectivity_index();
        auto name_seqs = G.getSeqVec();
        mmidx::MinimizerIndex GIndex(name_seqs.names, name_seqs.seqs, name_seqs.right_seqs, chainOpts, anchorOpts);
        GIndex.build_mm();
        GIndex.print_index_stats();

        // Align
        aligner::Alignmenter A(GIndex, name_seqs.names, name_seqs.seqs, extendOpts, alignOpts);
        A.align(cfg.map.reads, cfg.map.outFile, command_line);
    } else if (sub == "ressit") {
        AppConfig cfg = main_res_cut(argc, argv);

        rsite::Index idx;
        idx.build(cfg.ressit.genome_file, cfg.ressit.enzymes, cfg.ressit.scan_revcomp, cfg.global.threads);

        if (!cfg.ressit.out_bed.empty()) {
            idx.save_bed(cfg.ressit.out_bed);
        }
    } else if (sub == "split") {
        AppConfig cfg = main_split(argc, argv);

        GfaSplitter G;
        G.load_from_GFA(cfg.split.gfaFiles, cfg.split.gfaNames);

        G.split_by_components(
            cfg.split.prefix,
            /*write_paths=*/true,
            /*write_align=*/true,
            /*skip_comp=*/true
        );
    } else if (sub == "augment") {
        AppConfig cfg = main_augment(argc, argv);
        GfaAugmenter G(
            cfg.collapse.min_eq,
            cfg.collapse.min_match_ratio,
            cfg.collapse.min_ali_ratio,
            cfg.collapse.min_mapq,
            cfg.collapse.trim_min_len,
            cfg.collapse.trim_max_overlap,
            cfg.collapse.max_iters,
            cfg.collapse.max_abnormal_cut_len,
            cfg.collapse.min_abnormal_cut_count,
            cfg.collapse.min_trans_len,
            cfg.collapse.mm2_preset
        );
        G.set_opts(cfg.map.chainOpts, cfg.map.anchorOpts, cfg.map.extendOpts,
                   cfg.map.alignOpts, cfg.map.use_wfa);
        G.load_from_GFA(cfg.augment.gfa_files, cfg.augment.gfa_names);
        G.augment(cfg.augment.vcf_files, cfg.augment.prefix, cfg.augment.tag);
        G.save_to_disk(cfg.augment.prefix + ".augment.gfa", /*write_paths=*/false, /*write_align=*/false,
                       /*write_seq=*/true, command_line);
        G.save_to_disk(cfg.augment.prefix + ".augment.noseq.gfa", /*write_paths=*/false, /*write_align=*/false,
                       /*write_seq=*/false, command_line);
    } else if (sub == "clean") {
        AppConfig cfg = main_clean(argc, argv);
        GfaCleaner G({
            cfg.clean.min_reads,
            cfg.clean.min_purity,
            cfg.clean.min_comp_overlap
        });
        G.load_from_GFA({cfg.clean.utg_file});
        G.clean(cfg.clean.ctg_files);
        G.save_read_placements(cfg.clean.prefix + ".clean.reads.tsv");
        G.save_to_disk(cfg.clean.prefix + ".clean.gfa", /*write_paths=*/true, /*write_align=*/true, /*write_seq=*/true, command_line);
        G.save_to_disk(cfg.clean.prefix + ".clean.noseq.gfa", /*write_paths=*/true, /*write_align=*/true, /*write_seq=*/false, command_line);
    } else {
        error_stream() << "Unknown subcommand: " << sub << "\n";
        help(argv);
        return 1;
    }

    log_stream() 
        << "Real time: " << std::fixed << std::setprecision(3)
        << (realtime() - realtime0) << " sec; CPU: " << cputime()
        << " sec; Peak RSS: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB\n";

    return 0;
}
