#include <iostream>
#include <string>
#include <iomanip>
#include <filesystem>
#include <memory>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../include/OptionParser.hpp"
#include "../include/OptionParser_helper.hpp"
#include "../include/options.hpp"
#include "../include/gfa_parser.hpp"
#include "../include/gfa_depth.hpp"
#include "../include/gfa_seq.hpp"
#include "../include/gfa_deoverlapper.hpp"
#include "../include/gfa_collapser.hpp"
#include "../include/gfa_ctg_collapser.hpp"
#include "../include/gfa_bubble.hpp"
#include "../include/gfa_gfa2fa.hpp"
#include "../include/file2map.hpp"
#include "../include/liftover.hpp"
#include "../include/coordmap.hpp"
#include "../include/MapqBoost.hpp"
#include "../include/restriction_sites.hpp"
#include "../include/gfa_split.hpp"
#include "../include/gfa_augment.hpp"
#include "../include/gfa_clean.hpp"
#include "../include/gfa_gapfill.hpp"
#include "../include/aligner.hpp"
#include "../include/CommandRunner.hpp"

using namespace wfa;

namespace {

// Complete options shared with the global thread setting.
BubbleOpts bubble_options(const AppConfig& cfg)
{
    BubbleOpts options = cfg.bubble;
    options.threads = static_cast<uint32_t>(cfg.global.threads);
    return options;
}

opt::AlignmentOptions alignment_options(const AppConfig& cfg)
{
    opt::InitParams params;
    params.k = cfg.global.kmerLen;
    params.w = cfg.global.minimizerW;
    params.sec_pri_ratio = cfg.map.sec_pri_ratio;
    params.sec_pri_num = cfg.map.sec_pri_num;
    params.out_paf = cfg.map.outPAF;
    params.threads = cfg.global.threads;
    return opt::init_opts(params);
}

opt::AlignmentOptions graph_alignment_options(const AppConfig& cfg)
{
    opt::AlignmentOptions options;
    options.chain = cfg.map.chainOpts;
    options.anchor = cfg.map.anchorOpts;
    options.extend = cfg.map.extendOpts;
    options.align = cfg.map.alignOpts;
    options.use_wfa = cfg.map.use_wfa;
    return options;
}

// Select the values used by one collapse round.
CollapseOpts collapse_options(
    const AppConfig& cfg,
    double min_jaccard,
    int min_eq,
    double min_match_ratio,
    double min_ali_ratio
) {
    CollapseOpts options = cfg.collapse;
    options.min_eq = min_eq;
    options.min_match_ratio = min_match_ratio;
    options.min_ali_ratio = min_ali_ratio;
    options.min_jaccard = min_jaccard;
    options.homologous_k = cfg.bubble.homo_k;
    options.homologous_window = cfg.bubble.homo_w;
    return options;
}

void run_stat(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_stat(argc, argv);
    GfaGraph G;
    G.load_from_GFA(cfg.stat.gfaFiles);
}

void run_seq(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_seq(argc, argv);
    GfaSeq G;
    G.load_from_GFA(cfg.seq.gfaFiles);
    std::vector<std::string> paths, seqs;
    G.extract_from_file(cfg.seq.pathFile, paths, seqs);
    G.save_to_file(cfg.seq.outFile, paths, seqs);
}

void run_gfa2fa(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_gfa2fa(argc, argv);
    Gfa2fa G(cfg.gfa2fa);
    G.load_from_GFA(cfg.gfa2fa.gfaFiles);
    G.dump_to_file(cfg.gfa2fa.outFile);
}

void run_depth(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_depth(argc, argv);
    opt::AlignmentOptions options = alignment_options(cfg);

    // Load graph
    DepthOpts params = cfg.depth;
    params.forbid_overlap = !cfg.depth.reads.empty() || !cfg.depth.gafFile.empty();
    GfaDepth G(std::move(params));
    G.load_from_GFA(cfg.depth.gfaFiles);

    if (cfg.depth.reads.size() > 0) {
        auto name_seqs = G.getSeqVec(cfg.global.kmerLen);
        mmidx::MinimizerIndex GIndex(name_seqs, options);
        GIndex.build_mm(true);
        GIndex.print_index_stats();
        GIndex.count_depth(cfg.depth.reads);
        G.count_from_kmer(GIndex, cfg.depth.outFile);
    } else if (cfg.depth.gafFile.size() > 0) {
        G.count_from_gaf(cfg.depth.gafFile, cfg.depth.outFile);
    } else {
        G.count_from_A(cfg.depth.outFile);
    }
}

void run_bubble(int argc, char** argv, const std::string& command_line)
{
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

    GfaBubble::GfaBubbleFinder finder(G, bubble_options(cfg));

    finder.find_bubbles();
    finder.find_homologous_paths();

    if (DEBUG_ENABLED) finder.print_bubbles();
    if (DEBUG_ENABLED) finder.print_homologous_paths();

    finder.save_bubble_as_gfa(cfg.bubble.out_prefix + ".bubbles.gfa", cfg.bubble.min_len, cfg.bubble.min_num, /*write_seq=*/true, command_line);
    finder.save_bubble_as_gfa(cfg.bubble.out_prefix + ".bubbles.noseq.gfa", cfg.bubble.min_len, cfg.bubble.min_num, /*write_seq=*/false, command_line);
    if (cfg.bubble.write_vcf) {
        finder.save_bubble_as_vcf(cfg.bubble.out_prefix, cfg.bubble.paf_file, cfg.bubble.ref_file, cfg.bubble.ali_min_mapq, cfg.bubble.ali_min_len);
    }
}

void run_deoverlap(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_deoverlap(argc, argv);

    GfaDeoverlapper G(cfg.collapse);

    G.set_opts(graph_alignment_options(cfg));
    G.load_from_GFA(cfg.collapse.gfaFiles, cfg.collapse.gfaNames);
    G.print_graph_stats();

    G.deoverlap(cfg.collapse.prefix);

    G.save_to_disk(cfg.collapse.prefix + ".deoverlap.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/true, command_line);
    G.save_to_disk(cfg.collapse.prefix + ".deoverlap.noseq.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/false, command_line);
}

void run_collapse(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_collapse(argc, argv);
    const bool ctg_mode = !cfg.collapse.input.hap1Files.empty();
    std::vector<std::string> cur_gfa_files;
    std::vector<std::string> cur_gfa_names;

    if (ctg_mode) {
        CollapseOpts options = collapse_options(
            cfg,
            cfg.collapse.min_jaccards.front(),
            cfg.collapse.min_eqs.front(),
            cfg.collapse.min_match_ratios.front(),
            cfg.collapse.min_ali_ratios.front()
        );
        options.command_line = command_line;
        GfaCtgCollapser G(options);
        G.disable_repeat_mask();
        G.set_opts(graph_alignment_options(cfg));
        G.set_complex_marking(false);
        cur_gfa_files = G.collapse_samples(options);
    } else {
        cur_gfa_files = std::move(cfg.collapse.input.gfaFiles);
        cur_gfa_names = std::move(cfg.collapse.input.sampleNames);
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

        CollapseOpts round_options = collapse_options(
            cfg, min_jaccard, min_eq, min_match_ratio, min_ali_ratio
        );
        GfaCtgCollapser G(std::move(round_options));
        if (ctg_mode) G.disable_repeat_mask();

        G.set_opts(graph_alignment_options(cfg));
        G.set_complex_marking(!ctg_mode && cfg.bubble.check_complex);
        G.load_from_GFA(cur_gfa_files, cur_gfa_names);
        G.build_nodes_connectivity_index();
        G.build_vertex_topological_index();
        if (!ctg_mode && iter == 0 && cfg.bubble.check_complex) {
            G.detect_complex_regions(
                cfg.bubble.cx_branch_degree,
                cfg.bubble.cx_hub_degree,
                cfg.bubble.cx_min_nodes,
                cfg.bubble.cx_min_branches,
                /*skip_comp=*/false
            );
        }

        BubbleOpts finder_params = bubble_options(cfg);
        finder_params.keep_nested = true;
        finder_params.same_sim = same_sim;
        finder_params.same_min_len = same_min_len;
        finder_params.diff_min_src = diff_min_src;
        finder_params.diff_sim = diff_sim;
        finder_params.diff_min_len = diff_min_len;
        GfaBubble::GfaBubbleFinder finder(G, std::move(finder_params));

        finder.find_bubbles();
        const auto& bubbles = finder.get_bubbles();
        if (DEBUG_ENABLED) finder.print_bubbles();

        finder.find_homologous_paths();
        const auto& homologous_paths = finder.get_homologous_paths();
        if (DEBUG_ENABLED) finder.print_homologous_paths();

        // Mark complex
        if (!ctg_mode && cfg.bubble.check_complex) {
            std::vector<uint32_t> complex_segs = finder.collect_orientation_conflict_segments();
            G.mark_segments_complex(complex_segs);
        }

        // Repeat normalization
        const bool is_last_iter = (iter + 1 == cfg.collapse.iterations);
        if (!ctg_mode && is_last_iter && cfg.collapse.repeat_norm_len > 0) {
            G.normalize_homopolymer_bubbles(
                bubbles,
                cfg.collapse.repeat_norm_len
            );

            G.disable_repeat_mask();
        }

        const std::string iter_prefix = cfg.collapse.prefix + ".iter" + std::to_string(iter + 1);

        G.collapse_homologous_seq(bubbles, homologous_paths, iter_prefix);
        if (G.getNumPaths() > 0) G.rebuild_component_paths();

        if (ctg_mode && is_last_iter) {
            BubbleOpts repeat_bubble_options = bubble_options(cfg);
            repeat_bubble_options.keep_nested = false;
            G.normalize_repeat_paths(cfg.collapse, std::move(repeat_bubble_options));
        }

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
}

void run_gapfill(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_gapfill(argc, argv);
    GapfillOpts params = cfg.gapfill;
    params.alignment_options = graph_alignment_options(cfg);
    params.threads = static_cast<uint32_t>(cfg.global.threads);
    params.html_file = cfg.gapfill.prefix + ".gapfill.html";

    auto find_gapfill_bubbles = [&cfg](const GfaGraph& graph) {
        BubbleOpts params = bubble_options(cfg);
        params.max_depth = cfg.gapfill.max_depth;
        params.max_paths = cfg.gapfill.max_paths;
        params.DFS_guard = cfg.gapfill.DFS_guard;
        params.path_sim = cfg.gapfill.path_sim;
        params.stall_round_limit = cfg.gapfill.stall_round_limit;
        params.keep_nested = true;

        const bool debug = DEBUG_ENABLED;
        DEBUG_ENABLED = false;
        auto finder = std::make_unique<GfaBubble::GfaBubbleFinder>(graph, std::move(params));
        finder->find_bubbles();
        DEBUG_ENABLED = debug;
        return finder;
    };
    auto save_primary_vcfs = [&cfg](
        const GfaBubble::GfaBubbleFinder& finder,
        const GfaGapfill::PrimaryPaths& paths
    ) {
        for (uint8_t hap = 0; hap < 2; ++hap) {
            finder.save_bubble_as_vcf(
                cfg.gapfill.prefix + ".primary.hap" + std::to_string(hap + 1) + ".gapfill",
                "", "", cfg.bubble.ali_min_mapq, cfg.bubble.ali_min_len,
                &paths[hap]
            );
        }
    };
    std::vector<std::string> intermediate_gfas;
    auto protect_gapfill_input = [&cfg](const std::string& file) {
        std::error_code error;
        if (!std::filesystem::equivalent(file, cfg.gapfill.gfa_file, error)) return;
        error_stream() << "Refusing to overwrite the input GFA with the intermediate file '" << file << "'\n";
        std::exit(1);
    };
    auto cleanup_gapfill_intermediates = [&intermediate_gfas]() {
        for (const std::string& file : intermediate_gfas) {
            std::error_code error;
            std::filesystem::remove(file, error);
            if (error) warning_stream() << "  ! Cannot remove intermediate GFA '" << file << "': " << error.message() << '\n';
        }
    };

    std::string input = cfg.gapfill.gfa_file;
    std::vector<GfaGapfill::MisassemblyContig> misassemblies;
    std::unordered_set<std::string> relocated;
    bool gapfill_complete = false;
    {
        GfaGapfill preliminary(params);
        preliminary.load_from_GFA({input});
        preliminary.build_nodes_connectivity_index(false);
        preliminary.build_vertex_topological_index(false);
        misassemblies = preliminary.prepare_misassemblies();

        if (misassemblies.empty()) {
            auto finder = find_gapfill_bubbles(preliminary);
            preliminary.gapfill(finder->get_bubbles());
            const GfaGapfill::PrimaryPaths primary_paths = preliminary.save_samples(
                cfg.gapfill.prefix, command_line
            );
            save_primary_vcfs(*finder, primary_paths);
            gapfill_complete = true;
        } else {
            const std::string ms_file = cfg.gapfill.prefix + ".gapfill.ms.gfa";
            protect_gapfill_input(ms_file);
            preliminary.save_to_disk(ms_file, true, false, true, command_line);
            intermediate_gfas.push_back(ms_file);
        }
    }

    if (!misassemblies.empty()) {
        std::unordered_map<std::string, uint32_t> source_components;
        source_components.reserve(misassemblies.size());
        for (const GfaGapfill::MisassemblyContig& ms : misassemblies) {
            source_components.emplace(ms.name, ms.source_component);
        }

        CollapseOpts collapse_params = collapse_options(
            cfg,
            cfg.collapse.min_jaccards.front(),
            cfg.collapse.min_eqs.front(),
            cfg.gapfill.min_match,
            cfg.gapfill.min_ali_ratio
        );
        collapse_params.min_mapq = cfg.gapfill.min_mapq;
        GfaCtgCollapser collapser(std::move(collapse_params));
        collapser.set_opts(graph_alignment_options(cfg));
        collapser.set_complex_marking(false);
        collapser.load_from_GFA({cfg.gapfill.prefix + ".gapfill.ms.gfa"});
        collapser.build_nodes_connectivity_index(false);
        collapser.build_vertex_topological_index(false);
        collapser.rebuild_component_paths();
        const std::string collapse_prefix = cfg.gapfill.prefix + ".gapfill.recollapse";
        protect_gapfill_input(collapse_prefix + ".gfa");
        relocated = collapser.collapse_misassemblies(source_components, collapse_prefix);
        input = collapse_prefix + ".gfa";
        collapser.save_to_disk(input, true, false, true, command_line);
        intermediate_gfas.push_back(input);
    }

    if (!gapfill_complete) {
        GfaGapfill G(params);
        G.set_misassemblies(misassemblies, relocated);
        G.load_from_GFA({input});
        G.build_nodes_connectivity_index(false);
        G.build_vertex_topological_index(false);
        auto finder = find_gapfill_bubbles(G);
        G.gapfill(finder->get_bubbles());
        const GfaGapfill::PrimaryPaths primary_paths = G.save_samples(
            cfg.gapfill.prefix, command_line
        );
        save_primary_vcfs(*finder, primary_paths);
    }
    cleanup_gapfill_intermediates();
}

void run_file2map(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_file2map(argc, argv);
    mapconv::file_to_map_auto(cfg.file2map);
}

void run_liftover(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_liftover(argc, argv);
    cfg.liftover.cm_max_hops = cfg.coordmap.max_hops;
    cfg.liftover.cm_max_fanout = cfg.coordmap.max_fanout;
    cfg.liftover.cm_min_len = cfg.coordmap.min_len;
    cfg.liftover.cm_min_frac = cfg.coordmap.min_frac;
    cfg.liftover.cm_max_total_hits = cfg.coordmap.max_total_hits;
    cfg.liftover.threads = cfg.global.threads;
    std::vector<std::string> blocks = liftover::liftover(cfg.liftover);
    liftover::save_liftover_results(cfg.liftover.outFile, blocks);
}

void run_mapq_boost(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_mapq_boost(argc, argv);
    coordmap::CoordMap coormap_idx;
    coormap_idx.load(cfg.homq.mapFiles);

    // Restriction-site index (optional)
    rsite::Index rsite_idx;
    if (!cfg.ressit.genome_file.empty() || cfg.ressit.enzymes.size() > 0) {
        cfg.ressit.threads = cfg.global.threads;
        rsite_idx.build(cfg.ressit);
    }

    cfg.homq.cm_max_hops = cfg.coordmap.max_hops;
    cfg.homq.cm_max_fanout = cfg.coordmap.max_fanout;
    cfg.homq.cm_min_len = cfg.coordmap.min_len;
    cfg.homq.cm_min_frac = cfg.coordmap.min_frac;
    cfg.homq.cm_max_total_hits = cfg.coordmap.max_total_hits;
    cfg.homq.threads = cfg.global.threads;
    cfg.homq.io_threads = cfg.global.IOthreads;
    mapqboost::MapqBooster booster(coormap_idx, rsite_idx, cfg.homq);

    booster.run(cfg.homq.in_bam, cfg.homq.out_bam);
}

void run_align(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_align(argc, argv);
    opt::AlignmentOptions options = alignment_options(cfg);

    // Preset
    if      (cfg.map.preset == "ont")      opt::Preset::map_ont(options.chain, options.anchor, options.extend, options.align);
    else if (cfg.map.preset == "hifi")     opt::Preset::map_hifi(options.chain, options.anchor, options.extend, options.align);
    else if (cfg.map.preset == "illumina") opt::Preset::map_illumina(options.chain, options.anchor, options.extend, options.align);
    else if (cfg.map.preset == "asm5")     opt::Preset::map_asm_5(options.chain, options.anchor, options.extend, options.align);
    else                                   opt::Preset::map_other(options.chain, options.anchor, options.extend, options.align);

    // Build index
    GfaGraph G;
    G.load_from_GFA(cfg.map.gfaFiles);
    G.build_nodes_connectivity_index();
    auto name_seqs = G.getSeqVec();
    mmidx::MinimizerIndex GIndex(name_seqs, options);
    GIndex.build_mm();
    GIndex.print_index_stats();

    // Align
    aligner::Alignmenter A(GIndex, name_seqs, options);
    A.align(cfg.map.reads, cfg.map.outFile, command_line);
}

void run_ressit(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_res_cut(argc, argv);

    rsite::Index idx;
    cfg.ressit.threads = cfg.global.threads;
    idx.build(cfg.ressit);

    if (!cfg.ressit.out_bed.empty()) {
        idx.save_bed(cfg.ressit.out_bed);
    }
}

void run_split(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_split(argc, argv);

    GfaSplitter G;
    G.load_from_GFA(cfg.split.gfaFiles, cfg.split.gfaNames);

    G.split_by_components(
        cfg.split.prefix,
        /*write_paths=*/true,
        /*write_align=*/true,
        /*skip_comp=*/true
    );
}

void run_augment(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_augment(argc, argv);
    GfaAugmenter G(cfg.collapse);
    G.set_opts(graph_alignment_options(cfg));
    G.load_from_GFA(cfg.augment.gfa_files, cfg.augment.gfa_names);
    G.augment(cfg.augment.vcf_files, cfg.augment.prefix, cfg.augment.tag);
    G.save_to_disk(cfg.augment.prefix + ".augment.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/true, command_line);
    G.save_to_disk(cfg.augment.prefix + ".augment.noseq.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/false, command_line);
}

void run_clean(int argc, char** argv, const std::string& command_line)
{
    AppConfig cfg = main_clean(argc, argv);
    GfaCleaner G(cfg.clean);
    G.load_from_GFA({cfg.clean.utg_file});
    G.clean(cfg.clean.ctg_files);
    G.save_read_placements(cfg.clean.prefix + ".clean.reads.tsv");
    G.save_to_disk(cfg.clean.prefix + ".clean.gfa", /*write_paths=*/true, /*write_align=*/true, /*write_seq=*/true, command_line);
    G.save_to_disk(cfg.clean.prefix + ".clean.noseq.gfa", /*write_paths=*/true, /*write_align=*/true, /*write_seq=*/false, command_line);
}

}  // namespace

int run_command(
    const std::string& subcommand,
    int argc,
    char** argv,
    const std::string& command_line
) {
    if (subcommand == "stat") {
        run_stat(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "seq") {
        run_seq(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "gfa2fa") {
        run_gfa2fa(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "depth") {
        run_depth(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "bubble") {
        run_bubble(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "deoverlap") {
        run_deoverlap(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "collapse") {
        run_collapse(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "gapfill") {
        run_gapfill(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "file2map") {
        run_file2map(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "liftover") {
        run_liftover(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "mapq_boost") {
        run_mapq_boost(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "align") {
        run_align(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "ressit") {
        run_ressit(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "split") {
        run_split(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "augment") {
        run_augment(argc, argv, command_line);
        return 0;
    }
    if (subcommand == "clean") {
        run_clean(argc, argv, command_line);
        return 0;
    }

    error_stream() << "Unknown subcommand: " << subcommand << "\n";
    help(argv);
    return 1;
}
