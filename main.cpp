// g++ main.cpp src/*.cpp -o liftasm -std=c++20 -march=native -O3 -lhts -lwfa2cpp -lz -lbz2 -lcrypto -lssl -lpthread -fopenmp -lstdc++fs -I /home/zd233/zd233/00-software/00-project/liftasm/minimap2 /home/zd233/zd233/00-software/00-project/liftasm/minimap2/libminimap2.a
// g++ main.cpp src/*.cpp -o liftasm -std=c++20 -march=native -O3 -fopenmp -I /home/zd233/zd233/00-software/00-project/liftasm/minimap2 /home/zd233/zd233/00-software/00-project/liftasm/minimap2/libminimap2.a /home/zd233/zd233/00-software/WFA2-lib/build/libwfa2cpp.a /home/zd233/zd233/00-software/WFA2-lib/build/libwfa2.a /home/zd233/zd233/00-software/htslib-1.22.1/build/libhts.a -lz -lbz2 -lcrypto -lssl -lpthread -ldl -lstdc++fs

#include <iostream>
#include <string>
#include <iomanip>

#include "include/OptionParser.hpp"
#include "include/options.hpp"
#include "include/gfa_parser.hpp"
#include "include/gfa_depth.hpp"
#include "include/gfa_seq.hpp"
#include "include/gfa_deoverlapper.hpp"
#include "include/gfa_collapser.hpp"
#include "include/gfa_bubble.hpp"
#include "include/gfa_gfa2fa.hpp"
#include "include/file2map.hpp"
#include "include/liftover.hpp"
#include "include/coordmap.hpp"
#include "include/MapqBoost.hpp"
#include "include/aligner.hpp"
#include "include/sys.hpp"

using namespace wfa;

static inline bool is_top_help_flag(const char* s) {
    return (std::string(s) == "-h" || std::string(s) == "--help");
}
static inline bool is_top_ver_flag(const char* s) {
    return (std::string(s) == "-v" || std::string(s) == "--version");
}

int main(int argc, char** argv) {
    if (argc < 2) { help(argv); return 1; }
    if (is_top_help_flag(argv[1])) { help(argv); return 0; }
    if (is_top_ver_flag(argv[1]))  { std::cerr << program::version << "\n"; return 0; }

    // Dispatch by subcommand
    const std::string sub = argv[1];

    // timing
    double realtime0 = realtime();

    if (sub == "stat") {
        AppConfig cfg = main_stat(argc, argv);
        GfaGraph G;
        G.load_from_GFA(cfg.in.gfaFile);
        G.print_graph_stats();
    } else if (sub == "bubble") {
        AppConfig cfg = main_bubble(argc, argv);
        GfaGraph G;
        G.load_from_GFA(cfg.in.gfaFile);
        GfaBubble::GfaBubbleFinder finder(
            G, cfg.bubble.max_depth, cfg.bubble.max_paths, cfg.bubble.DFS_guard, 
            cfg.bubble.cx_depth, cfg.bubble.cx_nodes, cfg.bubble.cx_branches, cfg.bubble.cx_deg_branch, cfg.bubble.cx_deg_hub, cfg.bubble.cx_deg_cap, 
            cfg.bubble.path_diff, /*skip_comp=*/false, cfg.global.threads
        );
        finder.find_bubbles();
        const auto& bubbles = finder.get_bubbles();
        finder.save_bubble_as_gfa(cfg.out.gfaBubbles, cfg.bubble.min_len, cfg.bubble.min_num);
    } else if (sub == "seq") {
        AppConfig cfg = main_seq(argc, argv);
        GfaSeq G;
        G.load_from_GFA(cfg.in.gfaFile);
        std::vector<std::string> paths, seqs;
        G.extract_from_file(cfg.in.pathFile, paths, seqs);
        G.save_to_file(cfg.out.gfaSeq, paths, seqs);
    } else if (sub == "depth") {
        AppConfig cfg = main_depth(argc, argv);
        opt::ChainOpts  chainOpts;
        opt::AnchorOpts anchorOpts;
        opt::ExtendOpts extendOpts;
        opt::AlignOpts  alignOpts;
        init_opts(
            cfg.global.kmerLen, cfg.global.minimizerW,
            cfg.map.sec_pri_ratio, cfg.map.sec_pri_num,
            cfg.map.outPAF, cfg.global.threads,
            chainOpts, anchorOpts, extendOpts, alignOpts
        );

        // Load graph
        GfaDepth G(cfg.depth.min_mapq, cfg.depth.min_frac);
        G.load_from_GFA(cfg.in.gfaFile);
        
        if (cfg.in.reads.size() > 0) {
            auto name_seqs = G.getSeqVec(cfg.global.kmerLen);
            mmidx::MinimizerIndex GIndex(name_seqs.names, name_seqs.seqs, name_seqs.right_seqs, chainOpts, anchorOpts);
            GIndex.build_mm(true);
            GIndex.print_index_stats();
            GIndex.count_depth(cfg.in.reads);
            G.count_from_kmer(GIndex, cfg.out.gfaDepth);
        } else if (cfg.in.gafFile.size() > 0) {
            G.count_from_gaf(cfg.in.gafFile, cfg.out.gfaDepth);
        }
    } else if (sub == "deoverlap") {
        AppConfig cfg = main_deoverlap(argc, argv);
        opt::ChainOpts  chainOpts;
        opt::AnchorOpts anchorOpts;
        opt::ExtendOpts extendOpts;
        opt::AlignOpts  alignOpts;
        init_opts(
            cfg.global.kmerLen, cfg.global.minimizerW,
            cfg.map.sec_pri_ratio, cfg.map.sec_pri_num,
            cfg.map.outPAF, cfg.global.threads,
            chainOpts, anchorOpts, extendOpts, alignOpts
        );
        // Preset
        opt::Preset::map_asm_5(chainOpts, anchorOpts, extendOpts, alignOpts);

        GfaDeoverlapper G(cfg.collapse.min_eq, cfg.collapse.max_iters);
        G.set_opts(chainOpts, anchorOpts, extendOpts, alignOpts, cfg.global.use_wfa);
        G.load_from_GFA(cfg.in.gfaFile);
        G.print_graph_stats();
        G.deoverlap(cfg.out.prefix);
        G.save_to_disk(cfg.out.prefix + ".deoverlap.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/true);
        G.save_to_disk(cfg.out.prefix + ".deoverlap.noseq.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/false);
    } else if (sub == "collapse") {
        AppConfig cfg = main_collapse(argc, argv);
        opt::ChainOpts  chainOpts;
        opt::AnchorOpts anchorOpts;
        opt::ExtendOpts extendOpts;
        opt::AlignOpts  alignOpts;
        init_opts(
            cfg.global.kmerLen, cfg.global.minimizerW,
            cfg.map.sec_pri_ratio, cfg.map.sec_pri_num,
            cfg.map.outPAF, cfg.global.threads,
            chainOpts, anchorOpts, extendOpts, alignOpts
        );
        // Preset
        opt::Preset::map_asm_5(chainOpts, anchorOpts, extendOpts, alignOpts);

        GfaCollapser G(cfg.collapse.min_jaccard, cfg.collapse.min_new_frac, cfg.collapse.min_eq, cfg.collapse.max_iters);
        G.set_opts(chainOpts, anchorOpts, extendOpts, alignOpts, cfg.global.use_wfa);
        G.load_from_GFA(cfg.in.gfaFile);
        G.collapse_unitigs();
        G.save_to_disk(cfg.out.prefix + ".collapse.unitig.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/true);

        // Detect Fork
        GfaBubble::GfaBubbleFinder finder(
            G, cfg.bubble.max_depth, cfg.bubble.max_paths, cfg.bubble.DFS_guard, 
            cfg.bubble.cx_depth, cfg.bubble.cx_nodes, cfg.bubble.cx_branches, cfg.bubble.cx_deg_branch, cfg.bubble.cx_deg_hub, cfg.bubble.cx_deg_cap, 
            cfg.bubble.path_diff, /*skip_comp=*/false, cfg.global.threads
        );
        finder.find_forks();
        const auto& forks = finder.get_forks();

        finder.find_bubbles(/*filter_nonlocal=*/false);
        const auto& bubbles = finder.get_bubbles();
        G.collapse_bubbles(bubbles, forks, cfg.out.prefix, finder);
        G.collapse_unitigs(); // collapse again after bubble collapse

        // Output
        G.save_to_disk(cfg.out.prefix + ".collapse.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/true);
        G.save_to_disk(cfg.out.prefix + ".collapse.noseq.gfa", /*write_paths=*/false, /*write_align=*/false, /*write_seq=*/false);
    } else if (sub == "align") {
        AppConfig cfg = main_align(argc, argv);
        opt::ChainOpts  chainOpts;
        opt::AnchorOpts anchorOpts;
        opt::ExtendOpts extendOpts;
        opt::AlignOpts  alignOpts;
        init_opts(
            cfg.global.kmerLen, cfg.global.minimizerW,
            cfg.map.sec_pri_ratio, cfg.map.sec_pri_num,
            cfg.map.outPAF, cfg.global.threads,
            chainOpts, anchorOpts, extendOpts, alignOpts
        );

        // Preset
        if      (cfg.map.preset == "ont")      opt::Preset::map_ont(chainOpts, anchorOpts, extendOpts, alignOpts);
        else if (cfg.map.preset == "hifi")     opt::Preset::map_hifi(chainOpts, anchorOpts, extendOpts, alignOpts);
        else if (cfg.map.preset == "illumina") opt::Preset::map_illumina(chainOpts, anchorOpts, extendOpts, alignOpts);
        else if (cfg.map.preset == "asm5")     opt::Preset::map_asm_5(chainOpts, anchorOpts, extendOpts, alignOpts);
        else                                   opt::Preset::map_other(chainOpts, anchorOpts, extendOpts, alignOpts);

        // Build index
        GfaGraph G;
        G.load_from_GFA(cfg.in.gfaFile);
        auto name_seqs = G.getSeqVec();
        mmidx::MinimizerIndex GIndex(name_seqs.names, name_seqs.seqs, name_seqs.right_seqs, chainOpts, anchorOpts);
        GIndex.build_mm();
        GIndex.print_index_stats();

        // Align
        aligner::Alignmenter A(GIndex, name_seqs.names, name_seqs.seqs, extendOpts, alignOpts);
        A.align(cfg.in.reads, cfg.out.alignOut, program::cmdline(argc, argv));
    } else if (sub == "gfa2fa") {
        AppConfig cfg = main_gfa2fa(argc, argv);
        Gfa2fa G(cfg.gfa2fa.min_len_xbp, cfg.gfa2fa.extend_ybp, cfg.gfa2fa.wrap_width, cfg.gfa2fa.skip_unknown);
        G.load_from_GFA(cfg.in.gfaFile);
        G.dump_to_file(cfg.out.gfa2faOut);
    } else if (sub == "file2map") {
        AppConfig cfg = main_file2map(argc, argv);
        mapconv::file_to_map_auto(cfg.in.inputFile, cfg.out.outFile, cfg.file2map.paf_primary_only, cfg.file2map.min_len, cfg.file2map.min_mapq);
    } else if (sub == "liftover") {
        AppConfig cfg = main_liftover(argc, argv);
        liftover::RunOpts opt = liftover::set_opts(
            cfg.in.mapFile, cfg.in.bedFile, cfg.in.referenceFile, cfg.out.outFile,
            cfg.liftover.regex, cfg.liftover.min_frac, 
            cfg.liftover.flank_win, cfg.liftover.max_flank, cfg.liftover.max_gap, cfg.liftover.max_hit,
            cfg.in.pafFile, cfg.liftover.min_mapq, cfg.liftover.min_len,
            cfg.liftover.do_check, cfg.liftover.win, cfg.liftover.step, cfg.liftover.max_examples, 
            cfg.coordmap.max_hops, cfg.coordmap.max_fanout, cfg.coordmap.min_len, cfg.coordmap.min_frac, cfg.coordmap.max_total_hits, 
            cfg.global.threads
        );
        std::vector<std::string> blocks = liftover::liftover(opt);
        liftover::save_liftover_results(opt.out_file, blocks);
    } else if (sub == "mapq_boost") {
        AppConfig cfg = main_mapq_boost(argc, argv);
        coordmap::CoordMap map;
        map.load(cfg.in.mapFile);
        log_stream() << "Loaded records: " << map.num_records() << "\n";

        mapqboost::MapqBooster booster(
            map, 
            cfg.homq.batch_size, cfg.homq.mapq_low, cfg.homq.mapq_new, cfg.homq.min_frac, cfg.homq.min_equiv, cfg.homq.name_check,
            cfg.coordmap.max_hops, cfg.coordmap.max_fanout, cfg.coordmap.min_len, cfg.coordmap.min_frac, cfg.coordmap.max_total_hits,
            cfg.global.threads, cfg.global.IOthreads
        );
        booster.run(cfg.homq.in_bam, cfg.homq.out_bam);
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