/*------------------------------------------------------------------*
 * Implementation of Alignmenter
 *------------------------------------------------------------------*/
#include <cassert>
#include <iostream>
#include <zlib.h>

#include "../include/aligner.hpp"
#include "../include/kseq.h"

KSEQ_INIT(gzFile, gzread)

namespace aligner {

/*====================== 2. Alignmenter ============================*/
Alignmenter::Alignmenter(
    const mmidx::MinimizerIndex& GIndex,
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs,
    const opt::ExtendOpts& extendOpts,
    const opt::AlignOpts& alignParams
) : GIndex_(GIndex), names_(names), seqs_(seqs), extendOpts_(extendOpts), alignParams_(alignParams) {}

/*------------------------------------------------------------------*
 * flush SAM buffer (under out_mtx_)
 *------------------------------------------------------------------*/
static void flush_buffer(
    const opt::AlignOpts& alignParams,
    SAVE& saver,
    std::vector<Alignmenter::outBatch>& buf,
    const std::vector<std::string>& names,
    const std::vector<std::string_view>& seqs
) {
    for (auto& b : buf)
        if (!alignParams.outPAF) { saver.save(seedExtend::format_sam_record(b.aligns, b.seq, b.qual, alignParams.sec_pri_ratio)); }
        else { saver.save(seedExtend::format_paf_record(b.aligns, names, seqs, b.seq)); }
    buf.clear();
}

/*====================== 3. Main driver ============================*/
void Alignmenter::align(
    const std::vector<std::string>& reads, 
    const std::string& outputFile, 
    const std::string& cmd
) {
    // Log input files
    if (!reads.empty()) {
        std::string all = std::accumulate(
            std::next(reads.begin()), reads.end(), reads[0],
            [](const std::string &a, const std::string &b){
                return a + ' ' + b;
            }
        );
        log_stream() << "Aligning reads from: " << all << "\n";
    }

    // Open output file and write header if needed
    SAVE saver(outputFile);
    if (!alignParams_.outPAF) {
        saver.save(seedExtend::format_sam_header(names_, seqs_, cmd));
    }

    // ThreadPool with bounded internal queue.
    // If alignParams_.queue_cap == 0, ThreadPool defaults to 2 * threads.
    ThreadPool pool(alignParams_.threads, alignParams_.queue_cap);

    // We collect futures to ensure all tasks are finished before final flush.
    std::vector<std::future<void>> futs;
    futs.reserve(1024);

    // Producer: read sequencing reads and submit one task per record.
    // submit() will block when the pool's internal queue is full, providing back-pressure.
    auto tracker = ProgressTracker::Every(30000);  // print every 30,000 reads
    for (const auto& read : reads) {
        gzFile fq = gzopen(read.c_str(), "r");
        if (!fq) {
            error_stream() << read << ": No such file or directory\n";
            std::exit(1);
        }
        kseq_t* ks = kseq_init(fq);

        while (kseq_read(ks) >= 0) {
            tracker.hit();  // update progress

            // Materialize record strings and move them into the task.
            std::string name(ks->name.s, ks->name.l);
            std::string seq (ks->seq.s , ks->seq.l);
            std::string qual;
            if (ks->qual.l > 0) qual.assign(ks->qual.s, ks->qual.l);
            else                qual.assign(ks->seq.l, 'I');

            // Submit one task per read. This lambda owns its strings (moved in).
            futs.emplace_back(pool.submit(
                [this, &saver, name = std::move(name), seq = std::move(seq), qual = std::move(qual)]() mutable {
                    // 1) Run the per-read pipeline
                    auto aligns = produce_read(name, seq);

                    // 2) Buffer output (SAM/PAF) under out_mtx_ and periodic flush
                    {
                        std::lock_guard<std::mutex> lk(out_mtx_);
                        out_buf_.push_back({std::move(aligns), std::move(seq), std::move(qual)});
                        ++processed_reads_;
                        if (processed_reads_ % alignParams_.buffer_size == 0) {
                            flush_buffer(alignParams_, saver, out_buf_, names_, seqs_);
                        }
                    }
                }
            ));
        }

        kseq_destroy(ks);
        gzclose(fq);
    }
    tracker.finish();  // finish progress

    // Wait for all tasks to finish
    for (auto& f : futs) f.get();
    pool.stop();

    // Final flush
    {
        std::lock_guard<std::mutex> lk(out_mtx_);
        flush_buffer(alignParams_, saver, out_buf_, names_, seqs_);
    }
}
// void Alignmenter::align(
//     const std::vector<std::string>& reads, 
//     const std::string& outputFile, 
//     const std::string& cmd
// ) {
//     if (!reads.empty()) {
//         std::string all = std::accumulate(
//             std::next(reads.begin()), reads.end(), reads[0],
//             [](const std::string &a, const std::string &b){
//                 return a + ' ' + b;
//             }
//         );
//         log_stream() << "Aligning reads from: " << all << "\n";
//     }

//     /* ---------- open output file ---------- */
//     SAVE saver(outputFile);
//     if (!alignParams_.outPAF) saver.save(seedExtend::format_sam_header(names_, seqs_, cmd));

//     /* ---------- ThreadPool ---------- */
//     ThreadPool pool(alignParams_.threads);
//     std::vector<std::future<void>> futs;
//     for (uint16_t t = 0; t < alignParams_.threads; ++t)
//         futs.emplace_back(pool.submit([this, &saver]{ this->worker_loop(saver); }));

//     /* ---------- producer: read sequencing reads and push jobs ---------- */
//     auto tracker = ProgressTracker::Every(30000);  // print every 30,000 reads
//     for (const auto& read : reads) {
//         gzFile fq = gzopen(read.c_str(), "r");
//         if (!fq) {
//             error_stream() << "'" << read << "': No such file or directory.\n";
//             std::exit(1);
//         }
//         kseq_t* ks = kseq_init(fq);
//         while (kseq_read(ks) >= 0) {
//             tracker.hit();  // update progress
//             auto rec = std::make_unique<ReadRec>();
//             rec->name.assign(ks->name.s, ks->name.l);
//             rec->seq.assign(ks->seq.s , ks->seq.l);
//             if (ks->qual.l > 0) rec->qual.assign(ks->qual.s, ks->qual.l);
//             else rec->qual.assign(ks->seq.l, 'I');

//             std::unique_lock<std::mutex> lk(read_mtx_);
//             read_cv_prod_.wait(lk, [&]{ return read_q_.size() < alignParams_.queue_cap; });
//             read_q_.push(std::move(rec));
//             lk.unlock();
//             read_cv_cons_.notify_one();
//         }
//         kseq_destroy(ks); gzclose(fq);
//     }
//     tracker.finish();  // finish progress

//     // notify consumer threads that producer is done
//     producer_done_.store(true, std::memory_order_release);
//     read_cv_cons_.notify_all();

//     /* ---------- wait workers ---------- */
//     for (auto& f : futs) f.get();
//     pool.stop();

//     /* ---------- final flush ---------- */
//     {
//         std::lock_guard<std::mutex> lk(out_mtx_);
//         flush_buffer(alignParams_, saver, out_buf_, names_, seqs_);
//     }
// }

/*====================== 4. Worker loop ============================*/
void Alignmenter::worker_loop(SAVE& saver)
{
    for (;;) {
        std::unique_ptr<ReadRec> rec;
        /* ---- blocking pop ---- */
        {
            std::unique_lock<std::mutex> lk(read_mtx_);
            read_cv_cons_.wait(lk, [&]{
                return !read_q_.empty() || producer_done_.load();
            });
            if (read_q_.empty()) break;                // finished
            rec = std::move(read_q_.front());
            read_q_.pop();
            lk.unlock();
            read_cv_prod_.notify_one();
        }

        /* ---- pipeline ---- */
        auto aligns = produce_read(rec->name, rec->seq);

        /* ---- buffer ---- */
        {
            std::lock_guard<std::mutex> lk(out_mtx_);
            out_buf_.push_back({std::move(aligns), std::move(rec->seq), std::move(rec->qual)});
            ++processed_reads_;
            if (processed_reads_ % alignParams_.buffer_size == 0)
                flush_buffer(alignParams_, saver, out_buf_, names_, seqs_);
        }
    }
}

/*====================== 5. per-read pipeline ======================*/
std::vector<seedExtend::FragAlign> Alignmenter::produce_read(
    std::string_view read_name,
    std::string_view read_seq, 
    bool keep_same_strand_only
) {
    const uint32_t read_len = read_seq.size();
    
    /* 1. seeds */
    std::vector<mmidx::MM128> seeds = GIndex_.collect_seeds(read_seq, keep_same_strand_only);
    GIndex_.sort_seeds_by_ref(seeds);  // sort by reference coordinate
    if (DEBUG_ENABLED) GIndex_.print_seeds(read_name, read_len, seeds);

    /* 2. chain */
    std::vector<mmidx::MM128> b;
    std::vector<mmidx::Sc_Len> u;
    GIndex_.chain_dp(seeds, b, u);
    if (DEBUG_ENABLED) GIndex_.print_chains(read_name, read_len, b, u);

    /* 3.1 anchor */
    std::vector<mmidx::AnchorChain> anchors = GIndex_.build_anchor_from_chains(b, u, read_len);

    /* 3.2 merge anchors */
    anchors = GIndex_.merge_anchor_chains(anchors);

    /* 3.3 group anchors */
    std::vector<std::vector<uint32_t>> groups = GIndex_.group_anchor_chains_by_q_overlap(anchors);

    /* 3.4 filter anchors */
    std::vector<std::vector<uint32_t>> filtered_groups = GIndex_.filter_anchor_chains(groups, alignParams_.sec_pri_num);
    if (DEBUG_ENABLED) GIndex_.print_anchors(read_name, read_len, anchors, filtered_groups, "Filtered");

    /* 4. for each group: extend + score + filter */
    std::vector<seedExtend::FragAlign> fragAligns;
    for (auto& g : filtered_groups) {
        std::vector<mmidx::AnchorChain> sub_anchors;
        for (auto idx : g) sub_anchors.push_back(anchors[idx]);

        /* 4.1 extension */
        std::vector<seedExtend::FragAlign> sub_fragAligns = seedExtend::extend_chain_wfa(
            sub_anchors, names_, seqs_, read_name, read_seq, extendOpts_
        );
        /* 4.2 compute alignment scores */
        seedExtend::cal_align_scores(sub_fragAligns, alignParams_.sec_pri_ratio);
        /* 4.3 filter aligns */
        seedExtend::filter_aligns(sub_fragAligns, extendOpts_, alignParams_.sec_pri_num);

        fragAligns.insert(fragAligns.end(), sub_fragAligns.begin(), sub_fragAligns.end());
    }

    return fragAligns;
}

} // namespace aligner