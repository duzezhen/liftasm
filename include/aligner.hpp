#pragma once
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <numeric>

#include "gfa_parser.hpp"
#include "mmidx.hpp"
#include "seed_extend.hpp"
#include "ThreadPool.hpp"
#include "save.hpp"
#include "options.hpp"
#include "progress_tracker.hpp"
#include "logger.hpp"

namespace aligner {

/*------------------------------ 2. aligner --------------------------*/
struct ReadRec { std::string name; std::string seq; std::string qual; };

class Alignmenter {
public:
    Alignmenter(
        const mmidx::MinimizerIndex& GIndex,
        const std::vector<std::string>& names,
        const std::vector<std::string_view>& seqs,
        const opt::ExtendOpts& extendOpts,
        const opt::AlignOpts& alignParams
    );

    void align(
        const std::vector<std::string>& reads,
        const std::string& outputFile, 
        const std::string& cmd
    );

    /* single-read pipeline */
    // false = all seeds; true = only seeds on the same strand as the read
    std::vector<seedExtend::FragAlign> produce_read(std::string_view read_name, std::string_view read_seq, bool keep_same_strand_only = false);

    /* batch for buffered output */
    struct outBatch {
        std::vector<seedExtend::FragAlign> aligns;
        std::string                        seq;
        std::string                        qual;
    };

private:
    const opt::AlignOpts& alignParams_;
    const opt::ExtendOpts& extendOpts_;

    const mmidx::MinimizerIndex&   GIndex_;
    // reference sequences
    const std::vector<std::string>& names_;
    const std::vector<std::string_view>& seqs_;

    /* one worker thread */
    void worker_loop(SAVE& saver);

    /* ---------- bounded queue (producer / consumer) ---------- */
    std::queue<std::unique_ptr<ReadRec>> read_q_;
    std::mutex              read_mtx_;
    std::condition_variable read_cv_prod_;
    std::condition_variable read_cv_cons_;
    std::atomic<bool>       producer_done_{false};

    /* ---------- SAM buffer ---------- */
    std::vector<outBatch>   out_buf_;
    std::mutex              out_mtx_;
    size_t                  processed_reads_ = 0;

};

} // namespace aligner