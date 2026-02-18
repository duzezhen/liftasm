#pragma once
#include <cstddef>
#include <string>
#include <iomanip>
#include "logger.hpp"

/*
Usage:
------
1. Percentage mode (known total number of tasks)
   --------------------------------------------
   ProgressTracker tracker(total_tasks, 10); // print every 10%
   for (size_t i = 0; i < total_tasks; ++i) {
       // ... do work ...
       tracker.hit(); // or tracker.update(i+1);
   }
   tracker.finish(); // ensure the final 100% is printed

   Example output:
   [I::...::check_progress....] [Progress]  10% ( 50/492) done
   ...
   [I::...::check_progress....] [Progress] 100% (492/492) done

2. Step mode (unknown total number of tasks)
   -----------------------------------------
   auto tracker = ProgressTracker::Every(1000); // print every 1000 items
   size_t count = 0;
   while (has_more()) {
       // ... do work ...
       tracker.hit(); // or tracker.update(++count);
   }
   tracker.finish(); // print final processed count

   Example output:
   [I::...::check_progress....] [Progress] processed 1000 items
   [I::...::check_progress....] [Progress] processed 2000 items
   ...
   [I::...::check_progress....] [Progress] processed 5234 items (done)
*/

class ProgressTracker {
public:
    // ===== Percentage mode (known total) =====
    ProgressTracker(std::size_t total, std::size_t step_percent = 10)
        : mode_(Mode::Percent),
          total_(total),
          step_percent_(step_percent ? step_percent : 1),
          step_count_(0),
          next_step_(step_percent ? step_percent : 1),
          next_count_(0),
          processed_(0),
          finished_printed_(false),
          width_count_(total ? int(std::to_string(total).size()) : 1) {}

    // ===== Step mode (unknown total): print once every step_every items =====
    static ProgressTracker Every(std::size_t step_every) {
        ProgressTracker t;
        t.mode_ = Mode::Step;
        t.total_ = 0;
        t.step_percent_ = 0;
        t.step_count_ = (step_every ? step_every : 1);
        t.next_count_ = t.step_count_;
        t.next_step_ = 0;
        t.processed_ = 0;
        t.finished_printed_ = false;
        t.width_count_ = 1;
        return t;
    }

    void hit() {
        ++processed_;
        check_progress_();
    }

    void update(std::size_t done) {
        processed_ = done;
        check_progress_();
    }

    void finish() {
        if (mode_ == Mode::Percent) {  // Mode::Percent
            if (!finished_printed_ && total_ > 0) {
                log_stream()
                    << " [Progress] " << std::setw(3) << 100 << "% ("
                    << std::setw(width_count_) << processed_ << "/" << total_
                    << ") done\n";
                finished_printed_ = true;
            }
        } else {  // Mode::Step
            log_stream() << " [Progress] processed " << processed_ << " items (done)\n";
            finished_printed_ = true;
        }
    }

private:
    enum class Mode { Percent, Step };

    ProgressTracker() = default;

    void check_progress_() {
        if (mode_ == Mode::Percent) {  // Mode::Percent
            if (total_ == 0) return;

            double percent = (100.0 * processed_) / total_;

            if (percent >= 100.0 && !finished_printed_) {
                log_stream()
                    << " [Progress] " << std::setw(3) << 100 << "% ("
                    << std::setw(width_count_) << processed_ << "/" << total_
                    << ") done\n";
                finished_printed_ = true;
                return;
            }

            if (!finished_printed_ && percent >= next_step_) {
                log_stream()
                    << " [Progress] " << std::setw(3) << static_cast<int>(percent) << "% ("
                    << std::setw(width_count_) << processed_ << "/" << total_
                    << ") done\n";
                next_step_ += step_percent_;
            }
        } else {  // Mode::Step
            if (step_count_ == 0) return;
            if (processed_ >= next_count_ && !finished_printed_) {
                log_stream() << " [Progress] processed " << processed_ << " items\n";
                next_count_ += step_count_;
            }
        }
    }

    // === Mode & Status ===
    Mode mode_{Mode::Percent};

    std::size_t total_{0};         // Percentage mode
    std::size_t step_percent_{0};  // Percentage step size
    std::size_t step_count_{0};    // Step mode step size

    double      next_step_{0};     // Next percentage threshold
    std::size_t next_count_{0};    // Next count threshold

    std::size_t processed_{0};
    bool        finished_printed_{false};
    int         width_count_{1};
    std::string task_name_;        // Reserved, not used in logic
};