#include "../include/save.hpp"

/*------------------------------------------------------------*/
/*                       constructor                         */
/*------------------------------------------------------------*/
SAVE::SAVE(const std::string& outFileName, size_t cacheSize)
    : outputFileName_(outFileName), cache_size_(cacheSize) {
    is_gzip_ = (outputFileName_.find(".gz") != std::string::npos ||
                outputFileName_.find(".GZ") != std::string::npos);

    if (is_gzip_) {
        gzFile fp = gzopen(outputFileName_.c_str(), "wb");
        if (!fp) {
            error_stream() << outputFileName_ << ": No such file or directory\n";
            std::exit(1);
        }
        gzfpO_.reset(fp);
    } else if (!outputFileName_.empty()) {
        fpO.open(outputFileName_, std::ios::out);
        if (!fpO) {
            error_stream() << outputFileName_ << ": No such file or directory\n";
            std::exit(1);
        }
    }

    buffer_.reserve(cache_size_);
}

/*------------------------------------------------------------*/
/*                         destructor                         */
/*------------------------------------------------------------*/
SAVE::~SAVE() {
    flush();                                // write remaining bytes
}

/*------------------------------------------------------------*/
/*                          flush                             */
/*------------------------------------------------------------*/
void SAVE::flush() {
    if (buffer_.empty()) return;

    if (is_gzip_) {
        gzwrite(gzfpO_.get(), buffer_.c_str(),
                static_cast<unsigned int>(buffer_.size()));
    } else if (!outputFileName_.empty()) {
        fpO.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
    } else {
        std::cout.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
    }
    buffer_.clear();
}

/*------------------------------------------------------------*/
/*                           save                             */
/*------------------------------------------------------------*/
int SAVE::save(const std::string& outTxt) {
    if (outTxt.empty()) return 0;

    buffer_.append(outTxt);

    if (buffer_.size() >= cache_size_) flush();
    return 0;
}
