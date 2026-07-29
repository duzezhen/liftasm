#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <cstdlib>

#include "gfa_name.hpp"
#include "logger.hpp"

namespace gfa_name_check {

struct RenameStats {
    uint32_t invalid_reserved_total = 0;
    uint32_t trailing_sign_total    = 0;
    uint32_t duplicate_total        = 0;
    uint32_t namespace_total        = 0;

    std::vector<std::pair<std::string, std::string>> invalid_reserved_examples;
    std::vector<std::pair<std::string, std::string>> trailing_sign_examples;
    std::vector<std::pair<std::string, std::string>> duplicate_examples;
    std::vector<std::pair<std::string, std::string>> namespace_examples;
};

class FileRenamePlan {
public:
    explicit FileRenamePlan(std::size_t max_log_examples = 5, std::string name_suffix = {})
        : max_log_examples_(max_log_examples), name_suffix_(std::move(name_suffix)) {}

    std::string register_segment_name(
        const std::string& raw_name,
        const std::unordered_map<std::string, uint64_t>& global_used_names,
        const std::string& filename = ""
    ) {
        if (!seen_raw_segment_names_.insert(raw_name).second) {
            std::string msg = "Duplicate raw segment name '" + raw_name + "'";
            if (!filename.empty()) msg += " in file '" + filename + "'";
            error_stream() << msg << "\n";
            std::exit(1);
        }

        std::string base = raw_name;

        if (!base.empty() && (base.back() == '+' || base.back() == '-')) {
            std::string before = base;
            base.back() = '_';
            ++stats_.trailing_sign_total;
            push_example_(stats_.trailing_sign_examples, before, base);
        }

        if (contains_reserved_chars_(base) && !is_valid_coord_name_(base)) {
            std::string before = base;
            sanitize_reserved_chars_(base);
            ++stats_.invalid_reserved_total;
            push_example_(stats_.invalid_reserved_examples, before, base);
        }

        if (!name_suffix_.empty()) {
            const std::string before = base;
            base = append_root_suffix_(base, name_suffix_);
            ++stats_.namespace_total;
            push_example_(stats_.namespace_examples, before, base);
        }

        std::string final_name = base;
        if (is_name_used_(final_name, global_used_names)) {
            final_name = make_unique_name_(base, global_used_names);
            ++stats_.duplicate_total;
            push_example_(stats_.duplicate_examples, raw_name, final_name);
        }

        raw_to_norm_.emplace(raw_name, final_name);
        staged_norm_names_.insert(final_name);
        return final_name;
    }

    std::string rewrite_gfa_line(const std::string& line) const {
        if (line.empty()) return line;

        const char type = line[0];
        std::stringstream ss(line);

        if (type == 'S') {
            char c;
            std::string name, seq, rest;
            ss >> c >> name >> seq;
            std::getline(ss, rest);

            name = remap_seg_name_(name);

            std::ostringstream out;
            out << c << '\t' << name << '\t' << seq << rest;
            return out.str();
        }

        if (type == 'L') {
            char c, v_ori, w_ori;
            std::string v_name, w_name, ovlp, rest;
            ss >> c >> v_name >> v_ori >> w_name >> w_ori >> ovlp;
            std::getline(ss, rest);

            v_name = remap_seg_name_(v_name);
            w_name = remap_seg_name_(w_name);

            std::ostringstream out;
            out << c << '\t' << v_name << '\t' << v_ori
                << '\t' << w_name << '\t' << w_ori
                << '\t' << ovlp << rest;
            return out.str();
        }

        if (type == 'P') {
            char c;
            std::string pname, segs, cigar;
            ss >> c >> pname >> segs >> cigar;

            std::ostringstream out;
            out << c << '\t' << pname << '\t' << remap_p_segments_(segs) << '\t' << cigar;
            return out.str();
        }

        if (type == 'A') {
            char c, strand;
            std::string unitig, rname, rest;
            uint32_t pos = 0, rstart = 0, rend = 0;
            ss >> c >> unitig >> pos >> strand >> rname >> rstart >> rend;
            std::getline(ss, rest);

            unitig = remap_seg_name_(unitig);

            std::ostringstream out;
            out << c << '\t' << unitig << '\t' << pos << '\t' << strand
                << '\t' << rname << '\t' << rstart << '\t' << rend << rest;
            return out.str();
        }

        return line;
    }

    void print_log(const std::string& filename) const {
        if (stats_.invalid_reserved_total == 0 &&
            stats_.trailing_sign_total == 0 &&
            stats_.duplicate_total == 0 &&
            stats_.namespace_total == 0) {
            return;
        }

        log_stream() << "Renaming segments ...\n";
        print_one_log_("Fixed segment names: replaced invalid characters (':', ';', '-') with '_'", stats_.invalid_reserved_total, stats_.invalid_reserved_examples);
        print_one_log_("Fixed segment names: replaced trailing '+' or '-' with '_'", stats_.trailing_sign_total, stats_.trailing_sign_examples);
        print_one_log_("Fixed segment names: renamed duplicates by adding numeric suffix", stats_.duplicate_total, stats_.duplicate_examples);
        print_one_log_("Namespaced all segment names", stats_.namespace_total, stats_.namespace_examples);
        std::cerr << std::endl;
    }

private:
    std::size_t max_log_examples_{5};
    std::string name_suffix_;
    std::unordered_map<std::string, std::string> raw_to_norm_;
    std::unordered_set<std::string> staged_norm_names_;
    std::unordered_set<std::string> seen_raw_segment_names_;
    RenameStats stats_;

private:
    static bool contains_reserved_chars_(const std::string& s) {
        return s.find_first_of(":;-") != std::string::npos;
    }

    static void sanitize_reserved_chars_(std::string& s) {
        for (char& c : s) {
            if (c == ':' || c == ';' || c == '-') c = '_';
        }
        if (!s.empty() && (s.back() == '+' || s.back() == '-')) {
            s.back() = '_';
        }
    }

    static bool is_valid_coord_name_(const std::string& name) {
        if (name.empty()) return false;
        if (name.back() == '+' || name.back() == '-') return false;
        if (name.find('+') != std::string::npos) return false;

        gfaName gname(';');
        auto pieces = gname.parse_composite_with_dir(name);
        if (pieces.empty()) return false;

        const std::size_t token_cnt = 1 + static_cast<std::size_t>(std::count(name.begin(), name.end(), ';'));
        if (pieces.size() != token_cnt) return false;

        for (const auto& p : pieces) {
            if (p.root.empty()) return false;
            if (p.hi <= p.lo) return false;
            if (p.root.find_first_of(":;-+") != std::string::npos) return false;
            if (!p.root.empty() && (p.root.back() == '+' || p.root.back() == '-')) return false;
        }

        return true;
    }

    static std::string append_root_suffix_(const std::string& name, const std::string& suffix) {
        if (name.find(':') == std::string::npos) return name + suffix;

        gfaName gname(';');
        const auto pieces = gname.parse_composite_with_dir(name);
        const std::size_t token_count = 1 + static_cast<std::size_t>(std::count(name.begin(), name.end(), ';'));
        if (pieces.size() != token_count) return name + suffix;

        std::ostringstream out;
        for (size_t i = 0; i < pieces.size(); ++i) {
            if (i) out << ';';
            const auto& piece = pieces[i];
            out << piece.root << suffix << ':';
            if (piece.rev) out << piece.hi << '-' << piece.lo;
            else out << piece.lo << '-' << piece.hi;
        }
        return out.str();
    }

    bool is_name_used_(
        const std::string& s,
        const std::unordered_map<std::string, uint64_t>& global_used_names
    ) const {
        return global_used_names.find(s) != global_used_names.end() ||
               staged_norm_names_.find(s) != staged_norm_names_.end();
    }

    std::string make_unique_name_(
        const std::string& base,
        const std::unordered_map<std::string, uint64_t>& global_used_names
    ) const {
        if (!is_name_used_(base, global_used_names)) return base;

        for (uint32_t k = 1; ; ++k) {
            std::string cand = base + "_" + std::to_string(k);
            if (!is_name_used_(cand, global_used_names)) return cand;
        }
    }

    std::string remap_seg_name_(const std::string& name) const {
        auto it = raw_to_norm_.find(name);
        return (it == raw_to_norm_.end()) ? name : it->second;
    }

    std::string remap_p_segments_(const std::string& segs) const {
        if (segs.empty()) return segs;

        //   >seg1<seg2>seg3
        // or
        //   <u1tg034912l:27463-28103>u1tg034917l:256-5433>u1tg034917l:5433-5436
        if (segs.find('>') != std::string::npos || segs.find('<') != std::string::npos) {
            std::string out;
            out.reserve(segs.size() + 16);

            const std::size_t n = segs.size();
            std::size_t i = 0;

            while (i < n) {
                const char c = segs[i];

                if (c == '>' || c == '<') {
                    out.push_back(c);
                    ++i;

                    const std::size_t beg = i;
                    while (i < n && segs[i] != '>' && segs[i] != '<') ++i;

                    const std::string seg_name = segs.substr(beg, i - beg);
                    out += remap_seg_name_(seg_name);
                } else {
                    // Preserve any unexpected characters verbatim.
                    out.push_back(c);
                    ++i;
                }
            }

            return out;
        }

        //   seg1+,seg2-,seg3+
        if (segs.find(',') != std::string::npos) {
            std::stringstream in(segs);
            std::string tok;
            std::string out;
            bool first = true;

            while (std::getline(in, tok, ',')) {
                if (!first) out.push_back(',');
                first = false;

                if (tok.empty()) continue;

                std::size_t l = 0, r = tok.size();
                while (l < r && std::isspace(static_cast<unsigned char>(tok[l]))) ++l;
                while (r > l && std::isspace(static_cast<unsigned char>(tok[r - 1]))) --r;

                const std::string lead  = tok.substr(0, l);
                const std::string core  = tok.substr(l, r - l);
                const std::string trail = tok.substr(r);

                if (core.size() >= 2 && (core.back() == '+' || core.back() == '-')) {
                    const char ori = core.back();
                    const std::string seg_name = core.substr(0, core.size() - 1);
                    out += lead;
                    out += remap_seg_name_(seg_name);
                    out.push_back(ori);
                    out += trail;
                } else {
                    out += lead;
                    out += remap_seg_name_(core);
                    out += trail;
                }
            }

            return out;
        }

        //   seg1+
        if (segs.size() >= 2 && (segs.back() == '+' || segs.back() == '-')) {
            std::string out = remap_seg_name_(segs.substr(0, segs.size() - 1));
            out.push_back(segs.back());
            return out;
        }

        return remap_seg_name_(segs);
    }


    void push_example_(
        std::vector<std::pair<std::string, std::string>>& dst,
        const std::string& before,
        const std::string& after
    ) {
        if (dst.size() < max_log_examples_) {
            dst.emplace_back(before, after);
        }
    }

    void print_one_log_(
        const char* title,
        uint32_t total,
        const std::vector<std::pair<std::string, std::string>>& ex
    ) const {
        if (total == 0) return;

        warning_stream() << "    ! " << title << ": " << total << "\n";
        for (const auto& kv : ex) {
            warning_stream() << "      ! " << kv.first << " -> " << kv.second << "\n";
        }
        if (total > ex.size()) {
            warning_stream() << "      ! ... (" << (total - ex.size()) << " more)\n";
        }
    }
};

} // namespace gfa_name_check
