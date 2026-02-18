#include "../include/gfa_seq.hpp"
#include "../include/save.hpp"
#include "../include/logger.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>

void GfaSeq::trim_(std::string& s) {
    size_t i = 0, j = s.size();
    while (i < j && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    while (j > i && std::isspace(static_cast<unsigned char>(s[j-1]))) --j;
    s.assign(s.substr(i, j - i));
}

std::vector<std::string>
GfaSeq::extract_from_paths_(const std::vector<std::string>& paths) const {
    std::vector<std::string> results;
    results.reserve(paths.size());

    for (const auto& path : paths) {
        struct Item { uint32_t sid; bool rev; };
        std::vector<Item> items; items.reserve(16);

        const size_t n = path.size();
        size_t i = 0;

        // Automatically determine: if there is a comma, use comma-separated parsing; otherwise, fall back to prefix marker parsing
        const bool has_commas = (path.find(',') != std::string::npos);

        while (i < n) {
            // 1) Skip whitespace and (if using comma mode) commas
            while (i < n && (std::isspace(static_cast<unsigned char>(path[i])) || (has_commas && path[i] == ','))) ++i;
            if (i >= n) break;

            // 2) Optional prefix direction marker
            bool have_prefix_sign = false;
            bool prefix_rev = false; // '<' -> true(reverse), '>' -> false(forward)
            if (!has_commas && is_marker_(path[i])) {
                have_prefix_sign = true;
                prefix_rev = (path[i] == '<');
                ++i;
                // Skip whitespace after prefix
                while (i < n && std::isspace(static_cast<unsigned char>(path[i]))) ++i;
            } else if (has_commas && is_marker_(path[i])) {
                // Even in comma mode, allow prefix marker at the start of an item (tolerant to mixed input)
                have_prefix_sign = true;
                prefix_rev = (path[i] == '<');
                ++i;
                while (i < n && std::isspace(static_cast<unsigned char>(path[i]))) ++i;
            }

            // 3) Get token: comma mode -> read until next comma; otherwise -> read until next prefix marker
            size_t j = i;
            if (has_commas) {
                while (j < n && path[j] != ',') ++j;
            } else {
                while (j < n && !is_marker_(path[j])) ++j;
            }
            std::string token = path.substr(i, j - i);
            trim_(token);

            // 4) Skip empty token
            if (token.empty()) {
                i = j + (has_commas && j < n && path[j] == ',' ? 1 : 0);
                continue;
            }

            // 5) Analyze trailing '+' / '-' (if present, take precedence in direction)
            bool use_suffix = false;
            bool suffix_rev = false;
            if (!token.empty() && is_suffix_sign_(token.back())) {
                char last = token.back();
                token.pop_back();
                trim_(token);
                use_suffix = true;
                suffix_rev = (last == '-'); // '-' -> reverse
            }

            // 6) Now the token is the "pure segment name" (without direction markers, commas, or whitespace)
            //    Note: If the segment name itself ends with '+' or '-', use prefix marker syntax to disambiguate (here, trailing '+'/'-' is treated as direction)
            if (!token.empty()) {
                uint64_t sid64 = 0;
                try {
                    sid64 = getNodeInternalId(token); // If not exists, throw exception
                } catch (...) {
                    results.emplace_back("*");
                    goto next_path;
                }
                bool rev = false; // Default forward
                if (use_suffix) {
                    rev = suffix_rev;
                } else if (have_prefix_sign) {
                    rev = prefix_rev;
                } else {
                    rev = false;  // Default forward
                }
                items.push_back({ static_cast<uint32_t>(sid64), rev });
            }

            // 7) Continue
            i = j;
            // In comma mode, j falls on the comma position, need to skip in the top while
        }

        // ===== Concatenate output =====
        {
            // Estimate length
            size_t total_len = 0;
            for (const auto &it : items) {
                const GfaNode* node = getNode(it.sid);
                if (node && !node->sequence.empty() && node->sequence != "*")
                    total_len += node->sequence.size();
            }
            std::string out;
            out.reserve(total_len);

            bool unknown = false;
            for (const auto &it : items) {
                uint32_t vtx = (it.sid << 1) | (it.rev ? 1u : 0u);
                std::string s = get_oriented_sequence(vtx); // If "*", then empty
                if (s.empty()) { unknown = true; break; }
                out += s;
            }
            results.emplace_back(unknown ? "*" : std::move(out));
        }

        next_path: ;
    }

    return results;
}

void GfaSeq::extract_from_file(const std::string& path_file, std::vector<std::string>& paths, std::vector<std::string>& seqs) const {
    std::ifstream fin(path_file);
    if (!fin) {
        error_stream() << path_file << ": No such file or directory\n";
        std::exit(1);
    }

    paths.reserve(1024);

    std::string line;
    while (std::getline(fin, line)) {
        trim_(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue;
        paths.push_back(std::move(line));
    }
    seqs = extract_from_paths_(paths);
    return;
}

void GfaSeq::save_to_file(const std::string& out_file, const std::vector<std::string>& paths, const std::vector<std::string>& seqs) {
    if (paths.size() != seqs.size()) {
        error_stream() << "paths.size() != seqs.size() (" << paths.size() << " != " << seqs.size() << ")\n";
        std::exit(1);
    }

    SAVE saver(out_file);

    for (size_t i = 0; i < seqs.size(); ++i) {
        saver.save(">" + paths[i] + "\n" + seqs[i] + "\n");
    }
    return;
}