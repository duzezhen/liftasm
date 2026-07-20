#include "../include/gfa_seq.hpp"
#include "../include/save.hpp"
#include "../include/logger.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <unordered_set>

void GfaSeq::trim_(std::string& s) {
    size_t i = 0, j = s.size();
    while (i < j && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    while (j > i && std::isspace(static_cast<unsigned char>(s[j-1]))) --j;
    s.assign(s.substr(i, j - i));
}

bool GfaSeq::parse_open_walk_line_(const std::string& line, OpenWalkRequest& req) {
    std::istringstream iss(line);
    std::vector<std::string> cols;
    std::string col;
    while (iss >> col) cols.push_back(std::move(col));

    if (cols.size() != 3 || (cols[1] != "+" && cols[1] != "-")) return false;

    size_t used = 0;
    uint64_t len = 0;
    try {
        len = std::stoull(cols[2], &used);
    } catch (...) {
        return false;
    }
    if (used != cols[2].size()) return false;

    req.segment = cols[0];
    req.rev = (cols[1] == "-");
    req.len = len;
    return true;
}

std::string GfaSeq::extract_from_path_(const std::string& path) const {
    std::vector<uint32_t> vertices;
    vertices.reserve(16);

    const size_t n = path.size();
    size_t i = 0;
    const bool has_commas = (path.find(',') != std::string::npos);

    while (i < n) {
        while (i < n && (std::isspace(static_cast<unsigned char>(path[i])) || (has_commas && path[i] == ','))) ++i;
        if (i >= n) break;

        bool have_prefix_sign = false;
        bool prefix_rev = false;
        if (is_marker_(path[i])) {
            have_prefix_sign = true;
            prefix_rev = (path[i] == '<');
            ++i;
            while (i < n && std::isspace(static_cast<unsigned char>(path[i]))) ++i;
        }

        size_t j = i;
        if (has_commas) {
            while (j < n && path[j] != ',') ++j;
        } else {
            while (j < n && !is_marker_(path[j])) ++j;
        }

        std::string token = path.substr(i, j - i);
        trim_(token);

        if (!token.empty()) {
            bool use_suffix = false;
            bool suffix_rev = false;
            if (is_suffix_sign_(token.back())) {
                suffix_rev = (token.back() == '-');
                token.pop_back();
                trim_(token);
                use_suffix = true;
            }

            try {
                const uint32_t sid = static_cast<uint32_t>(getNodeInternalId(token));
                const bool rev = use_suffix ? suffix_rev : (have_prefix_sign ? prefix_rev : false);
                vertices.push_back(Vertex(sid, rev).vertex_id());
            } catch (...) {
                return "*";
            }
        }

        i = j;
    }

    return get_path_sequence(vertices);
}

std::vector<std::string> GfaSeq::extract_from_paths_(const std::vector<std::string>& paths) const {
    std::vector<std::string> results;
    results.reserve(paths.size());
    for (const auto& path : paths) {
        results.push_back(extract_from_path_(path));
    }
    return results;
}

std::string GfaSeq::extract_open_walk_(const OpenWalkRequest& req) const {
    uint32_t sid = UINT32_MAX;
    try {
        sid = static_cast<uint32_t>(getNodeInternalId(req.segment));
    } catch (...) {
        return "*";
    }

    const uint32_t src = Vertex(sid, req.rev).vertex_id();
    const std::unordered_set<uint32_t> region_set;
    const std::vector<std::vector<uint32_t>> walks = open_walk(
        src,
        region_set,
        /*max_depth=*/0,
        /*skip_comp=*/false,
        /*DFS_guard=*/1000000,
        req.len
    );

    if (walks.empty()) return "*";
    return get_path_sequence(walks.front());
}

void GfaSeq::extract_from_file(const std::string& path_file, std::vector<std::string>& paths, std::vector<std::string>& seqs) {
    std::ifstream fin(path_file);
    if (!fin) {
        error_stream() << path_file << ": No such file or directory\n";
        std::exit(1);
    }

    paths.reserve(1024);
    seqs.reserve(1024);

    bool topo_ready = false;
    std::string line;
    while (std::getline(fin, line)) {
        trim_(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue;

        OpenWalkRequest req;
        if (parse_open_walk_line_(line, req)) {
            if (!topo_ready) {
                build_vertex_topological_index();
                topo_ready = true;
            }
            seqs.push_back(extract_open_walk_(req));
        } else {
            seqs.push_back(extract_from_path_(line));
        }

        paths.push_back(std::move(line));
    }
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
