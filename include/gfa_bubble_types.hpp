#pragma once
#include <vector>
#include <string>
#include <utility>

#include "minimizer_dna.hpp"

namespace GfaBubble {

enum class Type : uint8_t {
    Tip,        // bubble formed by dead-end branches
    InDel,      // insertion/deletion bubble
    Normal,     // ordinary alternative path bubble
    SameSource, // Homologous paths from the same source
    DiffSource  // Homologous paths from different sources
};

struct Bubble {
    Bubble() = default;
    Bubble(uint32_t s, uint32_t t, std::vector<std::vector<uint32_t>> p, Type type, uint64_t len = 0)
        : source_(s), sink_(t), paths_(std::move(p)), type_(type), len_(len) {}

    // Setters
    void set_source(uint32_t s) { source_ = s; }
    void set_sink(uint32_t t) { sink_ = t; }
    void set_paths(const std::vector<std::vector<uint32_t>>& p) { paths_ = p; }
    void set_path_clusters(std::vector<uint32_t> clusters) { path_clusters_ = std::move(clusters); }
    void set_type(Type t) { type_ = t; }
    void set_len(uint64_t len) { len_ = len; }

    // Getters
    uint32_t get_source() const noexcept { return source_; }
    uint32_t get_sink()   const noexcept { return sink_;   }
    uint64_t get_len()    const noexcept { return len_;    }
    const std::vector<std::vector<uint32_t>>& get_paths() const noexcept { return paths_; }
    const std::vector<uint32_t>& get_path_clusters() const noexcept { return path_clusters_; }
    std::string get_type() const {
        switch (type_) {
            case Type::Tip:    return "Tip";
            case Type::InDel:  return "Indel";
            case Type::Normal: return "Normal";
        }
        return "unknown";
    }

private:
    uint32_t source_{UINT32_MAX};  // [31:1]=segment id, [0]=rev, source vertex id
    uint32_t sink_{UINT32_MAX};
    std::vector<std::vector<uint32_t>> paths_;
    std::vector<uint32_t> path_clusters_;
    Type type_{Type::Normal};  // "Tip", "InDel", "Normal"
    uint64_t len_{0};
};

struct ForkGroup {
public:
    ForkGroup() = default;

    explicit ForkGroup(uint32_t src, std::vector<uint32_t> brs, Type type = Type::Normal)
        : source_(src), branches_(std::move(brs)), type_(type) {}

    // Setters
    void set_source(uint32_t s) { source_ = s; }
    void set_branches(const std::vector<uint32_t>& brs) { branches_ = brs; }
    void set_type(Type t) { type_ = t; }
    void add_branch(uint32_t v) { branches_.push_back(v); }

    // Getters
    uint32_t get_source() const noexcept { return source_; }
    const std::vector<uint32_t>& get_branches() const noexcept { return branches_; }
    std::string get_type() const {
        switch (type_) {
            case Type::Tip:    return "Tip";
            case Type::InDel:  return "Indel";
            case Type::Normal: return "Normal";
        }
        return "unknown";
    }
    std::size_t size() const noexcept { return branches_.size(); }

    void clear() { branches_.clear(); }

private:
    uint32_t source_{UINT32_MAX};     // [31:1]=segment id, [0]=rev, source vertex id
    std::vector<uint32_t> branches_;  // first-hop target vertex ids (These are strictly branches connected to the same node, i.e., their incoming edges are unique and identical.)
    Type type_{Type::Normal};         // "Tip", "InDel", "Normal"
};


struct VcfRecord {
    std::string chrom;
    uint32_t pos{0};
    uint32_t end{0};
    std::string id;
    std::string ref;
    std::string alt;
    std::string src;
    std::string snk;
    uint32_t ns{0};
    uint64_t ln{0};
    uint32_t nc{0};
    std::string vt;
    std::vector<std::string> path_names;
    std::vector<std::string> gt_tokens;         // path-level allel ids
    std::vector<std::string> sample_gt_tokens;  // sample-level GT columns
    char gt_sep{'/'};
};


// For homologous path finding (v0.1.3-r11, 2026-06-13)
struct BubbleBranchPairs_ {
    std::unordered_map<uint32_t, std::vector<uint32_t>> seg_to_groups;  // seg_id -> bubble group ids

    bool empty() const {
        return seg_to_groups.empty();
    }
};


struct NodeSketchEntry {
    Vertex vtx{Vertex(UINT32_MAX)};
    uint32_t len{0};
    minimizerdna::Sketch sketch;
};
struct NodeGroupEntry {
    std::vector<Vertex> vtxs;
    double min_similarity{0.0};
};

struct HomologousParam {
    uint32_t min_source_len{100000};
    double   min_similarity{0.9};
    uint64_t min_path_bp{100000};  // The minimum total length of the path
    Type type{Type::SameSource};

    HomologousParam(uint32_t min_source_len, double min_similarity, uint64_t min_path_bp, Type type)
        : min_source_len(min_source_len), min_similarity(min_similarity), min_path_bp(min_path_bp), type(type) {}

    void set_preset_sameSource() {
        min_source_len = 1e5;
        min_similarity = 0.8;
        min_path_bp = 1e5;
        type = Type::SameSource;
    }
    void set_preset_diffSource() {
        min_source_len = 1e6;
        min_similarity = 0.9;
        min_path_bp = 3e6;
        type = Type::DiffSource;
    }
};
struct HomologousPath {
    std::vector<Vertex> sources;
    std::vector<Vertex> starts;
    std::vector<Vertex> ends;
    std::vector<std::uint64_t> lens;  // The size of lens is equal to sources/starts/ends.
    std::vector<std::vector<Vertex>> paths;  // <start1->end1, start2->end2, ...>
    std::vector<uint32_t> path_clusters;
    double min_hi_similarity = 0.0;
    double max_hi_similarity = 0.0;
    Type type = Type::SameSource;  // "SameSource" or "DiffSource"

    const std::vector<std::vector<Vertex>>& get_paths() const noexcept { return paths; }
    const std::vector<uint32_t>& get_path_clusters() const noexcept { return path_clusters; }
    std::vector<std::vector<uint32_t>> get_paths_u32() const {
        std::vector<std::vector<uint32_t>> paths_u32;
        paths_u32.reserve(paths.size());
        for (const auto& p : paths) {
            std::vector<uint32_t> pu32;
            pu32.reserve(p.size());
            for (const auto& v : p) pu32.push_back(v.vertex_id());
            paths_u32.push_back(std::move(pu32));
        }
        return paths_u32;
    }

    std::string get_type() const {
        switch (type) {
            case Type::SameSource: return "SameSource";
            case Type::DiffSource: return "DiffSource";
            default: return "unknown";
        }
    }
};

struct HomoTask_ {
    Type type{Type::SameSource};
    uint32_t topo_rank{UINT32_MAX};
    uint32_t tie_id{UINT32_MAX};
    std::vector<Vertex> srcs;
};

} // namespace GfaBubble
