#pragma once
#include <vector>
#include <string>
#include <utility>

namespace GfaBubble {

// Bubble type enumeration
enum class Type : uint8_t {
    Tip,        // bubble formed by dead-end branches
    InDel,      // insertion/deletion bubble
    Normal      // ordinary alternative path bubble
};

struct Bubble {
    Bubble() = default;
    Bubble(uint32_t s, uint32_t t, std::vector<std::vector<uint32_t>> p, Type type)
        : source_(s), sink_(t), paths_(std::move(p)), type_(type) {}

    // Setters
    void set_source(uint32_t s) { source_ = s; }
    void set_sink(uint32_t t) { sink_ = t; }
    void set_paths(const std::vector<std::vector<uint32_t>>& p) { paths_ = p; }
    void set_type(Type t) { type_ = t; }

    // Getters
    uint32_t get_source() const noexcept { return source_; }
    uint32_t get_sink()   const noexcept { return sink_;   }
    const std::vector<std::vector<uint32_t>>& get_paths() const noexcept { return paths_; }
    std::string get_type() const {
        switch (type_) {
            case Type::Tip:    return "Tip";
            case Type::InDel:  return "Indel";
            case Type::Normal: return "Normal";
        }
        return "unknown";
    }

private:
    uint32_t source_{0};
    uint32_t sink_{0};
    std::vector<std::vector<uint32_t>> paths_;
    Type type_{Type::Normal};  // "Tip", "InDel", "Normal"
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
    uint32_t source_{0};              // [31:1]=segment id, [0]=rev, source vertex id
    std::vector<uint32_t> branches_;  // first-hop target vertex ids (These are strictly branches connected to the same node, i.e., their incoming edges are unique and identical.)
    Type type_{Type::Normal};         // "Tip", "InDel", "Normal"
};

} // namespace GfaBubble