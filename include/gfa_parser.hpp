#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <functional>
#include "bindings/cpp/WFAligner.hpp"

#include <iterator>
#include <unordered_set>
#include <utility>

#include "CIGAR.hpp"
#include "seq_utils.hpp"
#include "gfa_parser_types.hpp"
#include "seg_replace.hpp"
#include "logger.hpp"
#include "GZ_chunk_reader.hpp"
#include "save.hpp"
#include "get_time.hpp"

#include <iomanip>
#include <cstring>
#include <stack>
#include <queue>


/*============================================================*/
/*                         GfaGraph                            */
/*============================================================*/
class GfaGraph {
public:
    friend class GfaBamLoader;

    GfaGraph();
    ~GfaGraph();

    void load_from_GFA(const std::string& filename);  // Load GFA file into the graph structure

    /* ---------- accessors ---------- */
    const GfaNode*                       getNode(uint64_t internal_id)              const { return internal_id < nodes_.size() ? &nodes_[internal_id] : nullptr; }
    uint64_t                             getNodeInternalId(const std::string& name) const;
    const std::string                    getNodeName(uint64_t internal_id)          const { return internal_id < nodes_.size() ? nodes_[internal_id].name : std::to_string(internal_id); }
    const uint32_t                       getNodeLength(uint64_t internal_id)        const { return internal_id < nodes_.size() ? nodes_[internal_id].length : 0; }
    const bool                           getNodeDeleted(uint64_t internal_id)       const { return internal_id < nodes_.size() ? nodes_[internal_id].deleted : true; }
    const std::vector<GfaArc>&           getAllArcs()                               const { return arcs_; }
    std::vector<size_t>                  getArcsIdxFromVertex(uint32_t v, bool skip_self = false)           const;  // indices of outgoing arcs per vertex, excluding deleted arcs
    std::vector<const GfaArc*>           getArcsFromVertex(uint32_t v, bool skip_self = false)              const;  // outgoing arcs per vertex, excluding deleted arcs
    std::vector<size_t>                  getArcsIdxToVertex(uint32_t v, bool skip_self = false)             const;  // indices of incoming arcs per vertex, excluding deleted arcs
    std::vector<const GfaArc*>           getArcsToVertex(uint32_t v, bool skip_self = false)                const;  // incoming arcs per vertex, excluding deleted arcs
    const std::vector<GfaPath>&          getPaths()                                 const { return paths_; }
    const std::vector<GfaAlignment>&     getAlignments()                            const { return alignments_; }
    size_t                               getNumNodes()                              const { return total_segments_; }
    size_t                               getNumUniqLinks()                          const { return total_uniq_links_; }
    size_t                               getNumPaths()                              const { return paths_.size(); }
    size_t                               getNumAlignments()                         const { return alignments_.size(); }
    size_t                               getTotalSegmentLength()                    const { return total_segment_length_; }
    double                               getAveSegmentsLength()                     const { return total_segments_ ? static_cast<double>(total_segment_length_) / total_segments_ : 0.0f; }
    uint32_t                             getMaxRank()                               const { return max_rank_; }
    uint64_t                             getTotalDeg()                              const { return total_deg_; }
    uint64_t                             getMaxDeg()                                const { return max_deg_; }
    double                               getAveDeg()                                const { return arc_indexs_.size() ? static_cast<double>(total_deg_) / arc_indexs_.size() : 0.0f; }
    ExpandedSeqs                         getSeqVec(uint32_t k = 0)                  const;
    std::vector<std::string>             getAllSegmentNames()                       const;

public:
    /* ---------- setters ---------- */
    void set_forbid_overlap(bool v) { forbid_overlap_ = v; }

    // Format overlap field ov/ow as in GFA
    std::string format_overlap_field(int32_t ov, int32_t ow) const;
    // Fetch oriented string (forward if ori=0; reverse-complement if ori=1)
    std::string get_oriented_sequence(uint32_t vertex) const;
    std::string get_path_sequence(const std::vector<uint32_t>& path) const;

    // traversal helpers
    // start_v -> vertex-id (segment<<1 | ori)
    // void dfs_from(uint32_t start_v, const std::function<void(uint32_t)>& visit) const;
    // void bfs_from(uint32_t start_v, const std::function<void(uint32_t)>& visit) const;
    std::vector<PathSequence> walk_dfs(
        uint32_t start_v,
        bool to_direction,      // false=From, true=To
        uint32_t max_paths=0,
        uint32_t max_steps=0,
        uint32_t max_bases=0
    ) const;

    std::vector<PathSequence> walk_bfs(
        uint32_t start_v,
        bool to_direction,      // false=From, true=To
        uint32_t max_paths=0,
        uint32_t max_steps=0,
        uint32_t max_bases=0
    ) const;
    
    // Enumerate paths from `src` to `sink` using DFS up to `max_depth` and `max_paths`.
    std::vector<std::vector<uint32_t>> enumerate_paths_DFS(
        const uint32_t src, const uint32_t sink, const std::unordered_set<uint32_t>& region_set, 
        const uint32_t max_depth, const uint32_t max_paths,
        const bool skip_comp, bool& hit_limits,
        const uint64_t DFS_guard
    ) const;

    bool is_connected_vertex(uint32_t v_from, uint32_t v_to, uint32_t step_cap = 200000) const;

    bool shortest_distance_between_offsets(
        uint32_t v_from, uint32_t off_from,
        uint32_t v_to,   uint32_t off_to,
        uint64_t& out_dist,
        uint32_t step_cap = 200000
    ) const;

    /* ---------- debug ---------- */
    void print_graph_stats() const;
    void printArcList()      const;
    void printArcIndex()     const;

    // Save the graph to a disk. Optionally write P-lines and A-lines.
    void save_to_disk(const std::string& filename, bool write_paths, bool write_align, bool write_seq) const;

protected:
    /* ---------- graph properties ---------- */
    bool forbid_overlap_{false};
    bool has_overlap_{false};

    /* ---------- segments ---------- */
    uint64_t                                   total_segments_{0};
    uint64_t                                   total_segment_length_{0};
    uint32_t                                   max_rank_{0};

    /* ---------- alignment ---------- */
    std::vector<GfaAlignment>                  alignments_;
    uint64_t                                   total_alignments_{0};

    /* ---------- nodes ---------- */
    std::vector<GfaNode>                       nodes_;
    std::unordered_map<std::string, uint64_t>  name_to_id_map_;

    /* ---------- arcs & links ---------- */
    std::vector<GfaArc>                        arcs_;  // e.g., 1->2, 1->3, 2->4, 3->4
    std::vector<uint64_t>                      arc_indexs_;  // arc_indexs_[v] = (start position in arcs_ << 32) | count
    std::vector<GfaAux>                        link_aux_;
    uint64_t                                   total_uniq_links_{0};
    uint64_t                                   total_deg_{0};
    uint64_t                                   max_deg_{0};

    /* ---------- paths ---------- */
    std::vector<GfaPath>                       paths_;


protected:
    /* ---------- helpers ---------- */
    // Get or add a segment by name, and return its ID
    uint32_t get_or_add_segment(const std::string& name);
    // Slice the sequence of a segment or return "*" if invalid
    static std::string slice_seq_or_star_(const std::vector<GfaNode>& nodes, uint32_t seg_id, uint32_t beg, uint32_t end, bool is_rev);

    // Add an arc to the graph, ensuring synchronization between arcs_ and link_aux_ vectors.
    GfaArc* add_arc(uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, bool comp_flag);

    void finalize_();  // Finalize the graph after loading, including sorting arcs, building index, etc.

    /* line parsers */
    bool parseSLine(std::stringstream& ss);
    bool parseLLine(std::stringstream& ss);
    bool parsePLine(std::stringstream& ss);
    bool parseALine(std::stringstream& ss);

    /* ---------- finalize steps ---------- */
    void fixNoSeg();     /* step 1 */
    void arcSort();      /* step 2 */
    void arcIndex();     /* step 3 */
    void fixSemiArc();   /* step 4 */
    void fixSymmAdd();   /* step 5 */
    void fixArcLen();    /* step 6 */
    void cleanup();      /* step 7 */
    void countLinks();   /* step 8 */
    void calculateDeg(); /* step 9 */

    uint32_t add_segment(const std::string& name, const std::string& seq, bool name_check=true);  // Add a segment to the graph (S-line)
    bool delete_segment(uint32_t seg_id);  // Delete a segment and all its incident arcs (both orientations).
    void rebuild_after_edits();  // Deletion or Insertion

    // Extend left from predecessors of incoming edges (suffix only), excluding start
    PathSequence extend_left_from(
        uint32_t start_vertex,
        uint32_t min_bases,
        uint32_t max_steps
    ) const;

    // Extend right from successors of outgoing edges (prefix only), excluding start
    PathSequence extend_right_from(
        uint32_t start_vertex,
        uint32_t min_bases,
        uint32_t max_steps
    ) const;
};