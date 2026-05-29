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

    void load_from_GFA(const std::vector<std::string>& filenames);

    /* ---------- accessors ---------- */
    const GfaNode*                       getNode(uint64_t internal_id)              const { return internal_id < nodes_.size() ? &nodes_[internal_id] : nullptr; }
    uint64_t                             getNodeInternalId(const std::string& name) const;
    const std::string                    getNodeName(uint64_t internal_id)          const { return internal_id < nodes_.size() ? nodes_[internal_id].name : "*"; }
    const uint32_t                       getNodeLength(uint64_t internal_id)        const { return internal_id < nodes_.size() ? nodes_[internal_id].length : 0; }
    const bool                           getNodeDeleted(uint64_t internal_id)       const { return internal_id < nodes_.size() ? nodes_[internal_id].deleted : true; }
    const GfaArc&                        getArc(size_t arc_idx)                     const { if (arc_idx < arcs_.size()) return arcs_[arc_idx]; throw std::out_of_range("Arc index out of bounds"); }
    const std::vector<GfaArc>&           getAllArcs()                               const { return arcs_; }
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

    std::string format_overlap_field(int32_t ov, int32_t ow) const;  // Format overlap field ov/ow as in GFA

    // Get outgoing/incoming arcs of a vertex, excluding deleted arcs. If skip_self is true, also exclude self-loop arcs.
    std::vector<size_t> getArcsIdxFromVertex(uint32_t v, bool skip_self = false) const;
    std::vector<const GfaArc*> getArcsFromVertex(uint32_t v, bool skip_self = false) const;
    std::vector<size_t> getArcsIdxToVertex(uint32_t v, bool skip_self = false) const;
    std::vector<const GfaArc*> getArcsToVertex(uint32_t v, bool skip_self = false) const;
    
    // Fetch oriented string (forward if ori=0; reverse-complement if ori=1)
    std::string get_oriented_sequence(Vertex vertex, uint32_t ow=0) const;

    // Get the assembled sequence of a path defined by a chain of vertex IDs, using the ow rule for overlaps. Return "*" if any segment is missing or has "*".
    std::string get_path_sequence(const std::vector<uint32_t>& path_vtx_ids) const;
    std::string get_path_sequence(const std::vector<Vertex>& path_vtxs) const;
    std::vector<std::uint32_t> get_path_overlaps(const std::vector<uint32_t>& path_vtx_ids) const;
    std::vector<std::uint32_t> get_path_overlaps(const std::vector<Vertex>& path_vtxs) const;

    // Get the overlap length (ow) of the edge from from_vtx to to_vtx, or 0 if no such edge exists or if it's deleted.
    uint32_t get_edge_ow(uint32_t from_vtx_id, uint32_t to_vtx_id) const;
    uint32_t get_edge_ow(Vertex from_vtx, Vertex to_vtx) const;

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
        const uint64_t DFS_guard, const uint32_t stall_round_limit
    ) const;
    std::vector<std::vector<uint32_t>> enumerate_paths_greedy_DFS(
        const uint32_t src, const uint32_t sink, const std::unordered_set<uint32_t>& region_set,
        const uint32_t max_depth, const uint32_t max_paths,
        const bool skip_comp, bool& hit_limits, const uint64_t DFS_guard, const uint32_t stall_round_limit
    ) const;

    /**
     * @brief Build a connectivity index for all nodes in the graph, where each node is assigned a component ID based on its connectivity.
     * 
     *   - Nodes that are reachable from each other (ignoring edge directions) will share the same component ID.
     * 
     */
    bool build_nodes_connectivity_index(bool skip_comp = false);
    const std::vector<uint32_t>& get_nodes_connectivity_index() const { return connectivity_index_; }
    uint32_t node_connectivity_id(Vertex v) const;
    bool nodes_connected(Vertex a, Vertex b) const;

    /**
     * @brief Build a topological index for all vertices in the graph, where each vertex is annotated with its strongly connected component (SCC) ID, topological rank, and cycle membership.
     * 
     */
    bool build_vertex_topological_index(bool skip_comp = false);
    const std::vector<GfaTopoIndex>& get_vertex_topological_index() const { return topo_index_; }
    uint32_t vertex_topo_rank(Vertex v) const;
    const GfaTopoIndex& vertex_topo_info(Vertex v) const;

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

    /* ---------- connectivity & topological index ---------- */
    std::vector<uint32_t> connectivity_index_;  // node_id -> component_id
    std::vector<GfaTopoIndex> topo_index_;      // vertex_id -> { scc_id, topo_rank, scc_size, in_cycle }


protected:
    /**
     * @brief Get the segment ID for a given name, or create a new segment entry if it does not exist.
     *
     * @param name    Segment name.
     * @param is_new  Optional output flag; set to true if a new segment was created, false otherwise.
     *
     * @return Segment ID corresponding to @p name.
     */
    uint32_t get_or_add_segment(const std::string& name, bool* is_new = nullptr);
    
    
    /**
     * @brief Extract a sequence slice from a segment, or return "*" if unavailable/invalid.
     *
     * The slice is taken from [beg, end) on the forward strand. If @p is_rev is true,
     * the extracted subsequence is reverse-complemented before returning.
     *
     * @param nodes    Node table.
     * @param seg_id   Segment ID.
     * @param beg      Begin position (inclusive).
     * @param end      End position (exclusive).
     * @param is_rev   Whether to return the reverse-complemented slice.
     *
     * @return Extracted subsequence, or "*" if the source sequence is missing or the range is invalid.
     */
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