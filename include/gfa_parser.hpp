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

    void load_from_GFA(
        const std::vector<std::string>& filenames,
        const std::vector<std::string>& sample_names = {},
        bool read_links = true,
        const std::vector<uint32_t>& duplicate_groups = {}
    );

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

    // Enumerate paths (DFS)
    std::vector<std::vector<uint32_t>> enumerate_paths_greedy_DFS(
        const uint32_t src, const uint32_t sink, const std::unordered_set<uint32_t>& region_set,
        const uint32_t max_depth, const uint32_t max_paths,
        const bool skip_comp, bool& hit_limits, const uint64_t DFS_guard, const uint32_t stall_round_limit,
        const uint32_t required_sample = UINT32_MAX
    ) const;
    std::vector<std::vector<uint32_t>> open_walk(
        const uint32_t src,
        const std::unordered_set<uint32_t>& region_set,
        const uint32_t max_depth,
        const bool skip_comp,
        const uint64_t DFS_guard,
        const uint64_t walk_bp,
        const std::unordered_set<uint32_t>* blocked_seg = nullptr,
        const uint32_t required_sample = UINT32_MAX
    ) const;

    // Connectivity index
    bool build_nodes_connectivity_index(bool skip_comp = false);
    const std::vector<uint32_t>& get_nodes_connectivity_index() const { return connectivity_index_; }
    uint32_t node_connectivity_id(Vertex v) const;
    uint32_t node_connectivity_size(Vertex v) const;
    bool nodes_connected(Vertex a, Vertex b) const;

    // Topological index
    bool build_vertex_topological_index(bool skip_comp = false);
    const std::vector<GfaTopoIndex>& get_vertex_topological_index() const { return topo_index_; }
    uint32_t vertex_topo_rank(Vertex v) const;
    bool vertex_in_cycle(Vertex v) const;
    const GfaTopoIndex& vertex_topo_info(Vertex v) const;

    // Complexity detection
    bool getNodeComplex(uint32_t sid) const { return sid < nodes_.size() && !nodes_[sid].deleted && nodes_[sid].is_complex; }
    void detect_complex_regions(
        uint16_t branch_degree,
        uint16_t hub_degree,
        uint32_t min_nodes,
        uint32_t min_branches,
        bool skip_comp = false
    );
    void mark_segments_complex(const std::vector<uint32_t>& seg_ids);

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
    void save_to_disk(
        const std::string& filename,
        bool write_paths,
        bool write_align,
        bool write_seq,
        const std::string& command_line = ""
    ) const;

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
    std::unordered_map<std::string, uint64_t>  name_to_id_map_;     // segment name -> internal id (index in nodes_)
    std::vector<std::string>                   sample_id_to_name_;  // sample id -> sample name (tags in S-lines and from -n/--name)
    std::unordered_map<std::string, uint32_t>  sample_name_to_id_;  // sample name -> sample id (index in sample_id_to_name_)

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
    std::vector<uint32_t> connectivity_sizes_;  // component_id -> segment count
    std::vector<GfaTopoIndex> topo_index_;      // vertex_id -> { scc_id, topo_rank, scc_size, in_cycle }


protected:
    static bool input_segment_names_need_namespace_(
        const std::vector<std::string>& filenames,
        const std::vector<uint32_t>& duplicate_groups
    );

    uint32_t get_or_add_segment(const std::string& name, bool* is_new = nullptr);

    static std::string slice_seq_or_star_(const std::vector<GfaNode>& nodes, uint32_t seg_id, uint32_t beg, uint32_t end, bool is_rev);

    // Add an arc to the graph, ensuring synchronization between arcs_ and link_aux_ vectors.
    GfaArc* add_arc(uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, bool comp_flag);

    void finalize_();  // Finalize the graph after loading, including sorting arcs, building index, etc.


    /* line parsers */
    bool parseSLine(std::stringstream& ss, const std::string& sample_name);
    bool parseLLine(std::stringstream& ss, const std::string& filename, uint64_t line_no);
    bool parsePLine(std::stringstream& ss, const std::string& filename, uint64_t line_no);
    bool parseALine(std::stringstream& ss, const std::string& filename, uint64_t line_no);


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

    uint32_t add_segment(const std::string& name, const std::string& seq, bool name_check=true, const std::vector<uint32_t>& sample_ids = {});  // Add a segment to the graph (S-line)
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


protected:  // Functions will be used in SN Tag.
    bool has_node_source_tag_{false};

    // Used to identify samples from S-line tags and -n/--name, and to manage sample ID mappings. (v0.1.3-r9 2026-06-03)
    uint32_t intern_sample_(const std::string& name);

    std::vector<uint32_t> parse_sample_csv_(const std::string& csv); 
    std::string sample_ids_to_csv_(const std::vector<uint32_t>& ids) const;

    bool hasNodeSourceTag() const { return has_node_source_tag_; }

public:  // Functions will be used in SN Tag.
    void merge_sample_ids_into(std::vector<uint32_t>& dst, const std::vector<uint32_t>& src) const;
    void merge_variant_tags_into(GfaAux& dst, const GfaAux& src) const;

    bool sample_ids_intersect(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) const;
    bool sample_ids_equal_or_empty(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) const;

    bool intersect_node_sample(uint32_t a_seg, uint32_t b_seg) const;
    bool intersect_vertex_sample(uint32_t a_vtx, uint32_t b_vtx) const;

    bool same_node_sample(uint32_t a_seg, uint32_t b_seg) const;
    bool same_vertex_sample(uint32_t a_vtx, uint32_t b_vtx) const;

    const std::vector<std::string>& getSampleNames() const;

    const std::vector<uint32_t>& getNodeSampleIds(uint32_t seg_id) const;
};
