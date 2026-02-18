#include "../include/gfa_gfa2fa.hpp"
#include "../include/save.hpp"
#include <fstream>
#include <algorithm>
#include <sstream>

/* ---------- static / private helpers ---------- */

uint32_t Gfa2fa::real_len_of_(const GfaNode& n) {
    if (n.sequence.empty() || n.sequence == "*") return 0u;
    if (n.length == 0) return static_cast<uint32_t>(n.sequence.size());
    return std::min<uint32_t>(n.length, static_cast<uint32_t>(n.sequence.size()));
}

void Gfa2fa::write_fasta_record_(
    std::ostream& out,
    const std::string& name,
    const std::string& seq,
    uint32_t wrap_width
) {
    out << '>' << name << '\n';
    if (wrap_width == 0) { out << seq << '\n'; return; }
    for (size_t i = 0; i < seq.size(); i += wrap_width) {
        const size_t take = std::min<size_t>(wrap_width, seq.size() - i);
        out.write(seq.data() + i, take);
        out.put('\n');
    }
}

bool Gfa2fa::node_good_(uint32_t seg) const {
    const GfaNode* nd = getNode(seg);
    return nd && !nd->deleted;
}

/* ---------- left/right extensions: return extension only ---------- */

std::string Gfa2fa::extend_left_seq_(uint32_t v_forward, uint32_t need_ybp) const {
    PathSequence p = extend_left_from(v_forward, need_ybp, /*max_steps=*/0);
    return std::move(p.sequence);
}

std::string Gfa2fa::extend_right_seq_(uint32_t v_forward, uint32_t need_ybp) const {
    PathSequence p = extend_right_from(v_forward, need_ybp, /*max_steps=*/0);
    return std::move(p.sequence);
}

/* ---------- public I/O ---------- */

void Gfa2fa::dump_to_stream(std::ostream& out) const {
    const uint32_t N = static_cast<uint32_t>(getNumNodes());
    for (uint32_t seg = 0; seg < N; ++seg) {
        const GfaNode* nd = getNode(seg);
        if (!nd || nd->deleted) continue;

        const uint32_t core_len = real_len_of_(*nd);
        if (skip_unknown_ && core_len == 0) continue;

        std::string core;
        if (!nd->sequence.empty() && nd->sequence != "*") {
            core.assign(nd->sequence.data(), core_len);
        }

        const bool do_extend = (extend_ybp_ > 0) && (core_len < static_cast<uint32_t>(min_len_xbp_));

        std::string left, right;
        if (do_extend) {
            const uint32_t v_fw = (seg << 1) | 0; // forward vertex of this segment
            left  = extend_left_seq_( v_fw, static_cast<uint32_t>(extend_ybp_) );
            right = extend_right_seq_(v_fw, static_cast<uint32_t>(extend_ybp_) );
        }

        std::string all; all.reserve(left.size() + core.size() + right.size());
        all.append(left);
        all.append(core);
        all.append(right);

        write_fasta_record_(out, nd->name, all, static_cast<uint32_t>(wrap_width_));
    }
}

bool Gfa2fa::dump_to_file(const std::string& fasta_path) const {
    std::ofstream ofs(fasta_path, std::ios::out | std::ios::binary);
    if (!ofs) return false;
    dump_to_stream(ofs);
    ofs.flush();
    return static_cast<bool>(ofs);
}