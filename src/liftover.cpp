#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <ostream>
#include <fstream>
#include <regex>
#include <algorithm>
#include <cctype>
#include <limits>
#include <zlib.h>

#include "../include/logger.hpp"
#include "../include/kio.hpp"
#include "../include/kseq.h"
#include "../include/coordmap.hpp"
#include "../include/liftover.hpp"
#include "../include/progress_tracker.hpp"
#include "../include/ThreadPool.hpp"

namespace liftover {

namespace detail {
    bool parse_u32(std::string_view sv, uint32_t& out) {
        if (sv.empty()) return false;
        uint64_t x = 0;
        for (char c : sv) {
            if (c < '0' || c > '9') return false;
            x = x * 10u + (uint64_t)(c - '0');
            if (x > std::numeric_limits<uint32_t>::max()) return false;
        }
        out = (uint32_t)x;
        return true;
    }

    bool split3_bed(std::string_view line, std::string_view& chrom, std::string_view& s0, std::string_view& s1, std::string_view& rest) {
        auto ws = [](char c){ return c==' ' || c=='\t'; };

        size_t i = 0;
        while (i < line.size() && ws(line[i])) ++i;
        if (i >= line.size() || line[i] == '#') return false;

        size_t a = i; while (i < line.size() && !ws(line[i])) ++i;
        chrom = line.substr(a, i - a);

        while (i < line.size() && ws(line[i])) ++i;
        size_t b = i; while (i < line.size() && !ws(line[i])) ++i;
        if (b == i) return false;
        s0 = line.substr(b, i - b);

        while (i < line.size() && ws(line[i])) ++i;
        size_t c = i; while (i < line.size() && !ws(line[i])) ++i;
        if (c == i) return false;
        s1 = line.substr(c, i - c);

        while (i < line.size() && ws(line[i])) ++i;
        rest = (i < line.size()) ? line.substr(i) : std::string_view{};
        return true;
    }

    std::string hit_label(const coordmap::CoordMap& M, const coordmap::Hit& h) {
        std::string s = M.contig_name(h.ctg);
        s.push_back(':');
        s += std::to_string(h.beg);
        s.push_back('-');
        s += std::to_string(h.end);
        s.push_back(h.rev ? '-' : '+');
        return s;
    }

    static inline char comp_base(char b) {
        switch (b) {
            case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
            case 'R': return 'Y'; case 'Y': return 'R'; case 'S': return 'S'; case 'W': return 'W';
            case 'K': return 'M'; case 'M': return 'K'; case 'B': return 'V'; case 'D': return 'H';
            case 'H': return 'D'; case 'V': return 'B';
            default:  return b;
        }
    }

    void rc_inplace(std::string& s) {
        std::reverse(s.begin(), s.end());
        for (char& c : s) c = comp_base((char)std::toupper((unsigned char)c));
    }
} // namespace detail


RunOpts set_opts(
    std::string map_file, std::string bed_file, std::string ref_file, std::string out_file, 
    std::string regex, double min_frac, 
    uint32_t flank_win, uint32_t max_flank, uint32_t max_gap, uint16_t max_hit,
    std::string paf_file, uint32_t min_mapq, uint32_t min_len,
    bool do_check, uint32_t win, uint32_t step, uint32_t max_examples, 
    int cm_max_hops, uint32_t cm_max_fanout, uint32_t cm_min_len, double cm_min_frac, uint32_t cm_max_total_hits,
    int threads
) {
    RunOpts o;
    o.map_file          = map_file;
    o.bed_file          = bed_file;
    o.ref_file          = ref_file;
    o.out_file          = out_file;
    o.regex             = regex;
    o.min_frac          = min_frac;
    o.flank_win         = flank_win;
    o.max_flank         = max_flank;
    o.max_gap           = max_gap;
    o.max_hit           = max_hit;

    o.paf_file          = paf_file;
    o.min_mapq          = min_mapq;
    o.min_len           = min_len;

    o.do_check          = do_check;
    o.win               = win;
    o.step              = step;
    o.max_examples      = max_examples;
    o.cm_max_hops       = cm_max_hops;
    o.cm_max_fanout     = cm_max_fanout;
    o.cm_min_len        = cm_min_len;
    o.cm_min_frac       = cm_min_frac;
    o.cm_max_total_hits = cm_max_total_hits;

    o.threads           = threads;

    if (o.map_file.empty()) { error_stream() << "--map is required"; std::exit(1); }
    if (o.do_check) {
        if (o.ref_file.empty()) { error_stream() << "--check requires --ref"; std::exit(1); }
    } else {
        if (o.bed_file.empty()) { error_stream() << "need --bed (or use --check)"; std::exit(1); }
    }

    return o;
}

int file_reader(
    const RunOpts& opt, 
    std::unordered_map<std::string, std::vector<BEDinfo>>& beds, 
    coordmap::CoordMap& M, 
    std::unordered_map<std::string, std::vector<Pafinfo>>& pafs,
    std::unordered_map<std::string, std::vector<uint32_t>>& paf_ends,
    std::unordered_map<std::string, std::vector<uint32_t>>& paf_idx_ends
) {
    // BED load
    {
        kio::LineReader lr(opt.bed_file);
        std::string line;
        while (lr.getline(line)) {
            std::string_view chrom, s0, s1, rest;
            if (!detail::split3_bed(line, chrom, s0, s1, rest)) continue;

            uint32_t b=0, e=0;
            if (!detail::parse_u32(s0, b) || !detail::parse_u32(s1, e)) continue;
            if (b > e) std::swap(b, e);
            if (b == e) continue;

            beds[std::string(chrom)].push_back(BEDinfo{b, e, std::string{rest}});
        }

        // Sort
        for (auto& kv : beds) {
            auto& v = kv.second;
            std::sort(v.begin(), v.end(), [](const BEDinfo& a, const BEDinfo& b){
                if (a.beg != b.beg) return a.beg < b.beg;
                return a.end < b.end;
            });
        }
    }

    // MAP load
    {
        M.load(opt.map_file);
        log_stream() << "Loaded records: " << M.num_records() << "\n";
    }

    // PAF load
    {
        auto cigar_check = [](const std::string& cg) ->bool {
            uint64_t num = 0;
            for (char c : cg) {
                if (c >= '0' && c <= '9') { num = num*10 + (uint64_t)(c - '0'); continue; }
                if (num == 0) continue;
                if (c!='M' && c!='I' && c!='D' && c!='=' && c!='X') return false;
                num = 0;
            }
            return true;
        };

        if (!opt.paf_file.empty()) {
            paf::Reader pr(opt.paf_file);
            paf::Record r;
            while (pr.next(r)) {
                if (r.mapq < opt.min_mapq || r.alen < opt.min_len) continue;
                if (r.strand != '+' && r.strand != '-') continue;
                if (r.tname == r.qname) {
                    error_stream() << "PAF record with identical query and target: " << r.qname << std::endl;
                    std::exit(1);
                }

                char tp_type=0; std::string_view tp_val;
                char cg_type=0; std::string_view cg_val;
                if (!paf::find_tag(r, "tp", tp_type, tp_val) || tp_type!='A' || tp_val.empty()) continue;
                if (!paf::find_tag(r, "cg", cg_type, cg_val) || cg_type!='Z' || cg_val.empty()) continue;

                const char tp = tp_val[0];
                if (!(tp=='P' || tp=='I')) continue;

                Pafinfo a;
                a.qlen   = r.qlen;
                a.qstart = r.qstart;
                a.qend   = r.qend;
                a.strand = r.strand;
                a.tname  = r.tname;
                a.tlen   = r.tlen;
                a.tstart = r.tstart;
                a.tend   = r.tend;
                a.nmatch = r.nmatch;
                a.alen   = r.alen;
                a.mapq   = r.mapq;
                a.tp     = tp;
                a.cg.assign(cg_val);
                if (!cigar_check(a.cg)) {
                    error_stream() << "Invalid CIGAR string in PAF record: " << a.cg << std::endl;
                    error_stream() << "Only M/I/D/=/X operations are allowed." << std::endl;
                    std::exit(1);
                }
                pafs[r.qname].push_back(std::move(a));

                // Swap ref and qry
                Pafinfo b;
                b.qlen   = r.tlen;
                b.qstart = r.tstart;
                b.qend   = r.tend;
                b.strand = r.strand;
                b.tname  = r.qname;
                b.tlen   = r.qlen;
                b.tstart = r.qstart;
                b.tend   = r.qend;
                b.nmatch = r.nmatch;
                b.alen   = r.alen;
                b.mapq   = r.mapq;
                b.tp     = tp;
                b.cg     = paf::cigar_swap(cg_val);
                pafs[r.tname].push_back(std::move(b));
            }
        }

        // Sort PAFs by query positions and build end-indexes
        paf_ends.reserve(pafs.size()); paf_idx_ends.reserve(pafs.size());

        for (auto& kv : pafs) {
            const std::string& qname = kv.first;
            auto& v = kv.second;
            if (v.empty()) continue;

            std::sort(v.begin(), v.end(), [](const Pafinfo& a, const Pafinfo& b){
                if (a.qstart != b.qstart) return a.qstart < b.qstart;
                if (a.qend   != b.qend)   return a.qend   < b.qend;
                return a.mapq > b.mapq;
            });

            auto& ends = paf_ends[qname];
            ends.resize(v.size());
            for (size_t i = 0; i < v.size(); ++i) ends[i] = v[i].qend;

            auto& idx = paf_idx_ends[qname];
            idx.resize(v.size());
            for (uint32_t i = 0; i < (uint32_t)v.size(); ++i) idx[i] = i;

            std::sort(idx.begin(), idx.end(), [&](uint32_t a, uint32_t b){
                if (ends[a] != ends[b]) return ends[a] < ends[b];
                if (v[a].qstart != v[b].qstart) return v[a].qstart < v[b].qstart;
                return v[a].mapq > v[b].mapq;
            });
        }
    }

    return 0;
}

std::vector<std::string> liftover(const RunOpts& opt) {
    std::vector<std::string> blocks;

    std::unordered_map<std::string, std::vector<BEDinfo>> beds;
    coordmap::CoordMap M;
    std::unordered_map<std::string, std::vector<Pafinfo>> pafs;
    std::unordered_map<std::string, std::vector<uint32_t>> paf_ends;
    std::unordered_map<std::string, std::vector<uint32_t>> paf_idx_ends;

    file_reader(opt, beds, M, pafs, paf_ends, paf_idx_ends);

    // Check mode
    if (opt.do_check) {
        check(M, opt);
        return blocks;
    }

    // Compile regex once (thread-safe for const usage)
    std::shared_ptr<const std::regex> rxp;
    if (!opt.regex.empty()) {
        rxp = std::make_shared<const std::regex>(opt.regex, std::regex::ECMAScript);
    }

    // ---------------- Thread pool ----------------
    ThreadPool pool(opt.threads);
    std::vector<std::future<std::string>> futs;
    futs.reserve(1024);

    auto submit_one =
        [&] (std::string qname, BEDinfo b)
        {
            futs.emplace_back(
                pool.submit([&, qname = std::move(qname), b = std::move(b)]() -> std::string {
                    const uint32_t qbeg = b.beg;
                    const uint32_t qend = b.end;
                    const std::string rest = b.rest;

                    std::vector<LIFTresult> results1 = paf_liftover(pafs, paf_ends, paf_idx_ends, qname, qbeg, qend);
                    std::vector<LIFTresult> results2 = map_liftover(opt, M, qname, qbeg, qend);
                    auto results  = dedup_keep_best_gap(results1, results2, qbeg, qend);

                    std::string buf;
                    buf.reserve(results.size() * 80);

                    for (const auto& r : results) {
                        if (rxp && !std::regex_search(r.tname, *rxp)) continue;

                        buf += r.tname;
                        buf += '\t';
                        buf += std::to_string(r.tbeg);
                        buf += '\t';
                        buf += std::to_string(r.tend);

                        if (!rest.empty()) {
                            buf += '\t';
                            buf += rest;
                        }

                        buf += '\t';
                        buf += r.rest;
                        buf += '\n';
                    }
                    return buf;
                })
            );
        };

    // Submit
    for (const auto& bed : beds) {
        const std::string& qname = bed.first;
        for (const auto& b : bed.second) {
            submit_one(qname, b);
        }
    }

    // Collect
    blocks.reserve(futs.size());
    ProgressTracker prog(futs.size());
    for (auto& f : futs) {
        prog.hit();
        blocks.emplace_back(f.get());
    }

    pool.stop();

    return blocks;
}

std::vector<LIFTresult> paf_liftover(
    const std::unordered_map<std::string, std::vector<Pafinfo>>& pafs,
    const std::unordered_map<std::string, std::vector<uint32_t>>& paf_ends,
    const std::unordered_map<std::string, std::vector<uint32_t>>& paf_idx_ends,
    std::string qname, uint32_t qbeg, uint32_t qend
) {
    auto find_ovlp_paf = [&](std::vector<const Pafinfo*>& out) ->bool {
        out.clear();
        auto it = pafs.find(qname);
        if (it == pafs.end() || it->second.empty()) return false;

        const auto& v    = it->second;
        const auto& ends = paf_ends.at(qname);
        const auto& idxE = paf_idx_ends.at(qname);

        size_t lo = 0, hi = idxE.size();
        while (lo < hi) {
            size_t mid = (lo + hi) >> 1;
            uint32_t id = idxE[mid];
            if (ends[id] <= qbeg) lo = mid + 1;
            else hi = mid;
        }

        for (size_t k = lo; k < idxE.size(); ++k) {
            uint32_t id = idxE[k];
            const auto& a = v[id];
            if (a.qstart < qend) out.push_back(&a);
        }
        if (out.empty()) return false;

        return true;
    };

    auto lift_by_cg = [](const Pafinfo& a, uint32_t qpos) ->int64_t {
        uint32_t qpos_oriented = qpos;
        if (a.strand == '-') {
            if (qpos >= a.qlen) return -1;
            qpos_oriented = a.qlen - qpos - 1;
        }
        
        // Start positions in oriented query coordinate
        uint32_t y = (a.strand == '+') ? a.qstart : (a.qlen - a.qend);
        uint32_t x = a.tstart;

        const uint32_t y_end = (a.strand == '+') ? a.qend : (a.qlen - a.qstart);
        if (qpos_oriented < y || qpos_oriented >= y_end) return -1;

        uint64_t num = 0;
        for (size_t i = 0; i < a.cg.size(); ++i) {
            char op = a.cg[i];
            if (op >= '0' && op <= '9') { num = num * 10 + (uint64_t)(op - '0'); continue; }
            if (num == 0) continue;

            const uint32_t len = (uint32_t)num;

            if (op == 'D') {
                x += len;
            } else if (op == 'I') {
                if (qpos_oriented >= y && qpos_oriented < y + len) return (int64_t)x;
                y += len;
            } else {
                if (qpos_oriented >= y && qpos_oriented < y + len) {
                    return (int64_t)x + (int64_t)(qpos_oriented - y);
                }
                y += len;
                x += len;
            }

            num = 0;
        }
        return -1;
    };

    std::vector<LIFTresult> results;

    std::vector<const Pafinfo*> hits;
    hits.reserve(16);

    if (pafs.find(qname) == pafs.end()) return results;

    if (qbeg >= qend) return results;

    if (!find_ovlp_paf(hits)) return results;

    for (const Pafinfo* a : hits) {
        if (qbeg >= qend || qend == 0) return results;

        uint32_t q0 = a->strand == '+' ? qbeg : qend - 1;
        uint32_t q1 = a->strand == '+' ? (qend - 1) : qbeg;

        int64_t tb0 = lift_by_cg(*a, q0);
        int64_t tb1 = lift_by_cg(*a, q1);

        std::string suffix;
        if (tb0 < 0) { suffix += "_t5"; tb0 = (int64_t)a->tstart; }
        if (tb1 < 0) { suffix += "_t3"; tb1 = (int64_t)a->tend - 1; }

        int64_t tb = tb0;
        int64_t te = tb1 + 1;

        if (tb < 0 || te < 0 || te <= tb) continue;

        if ((uint64_t)tb > (uint64_t)a->tlen) tb = (int64_t)a->tlen;
        if ((uint64_t)te > (uint64_t)a->tlen) te = (int64_t)a->tlen;
        std::string rest = qname + "_" + std::to_string(qbeg) + "_" + std::to_string(qend);
        rest += suffix;
        rest += "\t" + std::string(a->strand == '+' ? "+" : "-");

        results.emplace_back(a->tname, tb, te, qbeg, qend, a->strand, rest);
    }

    return results;
}

std::vector<LIFTresult> map_liftover(
    const RunOpts& opt, 
    const coordmap::CoordMap& M, 
    std::string qname, uint32_t qbeg, uint32_t qend
) {
    auto keep = [&](const coordmap::Hit& h){
        uint32_t hlen = (h.qend > h.qbeg) ? (h.qend - h.qbeg) : (h.qbeg - h.qend);
        uint32_t  len = (qend > qbeg) ? (qend - qbeg) : (qbeg - qend);
        return (len > 0 ? double(hlen) / double(len) : 0) >= opt.min_frac;
    };

    enum class Side : uint8_t { NORMAL=0, LEFT=1, RIGHT=2 };
    auto query_gap = [&](const coordmap::Hit& h, Side side)->int64_t {
        // NORMAL: gap = |h.qbeg - qbeg| + |h.qend - qend|
        //      (e.g. qbeg --- h.qbeg --- h.qend --- qend)
        //   LEFT: gap = |qbeg - h.qend|
        //      (e.g. h.qbeg --- h.qend --- qbeg --- qend)
        //  RIGHT: gap = |h.qbeg - qend|
        //      (e.g. qbeg --- qend --- h.qbeg --- h.qend)

        if (side == Side::LEFT)  return std::abs((int64_t)qbeg - (int64_t)h.qend);
        if (side == Side::RIGHT) return std::abs((int64_t)h.qbeg - (int64_t)qend);

        int64_t g = 0;
        g += std::abs((int64_t)h.qbeg - (int64_t)qbeg);
        g += std::abs((int64_t)qend - (int64_t)h.qend);
        return g;
    };

    auto take_top_hits = [&](uint32_t wb, uint32_t we, size_t K, Side side) {
        std::vector<coordmap::Hit> hits = M.map_range(
            qname, wb, we,
            opt.cm_max_hops, opt.cm_max_fanout, opt.cm_min_len, opt.cm_min_frac, opt.cm_max_total_hits,
            true
        );

        // pair.first.first=length, pair.first.second=query_gap, pair.second=hit*
        std::vector<std::pair<std::pair<uint32_t,int64_t>, coordmap::Hit*>> scored;
        scored.reserve(hits.size());

        for (auto& h : hits) {
            uint32_t len = (h.qend > h.qbeg) ? (h.qend - h.qbeg) : h.qbeg - h.qend;
            int64_t  gap = query_gap(h, side);
            scored.push_back({{len, gap}, &h});
        }

        auto better = [](const auto& a, const auto& b) {
            if (a.first.second != b.first.second) return a.first.second < b.first.second;
            return a.first.first > b.first.first;
        };

        if (scored.size() > K) {
            std::nth_element(scored.begin(), scored.begin() + K, scored.end(), better);
            scored.resize(K);
        }
        std::sort(scored.begin(), scored.end(), better);

        // pair.first=query_gap, pair.second=hit
        std::vector<std::pair<int64_t, coordmap::Hit>> out;
        out.reserve(scored.size());
        for (auto& x : scored) out.push_back({x.first.second, *x.second});
        return out;
    };

    std::vector<LIFTresult> results;

    if (qend == 0) return results;
    
    if (qbeg >= qend) return results;

    /* -------- normal liftover -------- */
    {
        std::vector<coordmap::Hit> hits = M.map_range(qname, qbeg, qend, opt.cm_max_hops, opt.cm_max_fanout, opt.cm_min_len, opt.cm_min_frac, opt.cm_max_total_hits, true);
        std::sort(hits.begin(), hits.end(), [&](const coordmap::Hit& a, const coordmap::Hit& b){ return query_gap(a, Side::NORMAL) < query_gap(b, Side::NORMAL); });

        for (auto& h : hits) {
            if (!keep(h)) continue;

            std::string suffix;
            if (h.qbeg > qbeg) { suffix += "_t5"; }
            if (h.qend < qend) { suffix += "_t3"; }

            if (h.beg == h.end || (!h.rev && h.qend < qend) || (h.rev && qbeg < h.qbeg)) h.end += 1;

            std::string rest = qname + "_" + std::to_string(qbeg) + "_" + std::to_string(qend);
            rest += suffix;
            rest += "\t" + std::string(h.rev ? "-" : "+");
            rest += "\t" + qname + ':' + std::to_string(h.qbeg) + '-' + std::to_string(h.qend);
            rest += "\t" + M.contig_name(h.ctg) + ':' + std::to_string(h.beg) + '-' + std::to_string(h.end);
            
            results.emplace_back(M.contig_name(h.ctg), h.beg, h.end, h.qbeg, h.qend, (h.rev ? '-' : '+'), rest);
        }
    }

    if (opt.max_flank == 0) return results;

    /* -------- extend liftover -------- */
    std::vector<LIFTresult> results_ext;
    {
        // pair.first=query_gap, pair.second=hit
        std::vector<std::pair<int64_t, coordmap::Hit>> Ltop, Rtop;

        uint32_t Loff=0, Roff=0;
        while (Ltop.empty() && Loff < opt.max_flank) {
            uint32_t l0 = (qbeg > Loff + opt.flank_win) ? (qbeg - Loff - opt.flank_win) : 0;
            uint32_t l1 = (qbeg > Loff) ? (qbeg - Loff) : 0;
            Loff += opt.flank_win;
            if (l0 >= l1) break;
            Ltop = take_top_hits(l0, l1, opt.max_hit, Side::LEFT);
        }
        while (Rtop.empty() && Roff < opt.max_flank) {
            uint32_t r0 = qend + Roff;
            uint32_t r1 = r0 + opt.flank_win;
            Roff += opt.flank_win;
            if (r0 >= r1) break;
            Rtop = take_top_hits(r0, r1, opt.max_hit, Side::RIGHT);
        }

        if (Ltop.empty() || Rtop.empty()) return results;

        // // print
        // for (size_t i = 0; i < Ltop.size(); ++i) {
        //     const auto& p = Ltop[i];
        //     const auto& h = p.second;

        //     std::cerr
        //         << "[L " << i << "] gap=" << p.first << "\n"
        //         << "  ctg=" << h.ctg << " strand=" << (h.rev ? '-' : '+')
        //         << "  T:[" << h.beg << "," << h.end << ")"
        //         << "  Q:[" << h.qbeg << "," << h.qend << ")\n";
        // }
        // for (size_t i = 0; i < Rtop.size(); ++i) {
        //     const auto& p = Rtop[i];
        //     const auto& h = p.second;

        //     std::cerr
        //         << "[R " << i << "] gap=" << p.first << "\n"
        //         << "  ctg=" << h.ctg << " strand=" << (h.rev ? '-' : '+')
        //         << "  T:[" << h.beg << "," << h.end << ")"
        //         << "  Q:[" << h.qbeg << "," << h.qend << ")\n";
        // }

        struct Pair { int64_t score; coordmap::Hit L, R; };
        std::vector<Pair> pairs;
        pairs.reserve(Ltop.size() * Rtop.size());

        for (auto& L : Ltop) for (auto& R : Rtop) {
            auto& l = L.second; auto& r = R.second;
            if (l.ctg != r.ctg || l.rev != r.rev) continue;

            int64_t tgap = l.rev ? (int64_t)l.beg - (int64_t)r.end : (int64_t)r.beg - (int64_t)l.end;
            int64_t qgap = (int64_t)r.qbeg - (int64_t)l.qend;
            if (tgap < 0 || qgap < 0) continue;
            int64_t gap  = (tgap > qgap) ? (tgap - qgap) : (qgap - tgap);
            if (gap> opt.max_gap) continue;

            // int64_t score = gap * 1000000LL + (L.first + R.first);
            int64_t score = L.first + R.first;
            pairs.push_back(Pair{score, l, r});
        }

        std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b){ return a.score < b.score; });

        // // print pairs information
        // for (size_t i = 0; i < pairs.size(); ++i) {
        //     const auto& p = pairs[i];
        //     const auto& l = p.L;
        //     const auto& r = p.R;

        //     std::cerr
        //         << "[pair " << i << "] score=" << p.score << "\n";

        //     std::cerr
        //         << "  [L] ctg=" << l.ctg << " strand=" << (l.rev ? '-' : '+')
        //         << "  T:[" << l.beg << "," << l.end << ")"
        //         << "  Q:[" << l.qbeg << "," << l.qend << ")\n";

        //     std::cerr
        //         << "  [R] ctg=" << r.ctg << " strand=" << (r.rev ? '-' : '+')
        //         << "  T:[" << r.beg << "," << r.end << ")"
        //         << "  Q:[" << r.qbeg << "," << r.qend << ")\n";
        // }

        size_t takeN = std::min<size_t>(pairs.size(), opt.max_hit);
        for (size_t i = 0; i < takeN; ++i) {
            auto& bestL = pairs[i].L;
            auto& bestR = pairs[i].R;

            int64_t t_beg = bestL.rev ? bestR.end : bestL.end;
            int64_t t_end = bestL.rev ? bestL.beg : bestR.beg;
            int64_t q_beg = bestL.qend;
            int64_t q_end = bestR.qbeg;

            // if (t_beg == t_end) t_end += 1;
            // t_end += 1;
            if ((t_beg == t_end) || (!bestL.rev && q_end < qend)) t_end += 1;
            if (t_beg >= t_end) continue;

            // if (q_beg == q_end) q_end += 1;
            q_end += 1;

            std::string rest = qname + "_" + std::to_string(qbeg) + "_" + std::to_string(qend);
            rest += "_ext";
            rest += "\t" + std::string(bestL.rev ? "-" : "+");
            rest += "\t" + qname + ':' + std::to_string(bestL.qbeg) + '-' + std::to_string(bestL.qend) + '|' + std::to_string(bestR.qbeg) + '-' + std::to_string(bestR.qend);
            rest += "\t" + M.contig_name(bestL.ctg) + ':' + std::to_string(bestL.beg) + '-' + std::to_string(bestL.end) + '|' + std::to_string(bestR.beg) + '-' + std::to_string(bestR.end);

            LIFTresult r(M.contig_name(bestL.ctg), (uint32_t)t_beg, (uint32_t)t_end, (uint32_t)q_beg, (uint32_t)q_end, bestL.rev ? '-' : '+', rest);

            results_ext.push_back(std::move(r));
        }
    }

    results = dedup_keep_best_gap(results, results_ext, qbeg, qend);

    return results;
}

std::vector<LIFTresult> dedup_keep_best_gap(
    std::vector<LIFTresult> fixed,
    std::vector<LIFTresult> cand,
    uint32_t qbeg, uint32_t qend
) {
    auto ovlp = [](const LIFTresult& a, const LIFTresult& b){
        return a.tname == b.tname && a.strand == b.strand && !(a.tend <= b.tbeg || b.tend <= a.tbeg);
    };

    // auto qgap = [&](const LIFTresult& h)->int64_t {
    //     int64_t g = 0;
    //     g += (h.qbeg < qbeg) ? (int64_t)(qbeg - h.qbeg) : (int64_t)(h.qbeg - qbeg);
    //     g += (h.qend < qend) ? (int64_t)(qend - h.qend) : (int64_t)(h.qend - qend);
    //     return g;
    // };

    // std::stable_sort(cand.begin(), cand.end(), [&](const LIFTresult& a, const LIFTresult& b){
    //     if (a.tname != b.tname) return a.tname < b.tname;
    //     if (a.strand != b.strand) return a.strand < b.strand;
    //     auto ga = qgap(a), gb = qgap(b);
    //     if (ga != gb) return ga < gb;
    //     auto la = (a.tend > a.tbeg) ? (a.tend - a.tbeg) : 0u;
    //     auto lb = (b.tend > b.tbeg) ? (b.tend - b.tbeg) : 0u;
    //     if (la != lb) return la > lb;
    //     return a.tbeg < b.tbeg;
    // });

    std::vector<LIFTresult> kept;
    kept.reserve(fixed.size() + cand.size());

    for (auto& r : fixed) kept.push_back(std::move(r));

    for (auto& r : cand) {
        bool bad = false;
        for (const auto& k : kept) {
            if (ovlp(k, r)) { bad = true; break; }
        }
        if (!bad) kept.push_back(std::move(r));
    }

    return kept;
}

void save_liftover_results(const std::string out_file, const std::vector<std::string>& blocks)
{
    std::ostream* out = &std::cout;
    std::ofstream fout;

    if (!out_file.empty() && out_file != "-") {
        fout.open(out_file);
        if (!fout) {
            error_stream() << out_file << ": No such file or directory";
            std::exit(1);
        }
        out = &fout;
    }

    for (const auto& s : blocks) {
        (*out) << s;
    }
}

// Function check
void check(const coordmap::CoordMap& M, const RunOpts& opt) {
    struct CheckStats {
        uint64_t windows_total = 0, windows_hit = 0;
        uint64_t mapped_bases = 0, equal_bases = 0, diff_bases = 0;
        uint64_t len_mismatch = 0, bounds_error = 0;
        uint32_t examples = 0;
    };

    SeqDB db;
    if (!db.load(opt.ref_file)) {
        error_stream() << "Load ref failed: " + opt.ref_file << std::endl;
        std::exit(1);
    }

    CheckStats st;
    log_stream() << "Win=" << opt.win << " step=" << opt.step << "\n";

    for (const auto& ctg : db.names) {
        const std::string* src_seq = db.get(ctg);
        if (!src_seq) continue;
        const uint32_t L = (uint32_t)src_seq->size();
        if (!L) continue;

        for (uint32_t qbeg = 0; qbeg < L; qbeg += opt.step) {
            uint32_t qend = std::min<uint32_t>(qbeg + opt.win, L);
            if (qbeg >= qend) continue;
            ++st.windows_total;

            auto hits = M.map_range(ctg, qbeg, qend, opt.cm_max_hops, opt.cm_max_fanout, opt.cm_min_len, opt.cm_min_frac, opt.cm_max_total_hits, false);
            if (hits.empty()) continue;
            ++st.windows_hit;

            for (const auto& h : hits) {
                const std::string tname = M.contig_name(h.ctg);
                const std::string* dst_seq = db.get(tname);

                const uint32_t qlen = h.qend - h.qbeg;
                const uint32_t dlen = h.end  - h.beg;

                if (!dst_seq || h.qend > src_seq->size() || h.end > dst_seq->size()) {
                    ++st.bounds_error;
                    if (st.examples < opt.max_examples) {
                        warning_stream() << "[BOUNDS] " << ctg << ":" << h.qbeg << "-" << h.qend
                                         << " -> " << tname << ":" << h.beg << "-" << h.end
                                         << (h.rev ? " (-)\n" : " (+)\n");
                        ++st.examples;
                    }
                    continue;
                }
                if (!qlen || qlen != dlen) {
                    ++st.len_mismatch;
                    if (st.examples < opt.max_examples) {
                        warning_stream() << "[LEN] " << ctg << ":" << h.qbeg << "-" << h.qend
                                         << " -> " << tname << ":" << h.beg << "-" << h.end
                                         << " qlen=" << qlen << " dlen=" << dlen
                                         << (h.rev ? " (-)\n" : " (+)\n");
                        ++st.examples;
                    }
                    continue;
                }

                std::string src = src_seq->substr(h.qbeg, qlen);
                std::string dst = dst_seq->substr(h.beg,  dlen);
                if (h.rev) detail::rc_inplace(dst);

                st.mapped_bases += qlen;

                uint32_t diffs = 0, first = qlen;
                for (uint32_t i = 0; i < qlen; ++i) {
                    if (src[i] != dst[i]) { ++diffs; if (first == qlen) first = i; }
                }
                st.equal_bases += (qlen - diffs);
                st.diff_bases  += diffs;

                if (diffs && st.examples < opt.max_examples) {
                    const uint32_t ctx0 = (first > 10 ? first - 10 : 0);
                    const uint32_t ctx1 = std::min<uint32_t>(first + 10, qlen);
                    warning_stream() << "[DIFF] " << ctg << ":" << h.qbeg << "-" << h.qend
                                     << " -> " << tname << ":" << h.beg << "-" << h.end
                                     << (h.rev ? " (-)" : " (+)")
                                     << " first=" << first << " diffs=" << diffs << "\n"
                                     << "       SRC(" << ctx0 << "-" << ctx1 << "): " << src.substr(ctx0, ctx1 - ctx0) << "\n"
                                     << "       DST(" << ctx0 << "-" << ctx1 << "): " << dst.substr(ctx0, ctx1 - ctx0) << "\n";
                    ++st.examples;
                }
            }
        }
    }

    log_stream()
        << "[check] windows=" << st.windows_total
        << " hit=" << st.windows_hit
        << " mapped=" << st.mapped_bases
        << " equal=" << st.equal_bases
        << " diff=" << st.diff_bases
        << " len_mis=" << st.len_mismatch
        << " bounds=" << st.bounds_error
        << "\n";

    return;
}

} // namespace liftover