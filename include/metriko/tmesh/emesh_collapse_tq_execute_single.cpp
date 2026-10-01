//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include <algorithm>
#include <format>
#include <set>
#include "emesh.h"
using namespace metriko;

void Emesh::collapse_tquad_execute(Tqaux& aux) {
    const auto& th_l = thalfs[aux.thid_l];
    const auto& th_r = thalfs[aux.thid_r];

    auto lateral_run = [&](int fr, int to) -> vec<int> {
        for (auto [a, b]: {std::pair{fr, to}, std::pair{to, fr}}) {
        for (auto* s: { &aux.thids_t, &aux.thids_b }) {
            vec<int> res;
            int c = a;
            while (c != b)
                if (auto it = rg::find_if(*s, [&](int t) { return thalfs[t].nid_fr() == c; }); it != s->end()) {
                    res.push_back(*it);
                    c = thalfs[*it].nid_to();
                }
            if (c == b) { if (a != fr) rg::reverse(res); return res; }
        }}
        METRIKO_FAIL("no thalf run on a lateral side from node {} to node {}", fr, to);
    };

    auto lateral_btn = [&](int n0, int n1) -> int {
        for (auto* s: { &aux.thids_t, &aux.thids_b })
        for (int thid: *s) {
            auto& th = thalfs[thid];
            if ((th.nid_fr() == n0 && th.nid_to() == n1) || (th.nid_fr() == n1 && th.nid_to() == n0)) return thid;
        }
        METRIKO_FAIL("no lateral thalf between nodes {} and {}", n0, n1);
    };

    // the band closes into a loop around a tube: a shortest path between two points would jump the seam, so the
    // merged line follows the band instead, along n0's side to its far end and over the zero side to n1. empty
    // when n0 is not at the end of its side (an inner merge) or the walk does not reach n1
    auto along_band = [&](int n0, int n1) -> vec<HmLoc> {
        const auto& tz  = thalfs[rg::contains(std::array{ th_l.nid_fr(), th_l.nid_to() }, n1) ? aux.thid_l : aux.thid_r];
        const int   end = tz.nid_fr() == n1 ? tz.nid_to() : tz.nid_fr();
        vec<int> thids;
        try { thids = lateral_run(n0, end); } catch (const std::runtime_error&) { return {}; }
        thids.push_back(tz.id);
        vec<int> nids = { n0 };
        for (int t: thids) {
            auto ns = tedges[thalfs[t].teid].nids;
            if (ns.back() == nids.back()) rg::reverse(ns);
            if (ns.front() != nids.back()) return {};
            nids.insert(nids.end(), ns.begin() + 1, ns.end());
        }
        if (nids.back() != n1) return {};
        vec<HmLoc> path;
        for (int n: nids) path.push_back(tnodes[n]);
        return path;
    };

    const auto region = allowed_range_tquads({ aux.tqid });
    const bool tube   = [&] {
        std::set<int> fs, vs, es;
        for (const auto& [eid, r0, r1]: region) { fs.insert(hm.edges[eid].face0().id); fs.insert(hm.edges[eid].face1().id); }
        for (int fid: fs) for (Half h: hm.faces[fid].adjHalfs()) { vs.insert(h.tail().id); es.insert(h.edge().id); }
        return (int)vs.size() - (int)es.size() + (int)fs.size() == 0;
    }();

    // 1: the pieces of the merged line between consecutive points: a lateral thalf when both points lie on the
    // same side, a new tedge across the band when the line switches sides
    struct Piece { int thid; bool t_fr; bool t_to; };
    vec<Piece> pieces;
    for (size_t i = 0; i + 1 < aux.pts.size(); ++i) {
        const auto& [n0, v0, ord0, s0] = aux.pts[i];
        const auto& [n1, v1, ord1, s1] = aux.pts[i + 1];
        if (s0 == s1) { pieces.emplace_back(lateral_btn(n0, n1), s0, s1); continue; }

        auto path = tube ? along_band(n0, n1) : vec<HmLoc>{};
        if (path.empty()) path = approx_shortest_path(30, hm, tnodes[n0], tnodes[n1], region);
        METRIKO_CHECK(path.size() >= 2, "no path within the allowed region: tqid {}", aux.tqid);

        const auto   nids  = add_new_path(path, n0, n1);
        const int    teid  = tedges.size();
        const int    thid0 = thalfs.size();
        const int    thid1 = thalfs.size() + 1;
        const double x     = std::abs(v1 - v0);
        const double r     = path_length(nids);
        tedges.emplace_back(teid, nids);
        thalfs.emplace_back(this, thid0, thid1, teid, true,  x, r);
        thalfs.emplace_back(this, thid1, thid0, teid, false, x, r);
        pieces.emplace_back(thid0, s0, s1);
        pieces.emplace_back(thid1, s1, s0);
    }

    // 2: group the pieces into runs, each replacing the thalfs of one outer tquad along the band. a run walks
    // from a node piece by piece and stops once it reaches the side it was told to stop at
    auto take_run = [&](int nid, bool side_to_stop, bool backwards) {
        vec<int> res;
        int cur = nid;
        while (true) {
            auto it = rg::find_if(pieces, [&](const Piece& p) { return (backwards ? thalfs[p.thid].nid_to() : thalfs[p.thid].nid_fr()) == cur; });
            if (it == pieces.end()) break;
            const bool stop = backwards ? it->t_fr == side_to_stop : it->t_to == side_to_stop;
            res.push_back(it->thid);
            cur = backwards ? thalfs[it->thid].nid_fr() : thalfs[it->thid].nid_to();
            pieces.erase(it);
            if (stop) break;
        }
        if (backwards) rg::reverse(res);
        return res;
    };

    const auto [nid_bgn, v_bgn, o_bgn, top_bgn] = aux.pts.front();
    const auto [nid_end, v_end, o_end, top_end] = aux.pts.back();
    vec<int> run_bgn = take_run(nid_bgn, !top_bgn,  top_bgn);
    vec<int> run_end = take_run(nid_end, !top_end, !top_end);
    vec<vec<int>> runs;
    auto is_head = [&](const Piece& q) { return rg::none_of(pieces, [&](const Piece& p) { return thalfs[p.thid].nid_to() == thalfs[q.thid].nid_fr(); }); };
    auto is_tail = [&](const Piece& q) { return rg::none_of(pieces, [&](const Piece& p) { return thalfs[p.thid].nid_fr() == thalfs[q.thid].nid_to(); }); };
    while (!pieces.empty()) {
        if (auto it = rg::find_if(pieces, is_head); it != pieces.end()) { runs.push_back(take_run(thalfs[it->thid].nid_fr(), it->t_fr, false)); continue; }
        if (auto it = rg::find_if(pieces, is_tail); it != pieces.end()) { runs.push_back(take_run(thalfs[it->thid].nid_to(), it->t_to, true));  continue; }
        METRIKO_FAIL("pieces of the merged line on one side only");
    }

    // 3: the outer tquads: each old thalf across a lateral side is replaced by the part of a run that spans it
    auto replace = [&](const vec<int>& olds, const vec<int>& run) {
        size_t ri = 0;
        for (int old: olds) {
            const int e0 = thalfs[old].nid_fr();
            const int e1 = thalfs[old].nid_to();
            vec<int> seg;
            while (ri < run.size()) {
                seg.push_back(run[ri]);
                const int to = thalfs[run[ri]].nid_to(); ++ri;
                if (ri >= run.size() || to == e0 || to == e1) break;
            }
            auto& [id, data] = tquads[thalfs[old].tqid];
            auto it = rg::find(data, old, &Edata::thid);
            METRIKO_CHECK(it != data.end(), "thalf {} not on its tquad", old);
            const auto si = it->side;
            it = data.erase(it);
            vec<Edata> repl;
            for (int t: seg) { thalfs[t].tqid = id; repl.push_back({ t, si }); }
            data.insert(it, repl.begin(), repl.end());
        }
    };
    auto outer = [&](int n0, int n1) { return lateral_run(n0, n1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>(); };

    // 4: the zero thalfs at the two ends. absorbed: the twin leaves its tquad and the tedge continuing past the
    // corner is extended along the zero path, so the corner moves to the merged line. kept: the zero thalf is
    // handed to the tquad across the corner as a zero-length thalf, to be dealt with later, when the corner is a
    // full junction or when absorbing would leave the tquad across the zero side without a side
    auto extend = [&](const Ehalf& th, bool ahd, bool keep) {
        auto step = [&](int thid) { return ahd ? step_next(thid) : step_prev(thid); };
        auto& th1 = thalfs[step(th.id)];
        auto& th2 = thalfs[step(th1.twid)];
        auto& [id, data] = tquads[th2.tqid];
        if (keep) {
            auto it = rg::find(data, th2.id, &Edata::thid);
            const auto si = it->side;
            thalfs[th.id].tqid = id;
            data.insert(ahd ? it : it + 1, Edata{ th.id, si });
            return;
        }
        auto& nids    = tedges[th.teid].nids;
        auto& th_twn  = thalfs[th.twid];
        auto& th2_twn = thalfs[th2.twid];
        std::erase_if(tquads[th_twn.tqid].data, [&](const Edata& d) { return d.thid == th_twn.id; });

        // the tquad across the corner also lies across the zero side, so the thalf following th1 is th's own
        // twin: there is no tedge to extend, dropping the twin is all that is needed
        if (th2.teid == th.teid) return;
        tedges[th2.teid].insert_locs(nids);

        // the zero thalf closed th2 into a loop (a band around a tube): both attachments were valid, so make
        // sure the loop's endpoint is the node where the other tedges meet, not the node being absorbed
        if (auto& loop = tedges[th2.teid].nids; loop.front() == loop.back()) {
            auto deg = [&](int n) {
                int c = 0;
                for (const auto& [teid, ns]: live_tedges()) if (teid != th2.teid && (ns.front() == n || ns.back() == n)) ++c;
                return c;
            };
            const int crnr = deg(nids.front()) >= deg(nids.back()) ? nids.front() : nids.back();
            if (loop.front() != crnr) {
                loop.pop_back();
                rg::rotate(loop, rg::find(loop, crnr));
                loop.push_back(crnr);
            }
        }

        if      (th.bgn)     { (th2.nid_fr() == th.nid_fr()     ? th2 : th2_twn).bgn = true; }
        else if (th_twn.bgn) { (th2.nid_fr() == th_twn.nid_fr() ? th2 : th2_twn).bgn = true; }
        else if (th.end)     { (th2.nid_to() == th.nid_to()     ? th2 : th2_twn).end = true; }
        else if (th_twn.end) { (th2.nid_to() == th_twn.nid_to() ? th2 : th2_twn).end = true; }
    };

    auto on_line = [&](int nid) { return rg::any_of(aux.pts, [&](const Tqpoint& p) { return p.nid == nid; }); };
    // absorbing a terminal drops its twin from the tquad across the zero side, which survives and must stay a valid
    // quad. both terminals are absorbed here, so both twins are taken away before checking
    auto absorb_breaks_across = [&](const Ehalf& th) {
        Equad copy = tquads[thalfs[th.twid].tqid];
        std::erase_if(copy.data, [&](const Edata& d) { return d.thid == th_l.twid || d.thid == th_r.twid; });
        return !copy.is_valid();
    };
    auto keep_of = [&](const Ehalf& th) {
        bool full_fr = count_adj_tquads(th.id)   == 4 && !on_line(th.nid_fr());
        bool full_to = count_adj_tquads(th.twid) == 4 && !on_line(th.nid_to());
        return full_fr || full_to || absorb_breaks_across(th);
    };
    bool keep_l = keep_of(th_l);
    bool keep_r = keep_of(th_r);

    if (!run_bgn.empty() && run_end.empty()) {
        // the whole line was taken from the left end: one outer stretch spans the band, both ends extend it
        METRIKO_CHECK(top_bgn == top_end, "the merged line ends on different sides");
        METRIKO_CHECK(run_bgn.size() == aux.pts.size() - 1, "the merged line is incomplete");
        bool ahd_l = th_l.nid_fr() == nid_bgn;
        bool ahd_r = th_r.nid_fr() == nid_end;
        extend(th_l, ahd_l, keep_l);
        extend(th_r, ahd_r, keep_r);
        replace(outer(ahd_l ? th_l.nid_to() : th_l.nid_fr(), ahd_r ? th_r.nid_to() : th_r.nid_fr()), run_bgn);
    } else {
        auto end_run = [&](const Ehalf& th, int nid_term, bool keep, const vec<int>& run) {
            bool ahd = th.nid_fr() == nid_term;
            extend(th, ahd, keep);
            int n0 = ahd ? th.nid_to() : th.nid_fr();
            int na = thalfs[run.back()].nid_to();
            int nb = thalfs[run.front()].nid_fr();
            int n1 = na == th.nid_fr() || na == th.nid_to() ? nb : na;
            replace(outer(n0, n1), run);
        };
        end_run(th_l, nid_bgn, keep_l, run_bgn);
        end_run(th_r, nid_end, keep_r, run_end);
        for (const auto& run: runs) {
            METRIKO_CHECK(run.size() >= 2, "a middle run of the merged line has a single piece");
            replace(outer(thalfs[run.front()].nid_fr(), thalfs[run.back()].nid_to()), run);
        }
    }

    // 5: the band itself is gone
    auto& [id, data] = tquads[aux.tqid];
    for (const auto& [thid, _]: data) {
        auto& th0 = thalfs[thid];
        if (th0.tqid != id) continue;
        auto& th1 = thalfs[th0.twid];
        auto& te  = tedges[th0.teid];
        th0 = {}; th1 = {}; te = {};
    }
    data.clear();
    id = -1;
}
