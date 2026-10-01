//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include <set>
#include "emesh.h"
using namespace metriko;

void Emesh::collapse_thalf(int thid) {
    auto& th_crr = thalfs[thid];
    auto& th_twn = thalfs[th_crr.twid];
    auto& tq_crr = tquads[th_crr.tqid];
    auto& tq_twn = tquads[th_twn.tqid];
    auto  it_crr = rg::find(tq_crr.data, thid     , &Edata::thid);
    auto  it_twn = rg::find(tq_twn.data, th_twn.id, &Edata::thid);
    auto  it_prv = circular_prev(tq_crr.data, it_crr);
    auto  it_nxt = circular_next(tq_crr.data, it_crr);
    auto& th_prv = thalfs[it_prv->thid];
    auto& th_nxt = thalfs[it_nxt->thid];
    auto& te_crr = tedges[th_crr.teid];
    auto& te_prv = tedges[th_prv.teid];
    auto& te_nxt = tedges[th_nxt.teid];

    METRIKO_CHECK(it_crr->side == it_prv->side || it_crr->side == it_nxt->side, "collapse thalf error");
    METRIKO_CHECK(it_crr->side != it_prv->side || it_crr->side != it_nxt->side, "collapse thalf error");

    int cout_fr = count_adj_tquads(th_crr.id);
    int cout_to = count_adj_tquads(th_twn.id);

    if      (cout_fr == 2) { te_prv.insert_locs(te_crr.nids); }
    else if (cout_to == 2) { te_nxt.insert_locs(te_crr.nids); }
    else {
        bool collapse_to_prev = it_crr->side != it_prv->side;

        auto [n_fr, n_to] = [&]{
            if ( collapse_to_prev &&  th_prv.cano) return std::pair{th_prv.nid_fr(), th_crr.nid_to()};
            if ( collapse_to_prev && !th_prv.cano) return std::pair{th_crr.nid_to(), th_prv.nid_fr()};
            if (!collapse_to_prev &&  th_nxt.cano) return std::pair{th_crr.nid_fr(), th_nxt.nid_to()};
            if (!collapse_to_prev && !th_nxt.cano) return std::pair{th_nxt.nid_to(), th_crr.nid_fr()};
            throw std::runtime_error("unreachable");
        }();

        auto region = allowed_range_tquads({tq_crr.id});
        auto path   = approx_shortest_path(30, hm, tnodes[n_fr], tnodes[n_to], region);
        auto nids   = add_new_path(path, n_fr, n_to);

        { // Count duplication of verts for checking loc-inj
            vec<HmLoc> locs;
            for (const auto& [id, side]: tq_crr.data) {
                auto ns = tedges[thalfs[id].teid].nids;
                if (!thalfs[id].cano) rg::reverse(ns);
                for (size_t k = 0; k + 1 < ns.size(); ++k) locs.push_back(tnodes[ns[k]]);
            }

            int dups = 0;
            for (size_t i = 0; i < locs.size(); ++i)
            for (size_t j = 0; j < i; ++j)
                if (locs[i] == locs[j]) { dups++; break; }

            // if dups >= 2, tquad is self intersected by one of its thalfs
            // if dups == 1, tquad is self intersected by one of its singular (now tempolary skip, and hope aother tquad is loc-inj)
            METRIKO_CHECK(dups < 2, "tquad {} is self-intersected by one of its thalfs", tq_crr.id);
            if (dups == 1) {
                auto& te_tgt = collapse_to_prev ? te_prv : te_nxt;
                if (path_length(nids) < 0.5 * (path_length(te_crr.nids) + path_length(te_tgt.nids))) return;
            }
        }

        auto  it_twn_adj = collapse_to_prev ? circular_next(tq_twn.data, it_twn) : circular_prev(tq_twn.data, it_twn);
        auto& th_twn_adj = thalfs[it_twn_adj->thid];
        auto& te_twn_adj = tedges[th_twn_adj.teid];
        METRIKO_CHECK(it_twn_adj->side == it_twn->side, "collapse thalf error");

        // 1: collapse to the prv/nxt edge
        // 2: insert missing segments to the adjacent edge.
        if (collapse_to_prev) { te_prv.nids = nids; te_twn_adj.insert_locs(te_crr.nids); }
        else                  { te_nxt.nids = nids; te_twn_adj.insert_locs(te_crr.nids); }
    }

    // remove the data from tquads which have collapsed thalfs
    tq_crr.data.erase(it_crr);
    tq_twn.data.erase(it_twn);
    th_crr = {};
    th_twn = {};
    te_crr = {};
}

// a point tquad: every side is a single zero thalf, so its four corners are one point of the integer grid. contract
// them into one tnode: the tedges ending at the other corners are extended to it (along the zero tedge to an adjacent
// corner, through the tquad from the opposite one), then the four zero tedges and the tquad are removed. every far
// side must keep another thalf (a zero one is fine: the side stays the ladder of a band that collapses afterwards),
// and at most one other tedge may end at each moved corner, since more would be extended along the same path
bool Emesh::collapse_point_tquad(int tqid) {
    auto& tq = tquads[tqid];
    if (tq.id == -1 || tq.data.size() != 4) return false;

    vec<int> zs;
    for (const Edata& d: tq.data) zs.push_back(d.thid);
    for (int k = 0; k < 4; ++k) {
        const auto& th = thalfs[zs[k]];
        if (th.x != 0 || th.nid_to() != thalfs[zs[(k + 1) % 4]].nid_fr()) return false;
        const auto& tw   = thalfs[th.twid];
        const auto& tq_f = tquads[tw.tqid];
        if (tq_f.thids(tq_f.side_of(tw)).size() < 2) return false;
    }

    vec<int> cs;
    std::set<int> zteids;
    for (int z: zs) { cs.push_back(thalfs[z].nid_fr()); zteids.insert(thalfs[z].teid); }
    auto ends_at = [&](int nid) {
        vec<int> res;
        for (const auto& te: live_tedges()) if (!zteids.contains(te.id) && (te.nids.front() == nid || te.nids.back() == nid)) res.push_back(te.id);
        return res;
    };
    vec<vec<int>> inc;
    for (int c: cs) inc.push_back(ends_at(c));
    const int keep = (int)std::distance(inc.begin(), rg::max_element(inc, {}, &vec<int>::size));
    for (int i = 0; i < 4; ++i) if (i != keep && inc[i].size() > 1) return false;

    for (int i = 0; i < 4; ++i) {
        if (i == keep || inc[i].empty()) continue;
        vec<int> route;
        if      ((i + 1) % 4 == keep) route = tedges[thalfs[zs[i]].teid].nids;    // zs[i] runs cs[i] -> cs[keep]
        else if ((keep + 1) % 4 == i) route = tedges[thalfs[zs[keep]].teid].nids; // zs[keep] runs cs[keep] -> cs[i]
        else route = add_new_path(approx_shortest_path(30, hm, tnodes[cs[i]], tnodes[cs[keep]], allowed_range_tquads({tqid})), cs[i], cs[keep]);
        METRIKO_CHECK(route.size() >= 2, "point tquad {}: no path from corner {} to corner {}", tqid, cs[i], cs[keep]);
        auto& te = tedges[inc[i].front()];
        te.insert_locs(route);
        const double r = path_length(te.nids);
        for (auto& th: thalfs) if (th.id != -1 && th.teid == te.id) th.r = r;
    }

    for (int z: zs) {
        const int twid = thalfs[z].twid;
        const int teid = thalfs[z].teid;
        std::erase_if(tquads[thalfs[twid].tqid].data, [&](const Edata& d) { return d.thid == twid; });
        thalfs[twid] = {};
        thalfs[z]    = {};
        tedges[teid] = {};
    }
    tq.data.clear();
    tq.id = -1;
    return true;
}
