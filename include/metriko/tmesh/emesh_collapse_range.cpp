//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include <set>
#include "emesh.h"

namespace metriko {
vec<std::tuple<int, double, double>> Emesh::allowed_range_thalfs(const vec<int>& thids) const {
    vec v_stop(hm.nV, false);
    umap<int, double> lo, hi;
    std::set<int> v_wall = {};

    for (int thid : thids) {
        const auto& th = thalfs[thid];
        const auto& te = tedges[th.teid];
        int n = te.nids.size();
        for (int i = 0; i < n - 1; ++i) {
            auto  j  = th.cano ? i : n - 1 - i;
            auto  k  = th.cano ? j + 1 : j - 1;
            auto& lj = tnodes[te.nids[j]];
            if (auto* v = std::get_if<HmLocOnV>(&lj)) {
                v_stop[v->id] = true;
                v_wall.insert(v->id);
                continue;
            }
            auto er = try_get_edge_ratio(hm, lj); if (!er) continue;
            auto [e, r] = *er;

            if (!lo.contains(e.id)) { lo[e.id] = 0.; hi[e.id] = 1.; }
            // Row3d, not auto: `auto` would keep an Eigen expression referencing the two temporaries
            Row3d dir = get_ptloc_pos(hm, tnodes[te.nids[k]]) - get_ptloc_pos(hm, lj);
            Row3d nrm = get_ptloc_nml(hm, lj);
            if (nrm.dot(dir.cross(e.vec())) > 0) lo[e.id] = std::max(lo[e.id], r); // inner side is [r, 1]
            else                                 hi[e.id] = std::min(hi[e.id], r); // inner side is [0, r]
        }
    }

    vec<std::tuple<int, double, double>> res;
    for (auto& [eid, l] : lo) if (l < hi[eid]) res.emplace_back(eid, l, hi[eid]);

    // flood the interior vertex by vertex: walls are the crossing edges and the OnV nodes;
    // every edge met on the way is fully admissible
    vec v_seen(hm.nV, false);
    vec e_done(hm.nE, false);
    for (const auto& eid: lo | std::views::keys) e_done[eid] = true;

    vec<int> stack;

    auto push_v = [&](int vid) {
        if (v_stop[vid] || v_seen[vid]) return;
        v_seen[vid] = true;
        stack.push_back(vid);
    };

    // 種：交差 edge の内側端点（lo==0 → tail 内側 / hi==1 → head 内側）
    for (auto& [eid, l] : lo) {
        if (!(l < hi[eid])) continue;
        Edge e = hm.edges[eid];
        if (l <= 0)       push_v(e.vert0().id);
        if (hi[eid] >= 1) push_v(e.vert1().id);
    }

    while (!stack.empty()) {
        int vid = stack.back(); stack.pop_back();
        for (Half h : hm.verts[vid].adjHalfs()) {
            int eid = h.edge().id;
            if (e_done[eid]) continue;
            e_done[eid] = true;
            res.emplace_back(eid, 0., 1.);
            push_v(h.head().id);
        }
    }

    for (int vid : v_wall)
    for (Half h : hm.verts[vid].adjHalfs()) {
        int eid = h.edge().id;
        int uid = h.head().id;
        if (!e_done[eid] && v_wall.contains(uid)) {
            e_done[eid] = true;
            res.emplace_back(eid, 0., 1.);
            push_v(h.head().id);
        }
    }

    return res;
}

vec<std::tuple<int, double, double>> Emesh::allowed_range_tquads(const vec<int>& tqids) const {
    vec<int> thids;
    for (int tqid: tqids)
    for (auto& [thid, _]: tquads[tqid].data)
        thids.push_back(thid);

    std::set ids(thids.begin(), thids.end());
    std::erase_if(thids, [&](int thid) { return ids.contains(thalfs[thid].twid); });
    return allowed_range_thalfs(thids);
}
}
