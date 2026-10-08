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
inline vec<Erng> Emesh::allowed_range_thalfs(const vec<int>& thids) const {
    vec<Erng> res;
    vec v_stop(hm.nV, false);
    umap<int, Erng> rngs;
    std::set<int> v_wall = {};

    for (int thid : thids) {
        const auto& th = thalfs[thid];
        const auto& te = th.tedge();
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

            auto& rng = rngs.try_emplace(e.id, Erng{.id = e.id, .fr = 0., .to = 1.}).first->second;
            Row3d dif = get_ptloc_dif(hm, lj, tnodes[te.nids[k]]);
            Row3d nrm = get_ptloc_nml(hm, lj);
            if (nrm.dot(dif.cross(e.vec())) > 0) rng.clip(r, 1.);
            else                                 rng.clip(0., r);
        }
    }

    vec v_seen(hm.nV, false);
    vec e_done(hm.nE, false);
    for (const auto& [eid, r] : rngs) if (!r.empty()) res.push_back(r);
    for (const auto& [eid, r] : rngs) e_done[eid] = true;

    vec<int> stack;

    auto push_v = [&](int vid) {
        if (v_stop[vid] || v_seen[vid]) return;
        v_seen[vid] = true;
        stack.push_back(vid);
    };

    for (const auto& [eid, rng] : rngs) {
        if (rng.empty()) continue;
        if (rng.fr <= 0.) push_v(hm.edges[eid].vert0().id);
        if (rng.to >= 1.) push_v(hm.edges[eid].vert1().id);
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

inline vec<Erng> Emesh::allowed_range_tquads(const vec<int>& tqids) const {
    vec<int> thids;
    for (int tqid: tqids)
    for (auto& [thid, _]: tquads[tqid].data)
        thids.push_back(thid);

    std::set ids(thids.begin(), thids.end());
    std::erase_if(thids, [&](int thid) { return ids.contains(thalfs[thid].twid); });
    return allowed_range_thalfs(thids);
}
}
