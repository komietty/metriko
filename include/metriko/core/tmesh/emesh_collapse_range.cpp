#include "emesh.h"
#include <set>

namespace metriko {

vec<std::tuple<int, double, double>> Emesh::allowed_range_thalfs(const vec<int>& thids) const {
    umap<int, double> lo, hi;
    vec v_stop(hm.nV, false);

    for (int thid : thids) {
        const auto& th = thalfs[thid];
        const auto& te = tedges[th.teid];
        int n = te.nids.size();
        for (int i = 0; i < n - 1; ++i) {
            auto  j  = th.cano ? i : n - 1 - i;
            auto  k  = th.cano ? j + 1 : j - 1;
            auto& lj = tnodes[te.nids[j]];
            if (auto* v = std::get_if<HmLocOnV>(&lj)) { v_stop[v->id] = true; continue; }
            auto er = try_get_edge_ratio(hm, lj); if (!er) continue;
            auto [e, r] = *er;

            if (!lo.contains(e.id)) { lo[e.id] = 0.; hi[e.id] = 1.; }
            auto dir = get_ptloc_pos(hm, tnodes[te.nids[k]]) - get_ptloc_pos(hm, lj);
            auto nrm = get_ptloc_nml(hm, lj);
            if (nrm.dot(dir.cross(e.vec())) > 0) lo[e.id] = std::max(lo[e.id], r);  // 内側 ⊂ [r,1]
            else                                 hi[e.id] = std::min(hi[e.id], r);  // 内側 ⊂ [0,r]
        }
    }

    vec<std::tuple<int, double, double>> res;
    for (auto& [eid, l] : lo) if (l < hi[eid]) res.emplace_back(eid, l, hi[eid]);

    // --- 隣接 edge を stack で flood し、内部 edge を range [0,1] で追加 ---
    // 壁 = 交差 edge（候補）と境界頂点（HmLocOnV ノード）。交差 edge の「内側端点」を種に、
    // 頂点づたいに edge を辿り、他の候補 edge か境界頂点に当たるまで内部 edge を集める。
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
