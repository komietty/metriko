#ifndef METRIKO_TMESH_MUT_COLLAPSE_RANGE_H
#define METRIKO_TMESH_MUT_COLLAPSE_RANGE_H
#include "tmesh_mut.h"
#include <queue>

namespace metriko {

vec<std::tuple<int, double, double>> TmeshMut::allowed_range(int tqid) const {

    auto lpos = [&](int nid) -> Row3d { return get_ptloc_pos(hm, tnodes[nid]); };
    auto edir = [&](int eid) -> Row3d { return hm.edges[eid].half().vec(); };

    const auto& data = tquads[tqid].data;
    // 各 edge 交点を左側判定し、edge ごとに内側区間 [lo,hi] を交差で絞る（2回横切れば [min,max]、
    // 1回なら [r,1] か [0,r]）。境界頂点（HmLocOnV）は flood の壁として v_stop に集める。
    umap<int, double> lo, hi;
    vec v_stop(hm.nV, false);
    for (auto& [thid, _] : data) {
        const auto& th = thalfs[thid];
        const auto& te = tedges[th.teid];
        int n = te.nids.size();
        for (int i = 0; i < n - 1; ++i) {
            int j   = th.cano ? i : n - 1 - i;
            int nid = te.nids[j];
            const HmLoc& loc = tnodes[nid];
            if (auto* v = std::get_if<HmLocOnV>(&loc)) { v_stop[v->id] = true; continue; }
            auto* e = std::get_if<HmLocOnE>(&loc);
            if (!e) continue;
            if (!lo.contains(e->id)) { lo[e->id] = 0.; hi[e->id] = 1.; }
            int jn = th.cano ? j + 1 : j - 1;
            Row3d dir = lpos(te.nids[jn]) - lpos(nid);
            Row3d nrm = get_ptloc_normal(hm, loc);
            if (nrm.dot(dir.cross(edir(e->id))) > 0) lo[e->id] = std::max(lo[e->id], e->r);  // 内側 ⊂ [r,1]
            else                                     hi[e->id] = std::min(hi[e->id], e->r);  // 内側 ⊂ [0,r]
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
}
#endif
