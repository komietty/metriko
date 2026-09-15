#include "./emesh.h"

namespace metriko {
void Emesh::collapse_tedge_snap(bool flag) {

    auto snap_valid = [&](Vert v) {
        for (auto& [id, nids]: live_tedges()) {
            if (rg::any_of(nids, [&](int nid) {
                auto* l = std::get_if<HmLocOnV>(&tnodes[nid]);
                return l && l->id == v.id;
            })) return false;
        }
        return true;
    };

    auto snap_joint = [&](int teid, Vert v) {
        auto  nid = tedges[teid].nids.back();
        auto* loc = std::get_if<HmLocOnF>(&tnodes[nid]);
        if (!loc || !snap_valid(v)) return;
        // moving the joint from p to v must not cross another tedge inside its face
        auto f = hm.faces[loc->id];
        auto a = f.to_local(get_ptloc_pos(hm, *loc));
        auto b = f.to_local(v.pos());
        for (auto& [id, nids]: live_tedges()) {
        for (int k = 0; k + 1 < nids.size(); ++k) {
            if (nids[k] == nid || nids[k + 1] == nid) continue;
            if (!is_in_face(f, tnodes[nids[k]]) || !is_in_face(f, tnodes[nids[k + 1]])) continue;
            auto c = f.to_local(get_ptloc_pos(hm, tnodes[nids[k]]));
            auto d = f.to_local(get_ptloc_pos(hm, tnodes[nids[k + 1]]));
            if (find_strict_intersection(a, b, c, d)) return;
        }}

        //for (auto& [id, nids]: live_tedges()) {
        //    int n = nids.size();
        //    if (nid == nids.front()) { int j = 0;  for (int k = 0; k < n; k++) if (in_ring(v, tnodes[nids[k]])) j = std::max(j, k); if (j >= 2)    nids.erase(nids.begin() + 1,     nids.begin() + j); }
        //    if (nid == nids.back())  { int j = n;  for (int k = 0; k < n; k++) if (in_ring(v, tnodes[nids[k]])) j = std::min(j, k); if (j + 2 < n) nids.erase(nids.begin() + j + 1, nids.end() - 1);   }
        //}

        // trim the incident chains inside v's one-ring and move the joint
        vec<std::pair<int, int>> te_tails; // teid, the biggest  nid
        vec<std::pair<int, int>> te_heads; // teid, the smallest nid
        for (auto& [id, nids]: live_tedges()) {
            if (nid == nids.front()) te_tails.emplace_back(id, 0);
            if (nid == nids.back())  te_heads.emplace_back(id, nids.size());
        }
        for (Face g: v.adjHalfs() | vw::transform(&Half::face)) {
            for (auto& [i, j]: te_tails) { auto& nids = tedges[i].nids; for (int k = 0; k < nids.size(); k++) { if (is_in_face(g, tnodes[nids[k]])) j = std::max(j, k); }}
            for (auto& [i, j]: te_heads) { auto& nids = tedges[i].nids; for (int k = 0; k < nids.size(); k++) { if (is_in_face(g, tnodes[nids[k]])) j = std::min(j, k); }}
        }
        for (auto& [i, j]: te_tails) { auto& nids = tedges[i].nids; if (j >= 2)              nids.erase(nids.begin() + 1, nids.begin() + j);   }
        for (auto& [i, j]: te_heads) { auto& nids = tedges[i].nids; if (j + 2 < nids.size()) nids.erase(nids.begin() + j + 1, nids.end() - 1); }
        tnodes[nid] = HmLocOnV{.id = v.id};
    };

    auto snap_inter = [&](int teid, int nid, Vert v) {
        auto& nids = tedges[teid].nids;
        auto it = rg::find(nids, nid); if (it == nids.end()) return;
        int i = it - nids.begin(), lo = i, hi = i;
        while (lo > 0               && is_in_ring(v, tnodes[nids[lo - 1]])) --lo;
        while (hi < nids.size() - 1 && is_in_ring(v, tnodes[nids[hi + 1]])) ++hi;
        if (hi - i >= 2) nids.erase(nids.begin() + i + 1, nids.begin() + hi);
        if (i - lo >= 2) nids.erase(nids.begin() + lo + 1, nids.begin() + i);
        tnodes[nid] = HmLocOnV{.id = v.id};
    };


    struct Cand { int nid; int eid; Vert v; double d; };

    // 1: snap joint tnodes. every (joint, face vertex) pair is a candidate,
    //    processed nearest first; the validation happens at apply time because
    //    earlier snaps change the occupancy and the geometry
    vec<Cand> cands;
    for (auto& [teid, nids]: live_tedges()) {
        auto* l = std::get_if<HmLocOnF>(&tnodes[nids.back()]); if (!l) continue;
        auto  p = get_ptloc_pos(hm, *l);
        for (Vert v: hm.faces[l->id].verts())
            cands.push_back({.eid=teid, .v=v, .d=(v.pos() - p).squaredNorm()});
    }
    rg::sort(cands, {}, &Cand::d);
    for (auto& [nid, teid, v, d]: cands) { snap_joint(teid, v); }
    vec candidates(tnodes.size(), vec<Cand>{});

    // 2: snap inter tnodes
    for (auto& [teid, nids]: live_tedges()) {
    for (int i = 1; i < nids.size() - 1; i++) {
        auto nid = nids[i];
        auto n = tnodes[nid];
        auto p = get_ptloc_pos(hm, n);

        if (std::holds_alternative<HmLocOnV>(n)) continue;
        for (int vid: get_ptloc_verts(hm, n)) {
            auto v = hm.verts[vid];
            auto d = (v.pos() - p).squaredNorm();
            candidates[nid].push_back({.nid = nid, .eid = teid, .v = v, .d = d});
        }
    }}

    // 3: sort candidates in inner/outer order
    for (auto& c: candidates) rg::sort(c, {}, &Cand::d);
    rg::sort(candidates, {}, [](const vec<Cand>& c) { return c.empty() ? 1e9 : c.front().d; });

    // 4: snap if it's valid
    //for (auto& c: candidates)
    //for (auto& [nid, eid, vrt, _]: c | vw::take(flag ? c.size() : 1))
    //    if (collapse_valid_snap_0(vrt)) { collapse_tedge_snap_inter(eid, nid, vrt); break; }

    if (flag) {
        for (auto& c: candidates) {
        for (auto& [nid, eid, vrt, _]: c) {
            if (snap_valid(vrt)) { snap_inter(eid, nid, vrt); break; }
        }}
    } else {
        for (auto& c: candidates) {
            if (c.size() == 0) continue;
            auto& [nid, eid, vrt, _] = c.front();
            if (snap_valid(vrt)) { snap_inter(eid, nid, vrt); }
        }
    }
}

void Emesh::collapse_tedge_snap_dedup(int teid) {
    auto& nids = tedges[teid].nids;

    // nodes shared with other chains (junctions / crossings) must stay: dropping
    // one from this chain only would let the straightened chord cross the other
    // chain inside a face without a shared vertex
    vec shared(tnodes.size(), false);
    for (auto& [id, nids_]: live_tedges()) {
        if (id == teid) continue;
        for (int nid: nids_) shared[nid] = true;
    }

    for (int i = 1; i < nids.size() - 1;) {
        int nid_prev = nids[i - 1];
        int nid_curr = nids[i];
        int nid_next = nids[i + 1];

        if (shared[nid_curr] || nid_prev == nid_next) { ++i; continue; }

        // a face carrying all three nodes: the chord prev-next stays inside it
        auto fids_curr = get_ptloc_faces(hm, tnodes[nid_curr]);
        auto fids_prev = get_ptloc_faces(hm, tnodes[nid_prev]);
        auto fids_next = get_ptloc_faces(hm, tnodes[nid_next]);

        bool same_face = rg::any_of(fids_curr, [&](int fid) { return rg::contains(fids_prev, fid) && rg::contains(fids_next, fid); });
        if (!same_face) { ++i; continue; }

        // drop the middle node; keep the removal only if the validation holds
        nids.erase(nids.begin() + i);
    }
}

// re-trace a tedge inside the union corridor of its two tquads. used to resolve
// tedge-tedge contacts created by snapping: the corridor walls are the other
// boundary tedges, so the new path cannot touch them by construction
bool Emesh::reroute_tedge(int teid) {
    auto ths = thalfs | vw::filter([&](const Ehalf& th) { return th.id != -1 && th.teid == teid; });
    auto tq0 = -1;
    auto tq1 = -1;
    for (auto& th: ths) (th.cano ? tq0 : tq1) = th.tqid;
    if (tq0 < 0 || tq1 < 0) return false;

    auto  allw = allowed_range_tquads({tq0, tq1});
    auto& nids = tedges[teid].nids;
    auto  path = approx_shortest_path(30, hm, tnodes[nids.front()], tnodes[nids.back()], allw);
    if (path.size() < 2) return false;
    nids = add_new_path(path, nids.front(), nids.back());

    auto r = path_length(nids);
    for (auto& th: ths) th.r = r;
    return true;
}
}
