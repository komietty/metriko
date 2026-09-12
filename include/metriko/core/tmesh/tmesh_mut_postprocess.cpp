#include <set>
#include "./tmesh_mut.h"

namespace metriko {
bool TmeshMut::collapse_valid_snap_0(Vert v) {
    for (auto& [id, nids]: live_tedges()) {
        if (rg::any_of(nids, [&](int nid) {
            auto* l = std::get_if<HmLocOnV>(&tnodes[nid]);
            return l && l->id == v.id;
        })) return false;
    }
    return true;
}

// validate and apply one joint candidate: move the OnF joint ending tedge
// `teid` onto vertex v, trimming the incident chains inside v's one-ring.
// returns false when the joint is already snapped, v is occupied, or the move
// would cross another tedge inside the joint's face
bool TmeshMut::collapse_tedge_snap_joint(int teid, Vert v) {
    auto  nid = tedges[teid].nids.back();
    auto* loc = std::get_if<HmLocOnF>(&tnodes[nid]);
    if (!loc || !collapse_valid_snap_0(v)) return false;

    // moving the joint from p to v must not cross another tedge inside its face
    Face f = hm.faces[loc->id];
    auto a = f.to_local(get_ptloc_pos(hm, *loc));
    auto b = f.to_local(v.pos());
    for (auto& [id, nids]: live_tedges()) {
    for (int k = 0; k + 1 < nids.size(); ++k) {
        if (nids[k] == nid || nids[k + 1] == nid) continue;
        if (!is_in_face(f, tnodes[nids[k]]) || !is_in_face(f, tnodes[nids[k + 1]])) continue;
        double rab;
        double rcd;
        auto c = f.to_local(get_ptloc_pos(hm, tnodes[nids[k]]));
        auto d = f.to_local(get_ptloc_pos(hm, tnodes[nids[k + 1]]));
        if (find_strict_intersection(a, b, c, d, rab, rcd)) return false;
    }}

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
    return true;
}

void TmeshMut::collapse_tedge_snap_dedup(int teid) {
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
        vec<int> fids_curr = get_ptloc_faces(hm, tnodes[nid_curr]);
        vec<int> fids_prev = get_ptloc_faces(hm, tnodes[nid_prev]);
        vec<int> fids_next = get_ptloc_faces(hm, tnodes[nid_next]);

        bool same_face = rg::any_of(fids_curr, [&](int fid) { return rg::contains(fids_prev, fid) && rg::contains(fids_next, fid); });
        if (!same_face) { ++i; continue; }

        // drop the middle node; keep the removal only if the validation holds
        nids.erase(nids.begin() + i);
    }
}

void TmeshMut::collapse_tedge_snap_inter(int teid, int nid, Vert v) {
    auto in_ring = [&](const HmLoc& l) { return rg::any_of(v.adjHalfs(), [&](Half h) { return is_in_face(h.face(), l); });};

    auto& nids = tedges[teid].nids;
    auto it = rg::find(nids, nid);
    if (it == nids.end()) return;   // trimmed away by an earlier snap
    int i_cur = it - nids.begin();
    int i_min = i_cur;
    int i_max = i_cur;
    while (i_min > 0               && in_ring(tnodes[nids[i_min - 1]])) --i_min;
    while (i_max < nids.size() - 1 && in_ring(tnodes[nids[i_max + 1]])) ++i_max;

    if (i_max - i_cur >= 2) nids.erase(nids.begin() + i_cur + 1, nids.begin() + i_max);
    if (i_cur - i_min >= 2) nids.erase(nids.begin() + i_min + 1, nids.begin() + i_cur);
    tnodes[nid] = HmLocOnV{.id = v.id};
}

// re-trace a tedge inside the union corridor of its two tquads. used to resolve
// tedge-tedge contacts created by snapping: the corridor walls are the other
// boundary tedges, so the new path cannot touch them by construction
bool TmeshMut::reroute_tedge(int teid) {
    auto ths = thalfs | vw::filter([&](const ThalfMut& th) { return th.id != -1 && th.teid == teid; });
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

void TmeshMut::collapse_tedge_snap(bool flag) {
    struct Cand { int nid; int eid; Vert v; double d; };

    // 1: snap joint tnodes. every (joint, face vertex) pair is a candidate,
    //    processed nearest first; the validation happens at apply time because
    //    earlier snaps change the occupancy and the geometry
    {
        vec<Cand> cands;
        for (auto& [teid, nids]: live_tedges()) {
            auto* loc = std::get_if<HmLocOnF>(&tnodes[nids.back()]);
            if (!loc) continue;
            auto p = get_ptloc_pos(hm, *loc);
            for (Vert v: hm.faces[loc->id].adjHalfs() | vw::transform(&Half::tail))
                cands.push_back({.eid=teid, .v=v, .d=(v.pos() - p).squaredNorm()});
        }
        rg::sort(cands, {}, &Cand::d);
        for (auto& [nid, teid, v, d]: cands) collapse_tedge_snap_joint(teid, v);
    }

    vec candidates(tnodes.size(), vec<Cand>{});

    // 2: snap inter tnodes
    for (auto& [teid, nids]: live_tedges()) {
    for (int i = 1; i < nids.size() - 1; i++) {
        auto nid = nids[i];
        auto pos = get_ptloc_pos(hm, tnodes[nid]);

        if (std::holds_alternative<HmLocOnV>(tnodes[nid])) continue;
        for (int vid: get_ptloc_verts(hm, tnodes[nid])) {
            auto v = hm.verts[vid];
            auto d = (v.pos() - pos).squaredNorm();
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
            if (collapse_valid_snap_0(vrt)) { collapse_tedge_snap_inter(eid, nid, vrt); break; }
        }}
    } else {
        for (auto& c: candidates) {
            if (c.size() == 0) continue;
            auto& [nid, eid, vrt, _] = c.front();
            if (collapse_valid_snap_0(vrt)) { collapse_tedge_snap_inter(eid, nid, vrt); }
        }
    }
}
}
