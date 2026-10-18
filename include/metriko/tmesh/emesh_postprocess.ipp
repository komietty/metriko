//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include "emesh.h"

namespace metriko {
inline void Emesh::collapse_tedge_snap(bool flag) {
    // the tedges with a segment in each face (a segment lies in the faces both its ends touch): an intersection test
    // in a face only needs these. a snap only adds segments next to the node it moves, so the faces of those are
    // added; a tedge that has left a face stays listed there, which costs a test but changes no answer
    vec<vec<int>> face_tes(hm.nF), node_tes(tnodes.size());
    for (auto& [id, ns]: live_tedges()) for (int n: ns) node_tes[n].push_back(id);
    auto add_segment = [&](int te, int k) {
        const auto& ns = tedges[te].nids;
        auto fb = get_ptloc_faces(hm, tnodes[ns[k + 1]]);
        for (int f: get_ptloc_faces(hm, tnodes[ns[k]]))
            if (rg::contains(fb, f) && !rg::contains(face_tes[f], te)) face_tes[f].push_back(te);
    };
    for (auto& [id, ns]: live_tedges()) for (int k = 0; k + 1 < std::ssize(ns); ++k) add_segment(id, k);
    auto reindex_node = [&](int nid) {  // the segments next to a moved node
        for (int te: node_tes[nid]) {
            const auto& ns = tedges[te].nids;
            for (int k = 0; k < std::ssize(ns); ++k) {
                if (ns[k] != nid) continue;
                if (k > 0) add_segment(te, k - 1);
                if (k + 1 < std::ssize(ns)) add_segment(te, k);
            }
        }
    };
    // the nodes on each vertex, and how often the live tedges pass each node
    vec<vec<int>> vert_nodes(hm.nV);
    vec<int> passes(tnodes.size(), 0);
    for (int n = 0; n < std::ssize(tnodes); ++n) if (auto* l = std::get_if<HmLocOnV>(&tnodes[n])) vert_nodes[l->id].push_back(n);
    for (auto& [id, ns]: live_tedges()) for (int n: ns) ++passes[n];
    auto erase_nodes = [&](vec<int>& ns, int b, int e) {
        for (int k = b; k < e; ++k) --passes[ns[k]];
        ns.erase(ns.begin() + b, ns.begin() + e);
    };
    auto move_node = [&](int nid, Vert v) {
        tnodes[nid] = HmLocOnV{.id = v.id};
        vert_nodes[v.id].push_back(nid);
        reindex_node(nid);
    };
    auto in_faces = [&](auto&& faces, auto&& test) {   // test(teid, k) on every segment of a tedge in one of the faces
        for (Face f: faces) for (int te: face_tes[f.id]) { const auto& ns = tedges[te].nids; for (int k = 0; k + 1 < ns.size(); ++k) if (test(f, te, k)) return true; }
        return false;
    };

    auto snap_valid = [&](Vert v) { return rg::none_of(vert_nodes[v.id], [&](int n) { return passes[n] > 0; }); };

    auto snap_joint = [&](int teid, Vert v) {
        auto  nid = tedges[teid].nids.back();
        auto* loc = std::get_if<HmLocOnF>(&tnodes[nid]);
        if (!loc || !snap_valid(v)) return;

        auto f = hm.faces[loc->id];
        if (in_faces(std::array{f}, [&](Face f, int te, int k) {
            const auto& nids = tedges[te].nids;
            if (nids[k] == nid || nids[k + 1] == nid) return false;
            return find_strict_intersection(f, *loc, HmLocOnV{v.id}, tnodes[nids[k]], tnodes[nids[k + 1]]);
        })) return;

        vec<std::pair<int, int>> te_tails; // teid, the biggest  nid
        vec<std::pair<int, int>> te_heads; // teid, the smallest nid
        for (auto& [id, nids]: live_tedges()) {
            if (nid == nids.front()) te_tails.emplace_back(id, 0);           // initial value
            if (nid == nids.back())  te_heads.emplace_back(id, nids.size()); // initial value
        }
        // only the consecutive run of in-star nodes at the joint end is trimmed: not leave and come back case...
        for (auto& [i, j]: te_tails) { auto& nids = tedges[i].nids; j = 0;               while (j + 1 < nids.size() && is_in_star(v, tnodes[nids[j + 1]])) ++j; }
        for (auto& [i, j]: te_heads) { auto& nids = tedges[i].nids; j = nids.size() - 1; while (j > 0               && is_in_star(v, tnodes[nids[j - 1]])) --j; }

        // todo: apply this
        // the legs the snap leaves: from v to the first node each joint tedge keeps. they must not cross or touch any
        // other segment (the trimmed part of their own tedge aside), or the contact is one the repair cannot undo
        // struct Leg { int teid; int keep; int k0; int k1; }; // the kept node, and the replaced segments [k0, k1)
        // vec<Leg> legs;
        // for (auto& [i, j]: te_tails) { int n = tedges[i].nids.size(); int p = std::min(std::max(j, 1), n - 1); legs.emplace_back(i, p, 0, p); }
        // for (auto& [i, j]: te_heads) { int n = tedges[i].nids.size(); int p = std::max(std::min(j, n - 2), 0); legs.emplace_back(i, p, p, n - 1); }
        // for (auto& [lt, keep, k0, k1]: legs) {
        //     const HmLoc& l = tnodes[tedges[lt].nids[keep]];
        //     for (Face g: v.adjHalfs() | vw::transform(&Half::face) | vw::filter([&](auto g) { return is_in_face(g, l); })) {
        //         for (auto& [id, ns]: live_tedges()) {
        //             for (int k = 0; k + 1 < ns.size(); ++k) {
        //                 if (rg::any_of(legs, [&](const Leg& o) { return o.teid == id && k >= o.k0 && k < o.k1; })) continue;
        //                 if (has_contact(g, l, HmLocOnV{v.id}, tnodes[ns[k]], tnodes[ns[k + 1]])) return;
        //             }}}
        // }

        for (auto& [i, j]: te_tails) { auto& nids = tedges[i].nids; if (j >= 2)              erase_nodes(nids, 1, j); }
        for (auto& [i, j]: te_heads) { auto& nids = tedges[i].nids; if (j + 2 < nids.size()) erase_nodes(nids, j + 1, nids.size() - 1); }
        move_node(nid, v);
    };

    auto snap_inter = [&](int teid, int nid, Vert v) -> bool {
        auto& nids = tedges[teid].nids;
        auto it = rg::find(nids, nid); if (it == nids.end()) return false;
        int i = it - nids.begin(), lo = i, hi = i;
        while (lo > 0               && is_in_star(v, tnodes[nids[lo - 1]])) --lo;
        while (hi < nids.size() - 1 && is_in_star(v, tnodes[nids[hi + 1]])) ++hi;
        int L = lo < i ? lo : i - 1;   // neighbours after the trim: the two new legs are L -> v -> R
        int R = hi > i ? hi : i + 1;

        for (auto& l: {tnodes[nids[L]], tnodes[nids[R]]}) {
            if (in_faces(v.adjHalfs() | vw::transform(&Half::face) | vw::filter([&](auto f) { return is_in_face(f, l); }), [&](Face f, int id, int k) {
                const auto& ns = tedges[id].nids;
                if (id == teid && k >= L && k < R) return false;
                return find_strict_intersection(f, l, HmLocOnV{v.id}, tnodes[ns[k]], tnodes[ns[k + 1]]);
            })) return false;
        }

        if (hi - i >= 2) erase_nodes(nids, i + 1, hi);
        if (i - lo >= 2) erase_nodes(nids, lo + 1, i);
        move_node(nid, v);
        return true;
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
            if (snap_valid(vrt) && snap_inter(eid, nid, vrt)) break;
        }}
    } else {
        for (auto& c: candidates) {
            if (c.size() == 0) continue;
            auto& [nid, eid, vrt, _] = c.front();
            if (snap_valid(vrt)) snap_inter(eid, nid, vrt);
        }
    }
}

inline void Emesh::collapse_tedge_snap_dedup(int teid) {
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

// retrace a tedge inside its two tquads together, to resolve a contact made by snapping
inline bool Emesh::reroute_tedge(int teid) {
    auto th = rg::find_if(thalfs, [&](const Ehalf& th) { return th.id != -1 && th.teid == teid && th.cano; });
    if (th == thalfs.end() || th->tqid < 0 || th->twin().tqid < 0) return false;
    auto& nids = tedges[teid].nids;
    if (nids.front() == nids.back()) return false;   // a closed tedge: no two ends to join

    // the boundary of the two tquads together, the tedge itself left out: from its head around to its tail
    auto rest = [&](const Equad& tq, int thid) {
        vec<int> res;
        for (auto it = circular_next(tq.data, tq.find_of(thid)); it->thid != thid; it = circular_next(tq.data, it)) res.push_back(it->thid);
        return res;
    };
    vec<int> loop = rest(th->tquad(), th->id);
    const int tail = loop.size();
    for (int t: rest(th->twin().tquad(), th->twid)) loop.push_back(t);
    const Unfolded chart = unfold(*this, loop);
    const auto  path  = retrace(*this, chart, chart.at(tail, nids.front()), chart.at(0, nids.back()));
    if (path.empty()) return false;
    nids = add_new_path(path, nids.front(), nids.back());

    auto r = path_length(nids);
    for (auto& t: thalfs) if (t.id != -1 && t.teid == teid) t.r = r;
    return true;
}
}
