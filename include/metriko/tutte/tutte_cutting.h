//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_TUTTE_CUTTING_H
#define METRIKO_TUTTE_CUTTING_H
#include <format>
#include <map>
#include <ranges>
#include <set>
#include <unordered_map>
#include "tutte.h"
#include "metriko/hmesh/hmesh.h"
#include "metriko/tmesh/emesh.h"

namespace metriko {
// split points on an original halfedge: (ratio, vid)
using EdgeSplits = std::set<std::pair<double, int>>;

struct SegRec {
    int i0;
    int i1;
    HalfData d0;
    HalfData d1;
};

inline bool is_on_side(const vec<vec<int>>& sides, int vid, int sid) {
    return vid < sides.size() && rg::contains(sides[vid], sid);
}

inline int common_side(const vec<vec<int>>& sides, int vid0, int vid1) {
    if (vid0 >= sides.size()) return -1;
    for (int l: sides[vid0])
        if (is_on_side(sides, vid1, l)) return l;
    return -1;
}

inline int common_side(const vec<vec<int>>& sides, int vid0, int vid1, int vid2) {
    if (vid0 >= sides.size()) return -1;
    for (int l: sides[vid0])
        if (is_on_side(sides, vid1, l) &&
            is_on_side(sides, vid2, l)) return l;
    return -1;
}

// vertex chain along an original halfedge, tail -> head, split points included
inline vec<int> chain_of(const vec<EdgeSplits>& splits, Half h) {
    vec chain = { h.tail().id };
    for (const auto& [r, vid]: vw::reverse(splits[h.id])) chain.push_back(vid);  // descending r = tail -> head
    chain.push_back(h.head().id);
    return chain;
}

// triangulate one face against the segments crossing it. the resulting pieces may
// be CONCAVE (traced tedges bend inside a face), so the trace allows reflex turns
// and the pieces are ear-clipped rather than fanned.
inline void face_cutting(
    const Face& f,                 // face to be cut
    const vec<Row2i>& sgs,         // segment endpoints (cut-mesh vertex indices)
    const vec<EdgeSplits>& splits, // split points per original halfedge
    const vec<vec<int>>& sides,    // vertex -> tquad * side
    const vec<Row3d>& vpos,        // vertex position
    vec<Row3i>& tris               // output triangles (appended, CCW wrt f)
) {
    /// 1: collect the directed edges bounding the pieces:
    ///    segments (both directions) + subdivided boundary halfedges (one direction)
    vec<Row2i> halfs;
    for (auto s: sgs) {
        halfs.emplace_back(s.x(), s.y());
        halfs.emplace_back(s.y(), s.x());
    }

    for (Half h: f.adjHalfs()) {
        auto ch = chain_of(splits, h);
        for (size_t k = 0; k + 1 < ch.size(); ++k) halfs.emplace_back(ch[k], ch[k + 1]);
    }

    auto cr   = [](complex u, complex v) { return (std::conj(u) * v).imag(); };

    constexpr double penalty = 100.; // dominates min_angle in (-pi, pi]

    /// 2: trace the pieces by the leftmost-turn rule (max signed CCW angle).
    ///    reflex turns are legal, and every chain must close.
    while (!halfs.empty()) {
        auto h = halfs.back();
        halfs.pop_back();

        vec poly = {h};
        int sta = h.x();
        int cur = h.y();

        while (cur != sta) {
            int prev_v = poly.back().x();
            int curr_v = poly.back().y();
            const auto p_curr = f.to_local(vpos[curr_v]);
            const auto d_prev = p_curr - f.to_local(vpos[prev_v]);

            auto   best_it = halfs.end();
            double best    = -10;   // signed turn angle in (-pi, pi], in the face plane like the ear clipping below
            for (auto it = halfs.begin(); it != halfs.end(); ++it) {
                if (it->x() != curr_v) continue;
                if (it->y() == prev_v) continue;   // no immediate u-turn
                double ang = std::arg((f.to_local(vpos[it->y()]) - p_curr) / d_prev);
                if (ang > best) { best = ang; best_it = it; }
            }

            METRIKO_CHECK(best_it != halfs.end(), "open chain in face {}", f.id);

            poly.push_back(*best_it);
            cur = best_it->y();
            halfs.erase(best_it);
        }

        /// 3: ear-clip the piece in the face plane (concave-capable, keeps CCW
        ///    orientation so no flipped faces can be produced)
        auto cyc = poly | vw::transform([](const Row2i& p) { return p.x(); }) | rg::to<vec<int>>();

        // interior angles of a CCW triangle in the face plane (all positive);
        // negative for a CW triangle, so it also ranks invalid choices last
        auto min_angle = [&](int va, int vb, int vc) {
            auto a = f.to_local(vpos[va]),
                 b = f.to_local(vpos[vb]),
                 c = f.to_local(vpos[vc]);
            return std::min({std::arg((c - a) / (b - a)),
                             std::arg((a - b) / (c - b)),
                             std::arg((b - c) / (a - c))});
        };

        METRIKO_CHECK(cyc.size() >= 3, "polygon with {} vertices", cyc.size());

        while (cyc.size() > 3) {
            size_t n = cyc.size();
            size_t best_k = n;
            double best_q = std::numeric_limits<double>::lowest();
            for (size_t k = 0; k < n; ++k) {
                int va = cyc[(k + n - 1) % n];
                int vb = cyc[k];
                int vc = cyc[(k + 1) % n];
                auto a = f.to_local(vpos[va]),
                     b = f.to_local(vpos[vb]),
                     c = f.to_local(vpos[vc]);
                if (cr(b - a, c - b) <= EPS * std::abs(b - a) * std::abs(c - b)) continue; // reflex or flat corner (sine below EPS): not an ear

                bool empty = true; // no other cycle vertex inside the ear
                for (size_t m = 0; m < n && empty; ++m) {
                    if (m == k || m == (k + 1) % n || m == (k + n - 1) % n) continue;
                    auto q = f.to_local(vpos[cyc[m]]);
                    empty = cr(b - a, q - a) < -EPS * std::abs(b - a) * std::abs(q - a)
                         || cr(c - b, q - b) < -EPS * std::abs(c - b) * std::abs(q - b)
                         || cr(a - c, q - c) < -EPS * std::abs(a - c) * std::abs(q - c);
                }
                if (!empty) continue;

                // quality of this clip: the worst interior angle it commits to.
                // a quad also fixes the remaining triangle, so score the pair —
                // this picks the better diagonal and avoids slivers whenever the
                // other choice is available
                double q = min_angle(va, vb, vc);
                if (n == 4) q = std::min(q, min_angle(vc, cyc[(k + 2) % n], va));
                if (common_side(sides, va, vb, vc) >= 0)                         q -= penalty;
                if (n == 4 && common_side(sides, vc, cyc[(k + 2) % n], va) >= 0) q -= penalty;
                if (q > best_q) { best_q = q; best_k = k; }
            }
            METRIKO_CHECK(best_k != n, "ear clipping failed in face {} (polygon of {})", f.id, cyc.size());
            tris.emplace_back(cyc[(best_k + n - 1) % n], cyc[best_k], cyc[(best_k + 1) % n]);
            cyc.erase(cyc.begin() + (long)best_k);
        }

        tris.emplace_back(cyc[0], cyc[1], cyc[2]);
    }
}

inline void split_degenerate_faces(
    const Hmesh& hm,
    const vec<bool>& seam0,
    const vec<SegRec>& sgs,
    vec<EdgeSplits>& h_aux,
    const vec<vec<int>>& sides,
    vec<Row3d>& vpos,
    vec<Row3i>& tris
) {
    std::set<std::pair<int, int>> sg_set;
    for (auto& s: sgs) sg_set.insert(std::minmax(s.i0, s.i1));

    // seam sub-edge -> (halfedge id, ratio of key.first, ratio of key.second):
    // lets a midpoint inserted on it be threaded back into h_aux so the
    // seam1 / matching1 propagation still finds every sub-halfedge
    std::map<std::pair<int, int>, std::tuple<int, double, double>> sm_info;

    for (Edge e: hm.edges) {
        if (!seam0[e.id]) continue;
        Half h  = e.half();
        auto ch = chain_of(h_aux, h);
        vec<double> rs = {1.}; // r = 1 at tail, 0 at head along h
        for (const auto& [r, _]: vw::reverse(h_aux[h.id])) rs.push_back(r);
        rs.push_back(0.);
        for (size_t k = 0; k + 1 < ch.size(); ++k) {
            auto key = std::minmax(ch[k], ch[k + 1]);
            sm_info[key] = key.first == ch[k]
                ? std::tuple{h.id, rs[k], rs[k + 1]}
                : std::tuple{h.id, rs[k + 1], rs[k]};
        }
    }

    std::map<std::pair<int, int>, int> tri_of; // directed edge -> tri on its left, kept up to date by split_at
    auto index   = [&](int i) { for (int j = 0; j < 3; ++j) tri_of[{tris[i][j], tris[i][(j + 1) % 3]}] = i; };
    auto unindex = [&](int i) { for (int j = 0; j < 3; ++j) tri_of.erase({tris[i][j], tris[i][(j + 1) % 3]}); };
    for (int i = 0; i < tris.size(); ++i) index(i);

    for (bool again = true; std::exchange(again, false);) {
        // split edge j of triangle i (and the neighbor across it) at the midpoint.
        // s: the side the degeneracy is on; rejected when the neighbor's opposite
        // vertex is also on it (the split would not resolve anything)
        auto split_at = [&](int i, int j, int s) -> bool {
            int a = tris[i][j];
            int b = tris[i][(j + 1) % 3];
            int c = tris[i][(j + 2) % 3];
            if (sg_set.contains(std::minmax(a, b))) return false;
            int k = tri_of.at({b, a});
            int d = tris[k][0] + tris[k][1] + tris[k][2] - a - b;
            if (is_on_side(sides, d, s)) return false;
            int m = vpos.size();

            if (auto it = sm_info.find(std::minmax(a, b)); it != sm_info.end()) {
                // the edge lies on a seam: thread the midpoint into h_aux
                auto [hid, r_lo, r_hi] = it->second;
                double rm = (r_lo + r_hi) / 2;
                h_aux[hid].emplace(rm, m);
                h_aux[hm.halfs[hid].twin().id].emplace(1 - rm, m);
                auto [lo, hi] = std::minmax(a, b);
                sm_info.erase(it);
                sm_info[{lo, m}] = {hid, r_lo, rm}; // m is always the largest id
                sm_info[{hi, m}] = {hid, r_hi, rm};
            }

            Row3d mid = (vpos[a] + vpos[b]) / 2;
            unindex(i);
            unindex(k);
            tris[i] = {a, m, c};
            tris[k] = {b, m, d};
            tris.emplace_back(m, b, c);
            tris.emplace_back(m, a, d);
            for (int t: { i, k, (int)tris.size() - 2, (int)tris.size() - 1 }) index(t);
            vpos.emplace_back(mid);
            return true;
        };

        // fully collinear triangles
        for (int i = 0; i < tris.size() && !again; ++i) {
            int l = common_side(sides, tris[i].x(), tris[i].y(), tris[i].z());
            if (l < 0) continue;
            for (int j = 0; j < 3 && !again; ++j) again = split_at(i, j, l);
        }

        // interior chords on a side
        for (int i = 0; i < tris.size() && !again; ++i)
        for (int j = 0; j < 3 && !again; ++j) {
            int a = tris[i][j];
            int b = tris[i][(j + 1) % 3];
            int c = tris[i][(j + 2) % 3];
            if (int s = common_side(sides, a, b); s >= 0 && !is_on_side(sides, c, s)) again = split_at(i, j, s);
        }
    }
}

inline std::unique_ptr<Hmesh> compute_embedding_cut_hmesh(
    const Hmesh& hm,
    const Emesh& tm,
    const vec<bool>& seam0,
    const VecXi& matching0,
    const VecXi& singular0,
          vec<bool>& seam1,
          VecXi& matching1,
          VecXi& singular1,
    vec<HalfData>& data
) {
    std::map<int, vec<Row2i>> cuts; // face id -> segment endpoints

    vec<SegRec> sgms; // flat annotation list; matched to cut halfedges at the end
    vec auxs(hm.nH, EdgeSplits{});

    vec<Row3d> vpos;
    for (auto p: hm.pos.rowwise()) { vpos.emplace_back(p); }

    // nid -> cut-mesh vertex index: tnode identity replaces epsilon-based position
    // matching. creates the vertex and registers the edge split exactly once.
    umap<int, int> vid_of_nid;
    vec<vec<int>> sides; // cut vertex -> tquad * side

    auto vid_of = [&](int nid) -> int {
        if (!vid_of_nid.contains(nid)) {
            vid_of_nid[nid] = std::visit(overloaded{
                [&](const HmLocOnV& v) { return v.id; },
                [&](const auto& l) {
                    vpos.emplace_back(get_ptloc_pos(hm, l));
                    int i = vpos.size() - 1;
                    if (auto hr = try_get_half_ratio(hm, l)) {
                        auto [h0, r] = hr.value();
                        auxs[h0.id].emplace(r, i);
                        auxs[h0.twin().id].emplace(1 - r, i);
                    }
                    return i;
                },
            }, tm.tnodes[nid]);
        }
        return vid_of_nid[nid];
    };

    for (const auto& th: tm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids = tm.tedges[th.teid].nids;
        const int   nsgs = nids.size() - 1;

        const auto& tw = tm.thalfs[th.twid];
        const int   sa = th.tqid * 4 + tm.tquads[th.tqid].side_of(th);
        const int   sb = tw.tqid * 4 + tm.tquads[tw.tqid].side_of(tw);

        double total = 0;
        vec<double> len(nsgs);
        for (int k = 0; k < nsgs; ++k) {
            len[k] = (get_ptloc_pos(hm, tm.tnodes[nids[k + 1]]) - get_ptloc_pos(hm, tm.tnodes[nids[k]])).norm();
            total += len[k];
        }
        METRIKO_CHECK(total >= EPS, "zero-length tedge");

        double sum = 0;
        for (int k = 0; k < nsgs; ++k) {
            auto& fr = tm.tnodes[nids[k]];
            auto& to = tm.tnodes[nids[k + 1]];
            auto  v0 = sum / total; sum += len[k];
            auto  v1 = sum / total;
            auto  i0 = vid_of(nids[k]);
            auto  i1 = vid_of(nids[k + 1]);
            for (int i: {i0, i1}) {
                if (i >= sides.size()) sides.resize(vpos.size());
                auto& ss = sides[i];
                if (!rg::contains(ss, sa)) ss.push_back(sa);
                if (!rg::contains(ss, sb)) ss.push_back(sb);
            }
            auto eo = try_get_edge(hm, fr, to);
            auto fo = try_get_face(hm, fr, to);
            auto t  = nsgs - k - 1;
            if (!eo.has_value()) cuts[fo.value().id].emplace_back(i0, i1);
            sgms.push_back({
                .i0 = i0, .i1 = i1,
                .d0 = HalfData{.v0 = v0,     .v1 = v1,     .thid = th.id,   .tqid = th.tqid,                 .order = k},
                .d1 = HalfData{.v0 = 1 - v1, .v1 = 1 - v0, .thid = th.twid, .tqid = tm.thalfs[th.twid].tqid, .order = t}
            });
        }
    }

    vec<Row3i> tris;
    for (Face f: hm.faces) {
        const auto it    = cuts.find(f.id);
        const bool cut   = it != cuts.end();
        const bool split = rg::any_of(f.halfs(), [&](Half h) { return !auxs[h.id].empty(); });
        if (!cut && !split) { auto [a, b, c] = f.verts(); tris.emplace_back(a.id, b.id, c.id); }
        else face_cutting(f, cut ? it->second : vec<Row2i>{}, auxs, sides, vpos, tris);
    }

    split_degenerate_faces(hm, seam0, sgms, auxs, sides, vpos, tris);

    MatXd vert_info = Eigen::Map<MatX3d>(vpos[0].data(), vpos.size(), 3);
    MatXi face_info = Eigen::Map<MatX3i>(tris[0].data(), tris.size(), 3);

    auto hm_cut = std::make_unique<Hmesh>(vert_info, face_info);
    seam1     = vec(hm_cut->nE, false);
    matching1 = VecXi::Zero(hm_cut->nE);
    singular1 = VecXi::Zero(hm_cut->nV);
    singular1.head(hm.nV) = singular0;

    std::map<std::pair<int, int>, Half> half_by_verts;
    for (Half h: hm_cut->halfs)
        half_by_verts.insert({{h.tail().id, h.head().id}, h});

    for (Edge e: hm.edges) {
        if (!seam0[e.id]) continue;
        auto ch = chain_of(auxs, e.half());
        for (int k = 0; k + 1 < ch.size(); ++k) {
            Half h = half_by_verts.at({ch[k], ch[k + 1]}); // same direction as h
            seam1[h.edge().id] = true;
            matching1(h.edge().id) = h.isCanonical() ? matching0(e.id) : -matching0(e.id);
        }
    }

    // annotate the cut halfedges: both directions pushed adjacently, twin = index ^ 1
    data.clear();
    for (auto& [i0, i1, d0, d1]: sgms) {
        auto it = half_by_verts.find({i0, i1});
        METRIKO_CHECK(it != half_by_verts.end(), "no half between vertices");
        d0.half = it->second;
        d1.half = it->second.twin();
        data.push_back(d0);
        data.push_back(d1);
    }

    std::sort(data.begin(), data.end());
    return hm_cut;
}
}
#endif
