#ifndef METRIKO_TUTTE_CUTTING_H
#define METRIKO_TUTTE_CUTTING_H
#include <set>
#include <unordered_map>
#include "./tutte.h"
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/tmesh/tmesh_mut.h"

namespace metriko {

// split points on an original halfedge: (ratio in the tail-weighted convention,
// cut-mesh vertex index), ordered along the halfedge
using EdgeSplits = std::set<std::pair<double, int>>;

// a directed edge of the subdivided face; exists only while cutting one face
struct DirEdge { int i0; int i1; };

inline std::optional<std::pair<Half, double>> half_ratio_of(const Hmesh& hm, const HmLoc& l) {
    return std::visit(overloaded{
        [&](const HmLocOnE& e) -> std::optional<std::pair<Half, double>> {
            Edge ed = hm.edges[e.id];
            Half h  = ed.half();
            double t = h.tail() == ed.vert0() ? e.r : 1. - e.r;  // t: 0 at h.tail
            return std::pair(h, 1. - t);
        },
        [&](const HmLocOnH& hl) -> std::optional<std::pair<Half, double>> {
            return std::pair(hm.halfs[hl.id], 1. - hl.r);
        },
        [&](const auto&) -> std::optional<std::pair<Half, double>> { return std::nullopt; },
    }, l);
}

// triangulate one face against the segments crossing it. the resulting pieces may
// be CONCAVE (traced tedges bend inside a face), so the trace allows reflex turns
// and the pieces are ear-clipped rather than fanned.
inline void face_cutting(
    const Face& f,                       // face to be cut
    const vec<std::pair<int, int>>& sgs, // segment endpoints (cut-mesh vertex indices)
    const vec<Row3d>& vpos,              // vertex position
    const vec<EdgeSplits>& splits,       // split points per original halfedge
    const std::unordered_map<int, vec<int>>& lines, // vertex -> (tquad,side) lines through it
    vec<std::array<int, 3>>& tris        // output triangles (appended, CCW wrt f)
) {
    if (sgs.empty() && splits[f.half().id].empty() && splits[f.half().next().id].empty() && splits[f.half().prev().id].empty()) {
        tris.push_back({f.half().tail().id, f.half().head().id, f.half().prev().tail().id});
        return;
    }

    /// 1: collect the directed edges bounding the pieces:
    ///    segments (both directions) + subdivided boundary halfedges (one direction)
    vec<DirEdge> halfs;
    for (auto [i0, i1]: sgs) {
        halfs.push_back({i0, i1});
        halfs.push_back({i1, i0});
    }
    for (Half h: f.adjHalfs()) {
        const auto& sp = splits[h.id];
        if (sp.empty()) { halfs.push_back({h.tail().id, h.head().id}); continue; }
        halfs.push_back({sp.begin()->second, h.head().id});
        auto prev = sp.begin();
        auto curr = std::next(sp.begin());
        for (; curr != sp.end(); ++curr, ++prev) { halfs.push_back({curr->second, prev->second}); }
        halfs.push_back({h.tail().id, sp.rbegin()->second});
    }

    const Row3d org = f.half().tail().pos();
    const Row3d bx  = f.basisX();
    const Row3d by  = f.basisY();
    auto pos2 = [&](int vid) { Row3d d = vpos[vid] - org; return complex(d.dot(bx), d.dot(by)); };
    auto cr   = [](complex u, complex v) { return (std::conj(u) * v).imag(); };

    // three vertices on one patch side are collinear in the tutte uv (a side maps
    // to a straight rectangle side, possibly spanning several tedges joined at
    // junctions), so such a triangle collapses in uv even when well-shaped in 3d
    auto on_one_line = [&](int va, int vb, int vc) {
        auto ia = lines.find(va); if (ia == lines.end()) return false;
        auto ib = lines.find(vb); if (ib == lines.end()) return false;
        auto ic = lines.find(vc); if (ic == lines.end()) return false;
        for (int x: ia->second) for (int y: ib->second) for (int z: ic->second)
            if (x == y && y == z) return true;
        return false;
    };
    constexpr double penalty = 100.; // dominates min_angle in (-pi, pi]

    /// 2: trace the pieces by the leftmost-turn rule (max signed CCW angle).
    ///    reflex turns are legal, and every chain must close.
    while (!halfs.empty()) {
        DirEdge h = halfs.back();
        halfs.pop_back();

        vec poly = {h};
        int sta = h.i0;
        int cur = h.i1;

        while (cur != sta) {
            int prev_v = poly.back().i0;
            int curr_v = poly.back().i1;
            const Row3d& p_prev = vpos[prev_v];
            const Row3d& p_curr = vpos[curr_v];
            const Row3d  d_prev = (p_curr - p_prev).normalized();

            auto   best_it = halfs.end();
            double best    = -10;   // signed turn angle in (-pi, pi]
            for (auto it = halfs.begin(); it != halfs.end(); ++it) {
                if (it->i0 != curr_v) continue;
                if (it->i1 == prev_v) continue;   // no immediate u-turn
                Row3d d_cand = (vpos[it->i1] - p_curr).normalized();
                double ang = std::atan2(f.normal().dot(d_prev.cross(d_cand)), d_prev.dot(d_cand));
                if (ang > best) { best = ang; best_it = it; }
            }

            if (best_it == halfs.end()) throw std::runtime_error(std::format("face_cutting: open chain (face {})", f.id));

            poly.push_back(*best_it);
            cur = best_it->i1;
            halfs.erase(best_it);
        }

        /// 3: ear-clip the piece in the face plane (concave-capable, keeps CCW
        ///    orientation so no flipped faces can be produced)
        vec<int> cyc;
        for (auto& [i0, i1]: poly) cyc.push_back(i0);
        if (cyc.size() < 3) continue;

        // interior angles of a CCW triangle in the face plane (all positive);
        // negative for a CW triangle, so it also ranks invalid choices last
        auto min_angle = [&](int va, int vb, int vc) {
            complex a = pos2(va), b = pos2(vb), c = pos2(vc);
            return std::min({std::arg((c - a) / (b - a)),
                             std::arg((a - b) / (c - b)),
                             std::arg((b - c) / (a - c))});
        };

        while (cyc.size() > 3) {
            const size_t n = cyc.size();
            size_t best_k = n;
            double best_q = std::numeric_limits<double>::lowest();
            for (size_t k = 0; k < n; ++k) {
                int va = cyc[(k + n - 1) % n];
                int vb = cyc[k];
                int vc = cyc[(k + 1) % n];
                complex a = pos2(va), b = pos2(vb), c = pos2(vc);
                if (cr(b - a, c - b) <= EPS) continue;   // reflex or flat corner: not an ear

                bool empty = true;                        // no other cycle vertex inside the ear
                for (size_t m = 0; m < n && empty; ++m) {
                    if (m == k || m == (k + 1) % n || m == (k + n - 1) % n) continue;
                    complex q = pos2(cyc[m]);
                    empty = cr(b - a, q - a) < -EPS || cr(c - b, q - b) < -EPS || cr(a - c, q - c) < -EPS;
                }
                if (!empty) continue;

                // quality of this clip: the worst interior angle it commits to.
                // a quad also fixes the remaining triangle, so score the pair —
                // this picks the better diagonal and avoids slivers whenever the
                // other choice is available
                double q = min_angle(va, vb, vc);
                if (n == 4) q = std::min(q, min_angle(vc, cyc[(k + 2) % n], va));
                if (on_one_line(va, vb, vc))                         q -= penalty;
                if (n == 4 && on_one_line(vc, cyc[(k + 2) % n], va)) q -= penalty;
                if (q > best_q) { best_q = q; best_k = k; }
            }
            if (best_k == n) throw std::runtime_error(std::format("face_cutting: ear clipping failed (face {}, size {})", f.id, cyc.size()));
            if (best_q < -4) std::println("[warn] face_cutting: uv-degenerate triangle unavoidable (face {})", f.id);
            tris.push_back({cyc[(best_k + n - 1) % n], cyc[best_k], cyc[(best_k + 1) % n]});
            cyc.erase(cyc.begin() + (long)best_k);
        }
        if (on_one_line(cyc[0], cyc[1], cyc[2]))
            std::println("[warn] face_cutting: uv-degenerate triangle unavoidable (face {})", f.id);
        tris.push_back({cyc[0], cyc[1], cyc[2]});
    }
}

inline std::unique_ptr<Hmesh> compute_embedding_cut_hmesh(
    const Hmesh& hm,
    const TmeshMut& tmm,
    const vec<bool>& seam0,
    const VecXi& matching0,
    const VecXi& singular0,
          vec<bool>& seam1,
          VecXi& matching1,
          VecXi& singular1,
    vec<HalfData>& data
) {
    std::map<int, vec<std::pair<int, int>>> cuts; // face id -> segment endpoints

    struct SegRec {
        int i0;
        int i1;
        HalfData d0;
        HalfData d1;
    };
    vec<SegRec> sgms; // flat annotation list; matched to cut halfedges at the end

    vec h_aux(hm.nH, EdgeSplits{});

    vec<Row3d> vpos;
    for (auto p: hm.pos.rowwise()) { vpos.emplace_back(p); }

    // nid -> cut-mesh vertex index: tnode identity replaces epsilon-based position
    // matching. creates the vertex and registers the edge split exactly once.
    std::unordered_map<int, int> vid_of_nid;
    std::unordered_map<int, vec<int>> vert_lines; // cut vertex -> (tquad,side) lines through it
    auto vid_of = [&](int nid) -> int {
        if (auto it = vid_of_nid.find(nid); it != vid_of_nid.end()) return it->second;
        const HmLoc& l = tmm.tnodes[nid];
        int vid;
        if (std::holds_alternative<HmLocOnV>(l)) { vid = std::get<HmLocOnV>(l).id; }  // reuse the original vertex
        else {
            vpos.emplace_back(get_ptloc_pos(hm, l));
            vid = (int)vpos.size() - 1;
            if (auto hr = half_ratio_of(hm, l)) {
                auto [h0, r] = hr.value();
                h_aux[h0.id].emplace(r, vid);
                h_aux[h0.twin().id].emplace(1 - r, vid);
            }
        }
        return vid_of_nid[nid] = vid;
    };

    for (const auto& th: tmm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids  = tmm.tedges[th.teid].nids;
        const int   n_sgs = (int)nids.size() - 1;

        // both patches bordered by this tedge see its nodes on one straight
        // rectangle side in the tutte uv; key the line as tqid * 4 + side
        const auto& tw = tmm.thalfs[th.twid];
        const int   la = th.tqid * 4 + tmm.tquads[th.tqid].side_of(th);
        const int   lb = tw.tqid * 4 + tmm.tquads[tw.tqid].side_of(tw);
        auto add_line = [&](int vid, int line) {
            auto& ls = vert_lines[vid];
            if (std::find(ls.begin(), ls.end(), line) == ls.end()) ls.push_back(line);
        };

        // boundary spacing by geometric arc length: uv lengths are unavailable for
        // edges created by collapse, and any monotone spacing is valid for tutte
        double total = 0;
        vec<double> len(n_sgs);
        for (int k = 0; k < n_sgs; ++k) {
            len[k] = (get_ptloc_pos(hm, tmm.tnodes[nids[k + 1]]) - get_ptloc_pos(hm, tmm.tnodes[nids[k]])).norm();
            total += len[k];
        }
        if (total < EPS) throw std::runtime_error("compute_embedding_cut_hmesh: zero-length tedge");

        double sum = 0;
        for (int k = 0; k < n_sgs; ++k) {
            const HmLoc& fr = tmm.tnodes[nids[k]];
            const HmLoc& to = tmm.tnodes[nids[k + 1]];
            double v0  = sum / total; sum += len[k];
            double v1  = sum / total;
            double v0i = 1 - v0;
            double v1i = 1 - v1;
            int i0 = vid_of(nids[k]);
            int i1 = vid_of(nids[k + 1]);
            for (int i: {i0, i1}) { add_line(i, la); add_line(i, lb); }
            int l  = n_sgs - k - 1;
            // on-edge segments cut nothing: the edge subdivision already realizes them
            auto eo = try_get_edge(hm, fr, to);
            auto fo = try_get_face(hm, fr, to);
            if (!eo.has_value()) cuts[fo.value().id].emplace_back(i0, i1);
            sgms.push_back({
                .i0 = i0, .i1 = i1,
                .d0 = HalfData{.half = Half(), .v0 = v0,  .v1 = v1,  .thid = th.id,   .tqid = th.tqid,                  .twin = -1, .order = k},
                .d1 = HalfData{.half = Half(), .v0 = v1i, .v1 = v0i, .thid = th.twid, .tqid = tmm.thalfs[th.twid].tqid, .twin = -1, .order = l}
            });
        }
    }

    vec<std::array<int, 3>> tris;
    tris.reserve(hm.nF * 2);
    for (Face f: hm.faces) { face_cutting(f, cuts[f.id], vpos, h_aux, vert_lines, tris); }

    MatXi face_info(tris.size(), 3);
    MatXd vert_info(vpos.size(), 3);
    for (int i = 0; i < vpos.size(); i++) { vert_info.row(i) = vpos[i]; }
    for (int i = 0; i < tris.size(); i++) { face_info.row(i) << tris[i][0], tris[i][1], tris[i][2]; }

    auto hm_cut = std::make_unique<Hmesh>(vert_info, face_info);
    seam1 = std::vector(hm_cut->nE, false);
    matching1 = VecXi::Zero(hm_cut->nE);

    struct EdgeKey {
        int tail;
        int head;
        bool operator==(const EdgeKey& o) const noexcept { return tail == o.tail && head == o.head; }
    };

    struct EdgeKeyHash {
        std::size_t operator()(const EdgeKey& k) const noexcept {
            return (static_cast<std::size_t>(k.tail) << 32) ^ static_cast<std::size_t>(k.head);
        }
    };

    std::unordered_map<EdgeKey, Half, EdgeKeyHash> half_by_verts;
    half_by_verts.reserve(hm_cut->nH * 2);

    for (Half h: hm_cut->halfs) { half_by_verts.insert({EdgeKey{.tail=h.tail().id, .head=h.head().id}, h}); }

    // propagate seam flags: walk the split chain along each seam edge of the original mesh
    for (Edge e: hm.edges) {
        if (!seam0[e.id]) continue;
        Half h = e.half();
        vec chain = { h.tail().id };
        for (auto it = h_aux[h.id].rbegin(); it != h_aux[h.id].rend(); ++it) chain.push_back(it->second);  // tail -> head (descending r)
        chain.push_back(h.head().id);
        for (size_t k = 0; k + 1 < chain.size(); ++k) {
            Half hh = half_by_verts.at({chain[k], chain[k + 1]});   // same direction as h
            seam1[hh.edge().id] = true;
            matching1(hh.edge().id) = hh.isCanonical() ? matching0(e.id) : -matching0(e.id);
        }
        //for (size_t k = 0; k + 1 < chain.size(); ++k)
        //    seam1[half_by_verts.at({chain[k], chain[k + 1]}).edge().id] = true;
    }

    // update singular
    singular1 = VecXi::Zero(hm_cut->nV);
    singular1.head(hm.nV) = singular0;

    // annotate the cut halfedges: both directions pushed adjacently, twin = index ^ 1
    data.clear();
    data.reserve(sgms.size() * 2);
    for (auto& [i0, i1, d0, d1]: sgms) {
        auto it = half_by_verts.find({.tail = i0, .head = i1});
        if (it == half_by_verts.end()) {
            std::println("[cut] no halfedge {} -> {} (thid {}, order {}, dist {})", i0, i1, d0.thid, d0.order, (vpos[i0] - vpos[i1]).norm());
            continue; // TEMP: skip to collect all offenders
        }
        Half h0 = it->second;

        d0.half = h0;        d0.twin = (int)data.size() + 1;
        d1.half = h0.twin(); d1.twin = (int)data.size();
        data.push_back(d0);
        data.push_back(d1);
    }

    return hm_cut;
}
}
#endif
