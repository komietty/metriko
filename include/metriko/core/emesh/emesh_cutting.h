#ifndef METRIKO_TUTTE_CUTTING_H
#define METRIKO_TUTTE_CUTTING_H
#include <set>
#include <unordered_map>
#include "./emesh.h"
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/tmesh/tmesh_mut.h"

namespace metriko::emesh {

// split points on an original halfedge: (ratio in the tail-weighted convention,
// cut-mesh vertex index), ordered along the halfedge
using EdgeSplits = std::set<std::pair<double, int>>;

// a directed edge of the subdivided face; exists only while cutting one face
struct DirEdge { int i0; int i1; };

// the face shared by two locations (a segment always lies inside one face)
inline int common_face(const Hmesh& hm, const HmLoc& a, const HmLoc& b) {
    auto faces_of = [&](const HmLoc& l, vec<int>& out) {
        std::visit(overloaded{
            [&](const HmLocOnV& v) { for (Half h : hm.verts[v.id].adjHalfs()) out.push_back(h.face().id); },
            [&](const HmLocOnE& e) { out.push_back(hm.edges[e.id].face0().id); out.push_back(hm.edges[e.id].face1().id); },
            [&](const HmLocOnH& h) { Half hh = hm.halfs[h.id]; out.push_back(hh.face().id); out.push_back(hh.twin().face().id); },
            [&](const HmLocOnF& f) { out.push_back(f.id); },
            [&](const auto&)       { throw std::runtime_error("common_face: unsupported loc"); },
        }, l);
    };
    vec<int> fa, fb;
    faces_of(a, fa); faces_of(b, fb);
    for (int x : fa) for (int y : fb) if (x == y) return x;
    throw std::runtime_error("common_face: endpoints share no face");
}

// (half, r) in the legacy convention (r weights the TAIL) when the point sits on an edge
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

    /// 2: trace the pieces by the leftmost-turn rule (max signed CCW angle).
    ///    reflex turns are legal, and every chain must close.
    while (!halfs.empty()) {
        DirEdge h = halfs.back();
        halfs.pop_back();

        vec<DirEdge> poly = {h};
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
            if (best_it == halfs.end())
                throw std::runtime_error(std::format("face_cutting: open chain (face {})", f.id));

            poly.push_back(*best_it);
            cur = best_it->i1;
            halfs.erase(best_it);
        }

        /// 3: ear-clip the piece in the face plane (concave-capable, keeps CCW
        ///    orientation so no flipped faces can be produced)
        vec<int> cyc;
        for (auto& e: poly) cyc.push_back(e.i0);
        if (cyc.size() < 3) continue;   // zero-area piece (e.g. a segment lying on a mesh edge)

        while (cyc.size() > 3) {
            bool clipped = false;
            const size_t n = cyc.size();
            for (size_t k = 0; k < n && !clipped; ++k) {
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

                tris.push_back({va, vb, vc});
                cyc.erase(cyc.begin() + (long)k);
                clipped = true;
            }
            if (!clipped)
                throw std::runtime_error(std::format("face_cutting: ear clipping failed (face {}, size {})", f.id, cyc.size()));
        }
        tris.push_back({cyc[0], cyc[1], cyc[2]});
    }
}

// Cut the original mesh for tutte parameterization.
// `data` receives one HalfData per segment direction; the two directions of a
// segment are adjacent, so the twin of data[i] is always data[i ^ 1] (unsorted).
inline std::unique_ptr<Hmesh> compute_embedding_cut_hmesh(
    const Hmesh& hm,         // input hmesh
    const TmeshMut& tmm,     // input tmesh (post-collapse; new traced edges included)
    const vec<bool>& seam0,  //
          vec<bool>& seam1,  //
    vec<HalfData>& data      //
) {
    std::map<int, vec<std::pair<int, int>>> cuts; // face id -> segment endpoints

    struct SegRec { int i0; int i1; HalfData d0; HalfData d1; };
    vec<SegRec> sgms; // flat annotation list; matched to cut halfedges at the end

    vec h_aux(hm.nH, EdgeSplits{});

    vec<Row3d> vpos;
    for (auto p: hm.pos.rowwise()) { vpos.emplace_back(p); }

    // nid -> cut-mesh vertex index: tnode identity replaces epsilon-based position
    // matching. creates the vertex and registers the edge split exactly once.
    std::unordered_map<int, int> vid_of_nid;
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
        const auto& [nids] = tmm.tedges[th.teid];
        const int   n_sgs  = (int)nids.size() - 1;

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
            double v0 = sum / total; sum += len[k];
            double v1 = sum / total;
            int i0 = vid_of(nids[k]);
            int i1 = vid_of(nids[k + 1]);
            int l = n_sgs - k - 1;
            cuts[common_face(hm, fr, to)].emplace_back(i0, i1);
            sgms.push_back({
                i0, i1,
                HalfData{ Half(),     v0,     v1, th.id,   th.tqid,                  -1, k },
                HalfData{ Half(), 1 - v1, 1 - v0, th.twid, tmm.thalfs[th.twid].tqid, -1, l }
            });
        }
    }

    vec<std::array<int, 3>> tris;
    tris.reserve(hm.nF * 2);
    for (Face f: hm.faces) { face_cutting(f, cuts[f.id], vpos, h_aux, tris); }

    MatXi face_info((int)tris.size(), 3);
    MatXd vert_info((int)vpos.size(), 3);
    for (int i = 0; i < (int)vpos.size(); i++) { vert_info.row(i) = vpos[i]; }
    for (int i = 0; i < (int)tris.size(); i++) { face_info.row(i) << tris[i][0], tris[i][1], tris[i][2]; }

    auto hm_cut = std::make_unique<Hmesh>(vert_info, face_info);
    seam1 = std::vector(hm_cut->nE, false);

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

    for (Half h: hm_cut->halfs) {
        half_by_verts.insert({EdgeKey{h.tail().id, h.head().id}, h});
    }

    // propagate seam flags: walk the split chain along each seam edge of the original mesh
    for (Edge e: hm.edges) {
        if (!seam0[e.id]) continue;
        Half h = e.half();
        vec<int> chain = { h.tail().id };
        for (auto it = h_aux[h.id].rbegin(); it != h_aux[h.id].rend(); ++it) chain.push_back(it->second);  // tail -> head (descending r)
        chain.push_back(h.head().id);
        for (size_t k = 0; k + 1 < chain.size(); ++k)
            seam1[half_by_verts.at({chain[k], chain[k + 1]}).edge().id] = true;
    }

    // annotate the cut halfedges: both directions pushed adjacently, twin = index ^ 1
    data.clear();
    data.reserve(sgms.size() * 2);
    for (auto& s: sgms) {
        Half h0 = half_by_verts.at({s.i0, s.i1});
        s.d0.half = h0;        s.d0.twin = (int)data.size() + 1;
        s.d1.half = h0.twin(); s.d1.twin = (int)data.size();
        data.push_back(s.d0);
        data.push_back(s.d1);
    }

    return hm_cut;
}
}
#endif
