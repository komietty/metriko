//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// quad patching: build the quad mesh from the extracted q-faces and replay the
// t-mesh on it, so that every quad knows which tquad (patch) it belongs to.
#ifndef METRIKO_PATCHING_H
#define METRIKO_PATCHING_H
#include <set>
#include <map>
#include <queue>
#include <Eigen/Geometry>
#include <igl/AABB.h>
#include <igl/per_face_normals.h>
#include <igl/remove_duplicate_vertices.h>
#include "qex/common.h"
#include "tmesh/emesh.h"

namespace metriko {

//------------------------------------------------------------------------------
// quad mesh from the q-faces
//------------------------------------------------------------------------------

// quad soup of the extracted q-faces: 4 corners per quad, corner (i, j) is qfaces[i].qhalfs[j].port1()
inline void quad_soup(const vec<qex::Qface>& qfaces, MatXd& pos, MatXi& idx) {
    const int l = (int)qfaces.size();
    pos.resize(l * 4, 3);
    idx.resize(l, 4);
    for (int i = 0; i < l; i++)
    for (int j = 0; j < 4; j++) {
        pos.row(i * 4 + j) = qfaces[i].qhalfs[j].port1().pos;
        idx(i, j) = i * 4 + j;
    }
}

// weld the per-quad duplicated corners into a connected quad mesh
inline void weld_quad_soup(const MatXd& pos, const MatXi& idx, MatXd& pos_w, MatXi& idx_w) {
    VecXi SVI, SVJ;
    const double eps = 1e-7 * (pos.colwise().maxCoeff() - pos.colwise().minCoeff()).norm();
    igl::remove_duplicate_vertices(pos, eps, pos_w, SVI, SVJ);
    idx_w = idx.unaryExpr([&](int i) { return SVJ(i); });
}

// conformal relaxation of a quad mesh on the input surface: shape-up (bouaziz et al. 2012) with the set of squares as
// the shape constraint. local step: every quad is projected onto its closest square, a similarity of the unit square
// fitted to its four corners (umeyama 1991). global step: every vertex takes the least-squares compromise of the
// corners its quads ask for, which for this constraint set alone is their mean. the surface closeness is enforced as
// a projection onto the input surface after every global step (the limit of an infinite weight on that constraint),
// so the vertices never leave it. a square has right angles and equal sides, and its size is free per quad: the
// iteration drives every corner towards 90 degrees while the quad size varies smoothly over the mesh.
// not part of shape-up: a step must not fold a quad. a quad convex with respect to the surface normal (summed over the
// faces closest to its corners) before a step keeps its corners where they were if the step would make it non-convex;
// reverting corners can fold a neighbour in turn, so this repeats until nothing changes. quads that came out of the
// quad extraction folded move freely, and are guarded once they come out convex
inline void conformalize_on_surface(MatXd& pos, const MatXi& idx, const MatXd& V, const MatXi& F, const int max_iter = 200, const double tol = 1e-6) {
    Eigen::Matrix<double, 3, 4> square;   // the unit square, corners in the cyclic order of a quad
    square << -1,  1, 1, -1,
              -1, -1, 1,  1,
               0,  0, 0,  0;
    VecXi n_quads = VecXi::Zero(pos.rows());
    for (int i = 0; i < idx.rows(); ++i) for (int j = 0; j < 4; ++j) ++n_quads(idx(i, j));

    MatXd N;
    igl::per_face_normals(V, F, N);
    igl::AABB<MatXd, 3> tree;
    tree.init(V, F);
    VecXd sqrD; VecXi I; MatXd C;
    tree.squared_distance(V, F, pos, sqrD, I, C);

    auto convex = [&](const MatXd& P, const VecXi& face_of, int q) {
        Row3d n = Row3d::Zero();
        for (int j = 0; j < 4; ++j) n += N.row(face_of(idx(q, j)));
        for (int j = 0; j < 4; ++j) {
            const Row3d c = P.row(idx(q, j)), a = P.row(idx(q, (j + 3) % 4)), b = P.row(idx(q, (j + 1) % 4));
            if ((b - c).cross(a - c).dot(n) <= 0) return false;
        }
        return true;
    };
    vec<char> guarded(idx.rows());
    for (int q = 0; q < idx.rows(); ++q) guarded[q] = convex(pos, I, q);

    vec<Eigen::Matrix<double, 3, 4>> fit(idx.rows());
    for (int it = 0; it < max_iter; ++it) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < idx.rows(); ++i) {
            Eigen::Matrix<double, 3, 4> p;
            for (int j = 0; j < 4; ++j) p.col(j) = pos.row(idx(i, j)).transpose();
            const Eigen::Matrix4d T = Eigen::umeyama(square, p, true);
            fit[i] = (T.topLeftCorner<3, 3>() * square).colwise() + T.topRightCorner<3, 1>();
        }
        MatXd next = MatXd::Zero(pos.rows(), 3);
        for (int i = 0; i < idx.rows(); ++i) for (int j = 0; j < 4; ++j) next.row(idx(i, j)) += fit[i].col(j).transpose();
        for (int v = 0; v < pos.rows(); ++v) next.row(v) = n_quads(v) > 0 ? Row3d(next.row(v) / n_quads(v)) : Row3d(pos.row(v));

        VecXi I_next;
        tree.squared_distance(V, F, next, sqrD, I_next, C);
        for (bool changed = true; changed; ) {
            changed = false;
            for (int q = 0; q < idx.rows(); ++q) {
                if (!guarded[q] || convex(C, I_next, q)) continue;
                for (int j = 0; j < 4; ++j) {
                    const int v = idx(q, j);
                    if (C.row(v) == pos.row(v)) continue;
                    C.row(v) = pos.row(v);
                    I_next(v) = I(v);
                    changed = true;
                }
            }
        }
        for (int q = 0; q < idx.rows(); ++q) if (!guarded[q]) guarded[q] = convex(C, I_next, q);

        const double move = (C - pos).rowwise().norm().maxCoeff();   // in grid units: the pipeline runs on V / grid_unit
        pos = C;
        I   = I_next;
        if (move < tol) break;
    }
}

// welded (and optionally smoothed) quad mesh of the q-faces. face (i, j) still corresponds to qfaces[i].qhalfs[j].port1()
inline std::pair<MatXd, MatXi> extract_quad_mesh(const Hmesh& hm, const vec<qex::Qface>& qfaces, const bool refine = true) {
    MatXd pos, pos_w;
    MatXi idx, idx_w;
    quad_soup(qfaces, pos, idx);
    weld_quad_soup(pos, idx, pos_w, idx_w);
    if (refine) conformalize_on_surface(pos_w, idx_w, hm.pos, hm.idx);
    return {pos_w, idx_w};
}

//------------------------------------------------------------------------------
// patch labeling: replay the t-mesh on the quad mesh
//------------------------------------------------------------------------------

// the quad mesh as directed edges
struct QuadGraph {
    const MatXi& idx;
    std::map<std::pair<int, int>, std::pair<int, int>> corner;   // directed edge a -> b: (quad on its left, corner a)
    vec<vec<int>> nbrs;                                           // quad vertex -> neighbours
    umap<int, int> qv_of_vid;                                     // mesh vertex -> quad vertex (q-vertices on a vertex)

    QuadGraph(const vec<qex::Qface>& qfaces, const MatXi& idx): idx(idx), nbrs(idx.size() ? idx.maxCoeff() + 1 : 0) {
        for (int i = 0; i < idx.rows(); ++i)
        for (int j = 0; j < 4; ++j) {
            int a = idx(i, j), b = idx(i, (j + 1) % 4);
            corner[{a, b}] = {i, j};
            nbrs[a].push_back(b);
            if (int v = qfaces[i].qhalfs[j].port1().vid; v >= 0) qv_of_vid[v] = a;
        }
    }
    int quad(int a, int b)     const { return corner.at({a, b}).first; }                                // quad on the left of a -> b
    int valence(int v)         const { return (int)nbrs[v].size(); }
    int rot(int v, int w)      const { auto [q, j] = corner.at({v, w}); return idx(q, (j + 3) % 4); }   // neighbour of v after w, CCW
    int straight(int t, int v) const { return rot(v, rot(v, t)); }                                      // arrived t -> v: keep going
    vec<int> ring(int v, int first) const {                                                             // neighbours of v, CCW from `first`
        vec<int> r = {first};
        while ((int)r.size() < valence(v)) r.push_back(rot(v, r.back()));
        return r;
    }
};

// the t-mesh as seen from the quad mesh, built once and never changed
struct TmeshView {
    struct Tedge  { int cano = -1; int steps = 0; int left = -1; int right = -1; };   // canonical thalf, quantized length, tquads along nids
    struct Branch { int teid; bool fwd; int slot; };   // a tedge leaving a node, forward along its nids or not, and its slot offset around the node

    const Emesh& tm;
    const VecXi& singular;
    const QuadGraph& g;
    vec<Tedge> te;                                       // per tedge
    umap<int, vec<std::pair<int, bool>>> branches;       // node -> (teid, leaves the node forward)

    TmeshView(const Emesh& tm, const VecXi& singular, const QuadGraph& g): tm(tm), singular(singular), g(g), te(tm.tedges.size()) {
        for (auto& th: tm.thalfs) if (th.id != -1 && th.cano) te[th.teid] = {th.id, (int)std::round(th.x), th.tqid, th.twin().tqid};
        for (auto& [teid, nids]: tm.live_tedges()) {
            branches[nids.front()].emplace_back(teid, true);
            branches[nids.back()].emplace_back(teid, false);
        }
    }
    int far_node(int teid, bool fwd) const { auto& nids = tm.tedges[teid].nids; return fwd ? nids.back() : nids.front(); }
    const HmLocOnV* singular_vert(int nid) const {   // the singular vertex a node sits on, or null
        auto* lv = std::get_if<HmLocOnV>(&tm.tnodes[nid]);
        return lv && singular(lv->id) ? lv : nullptr;
    }
    int singular_qv(int nid) const {                 // quad vertex of a singular node, -1 otherwise
        auto* lv = singular_vert(nid);
        return lv && g.qv_of_vid.contains(lv->id) ? g.qv_of_vid.at(lv->id) : -1;
    }
    // all tedges at a node, CCW from (t, fwd), each with its slot offset in quad edges. rotating CCW around the node
    // from an outgoing tedge sweeps through the tquad P on its left and reaches P's previous boundary thalf, which
    // ends at the node. this follows P's boundary order, so a tquad touching itself at the node is handled; the slot
    // advances by 2 when the node lies inside a side of P
    vec<Branch> chain(int nid, int t, bool fwd) const {
        vec<Branch> ch = {{t, fwd, 0}};
        for (int k = 1, slot = 0; k < (int)branches.at(nid).size(); ++k) {
            const int c    = te[t].cano;
            const auto& th = fwd ? tm.thalfs[c] : tm.thalfs[c].twin(); // P on its left, pointing away from the node
            const auto& tq = tm.tquads[th.tqid];
            auto cur = rg::find(tq.data, th.id, &Edata::thid);
            auto prv = circular_prev(tq.data, cur);
            const auto& th2 = tm.thalfs[prv->thid];                     // ends at the node
            slot += prv->side == cur->side ? 2 : 1;
            t   = th2.teid;
            fwd = !th2.cano;                                            // outgoing from the node: backwards along nids if cano
            ch.push_back({t, fwd, slot});
        }
        return ch;
    }
};

struct QuadPatch {
    vec<double> tqid_of_quad;                  // per quad, -1 when unlabeled
    std::set<std::pair<int, int>> track;       // quad edges (min, max) lying on a tedge
    umap<int, int> node_qv;                    // t-node -> quad vertex
    // singular nodes whose rotation was fixed / with no consistent rotation / with several, tedges never walked, quads with no tquad
    int anchors = 0, no_rotation = 0, several = 0, unreached = 0, unlabeled = 0;
    bool ok() const { return no_rotation == 0 && several == 0 && unreached == 0 && unlabeled == 0; }
};

// walk every tedge on the quad mesh, labeling the quads on its two sides. the state of the walk only: a replay is
// copied to try the rotations at a singular
struct Replay {
    struct Job { int teid; bool fwd; int from; int to; };   // walk teid starting with the quad edge from -> to

    const TmeshView* tv;                   // pointer, so that a replay can be copied and assigned
    QuadPatch out;                         // the labels so far
    vec<bool> done;                        // per tedge
    std::queue<Job> jobs;
    bool ok = true;

    explicit Replay(const TmeshView& tv): tv(&tv), done(tv.tm.tedges.size(), false) { out.tqid_of_quad.assign(tv.g.idx.rows(), -1); }

    // anchor a node at a quad vertex and queue the tedges leaving it, the first at slot `rot`
    void anchor(int nid, int qv, const vec<TmeshView::Branch>& branches, const vec<int>& ring, int rot) {
        out.node_qv[nid] = qv;
        for (auto [t, fwd, slot]: branches) if (!done[t]) jobs.push({t, fwd, qv, ring[(slot + rot) % ring.size()]});
    }

    // walk a tedge straight along the quad edges; false at the first contradiction
    bool walk(const Job& jb) {
        auto label = [&](int q, int tqid) { auto& l = out.tqid_of_quad[q]; if (l >= 0 && l != tqid) return false; l = tqid; return true; };
        const auto& e = tv->te[jb.teid];
        const int lq = jb.fwd ? e.left : e.right, rq = jb.fwd ? e.right : e.left;
        int a = jb.from, b = jb.to;
        for (int k = 0; k < e.steps; ++k) {
            if (!label(tv->g.quad(a, b), lq) || !label(tv->g.quad(b, a), rq)) return false;
            out.track.insert(std::minmax(a, b));
            if (k + 1 == e.steps) break;
            if (tv->g.valence(b) != 4) return false;   // irregular vertex before the far node
            std::tie(a, b) = std::pair(b, tv->g.straight(a, b));
        }
        const int nid = tv->far_node(jb.teid, jb.fwd);
        if (int s = tv->singular_qv(nid); s >= 0 && s != b) return false;               // landed away from the singular
        if (auto it = out.node_qv.find(nid); it != out.node_qv.end()) return it->second == b;   // anchored elsewhere
        anchor(nid, b, tv->chain(nid, jb.teid, !jb.fwd), tv->g.ring(b, a), 0);         // slots CCW from the reverse of arrival
        return true;
    }

    void run() {
        while (ok && !jobs.empty()) {
            Job jb = jobs.front(); jobs.pop();
            if (done[jb.teid]) continue;
            done[jb.teid] = true;
            ok = walk(jb);
        }
    }
};

// replay the t-mesh on the quad mesh and label every quad by its tquad, with no geometry.
//
// singular vertices are q-vertices, so their quad vertex is exact. every tedge walks exactly x quad edges straight
// ahead, and its landing vertex IS its far node, which anchors the tedges there in CCW order. the only unknown is
// the rotation of the chain at the first singular of each component: every rotation is replayed and the one
// without contradictions is kept. the labels are then flooded into the patch interiors without crossing a track
inline QuadPatch label_quad_patches(
    const Emesh& tm,
    const VecXi& singular,
    const vec<qex::Qface>& qfaces,
    const MatXi& qidx   // welded quad mesh faces, corner (i, j) <-> qfaces[i].qhalfs[j].port1()
) {
    const QuadGraph g(qfaces, qidx);
    const TmeshView tv(tm, singular, g);
    Replay rp(tv);
    for (auto& [nid, brs]: tv.branches) {
        const int s = tv.singular_qv(nid);
        if (s < 0 || rp.out.node_qv.contains(nid)) continue;
        const auto ring  = g.ring(s, g.nbrs[s].front());
        const auto chain = tv.chain(nid, brs[0].first, brs[0].second);
        vec<Replay> good;
        for (int rot = 0; rot < ring.size(); ++rot) {
            Replay t = rp;
            t.anchor(nid, s, chain, ring, rot);
            t.run();
            if (t.ok) good.push_back(std::move(t));
        }
        if (good.empty()) { ++rp.out.no_rotation; continue; }
        const bool several = good.size() > 1;
        rp = std::move(good.front());
        rp.out.several += several;
        ++rp.out.anchors;
    }
    for (auto& [teid, nids]: tm.live_tedges()) if (tv.te[teid].steps > 0 && !rp.done[teid]) ++rp.out.unreached;

    QuadPatch res = std::move(rp.out);
    std::queue<int> que;   // flood
    for (int i = 0; i < (int)res.tqid_of_quad.size(); ++i) if (res.tqid_of_quad[i] >= 0) que.push(i);
    while (!que.empty()) {
        int q = que.front(); que.pop();
        for (int j = 0; j < 4; ++j) {
            int a = qidx(q, j), b = qidx(q, (j + 1) % 4);
            if (res.track.contains(std::minmax(a, b))) continue;
            if (int nb = g.quad(b, a); res.tqid_of_quad[nb] < 0) { res.tqid_of_quad[nb] = res.tqid_of_quad[q]; que.push(nb); }
        }
    }
    res.unlabeled = rg::count(res.tqid_of_quad, -1.);
    return res;
}

}
#endif
