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
#include <igl/AABB.h>
#include <igl/remove_duplicate_vertices.h>
#include "core/qex/common.h"
#include "core/tmesh/emesh.h"

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
    idx_w.resize(idx.rows(), idx.cols());
    for (int i = 0; i < idx.rows(); ++i)
    for (int j = 0; j < idx.cols(); ++j) idx_w(i, j) = SVJ(idx(i, j));
}

// laplacian smoothing of a quad mesh, projected back onto the input surface after every step
inline void smooth_on_surface(MatXd& pos, const MatXi& idx, const MatXd& V, const MatXi& F, const int iters = 300, const double lambda = 0.01) {
    vec<std::set<int>> adj(pos.rows());
    for (int i = 0; i < idx.rows(); ++i)
    for (int j = 0; j < 4; ++j) {
        int a = idx(i, j), b = idx(i, (j + 1) % 4);
        adj[a].insert(b);
        adj[b].insert(a);
    }
    igl::AABB<MatXd, 3> tree;
    tree.init(V, F);
    for (int it = 0; it < iters; ++it) {
        MatXd next = pos;
        for (int v = 0; v < pos.rows(); ++v) {
            if (adj[v].empty()) continue;
            Row3d c = Row3d::Zero();
            for (int n: adj[v]) c += pos.row(n);
            next.row(v) = (1 - lambda) * pos.row(v) + lambda * c / (double)adj[v].size();
        }
        VecXd sqrD; VecXi I; MatXd C;
        tree.squared_distance(V, F, next, sqrD, I, C);
        pos = C;
    }
}

// welded (and optionally smoothed) quad mesh of the q-faces. face (i, j) still corresponds to qfaces[i].qhalfs[j].port1()
inline std::pair<MatXd, MatXi> extract_quad_mesh(const Hmesh& hm, const vec<qex::Qface>& qfaces, const bool refine = true) {
    MatXd pos, pos_w;
    MatXi idx, idx_w;
    quad_soup(qfaces, pos, idx);
    weld_quad_soup(pos, idx, pos_w, idx_w);
    if (refine) smooth_on_surface(pos_w, idx_w, hm.pos, hm.idx);
    return {pos_w, idx_w};
}

//------------------------------------------------------------------------------
// navigation on the quad mesh
//------------------------------------------------------------------------------

struct QuadGraph {
    const MatXi& idx;
    std::map<std::pair<int, int>, std::pair<int, int>> corner;   // directed edge a -> b: (quad on its left, corner a)
    umap<int, vec<int>> nbrs;                                     // vertex -> neighbours
    umap<int, int> qv_of_vid;                                     // mesh vertex -> quad vertex (q-vertices on a vertex)

    QuadGraph(const vec<qex::Qface>& qfaces, const MatXi& idx): idx(idx) {
        for (int i = 0; i < idx.rows(); ++i)
        for (int j = 0; j < 4; ++j) {
            int a = idx(i, j), b = idx(i, (j + 1) % 4);
            corner[{a, b}] = {i, j};
            nbrs[a].push_back(b);
            if (int v = qfaces[i].qhalfs[j].port1().vid; v >= 0) qv_of_vid[v] = a;
        }
    }
    int quad(int a, int b)     const { return corner.at({a, b}).first; }                        // quad on the left of a -> b
    int valence(int v)         const { return (int)nbrs.at(v).size(); }
    int rot(int v, int w)      const { auto [q, j] = corner.at({v, w}); return idx(q, (j + 3) % 4); }   // neighbour of v after w, CCW
    int straight(int t, int v) const { return rot(v, rot(v, t)); }                              // arrived t -> v: keep going
    vec<int> ring(int v, int first) const {                                                     // neighbours of v, CCW from `first`
        vec<int> r = {first};
        while ((int)r.size() < valence(v)) r.push_back(rot(v, r.back()));
        return r;
    }
};

//------------------------------------------------------------------------------
// the t-mesh as seen from the quad mesh
//------------------------------------------------------------------------------

struct Branch { int teid; bool fwd; int slot; };   // a tedge leaving a node, forward along its nids or not, and its slot offset around the node

struct TmeshView {
    const Emesh& tm;
    const VecXi& singular;
    const QuadGraph& g;
    umap<int, int> steps;                                // teid -> quantized length
    umap<int, int> cano;                                 // teid -> canonical thalf
    umap<int, std::pair<int, int>> sides;                // teid -> (left, right) tquad along nids
    umap<int, vec<std::pair<int, bool>>> branches;       // node -> (teid, leaves the node forward)

    TmeshView(const Emesh& tm, const VecXi& singular, const QuadGraph& g): tm(tm), singular(singular), g(g) {
        for (auto& th: tm.thalfs) {
            if (th.id == -1 || !th.cano) continue;
            steps[th.teid] = (int)std::round(th.x);
            cano[th.teid]  = th.id;
            sides[th.teid] = {th.tqid, tm.thalfs[th.twid].tqid};
        }
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
    // all tedges at a node, CCW from (te, fwd), each with its slot offset in quad edges. rotating CCW around the
    // node from an outgoing tedge sweeps through the tquad P on its left and reaches P's previous boundary thalf,
    // which ends at the node. this follows P's boundary order, so a tquad touching itself at the node is handled;
    // the slot advances by 2 when the node lies inside a side of P
    vec<Branch> chain(int nid, int te, bool fwd) const {
        vec<Branch> ch = {{te, fwd, 0}};
        for (int k = 1, slot = 0; k < (int)branches.at(nid).size(); ++k) {
            const int c = cano.at(te);
            const auto& th = tm.thalfs[fwd ? c : tm.thalfs[c].twid];   // P on its left, pointing away from the node
            const auto& tq = tm.tquads[th.tqid];
            auto cur = rg::find(tq.data, th.id, &Edata::thid);
            auto prv = circular_prev(tq.data, cur);
            const auto& th2 = tm.thalfs[prv->thid];                    // ends at the node
            slot += prv->side == cur->side ? 2 : 1;
            te  = th2.teid;
            fwd = !th2.cano;                                           // outgoing from the node: backwards along nids if cano
            ch.push_back({te, fwd, slot});
        }
        return ch;
    }
};

//------------------------------------------------------------------------------
// replay: walk every tedge on the quad mesh, labeling the quads on its two sides
//------------------------------------------------------------------------------

struct Replay {
    struct Job { int teid; bool fwd; int from; int to; };   // walk teid starting with the quad edge from -> to

    const TmeshView* tv;   // pointer, so that a replay can be copied and assigned
    vec<double> tqid_of_quad;
    std::set<std::pair<int, int>> track;   // quad edges (min, max) lying on a tedge
    vec<bool> done;                        // per tedge
    umap<int, int> node_qv;                // node -> quad vertex
    std::queue<Job> jobs;
    bool ok = true;

    explicit Replay(const TmeshView& tv): tv(&tv), tqid_of_quad(tv.g.idx.rows(), -1), done(tv.tm.tedges.size(), false) {}

    // anchor a node at a quad vertex and queue the tedges leaving it, `first` at slot `rot`
    void anchor(int nid, int qv, const vec<Branch>& branches, const vec<int>& ring, int rot) {
        node_qv[nid] = qv;
        for (auto [te, fwd, slot]: branches) if (!done[te]) jobs.push({te, fwd, qv, ring[(slot + rot) % ring.size()]});
    }

    // label one quad: false when it already carries another patch
    bool label(int q, int tqid) {
        if (tqid_of_quad[q] >= 0 && tqid_of_quad[q] != tqid) return false;
        tqid_of_quad[q] = tqid;
        return true;
    }

    // walk a tedge straight along the quad edges; false at the first contradiction
    bool walk(const Job& jb) {
        auto [lq, rq] = tv->sides.at(jb.teid);
        if (!jb.fwd) std::swap(lq, rq);
        int a = jb.from, b = jb.to;
        for (int k = 0, n = tv->steps.at(jb.teid); k < n; ++k) {
            if (!label(tv->g.quad(a, b), lq) || !label(tv->g.quad(b, a), rq)) return false;
            track.insert(std::minmax(a, b));
            if (k + 1 == n) break;
            if (tv->g.valence(b) != 4) return false;   // irregular vertex before the far node
            std::tie(a, b) = std::pair(b, tv->g.straight(a, b));
        }
        const int nid = tv->far_node(jb.teid, jb.fwd);
        if (int s = tv->singular_qv(nid); s >= 0 && s != b) return false;   // landed away from the singular
        if (auto it = node_qv.find(nid); it != node_qv.end()) return it->second == b;   // anchored elsewhere
        anchor(nid, b, tv->chain(nid, jb.teid, !jb.fwd), tv->g.ring(b, a), 0);       // slots CCW from the reverse of arrival
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

//------------------------------------------------------------------------------
// result
//------------------------------------------------------------------------------

struct QuadPatch {
    struct Singular { int vid; int qv; int tedges; int tracks; int valence; };   // vertex, quad vertex, tedges at it, quad edges on a track, quad valence
    vec<double> tqid_of_quad;                  // per quad, -1 when unlabeled
    vec<double> colour;                        // per quad in [0, 1), neighbouring patches far apart; -1 when unlabeled
    std::set<std::pair<int, int>> track;       // quad edges (min, max) lying on a tedge
    umap<int, int> node_qv;                    // t-node -> quad vertex
    vec<Singular> singulars;
    vec<int> junctions;                        // quad vertices of the anchored non-singular nodes
    int anchors = 0;                           // singular nodes whose rotation was fixed
    int no_rotation = 0;                       // singular nodes with no consistent rotation
    int several = 0;                           // singular nodes with several consistent rotations
    int unreached = 0;                         // tedges never walked
    int unlabeled = 0;                         // quads with no tquad
    bool ok() const { return no_rotation == 0 && several == 0 && unreached == 0 && unlabeled == 0; }
};

// spread the labels into the patch interiors without crossing a track
inline void flood_patches(const QuadGraph& g, const std::set<std::pair<int, int>>& track, vec<double>& tqid_of_quad) {
    std::queue<int> que;
    for (int i = 0; i < (int)tqid_of_quad.size(); ++i) if (tqid_of_quad[i] >= 0) que.push(i);
    while (!que.empty()) {
        int q = que.front(); que.pop();
        for (int j = 0; j < 4; ++j) {
            int a = g.idx(q, j), b = g.idx(q, (j + 1) % 4);
            if (track.contains(std::minmax(a, b))) continue;
            int nb = g.quad(b, a);
            if (tqid_of_quad[nb] < 0) { tqid_of_quad[nb] = tqid_of_quad[q]; que.push(nb); }
        }
    }
}

// greedy colouring of the patch adjacency, then the colour index spread by the golden ratio so that neighbours map far apart
inline vec<double> colour_patches(const QuadGraph& g, const vec<double>& tqid_of_quad) {
    std::map<int, std::set<int>> adj;
    for (int i = 0; i < (int)tqid_of_quad.size(); ++i)
    for (int j = 0; j < 4; ++j) {
        int ta = (int)tqid_of_quad[i], tb = (int)tqid_of_quad[g.quad(g.idx(i, (j + 1) % 4), g.idx(i, j))];
        if (ta >= 0 && tb >= 0 && ta != tb) { adj[ta].insert(tb); adj[tb].insert(ta); }
    }
    vec<std::pair<int, int>> order;   // (degree, tqid), most constrained first
    for (auto& [t, ns]: adj) order.emplace_back((int)ns.size(), t);
    rg::sort(order, std::greater<>{});
    umap<int, int> colour_of;
    for (auto [deg, t]: order) {
        std::set<int> used;
        for (int u: adj[t]) if (colour_of.contains(u)) used.insert(colour_of.at(u));
        int c = 0;
        while (used.contains(c)) ++c;
        colour_of[t] = c;
    }
    vec<double> colour(tqid_of_quad.size(), -1);
    for (int i = 0; i < (int)tqid_of_quad.size(); ++i)
        if (int t = (int)tqid_of_quad[i]; t >= 0 && colour_of.contains(t)) colour[i] = std::fmod(colour_of.at(t) * 0.618033988749895, 1.);
    return colour;
}

// replay the t-mesh on the quad mesh and label every quad by its tquad, with no geometry.
//
// singular vertices are q-vertices, so their quad vertex is exact. every tedge walks exactly x quad edges straight
// ahead, and its landing vertex IS its far node, which anchors the tedges there in CCW order. the only unknown is
// the rotation of the chain at the first singular of each component: every rotation is replayed and the one
// without contradictions is kept.
inline QuadPatch label_quad_patches(
    const Emesh& tm,
    const VecXi& singular,
    const vec<qex::Qface>& qfaces,
    const MatXi& qidx   // welded quad mesh faces, corner (i, j) <-> qfaces[i].qhalfs[j].port1()
) {
    const QuadGraph g(qfaces, qidx);
    const TmeshView tv(tm, singular, g);
    QuadPatch res;

    Replay rp(tv);
    for (auto& [nid, brs]: tv.branches) {
        const int s = tv.singular_qv(nid);
        if (s < 0 || rp.node_qv.contains(nid)) continue;
        const auto ring  = g.ring(s, g.nbrs.at(s).front());
        const auto chain = tv.chain(nid, brs[0].first, brs[0].second);
        vec<Replay> good;
        for (int rot = 0; rot < (int)ring.size(); ++rot) {
            Replay t = rp;
            t.anchor(nid, s, chain, ring, rot);
            t.run();
            if (t.ok) good.push_back(std::move(t));
        }
        if (good.empty()) { ++res.no_rotation; continue; }
        if (good.size() > 1) ++res.several;
        rp = std::move(good.front());
        ++res.anchors;
    }
    for (auto& [teid, nids]: tm.live_tedges()) if (tv.steps.at(teid) > 0 && !rp.done[teid]) ++res.unreached;

    res.tqid_of_quad = std::move(rp.tqid_of_quad);
    res.track        = std::move(rp.track);
    res.node_qv      = std::move(rp.node_qv);
    flood_patches(g, res.track, res.tqid_of_quad);
    res.unlabeled = (int)rg::count(res.tqid_of_quad, -1.);
    res.colour    = colour_patches(g, res.tqid_of_quad);

    for (auto& [nid, v]: res.node_qv) {
        auto* lv = tv.singular_vert(nid);
        if (!lv) { res.junctions.push_back(v); continue; }
        int ntrack = 0;
        for (int w: g.nbrs.at(v)) if (res.track.contains(std::minmax(v, w))) ++ntrack;
        res.singulars.push_back({lv->id, v, (int)tv.branches.at(nid).size(), ntrack, g.valence(v)});
    }
    return res;
}

}
#endif
