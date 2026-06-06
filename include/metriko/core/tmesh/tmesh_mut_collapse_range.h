#ifndef METRIKO_TMESH_MUT_COLLAPSE_RANGE_H
#define METRIKO_TMESH_MUT_COLLAPSE_RANGE_H
#include "tmesh_mut.h"
#include <queue>

namespace metriko {

// Port of allowed_ranges_in_tquad for TmeshMut.
// TmeshMut has no uv (cf), so all geometry is done in 3D: tnodes (HmLoc) -> 3D via
// get_ptloc_pos, and the CCW "left = inside" test uses the face normal:
//   left(d, v)  <=>  dot(normal_f, cross(d, v)) > 0
//
// Returns: eid -> Row2d(r0, r1)  (allowed sub-range along edge.half(), tail->head).
inline umap<int, Row2d> allowed_ranges_in_tquad(const TquadMut& tq, const TmeshMut& tm) {
    const Hmesh& hm = tm.hm;

    auto pos = [&](int nid) -> Row3d { return get_ptloc_pos(hm, tm.tnodes[nid]); };

    // incident faces of a tnode
    auto node_faces = [&](int nid) -> vec<int> {
        return std::visit(overloaded{
            [&](const HmLocOnV& v) -> vec<int> {
                vec<int> fs;
                for (Half h : hm.verts[v.id].adjHalfs()) if (h.face().id != -1) fs.push_back(h.face().id);
                return fs;
            },
            [&](const HmLocOnE& e) -> vec<int> {
                vec<int> fs; Edge ed = hm.edges[e.id];
                if (ed.face0().id != -1) fs.push_back(ed.face0().id);
                if (ed.face1().id != -1) fs.push_back(ed.face1().id);
                return fs;
            },
            [&](const HmLocOnF& f) -> vec<int> { return { f.id }; },
            [&](const auto&)       -> vec<int> { return {}; },
        }, tm.tnodes[nid]);
    };

    // the face a segment (between two consecutive tnodes) lies in = their common incident face
    auto seg_face = [&](int a, int b) -> int {
        vec<int> fa = node_faces(a), fb = node_faces(b);
        for (int x : fa) for (int y : fb) if (x == y) return x;
        return -1;
    };

    // canonical (tail->head) 3D direction / midpoint of edge eid
    auto edge_dir = [&](int eid) -> Row3d { Half h = hm.edges[eid].half(); return h.head().pos() - h.tail().pos(); };
    auto edge_mid = [&](int eid) -> Row3d { Half h = hm.edges[eid].half(); return (h.head().pos() + h.tail().pos()) * 0.5; };

    // directed node sequence of a thalf in CCW boundary order (reversed for non-canonical)
    auto thalf_seq = [&](const ThalfMut& th) -> vec<int> {
        vec<int> seq = tm.tedges[th.teid].nids;
        if (!th.cano) rg::reverse(seq);
        return seq;
    };

    // --- 1. inside sub-range of crossed edges (CCW left = inside, via face normal) ---
    vec<bool> wall_faces(hm.nF, false);
    vec<bool> is_crossed(hm.nE, false);
    umap<int, double> lo, hi;

    auto handle_crossing = [&](int nid, const Row3d& d, int fid) {
        if (auto* e = std::get_if<HmLocOnE>(&tm.tnodes[nid])) {
            Row3d n = hm.faces[fid].normal();
            bool toward_head = n.dot(d.cross(edge_dir(e->id))) > 0;   // inside toward head vertex?
            if (!is_crossed[e->id]) { is_crossed[e->id] = true; lo[e->id] = 0.0; hi[e->id] = 1.0; }
            if (toward_head) lo[e->id] = std::max(lo[e->id], e->r);   // inside subset of [r, 1]
            else             hi[e->id] = std::min(hi[e->id], e->r);   // inside subset of [0, r]
        }
    };

    for (const TdataMut& td : tq.data) {
        vec<int> seq = thalf_seq(tm.thalfs[td.thid]);
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            int fr = seq[i], to = seq[i + 1];
            int fid = seg_face(fr, to);
            if (fid == -1) continue;
            wall_faces[fid] = true;
            Row3d d = pos(to) - pos(fr);
            handle_crossing(fr, d, fid);
            handle_crossing(to, d, fid);
        }
    }

    umap<int, Row2d> allowed;
    for (int eid = 0; eid < hm.nE; ++eid)
        if (is_crossed[eid] && lo[eid] < hi[eid])
            allowed[eid] = Row2d(lo[eid], hi[eid]);

    // --- 1b. fill uncrossed edges trapped inside the wall band (both faces are wall faces) ---
    for (const TdataMut& td : tq.data) {
        vec<int> seq = thalf_seq(tm.thalfs[td.thid]);
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            int fr = seq[i], to = seq[i + 1];
            int fid = seg_face(fr, to);
            if (fid == -1) continue;
            Row3d frp = pos(fr);
            Row3d d = pos(to) - frp;
            Row3d n = hm.faces[fid].normal();
            for (Half h : hm.faces[fid].adjHalfs()) {
                int eid = h.edge().id;
                if (is_crossed[eid] || allowed.count(eid)) continue;
                Edge e = hm.edges[eid];
                int f0 = e.face0().id, f1 = e.face1().id;
                if (f0 == -1 || f1 == -1 || !wall_faces[f0] || !wall_faces[f1]) continue;
                if (n.dot(d.cross(edge_mid(eid) - frp)) > 0) allowed[eid] = Row2d(0.0, 1.0);  // left = inside
            }
        }
    }

    // --- 2. fully-interior edges: components split by wall band, interior chosen by CCW vote ---
    vec<int> comp(hm.nF, -1);
    int ncomp = 0;
    for (Face f0 : hm.faces) {
        if (wall_faces[f0.id] || comp[f0.id] != -1) continue;
        int c = ncomp++;
        std::queue<int> q; q.push(f0.id); comp[f0.id] = c;
        while (!q.empty()) {
            int fid = q.front(); q.pop();
            for (Half h : hm.faces[fid].adjHalfs()) {
                int nf = h.twin().face().id;
                if (nf != -1 && !wall_faces[nf] && comp[nf] == -1) { comp[nf] = c; q.push(nf); }
            }
        }
    }

    vec<int> votes(ncomp, 0);
    for (const TdataMut& td : tq.data) {
        vec<int> seq = thalf_seq(tm.thalfs[td.thid]);
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            int fr = seq[i], to = seq[i + 1];
            int fid = seg_face(fr, to);
            if (fid == -1) continue;
            Row3d frp = pos(fr);
            Row3d d = pos(to) - frp;
            Row3d n = hm.faces[fid].normal();
            for (Half h : hm.faces[fid].adjHalfs()) {
                int nf = h.twin().face().id;
                if (nf == -1 || wall_faces[nf]) continue;
                votes[comp[nf]] += (n.dot(d.cross(edge_mid(h.edge().id) - frp)) > 0) ? 1 : -1;  // left = inside
            }
        }
    }

    int inner = -1;
    for (int c = 0; c < ncomp; ++c) if (votes[c] > 0 && (inner == -1 || votes[c] > votes[inner])) inner = c;

    if (inner != -1)
        for (int fid = 0; fid < hm.nF; ++fid) {
            if (wall_faces[fid] || comp[fid] != inner) continue;
            for (Half h : hm.faces[fid].adjHalfs())
                if (!is_crossed[h.edge().id]) allowed[h.edge().id] = Row2d(0.0, 1.0);
        }

    return allowed;
}

}
#endif
