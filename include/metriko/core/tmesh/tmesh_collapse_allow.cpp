#include "./tmesh.h"
#include <queue>
using namespace metriko;

/*
namespace metriko {

// Extract the "allowed ranges" of hmesh edges that lie inside a Tquad.
// Returns: eid -> Row2d(r0, r1)  (the `allowed` format for approx_shortest_path;
//          r is measured along edge.half(), i.e. tail->head).
//
// The Tquad boundary (thalf -> segments) is CCW-oriented, so the interior is
// always on the left of the travel direction. Using this we clip each edge's
// inside sub-range directly, instead of relying on whole-face classification
// (a Tquad is simply connected, so one interval per edge is enough; this also
//  works for very thin Tquads that have no fully-interior face).
//
//   1.  Crossed edges: classify each boundary OnEdge crossing as inside-on-the
//       -head-side / tail-side by the sign of the 2D cross product, then clip.
//   1b. Fill uncrossed edges trapped between two boundary ("wall") faces.
//   2.  Fully-interior edges: split faces by the wall band into components and
//       pick the interior one by CCW voting; emit its edges as [0,1].
umap<int, Row2d> allowed_ranges_in_tquad(
    const Tquad& tq,
    const Tmesh& tm,
    const mc::Mgrph& mg
) {
    const Hmesh& hm = tm.hm;
    const VecXc& cf = mg.cf;

    // Canonical (tail->head) position of edge `eid` at parameter r, in face `fid`'s uv.
    auto edge_uv = [&](int eid, double r, int fid) -> complex {
        for (Half h : hm.faces[fid].adjHalfs())
            if (h.edge().id == eid) {
                complex uv0 = cf[h.next().crnr().id];
                complex uv1 = cf[h.prev().crnr().id];
                return h.isCanonical() ? lerp(uv0, uv1, r) : lerp(uv1, uv0, r);
            }
        return {};
    };

    // --- 1. Inside sub-range of crossed edges (CCW left = inside) ---
    vec wall_faces(hm.nF, false); // faces the boundary passes through
    vec is_crossed(hm.nE, false); //
    umap<int, double> lo, hi;     // inside interval [lo, hi] per crossed edge (intersection of half-ranges)

    auto handle_crossing = [&](int nid, complex d, int fid) {
        const mc::Mnode& mn = mg.mnodes[nid];
        if (auto* e = std::get_if<mc::OnEdge>(&mn.loc)) {
            complex edge_dir = edge_uv(e->eid, 1.0, fid) - edge_uv(e->eid, 0.0, fid);
            bool toward_head = cross(d, edge_dir) > 0;   // is the left (inside) toward the head vertex?
            if (!is_crossed[e->eid]) { is_crossed[e->eid] = true; lo[e->eid] = 0.0; hi[e->eid] = 1.0; }
            if (toward_head) lo[e->eid] = std::max(lo[e->eid], e->r);   // inside subset of [r, 1]
            else             hi[e->eid] = std::min(hi[e->eid], e->r);   // inside subset of [0, r]
        }
    };

    for (const auto& [thid, _] : tq.data) {
        const Thalf& th = tm.thalfs[thid];
        const Tedge& te = th.edge();
        for (const mc::Msgmt& sg : te.segs) {
            if (sg.face_id != -1) wall_faces[sg.face_id] = true;
            // Directed segment aligned with the CCW travel direction (reversed for non-canonical thalf).
            int fr = th.cano ? sg.fr_nid : sg.to_nid;
            int to = th.cano ? sg.to_nid : sg.fr_nid;
            complex fr_uv = mc::get_face_uv(mg.mnodes[fr], sg.face_id, hm, cf);
            complex to_uv = mc::get_face_uv(mg.mnodes[to], sg.face_id, hm, cf);
            complex d = to_uv - fr_uv;
            handle_crossing(fr, d, sg.face_id);
            handle_crossing(to, d, sg.face_id);
        }
    }

    umap<int, Row2d> allowed;
    for (int eid = 0; eid < hm.nE; ++eid)
        if (is_crossed[eid] && lo[eid] < hi[eid])
            allowed[eid] = Row2d(lo[eid], hi[eid]);

    // --- 1b. Fill uncrossed edges trapped inside the wall band (between two wall faces) ---
    for (const Tdata& td : tq.data) {
        const Thalf& th = tm.thalfs[td.thid];
        for (const mc::Msgmt& sg : th.edge().segs) {
            int fid = sg.face_id;
            if (fid == -1) continue;
            int fr = th.cano ? sg.fr_nid : sg.to_nid;
            int to = th.cano ? sg.to_nid : sg.fr_nid;
            complex fr_uv = mc::get_face_uv(mg.mnodes[fr], fid, hm, cf);
            complex d = mc::get_face_uv(mg.mnodes[to], fid, hm, cf) - fr_uv;
            for (Half h : hm.faces[fid].adjHalfs()) {
                int eid = h.edge().id;
                if (is_crossed[eid] || allowed.count(eid)) continue;
                // Only edges interior to the wall band (both incident faces are wall faces);
                // the band's outer-boundary edges are intentionally excluded.
                Edge e = hm.edges[eid];
                int f0 = e.face0().id, f1 = e.face1().id;
                if (f0 == -1 || f1 == -1 || !wall_faces[f0] || !wall_faces[f1]) continue;
                complex mid = (edge_uv(eid, 0.0, fid) + edge_uv(eid, 1.0, fid)) * 0.5;
                if (cross(d, mid - fr_uv) > 0) allowed[eid] = Row2d(0.0, 1.0);  // left = inside
            }
        }
    }

    // --- 2. Fully-interior edges: components split by the wall band, interior chosen by CCW vote ---
    // Flood-fill the non-wall faces into connected components.
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

    // Vote for each component as interior (left, +1) or exterior (right, -1) of the boundary chords.
    vec<int> votes(ncomp, 0);
    for (const Tdata& td : tq.data) {
        const Thalf& th = tm.thalfs[td.thid];
        for (const mc::Msgmt& sg : th.edge().segs) {
            int fid = sg.face_id;
            if (fid == -1) continue;
            int fr = th.cano ? sg.fr_nid : sg.to_nid;
            int to = th.cano ? sg.to_nid : sg.fr_nid;
            complex fr_uv = mc::get_face_uv(mg.mnodes[fr], fid, hm, cf);
            complex d = mc::get_face_uv(mg.mnodes[to], fid, hm, cf) - fr_uv;
            for (Half h : hm.faces[fid].adjHalfs()) {
                int nf = h.twin().face().id;
                if (nf == -1 || wall_faces[nf]) continue;
                complex mid = (edge_uv(h.edge().id, 0.0, fid) + edge_uv(h.edge().id, 1.0, fid)) * 0.5;
                votes[comp[nf]] += (cross(d, mid - fr_uv) > 0) ? 1 : -1;   // left = inside
            }
        }
    }

    // The interior component is the one with the most "inside" votes.
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
*/
