
#include "./tmesh_mut.h"

namespace metriko {
void TmeshMut::collapse_valid() {
    // 1: two segments from two adjacent tedge shold not be too close to each other
    // 2: three nodes from same tedge should not belong to on the vertices of same triangle
    // report tedge pairs leaving a shared node (singularity / junction) at an angle
    // narrower than min_angle: such wedges become slivers in the tutte cutting
    constexpr double min_angle = PI / 12.;   // 15 deg

    struct Port { int teid; bool at_front; double ang; };

    umap<int, vec<Port>> node_ports;
    for (const auto& th: thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids = tedges[th.teid].nids;
        node_ports[nids.front()].push_back({.teid = th.teid, .at_front = true,  .ang = 0});
        node_ports[nids.back() ].push_back({.teid = th.teid, .at_front = false, .ang = 0});
    }

    for (auto& [nid, ports]: node_ports) {
        if (ports.size() < 2) continue;

        // tangent-plane basis at the node
        Row3d o   = get_ptloc_pos(hm, tnodes[nid]);
        Row3d nrm = get_ptloc_normal(hm, tnodes[nid]);
        Row3d t   = std::abs(nrm.x()) < 0.9 ? Row3d(1, 0, 0) : Row3d(0, 1, 0);
        Row3d u   = nrm.cross(t).normalized();
        Row3d w   = nrm.cross(u);

        // outgoing direction of the first segment of each incident tedge
        for (auto& p: ports) {
            const auto& nids = tedges[p.teid].nids;
            int nb = p.at_front ? nids[1] : nids[nids.size() - 2];
            Row3d d = get_ptloc_pos(hm, tnodes[nb]) - o;
            p.ang = std::atan2(d.dot(w), d.dot(u));
        }
        rg::sort(ports, {}, &Port::ang);

        // adjacent angular gaps (with wrap-around)
        for (size_t i = 0; i < ports.size(); ++i) {
            const auto& a = ports[i];
            const auto& b = ports[(i + 1) % ports.size()];
            double gap = b.ang - a.ang;
            if (i + 1 == ports.size()) gap += TwoPI;
            if (gap < min_angle)
                std::println("[collapse_valid] node {} {}: tedge {} / {} meet at {:.2f} deg", nid, loc_str(tnodes[nid]), a.teid, b.teid, gap * 180. / PI);
        }
    }
}

void TmeshMut::collapse_tedge_edge_snapping(int teid) {
    constexpr double delta = 0.05;
    const auto& nids = tedges[teid].nids;

    const int nid = nids.back();
    Row3d p = get_ptloc_pos(hm, tnodes[nid]);

    std::visit(overloaded{
        [&](const HmLocOnF& f) {
            for (Half h: hm.faces[f.id].adjHalfs()) {
                Row3d a  = h.tail().pos();
                Row3d d1 = h.vec();
                Row3d d2 = p - a;
                double t = d2.dot(d1) / d1.squaredNorm();
                if ((d2 - t * d1).norm() > delta * d1.norm()) continue;
                double r = h.isCanonical() ? t : 1 - t;
                tnodes[nid] = HmLocOnE{.id = h.edge().id, .r = r};
                return;
            }
        },
        [&](const auto&) {},   // OnV (singularity etc.): nothing to do
    }, tnodes[nid]);
}

void TmeshMut::collapse_tedge_vert_snapping(int teid) {
    constexpr double delta = 0.05;   // snap when within 10% of the carrier edge length
    auto& nids = tedges[teid].nids;

    // shared nodes (junctions / crossings / other tedges' endpoints) must not move
    vec shared(tnodes.size(), false);
    for (int i = 0; i < tedges.size(); ++i) {
        if (i == teid) continue;
        for (int nid: tedges[i].nids) shared[nid] = true;
    }

    // vertices already occupied by any tnode: snapping onto them would create an
    // accidental junction between unrelated tedges
    vec v_used(hm.nV, false);
    for (const auto& te: tedges)
        for (int nid: te.nids)
            if (auto* v = std::get_if<HmLocOnV>(&tnodes[nid])) v_used[v->id] = true;

    /// 1: snap exclusive interior nodes onto a nearby mesh vertex
    for (size_t i = 1; i + 1 < nids.size(); ++i) {
        int nid = nids[i];
        if (shared[nid]) continue;

        int vid = -1;
        std::visit(overloaded{
            [&](const HmLocOnE& e) {
                Edge ed = hm.edges[e.id];
                if      (e.r     < delta) vid = ed.vert0().id;
                else if (1 - e.r < delta) vid = ed.vert1().id;
            },
            [&](const HmLocOnH& hh) {
                Half h = hm.halfs[hh.id];
                if      (hh.r     < delta) vid = h.tail().id;
                else if (1 - hh.r < delta) vid = h.head().id;
            },
            [&](const HmLocOnF& f) {
                Row3d p = get_ptloc_pos(hm, tnodes[nid]);
                for (Half h: hm.faces[f.id].adjHalfs())
                    if ((h.tail().pos() - p).norm() < delta * h.vec().norm()) { vid = h.tail().id; break; }
            },
            [&](const auto&) {},   // OnV: nothing to do
        }, tnodes[nid]);
        if (vid < 0) continue;

        // allow the snap only if the vertex is free, or occupied by our direct chain
        // neighbor (that case merges in the dedup pass below)
        auto at_v = [&](int n) { auto* v = std::get_if<HmLocOnV>(&tnodes[n]); return v && v->id == vid; };
        if (v_used[vid] && !at_v(nids[i - 1]) && !at_v(nids[i + 1])) continue;

        tnodes[nid] = HmLocOnV{vid};
        v_used[vid] = true;
    }

    /// 2: dedup — snapping can land neighboring nodes on the same vertex
    for (size_t i = 0; i + 1 < nids.size();) {
        auto* a = std::get_if<HmLocOnV>(&tnodes[nids[i]]);
        auto* b = std::get_if<HmLocOnV>(&tnodes[nids[i + 1]]);
        bool dup = nids[i] == nids[i + 1] || (a && b && a->id == b->id);
        if (!dup) { ++i; continue; }
        bool can_drop_b = !shared[nids[i + 1]] && i + 1 < nids.size() - 1;   // interior, non-shared
        bool can_drop_a = !shared[nids[i]]     && i > 0;
        if      (can_drop_b) nids.erase(nids.begin() + (long)i + 1);   // stay at i:
        else if (can_drop_a) nids.erase(nids.begin() + (long)i);      // 3+ nodes may coincide
        else ++i;   // both are endpoints/shared: leave for validation to flag
    }
}

void TmeshMut::collapse_tedge_short_segment(int teid) {
    auto& nids = tedges[teid].nids;
    vec shared(tnodes.size(), false);

    for (int i = 0; i < tedges.size(); ++i) {
        if (i == teid) continue;
        for (int nid: tedges[i].nids) shared[nid] = true;
    }

    auto faces_of = [&](const HmLoc& l, vec<int>& out) {
        std::visit(overloaded{
            [&](const HmLocOnV& v) { for (Half h : hm.verts[v.id].adjHalfs()) out.push_back(h.face().id); },
            [&](const HmLocOnE& e) { out.push_back(hm.edges[e.id].face0().id); out.push_back(hm.edges[e.id].face1().id); },
            [&](const HmLocOnH& h) { Half hh = hm.halfs[h.id]; out.push_back(hh.face().id); out.push_back(hh.twin().face().id); },
            [&](const HmLocOnF& f) { out.push_back(f.id); },
            [&](const auto&)       {},
        }, l);
    };

    // drop an interior node when both incident segments lie in one common face:
    // the face is planar and convex, so the straightened chord stays inside it
    for (size_t i = 1; i + 1 < nids.size();) {
        int n1 = nids[i];
        if (shared[n1]) { ++i; continue; }

        vec<int> f0, f1, f2;
        faces_of(tnodes[nids[i - 1]], f0);
        faces_of(tnodes[n1],          f1);
        faces_of(tnodes[nids[i + 1]], f2);

        bool same_face = false;
        for (int a: f0)
        for (int b: f1)
        for (int c: f2) if (a == b && b == c) { same_face = true; goto done; }
        done:;

        if (same_face) { nids.erase(nids.begin() + i); }
        else ++i;
    }
}

// simplification objective:
// 1: straighter is better
// 2: fewer nodes are better
// 3: snapped to vertex or edge is better not to create sliver with tutte cutting
}
