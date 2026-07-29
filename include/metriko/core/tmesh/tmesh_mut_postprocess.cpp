#include <set>
#include "./tmesh_mut.h"

namespace metriko {
bool TmeshMut::collapse_valid_snap(Vert v) {

    for (auto& [id, nids]: tedges) {
        if (id == -1) continue;
        if (rg::any_of(nids, [&](int nid) {
            auto* l = std::get_if<HmLocOnV>(&tnodes[nid]);
            return l && l->id == v.id;
        })) return false;
    }



    return true;
}

bool TmeshMut::collapse_valid() {

    auto faces_of = [&](const HmLoc& l, vec<int>& out) {
        std::visit(overloaded{
            [&](const HmLocOnV& v) { for (Half h: hm.verts[v.id].adjHalfs()) out.push_back(h.face().id); },
            //[&](const HmLocOnE& e) { out.push_back(hm.edges[e.id].face0().id); out.push_back(hm.edges[e.id].face1().id); },
            //[&](const HmLocOnH& h) { Half hh = hm.halfs[h.id]; out.push_back(hh.face().id); out.push_back(hh.twin().face().id); },
            //[&](const HmLocOnF& f) { out.push_back(f.id); },
            [&](const auto&)       {},
        }, l);
    };

    bool ok = true;
    for (auto& tq: tquads) {
        if (tq.id == -1) continue;
        for (int side = 0; side < 4; side++) {
            // tnodes in same side must not pass three vertices of the same face
            umap<int, int> count;
            std::set<int>  seen;
            for (int thid: tq.thids(side)) {
            for (int nid: tedges[thalfs[thid].teid].nids) {
                if (!seen.insert(nid).second) continue;
                vec<int> fids;
                faces_of(tnodes[nid], fids);
                for (int fid: fids)
                    if (++count[fid] == 3) {
                        std::println("[collapse_valid] tquad {} side {}: face {} holds 3+ nodes", tq.id, side, fid);
                        ok = false;
                    }
            }}
        }
    }
    return ok;

    for (auto& [id, nids]: tedges) {
        if (id == -1) continue;
    }

    /*
    // 1: two segments from two adjacent tedge shold not be too close to each other
    // 2: three nodes from same tedge should not belong to on the vertices of same triangle
    // report tedge pairs leaving a shared node (singularity / junction) at an angle
    // narrower than min_angle: such wedges become slivers in the tutte cutting
    constexpr double min_angle = PI / 12.;   // 15 deg
    struct Port { int teid; bool at_front; double ang; };
    umap<int, vec<Port>> node_ports;
    for (const auto& [id, nids]: tedges) {
        if (id == -1) continue;
        node_ports[nids.front()].push_back({.teid = id, .at_front = true,  .ang = 0});
        node_ports[nids.back() ].push_back({.teid = id, .at_front = false, .ang = 0});
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
    */
}

static bool is_in_face(Face face, const HmLoc& l) {
    return std::visit(overloaded{
        [&](const HmLocOnV& v) { for (Half h_: face.adjHalfs()) { if (h_.tail().id == v.id) return true; } return false; },
        [&](const HmLocOnE& e) { for (Half h_: face.adjHalfs()) { if (h_.edge().id == e.id) return true; } return false; },
        [&](const HmLocOnH& h) { for (Half h_: face.adjHalfs()) { if (h_.id == h.id)        return true; } return false; },
        [&](const HmLocOnF& f) { return f.id == face.id; },
        [&](const auto&) -> bool { throw std::runtime_error("not implemented"); },
    }, l);
};

void TmeshMut::collapse_tedge_0(int teid) {
    int last_nid = tedges[teid].nids.back();
    vec<std::pair<int, int>> te_tails; // teid, the latest nid
    vec<std::pair<int, int>> te_heads; // teid, the smallest nid

    for (auto& [id, nids]: tedges) {
        if (id == -1) continue;
        if (last_nid == nids.front()) te_tails.emplace_back(id, 0);
        if (last_nid == nids.back())  te_heads.emplace_back(id, nids.size());
    }

    if (auto* l = std::get_if<HmLocOnF>(&tnodes[last_nid])) {
        Vert  v = {};
        Face  f = hm.faces[l->id];
        Row3d p = get_ptloc_pos(hm, *l);
        auto  d = 1e6;

        for (auto h: f.adjHalfs()) {
            Vert v1 = h.tail();
            auto d1 = (v1.pos() - p).norm();
            if (d1 < d && collapse_valid_snap(v1)) { d = d1; v = v1; }
        }

        tnodes[last_nid] = HmLocOnV{.id = v.id};

        for (auto h: v.adjHalfs()) {
            Face f1 = h.face();
            for (auto& [i_, j_]: te_tails) { auto& nids = tedges[i_].nids; for (int k = 0; k < nids.size(); k++) { if (is_in_face(f1, tnodes[nids[k]])) j_ = std::max(j_, k); } }
            for (auto& [i_, j_]: te_heads) { auto& nids = tedges[i_].nids; for (int k = 0; k < nids.size(); k++) { if (is_in_face(f1, tnodes[nids[k]])) j_ = std::min(j_, k); } }
        }

        for (auto& [i_, j_]: te_tails) { auto& nids = tedges[i_].nids; if (j_ >= 2)              nids.erase(nids.begin() + 1, nids.begin() + j_);   }
        for (auto& [i_, j_]: te_heads) { auto& nids = tedges[i_].nids; if (j_ + 2 < nids.size()) nids.erase(nids.begin() + j_ + 1, nids.end() - 1); }
    }
}

void TmeshMut::collapse_tedge_1() {
    Vert v_min = {};
    auto d_min = 1e6;
    auto i_min = -1;
    auto e_min = -1;

    auto cb = [&](const HmLoc& l, int teid,  int iter, Vert v) {
        if (!collapse_valid_snap(v) || !collapse_valid()) return;
        Row3d p = get_ptloc_pos(hm, l);
        auto  d = (v.pos() - p).norm();
        if (d < d_min) {
            e_min = teid;
            i_min = iter;
            v_min = v;
            d_min = d;
        }
    };

    for (auto& [teid, nids]: tedges) {
        if (teid == -1) continue;
        for (int nid: nids) {
            std::visit(overloaded{
                [&](const HmLocOnE& l) {
                    cb(l, teid, nid, hm.edges[l.id].vert0());
                    cb(l, teid, nid, hm.edges[l.id].vert1());
                },
                [&](const HmLocOnH& l) {
                    cb(l, teid, nid, hm.halfs[l.id].edge().vert0());
                    cb(l, teid, nid, hm.halfs[l.id].edge().vert1());
                },
                [&](const HmLocOnF& l) {
                    cb(l, teid, nid, hm.faces[l.id].half().tail());
                    cb(l, teid, nid, hm.faces[l.id].half().head());
                    cb(l, teid, nid, hm.faces[l.id].half().crnr().vert());
                },
                [&](const auto&) {},
            }, tnodes[nid]);
        }
    }

    if (i_min != -1) {
        auto& [teid, nids] = tedges[e_min];
        int ii = rg::find(nids, i_min) - nids.begin();

        auto in_ring = [&](int k) {
            for (auto h: v_min.adjHalfs()) if (is_in_face(h.face(), tnodes[nids[k]])) return true;
            return false;
        };
        int k_min = ii;
        int k_max = ii;
        while (k_min > 0               && in_ring(k_min - 1)) --k_min;
        while (k_max + 1 < nids.size() && in_ring(k_max + 1)) ++k_max;

        if (k_max - ii >= 2) nids.erase(nids.begin() + ii + 1,    nids.begin() + k_max);
        if (ii - k_min >= 2) nids.erase(nids.begin() + k_min + 1, nids.begin() + ii);
        tnodes[i_min] = HmLocOnV{.id = v_min.id};
    } else {
        std::println("not found");
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
        [&](const auto&) {},
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
