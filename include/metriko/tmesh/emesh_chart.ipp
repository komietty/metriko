//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// a region of the t-mesh (a tquad, or two tquads sharing a tedge) is a disk, so it unfolds into one chart of the
// integer-grid map, where its boundary is a simple polygon. a new tedge between two boundary nodes is the shortest
// path inside that polygon, traced back onto the mesh face by face
#include <queue>
#include "emesh.h"
namespace metriko {

constexpr double CHART_HAIR = 1e-4;   // how far a path bend is pushed off the boundary, in grid units

// a transition of the integer-grid map between two charts: uv' = rot * uv + off, rot a power of i
struct Xf {
    complex rot = 1;
    complex off = 0;
    complex operator()(complex uv) const { return rot * uv + off; }
    Xf operator*(const Xf& x) const { return {rot * x.rot, rot * x.off + off}; }   // this after x
    Xf inv() const { return {1. / rot, -off / rot}; }
};

// uv_f = T(uv_g) for the faces f, g across the edge e
inline Xf transition(const Hmesh& hm, Face f, Face g, Edge e) {
    auto uf = [&](Vert v) { return try_get_crnr(hm, v.id, f.id).value().uv(); };
    auto ug = [&](Vert v) { return try_get_crnr(hm, v.id, g.id).value().uv(); };
    const complex d   = (uf(e.vert1()) - uf(e.vert0())) / (ug(e.vert1()) - ug(e.vert0()));
    const complex rot = std::polar(1., PI / 2 * std::lround(std::arg(d) / (PI / 2)));
    return {rot, uf(e.vert0()) - rot * ug(e.vert0())};
}

// the uv of a location in the chart of face fid
inline complex face_uv(const Hmesh& hm, const HmLoc& loc, int fid) {
    if (auto* l = std::get_if<HmLocOnF>(&loc)) {
        Face f = hm.faces[fid];
        auto [a, b, c] = f.verts();
        auto uv = [&](Vert v) { return try_get_crnr(hm, v.id, fid).value().uv(); };
        const complex w = calc_coefficient(f.to_local(a.pos()), f.to_local(b.pos()), f.to_local(c.pos()), l->xy);
        return uv(a) + (uv(b) - uv(a)) * w.real() + (uv(c) - uv(a)) * w.imag();
    }
    return get_face_uv(loc, fid, hm);
}

// the boundary of a region as a polygon in one chart. the region lies on the left of the boundary
struct Chart {
    vec<int>     nids;   // the boundary nodes in order
    vec<complex> uvs;    // their uv in the chart
    vec<int>     fids;   // the face carrying the segment nids[k] -> nids[k + 1]
    vec<Xf>      xfs;    // uv of that face -> uv of the chart
    vec<int>     heads;  // the k-th thalf of the boundary starts at nids[heads[k]]

    int n()         const { return nids.size(); }
    int next(int i) const { return (i + 1) % n(); }
    int prev(int i) const { return (i + n() - 1) % n(); }
    int at(int k, int nid) const {   // the position of a node of the k-th thalf (its ends included)
        for (int i = heads[k]; ; i = next(i)) {
            if (nids[i] == nid) return i;
            METRIKO_CHECK(i != heads[(k + 1) % heads.size()], "node {} is not on the {}-th thalf of the region", nid, k);
        }
    }
};

// the faces around a node with their transforms: counter-clockwise from f up to stop (all of them when stop is f).
// at a boundary node, from the segment leaving it to the segment entering it is the inside wedge; a singular
// vertex is only ever met there, and the other way round would pick up its holonomy
inline vec<std::pair<Face, Xf>> Emesh::fan(const HmLoc& loc, Face f, Xf xf, Face stop) const {
    vec<std::pair<Face, Xf>> res = {{f, xf}};
    if (auto* e = std::get_if<HmLocOnE>(&loc)) {
        Edge ed = hm.edges[e->id];
        Face g  = ed.face0() == f ? ed.face1() : ed.face0();
        res.emplace_back(g, xf * transition(hm, f, g, ed));
    }
    if (auto* v = std::get_if<HmLocOnV>(&loc)) {
        Half h = try_get_crnr(hm, v->id, f.id).value().half().prev();   // the half of f leaving the vertex
        for (Half g: hm.verts[v->id].adjHalfs(h)) {                      // g.edge() leads from the previous face into g.face()
            if (g.face() == f) continue;
            res.emplace_back(g.face(), res.back().second * transition(hm, res.back().first, g.face(), g.edge()));
            if (g.face() == stop) break;
        }
    }
    return res;
}

inline Chart Emesh::unfold(const vec<int>& thids) const {
    Chart c;
    for (int thid: thids) {   // the nodes of each thalf in its direction, its last node left to the next thalf
        const Ehalf& th = thalfs[thid];
        vec<int> ns = th.tedge().nids;
        if (!th.cano) rg::reverse(ns);
        c.heads.push_back(c.nids.size());
        c.nids.insert(c.nids.end(), ns.begin(), ns.end() - 1);
    }
    const int n = c.n();

    // the face of every segment. a segment along a mesh edge takes the face on its left
    for (int k = 0; k < n; ++k) {
        const HmLoc& a = tnodes[c.nids[k]];
        const HmLoc& b = tnodes[c.nids[c.next(k)]];
        if (auto e = try_get_edge(hm, a, b)) {
            auto t = [&](const HmLoc& l) { auto* v = std::get_if<HmLocOnV>(&l); return v ? (v->id == e->vert0().id ? 0. : 1.) : std::get<HmLocOnE>(l).r; };
            c.fids.push_back((t(a) < t(b) ? e->half() : e->half().twin()).face().id);
        } else {
            auto f = try_get_face(hm, a, b);
            METRIKO_CHECK(f, "boundary segment {} -> {} lies in no face", loc_str(a), loc_str(b));
            c.fids.push_back(f->id);
        }
    }

    // the transform of every face: the chart is the first face, the next ones are reached through the shared node
    c.xfs.push_back({});
    for (int k = 1; k < n; ++k) {
        const Face f = hm.faces[c.fids[k - 1]];
        const Face g = hm.faces[c.fids[k]];
        const HmLoc& l = tnodes[c.nids[k]];
        Xf x = c.xfs.back();
        if (f != g) x = std::holds_alternative<HmLocOnV>(l) ? x * fan(l, g, {}, f).back().second.inv()   // around the inside wedge, g to f
                                                            : x * transition(hm, f, g, try_get_edge(hm, l, l).value());
        c.xfs.push_back(x);
    }
    for (int k = 0; k < n; ++k) c.uvs.push_back(c.xfs[k](face_uv(hm, tnodes[c.nids[k]], c.fids[k])));

    // a disk closes: the first node seen from the last segment is where it started
    const complex back = c.xfs[n - 1](face_uv(hm, tnodes[c.nids[0]], c.fids[n - 1]));
    METRIKO_CHECK(std::abs(back - c.uvs[0]) < 1e-6, "the region does not unfold into one chart (closure error {})", std::abs(back - c.uvs[0]));
    return c;
}

// the shortest path from boundary node a to boundary node b inside the polygon, as uv points (none when there is
// no way through). it bends only at corners that are not convex, and there it is pushed off the boundary by a hair
// so that the new tedge touches no other tedge (a wall traced as a straight line is a run of such corners)
inline vec<complex> route(const Chart& c, int a, int b) {
    const int n = c.n();
    auto cross = [](complex p, complex q) { return p.real() * q.imag() - p.imag() * q.real(); };
    auto turn  = [&](int i) { return cross(c.uvs[i] - c.uvs[c.prev(i)], c.uvs[c.next(i)] - c.uvs[i]); };   // > 0 at a convex corner
    auto hair  = [&](int i) {   // into the region, along the bisector (or the left normal where the boundary is straight)
        const complex o = c.uvs[c.next(i)] - c.uvs[i], e = c.uvs[c.prev(i)] - c.uvs[i];
        complex d = -(o / std::abs(o) + e / std::abs(e));
        if (std::abs(d) < 1e-6) d = o / std::abs(o) * complex(0, 1);
        return c.uvs[i] + d / std::abs(d) * CHART_HAIR;
    };
    auto into = [&](int i, complex d) {   // d leaves corner i strictly into the region
        const complex o = c.uvs[c.next(i)] - c.uvs[i], e = c.uvs[c.prev(i)] - c.uvs[i];
        return cross(o, e) > 0 ? cross(o, d) > EPS && cross(d, e) > EPS : cross(o, d) > EPS || cross(d, e) > EPS;
    };

    // the path nodes: the two ends on the boundary, and the bends pushed inside
    vec<int>     idx = {a, b};
    vec<complex> pts = {c.uvs[a], c.uvs[b]};
    auto visible = [&](int u, int v) {
        const complex p = pts[u], q = pts[v];
        if (u < 2 && !into(idx[u], q - p)) return false;
        if (v < 2 && !into(idx[v], p - q)) return false;
        for (int k = 0; k < n; ++k) {
            if (find_strict_intersection(p, q, c.uvs[k], c.uvs[c.next(k)])) return false;
            if (k != a && k != b && is_inside_segment(p, q, c.uvs[k])) return false;
        }
        return true;
    };
    if (visible(0, 1)) return pts;
    for (int i = 0; i < n; ++i) if (i != a && i != b && turn(i) < 1e-6) { idx.push_back(i); pts.push_back(hair(i)); }

    const int m = idx.size();
    vec dist(m, std::numeric_limits<double>::infinity());
    vec from(m, -1);
    std::priority_queue<std::pair<double, int>, vec<std::pair<double, int>>, std::greater<>> pq;
    dist[0] = 0; pq.emplace(0., 0);
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        if (u == 1) break;
        for (int v = 0; v < m; ++v)
            if (const double w = d + std::abs(pts[v] - pts[u]); w < dist[v] && visible(u, v)) { dist[v] = w; from[v] = u; pq.emplace(w, v); }
    }
    if (from[1] == -1) return {};   // no way through: the region is pinched between the two nodes

    vec<complex> path;
    for (int v = 1; v != -1; v = from[v]) path.push_back(pts[v]);
    return path | vw::reverse | rg::to<vec<complex>>();
}

// the path from boundary node a to boundary node b, traced on the mesh. its ends are the nodes themselves
inline vec<HmLoc> Emesh::retrace(const Chart& c, int a, int b) const {
    const vec<complex> pts = route(c, a, b);
    if (pts.empty()) return {};
    vec<HmLoc> path = {tnodes[c.nids[a]]};

    auto corners = [&](Face f, const Xf& x) {
        std::array<Half, 3> h = f.halfs();
        std::array<complex, 3> p;
        for (int k = 0; k < 3; ++k) p[k] = x(try_get_crnr(hm, h[k].tail().id, f.id).value().uv());
        return std::pair{p, h};
    };
    auto depth = [&](Face f, const Xf& x, complex q) {   // how far inside f (negative outside)
        auto [p, h] = corners(f, x);
        return std::min({side_distance(p[0], p[1], q), side_distance(p[1], p[2], q), side_distance(p[2], p[0], q)});
    };
    // the face among those around the current node that the segment p -> q runs into: the first one when a step
    // along the segment stays in it, otherwise the one that step lies deepest in
    auto enter = [&](const vec<std::pair<Face, Xf>>& fan, complex p, complex q) {
        const complex step = p + (q - p) * 1e-6;
        if (depth(fan[0].first, fan[0].second, step) > 0) return fan[0];
        return *rg::max_element(fan, {}, [&](auto& fx) { return depth(fx.first, fx.second, step); });
    };
    // a point of f as a location: on a vertex or an edge of f when it lies there, inside otherwise
    auto locate = [&](Face f, const Xf& x, complex q) -> HmLoc {
        auto [p, h] = corners(f, x);
        for (int k = 0; k < 3; ++k) {
            if (std::abs(side_distance(p[k], p[(k + 1) % 3], q)) > EPS) continue;
            const double r = std::abs(q - p[k]) / std::abs(p[(k + 1) % 3] - p[k]);
            if (r < EPS)     return HmLocOnV{h[k].tail().id};
            if (r > 1 - EPS) return HmLocOnV{h[k].head().id};
            return HmLocOnE{h[k].edge().id, h[k].isCanonical() ? r : 1 - r};
        }
        return HmLocOnF{f.id, f.to_local(f.uv2pos(x.inv()(q)))};
    };

    const Face f0 = hm.faces[c.fids[a]];
    auto [f, x] = enter(fan(path.back(), f0, c.xfs[a], hm.faces[c.fids[c.prev(a)]]), pts[0], pts[1]);
    for (size_t s = 0; s + 1 < pts.size(); ++s) {
        complex p = pts[s];
        const complex q = pts[s + 1];
        if (s > 0) std::tie(f, x) = enter(fan(path.back(), f, x, f), p, q);
        while (depth(f, x, q) < -EPS) {   // leave f through the edge the line crosses first
            auto [cs, hs] = corners(f, x);
            double best = 2, r = 0; int k = -1;
            for (int j = 0; j < 3; ++j) {
                double se, tp;
                if (!find_extended_intersection(cs[j], cs[(j + 1) % 3], p, q, se, tp)) continue;
                if (tp > EPS && se > -EPS && se < 1 + EPS && tp < best) { best = tp; k = j; r = se; }
            }
            METRIKO_CHECK(k != -1, "the traced line leaves face {} through no edge", f.id);
            p = cs[k] + (cs[(k + 1) % 3] - cs[k]) * r;
            path.push_back(locate(f, x, p));
            if (std::holds_alternative<HmLocOnV>(path.back())) {   // through a vertex: pick the face it continues into
                std::tie(f, x) = enter(fan(path.back(), f, x, f), p, q);
            } else {
                x = x * transition(hm, f, hs[k].twin().face(), hs[k].edge());
                f = hs[k].twin().face();
            }
        }
        if (s + 2 < pts.size()) path.push_back(locate(f, x, q));
    }
    path.push_back(tnodes[c.nids[b]]);
    return path;
}
}
