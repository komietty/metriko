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

// a transition between two charts of the integer-grid map: uv' = rot * uv + off, rot a power of i
struct Xf {
    complex rot = 1;
    complex off = 0;
    complex operator()(complex uv) const { return rot * uv + off; }
    Xf operator*(const Xf& x) const { return {rot * x.rot, rot * x.off + off}; }   // this after x
    Xf inv() const { return {1. / rot, -off / rot}; }
};

// a face placed in the chart: uv in the chart = x(uv in the face)
struct Pface {
    Face f;
    Xf   x = {};

    complex uv(Vert v) const { return x(try_get_crnr(*f.m, v.id, f.id).value().uv()); }
    complex uv(const HmLoc& loc) const {
        return std::visit(overloaded{
            [&](const HmLocOnV& l) { return uv(f.m->verts[l.id]); },
            [&](const HmLocOnE& l) { Edge e = f.m->edges[l.id]; return lerp(uv(e.vert0()), uv(e.vert1()), l.r); },
            [&](const HmLocOnF& l) {
                auto [a, b, c] = f.verts();
                const complex w = calc_coefficient(f.to_local(a.pos()), f.to_local(b.pos()), f.to_local(c.pos()), l.xy);
                return uv(a) + (uv(b) - uv(a)) * w.real() + (uv(c) - uv(a)) * w.imag();
            },
            [&](const auto& l) -> complex { METRIKO_FAIL("{} has no uv", loc_str(l)); },
        }, loc);
    }
    std::array<complex, 3> corners() const { auto [a, b, c] = f.verts(); return {uv(a), uv(b), uv(c)}; }   // corner k is the tail of halfs()[k]
    double depth(complex q) const {   // how far inside the face (negative outside)
        auto [a, b, c] = corners();
        return std::min({side_distance(a, b, q), side_distance(b, c, q), side_distance(c, a, q)});
    }

    // the face across the edge e, placed so that the two agree along e
    Pface across(Edge e) const {
        const Face g  = e.face0() == f ? e.face1() : e.face0();
        const Pface p = {g};
        const complex d   = (uv(e.vert1()) - uv(e.vert0())) / (p.uv(e.vert1()) - p.uv(e.vert0()));
        const complex rot = std::polar(1., PI / 2 * std::lround(std::arg(d) / (PI / 2)));
        return {g, {rot, uv(e.vert0()) - rot * p.uv(e.vert0())}};
    }
    // the faces around the vertex v, counter-clockwise from this one up to stop (all of them when stop is f)
    vec<Pface> around(Vert v, Face stop) const {
        vec<Pface> res = {*this};
        const Half h = try_get_crnr(*f.m, v.id, f.id).value().half().prev();   // the half of f leaving v
        for (Half g: v.adjHalfs(h)) {   // g.edge() leads from the previous face into g.face()
            if (g.face() == f) continue;
            res.push_back(res.back().across(g.edge()));
            if (g.face() == stop) break;
        }
        return res;
    }
};

// the faces a line may run into from a node of p: both sides of an edge node, the fan of a vertex node up to stop
inline vec<Pface> fan(const Hmesh& hm, const HmLoc& loc, const Pface& p, Face stop) {
    if (auto* e = std::get_if<HmLocOnE>(&loc)) return {p, p.across(hm.edges[e->id])};
    if (auto* v = std::get_if<HmLocOnV>(&loc)) return p.around(hm.verts[v->id], stop);
    return {p};
}

// the boundary of a region as a polygon in one chart. the region lies on the left of the boundary
struct Chart {
    vec<int>     nids;    // the boundary nodes in order
    vec<complex> uvs;     // their uv in the chart
    vec<Pface>   faces;   // faces[k] carries the segment nids[k] -> nids[k + 1]
    vec<int>     heads;   // the k-th thalf of the boundary starts at nids[heads[k]]

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
    auto face_of = [&](int k) {   // a face carrying segment k (either one when it runs along a mesh edge)
        auto f = try_get_face(hm, tnodes[c.nids[k]], tnodes[c.nids[c.next(k)]]);
        METRIKO_CHECK(f, "boundary segment {} -> {} lies in no face", loc_str(tnodes[c.nids[k]]), loc_str(tnodes[c.nids[c.next(k)]]));
        return *f;
    };

    // the faces are placed one after the other through the node they share. around a vertex the way goes through
    // the inside wedge, from the segment leaving it to the one entering it: a singular vertex is only met on the
    // boundary, and the other way round would pick up its holonomy
    c.faces.push_back({face_of(0)});
    for (int k = 1; k < n; ++k) {
        const Pface& p = c.faces.back();
        const Face   g = face_of(k);
        const HmLoc& l = tnodes[c.nids[k]];
        if (g == p.f)                                 c.faces.push_back(p);
        else if (auto* v = std::get_if<HmLocOnV>(&l)) c.faces.push_back({g, p.x * Pface{g}.around(hm.verts[v->id], p.f).back().x.inv()});
        else                                          c.faces.push_back(p.across(try_get_edge(hm, l, l).value()));
    }
    for (int k = 0; k < n; ++k) c.uvs.push_back(c.faces[k].uv(tnodes[c.nids[k]]));

    // a disk closes: the first node seen from the last segment is where it started
    const double gap = std::abs(c.faces[n - 1].uv(tnodes[c.nids[0]]) - c.uvs[0]);
    METRIKO_CHECK(gap < 1e-6, "the region does not unfold into one chart (closure error {})", gap);
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

    // the placed face a line from p toward q runs into, among those around the current node: the first one when a
    // step along the line stays in it, otherwise the one that step lies deepest in
    auto enter = [&](const vec<Pface>& around, complex p, complex q) {
        const complex step = p + (q - p) * 1e-6;
        if (around[0].depth(step) > 0) return around[0];
        return *rg::max_element(around, {}, [&](const Pface& pf) { return pf.depth(step); });
    };
    // a point of the placed face as a location: on a vertex or an edge when it lies there, inside otherwise
    auto locate = [&](const Pface& p, complex q) -> HmLoc {
        const auto cs = p.corners();
        const auto hs = p.f.halfs();
        for (int k = 0; k < 3; ++k) {
            if (std::abs(side_distance(cs[k], cs[(k + 1) % 3], q)) > EPS) continue;
            const double r = std::abs(q - cs[k]) / std::abs(cs[(k + 1) % 3] - cs[k]);
            if (r < EPS)     return HmLocOnV{hs[k].tail().id};
            if (r > 1 - EPS) return HmLocOnV{hs[k].head().id};
            return HmLocOnE{hs[k].edge().id, hs[k].isCanonical() ? r : 1 - r};
        }
        return HmLocOnF{p.f.id, p.f.to_local(p.f.uv2pos(p.x.inv()(q)))};
    };

    Pface p = enter(fan(hm, path.back(), c.faces[a], c.faces[c.prev(a)].f), pts[0], pts[1]);   // the inside wedge at a
    for (size_t s = 0; s + 1 < pts.size(); ++s) {
        complex u = pts[s];
        const complex q = pts[s + 1];
        if (s > 0) p = enter(fan(hm, path.back(), p, p.f), u, q);
        while (p.depth(q) < -EPS) {   // leave the face through the edge the line crosses first
            const auto cs = p.corners();
            double best = 2, r = 0; int k = -1;
            for (int j = 0; j < 3; ++j) {
                double se, tp;
                if (!find_extended_intersection(cs[j], cs[(j + 1) % 3], u, q, se, tp)) continue;
                if (tp > EPS && se > -EPS && se < 1 + EPS && tp < best) { best = tp; k = j; r = se; }
            }
            METRIKO_CHECK(k != -1, "the traced line leaves face {} through no edge", p.f.id);
            u = lerp(cs[k], cs[(k + 1) % 3], r);
            path.push_back(locate(p, u));
            p = std::holds_alternative<HmLocOnV>(path.back()) ? enter(fan(hm, path.back(), p, p.f), u, q)   // through a vertex
                                                              : p.across(p.f.halfs()[k].edge());
        }
        if (s + 2 < pts.size()) path.push_back(locate(p, q));
    }
    path.push_back(tnodes[c.nids[b]]);
    return path;
}
}
