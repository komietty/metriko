//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_HMFLAT_H
#define METRIKO_HMFLAT_H
#include "metriko/hmesh/hmloc.h"
#include "metriko/hmesh/utilities.h"

namespace metriko {

struct HmFlat {
    Face f;
    arr<complex, 3> uv;
};

inline HmFlat flat_own(Face f) {
    auto [h0, h1, h2] = f.halfs();
    return {f, {h0.cr_t().uv(), h1.cr_t().uv(), h2.cr_t().uv()}};
}

inline complex flat_uv(const HmFlat& c, Vert v) {
    auto vs = c.f.verts();
    return c.uv[rg::find(vs, v) - vs.begin()];
}

inline complex flat_uv(const HmFlat& c, const HmLoc& loc) {
    const Face f = c.f;
    return std::visit(overloaded{
        [&](const HmLocOnV& l) { return flat_uv(c, f.m->verts[l.id]); },
        [&](const HmLocOnE& l) { Edge e = f.m->edges[l.id]; return lerp(flat_uv(c, e.vert0()), flat_uv(c, e.vert1()), l.r); },
        [&](const HmLocOnF& l) {
            auto [a, b, d] = f.verts();
            auto w = calc_coefficient(f.to_local(a.pos()), f.to_local(b.pos()), f.to_local(d.pos()), l.xy);
            return c.uv[0] + (c.uv[1] - c.uv[0]) * w.real() + (c.uv[2] - c.uv[0]) * w.imag();
        },
        [&](const auto& l) -> complex { METRIKO_FAIL("{} has no uv", loc_str(l)); },
    }, loc);
}

inline HmFlat flat_across(const HmFlat& c, Edge e) {
    auto g = flat_own(e.face0() == c.f ? e.face1() : e.face0());
    auto d = (flat_uv(c, e.vert1()) - flat_uv(c, e.vert0())) / (flat_uv(g, e.vert1()) - flat_uv(g, e.vert0()));
    auto rot = std::polar(1., PI / 2 * std::lround(std::arg(d) / (PI / 2)));
    auto off = flat_uv(c, e.vert0()) - rot * flat_uv(g, e.vert0());
    for (auto& uv: g.uv) uv = rot * uv + off;
    return g;
}


inline vec<HmFlat> flat_around(const HmFlat& c, Vert v, Face stop, bool ccw = true) {
    vec res = {c};
    Half h = try_get_crnr(*c.f.m, v.id, c.f.id).value().half().prev();   // the half of f leaving v
    Half p = h;
    for (Half g: v.adjHalfs(h, ccw)) {
        if (g.face() == c.f) continue;
        res.push_back(flat_across(res.back(), (ccw ? g : p).edge()));
        p = g;
        if (g.face() == stop) break;
    }
    return res;
}

inline HmLoc flat_loc(const HmFlat& c, complex q) {
    auto hs = c.f.halfs();
    for (int k = 0; k < 3; ++k) {
        auto a = c.uv[k], b = c.uv[(k + 1) % 3];
        if (std::abs(side_distance(a, b, q)) > EPS) continue;
        auto r = std::abs(q - a) / std::abs(b - a);
        if (r < EPS)     return HmLocOnV{hs[k].tail().id};
        if (r > 1 - EPS) return HmLocOnV{hs[k].head().id};
        return HmLocOnE{hs[k].edge().id, hs[k].isCanonical() ? r : 1 - r};
    }
    auto [a, b, d] = c.f.verts();
    return HmLocOnF{c.f.id, c.f.to_local(conversion_2d_3d(c.uv[0], c.uv[1], c.uv[2], a.pos(), b.pos(), d.pos(), q))};
}

inline vec<HmLoc> flat_trace(HmFlat& c, const HmLoc& from, complex q, Face stop) {

    auto dist = [](const HmFlat& c, complex q) {
        return std::min({
            side_distance(c.uv[0], c.uv[1], q),
            side_distance(c.uv[1], c.uv[2], q),
            side_distance(c.uv[2], c.uv[0], q)
        });
    };

    auto exit = [&](complex u) {
        int k = -1;
        double r = 0;
        double t = 2;
        double se;
        double tp;
        for (int j = 0; j < 3; ++j) {
            if (!find_extended_intersection(c.uv[j], c.uv[(j + 1) % 3], u, q, se, tp)) continue;
            if (tp > EPS && se > -EPS && se < 1 + EPS && tp < t) { t = tp; k = j; r = se; }
        }
        return std::pair{k, r};
    };

    auto enter = [&](const HmLoc& loc, complex u, Face stop) {
        vec fan = {c};
        if (auto* e = std::get_if<HmLocOnE>(&loc)) fan.push_back(flat_across(c, c.f.m->edges[e->id]));
        if (auto* v = std::get_if<HmLocOnV>(&loc)) fan = flat_around(c, c.f.m->verts[v->id], stop);
        auto step = u + (q - u) * 1e-6;
        if (dist(fan[0], step) > 0) return fan[0];
        return *rg::max_element(fan, {}, [&](const HmFlat& c) { return dist(c, step); });
    };

    vec<HmLoc> path;
    auto u = flat_uv(c, from);
    c = enter(from, u, stop);
    while (dist(c, q) < -EPS) {
        auto [k, r] = exit(u);
        u = lerp(c.uv[k], c.uv[(k + 1) % 3], r);
        path.push_back(flat_loc(c, u));
        c = std::holds_alternative<HmLocOnV>(path.back()) ? enter(path.back(), u, c.f) : flat_across(c, c.f.halfs()[k].edge());
    }
    path.push_back(flat_loc(c, q));
    return path;
}
}
#endif
