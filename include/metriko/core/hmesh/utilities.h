//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//
#ifndef METRIKO_HMESH_UTILITIES_H
#define METRIKO_HMESH_UTILITIES_H
#include "../common/utilities.h"
#include "../common/predicates.h"
#include "hmesh.h"
#include "hmloc.h"

namespace metriko {
inline complex calc_coefficient(const Face f, const VecXc& cf, const complex uv) {
    auto cs = f.crnrs();
    return calc_coefficient(cf(cs[0].id), cf(cs[1].id), cf(cs[2].id), uv);
}

inline Row3d conversion_2d_3d(const Face& f, const VecXc& cf, const complex uv) {
    auto [c1, c2, c3] = f.crnrs();
    return conversion_2d_3d(
        cf(c1.id),
        cf(c2.id),
        cf(c3.id),
        c1.vert().pos(),
        c2.vert().pos(),
        c3.vert().pos(),
        uv);
}

inline bool is_inside_face(const Face f, const VecXc& cf, const complex uv) {
    auto uv1 = cf(f.id * 3 + 0);
    auto uv2 = cf(f.id * 3 + 1);
    auto uv3 = cf(f.id * 3 + 2);
    return is_inside_triangle(uv1, uv2, uv3, uv);
}

inline std::optional<Crnr> try_get_crnr(const Hmesh& hm, int vid, int fid) {
    Face f = hm.faces[fid];
    Vert v = hm.verts[vid];
    for (auto h: f.adjHalfs())
        if (h.crnr().vert() == v) return h.crnr();
    return std::nullopt;
}

inline std::optional<Half> try_get_half(const Hmesh& hm, int eid, int fid) {
    Face f = hm.faces[fid];
    Edge e = hm.edges[eid];
    for (auto h: f.adjHalfs())
        if (h.edge() == e) return h;
    return std::nullopt;
}

inline std::optional<Half> try_get_opposite_half_from_half(
    const Half h0,   // half came in
    const VecXc &cf, // corner function
    const complex o, // origin in 2d
    const complex d, // direction in 2d
    const double tol // tolerance to eliminate
) {
    for (Half h: h0.face().adjHalfs()) {
        if (h.id == h0.id) continue;
        auto a = cf(h.next().crnr().id);
        auto b = cf(h.prev().crnr().id);
        if (is_points_into(o, a, b, o + d, tol)) return h;
    }
    return std::nullopt;
}

inline std::optional<Crnr> try_get_opposite_crnr_from_half(
    const Half h0,   // half came in
    const VecXc &cf, // corner function
    const complex o, // origin in 2d
    const complex d, // direction in 2d
    const double tol // tolerance to include
) {
    Crnr c = h0.crnr();
    auto b = cf(c.id);
    auto v = b - o;
    if (abs(v) > tol && dot(v, d) > 0 && is_collinear(o, o + d, b, tol)) return c;
    return std::nullopt;
}

inline std::optional<Half> try_get_opposite_half_from_crnr(
    const Crnr c,    // crnr came in
    const VecXc &cf, // corner function
    const complex d, // direction in 2d
    const double tol // tolerance to include
) {
    auto o = cf(c.id);
    auto h = c.half();
    auto a = cf(h.next().crnr().id);
    auto b = cf(h.prev().crnr().id);
    if (is_points_into(o, a, b, o + d, tol)) return h;
    return std::nullopt;
}

inline std::optional<Crnr> try_get_opposite_crnr_from_crnr(
    const Crnr c,    // crnr came in
    const VecXc &cf, // corner function
    const complex d, // direction in 2d
    const double tol // tolerance to include
) {
    auto o = cf(c.id);
    auto h  = c.half();
    auto c1 = h.next().crnr();
    auto c2 = h.prev().crnr();
    if (is_collinear(o, o + d, cf(c1.id), tol)) return c1;
    if (is_collinear(o, o + d, cf(c2.id), tol)) return c2;
    return std::nullopt;
}

inline complex get_face_uv(const HmLoc& loc, int fid, const Hmesh& hm, const VecXc& cf) {
    return std::visit(overloaded {
        [&](const HmLocOnP& f) -> complex { return f.uv; },
        [&](const HmLocOnV& v) -> complex { return cf[try_get_crnr(hm, v.id, fid).value().id]; },
        [&](const HmLocOnE& e) -> complex {
            auto h  = try_get_half(hm, e.id, fid).value();
            auto p0 = cf[h.next().crnr().id];
            auto p1 = cf[h.prev().crnr().id];
            return lerp(p0, p1, h.isCanonical() ? e.r : 1 - e.r);
        },
        [&](const auto& _) -> complex { throw std::runtime_error("invalid arguments"); },
    }, loc);
}
}
#endif
