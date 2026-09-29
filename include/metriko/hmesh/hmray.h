//
// Created by saki on 2026/09/28.
//

#ifndef HMESH_H_HMRAY_H
#define HMESH_H_HMRAY_H
#include "metriko/common/predicates.h"
#include "metriko/hmesh/hmesh.h"
#include "metriko/hmesh/hmloc.h"
namespace metriko {
inline std::optional<HmLoc> find_ray_intersection(
    const Half h,
    const complex o,
    const complex d,
    const double tol
) {
    auto a = h.cr_t().uv();
    auto b = h.cr_h().uv();
    auto s = 0.;
    auto t = 0.;
    if (!find_extended_intersection(a, b, o, o + d, s, t)) return std::nullopt;
    if (t <= tol || s < -tol || s > 1 + tol) return std::nullopt;
    if (s < tol)     return HmLocOnC{.id = h.cr_t().id};
    if (s > 1 - tol) return HmLocOnC{.id = h.cr_h().id};
    return HmLocOnH{.id = h.id, .r = s};
}

inline HmLoc find_ray_intersection(
    const Hmesh& hm,
    const HmLoc& loc,
    const complex dir
) {
    double tol = EPS;

    return std::visit(overloaded{
        [&](const HmLocOnC& l) -> HmLoc {
            if (auto x = find_ray_intersection(hm.crnrs[l.id].half(), hm.crnrs[l.id].uv(), dir, tol)) return *x;
            METRIKO_FAIL("no intersection");
        },
        [&](const HmLocOnH& l) -> HmLoc {
            auto h = hm.halfs[l.id];
            if (auto x = find_ray_intersection(h.next(), h.uv(l.r), dir, tol)) return *x;
            if (auto x = find_ray_intersection(h.prev(), h.uv(l.r), dir, tol)) return *x;
            METRIKO_FAIL("no intersection");
        },
        [&](const HmLocOnP& l) -> HmLoc {
            auto hs = hm.faces[l.id].halfs();
            if (auto x = find_ray_intersection(hs[0], l.uv, dir, tol)) return *x;
            if (auto x = find_ray_intersection(hs[1], l.uv, dir, tol)) return *x;
            if (auto x = find_ray_intersection(hs[2], l.uv, dir, tol)) return *x;
            METRIKO_FAIL("no intersection");
        },
        [&](const auto& _) -> HmLoc { METRIKO_FAIL("no impl"); },
    }, loc);
}

inline std::pair<HmLoc, complex> cross_to_twin(
    const Hmesh& hm,
    const HmLoc& loc,
    const VecXi& matching,
    const complex dir
) {
    auto get_m = [&](Half h) { return (h.isCanonical() ? -1 : 1) * matching[h.edge().id]; };

    return std::visit(overloaded{
        [&](const HmLocOnC& l) -> std::pair<HmLoc, complex> {
            auto c = hm.crnrs[l.id];
            auto v = c.vert();
            auto d = dir;
            for (Half h: v.adjHalfs(c.half().next().twin())) {
                d *= std::polar(1., PI / 2 * get_m(h));
                auto c1  = h.next().crnr();
                auto uv0 = c1.uv();
                auto uv1 = c1.half().cr_t().uv();
                auto uv2 = c1.half().cr_h().uv();
                if (is_points_into(uv0, uv1, uv2, uv0 + d, 0) && c1 != c) return {HmLocOnC{c1.id}, d};
            }
            METRIKO_FAIL("no face around vert {} admits the ray direction", v.id);
        },
        [&](const HmLocOnH& l) -> std::pair<HmLoc, complex> {
            auto h = hm.halfs[l.id].twin();
            return {HmLocOnH{h.id, 1. - l.r}, std::polar(1., PI / 2 * get_m(h)) * dir};
        },
        [&](const auto& _) -> std::pair<HmLoc, complex> { METRIKO_FAIL("no impl"); },
    }, loc);
}
}
#endif
