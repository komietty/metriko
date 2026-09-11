//
// Created by saki on 2026/05/08.
//

#ifndef EXAMPLE_EBD_CPP_COMMON_H
#define EXAMPLE_EBD_CPP_COMMON_H

#include "../common/predicates.h"
#include "../common/utilities.h"
#include "../hmesh/hmesh.h"

namespace metriko::mc {
// case 1
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

// case 2
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

// case 3
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

// case 4
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
}

#endif //EXAMPLE_EBD_CPP_COMMON_H
