//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_EXAMPLE_PREDICATES_H
#define METRIKO_EXAMPLE_PREDICATES_H
#include "typedef.h"

namespace metriko {
// Todo: still using epsilon predicates... consider using strict predicates
inline double orientation(
    const complex pa,
    const complex pb,
    const complex pc
) {
    double acx = pa.real() - pc.real();
    double bcx = pb.real() - pc.real();
    double acy = pa.imag() - pc.imag();
    double bcy = pb.imag() - pc.imag();
    return acx * bcy - acy * bcx;
}

// Todo: consider using side_distance or side_ratio
inline bool is_collinear(
    const complex pa,
    const complex pb,
    const complex pc,
    const double eps = EPS
) {
    double o = orientation(pa, pb, pc);
    return abs(o) < eps;
}

// Todo: consider using side_distance or side_ratio
inline bool is_points_into(
    const complex p1,
    const complex p2,
    const complex p3,
    const complex uv,
    const double eps = EPS
) {
    return orientation(p1, p2, uv) > eps && orientation(p1, p3, uv) < -eps;
}

// Todo: consider using side_distance or side_ratio
inline bool is_inside_triangle(
    const complex pa,
    const complex pb,
    const complex pc,
    const complex uv
) {
    return orientation(pa, pb, uv) > EPS &&
           orientation(pb, pc, uv) > EPS &&
           orientation(pc, pa, uv) > EPS;
}

inline double side_distance(
    const complex a,
    const complex b,
    const complex p
) {
    return orientation(a, b, p) / std::abs(b - a);
}

inline double side_ratio(
    const complex a,
    const complex b,
    const complex p
) {
    return orientation(a, b, p) / std::norm(b - a);
}

inline bool is_inside_segment(
    const complex a,
    const complex b,
    const complex p
) {
    auto ab = b - a;
    auto ap = p - a;
    auto l  = std::abs(ab); if (l < EPS) return false;
    auto d  = (ab.real() * ap.real() + ab.imag() * ap.imag()) / (l * l);
    return std::abs(side_ratio(a, b, p)) < EPS && d > EPS && d < 1 - EPS;
}
}
#endif
