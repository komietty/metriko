//
// Created by saki on 2025/11/17.
//

#ifndef METRIKO_EXAMPLE_PREDICATES_H
#define METRIKO_EXAMPLE_PREDICATES_H
#include "typedef.h"

namespace metriko {

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

inline bool is_collinear(
    const complex pa,
    const complex pb,
    const complex pc,
    const double eps = EPS
) {
    double o = orientation(pa, pb, pc);
    return abs(o) < eps;
}

inline bool is_points_into(
    const complex p1,
    const complex p2,
    const complex p3,
    const complex uv,
    const double eps = EPS
) {
    return orientation(p1, p2, uv) > eps && orientation(p1, p3, uv) < -eps;
}

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
}

#endif