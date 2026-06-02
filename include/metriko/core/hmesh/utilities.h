//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_HMESH_UTILITIES_H
#define METRIKO_HMESH_UTILITIES_H
#include "../common/utilities.h"
#include "../common/predicates.h"
#include "hmesh.h"

namespace metriko {
inline complex calc_coefficient(
    const Face f,
    const VecXc &cf,
    const complex uv
) {
    return calc_coefficient(
        cf(f.half().crnr().id),
        cf(f.half().next().crnr().id),
        cf(f.half().prev().crnr().id),
        uv);
}

inline Row3d conversion_2d_3d(
    const complex o2, // origin of 2d
    const complex a2, //
    const complex b2, //
    const Row3d &o3,  // origin of 3d
    const Row3d &a3,  //
    const Row3d &b3,  //
    const complex uv  // target uv value
) {
    const complex c = calc_coefficient(o2, a2, b2, uv);
    return o3 + (a3 - o3) * c.real() + (b3 - o3) * c.imag();
}

inline Row3d conversion_2d_3d(
    const Face &f,
    const VecXc &cf,
    const complex uv
) {
    const Crnr c1 = f.half().crnr();
    const Crnr c2 = f.half().next().crnr();
    const Crnr c3 = f.half().prev().crnr();
    return conversion_2d_3d(
        cf(c1.id),
        cf(c2.id),
        cf(c3.id),
        c1.vert().pos(),
        c2.vert().pos(),
        c3.vert().pos(),
        uv);
}

inline bool is_inside_face(
    const Face f,
    const VecXc &cf,
    const complex uv
) {
    auto uv1 = cf(f.id * 3 + 0);
    auto uv2 = cf(f.id * 3 + 1);
    auto uv3 = cf(f.id * 3 + 2);
    return is_inside_triangle(uv1, uv2, uv3, uv);
}

}

#endif
