//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_HMESH_UTILITIES_H
#define METRIKO_HMESH_UTILITIES_H
#include "../common/utilities.h"
#include "hmesh.h"

namespace metriko {
    inline complex calc_coefficient(
        const Face face,
        const VecXc& cfn,
        const complex uv
    ) {
        return calc_coefficient(
            cfn(face.half().crnr().id),
            cfn(face.half().next().crnr().id),
            cfn(face.half().prev().crnr().id),
            uv);
    }

    inline Row3d conversion_2d_3d(
        const complex origin2d,
        const complex p1_2d,
        const complex p2_2d,
        const Row3d& origin3d,
        const Row3d& p1_3d,
        const Row3d& p2_3d,
        const complex uv
    ) {
        const complex c = calc_coefficient(origin2d, p1_2d, p2_2d, uv);
        return origin3d + (p1_3d - origin3d) * c.real() + (p2_3d - origin3d) * c.imag();
    }

    inline Row3d conversion_2d_3d(
        const Face& face,
        const VecXc& cfn,
        const complex uv
    ) {
        const Crnr c1 = face.half().crnr();
        const Crnr c2 = face.half().next().crnr();
        const Crnr c3 = face.half().prev().crnr();
        return conversion_2d_3d(
            cfn(c1.id),
            cfn(c2.id),
            cfn(c3.id),
            c1.vert().pos(),
            c2.vert().pos(),
            c3.vert().pos(),
            uv);
    }

    inline Half get_oppsite_half(
        const VecXc& cfn,
        const complex& origin,
        const complex& dir,
        const Half fr
    ) {
        for (Half h: fr.face().adjHalfs()) {
            if (h.id == fr.id) continue;
            auto a = cfn(h.prev().crnr().id) - origin;
            auto b = cfn(h.next().crnr().id) - origin;
            //if (abs(a) < 1e-3 || abs(b) < 1e-3) { return h; }
            if (cross(a, dir) * cross(b, dir) < 0) return h;
        }

        throw std::invalid_argument(
            "Consider two cases below!\n"
            "- The input uv is not strictly within uv-space of the face\n"
            "- The direction points to the joint of two halfedges\n"
        );
    }

}

#endif
