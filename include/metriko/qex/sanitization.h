//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QEX_SANITIZATION_H
#define METRIKO_QEX_SANITIZATION_H
#include "common.h"

namespace metriko::qex {
    inline void sanitization(
        Hmesh& mesh,
        const VecXi& matching,
        const VecXi& singular,
        const int rosyN
    ) {
        VecXc& cfn = mesh.cfn;
        VecXc heR;
        VecXc heT;
        compute_trs_matrix(mesh, matching, rosyN, heR, heT);
        for (Vert v: mesh.verts) {
            double max = 0;
            for (Crnr cc: v.adjCrnrs()) {
                auto uv = cfn(cc.id);
                max = std::max(max, abs(uv.real()));
                max = std::max(max, abs(uv.imag()));
            }

            bool init = false;
            for (Half h: v.adjHalfs()) {
                if(!init) {
                    if (singular[v.id]) {
                        auto cc = h.next().crnr();
                        auto uv = cfn(cc.id);
                        cfn(cc.id) = complex(std::round(uv.real()), std::round(uv.imag()));
                    } else {
                        auto cc = h.next().crnr();
                        auto uv = cfn(cc.id);
                        double delta = std::pow(2, log2(max));
                        complex sign((0 < uv.real()) - (uv.real() < 0), (0 < uv.imag()) - (uv.imag() < 0));

                        uv = (uv + delta * sign) - delta * sign;
                        // every coordinate within tol of an integer is that integer: a vertex near a grid point
                        // becomes the point, otherwise the same point is found as an edge or face q-vertex in the
                        // neighbouring faces and never pairs. a vertex near a single isoline lands exactly on it,
                        // which the tracer passes by stepping around the vertex (pick_next_half)
                        constexpr double tol = 1e-4;
                        complex g = nearby_grid(uv);
                        if (std::abs(uv.real() - g.real()) < tol) uv.real(g.real());
                        if (std::abs(uv.imag() - g.imag()) < tol) uv.imag(g.imag());
                        cfn(cc.id) = uv;
                    }
                } else {
                    Crnr cp = h.twin().prev().crnr();
                    Crnr cc = h.next().crnr();
                    complex r = heR(h.twin().id);
                    complex t = heT(h.twin().id);
                    cfn(cc.id) = r * cfn(cp.id) + t;
                }
                init = true;
            }
        }
    }
}
#endif
