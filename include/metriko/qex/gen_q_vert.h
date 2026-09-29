//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QEX_GEN_Q_VERT_H
#define METRIKO_QEX_GEN_Q_VERT_H
#include "common.h"
#include "metriko/hmesh/utilities.h"

namespace metriko::qex {
    inline void generate_q_vert(
        const Hmesh &mesh,
        vec<Qvert> &vqvs,
        vec<Qvert> &eqvs,
        vec<Qvert> &fqvs
    ) {
        // vert_q_vert
        for (Vert v: mesh.verts) {
            auto uv = v.half().next().crnr().uv(); // checks only one corner
            auto x = std::fmod(std::abs(uv.real()), 1.);
            auto y = std::fmod(std::abs(uv.imag()), 1.);
            if ((x < EPS || 1 - x < EPS) && (y < EPS || 1 - y < EPS)) vqvs.emplace_back(complex(x, y), v.pos(), v.id);
        }

        // edge_q_vert
        for (Edge e: mesh.edges) {
            auto uv1 = e.half().cr_t().uv();
            auto uv2 = e.half().cr_h().uv();
            auto [minX, minY, maxX, maxY] = get_minmax_int({uv1, uv2});
            for (int x = minX; x <= maxX; x++) {
            for (int y = minY; y <= maxY; y++) {
                auto xy = complex(x, y);
                auto a  = ((xy - uv1) / (uv2 - uv1)).real(); // imag part is evaluated by is_collinear
                if (is_collinear(uv1, uv2, xy) && a > EPS && a < 1 - EPS) eqvs.emplace_back(xy, e.lerp(a), e.id);
            }}
        }

        // face_q_vert
        for (Face f: mesh.faces) {
            auto [c1, c2, c3] = f.crnrs();
            auto [minX, minY, maxX, maxY] = get_minmax_int({c1.uv(), c2.uv(), c3.uv()});
            for (int x = minX; x <= maxX; x++) {
            for (int y = minY; y <= maxY; y++) {
                auto c = complex(x, y);
                if (is_inside_triangle(c1.uv(), c2.uv(), c3.uv(), c)) fqvs.emplace_back(c, f.uv2pos(c), f.id);
            }}
        }
    }
}
#endif
