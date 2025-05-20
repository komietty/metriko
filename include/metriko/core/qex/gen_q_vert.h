#ifndef GEN_Q_VERT_H
#define GEN_Q_VERT_H
#include "common.h"
#include "gen_q_edge.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::qex {
    inline void generate_q_vert(
        const Hmesh &mesh,
        const VecXc &cfn,
        std::vector<Qvert> &vqvs,
        std::vector<Qvert> &eqvs,
        std::vector<Qvert> &fqvs
    ) {

        // vert_q_vert
        for (Vert v: mesh.verts) {
            auto c = v.half().next().crnr(); // checks only one corner
            auto uv = cfn(c.id);
            auto x = std::fmod(std::abs(uv.real()), 1.);
            auto y = std::fmod(std::abs(uv.imag()), 1.);
            if ((x < ACCURACY || 1 - x < ACCURACY) && (y < ACCURACY || 1 - y < ACCURACY))
                vqvs.emplace_back(complex(x, y), v.pos(), v.id);
        }

        // edge_q_vert
        for (Edge e: mesh.edges) {
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            auto uv1 = cfn(e.half().next().crnr().id);
            auto uv2 = cfn(e.half().prev().crnr().id);
            auto [minX, minY, maxX, maxY] = get_minmax_int({uv1, uv2});
            for (int x = minX; x <= maxX; x++) {
            for (int y = minY; y <= maxY; y++) {
                auto xy = complex(x, y);
                auto a = abs(xy - uv1) / abs(uv2 - uv1);
                if (is_collinear(uv1, uv2, xy) && a > 0 && a < 1) {
                    eqvs.emplace_back(xy, p1 + (p2 - p1) * a, e.id);
                }
            }}
        }

        // face_q_vert
        auto nested = std::array{vqvs, eqvs} | vw::join;
        for (Face f: mesh.faces) {
            Row3d p1 = mesh.verts[mesh.idx(f.id, 0)].pos();
            Row3d p2 = mesh.verts[mesh.idx(f.id, 1)].pos();
            Row3d p3 = mesh.verts[mesh.idx(f.id, 2)].pos();
            auto uv1 = cfn(f.id * 3 + 0);
            auto uv2 = cfn(f.id * 3 + 1);
            auto uv3 = cfn(f.id * 3 + 2);
            auto [minX, minY, maxX, maxY] = get_minmax_int({uv1, uv2, uv3});
            for (int x = minX; x <= maxX; x++) {
            for (int y = minY; y <= maxY; y++) {
                auto xy = complex(x, y);
                auto p = conversion_2d_3d(uv1, uv2, uv3, p1, p2, p3, xy);
                if (is_inside_triangle(uv1, uv2, uv3, xy)) fqvs.emplace_back(xy, p, f.id);
            }}
        }
    }
}

#endif
