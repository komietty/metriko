#ifndef GEN_Q_VERT_H
#define GEN_Q_VERT_H
#include "common.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::qex {
    inline void generate_q_vert(
        const Hmesh &mesh,
        const VecXc &cfn,
        std::vector<Qvert> &vqvs,
        std::vector<Qvert> &eqvs,
        std::vector<Qvert> &fqvs
    ) {
        constexpr double prc = 1e-9;
        constexpr double prc2 = 1e-5;

        // vert_q_vert
        int vqv_counter = 0;
        for (Vert v: mesh.verts) {
            Half h = v.half(); // sanitization required
            Crnr c = h.next().crnr();
            complex uv = cfn(c.id);
            double x = std::fmod(std::abs(uv.real()), 1.);
            double y = std::fmod(std::abs(uv.imag()), 1.);
            if ((x < prc || 1 - x < prc) && (y < prc || 1 - y < prc)) {
                vqvs.emplace_back(x, y, v.pos(), vqv_counter, v.id);
                vqv_counter++;
            }
        }

        // edge_q_vert
        int eqv_counter = 0;
        for (Edge e: mesh.edges) {
            Vert v1 = e.half().tail();
            Vert v2 = e.half().head();
            complex uv1 = cfn(e.half().next().crnr().id);
            complex uv2 = cfn(e.half().prev().crnr().id);
            int minX = std::floor(std::min(uv1.real(), uv2.real()));
            int minY = std::floor(std::min(uv1.imag(), uv2.imag()));
            int maxX = std::ceil(std::max(uv1.real(), uv2.real()));
            int maxY = std::ceil(std::max(uv1.imag(), uv2.imag()));
            for (int x = minX + 1; x < maxX; x++) {
                for (int y = minY + 1; y < maxY; y++) {
                    complex vec1 = uv2 - uv1;
                    complex vec2 = complex(x, y) - uv1;
                    complex vec3 = complex(x, y) - uv2;
                    double a = abs(vec2) / abs(vec1);
                    if (abs(vec2) < prc) continue;
                    if (abs(vec3) < prc) continue;
                    complex uv = uv1 + (uv2 - uv1) * a;
                    Row3d pos = v1.pos() + (v2.pos() - v1.pos()) * a;
                    complex n1 = normalize(vec1);
                    complex n2 = normalize(vec2);
                    if (abs(1. - dot(n1, n2)) < prc2) {
                        if (rg::any_of(vqvs, [&](const Qvert &qv) { return (pos - qv.pos).norm() < prc2; })) continue;
                        eqvs.emplace_back(uv, pos, eqv_counter, e.id);
                        eqv_counter++;
                    }
                }
            }
        }

        // face_q_vert
        int fqv_counter = 0;
        for (Face f: mesh.faces) {
            Vert v1 = mesh.verts[mesh.idx(f.id, 0)];
            Vert v2 = mesh.verts[mesh.idx(f.id, 1)];
            Vert v3 = mesh.verts[mesh.idx(f.id, 2)];
            complex uv1 = cfn(f.id * 3 + 0);
            complex uv2 = cfn(f.id * 3 + 1);
            complex uv3 = cfn(f.id * 3 + 2);
            int minX = std::floor(std::min({uv1.real(), uv2.real(), uv3.real()}));
            int minY = std::floor(std::min({uv1.imag(), uv2.imag(), uv3.imag()}));
            int maxX = std::ceil(std::max({uv1.real(), uv2.real(), uv3.real()}));
            int maxY = std::ceil(std::max({uv1.imag(), uv2.imag(), uv3.imag()}));
            for (int x = minX; x < maxX; x++) {
            for (int y = minY; y < maxY; y++) {
                auto xy = complex(x, y);
                double s1 = cross(uv2 - uv1, xy - uv1);
                double s2 = cross(uv3 - uv2, xy - uv2);
                double s3 = cross(uv1 - uv3, xy - uv3);
                Row3d pos = conversion_2d_3d(uv1, uv2, uv3, v1.pos(), v2.pos(), v3.pos(), xy);
                if ((s1 > prc2 && s2 > prc2 && s3 > prc2) || (s1 < prc2 && s2 < prc2 && s3 < prc2)) {
                    if (
                        rg::any_of(vqvs, [&](const Qvert &qv) { return (pos - qv.pos).squaredNorm() < prc2; }) ||
                        rg::any_of(eqvs, [&](const Qvert &qv) { return (pos - qv.pos).squaredNorm() < prc2; }))
                        continue;

                    fqvs.emplace_back(x, y, pos, fqv_counter, f.id);
                    fqv_counter++;
                }
            }
            }
        }
    }
}

#endif
