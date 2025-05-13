#ifndef GEN_Q_VERT_H
#define GEN_Q_VERT_H
#include "common.h"

namespace metriko::qex {
    inline void generate_q_vert(
        const Hmesh& mesh,
        const MatXd& cfn,
        std::vector<Qvert>& vert_q_verts,
        std::vector<Qvert>& edge_q_verts,
        std::vector<Qvert>& face_q_verts
        ) {
        constexpr double prc = 1e-9;
        constexpr double prc2 = 1e-5;

        //--- vert_q_vert ----//
        int vqv_counter = 0;
        for (Vert v: mesh.verts) {
            Half h = v.half(); // sanitization required
            Crnr c = h.next().crnr();
            Row2d uv = cfn.row(c.id);
            double u_ = std::fmod(std::abs(uv.x()), 1.);
            double v_ = std::fmod(std::abs(uv.y()), 1.);
            if ((u_ < prc || 1 - u_ < prc) && (v_ < prc || 1 - v_ < prc)) {
                vert_q_verts.emplace_back(Qvert{
                    Row2d(u_, v_),
                    v.pos(),
                    vqv_counter,
                    v.id
                });
                vqv_counter++;
            }
        }

        //--- edge_q_vert ----//
        int eqv_counter = 0;
        for (Edge e: mesh.edges) {
            Vert v1 = e.half().tail();
            Vert v2 = e.half().head();
            Row2d uv1 = cfn.row(e.half().next().crnr().id); // sanitization required
            Row2d uv2 = cfn.row(e.half().prev().crnr().id); // sanitization required
            int minX = std::floor(std::min(uv1.x(), uv2.x()));
            int maxX = std::ceil( std::max(uv1.x(), uv2.x()));
            int minY = std::floor(std::min(uv1.y(), uv2.y()));
            int maxY = std::ceil( std::max(uv1.y(), uv2.y()));
            for (int x = minX + 1; x < maxX; x++) {
            for (int y = minY + 1; y < maxY; y++) {
                Row2d vec1 = uv2 - uv1;
                Row2d vec2 = Row2d{x, y} - uv1;
                Row2d vec3 = Row2d{x, y} - uv2;
                double a = vec2.norm() / vec1.norm();
                if (vec2.norm() < prc) continue;
                if (vec3.norm() < prc) continue;
                Row2d uvw = uv1 + (uv2 - uv1) * a;
                Row3d pos = v1.pos() + (v2.pos() - v1.pos()) * a;
                if (abs(1. - vec1.normalized().dot(vec2.normalized())) < prc2) {

                    if (std::ranges::any_of(vert_q_verts.begin(), vert_q_verts.end(), [&](const Qvert& qv) {
                        return (pos - qv.pos).norm() < prc2;
                    })) continue;

                    edge_q_verts.emplace_back(
                        Qvert{
                            uvw,
                            pos,
                            eqv_counter,
                            e.id
                        }
                        );
                    eqv_counter++;
                }
            }}
        }

        //--- face_q_vert ----//
        int fqv_counter = 0;
        for (Face f: mesh.faces) {
            Vert v1 = mesh.verts[mesh.idx(f.id, 0)];
            Vert v2 = mesh.verts[mesh.idx(f.id, 1)];
            Vert v3 = mesh.verts[mesh.idx(f.id, 2)];
            Row2d uv1 = cfn.row(f.id * 3 + 0); // sanitization required
            Row2d uv2 = cfn.row(f.id * 3 + 1); // sanitization required
            Row2d uv3 = cfn.row(f.id * 3 + 2); // sanitization required
            int minX = std::floor(std::min({uv1.x(), uv2.x(), uv3.x()}));
            int minY = std::floor(std::min({uv1.y(), uv2.y(), uv3.y()}));
            int maxX = std::ceil(std::max({uv1.x(), uv2.x(), uv3.x()}));
            int maxY = std::ceil(std::max({uv1.y(), uv2.y(), uv3.y()}));
            for (int x = minX; x < maxX; x++) {
            for (int y = minY; y < maxY; y++) {
                Row2d xy(x, y);
                double s1 = cross(uv2 - uv1, xy - uv1);
                double s2 = cross(uv3 - uv2, xy - uv2);
                double s3 = cross(uv1 - uv3, xy - uv3);
                Row3d pos = uv2pos(uv1, uv2, uv3, v1.pos(), v2.pos(), v3.pos(), xy);
                if ((s1 > prc2 && s2 > prc2 && s3 > prc2) || (s1 < prc2 &&  s2 < prc2 && s3 < prc2)) {

                    if (
                        std::ranges::any_of(vert_q_verts.begin(), vert_q_verts.end(), [&](const Qvert& qv) { return (pos - qv.pos).squaredNorm() < prc2; }) ||
                        std::ranges::any_of(edge_q_verts.begin(), edge_q_verts.end(), [&](const Qvert& qv) { return (pos - qv.pos).squaredNorm() < prc2; })
                    )
                        continue;

                    face_q_verts.emplace_back(
                        Qvert{
                            Row2d(x, y),
                            pos,
                            fqv_counter,
                            f.id
                        });
                    fqv_counter++;
                }
            }}
        }
    }
}

#endif
