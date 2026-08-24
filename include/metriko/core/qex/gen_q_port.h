#ifndef METRIKO_QEX_GEN_Q_PORT_H
#define METRIKO_QEX_GEN_Q_PORT_H
#include "common.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::qex {
    inline void generate_eqvert_qport(
        const Hmesh &mesh,
        const VecXc &cfn,
        const std::vector<Qvert> &eqverts,
        std::vector<Qport> &qports
    ) {
        for (const Qvert &qv: eqverts) {
            vec<std::pair<double, Qport>> ps;
            Edge e = mesh.edges[qv.sid];
            for (const Half h: std::vector{e.half(), e.half().twin()}) {
                Crnr c1  = h.next().crnr();
                Crnr c2  = h.prev().crnr();
                auto uv1 = cfn(c1.id);
                auto uv2 = cfn(c2.id);
                Row3d p1 = c1.vert().pos();
                Row3d p2 = c2.vert().pos();
                auto uv  = lerp(uv1, uv2, (qv.pos - p1).norm() / (p2 - p1).norm());

                for (int i = 0; i < 4; i++) {
                    auto d  = get_quater_rot(i);

                    if (is_collinear(uv1, uv2, uv + d) ?
                        !h.isCanonical() :                // if collinear, always non cano face will be skipped
                        orientation(uv1, uv2, uv + d) < 0 // if not collinear dir must be inside face otherwise skipped
                    ) continue;

                    double k = std::arg(d / (uv2 - uv1)); // angle from the edge: [0, pi] in this chart
                    if (k < 0) k += TwoPI;                // exact anti-parallel may yield -pi
                    ps.emplace_back((h.isCanonical() ? 0 : PI) + k, Qport(-1, -1, e.id, h.face().id, nearby_grid(uv), d, qv.pos));
                }
            }

            assert(ps.size() == 4);
            rg::sort(ps, {}, [](const auto& p) { return p.first; });

            const int s = ps.size();
            for (int i = 0; i < s; i++) { ps[i].second.idx = qports.size() + i; }
            for (int i = 0; i < s; i++) {
                ps[i].second.prev_id = ps[(i - 1 + s) % s].second.idx;
                ps[i].second.next_id = ps[(i + 1    ) % s].second.idx;
                qports.push_back(ps[i].second);
            }
        }
    }

    inline void generate_fqvert_qport(
        const Hmesh &mesh,
        const std::vector<Qvert> &fqverts,
        std::vector<Qport> &qports
    ) {
        for (const Qvert &qv: fqverts) {
            Face f = mesh.faces[qv.sid];
            for (int i = 0; i < 4; i++) {
                qports.emplace_back(qports.size(), -1, -1, f.id, qv.uv, get_quater_rot(i), qv.pos);
            }
            for (int i = 0; i < 4; i++) {
                int l = qports.size();
                qports[l - i - 1].next_id = qports[l - (i - 1 + 4) % 4 - 1].idx;
                qports[l - i - 1].prev_id = qports[l - (i + 1 + 4) % 4 - 1].idx;
            }
        }
    }

    inline void generate_vqvert_qport(
        const Hmesh &mesh,
        const VecXc &cfn,
        const std::vector<Qvert> &vqverts,
        std::vector<Qport> &qports
    ) {
        for (const Qvert &qv: vqverts) {
            vec<std::pair<double, Qport>> ps;
            double acc = 0.;
            Vert v = mesh.verts[qv.sid];
            for (Half h: v.adjHalfs()) {
                auto uv1 = cfn(h.next().crnr().id);
                auto uv2 = cfn(h.prev().crnr().id);
                auto uv3 = cfn(h.crnr().id);
                for (int i = 0; i < 4; i++) {
                    auto d  = get_quater_rot(i);
                    bool f1 = is_points_into(uv1, uv2, uv3, uv1 + d);
                    bool f2 = is_collinear(uv1, uv2, uv1 + d);
                    bool f3 = dot(uv2 - uv1, d) > 0;
                    if (f1 || (f2 && f3)) {
                        ps.emplace_back(
                            acc + std::arg(d / (uv2 - uv1)),
                            Qport{-1, v.id, -1, h.face().id, uv1, d, qv.pos}
                        );
                    }
                }
                acc += std::arg((uv3 - uv1) / (uv2 - uv1));
            }

            rg::sort(ps, {}, [](const auto& p) { return p.first; });

            const int s = ps.size();
            for (int i = 0; i < s; i++) { ps[i].second.idx = qports.size() + i; }
            for (int i = 0; i < s; i++) {
                ps[i].second.prev_id = ps[(i - 1 + s) % s].second.idx;
                ps[i].second.next_id = ps[(i + 1    ) % s].second.idx;
                qports.push_back(ps[i].second);
            }
        }
    }
}

#endif
