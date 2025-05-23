#ifndef METRIKO_QEX_GEN_Q_PORT_H
#define METRIKO_QEX_GEN_Q_PORT_H
#include "common.h"

namespace metriko::qex {
    inline void generate_eqvert_qport(
        const Hmesh &mesh,
        const VecXc &cfn,
        const std::vector<Qvert> &eqverts,
        std::vector<Qport> &qports
    ) {
        std::vector<Qport> tmp;
        std::vector<Row3d> visit;

        for (const Qvert &qv: eqverts) {
            tmp.clear();
            Edge e = mesh.edges[qv.sid];
            for (const Half h: std::vector{e.half(), e.half().twin()}) {
                Crnr c1 = h.next().crnr();
                Crnr c2 = h.prev().crnr();
                auto uv1 = cfn(c1.id);
                auto uv2 = cfn(c2.id);
                auto uv3 = cfn(h.crnr().id);
                Row3d p1 = c1.vert().pos();
                Row3d p2 = c2.vert().pos();
                auto uv = lerp(uv1, uv2, (qv.pos - p1).norm() / (p2 - p1).norm());

                int r;
                for (r = 0; r < 4; r++) {
                    if (is_points_into(uv1, uv2, uv3, uv1 + get_quater_rot(r))) break;
                }
                for (int i = 0; i < 4; i++) {
                    auto dir = get_quater_rot((r + i + 2) % 4); // considering r = 0 and boundary case
                    bool f1 = orientation(uv1, uv2, uv + dir) >= 0;
                    auto ev = conversion_2d_3d(h.face(), cfn, uv + dir).normalized();
                    auto it = rg::find_if(visit, [&](const Row3d &v) { return (ev - v).norm() < EPS; });
                    if (f1 && it == visit.end()) {
                        tmp.emplace_back(-1, -1, e.id, h.face().id, nearby_grid(uv), dir, qv.pos);
                        visit.emplace_back(ev);
                    }
                }
            }

            assert(tmp.size() == 4);
            for (int i = 0; i < 4; i++) { tmp[i].idx = (int) qports.size() + i; }
            for (int i = 0; i < 4; i++) {
                tmp[i].prev_id = tmp[(i - 1 + 4) % 4].idx;
                tmp[i].next_id = tmp[(i + 1 + 4) % 4].idx;
            }
            qports.insert(qports.end(), tmp.begin(), tmp.end());
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
                int l = (int) qports.size();
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
        std::vector<Qport> tmp;

        for (const Qvert &qv: vqverts) {
            tmp.clear();
            Vert v = mesh.verts[qv.sid];
            for (Half h: v.adjHalfs()) {
                auto uv1 = cfn(h.next().crnr().id);
                auto uv2 = cfn(h.prev().crnr().id);
                auto uv3 = cfn(h.crnr().id);
                int r;
                for (r = 0; r < 4; r++) {
                    if (is_points_into(uv1, uv2, uv3, uv1 + get_quater_rot(r))) break;
                }
                for (int i = 0; i < 4; i++) {
                    auto d = get_quater_rot((r + i + 3) % 4); // considering r = 0
                    bool f1 = is_points_into(uv1, uv2, uv3, uv1 + d);
                    bool f2 = is_collinear(uv1, uv2, uv1 + d);
                    bool f3 = dot(uv2 - uv1, d) > 0;
                    if (f1 || (f2 && f3)) tmp.emplace_back(-1, v.id, -1, h.face().id, uv1, d, qv.pos);
                }
            }
            const int s = (int) tmp.size();
            for (int i = 0; i < s; i++) { tmp[i].idx = (int) qports.size() + i; }
            for (int i = 0; i < s; i++) {
                tmp[i].prev_id = tmp[(i - 1 + s) % s].idx;
                tmp[i].next_id = tmp[(i + 1 + s) % s].idx;
            }

            qports.insert(qports.end(), tmp.begin(), tmp.end());
        }
    }
}

#endif
