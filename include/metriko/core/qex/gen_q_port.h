#ifndef GEN_Q_PORT_H
#define GEN_Q_PORT_H
#include "common.h"

namespace metriko::qex {
    inline bool is_points_into(complex p1, complex p2, complex p3, complex uv) {
        return orientation(p1, p2, uv) > 0 && orientation(p1, p3, uv) < 0;
    }

    inline int dir_to_int(complex dir) {
        if (equal(dir, complex(1, 0))) return 0;
        if (equal(dir, complex(0, 1))) return 1;
        if (equal(dir, complex(-1, 0))) return 2;
        if (equal(dir, complex(0, -1))) return 3;
        return -1;
    }

    inline void generate_eqvert_qport(
        const Hmesh &mesh,
        const VecXc &cfn,
        const std::vector<Qvert> &eqvs,
        std::vector<Qport> &q_ports
    ) {
        std::vector<Qport> ports_per_qv;
        std::vector<Row3d> visit;

        for (const Qvert &qv: eqvs) {
            ports_per_qv.clear();
            Edge e = mesh.edges[qv.sid];
            std::vector<Half> hs;
            hs.emplace_back(e.half());
            hs.emplace_back(e.half().twin());

            std::cout << "qv.sid: " << qv.sid << std::endl;
            for (const Half h: hs) {
                std::cout << "fid: " << h.face().id << std::endl;
                std::vector<Qport> ports_per_he;
                ports_per_he.clear();
                Crnr c1 = h.next().crnr();
                Crnr c2 = h.prev().crnr();
                auto uv1 = cfn(c1.id);
                auto uv2 = cfn(c2.id);
                auto uv3 = cfn(h.crnr().id);
                Row3d p1 = c1.vert().pos();
                Row3d p2 = c2.vert().pos();
                double alpha = (qv.pos - p1).norm() / (p2 - p1).norm();
                auto uv = lerp(uv1, uv2, alpha);
                std::cout << "ori: " << orientation(uv1, uv2, uv) << std::endl;
                for (int i = 0; i < 4; i++) {
                    complex dir = get_quater_rot(i);
                    // extrinsic... better solution??
                    assert(orientation(uv1, uv2, uv3) > 0);
                    bool f1 = orientation(uv1, uv2, uv + dir) >= 0;
                    Row3d ev = conversion_2d_3d(h.face(), cfn, uv + dir).normalized();
                    auto it = rg::find_if(visit, [&](const Row3d &v) { return (ev - v).norm() < ACCURACY; });
                    if (f1 && it == visit.end()) {
                        auto l = lerp(uv1, uv2, alpha);
                        auto x = std::round(l.real());
                        auto y = std::round(l.imag());
                        ports_per_he.emplace_back(-1, -1, e.id, h.face().id, complex(x, y), dir, qv.pos);
                        visit.emplace_back(ev);
                    }
                }

                // sort ports inside a face
                rg::sort(ports_per_he, [&](const Qport &pa, const Qport &pb) {
                    const complex dir = uv2 - uv1;
                    return dot(pa.dir, dir) > dot(pb.dir, dir);
                });

                for (auto &p: ports_per_he) { ports_per_qv.emplace_back(p); }
            }

            const int s = ports_per_qv.size();
            assert(s == 4);

            for (int i = 0; i < s; i++) { ports_per_qv[i].idx = (int) q_ports.size() + i; }
            for (int i = 0; i < s; i++) {
                ports_per_qv[i].prev_id = ports_per_qv[(i - 1 + s) % s].idx;
                ports_per_qv[i].next_id = ports_per_qv[(i + 1 + s) % s].idx;
            }
            q_ports.insert(q_ports.end(), ports_per_qv.begin(), ports_per_qv.end());
        }
    }

    inline void generate_fqvert_qport(
        const Hmesh &mesh,
        const std::vector<Qvert> &fqvs,
        std::vector<Qport> &q_ports
    ) {
        for (const Qvert &qv: fqvs) {
            const int fid = mesh.faces[qv.sid].id;
            for (int i = 0; i < 4; i++) {
                q_ports.emplace_back(q_ports.size(), -1, -1, fid, qv.uv, get_quater_rot(i), qv.pos);
            }
            for (int i = 0; i < 4; i++) {
                const int l = q_ports.size();
                q_ports[l - i - 1].next_id = q_ports[l - (i - 1 + 4) % 4 - 1].idx;
                q_ports[l - i - 1].prev_id = q_ports[l - (i + 1 + 4) % 4 - 1].idx;
            }
        }
    }

    inline void generate_vqvert_qport(
        const Hmesh &mesh,
        const VecXc &cfn,
        const std::vector<Qvert> &vqvs,
        std::vector<Qport> &q_ports
    ) {
        std::vector<Qport> ports_per_qv;

        for (const Qvert &qv: vqvs) {
            ports_per_qv.clear();
            Vert vert = mesh.verts[qv.sid];
            for (Half h: vert.adjHalfs()) {
                std::vector<Qport> ports_per_he;
                ports_per_he.clear();
                auto u = cfn(h.next().crnr().id);
                auto v = cfn(h.prev().crnr().id);
                auto w = cfn(h.crnr().id);
                const double orient = orientation(u, v, w);
                if (orient <= 0) std::cerr << "not locally injective!" << std::endl;

                int r;
                for (r = 0; r < 4; r++) {
                    if (is_points_into(u, v, w, u + get_quater_rot(r))) break;
                }
                for (int i = 0; i < 4; i++) {
                    auto d = get_quater_rot((r - i + 4) % 4);
                    bool f1 = is_points_into(u, v, w, u + d);
                    bool f2 = is_collinear(u, v, u + d);
                    bool f3 = dot(v - u, d) > 0;
                    if (f1 || (f2 && f3)) ports_per_he.emplace_back(-1, vert.id, -1, h.face().id, u, d, qv.pos);
                }

                // sort ports inside a face
                rg::sort(ports_per_he, [&](const Qport &pa, const Qport &pb) {
                    complex dir = v - u;
                    return dot(pa.dir, dir) > dot(pb.dir, dir);
                });

                for (auto &p: ports_per_he) ports_per_qv.emplace_back(p);
            }
            const int s = (int) ports_per_qv.size();
            for (int i = 0; i < s; i++) { ports_per_qv[i].idx = (int) q_ports.size() + i; }
            for (int i = 0; i < s; i++) {
                ports_per_qv[i].prev_id = ports_per_qv[(i - 1 + s) % s].idx;
                ports_per_qv[i].next_id = ports_per_qv[(i + 1 + s) % s].idx;
            }

            q_ports.insert(q_ports.end(), ports_per_qv.begin(), ports_per_qv.end());
        }
    }
}

#endif
