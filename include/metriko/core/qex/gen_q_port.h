#ifndef GEN_Q_PORT_H
#define GEN_Q_PORT_H
#include "common.h"

namespace metriko::qex {

    inline void generate_eqvert_qport(
        const Hmesh &mesh,
        const VecXc &cfn,
        const std::vector<Qvert> &edge_qverts,
        std::vector<Qport> &q_ports
    ) {
        std::vector<Qport> ports_per_qv;

        for (const Qvert &qv: edge_qverts) {
            ports_per_qv.clear();
            Edge edge = mesh.edges[qv.sid];
            auto hs = std::vector{edge.half(), edge.half().twin()};
            for (Half h: hs) {
                std::vector<Qport> ports_per_he;
                ports_per_he.clear();
                Crnr c1 = h.next().crnr();
                Crnr c2 = h.prev().crnr();
                complex uv1 = cfn(c1.id);
                complex uv2 = cfn(c2.id);
                complex uv3 = cfn(h.crnr().id);
                Row3d p1 = c1.vert().pos();
                Row3d p2 = c2.vert().pos();
                double alpha = (qv.pos - p1).norm() / (p2 - p1).norm();
                for (int i = 0; i < 4; i++) {
                    complex dir = get_quater_rot(i);
                    if (cross(dir, uv2 - uv1) * cross(uv3 - uv1, uv2 - uv1) > 0) {
                        complex xy = uv1 + (uv2 - uv1) * alpha;
                        complex xy_ = complex(std::round(xy.real()), std::round(xy.imag()));
                        ports_per_he.emplace_back(-1, -1, edge.id, h.face().id, convert(xy_), convert(dir), qv.pos);
                    }
                }

                // sort ports inside a face
                rg::sort(ports_per_he, [&](const Qport &pa, const Qport &pb) {
                    const complex dir = uv2 - uv1;
                    return dot(convert(pa.dir), dir) > dot(convert(pb.dir), dir);
                });

                for (auto &p: ports_per_he) { ports_per_qv.emplace_back(p); }
            }

            const int s = ports_per_qv.size();
            assert(s == 4);

            for (int i = 0; i < s; i++) { ports_per_qv[i].idx = (int)q_ports.size() + i; }
            for (int i = 0; i < s; i++) {
                ports_per_qv[i].prev_id = ports_per_qv[(i - 1 + s) % s].idx;
                ports_per_qv[i].next_id = ports_per_qv[(i + 1 + s) % s].idx;
            }
            q_ports.insert(q_ports.end(), ports_per_qv.begin(), ports_per_qv.end());
        }
    }

    inline void generate_fqvert_qport(
        const Hmesh &mesh,
        const std::vector<Qvert> &face_qverts,
        std::vector<Qport> &q_ports
    ) {
        for (const Qvert &qv: face_qverts) {
            const int fid = mesh.faces[qv.sid].id;
            for (int i = 0; i < 4; i++) {
                q_ports.emplace_back(q_ports.size(), -1, -1, fid, convert(qv.uv), convert(get_quater_rot(i)), qv.pos);
            }
            for (int i = 0; i < 4; i++) {
                const int l = q_ports.size();
                q_ports[l - i - 1].next_id = q_ports[l - (i - 1) % 4 - 1].idx;
                q_ports[l - i - 1].prev_id = q_ports[l - (i + 1) % 4 - 1].idx;
            }
        }
    }

    inline void generate_vqvert_qport(
        const Hmesh &mesh,
        const VecXc &cfn,
        const std::vector<Qvert> &vert_qverts,
        std::vector<Qport> &q_ports
    ) {
        std::vector<Qport> ports_per_qv;

        for (const Qvert &qv: vert_qverts) {
            ports_per_qv.clear();
            Vert vert = mesh.verts[qv.sid];
            for (Half h: vert.adjHalfs()) {
                std::vector<Qport> ports_per_he;
                ports_per_he.clear();
                complex u = cfn(h.next().crnr().id);
                complex v = cfn(h.prev().crnr().id);
                complex w = cfn(h.crnr().id);
                const double orient = orientation(u, v, w);
                if (orient <= 0) std::cerr << "not locally injective!" << std::endl;

                int r;
                for (r = 0; r < 4; r++) {
                    if (points_into(get_quater_rot(r), u, v, w)) break;
                }
                for (int i = 0; i < 4; i++) {
                    complex dir = get_quater_rot((r - i + 4) % 4);
                    if (points_into(dir, u, v, w) || dot(dir, normalize(v - u)) == 1.) {
                        ports_per_he.emplace_back(-1, vert.id, -1, h.face().id, convert(u), convert(dir), qv.pos);
                    }
                }

                // sort ports inside a face
                rg::sort(ports_per_he, [&](const Qport &pa, const Qport &pb) {
                    const complex dir = v - u;
                    return dot(convert(pa.dir), dir) > dot(convert(pb.dir), dir);
                });

                for (auto & p : ports_per_he) { ports_per_qv.emplace_back(p); }
            }
            const int s = (int)ports_per_qv.size();
            for (int i = 0; i < s; i++) { ports_per_qv[i].idx = (int)q_ports.size() + i; }
            for (int i = 0; i < s; i++) {
                ports_per_qv[i].prev_id = ports_per_qv[(i - 1 + s) % s].idx;
                ports_per_qv[i].next_id = ports_per_qv[(i + 1 + s) % s].idx;
            }

            q_ports.insert(q_ports.end(), ports_per_qv.begin(), ports_per_qv.end());
        }
    }
}

#endif
