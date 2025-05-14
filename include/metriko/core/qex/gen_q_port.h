#ifndef GEN_Q_PORT_H
#define GEN_Q_PORT_H
#include "common.h"

namespace metriko::qex {
    template<typename T>
    double orient2d(const T &u, const T &v, const T &w) { return cross(v - u, w - u); }

    template<typename T>
    bool points_into(const T &d, const T &u, const T &v, const T &w) {
        return isccw(u, v, u + d) && isccw(u, u + d, w);
    }

    inline Row2d rot(int i) {
        switch (i) {
            case 0: return Row2d(1, 0);
            case 1: return Row2d(0, 1);
            case 2: return Row2d(-1, 0);
            case 3: return Row2d(0, -1);
            default: throw std::runtime_error("invalid rot idx");
        }
    }

    inline double dot(
        const Row2d &a,
        const Row2d &b
    ) {
        return a.x() * b.x() + a.y() * b.y();
    }

    inline void generate_eqvert_qport(
        const Hmesh &mesh,
        const MatXd &cfn,
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
                Row2d uv1 = cfn.row(c1.id);
                Row2d uv2 = cfn.row(c2.id);
                Row2d uv3 = cfn.row(h.crnr().id);
                Row3d p1 = c1.vert().pos();
                Row3d p2 = c2.vert().pos();
                double alpha = (qv.pos - p1).norm() / (p2 - p1).norm();
                for (int i = 0; i < 4; i++) {
                    Row2d dir = rot(i);
                    if (cross(dir, uv2 - uv1) * cross(uv3 - uv1, uv2 - uv1) > 0) {
                        Row2d xy = uv1 + (uv2 - uv1) * alpha;
                        Row2d xy_ = Row2d(std::round(xy.x()), std::round(xy.y()));
                        ports_per_he.emplace_back(-1, qv.idx, -1, edge.id, h.face().id, xy_, dir, qv.pos);
                    }
                }

                // sort ports inside a face
                rg::sort(ports_per_he, [&](const Qport &pa, const Qport &pb) {
                    const Row2d dir = uv2 - uv1;
                    return dot(pa.dir, dir) > dot(pb.dir, dir);
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
        const MatXd &cfn,
        const std::vector<Qvert> &face_qverts,
        std::vector<Qport> &q_ports
    ) {
        for (const Qvert &qv: face_qverts) {
            Face face = mesh.faces[qv.sid];
            for (int i = 0; i < 4; i++) {
                q_ports.emplace_back(q_ports.size(), qv.idx, -1, -1, face.id, qv.uvw, rot(i), qv.pos);
            }
            const int l = q_ports.size();
            q_ports[l - 4].next_id = q_ports[l - 3].idx;
            q_ports[l - 4].prev_id = q_ports[l - 1].idx;
            q_ports[l - 3].next_id = q_ports[l - 2].idx;
            q_ports[l - 3].prev_id = q_ports[l - 4].idx;
            q_ports[l - 2].next_id = q_ports[l - 1].idx;
            q_ports[l - 2].prev_id = q_ports[l - 3].idx;
            q_ports[l - 1].next_id = q_ports[l - 4].idx;
            q_ports[l - 1].prev_id = q_ports[l - 2].idx;
        }
    }

    inline void generate_vqvert_qport(
        const Hmesh &mesh,
        const MatXd &cfn,
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
                Face f = h.face();
                Row2d u = cfn.row(h.next().crnr().id);
                Row2d v = cfn.row(h.prev().crnr().id);
                Row2d w = cfn.row(h.crnr().id);
                const double orient = orient2d(u, v, w);
                if (orient <= 0) std::cerr << "orientation is not locally injective!" << std::endl;

                int r;
                for (r = 0; r < 4; r++) { if (points_into(rot(r), u, v, w)) break; }
                for (int i = 0; i < 4; i++) {
                    Row2d dir = rot((r - i + 4) % 4);
                    if (points_into(dir, u, v, w) || dir.dot((v - u).normalized()) == 1.) {
                        ports_per_he.emplace_back(-1, vert.id, qv.idx, -1, f.id, u, dir, qv.pos);
                    }
                }

                // sort ports inside a face
                rg::sort(ports_per_he, [&](const Qport &pa, const Qport &pb) {
                    const Row2d dir = v - u;
                    return dot(pa.dir, dir) > dot(pb.dir, dir);
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
