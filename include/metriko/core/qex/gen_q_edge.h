#ifndef GEN_Q_EDGE_H
#define GEN_Q_EDGE_H
#include "common.h"

namespace metriko::qex {
    inline bool isInside(
        const Face f,
        const MatXd &cfn,
        const Row2d &target_uv
    ) {
        const double prc = 1e-3;
        const Crnr c1 = f.half().crnr();
        const Crnr c2 = f.half().next().crnr();
        const Crnr c3 = f.half().prev().crnr();
        Row2d coef1 = uv2coef(cfn.row(c1.id), cfn.row(c2.id), cfn.row(c3.id), target_uv);
        Row2d coef2 = uv2coef(cfn.row(c2.id), cfn.row(c3.id), cfn.row(c1.id), target_uv);
        return coef1.x() >= -prc && coef1.x() <= 1 + prc &&
               coef1.y() >= -prc && coef1.y() <= 1 + prc &&
               coef2.x() >= -prc && coef2.x() <= 1 + prc &&
               coef2.y() >= -prc && coef2.y() <= 1 + prc;
    }

    inline int pickNextHe(
        Face f,
        const int vFrId,
        const int eFrId,
        const MatXd &cfn,
        const Row2d &uvFr,
        const Row2d &uvTo
    ) {
        for (Half h: f.adjHalfs()) {
            if (h.edge().id == eFrId) continue;
            Row2d uv1 = cfn.row(h.next().crnr().id);
            Row2d uv2 = cfn.row(h.prev().crnr().id);
            Row2d d = uv2 - uv1;
            if (
                cross(uvFr - uv1, d) * cross(uvTo - uv1, d) < 0 &&
                cross(uv1 - uvFr, uvTo - uvFr) * cross(uv2 - uvFr, uvTo - uvFr) < 0
            ) return h.id;
        }
        return -1;
    }

    inline void generate_q_edge(
        const Hmesh &mesh,
        const MatXd &cfn,
        const VecXi &matching,
        std::vector<Qport> &q_ports,
        std::vector<Qedge> &q_edges
    ) {
        VecXi heMatching;
        MatXd heTranslation;
        compute_trs_matrix(mesh, cfn, matching, 4, heMatching, heTranslation); // todo: tmp rosyN = 4

        for (Qport& port: q_ports) {
            if (port.isConnected) continue;
            Row2d uv1 = port.uvw;
            Row2d uv2 = port.uvw + port.dir;
            Face f = mesh.faces[port.fid];
            int eFrId = port.eid;
            int vFrId = port.vid;
            int hFrId = -1;

            int counter = 0;
            Vec2d dir = port.dir;
            do {
                hFrId = pickNextHe(f, vFrId, eFrId, cfn, uv1, uv2);
                if (hFrId == -1) break;
                Half h = mesh.halfs[hFrId];
                if (h.twin().isBoundary()) break;
                eFrId = mesh.halfs[hFrId].edge().id;
                f = h.twin().face();
                Vec2d d = heTranslation.row(h.id).transpose();
                Mat2d r = matching2rot(heMatching[h.id]);
                uv2 = r * uv2.transpose() + d;
                uv1 = r * uv1.transpose() + d;
                dir = r * dir;
                //std::cout << "fid: " << f.id << std::endl;
                //std::cout << "uv2: " <<  << std::endl;
                // todo: sign of area condition
                counter++;
            } while (!isInside(f, cfn, uv2) && counter < 10);

            auto pair = std::ranges::find_if(q_ports, [&](const Qport &other) {
                if (other.idx == port.idx) return false;
                if (other.isConnected) return false;
                // todo without this there is a bug, which means this depends the order of ports
                return other.dir == -dir.transpose() && (other.uvw - uv2).norm() < 1e-5 && other.fid == f.id;
            });

            //assert(pair != q_ports.end());
            if (pair == q_ports.end()) {
                std::cout << "port idx: " << port.idx << std::endl;
                continue;
            }
            port.isConnected = true;
            pair->isConnected = true;
            q_edges.emplace_back(port, *pair, Eigen::Matrix2d::Identity(), Row2d::Identity());
        }
    }
}

#endif
