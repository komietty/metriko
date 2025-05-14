#ifndef GEN_Q_EDGE_H
#define GEN_Q_EDGE_H
#include "common.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::qex {
    inline bool is_inside(
        const Face f,
        const VecXc &cfn,
        const Row2d &target
    ) {
        auto p = 1e-3;
        auto c = calc_coefficient(f, cfn, convert(target));
        return c.real() > -p && c.imag() > -p && c.real() + c.imag() < 1 + p;
    }

    inline int pickNextHe(
        Face f,
        const int eFrId,
        const VecXc &cfn,
        const complex &uvFr,
        const complex &uvTo
    ) {
        for (Half h: f.adjHalfs()) {
            if (h.edge().id == eFrId) continue;
            complex uv1 = cfn(h.next().crnr().id);
            complex uv2 = cfn(h.prev().crnr().id);
            complex d = uv2 - uv1;
            if (
                cross(uvFr - uv1, d) * cross(uvTo - uv1, d) < 0 &&
                cross(uv1 - uvFr, uvTo - uvFr) * cross(uv2 - uvFr, uvTo - uvFr) < 0
            ) return h.id;
        }
        return -1;
    }

    inline std::vector<Qedge> generate_q_edge(
        const Hmesh &mesh,
        const VecXc &cfn,
        const VecXi &matching,
        std::vector<Qport> &qps
    ) {
        std::vector<Qedge> qes;
        VecXi heMatching;
        MatXd heTranslation;
        compute_trs_matrix_tmp(mesh, cfn, matching, 4, heMatching, heTranslation);

        for (Qport& port: qps) {
            if (port.isConnected) continue;
            Row2d uv1 = port.uvw;
            Row2d uv2 = port.uvw + port.dir;
            Face f = mesh.faces[port.fid];
            int eFrId = port.eid;
            int hFrId = -1;

            int counter = 0;
            Vec2d dir = port.dir;
            do {
                hFrId = pickNextHe(f, eFrId, cfn, convert(uv1), convert(uv2));
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
                counter++;
            } while (!is_inside(f, cfn, uv2) && counter < 10);

            auto pair = rg::find_if(qps, [&](const Qport &other) {
                if (other.idx == port.idx) return false;
                if (other.isConnected) return false;
                // todo without this there is a bug, which means this depends the order of ports
                return other.dir == -dir.transpose() && (other.uvw - uv2).norm() < 1e-5 && other.fid == f.id;
            });

            //assert(pair != qps.end());
            if (pair == qps.end()) {
                std::cout << "Error: no pair found for port: " << port.idx << std::endl;
                continue;
            }
            port.isConnected = true;
            pair->isConnected = true;
            qes.emplace_back(port, *pair, Mat2d::Identity(), Row2d::Identity());
        }
        return qes;
    }
}

#endif
