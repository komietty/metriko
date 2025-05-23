#ifndef METRIKO_QEX_GEN_Q_EDGE_H
#define METRIKO_QEX_GEN_Q_EDGE_H
#include "common.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::qex {
    inline bool predict_extrinsic_collinear(
        const Hmesh &mesh,
        const VecXc &cfn,
        complex ori,
        complex dir,
        Face face,
        const Qport &pair
    ) {
        Row3d pb1 = conversion_2d_3d(face, cfn, ori);
        Row3d pb2 = conversion_2d_3d(face, cfn, ori + dir);
        Row3d pa1 = conversion_2d_3d(mesh.faces[pair.fid], cfn, pair.uv);
        Row3d pa2 = conversion_2d_3d(mesh.faces[pair.fid], cfn, pair.uv + pair.dir);
        Row3d da = (pa2 - pa1).normalized();
        Row3d db = (pb2 - pb1).normalized();
        Row3d dc = (pa1 - pb1).normalized();
        return abs(1. - db.dot(dc)) < EPS && abs(1. + da.dot(db)) < EPS;
    }

    inline std::pair<Half, complex> pick_next_half(
        const VecXc &cfn,
        complex ori,
        complex dir,
        Face face
    ) {
        for (Half h: face.adjHalfs()) {
            auto uv1 = cfn(h.next().crnr().id);
            auto uv2 = cfn(h.prev().crnr().id);
            double rab, rcd;
            if (find_strict_intersection(ori, ori + dir * 1e2, uv1, uv2, rab, rcd) && rab > EPS)
                return std::make_pair(h, lerp(uv1, uv2, rcd));
        }
        throw std::runtime_error("no next half found");
    }

    inline std::vector<Qedge> generate_q_edge(
        const Hmesh &mesh,
        const VecXc &cfn,
        const VecXi &matching,
        std::vector<Qport> &qports
    ) {
        std::vector<Qedge> qedges;
        VecXc heR;
        VecXc heT;
        compute_trs_matrix(mesh, cfn, matching, 4, heR, heT);

        auto eqports = vw::filter(qports, [&](const Qport &qp) { return qp.eid >= 0; });
        auto vqports = vw::filter(qports, [&](const Qport &qp) { return qp.vid >= 0; });
        auto fqports = vw::filter(qports, [&](const Qport &qp) { return qp.fid >= 0; });

        for (Qport &pfr: qports) {
            if (pfr.isConnected) continue;
            auto ori = pfr.uv;
            auto dir = pfr.dir;
            auto gri = nearby_grid(ori, dir);
            Face f = mesh.faces[pfr.fid];

            while (true) {
                // face-qport case
                if (is_inside_triangle(f, cfn, gri)) {
                    auto it = rg::find_if(fqports, [&](const Qport &p) {
                        if (pfr.isConnected || p.idx == pfr.idx || p.fid != f.id) return false;
                        return equal(p.dir, -dir) && abs(p.uv - gri) < EPS;
                    });
                    assert(it != fqports.end());
                    pfr.isConnected = true;
                    it->isConnected = true;
                    qedges.emplace_back(pfr, *it);
                    goto loop_end;
                }

                // edge-qport case
                for (Half h: f.adjHalfs()) {
                    auto it = rg::find_if(eqports, [&](const Qport &p) {
                        if (p.isConnected || p.eid == pfr.eid || p.eid != h.edge().id) return false;
                        return predict_extrinsic_collinear(mesh, cfn, ori, dir, f, p);
                    });
                    if (it != eqports.end()) {
                        pfr.isConnected = true;
                        it->isConnected = true;
                        qedges.emplace_back(pfr, *it);
                        goto loop_end;
                    }
                }

                // vert-qport case
                for (Half h: f.adjHalfs()) {
                    auto it = rg::find_if(vqports, [&](const Qport &p) {
                        if (p.isConnected || p.vid == pfr.vid || p.vid != h.tail().id) return false;
                        return predict_extrinsic_collinear(mesh, cfn, ori, dir, f, p);
                    });
                    if (it != vqports.end()) {
                        pfr.isConnected = true;
                        it->isConnected = true;
                        qedges.emplace_back(pfr, *it);
                        goto loop_end;
                    }
                }

                // cannot find the pair. move to the next face
                auto [nh, hit] = pick_next_half(cfn, ori, dir, f);
                if (nh.twin().isBoundary()) throw std::runtime_error("not implemented yet");
                f = nh.twin().face();
                ori = hit;
                complex t = heT(nh.id);
                complex r = heR(nh.id);
                ori = r * ori + t;
                dir = r * dir;
                gri = nearby_grid(ori, dir);
            }
        loop_end:

        }
        return qedges;
    }
}

#endif
