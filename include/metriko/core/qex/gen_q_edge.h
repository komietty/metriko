#ifndef METRIKO_QEX_GEN_Q_EDGE_H
#define METRIKO_QEX_GEN_Q_EDGE_H
#include "common.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::qex {
    inline bool predict_extrinsic_collinear(
        const Hmesh &hm,  //
        const VecXc &cf,  // corner function
        const complex o,  // origin
        const complex d,  // direction
        const Face f,     // face
        const Qport &pair //
    ) {
        Row3d pb1 = conversion_2d_3d(f, cf, o);
        Row3d pb2 = conversion_2d_3d(f, cf, o + d);
        Row3d pa1 = conversion_2d_3d(hm.faces[pair.fid], cf, pair.uv);
        Row3d pa2 = conversion_2d_3d(hm.faces[pair.fid], cf, pair.uv + pair.dir);
        Row3d da = (pa2 - pa1).normalized();
        Row3d db = (pb2 - pb1).normalized();
        Row3d dc = (pa1 - pb1).normalized();
        return abs(1. - db.dot(dc)) < EPS && abs(1. + da.dot(db)) < EPS;
    }

    inline vec<std::pair<Half, complex>> pick_next_half(
        const VecXc &cf, // corner function
        const complex o, // origin
        const complex d, // direction
        const Face f     // face
    ) {
        //std::println("pick_next_half, fid: {}", f.id);
        for (Half h: f.adjHalfs()) {
            auto uv1 = cf(h.next().crnr().id);
            auto uv2 = cf(h.prev().crnr().id);
            double rab, rcd;

            //if (abs(uv1 - o) > EPS && 1 - dot(normalize(uv1 - o), d) < EPS) {
            //    return vec{
            //        std::make_pair(h, uv1),
            //        std::make_pair(h.prev(), uv1),
            //    };
            //}

            if (find_strict_intersection(o, o + d * 1e2, uv1, uv2, rab, rcd) && rab > EPS)
                return vec{std::make_pair(h, lerp(uv1, uv2, rcd))};
        }

        return {};
        //throw std::runtime_error("no next half found");
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

        struct Cache {
            complex ori;
            complex dir;
            complex gri;
            int fid;
        };



        for (Qport &pfr: qports) {
            if (pfr.isConnected) continue;
            vec<Cache> caches = {{.ori=pfr.uv, .dir=pfr.dir, .gri=nearby_grid(pfr.uv, pfr.dir), .fid=pfr.fid}};
            std::set<int> pushed;
            //auto ori = pfr.uv;
            //auto dir = pfr.dir;
            //auto gri = nearby_grid(ori, dir);
            //Face f = mesh.faces[pfr.fid];

            while (!caches.empty()) {
                auto [ori, dir, gri, fid] = caches.back();
                caches.pop_back();
                Face f = mesh.faces[fid];

                // face-qport case
                if (is_inside_face(f, cfn, gri)) {
                    auto it = rg::find_if(fqports, [&](const Qport &p) {
                        if (pfr.isConnected || p.idx == pfr.idx || p.fid != f.id) return false;
                        return equal(p.dir, -dir) && abs(p.uv - gri) < EPS;
                    });

                    if (it != fqports.end()) {
                        pfr.isConnected = true;
                        it->isConnected = true;
                        qedges.emplace_back(pfr, *it);
                        goto loop_end;
                    }
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
                //auto [nh, hit] = pick_next_half(cfn, ori, dir, f);
                //if (nh.twin().isBoundary()) throw std::runtime_error("not implemented yet");
                //f = nh.twin().face();
                //ori = hit;
                //complex t = heT(nh.id);
                //complex r = heR(nh.id);
                //ori = r * ori + t;
                //dir = r * dir;
                //gri = nearby_grid(ori, dir);
                for (auto& [nh, hit]: pick_next_half(cfn, ori, dir, f)) {
                    if (nh.twin().isBoundary()) throw std::runtime_error("not implemented yet");
                    if (!pushed.insert(nh.id).second) continue;   // already explored this crossing
                    complex r = heR(nh.id);
                    complex t = heT(nh.id);
                    complex o2 = r * hit + t;
                    complex d2 = r * dir;
                    caches.push_back({.ori=o2, .dir=d2, .gri=nearby_grid(o2, d2), .fid=nh.twin().face().id});
                }
            }
        loop_end:
        }
        return qedges;
    }
}

#endif
