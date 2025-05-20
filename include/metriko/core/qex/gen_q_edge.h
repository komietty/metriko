#ifndef GEN_Q_EDGE_H
#define GEN_Q_EDGE_H
#include "common.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::qex {

    inline complex closest_semi_integer(complex uv, complex dir) {
        double x = uv.real();
        double y = uv.imag();
        if (dir.real() != 0) {
            if (abs(fmod(x, 1.)) < ACCURACY) x += dir.real() * ACCURACY;
            x = dir.real() > 0 ? std::ceil(x) : std::floor(x);
        }
        if (dir.imag() != 0) {
            if (abs(fmod(y, 1.)) < ACCURACY) y += dir.imag() * ACCURACY;
            y = dir.imag() > 0 ? std::ceil(y) : std::floor(y);
        }
        return {x, y};
    }

    inline std::optional<std::pair<Half, complex> > pick_next_half(
        const VecXc &cfn,
        complex ori,
        complex dir,
        Face face
    ) {
        for (Half h: face.adjHalfs()) {
            auto uv1 = cfn(h.next().crnr().id);
            auto uv2 = cfn(h.prev().crnr().id);
            double rab, rcd;
            if (find_strict_intersection(ori, ori + dir * 1e2, uv1, uv2, rab, rcd) && rab > ACCURACY)
                return std::pair<Half, complex>({h, lerp(uv1, uv2, rcd)});
        }

        return std::nullopt;
    }

    /*
    inline void try_find_vqp(const Hmesh& mesh, const VecXc& cfn, Face f, ) {
        for (Half h: f.adjHalfs()) {
            Vert v = h.tail();
            //std::cout << "vid : " << v.id << std::endl;
            for (Qport &vq: vw::filter(vqports, [&](const Qport &qp) {
                return qp.vid != port.vid && qp.vid == v.id;
            })) {
                //std::cout << "pair port idx: " << vq.idx << std::endl;
                if (vq.isConnected) continue; // vq.isconnected is not alterd in the view
                Row3d pb1 = conversion_2d_3d(f, cfn, ori);
                Row3d pb2 = conversion_2d_3d(f, cfn, ori + dir);
                Row3d pa1 = conversion_2d_3d(mesh.faces[vq.fid], cfn, vq.uv);
                Row3d pa2 = conversion_2d_3d(mesh.faces[vq.fid], cfn, vq.uv + vq.dir);
                Row3d da = (pa2 - pa1).normalized();
                Row3d db = (pb2 - pb1).normalized();
                Row3d dc = (pa1 - pb1).normalized();
                if (abs(1. - db.dot(dc)) < 1e-9 && abs(1. + da.dot(db)) < 1e-9) {
                    std::cout << "pair pid: " << vq.idx << std::endl;
                    qes.emplace_back(port, vq);
                    port.isConnected = true;
                    vq.isConnected = true;
                    break;
                }
            }
        }
    }
    */

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

        auto eqports = vw::filter(qps, [&](const Qport &qp) { return qp.eid >= 0; });
        auto vqports = vw::filter(qps, [&](const Qport &qp) { return qp.vid >= 0; });
        auto fqports = vw::filter(qps, [&](const Qport &qp) { return qp.fid >= 0; });

        for (Qport &port: qps) {
            if (port.isConnected) continue;
            //if (port.idx != 184 && port.idx != 85) continue;
            //if (port.idx != 845) continue;
            //if (port.idx != 184) continue;
            complex ori = port.uv;
            complex dir = port.dir;
            Face f = mesh.faces[port.fid];
            //std::cout << "-----pid : " << port.idx << std::endl;
            //std::cout << "-----fid : " << f.id << std::endl;
            bool flag = false;

            while (!flag) {
                // fqport case
                complex gri = closest_semi_integer(ori, dir);
                if (is_inside_triangle(f, cfn, gri)) {
                    auto pair = rg::find_if(qps, [&](const Qport &op) {
                        if (op.idx == port.idx || port.isConnected || op.fid != f.id) return false;
                        return equal(op.dir, -dir) && abs(op.uv - gri) < ACCURACY;
                    });
                    if (pair == qps.end()) std::cout << "Error: no pair found for port: " << port.idx << std::endl;
                    assert(pair != qps.end());
                    //std::cout << "pair port idx: " << pair->idx << std::endl;
                    port.isConnected = true;
                    pair->isConnected = true;
                    //std::cout << "curr port fid: " << port.fid << std::endl;
                    //std::cout << "pair port fid: " << pair->fid << std::endl;
                    qes.emplace_back(port, *pair);
                    break;
                }

                // eqport case.
                for (Half h: f.adjHalfs()) {
                    Edge e = h.edge();
                    //std::cout << "vid : " << v.id << std::endl;
                    for (Qport &eq: vw::filter(eqports, [&](const Qport &qp) {
                        return qp.eid != port.eid && qp.eid == e.id;
                    })) {
                        //std::cout << "pair port idx: " << vq.idx << std::endl;
                        //std::cout << "eq idx: " << eq.idx << std::endl;
                        if (eq.isConnected) continue; // vq.isconnected is not alterd in the view
                        Row3d pb1 = conversion_2d_3d(f, cfn, ori);
                        Row3d pb2 = conversion_2d_3d(f, cfn, ori + dir);
                        Row3d pa1 = conversion_2d_3d(mesh.faces[eq.fid], cfn, eq.uv);
                        Row3d pa2 = conversion_2d_3d(mesh.faces[eq.fid], cfn, eq.uv + eq.dir);
                        Row3d da = (pa2 - pa1).normalized();
                        Row3d db = (pb2 - pb1).normalized();
                        Row3d dc = (pa1 - pb1).normalized();
                        //std::cout << "abs 1: " << abs(1. - db.dot(dc)) << std::endl;
                        //std::cout << "abs 2: " << abs(1. + da.dot(db)) << std::endl;
                        if (abs(1. - db.dot(dc)) < 1e-9 && abs(1. + da.dot(db)) < 1e-9) {
                            //std::cout << "pair port pid: " << eq.idx << std::endl;
                            qes.emplace_back(port, eq);
                            port.isConnected = true;
                            eq.isConnected = true;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) break;
                }
                if (flag) break;

                // vqport case. extrinsic check
                for (Half h: f.adjHalfs()) {
                    Vert v = h.tail();
                    //std::cout << "vid : " << v.id << std::endl;
                    for (Qport &vq: vw::filter(vqports, [&](const Qport &qp) {
                        return qp.vid != port.vid && qp.vid == v.id;
                    })) {
                        //std::cout << "pair port idx: " << vq.idx << std::endl;
                        if (vq.isConnected) continue; // vq.isconnected is not alterd in the view
                        Row3d pb1 = conversion_2d_3d(f, cfn, ori);
                        Row3d pb2 = conversion_2d_3d(f, cfn, ori + dir);
                        Row3d pa1 = conversion_2d_3d(mesh.faces[vq.fid], cfn, vq.uv);
                        Row3d pa2 = conversion_2d_3d(mesh.faces[vq.fid], cfn, vq.uv + vq.dir);
                        Row3d da = (pa2 - pa1).normalized();
                        Row3d db = (pb2 - pb1).normalized();
                        Row3d dc = (pa1 - pb1).normalized();
                        if (abs(1. - db.dot(dc)) < 1e-9 && abs(1. + da.dot(db)) < 1e-9) {
                            //std::cout << "pair port pid: " << vq.idx << std::endl;
                            qes.emplace_back(port, vq);
                            port.isConnected = true;
                            vq.isConnected = true;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) break;
                }
                if (flag) break;


                // cannot find. move to the next face
                auto res = pick_next_half(cfn, ori, dir, f);
                if (res == std::nullopt) {
                    std::cout << "Error: no next half found for port: " << port.idx << ", eid: " << port.eid << std::endl;
                    break;
                }
                auto [nh, hit] = *res;
                if (nh.twin().isBoundary()) throw std::runtime_error("not implemented yet");
                f = nh.twin().face();
                //std::cout << "-----fid : " << f.id << std::endl;
                ori = hit;
                complex t = convert(heTranslation.row(nh.id));
                complex r = get_quater_rot(heMatching[nh.id]);
                ori = r * ori + t;
                dir = r * dir;
                gri = closest_semi_integer(ori, dir);
            }

            //int counter = 0;
            //while (!is_inside(f, cfn, gri) && counter++ < 10) {
            //    auto res = pick_next_half(cfn, ori, dir, f);
            //    if (res == std::nullopt) break;
            //    auto [nh, hit] = *res;
            //    //auto [nh, hit] = pick_next_half(cfn, ori, dir, f, hid);
            //    if (nh.twin().isBoundary()) break;
            //    f = nh.twin().face();
            //    ori = hit;
            //    std::cout << "-----fid : " << f.id << std::endl;
            //    complex t = convert(heTranslation.row(nh.id));
            //    complex r = get_quater_rot(heMatching[nh.id]);
            //    ori = r * ori + t;
            //    dir = r * dir;
            //    gri = closest_semi_integer(ori, dir);
            //    std::cout << "r: " << r << std::endl;
            //    std::cout << "t: " << t << std::endl;
            //    std::cout << "ori: " << ori << std::endl;
            //    std::cout << "dir: " << dir << std::endl;
            //    std::cout << "gri: " << gri << std::endl;
            //} ;

            //auto pair = rg::find_if(qps, [&](const Qport &op) {
            //    if (op.idx == port.idx || port.isConnected || op.fid != f.id) return false;
            //    return equal(op.dir, -dir) && abs(op.uv - gri) < ACCURACY;
            //    // without this there is a bug, which means this depends the order of ports
            //});

            //if (pair == qps.end()) {
            //    std::cout << "Error: no pair found for port: " << port.idx << std::endl;
            //    continue;
            //} else {
            //    //std::cout << "pair: " << (*pair).idx << std::endl;
            //}
            //port.isConnected = true;
            //pair->isConnected = true;
            //qes.emplace_back(port, *pair);
        }
        return qes;
    }
}

#endif
