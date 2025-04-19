//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_VISUALIZE_QGP_EDGE_H
#define METRIKO_VISUALIZE_QGP_EDGE_H
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "gen_qgp_mesh.h"

namespace metriko::visualizer {
inline void visualize_split_vert(
    const Tmesh& tmesh,
    const VecXc& cfn,
    const std::vector<std::vector<SplitVert>>& split_verts
) {
    std::vector<glm::vec3> ps;
    for (int ei = 0; ei < tmesh.nTE; ei++) {
        for (auto &s: split_verts[ei]) {
            Row3d p = conversion_2d_3d(s.segment.face, cfn, s.uv);
            ps.emplace_back(p.x(), p.y(), p.z());
        }
    }
    auto p = polyscope::registerPointCloud("tedge split verts", ps);
    p->setEnabled(false);
    p->setPointRadius(0.002);
    p->resetTransform();
}

inline void visualize_split_table(
    const VecXc& cfn,
    const std::vector<std::vector<SplitElem>> &sv_quad0,
    const int tqid
) {
    std::vector<glm::vec3> ps;
    std::vector<double> sides;
    std::vector<double> order;
    std::vector<double> dists;
    for (int i = 0; i < 4; i++) {
        int c = 0;
        for (auto &sv: sv_quad0[i]) {
            Row3d p = conversion_2d_3d(sv.segment.face, cfn, sv.uv);
            ps.emplace_back(p.x(), p.y(), p.z());
            sides.emplace_back(i);
            order.emplace_back(c);
            dists.emplace_back(sv.distance);
            c++;
        }
    }
    auto p = polyscope::registerPointCloud("tquad split verts-" + std::to_string(tqid), ps);
    p->setEnabled(true);
    p->addScalarQuantity("sides", sides);
    p->addScalarQuantity("order", order);
    p->addScalarQuantity("distance", dists);
    p->setPointRadius(0.002);
    p->resetTransform();
}

inline void visualize_split_edge(
    const Tmesh& tmesh,
    const VecXd& R,
    const VecXc& cfn,
    const VecXi& matching,
    const std::vector<std::vector<SplitElem>> &sv_quad0,
    const int tqid,
    std::vector<SplitArc> &arcs1,
    std::vector<SplitArc> &arcs2
) {
    std::vector<glm::vec3> ps;
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2> > es;
    std::vector<double> be;
    size_t counter = 0;

    for (int s = 0; s < 2; s++) {
        auto &svs0 = sv_quad0[s];
        auto &svs2 = sv_quad0[(s + 2) % 4];

        assert(svs0.size() == svs2.size());
        for (int i = 1; i < svs0.size() - 1; i++) {
            auto sv_a = svs0[i];
            auto sv_b = svs2[svs2.size() - i - 1];

            Row3d pos_fr = conversion_2d_3d(sv_a.segment.face, cfn, sv_a.uv);
            Row3d pos_to = conversion_2d_3d(sv_b.segment.face, cfn, sv_b.uv);
            ps.emplace_back(pos_fr.x(), pos_fr.y(), pos_fr.z());
            ps.emplace_back(pos_to.x(), pos_to.y(), pos_to.z());
            be.emplace_back(0.);
            be.emplace_back(1.);


            const Tquad &tq_draw = tmesh.tquads[tqid];
            auto res = gen_split_arc(tmesh, tq_draw, R, cfn, matching, s, sv_a, sv_b);

            for (auto &re: res.edges) {
                Row3d p1 = conversion_2d_3d(re.face, cfn, re.uv1);
                Row3d p2 = conversion_2d_3d(re.face, cfn, re.uv2);
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array{counter, counter + 1});
                counter += 2;
            }

            if (s == 0) arcs1.emplace_back(res);
            else arcs2.emplace_back(res);
        }
    }

    auto c = polyscope::registerCurveNetwork("split edges-" + std::to_string(tqid), ns, es);
    c->setEnabled(false);
    c->resetTransform();
    c->setRadius(0.0005);
    c->setMaterial("flat");

    auto p = polyscope::registerPointCloud("tedge split edge fr to-" + std::to_string(tqid), ps);
    p->setEnabled(false);
    p->addScalarQuantity("bgn 0 end 1", be);
    p->setPointRadius(0.001);
    p->resetTransform();
}

inline void visualize_split_edge(
    const Tmesh& tmesh,
    const VecXd& X,
    const VecXd& R,
    const VecXc& cfn,
    const VecXi& matching,
    std::vector<SplitArc> &arcs1,
    std::vector<SplitArc> &arcs2
) {
    std::vector<std::vector<SplitVert>> split_verts(tmesh.nTE);
    for (int i = 0; i < tmesh.nTE; i++)
        split_verts[i] = construct_verts_on_tedge(tmesh, X, R, i);

    visualize_split_vert(tmesh, cfn, split_verts);

    //{
    //    const int tqid = 0;
    //    const Tquad &tq_draw = tmesh.tquads[tqid];
    //    std::vector<std::vector<SplitElem> > sv_quad0;
    //    gen_split_table(tmesh, tq_draw, split_verts, sv_quad0);
    //    visualize_split_table(cfn, sv_quad0, tqid);
    //    for (int s = 0; s < 2; s++) {
    //        auto &svs0 = sv_quad0[s];
    //        auto &svs2 = sv_quad0[(s + 2) % 4];
    //        std::cout << "size of "<< s << "is: " << svs0.size() << std::endl;
    //        std::cout << "size of "<< (s + 2) % 4 << "is: " << svs2.size() << std::endl;
    //    }
    //    //visualize_split_edge(tmesh, R, cfn, matching, sv_quad0, tqid);
    //}

    for (int ii = 0; ii < tmesh.tquads.size(); ii++) {
        const Tquad &tq_draw = tmesh.tquads[ii];
        std::vector<std::vector<SplitElem> > sv_quad0;
        gen_split_table(tmesh, tq_draw, split_verts, sv_quad0);
        //visualize_split_table(cfn, sv_quad0, ii);
        for (int s = 0; s < 2; s++) {
            auto &svs0 = sv_quad0[s];
            auto &svs2 = sv_quad0[s + 2];

            if (svs0.size() != svs2.size()) {
                std::cout << "tqid " << ii << ",  side " << s << " is: " << svs0.size() << ", " << svs2.size() << std::endl;

                auto thidsA = tq_draw.thids_by_side(s);
                auto thidsB = tq_draw.thids_by_side(s + 2);
                int sumA = 0;
                int sumB = 0;
                for (int thid: thidsA) {
                    int teid = tmesh.thalfs[thid].edge().id;
                    sumA += (int) X[teid];
                    std::cout << "val: " << X[teid] << ", teid: " << teid << std::endl;
                }
                std::cout << "sumA: " << sumA << std::endl;

                for (int thid: thidsB) {
                    int teid = tmesh.thalfs[thid].edge().id;
                    sumB += (int) X[teid];
                    std::cout << "val: " << X[teid] << ", teid: " << teid << std::endl;
                }
                std::cout << "sumB: " << sumB << std::endl;
            }
        }
        visualize_split_edge(tmesh, R, cfn, matching, sv_quad0, ii, arcs1, arcs2);
    }
}
}

#endif
