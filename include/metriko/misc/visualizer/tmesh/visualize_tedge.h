//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_VISUALIZE_TEDGE_H
#define METRIKO_VISUALIZE_TEDGE_H
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "gen_qgp_mesh.h"

namespace metriko::visualizer {
    inline void visualize_tedge(
        const Tmesh& tm,
        const VecXc& uv,
        const VecXd* X = nullptr,
        const VecXd* R = nullptr,
        const std::vector<int> &selector = std::vector<int>(),
        const std::string& prefix = std::string("")
    ) {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> teids;
        std::vector<double> randoms;
        size_t counter = 0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr(1, 30);

        std::vector<double> vecR;
        std::vector<double> vecX;
        std::vector<double> difx;
        std::vector<double> dify;

        for (int i = 0; i < tm.nTE; i++) {
            const auto& te = tm.tedges[i];
            int random_value = distr(gen);
            if (!selector.empty() && rg::find(selector, i) == selector.end()) continue;
            for (const Msgmt &ts: te.segments()) {
                Row3d p1 = conversion_2d_3d(ts.face, uv, ts.fr.uv);
                Row3d p2 = conversion_2d_3d(ts.face, uv, ts.to.uv);
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array{counter, counter + 1});
                teids.emplace_back(i);
                difx.emplace_back(ts.to.uv.real() - ts.fr.uv.real());
                dify.emplace_back(ts.to.uv.imag() - ts.fr.uv.imag());
                if (R != nullptr) vecR.emplace_back((*R)[i]);
                if (X != nullptr) vecX.emplace_back((*X)[i]);
                randoms.emplace_back(random_value);
                counter += 2;
            }
        }

        auto c = polyscope::registerCurveNetwork(prefix + "tedges", ns, es);
        c->setColor(glm::vec4(.0, .0, .0, 1.));
        c->addEdgeScalarQuantity("teid", teids);
        if (R != nullptr) c->addEdgeScalarQuantity("R", vecR);
        if (X != nullptr) c->addEdgeScalarQuantity("X", vecX);
        c->addEdgeScalarQuantity("difx", difx);
        c->addEdgeScalarQuantity("dify", dify);
        c->addEdgeScalarQuantity("random", randoms);
        c->setEnabled(false);
        c->resetTransform();
        c->setRadius(0.0005);
        c->setMaterial("flat");
    }

    /*
    inline void visualize_next_thalfs_on_joint(
        const Tmesh& tm,
        const VecXc &cfn,
        int thid_fr
    ) {
        auto pair = tm.thalfs[thid_fr].adj_thalfs();
        if (pair.size()) return;
        int pair_thid1 = pair[0].id;
        int pair_thid2 = pair[1].id;
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2> > es;
        std::vector<double> frtos;
        size_t counter = 0;
        std::cout << "thids, " << thid_fr << ", " << pair_thid1 << ", " << pair_thid2 << std::endl;

        for (int thid: std::vector{thid_fr, pair_thid1, pair_thid2}) {
            const Thalf &th = tm.thalfs[thid];
            const Tedge &te = th.edge();
            Msgmt seg = te.seg_fr;

            while (true) {
                Row3d p1 = conversion_2d_3d(seg.face, cfn, seg.fr.uv);
                Row3d p2 = conversion_2d_3d(seg.face, cfn, seg.to.uv);
                Row3d nor = seg.face.normal();
                Row3d bnr = nor.cross(p2 - p1).normalized();
                Row3d p1o = p1 + bnr * 0.002 * (th.cannonical ? 1 : -1);
                Row3d p2o = p2 + bnr * 0.002 * (th.cannonical ? 1 : -1);
                ns.emplace_back(p1o.x(), p1o.y(), p1o.z());
                ns.emplace_back(p2o.x(), p2o.y(), p2o.z());
                es.emplace_back(std::array{counter, counter + 1});
                if      (thid == thid_fr)    frtos.emplace_back(0);
                else if (thid == pair_thid1) frtos.emplace_back(1);
                else                         frtos.emplace_back(2);
                counter += 2;
                if (seg == te.seg_to) break;
                seg = seg.curv->sgmts[seg.next_id];
            }
        }

        auto c = polyscope::registerCurveNetwork("next thalfs " + std::to_string(thid_fr), ns, es);
        c->addEdgeScalarQuantity("frto", frtos)->setEnabled(true);
        c->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.001);
        c->setMaterial("flat");
    }
    */
}

#endif
