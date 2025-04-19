//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_GEN_QGP_TABLE_H
#define METRIKO_GEN_QGP_TABLE_H
#include "gen_qgp_edge.h"

namespace metriko::visualizer {
    inline double compute_accum_value(
        const SplitVert& sv0,
        const SplitVert& sv1
    ) {
        double result = 0;
        assert(sv0.segment.curv == sv1.segment.curv);
        if (sv0.segment.id == sv1.segment.id)  return abs(sv1.uv - sv0.uv);

        bool c = sv0.segment.id < sv1.segment.id;
        const SplitVert& fr = c ? sv0 : sv1;
        const SplitVert& to = c ? sv1 : sv0;
        Mcurv* spline = fr.segment.curv;

        bool bgn = false;
        for (const Msgmt &sg: spline->sgmts) {
            if (bgn) {
                if (sg == to.segment) {
                    result += abs(to.uv - sg.fr.uv);
                    break;
                }
                result += abs(sg.to.uv - sg.fr.uv);
            }
            if (sg == fr.segment) {
                result += abs(sg.to.uv - fr.uv);
                bgn = true;
            }
        }
        return result;
    }


    inline void gen_split_table(
        const Tmesh& tmesh,
        const Tquad& tquad,
        const std::vector<std::vector<SplitVert>>& sverts_per_tedge,
        std::vector<std::vector<SplitElem>>& sverts_per_tquad
    ) {
        for (int i = 0; i < 4; i++) {
            std::vector<SplitElem> sverts_per_side;
            auto& th = tmesh.thalfs[tquad.find_first_thid(i)];
            auto uv0 = th.uv_fr();
            auto sg0 = th.sg_fr();
            auto sv0 = SplitVert(sg0, uv0);
            for (int thid: tquad.thids_by_side(i)) {
                auto &th = tmesh.thalfs[thid];
                auto &vs = sverts_per_tedge[th.edge().id];
                for (int j = 0; j < vs.size(); j++) {
                    if (!sverts_per_side.empty() && j == 0) continue;
                    auto v = vs[th.cannonical ? j : vs.size() - 1 - j];
                    auto d = compute_accum_value(sv0, v);
                    sverts_per_side.emplace_back(v.segment, v.uv, d);
                }
            }

            sverts_per_tquad.emplace_back(sverts_per_side);
        }
    }
}

#endif
