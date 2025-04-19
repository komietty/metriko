//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_GEN_QGP_MESH_H
#define METRIKO_GEN_QGP_MESH_H
#include "gen_qgp_table.h"

namespace metriko::visualizer {
    struct SplitTable {
        std::vector<Row3d> pos;
        MatXi table;
    };

    inline Row3d get_intersection_point(
        const SplitArc &arc1,
        const SplitArc &arc2,
        const VecXc& cfn
    ) {
        for (auto& ea : arc1.edges) {
        for (auto& eb : arc2.edges) {
            double r_a2b, r_c2d;
            if (ea.face.id == eb.face.id && find_strict_intersection(ea.uv1, ea.uv2, eb.uv1, eb.uv2, r_a2b, r_c2d)) {
                return conversion_2d_3d(ea.face, cfn, lerp(ea.uv1, ea.uv2, r_a2b));
            }
        }}
        throw std::runtime_error("no intersection");
    }

    inline SplitTable gen_intersection_table(
        const Tmesh &tmesh,
        const Tquad &tquad,
        const std::vector<SplitArc> &arcs1,
        const std::vector<SplitArc> &arcs2,
        const VecXc& cfn
    ) {
        std::vector<Row3d> pos;
        MatXi table;
        int nx = arcs1.size() + 2;
        int ny = arcs2.size() + 2;
        table.resize(nx, ny);

        for (int i = 0; i < 4; i++) {
            int thid = tquad.find_first_thid(i);
            Thalf th = tmesh.thalfs[thid];
            complex uv = th.uv_fr();
            Row3d p = conversion_2d_3d(th.sg_fr().face, cfn, uv);
            pos.push_back(p);
            int x = 0, y = 0;
            if      (i == 0) { x = 0; y = 0; }
            else if (i == 1) { x = nx - 1; y = 0; }
            else if (i == 2) { x = nx - 1; y = ny - 1; }
            else if (i == 3) { x = 0; y = ny - 1; }
            table(x, y) = pos.size() - 1;
        }

        // suppose arcs are in collect order...
        for (int i = 0; i < arcs1.size(); i++) {
            auto& arc = arcs1[i];
            auto& ef = arc.edges.front();
            auto& eb = arc.edges.back();
            pos.push_back(conversion_2d_3d(ef.face, cfn, ef.uv1));
            table(i + 1, 0) = pos.size() - 1;
            pos.push_back(conversion_2d_3d(eb.face, cfn, eb.uv2));
            table(i + 1, ny - 1) = pos.size() - 1;
        }

        for (int i = 0; i < arcs2.size(); i++) {
            auto& arc = arcs2[i];
            auto& ef = arc.edges.front();
            auto& eb = arc.edges.back();
            pos.push_back(conversion_2d_3d(ef.face, cfn, ef.uv1));
            table(nx - 1, i + 1) = pos.size() - 1;
            pos.push_back(conversion_2d_3d(eb.face, cfn, eb.uv2));
            table(0, i + 1) = pos.size() - 1;
        }

        // intersection might be flip
        for (int i = 0; i < arcs1.size(); i++) {
        for (int j = 0; j < arcs2.size(); j++) {
            auto& arc1 = arcs1[i];
            auto& arc2 = arcs2[j];
            pos.push_back(get_intersection_point(arc1, arc2, cfn));
            table(i + 1, j + 1) = pos.size() - 1;
        }}

        return {pos, table};
    }
}

#endif
