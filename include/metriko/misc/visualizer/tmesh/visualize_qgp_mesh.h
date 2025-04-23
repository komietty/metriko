//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_VISUALIZE_QGP_MESH_H
#define METRIKO_VISUALIZE_QGP_MESH_H
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include <polyscope/surface_mesh.h>

#include "gen_qgp_mesh.h"

namespace metriko::visualizer {

    inline glm::vec3 hsv2rgb(float h, float s, float v) {
        float r, g, b;
        int i = int(h * 6);
        float f = h * 6 - i;
        float p = v * (1 - s);
        float q = v * (1 - f * s);
        float t = v * (1 - (1 - f) * s);
        switch(i % 6) {
            case 0: r = v, g = t, b = p; break;
            case 1: r = q, g = v, b = p; break;
            case 2: r = p, g = v, b = t; break;
            case 3: r = p, g = q, b = v; break;
            case 4: r = t, g = p, b = v; break;
            case 5: r = v, g = p, b = q; break;
            default: r = 0, g = 0, b = 0; break;
        }
        return {r, g, b};
    }

    inline void visualize_split_mesh(
    const std::vector<SplitArc> &arcs1,
    const std::vector<SplitArc> &arcs2,
    const Tmesh& tmesh,
    const VecXc& cfn,
    const int idx,
    MatXd& VQ,
    MatXi& FQ,
    std::vector<glm::vec3>& CQ
) {
        int nVQ = VQ.rows();
        int nFQ = FQ.rows();
        auto [pos, table] = gen_intersection_table(tmesh, tmesh.tquads[idx], arcs1, arcs2, cfn);
        assert(pos.size() == table.rows() * table.cols());

        MatXd V(table.rows() * table.cols(), 3);
        MatXi F((table.rows() - 1) * (table.cols() - 1), 4);

        for (int i = 0; i < table.rows(); ++i) {
        for (int j = 0; j < table.cols(); ++j) {
            auto i0 = table(i, j);
            auto p0 = pos[i0];
            V.row(i0) = p0;
        }}

        int k = 0;
        for (int i = 0; i < table.rows() - 1; ++i) {
        for (int j = 0; j < table.cols() - 1; ++j) {
            int i0 = table(i, j);
            int i1 = table(i + 1, j);
            int i2 = table(i, j + 1);
            int i3 = table(i + 1, j + 1);
            F.row(k) << nVQ + i0, nVQ + i2, nVQ + i3, nVQ + i1;
            k++;
        }}

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution dist(0.15, 0.55);
        double randomValue = dist(gen);

        for (int i = 0; i < F.rows(); ++i) {
            CQ.emplace_back(hsv2rgb(randomValue, 0.9, 1.0));
        }

        VQ.conservativeResize(nVQ + V.rows(), 3);
        FQ.conservativeResize(nFQ + F.rows(), 4);
        VQ.block(nVQ, 0, V.rows(), 3) << V;
        FQ.block(nFQ, 0, F.rows(), 4) << F;


        //auto c = polyscope::registerCurveNetwork("boundary of qgp mesh" + std::to_string(idx), ns, es);
        //c->setEnabled(false);
        //c->resetTransform();
        //c->setRadius(0.0005);
        //c->setMaterial("flat");
    }
}

#endif
