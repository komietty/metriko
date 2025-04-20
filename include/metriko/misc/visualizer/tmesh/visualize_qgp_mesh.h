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

    // Get RGB color from colormap at normalized position t ∈ [0,1]
    inline glm::vec3 get_color_from_color_map(const std::string& colormapName, float t) {
        const polyscope::render::ValueColorMap& cmap =
            polyscope::render::engine->getColorMap(colormapName);
        return cmap.getValue(t);
    }

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
    const int num
) {
        auto intersections = gen_intersection_table(tmesh, tmesh.tquads[idx], arcs1, arcs2, cfn);

        auto& table = intersections.table;
        assert(intersections.pos.size() == table.rows() * table.cols());

        std::vector<glm::vec3> P;
        std::vector<double> I;
        std::vector<double> J;
        for (int i = 0; i < table.rows(); ++i) {
            for (int j = 0; j < table.cols(); ++j) {
                auto p = intersections.pos[table(i, j)];
                P.emplace_back(p.x(), p.y(), p.z());
                I.emplace_back(i);
                J.emplace_back(j);
            }}

        //auto pc = polyscope::registerPointCloud("intersection point", P);
        //pc->setEnabled(false);
        //pc->addScalarQuantity("i", I);
        //pc->addScalarQuantity("j", J);
        //pc->resetTransform();
        //pc->setPointRadius(0.001);

        MatXd V(table.rows() * table.cols(), 3);
        MatXi F((table.rows() - 1) * (table.cols() - 1), 4);

        for (int i = 0; i < table.rows(); ++i) {
            for (int j = 0; j < table.cols(); ++j) {
                auto i0 = table(i, j);
                auto p0 = intersections.pos[i0];
                V.row(i0) = p0;
            }}

        int k = 0;
        for (int i = 0; i < table.rows() - 1; ++i) {
            for (int j = 0; j < table.cols() - 1; ++j) {
                auto i0 = table(i, j);
                auto i1 = table(i + 1, j);
                auto i2 = table(i, j + 1);
                auto i3 = table(i + 1, j + 1);
                F.row(k / 2) << i0, i2, i3, i1;
                k += 2;
            }}


        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution dist(0.1, 0.6);
        double randomValue = dist(gen);

        std::vector<glm::vec3> colors;
        for (int i = 0; i < F.rows(); ++i) {
            //colors.emplace_back(get_color_from_color_map("viridis", randomValue));
            colors.emplace_back(hsv2rgb(randomValue, 0.8, 0.9));
        }

        auto surf = polyscope::registerSurfaceMesh("patch mesh-" + std::to_string(idx), V, F);
        surf->addFaceColorQuantity("color", colors)->setEnabled(true);
        surf->setShadeStyle(polyscope::MeshShadeStyle::Smooth);
        surf->setEdgeWidth(1.);
    }

}

#endif
