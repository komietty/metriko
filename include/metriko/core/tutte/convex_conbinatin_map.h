//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_CONVEX_CONBINATIN_MAP_H
#define METRIKO_CONVEX_CONBINATIN_MAP_H

namespace metriko {

    inline SprsD cotan_laplacian(const Hmesh &mesh) {
        SprsD S(mesh.nV, mesh.nV);
        std::vector<TripD> T;

        for (Vert v: mesh.verts) {
            double sum = 0.;
            for (Half h: v.adjHalfs()) {
                const double c = h.edge().cot();
                sum += c;
                T.emplace_back(v.id, h.head().id, -c);
            }
            T.emplace_back(v.id, v.id, sum);
        }
        S.setFromTriplets(T.begin(), T.end());
        return S;
    }

    inline SprsD mass_matrix(const Hmesh &mesh) {
        SprsD M(mesh.nV, mesh.nV);
        std::vector<TripD> T;

        for (Vert v: mesh.verts) {
            T.emplace_back(v.id, v.id, v.baryArea());
        }
        M.setFromTriplets(T.begin(), T.end());
        return M;
    }

    inline SprsD boundary_snap_laplacian(const Hmesh &mesh) {
        SprsD S(mesh.nV, mesh.nV);
        std::vector<TripD> T;

        for (Vert v: mesh.verts) {
            if (v.isBoundary()) T.emplace_back(v.id, v.id, 1);
            else {
                double sum = 0.;
                for (Half h: v.adjHalfs()) {
                    const double c = 1.;
                    sum += c;
                    T.emplace_back(v.id, h.head().id, c);
                }
                T.emplace_back(v.id, v.id, -sum);
            }
        }
        S.setFromTriplets(T.begin(), T.end());
        return S;
    }
}

#endif
