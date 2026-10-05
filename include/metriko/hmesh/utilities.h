//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_HMESH_UTILITIES_H
#define METRIKO_HMESH_UTILITIES_H
#include <queue>
#include "../common/utilities.h"
#include "../common/predicates.h"
#include "hmesh.h"
#include "hmloc.h"

namespace metriko {
inline std::optional<Crnr> try_get_crnr(const Hmesh& hm, int vid, int fid) {
    Face f = hm.faces[fid];
    Vert v = hm.verts[vid];
    for (auto h: f.adjHalfs())
        if (h.crnr().vert() == v) return h.crnr();
    return std::nullopt;
}

inline std::optional<Half> try_get_half(const Hmesh& hm, int eid, int fid) {
    Face f = hm.faces[fid];
    Edge e = hm.edges[eid];
    for (auto h: f.adjHalfs())
        if (h.edge() == e) return h;
    return std::nullopt;
}

inline complex get_face_uv(const HmLoc& loc, int fid, const Hmesh& hm) {
    return std::visit(overloaded {
        [&](const HmLocOnP& f) -> complex { return f.uv; },
        [&](const HmLocOnV& v) -> complex { return try_get_crnr(hm, v.id, fid).value().uv(); },
        [&](const HmLocOnE& e) -> complex {
            auto h = try_get_half(hm, e.id, fid).value();
            return h.uv(h.isCanonical() ? e.r : 1 - e.r);
        },
        [&](const auto& _) -> complex { METRIKO_FAIL("no impl"); },
    }, loc);
}

// per face weight of the poisson energy. around a cone the least squares integration admits a mode whose gradient
// grows like r^{-3/4} (z^{1/4} at index -1, conj(z)^{1/4} at +1): its dirichlet energy is finite, and once the mesh
// resolves it the cone is wound by -2pi. the weight (R / d)^beta, d the distance of the face to the nearest
// singular, makes the energy of that mode infinite for beta >= 1/2 and leaves the right modes (z^{5/4}, z^{3/4})
// finite for beta < 3/2. beta = 1; R = 8 grid cells is measured (4 still leaves wrong cones on bimba, rolling_stage).
// returns the diagonal matrix with the weight of each face repeated over its k rows (the 2N field components of the
// integration). gridscale: the grid spacing as a fraction of the bounding box diagonal
inline SprsD compute_poisson_weight_matrix(
    const Hmesh& hm,
    const VecXi& singular,
    const double gridscale,
    const int k
) {
    using DI = std::pair<double, int>;
    constexpr double beta  = 1.;
    constexpr double cells = 8.;
    double R = (hm.pos.colwise().maxCoeff() - hm.pos.colwise().minCoeff()).norm() * gridscale * cells;

    // edge path distance to the nearest inner singular, up to R
    VecXd d = VecXd::Constant(hm.nV, std::numeric_limits<double>::infinity());
    std::priority_queue<DI, vec<DI>, std::greater<>> pq;
    for (Vert v: hm.verts) if (!v.isBoundary() && singular[v.id] != 0) { d(v.id) = 0; pq.emplace(0., v.id); }
    while (!pq.empty()) {
        auto [dv, v] = pq.top(); pq.pop();
        if (dv > d(v) || dv > R) continue;
        for (Half h: hm.verts[v].adjHalfs()) {
            int w = h.head().id;
            if (double nd = dv + h.len(); nd < d(w)) { d(w) = nd; pq.emplace(nd, w); }
        }
    }

    VecXd w(hm.nF);
    for (Face f: hm.faces) {
        double df = 0;
        double el = 0;
        for (Half h: f.adjHalfs()) { df += std::min(d(h.tail().id), R); el += h.len(); }
        w(f.id) = std::pow(R / std::max(df / 3, el / 9), beta);
    }
    return SprsD(w.replicate(1, k).transpose().reshaped().asDiagonal());
}
}
#endif
