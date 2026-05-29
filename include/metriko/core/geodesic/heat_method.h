#ifndef EXAMPLE_EBD_CPP_HEAT_METHOD_H
#define EXAMPLE_EBD_CPP_HEAT_METHOD_H
#include "metriko/core/hmesh/hmesh.h"

namespace metriko {
inline VecXd compute_heatflow(
    const Hmesh& hm,   //
    const SprsD& L,    // laplacian matrix
    const SprsD& M,    // mass matrix
    const VecXd& delta //
) {
    double mean = rg::fold_left(hm.edges | vw::transform([](auto e) { return e.len(); }), 0., std::plus()) / hm.nE;

    SprsD F = M + L * pow(mean, 2) * 10;
    Eigen::SimplicialLDLT solverF(F);
    VecXd u = solverF.solve(delta);

    MatXd X = MatXd::Zero(hm.nF, 3);
    VecXd D = VecXd::Zero(hm.nV);

    for (Face f: hm.faces) {
        Row3d d = Row3d::Zero();
        Row3d n = f.normal();
        auto a2 = f.area() * 2;
        for (Half h: f.adjHalfs()) { d += n.cross(h.vec()) * u(h.crnr().vert().id) / a2; }
        if (d.norm() > EPS) X.row(f.id) = -d.normalized();
    }

    for (Vert v: hm.verts) {
        double sum = 0;
        for (Half h: v.adjHalfs()) {
            if (h.isBoundary()) continue;
            auto xj = X.row(h.face().id);
            auto e1 = h.vec();
            auto e2 = h.prev().twin().vec();
            auto c1 = h.cot();
            auto c2 = h.prev().cot();
            sum += c1 * e1.dot(xj) + c2 * e2.dot(xj);
        }
        D(v.id) = sum * 0.5;
    }

    Eigen::SimplicialLDLT solverL(L);
    VecXd phi = solverL.solve(-D);
    phi.array() -= phi.minCoeff();
    return phi;
}
}

#endif
