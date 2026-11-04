//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_SOLVE_QP_H
#define METRIKO_SOLVE_QP_H
#include "sparse_solver.h"

namespace metriko {
// min 1/2 x^T H x + g^T x  s.t.  A x >= b, H positive definite, by the predictor-corrector interior point method of
// nocedal and wright, numerical optimization (2nd ed.), algorithm 16.4: slacks y = A x - b >= 0 and multipliers l >= 0,
// every newton step solved on the normal equations (H + A^T Y^-1 L A) dx = ... (16.62), the starting point by the
// heuristic of p. 484. x is the starting point and the result. returns whether the primal residual and the
// complementarity reached tol (the dual residual too, unless the factorization fails or the iterations run out
// first), false e.g. on an infeasible problem
inline bool solve_convex_qp(const SprsD& H, const VecXd& g, const SprsD& A, const VecXd& b, VecXd& x, const int maxIter = 200, const double tol = 1e-9) {
    const int m = A.rows();
    const SprsD At = A.transpose();
    SparseLDLT<SprsD> ldlt;
    bool analyzed = false;
    VecXd y, l;

    auto factorize = [&] {
        const SprsD M = H + At * l.cwiseQuotient(y).asDiagonal() * A;
        if (!analyzed) { ldlt.analyzePattern(M); analyzed = true; }
        ldlt.factorize(M);
        return ldlt.info() == Eigen::Success;
    };
    // the newton step (16.58) whose complementarity rows L dy + Y dl read r3
    auto step = [&](const VecXd& rd, const VecXd& rp, const VecXd& r3, VecXd& dx, VecXd& dy, VecXd& dl) {
        dx = ldlt.solve(VecXd(-rd + At * (r3 - l.cwiseProduct(rp)).cwiseQuotient(y)));
        dy = A * dx + rp;
        dl = (r3 - l.cwiseProduct(dy)).cwiseQuotient(y);
    };
    // the largest a in (0, 1] keeping v + a dv >= 0
    auto max_step = [](const VecXd& v, const VecXd& dv) {
        double a = 1;
        for (int i = 0; i < v.size(); ++i) if (dv(i) < 0) a = std::min(a, -v(i) / dv(i));
        return a;
    };

    // starting point: one affine step from (x, |A x - b|, 1), then the slacks and multipliers lifted to at least 1
    VecXd dx, dy, dl, dxa, dya, dla;
    y = (A * x - b).cwiseAbs().cwiseMax(1.);
    l = VecXd::Ones(m);
    if (!factorize()) return false;
    step(H * x + g - At * l, A * x - y - b, -l.cwiseProduct(y), dxa, dya, dla);
    y = (y + dya).cwiseAbs().cwiseMax(1.);
    l = (l + dla).cwiseAbs().cwiseMax(1.);

    double gn = 1 + g.lpNorm<Eigen::Infinity>(), bn = 1 + b.lpNorm<Eigen::Infinity>();
    for (int it = 0; it < maxIter; ++it) {
        VecXd rd = H * x + g - At * l;
        VecXd rp = A * x - y - b;
        double mu = y.dot(l) / m;
        // acceptable: primal feasible and complementary. the dual residual may stall above tol once the normal
        // equations lose their definiteness near the solution, where the slacks of the active rows vanish
        bool acceptable = rp.lpNorm<Eigen::Infinity>() <= tol * bn && mu <= tol;
        if (acceptable && rd.lpNorm<Eigen::Infinity>() <= tol * gn) return true;
        if (!factorize()) return acceptable;

        // predictor: the affine scaling step, giving the centering sigma = (mu_aff / mu)^3
        step(rd, rp, -l.cwiseProduct(y), dxa, dya, dla);
        double aa = std::min(max_step(y, dya), max_step(l, dla));
        double sigma = std::pow((y + aa * dya).dot(l + aa * dla) / m / mu, 3);

        // corrector, with the second order term of the affine step
        step(rd, rp, VecXd(-l.cwiseProduct(y) - dla.cwiseProduct(dya) + VecXd::Constant(m, sigma * mu)), dx, dy, dl);
        double a = std::min(1., 0.995 * std::min(max_step(y, dy), max_step(l, dl)));
        x += a * dx;
        y += a * dy;
        l += a * dl;
    }
    return (A * x - y - b).lpNorm<Eigen::Infinity>() <= tol * bn && y.dot(l) / m <= tol;
}
}
#endif
