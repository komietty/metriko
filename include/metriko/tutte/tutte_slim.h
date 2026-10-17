//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_TUTTE_SLIM_H
#define METRIKO_TUTTE_SLIM_H
#include <igl/doublearea.h>
#include <igl/flip_avoiding_line_search.h>
#include <igl/grad.h>
#include <igl/local_basis.h>
#include <igl/polar_svd.h>
#include "metriko/hmesh/hmesh.h"
#include "metriko/solver/sparse_solver.h"

namespace metriko {

// slim (rabinovich et al. 2017) with the symmetric dirichlet energy on the mesh cut along its seam, kept seamless:
// the copies of a vertex on the seam move together under the transitions between their charts (a rotation by a
// multiple of 90 degrees and an integer translation), so the result stays an integer-grid map while the seam is
// free to move. only the cones are fixed
inline MatXd seamless_slim(
    const Hmesh&     hm,       // the mesh with its seam
    const Hmesh&     hc,       // hm cut along the seam: the same faces, corner by corner
    const vec<bool>& seam,     // per edge of hm
    const VecXi&     singular, // per vertex of hm
    const MatXd&     uv0,      // per vertex of hc: a seamless map to start from
    const int        max_iter,
    const double     rel_tol
) {
    const int n  = hc.nV;
    const int nf = hc.nF;

    // 1: the transitions across the seam: z_b = r z_a + t from a copy a of a vertex to the copy b on the other side
    struct Link { int to; complex r; complex t; };
    vec<vec<Link>> links(n);
    vec<int>       org(n);
    for (int f = 0; f < nf; ++f) for (int j = 0; j < 3; ++j) org[hc.idx(f, j)] = hm.idx(f, j);

    auto copy_at = [&](Face f, Vert v) { for (int j = 0; j < 3; ++j) if (hm.idx(f.id, j) == v.id) return hc.idx(f.id, j); METRIKO_FAIL("vertex {} is not on face {}", v.id, f.id); };
    auto uv_at   = [&](int k) { return complex(uv0(k, 0), uv0(k, 1)); };

    double gap = 0;
    for (Edge e: hm.edges) {
        if (!seam[e.id] || e.isBoundary()) continue;
        Face f0 = e.face0();
        Face f1 = e.face1();
        Vert p  = e.vert0();
        Vert q  = e.vert1();
        auto d0 = uv_at(copy_at(f0, q)) - uv_at(copy_at(f0, p));
        auto d1 = uv_at(copy_at(f1, q)) - uv_at(copy_at(f1, p));
        auto r  = get_quater_rot(std::lround(std::arg(d1 / d0) / (PI / 2)));
        auto t0 = uv_at(copy_at(f1, p)) - r * uv_at(copy_at(f0, p));
        auto t  = complex(std::round(t0.real()), std::round(t0.imag()));
        gap = std::max(gap, std::abs(t0 - t));
        for (Vert v: {p, q}) {
            int a = copy_at(f0, v);
            int b = copy_at(f1, v);
            links[a].push_back({b, r, t});
            links[b].push_back({a, std::conj(r), -std::conj(r) * t});
        }
    }
    METRIKO_CHECK(gap < 1e-6, "the seam transitions are not integer (off by {})", gap);

    // 2: one free point per vertex of hm: its first copy, the others follow it. the cones stay where they are, and so
    //    does any vertex whose transitions do not close around it
    vec<int>     var(n, -1); // the free point a copy follows, -1 when fixed
    vec<complex> rk(n, 1);   // the copy is rk z + tk of its free point z
    vec<complex> tk(n, 0);
    vec<int>     first;      // the first copy of each free point
    vec<int>     copies;     // the number of copies of each free point
    vec<bool>    seen(n, false);
    for (int k = 0; k < n; ++k) {
        if (seen[k]) continue;
        vec<int> grp = {k};
        seen[k] = true;
        bool ok = singular[org[k]] == 0;
        for (size_t i = 0; i < grp.size(); ++i) {
            int a = grp[i];
            for (auto [b, r, t]: links[a]) {
                auto rb = r * rk[a];
                auto tb = r * tk[a] + t;
                if (!seen[b]) { seen[b] = true; rk[b] = rb; tk[b] = tb; grp.push_back(b); }
                else if (std::abs(rk[b] - rb) > 1e-9 || std::abs(tk[b] - tb) > 1e-6) ok = false;
            }
        }
        if (!ok) continue;
        for (int a: grp) var[a] = first.size();
        first.push_back(k);
        copies.push_back(grp.size());
    }
    const int m = first.size();

    // the copy k at z: rk z + tk, or where it is when fixed
    auto expand = [&](const VecXd& z) {
        MatXd uv(n, 2);
        for (int k = 0; k < n; ++k) {
            auto p = var[k] < 0 ? uv_at(k) : rk[k] * complex(z(2 * var[k]), z(2 * var[k] + 1)) + tk[k];
            uv(k, 0) = p.real();
            uv(k, 1) = p.imag();
        }
        return uv;
    };
    auto shrink = [&](const MatXd& uv) {
        VecXd z(2 * m);
        for (int g = 0; g < m; ++g) { z(2 * g) = uv(first[g], 0); z(2 * g + 1) = uv(first[g], 1); }
        return z;
    };

    // 3: per face, the gradient of the three corner values in a frame of the face, and the face area
    MatXd F1, F2, F3;
    igl::local_basis(hc.pos, hc.idx, F1, F2, F3);
    SprsD G;
    igl::grad(hc.pos, hc.idx, G);
    Eigen::SparseMatrix<double, Eigen::RowMajor> Gr = G;
    MatXd dx = MatXd::Zero(nf, 3);
    MatXd dy = MatXd::Zero(nf, 3);
    for (int f = 0; f < nf; ++f) {
    for (int d = 0; d < 3; ++d) {
        for (decltype(Gr)::InnerIterator it(Gr, f + d * nf); it; ++it) {
            int j = it.col() == hc.idx(f, 0) ? 0 : it.col() == hc.idx(f, 1) ? 1 : 2;
            dx(f, j) += F1(f, d) * it.value();
            dy(f, j) += F2(f, d) * it.value();
        }
    }}
    VecXd area;
    igl::doublearea(hc.pos, hc.idx, area);

    // 4: the system over the free points. its pattern is fixed: where each face adds its 6 x 6 block is found once
    SprsD K(2 * m, 2 * m);
    {
        vec<TripD> T;
        for (int g = 0; g < 2 * m; ++g) T.emplace_back(g, g, 0.);
        for (int f = 0; f < nf; ++f)
        for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            int ga = var[hc.idx(f, a)];
            int gb = var[hc.idx(f, b)];
            if (ga < 0 || gb < 0) continue;
            for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) T.emplace_back(2 * ga + i, 2 * gb + j, 0.);
        }
        K.setFromTriplets(T.begin(), T.end());
        K.makeCompressed();
    }
    auto slot = [&](int r, int c) {
        auto bgn = K.innerIndexPtr() + K.outerIndexPtr()[c];
        auto end = K.innerIndexPtr() + K.outerIndexPtr()[c + 1];
        return static_cast<int>(std::lower_bound(bgn, end, r) - K.innerIndexPtr());
    };
    vec<arr<int, 36>> slots(nf);
    for (int f = 0; f < nf; ++f)
    for (int a = 0; a < 6; ++a)
    for (int b = 0; b < 6; ++b) {
        int ga = var[hc.idx(f, a / 2)];
        int gb = var[hc.idx(f, b / 2)];
        slots[f][a * 6 + b] = ga < 0 || gb < 0 ? -1 : slot(2 * ga + a % 2, 2 * gb + b % 2);
    }
    vec<int> diag(2 * m);
    for (int g = 0; g < 2 * m; ++g) diag[g] = slot(g, g);

    // the jacobian of a face, rows u, v and columns along the frame
    auto jacobian = [&](const MatXd& uv, int f) {
        Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
        for (int j = 0; j < 3; ++j) J += uv.row(hc.idx(f, j)).transpose() * Eigen::RowVector2d(dx(f, j), dy(f, j));
        return J;
    };
    std::function<double(MatXd&)> energy = [&](MatXd& uv) {
        double e = 0;
        #pragma omp parallel for reduction(+:e)
        for (int f = 0; f < nf; ++f) {
            auto   J   = jacobian(uv, f);
            double det = J.determinant();
            e += det > 0 ? area(f) * J.squaredNorm() * (1 + 1 / (det * det)) : std::numeric_limits<double>::infinity();
        }
        return e;
    };

    constexpr double proximal = 1e-4; // as in libigl
    SparseLLT<SprsD>  llt;   // fast, but a nearly degenerate face can make it fail
    SparseLDLT<SprsD> ldlt;  // then this one
    VecXd  z  = shrink(uv0);
    MatXd  uv = expand(z);
    double e  = energy(uv);
    int    small = 0;

    for (int it = 0; it < max_iter; ++it) {
        // 5: per face, the closest rotation R and the weights W of the quadratic proxy (the local step), and its
        //    block of the weighted least squares || W (J - R) ||^2 over the free points of the face (the global step)
        std::fill(K.valuePtr(), K.valuePtr() + K.nonZeros(), 0.);
        VecXd rhs = VecXd::Zero(2 * m);
        #pragma omp parallel for
        for (int f = 0; f < nf; ++f) {
            Eigen::Matrix2d C = Eigen::Matrix2d::Zero();          // the part of J that does not move
            Eigen::Matrix<double, 4, 6> B = Eigen::Matrix<double, 4, 6>::Zero(); // J00, J01, J10, J11 by the free points
            for (int j = 0; j < 3; ++j) {
                int k = hc.idx(f, j);
                Eigen::RowVector2d d(dx(f, j), dy(f, j));
                if (var[k] < 0) { C += uv.row(k).transpose() * d; continue; }
                C += Eigen::Vector2d(tk[k].real(), tk[k].imag()) * d;
                Eigen::Matrix2d Rk;
                Rk << rk[k].real(), -rk[k].imag(), rk[k].imag(), rk[k].real();
                for (int c = 0; c < 2; ++c) {
                    B(0, 2 * j + c) = d(0) * Rk(0, c);
                    B(1, 2 * j + c) = d(1) * Rk(0, c);
                    B(2, 2 * j + c) = d(0) * Rk(1, c);
                    B(3, 2 * j + c) = d(1) * Rk(1, c);
                }
            }

            Eigen::Matrix2d R, T, U, V;
            Eigen::Vector2d s, w;
            igl::polar_svd(jacobian(uv, f), R, T, U, s, V);
            for (int i = 0; i < 2; ++i) w(i) = std::abs(s(i) - 1) < 1e-8 ? 2. : std::sqrt((s(i) - std::pow(s(i), -3)) / (s(i) - 1));
            Eigen::Matrix2d W = U * w.asDiagonal() * U.transpose();

            // the residual W (J - R) = Wb B z - vec(W (R - C)), row by row
            Eigen::Matrix4d Wb = Eigen::Matrix4d::Zero();
            Wb(0, 0) = W(0, 0); Wb(0, 2) = W(0, 1); Wb(1, 1) = W(0, 0); Wb(1, 3) = W(0, 1);
            Wb(2, 0) = W(1, 0); Wb(2, 2) = W(1, 1); Wb(3, 1) = W(1, 0); Wb(3, 3) = W(1, 1);
            Eigen::Matrix2d D = W * (R - C);
            Eigen::Vector4d b(D(0, 0), D(0, 1), D(1, 0), D(1, 1));
            Eigen::Matrix<double, 4, 6> A = Wb * B;
            Eigen::Matrix<double, 6, 6> Kf = area(f) * A.transpose() * A;
            Eigen::Matrix<double, 6, 1> rf = area(f) * A.transpose() * b;

            for (int a = 0; a < 6; ++a) {
                int ga = var[hc.idx(f, a / 2)];
                if (ga < 0) continue;
                #pragma omp atomic
                rhs(2 * ga + a % 2) += rf(a);
                for (int c = 0; c < 6; ++c) {
                    int p = slots[f][a * 6 + c];
                    if (p < 0) continue;
                    #pragma omp atomic
                    K.valuePtr()[p] += Kf(a, c);
                }
            }
        }
        for (int g = 0; g < m; ++g)
        for (int i = 0; i < 2; ++i) {
            K.valuePtr()[diag[2 * g + i]] += proximal * copies[g];
            rhs(2 * g + i) += proximal * copies[g] * z(2 * g + i);
        }
        if (it == 0) llt.analyzePattern(K);
        llt.factorize(K);
        VecXd zd;
        if (llt.info() == Eigen::Success) zd = llt.solve(rhs);
        else {
            ldlt.compute(K);
            METRIKO_CHECK(ldlt.info() == Eigen::Success, "the slim system is singular");
            zd = ldlt.solve(rhs);
        }
        MatXd dst = expand(zd);

        // 6: the step toward it that flips no face, and the energy there
        double e_new = igl::flip_avoiding_line_search(hc.idx, uv, dst, energy, e);
        z = shrink(uv);
        small = std::abs(e - e_new) < rel_tol * std::abs(e_new) ? small + 1 : 0; // a single short step is the flip barrier, not convergence
        e = e_new;
        if (small == 2) break;
    }
    return uv;
}
}
#endif
