//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_SPARSE_SOLVER_H
#define METRIKO_SPARSE_SOLVER_H
#include <random>
#include <Eigen/Sparse>
#include "metriko/common/typedef.h"
#ifdef GC_HAVE_SUITESPARSE
#include <Eigen/CholmodSupport>
#endif

namespace metriko {
// CHOLMOD when SuiteSparse is available, eigen's built-in solvers otherwise
#ifdef GC_HAVE_SUITESPARSE
template <class M> using SparseLDLT = Eigen::CholmodSimplicialLDLT<M>;
template <class M> using SparseLLT  = Eigen::CholmodSupernodalLLT<M>;
#else
template <class M> using SparseLDLT = Eigen::SimplicialLDLT<M>;
template <class M> using SparseLLT  = Eigen::SimplicialLLT<M>;
#endif
// eigen 3.4's SPQR wrapper does not compile against SuiteSparse 7 (SuiteSparse_long is int64_t), so rank
// revealing QR stays on eigen in both configurations
using SparseRankQR = Eigen::SparseQR<SprsD, Eigen::COLAMDOrdering<int>>;

// A is hermitian; only its lower triangle is read
inline VecXc solve_hermitian(const SprsC& A, const VecXc& b) {
    SparseLDLT<SprsC> ldlt(A);
    if (ldlt.info() != Eigen::Success) {
        METRIKO_FAIL("failed to factorize a hermitian system (info {}): the mesh is likely degenerate, too coarsely tessellated, or contains sharp features / bad triangles", static_cast<int>(ldlt.info()));
    }
    return ldlt.solve(b);
}

// inverse power iteration for the smallest generalized eigenvector of L x = s M x. iterates until the rayleigh
// quotient s = x* L x / x* M x stops changing (relative change below tol), or maxIter is reached
inline VecXc solve_smallest_eig(const SprsC& L, const SprsC& M, const int maxIter = 500, const double tol = 1e-8) {
    SparseLDLT<SprsC> ldlt(L);
    if (ldlt.info() != Eigen::Success) {
        METRIKO_FAIL("failed to factorize the connection Laplacian (info {}): the mesh is likely degenerate, too coarsely tessellated, or contains sharp features / bad triangles", static_cast<int>(ldlt.info()));
    }

    // fixed-seed start vector: the result must not depend on what else consumed std::rand in the process.
    // real and imaginary parts are drawn in separate statements to keep the draw order defined
    std::mt19937 rng(0);
    std::uniform_real_distribution<double> uni(-1., 1.);
    VecXc x(L.rows());
    for (int i = 0; i < x.size(); i++) {
        double re = uni(rng);
        double im = uni(rng);
        x(i) = complex(re, im);
    }

    double s_prev = std::numeric_limits<double>::infinity();
    for (int i = 0; i < maxIter; i++) {
        x = ldlt.solve(M * x);
        x /= std::sqrt(std::abs(x.dot(M * x)));   // x* M x = 1
        const double s = std::abs(x.dot(L * x));  // rayleigh quotient
        if (std::abs(s - s_prev) <= tol * std::abs(s)) break;
        s_prev = s;
    }
    return x;
}
}
#endif
