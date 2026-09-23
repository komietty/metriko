//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_SUITESPARSE_INTERFACE_H
#define METRIKO_SUITESPARSE_INTERFACE_H
#include <iostream>
#include <stdexcept>
#include <string>
#include "metriko/core/common/typedef.h"
#include <Eigen/Sparse>
#ifdef GC_HAVE_SUITESPARSE
#include "suitesparse/suitesparse_pdefinite.h"
#include "suitesparse/suitesparse_square.h"
#endif

namespace metriko {
    inline void check_solver_info(const Eigen::SparseLU<SprsC>& solver, const char *what) {
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error(
                std::string("Metriko SparseLU failed to factorize the ") + what +
                " (info=" + std::to_string(static_cast<int>(solver.info())) +
                "). The mesh is likely degenerate, too coarsely tessellated, or "
                "contains sharp features / bad triangles.");
        }
    }

    inline VecXc solveSquare(SprsC& lhs, VecXc& rhs) {
        #ifdef GC_HAVE_SUITESPARSE
        VecXc x;
        auto solver = std::make_unique<suitesparse::SquareSolver<complex>>(lhs);
        solver->solve(x, rhs);
        return x;
        #else
        Eigen::SparseLU<SprsC> solver;
        solver.compute(lhs);
        check_solver_info(solver, "system matrix");
        return solver.solve(rhs);
        #endif
    }

    inline VecXc solveSmallestEig(SprsC& L, SprsC& M, int nIter = 50) {
        #ifdef GC_HAVE_SUITESPARSE
        auto solver = std::make_unique<suitesparse::PositiveDefiniteSolver<complex>>(L);
        VecXc u = VecXc::Random(L.rows());
        VecXc x = u;
        for (size_t i = 0; i < nIter; i++) {
            solver->solve(x, M * u);
            double s = std::sqrt(std::abs(x.dot(M * x)));
            x /= s;
            u = x;
        }
        return x;
        #else
        Eigen::SparseLU<SprsC> solver;
        solver.compute(L);
        check_solver_info(solver, "connection Laplacian");
        VecXc u = VecXc::Random(L.rows());
        VecXc x = u;
        for (size_t i = 0; i < nIter; i++) {
            x = solver.solve(M * u);
            double s = std::sqrt(std::abs(x.dot(M * x)));
            x /= s;
            u = x;
        }
        return x;
        #endif
    }

    inline VecXd solveSmallestEig(const SprsD& L, const SprsD& M, int nIter = 50) {
        Eigen::SimplicialLDLT<SprsD> solver;
        solver.compute(L);
        VecXd u = VecXd::Random(L.rows());
        VecXd x = u;
        for (size_t i = 0; i < nIter; i++) {
            std::cout << "iter " << i << std::endl;
            x = solver.solve(M * u);
            //double s = std::sqrt(std::abs(x.dot(M * x)));
            //x /= s;
            //u = x;
            u = x / std::sqrt(x.dot(M * x));
        }
        return x;
    }
}
#endif
