//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_ITER_ROUNDING_H
#define METRIKO_ITER_ROUNDING_H
#include "metriko/core/solver/levenberg_marquardt.h"
#include "iter_rounding_common.h"
#include "injective_barrier.h"
#include "iter_rounding_init.h"
#include "iter_rounding_loop.h"

namespace metriko {
    class CholeskyWrapper {
    public:
        Eigen::SimplicialLLT<SprsD> llt;

        bool factorize(const SprsD &A) { llt.compute(A); return llt.info() == Eigen::Success; }
        bool solve(const VecXd &rhs, VecXd &x) const { x = llt.solve(rhs); return true; }

        // non-const lvalue reference to type 'SprsD' (aka 'SparseMatrix<double>') cannot bind to a temporary of type
        //std::unique_ptr<solver::PositiveDefiniteSolver<double>> llt;
        //bool factorize(SprsD &A) {
        //    llt = std::make_unique<solver::PositiveDefiniteSolver<double>>(A);
        //    return true;
        //}
        //bool solve(const VecXd &rhs, VecXd &x) const {
        //    llt->solve(x, rhs);
        //    return true;
        //}
    };

    inline bool iterative_rounding(
        const VecXi &fixedIdcs,
        const VecXd &fixedVals,
        const VecXi &singularIdcs,
        const VecXi &integerIdcs,
        const double length,
        const SprsD &Cfull,
        const SprsD &G2,
        const int N,
        const int n,
        const int nF,
        const bool seamless,
        const bool roundSeams,
        const bool localInjectivity,
        const bool verbose,
        const VecXd& intrinsicField,
        const std::function<bool(const VecXd& x)> &iter_cb,
        VecXd& fullx
    ) {
        using NI = NaiveIntegration;
        using SI = SeamlessIntegration<NaiveIntegration>;

        CholeskyWrapper llt_ni, llt_si;

        DiagonalDamping<NI> dd_ni(localInjectivity ? 0.01 : 0);
        LMSolver<CholeskyWrapper, NI, DiagonalDamping<NI> > isLM;

        InjectiveBarrier barrier(N, nF, intrinsicField);

        NI ni(
            G2,
            Cfull,
            intrinsicField,
            fixedIdcs,
            fixedVals,
            length,
            localInjectivity,
            nF,
            N,
            n,
            iter_cb,
            &barrier
        );

        //----- Naive Solution -----//
        if (!seamless && !localInjectivity) {
            fullx = ni.x0;
            return true;
        }
        isLM.init(&llt_ni, &ni, &dd_ni);
        isLM.solve(verbose);

        // todo need to assign minix!
        fullx = (ni.UExt * isLM.x).head(ni.UFull.rows());
        if (!seamless) return true;

        //----- Seamless Solution -----//
        VecXd X_  = isLM.x.head(ni.UFull.cols());
        VecXd F2_ = (ni.UExt * isLM.x).tail(2 * N * nF); // x is changed from F2, but how is this important??
        std::cout << "F2 diff: " << (F2_ - intrinsicField).norm() << std::endl;
        SI si(ni, X_, F2_, integerIdcs, singularIdcs, roundSeams);

        bool success = true;
        bool rounded = false;
        DiagonalDamping<SI> dd_si(localInjectivity ? 0.01 : 0);
        LMSolver<CholeskyWrapper, SI, DiagonalDamping<SI> > irLM;

        while (!si.leftIdcs.empty()) {
            if (!si.initFixedIndices()) continue;
            rounded = true;
            dd_si.currLambda = localInjectivity ? 0.01 : 0.;
            irLM.init(&llt_si, &si, &dd_si, 100, 1e-7, 1e-7);
            irLM.solve(verbose);
            if (!si.post_checking(irLM.x)) {
                success = false;
                break;
            }
        }

        fullx = rounded ? ni.UFull * irLM.x : ni.x0;
        return success;
    }
}

#endif
