//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_ITER_ROUNDING_H
#define METRIKO_ITER_ROUNDING_H
#include <algorithm>
#include "metriko/solver/levenberg_marquardt.h"
#include "metriko/solver/sparse_solver.h"
#include "iter_rounding_common.h"
#include "injective_barrier.h"
#include "iter_rounding_init.h"
#include "iter_rounding_loop.h"

namespace metriko {
class CholeskyWrapper {
    vec<int> outer;
    vec<int> inner;

    template <class M>
    bool same_pattern(const M &A) const {
        return A.isCompressed()
            && std::ssize(outer) == A.outerSize() + 1
            && std::ssize(inner) == A.nonZeros()
            && std::equal(outer.begin(), outer.end(), A.outerIndexPtr())
            && std::equal(inner.begin(), inner.end(), A.innerIndexPtr());
    }


#ifdef GC_HAVE_SUITESPARSE
    cholmod_common c;
    cholmod_factor* L = nullptr;
    Eigen::SparseMatrix<double, Eigen::RowMajor> Jt;   // J row-major: the column-compressed storage of J'

    cholmod_sparse view() {   // J' as a cholmod matrix (unknowns x residuals), unsymmetric: cholmod factorizes A A'
        cholmod_sparse A{};
        A.nrow = Jt.cols(); A.ncol = Jt.rows(); A.nzmax = Jt.nonZeros();
        A.p = Jt.outerIndexPtr(); A.i = Jt.innerIndexPtr(); A.x = Jt.valuePtr();
        A.stype = 0; A.itype = CHOLMOD_INT; A.xtype = CHOLMOD_REAL; A.dtype = CHOLMOD_DOUBLE; A.sorted = 1; A.packed = 1;
        return A;
    }

public:
    CholeskyWrapper() { cholmod_start(&c); c.supernodal = CHOLMOD_SUPERNODAL; }
    ~CholeskyWrapper() { if (L) cholmod_free_factor(&L, &c); cholmod_finish(&c); }
    CholeskyWrapper(const CholeskyWrapper&) = delete;
    CholeskyWrapper& operator=(const CholeskyWrapper&) = delete;

    bool factorize(const SprsD &J) {
        Jt = J;
        Jt.makeCompressed();
        cholmod_sparse A = view();
        if (!L || !same_pattern(Jt)) {
            if (L) cholmod_free_factor(&L, &c);
            L = cholmod_analyze(&A, &c);
            if (!L) { outer.clear(); inner.clear(); return false; }
            outer.assign(Jt.outerIndexPtr(), Jt.outerIndexPtr() + Jt.outerSize() + 1);
            inner.assign(Jt.innerIndexPtr(), Jt.innerIndexPtr() + Jt.nonZeros());
        }
        cholmod_factorize(&A, L, &c);
        return c.status == CHOLMOD_OK && L->minor == L->n;   // minor < n: not positive definite
    }

    bool solve(const VecXd &rhs, VecXd &x) {
        cholmod_dense b{};
        b.nrow = rhs.size(); b.ncol = 1; b.nzmax = rhs.size(); b.d = rhs.size();
        b.x = const_cast<double*>(rhs.data()); b.xtype = CHOLMOD_REAL; b.dtype = CHOLMOD_DOUBLE;
        cholmod_dense* y = cholmod_solve(CHOLMOD_A, L, &b, &c);
        if (!y) return false;
        x = Eigen::Map<VecXd>((double*)y->x, y->nrow);
        cholmod_free_dense(&y, &c);
        return true;
    }
#else
    SparseLLT<SprsD> llt;

public:
    bool factorize(const SprsD &J) {
        SprsD A = J.transpose() * J;
        if (!same_pattern(A)) {
            llt.analyzePattern(A);
            if (llt.info() != Eigen::Success) { outer.clear(); inner.clear(); return false; }
            outer.assign(A.outerIndexPtr(), A.outerIndexPtr() + A.outerSize() + 1);
            inner.assign(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros());
        }
        llt.factorize(A);
        return llt.info() == Eigen::Success;
    }

    bool solve(const VecXd &rhs, VecXd &x) const { x = llt.solve(rhs); return true; }
#endif
};

    inline bool iterative_rounding(
        const VecXi &fixedIdcs,
        const VecXd &fixedVals,
        const SprsD &weightMatrix,
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
            weightMatrix,
            length,
            localInjectivity,
            N,
            n,
            iter_cb,
            &barrier
        );

        //----- Naive Solution -----//
        if (!seamless) {
            // Either it doesn't need loc_inj or it is already loc_inj (because post_iteration has loc_inj check)
            if (!localInjectivity || ni.post_iteration(ni.XF_Small)) { fullx = ni.x0; return true; }
        }

        isLM.init(&llt_ni, &ni, &dd_ni);
        isLM.solve(verbose);

        // todo need to assign minix!
        fullx = (ni.UExt * isLM.x).head(ni.UFull.rows());
        if (!seamless) return true;

        //----- Seamless Solution -----//
        VecXd X_  = isLM.x.head(ni.UFull.cols());
        VecXd F2_ = (ni.UExt * isLM.x).tail(2 * N * nF); // x is changed from F2, but how is this important??
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
