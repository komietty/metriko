//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_ITER_ROUNDING_INIT_H
#define METRIKO_ITER_ROUNDING_INIT_H
#include <stdexcept>
#include <igl/local_basis.h>
#include <igl/slice.h>
#include <igl/speye.h>
#include <igl/unique.h>
#include "iter_rounding_common.h"

namespace metriko {
    class NaiveIntegration : public Integration {
    public:
        const VecXi& fixedIdcs;
        const VecXd& fixedVals;
        SprsD UFull; // parmMat * URaw.
        SprsD UExt;  // direct sum of UFull and rawF2. Just pass through rawF2.
        VecXd XF_Small;
        VecXd x0;
        std::function<bool(const VecXd &x)> iter_cb;

        void initial_solution(VecXd &xf_) const { xf_ = XF_Small; }

        bool post_iteration(const VecXd &xf_) const { return iter_cb((UExt * xf_).head(x0.size())); }

        void prepare_jacobian_component() {
            const int l1 = F2.size();
            const int l2 = UExt.rows();
            const int l3 = fixedIdcs.size();

            gInteg.resize(G2.rows(), G2.cols() + l1);
            gClose.resize(l1, l2);
            gConst.resize(l3, l2);

            vec<TripD> T;

            for (int k = 0; k < G2.outerSize(); ++k)
                for (SprsD::InnerIterator it(G2, k); it; ++it)
                    T.emplace_back(it.row(), it.col(), -length * it.value());
            for (int i = 0; i < l1; i++) T.emplace_back(i, G2.cols() + i, 1.);
            gInteg.setFromTriplets(T.begin(), T.end());

            T.clear();
            for (int i = 0; i < l1; i++) T.emplace_back(i, x0.size() + i, 1.);
            gClose.setFromTriplets(T.begin(), T.end());

            T.clear();
            for (int i = 0; i < l3; i++) T.emplace_back(i, fixedIdcs(i), 1.);
            gConst.setFromTriplets(T.begin(), T.end());

            gInteg = gInteg * UExt * wInteg;
            gClose = gClose * UExt * wClose;
            gConst = gConst * UExt * wConst;
        }

        void objective_jacobian(const VecXd &xf_, VecXd &E, SprsD &J, const bool updateJ) {
            const VecXd xf = UExt * xf_; // here UExt...
            const VecXd x = xf.head(x0.size());
            const VecXd f = xf.tail(F2.size());
            const VecXd fInteg = f - length * G2 * x;
            const VecXd fClose = f - F2;
            //const VecXd fInteg = F2 - length * G2 * xCurr;
            //const VecXd fClose = VecXd::Zero(F2.size());
            VecXd fConst(fixedIdcs.size());
            for (int i = 0; i < fixedIdcs.size(); i++)
                fConst(i) = x(fixedIdcs(i)) - fixedVals(i);
            compute_jacobian(f, fInteg, fClose, fConst, updateJ, E, J, false);
        }

        NaiveIntegration(
            const SprsD &G2_,
            const SprsD &Cfull,
            const VecXd &F2_,
            const VecXi &fixedIdcs,
            const VecXd &fixedVals,
            const SprsD &weightMatrix, // per face weight of the poisson energy, repeated over the 2N rows of the face
            const double length_,
            const bool locinj_,
            const int N,
            const int n,
            std::function<bool(const VecXd &x)> iter_cb,
            InjectiveBarrier *barrier
        ): Integration(G2_, F2_, N, n, locinj_, length_, 10e3, 10e3, 1, 1e-4, barrier),
           fixedIdcs(fixedIdcs),
           fixedVals(fixedVals),
           iter_cb(iter_cb)
        {
            //----- Generating naive poisson solution -----
            // Compute poisson eq. so that the gradient of enegy function equal to zero.
            // Conceptually it computes uv to follow the given nvec with the constraint.
            // See eq. (6) in the report by Bommes(2012)
            // if isometricity is not important, then the solver below might have room for optimization
            // e.g. consider only the conformality... use conformal optimization

            UFull = sparse_null_space(Cfull);
            X2F = (G2 * UFull).pruned();
            SprsD E = X2F.transpose() * weightMatrix * X2F * length;
            VecXd f = X2F.transpose() * weightMatrix * F2;
            SprsD constMat(fixedIdcs.size(), UFull.cols());

            igl::slice(UFull, fixedIdcs, 1, constMat);

            vec<TripD> bmT;

            for (int k = 0; k < E.outerSize(); ++k)
            for (SprsD::InnerIterator it(E, k); it; ++it)
                bmT.emplace_back(it.row(), it.col(), it.value());

            for (int k = 0; k < constMat.outerSize(); ++k) {
            for (SprsD::InnerIterator it(constMat, k); it; ++it) {
                bmT.emplace_back(it.row() + E.rows(), it.col(), it.value());
                bmT.emplace_back(it.col(), it.row() + E.rows(), it.value());
            }}

            // the fixed values are constraints C x = v on a few variables only: write x = x0 + Z y, with Z the identity
            // on the untouched variables and the (dense, small) kernel of C on the touched ones. E restricted to y is
            // symmetric positive definite and is solved by cholesky; the LU of the full KKT system is the fallback
            VecXd XSmall;
            {
                VecXd xp;
                const SprsD Z   = sparse_null_space(constMat, fixedVals, xp);
                const SprsD Ey  = Z.transpose() * E * Z;
                const VecXd rhs = Z.transpose() * (f - E * xp);
                SparseLLT llt(Ey);
                const VecXd y = llt.solve(rhs);
                // a (nearly) singular E passes the factorization but not the solve
                if (llt.info() == Eigen::Success && y.allFinite() && (Ey * y - rhs).norm() <= 1e-8 * std::max(rhs.norm(), 1.)) {
                    XSmall = xp + Z * y;
                } else {
                    SprsD bigMat(E.rows() + constMat.rows(), E.rows() + constMat.rows());
                    bigMat.setFromTriplets(bmT.begin(), bmT.end());
                    VecXd bigRhs(f.size() + fixedVals.size());
                    bigRhs << f, fixedVals;
                    Eigen::SparseLU solver(bigMat);
                    if (solver.info() != Eigen::Success) METRIKO_FAIL("initial Poisson solve (SparseLU) failed");
                    VecXd XSmallFull = solver.solve(bigRhs);
                    XSmall = XSmallFull.head(UFull.cols());
                }
            }

            x0 = UFull * XSmall;
            XF_Small.resize(XSmall.size() + F2.size());
            XF_Small << XSmall, F2;

            vec<TripD> ueT;
            for (int k = 0; k < UFull.outerSize(); ++k)
            for (SprsD::InnerIterator it(UFull, k); it; ++it)
                ueT.emplace_back(it.row(), it.col(), it.value());

            for (int k = 0; k < F2.size(); k++)
                ueT.emplace_back(UFull.rows() + k, UFull.cols() + k, 1.);

            UExt.resize(UFull.rows() + F2.size(), UFull.cols() + F2.size());
            UExt.setFromTriplets(ueT.begin(), ueT.end());

            VecXd E_;
            SprsD J_;
            objective_jacobian(XF_Small, E_, J_, false);
            ESize = E_.size();
            xSize = UExt.cols();

            prepare_jacobian_component();
        }
    };
}
#endif
