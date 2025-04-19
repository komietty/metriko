//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_ITER_ROUNDING_INIT_H
#define METRIKO_ITER_ROUNDING_INIT_H
#include <set>
#include <igl/local_basis.h>
#include <igl/unique.h>
#include <igl/setdiff.h>
#include <igl/speye.h>
#include <igl/slice.h>
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

            std::vector<TripD> T;

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
            const double length_,
            const bool locinj_,
            const int nF,
            const int N,
            const int n,
            std::function<bool(const VecXd &x)> iter_cb,
            InjectiveBarrier *barrier
        ): Integration(G2_, F2_, N, n, locinj_, length_, 10e3, 10e3, 1, 1e-4, barrier),
           fixedIdcs(fixedIdcs),
           fixedVals(fixedVals),
           iter_cb(iter_cb)
        {

            // ---
            // Reducing constraint matrix:
            // reducing duplication of Cfull
            // all non-zero columns
            // ---
            VecXi I(Cfull.nonZeros());
            VecXi J(Cfull.nonZeros());
            VecXd S(Cfull.nonZeros());
            std::set<int> uniqueJ;
            int counter = 0;
            for (int k = 0; k < Cfull.outerSize(); ++k) {
                for (SprsD::InnerIterator it(Cfull, k); it; ++it) {
                    I(counter) = it.row();
                    J(counter) = it.col();
                    uniqueJ.insert(it.col());
                    S(counter++) = it.value();
                }
            }

            // creating small dense matrix with all non-zero columns
            VecXi JMask = VecXi::Constant(Cfull.cols(), -1);
            std::vector uj(uniqueJ.begin(), uniqueJ.end());
            VecXi uniqueJVec = Eigen::Map<VecXi>(uj.data(), uj.size());
            for (int i = 0; i < uj.size(); i++) JMask(uj[i]) = i;

            MatXd CSmall = MatXd::Zero(Cfull.rows(), JMask.maxCoeff() + 1);
            for (int i = 0; i < I.size(); i++) CSmall(I(i), JMask(J(i))) = S(i);

            // converting into the big matrix
            VecXi nonPartIndices, stub;
            VecXi n_vtrans_ids(Cfull.cols());
            for (int i = 0; i < Cfull.cols(); i++) n_vtrans_ids(i) = i;
            igl::setdiff(n_vtrans_ids, uniqueJVec, nonPartIndices, stub);

            MatXd USmall(0, 0);
            if (CSmall.rows() != 0) {
                Eigen::FullPivLU<MatXd> lu(CSmall);
                USmall = lu.kernel(); // lu.rank(): num of the constraints
            } else
                nonPartIndices = n_vtrans_ids;

            SprsD URaw(
                nonPartIndices.size() + USmall.rows(),
                nonPartIndices.size() + USmall.cols()
            );
            std::vector<TripD> urT;
            for (int i = 0; i < nonPartIndices.size(); i++) urT.emplace_back(i, i, 1.);

            for (int i = 0; i < USmall.rows(); i++)
            for (int j = 0; j < USmall.cols(); j++)
                urT.emplace_back(nonPartIndices.size() + i, nonPartIndices.size() + j, USmall(i, j));

            URaw.setFromTriplets(urT.begin(), urT.end());

            SprsD permMat(URaw.rows(), URaw.rows());
            std::vector<TripD> pmT;
            for (int ci = 0; ci < nonPartIndices.size(); ci++) pmT.emplace_back(nonPartIndices(ci), ci, 1.);
            for (int ci = 0; ci < uniqueJVec.size(); ci++) pmT.emplace_back(
                uniqueJVec(ci), nonPartIndices.size() + ci, 1.);
            permMat.setFromTriplets(pmT.begin(), pmT.end());

            UFull = permMat * URaw;

            //----- Generating naive poisson solution -----
            // Compute poisson eq. so that the gradient of enegy function equal to zero.
            // Conceptually it computes uv to follow the given vectorfield with the constraint.
            // See eq. (6) in the report by Bommes(2012)
            //
            // if isometricity is not important, then the solver below might have room for optimization
            // e.g. consider only the conformality... use conformal optimization
            //
            // maybe it is better to use correct mass matrix ...
            // double max_mass = 0;
            // for (Face f: mesh.faces) { max_mass = std::max(max_mass, f.area() * 0.5); }
            // std::vector<TripD> T;
            // for (Face f: mesh.faces) {
            // for (int i = 0; i < N * 2; i++) {
            //     T.emplace_back(f.id * N * 2 + i, f.id * N * 2 + i, f.area() / max_mass);
            // }}

            SprsD Mass;
            igl::speye(2 * N * nF, Mass);
            X2F = (G2 * UFull).pruned();
            SprsD E = X2F.transpose() * Mass * X2F * length;
            VecXd f = X2F.transpose() * Mass * F2;
            SprsD constMat(fixedIdcs.size(), UFull.cols());

            igl::slice(UFull, fixedIdcs, 1, constMat);

            std::vector<TripD> bmT;

            for (int k = 0; k < E.outerSize(); ++k)
            for (SprsD::InnerIterator it(E, k); it; ++it)
                bmT.emplace_back(it.row(), it.col(), it.value());

            for (int k = 0; k < constMat.outerSize(); ++k) {
            for (SprsD::InnerIterator it(constMat, k); it; ++it) {
                bmT.emplace_back(it.row() + E.rows(), it.col(), it.value());
                bmT.emplace_back(it.col(), it.row() + E.rows(), it.value());
            }}

            SprsD bigMat(E.rows() + constMat.rows(), E.rows() + constMat.rows());
            bigMat.setFromTriplets(bmT.begin(), bmT.end());

            VecXd bigRhs(f.size() + fixedVals.size());
            bigRhs << f, fixedVals;

            Eigen::SparseLU<SprsD> solver;
            solver.compute(bigMat);
            if (solver.info() != Eigen::Success) {
                std::cout << "Solver failed..." << std::endl;
                return;
            }

            VecXd XSmallFull = solver.solve(bigRhs);
            VecXd XSmall = XSmallFull.head(UFull.cols());

            x0 = UFull * XSmall;
            XF_Small.resize(XSmall.size() + F2.size());
            XF_Small << XSmall, F2;

            std::vector<TripD> ueT;
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
