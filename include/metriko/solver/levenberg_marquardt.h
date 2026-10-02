//
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_LEVENBERG_MARQUARDT_H
#define METRIKO_LEVENBERG_MARQUARDT_H
#include <iostream>
#include "matrix_ops.h"

namespace metriko {
    template<class SolverTraits>
    class DiagonalDamping {
    public:
        double currLambda;


        void build_damped_matrix(const SprsD& J, SprsD& dampJ) const {
            VecXd dampVec = VecXd::Zero(J.cols());
            vec<TripD> dampJTris;
            for (int k = 0; k < J.outerSize(); ++k)
                for (SprsD::InnerIterator it(J, k); it; ++it) {
                    dampVec(it.col()) += currLambda * it.value() * it.value();
                    dampJTris.emplace_back(it.row(), it.col(), it.value());
                }
            for (int i = 0; i < dampVec.size(); i++)
                dampJTris.emplace_back(J.rows() + i, i, std::sqrt(dampVec(i)));

            dampJ.resize(J.rows() + J.cols(), J.cols());
            dampJ.setFromTriplets(dampJTris.begin(), dampJTris.end());
        }

        void init( const SprsD& J, SprsD& dampJ) const { build_damped_matrix(J, dampJ); }

        void update(
            const SprsD& J,
            double prvEnergy2,
            double newEnergy2,
            SprsD& dampJ
        ) {
            bool f = prvEnergy2 > newEnergy2 && newEnergy2 != std::numeric_limits<double>::infinity();
            if (f) currLambda /= 10.;
            else   currLambda *= 10.;
            build_damped_matrix(J, dampJ);
        }

        DiagonalDamping(double _currLambda = 0.01) : currLambda(_currLambda) { }
    };

    template<class LinearSolver, class SolverTraits, class DampingTraits>
    class LMSolver {
    public:
        Eigen::VectorXd x; //current solution; always updated
        Eigen::VectorXd prevx; //the solution of the previous iteration
        Eigen::VectorXd x0; //the initial solution to the system
        Eigen::VectorXd d; //the direction taken.
        Eigen::VectorXd currObjective; //the current value of the energy
        Eigen::VectorXd prevObjective; //the previous value of the energy

        LinearSolver *LS;
        SolverTraits *ST;
        DampingTraits *DT;

        double funcTolerance;
        double fooTolerance;

        //always updated to the current iteration
        double energy;
        double fooOptimality;
        int maxIter;
        int curIter;

        LMSolver() { }

        void init(
            LinearSolver *_LS,
            SolverTraits *_ST,
            DampingTraits *_DT,
            int _maxIterations = 100,
            double _funcTolerance = 10e-10,
            double _fooTolerance = 10e-10
        ) {
            LS = _LS;
            ST = _ST;
            DT = _DT;
            maxIter = _maxIterations;
            funcTolerance = _funcTolerance;
            fooTolerance  = _fooTolerance;

            d.resize(ST->xSize);
            x.resize(ST->xSize);
            x0.resize(ST->xSize);
            prevx.resize(ST->xSize);
            currObjective.resize(ST->ESize);
            currObjective.resize(ST->ESize);
        }

        bool solve(const bool verbose) {
            using namespace Eigen;
            ST->initial_solution(x0);
            prevx << x0;

            VectorXd rhs(ST->xSize);
            VectorXd dir;
            if (verbose) std::cout << "******Beginning Optimization******" << std::endl;

            //estimating initial miu
            SprsD dampJ;
            VecXd ECur;
            VecXd Eprv;
            SprsD J;

            curIter = 0;
            ST->objective_jacobian(prevx, ECur, J, true);
            DT->init(J, dampJ);

            do {
                ST->pre_iteration(prevx);
                double prvEnergy2 = ECur.squaredNorm();
                rhs = -(J.transpose() * ECur);

                fooOptimality = rhs.lpNorm<Infinity>();
                if (fooOptimality < fooTolerance) { x = prevx; break; }

                //solving to get the LM direction
                if (!LS->factorize(dampJ)) { std::cout << "Solver Failed to factorize! " << std::endl; return false; }

                LS->solve(rhs, dir);

                if (dir.norm() < funcTolerance) { x = prevx; return true; }

                Eprv = ECur;

                ST->objective_jacobian(prevx + dir, ECur, J, false);
                double newEnergy2 = ECur.squaredNorm();
                energy = newEnergy2;
                bool accepted = prvEnergy2 > newEnergy2;
                x = accepted ? prevx + dir : prevx;

                if (std::abs(prvEnergy2 - newEnergy2) < funcTolerance) { break; }
                if (accepted) ST->objective_jacobian(x, ECur, J, true);
                else ECur = Eprv;

                energy = ECur.squaredNorm();

                DT->update(J, prvEnergy2, newEnergy2, dampJ);

                if (ST->post_iteration(x)) { return true; }

                curIter++;
                prevx = x;
            } while (curIter <= maxIter);

            return false;
        }
    };
}
#endif
