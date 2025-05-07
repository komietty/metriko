//
//--- Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_LEVENBERG_MARQUARDT_H
#define METRIKO_LEVENBERG_MARQUARDT_H
#include <iostream>
#include "matrix_ops.h"

namespace metriko {

    template<class SolverTraits>
    class DiagonalDamping {
    public:
        double currLambda;

        void init(
            const Eigen::SparseMatrix<double> &J,
            const Eigen::VectorXd &initSolution,
            const bool verbose,
            Eigen::SparseMatrix<double> &dampJ
        ) {
            //collecting the diagonal values
            Eigen::VectorXd dampVector = Eigen::VectorXd::Zero(initSolution.size());
            std::vector<Eigen::Triplet<double> > dampJTris;
            for (int k = 0; k < J.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(J, k); it; ++it) {
                    dampVector(it.col()) += currLambda * it.value() * it.value();
                    dampJTris.emplace_back(it.row(), it.col(), it.value());
                }
            }
            for (int i = 0; i < dampVector.size(); i++)
                dampJTris.emplace_back(J.rows() + i, i, sqrt(dampVector(i)));

            dampJ.conservativeResize(J.rows() + dampVector.size(), J.cols());
            dampJ.setFromTriplets(dampJTris.begin(), dampJTris.end());

            if (verbose)
                std::cout << "Initial Lambda: " << currLambda << std::endl;
        }

        bool update(
            SolverTraits &ST,
            const Eigen::SparseMatrix<double> &J,
            const Eigen::VectorXd &currSolution,
            const Eigen::VectorXd &direction,
            const bool verbose,
            Eigen::SparseMatrix<double> &dampJ
        ) {
            Eigen::VectorXd EVec;
            Eigen::SparseMatrix<double> stubJ;
            ST.objective_jacobian(currSolution, EVec, stubJ, false);
            double prevEnergy2 = EVec.squaredNorm();
            ST.objective_jacobian(currSolution + direction, EVec, stubJ, false);
            double newEnergy2 = EVec.squaredNorm();

            if ((prevEnergy2 > newEnergy2) &&
                (newEnergy2 != std::numeric_limits<double>::infinity())) //progress; making it more gradient descent
                currLambda /= 10.;
            else
                currLambda *= 10.;

            if (verbose)
                std::cout << "Current Lambda: " << currLambda << std::endl;
            //collecting the diagonal values
            Eigen::VectorXd dampVector = Eigen::VectorXd::Zero(currSolution.size());
            std::vector<Eigen::Triplet<double> > dampJTris;
            for (int k = 0; k < J.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(J, k); it; ++it) {
                    dampVector(it.col()) += currLambda * it.value() * it.value();
                    dampJTris.emplace_back(it.row(), it.col(), it.value());
                }
            }
            for (int i = 0; i < dampVector.size(); i++)
                dampJTris.emplace_back(J.rows() + i, i, sqrt(dampVector(i)));

            dampJ.conservativeResize(J.rows() + dampVector.size(), J.cols());
            dampJ.setFromTriplets(dampJTris.begin(), dampJTris.end());

            return prevEnergy2 > newEnergy2;
            //this preconditioner always approves new direction
        }

        DiagonalDamping(double _currLambda = 0.01) : currLambda(_currLambda) { }
        ~DiagonalDamping() {
        }
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

        int maxIterations;
        double funcTolerance;
        double fooTolerance;

        //always updated to the current iteration
        double energy;
        double fooOptimality;
        int currIter;

    public:
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
            maxIterations = _maxIterations;
            funcTolerance = _funcTolerance;
            fooTolerance = _fooTolerance;

            d.resize(ST->xSize);
            x.resize(ST->xSize);
            x0.resize(ST->xSize);
            prevx.resize(ST->xSize);
            currObjective.resize(ST->ESize);
            currObjective.resize(ST->ESize);
        }


        bool solve(const bool verbose) {
            using namespace Eigen;
            using namespace std;
            ST->initial_solution(x0);
            prevx << x0;

            VectorXd rhs(ST->xSize);
            VectorXd direction;
            if (verbose)
                cout << "******Beginning Optimization******" << endl;

            //estimating initial miu
            SparseMatrix<double> dampJ;
            VectorXd EVec;
            SparseMatrix<double> J;

            currIter = 0;
            ST->objective_jacobian(prevx, EVec, J, true);
            DT->init(J, prevx, verbose, dampJ);

            do {
                ST->pre_iteration(prevx);

                if (verbose) cout << "Initial objective for Iteration " << currIter << ": " << EVec.squaredNorm() << endl;

                //multiply_adjoint_vector(ST->JRows, ST->JCols, JVals, -EVec, rhs);
                rhs = -(J.transpose() * EVec);

                fooOptimality = rhs.template lpNorm<Infinity>();
                if (verbose) cout << "firstOrderOptimality: " << fooOptimality << endl;

                if (fooOptimality < fooTolerance) {
                    x = prevx;
                    if (verbose) {
                        cout << "First-order optimality has been reached" << endl;
                        break;
                    }
                }

                //trying to do A'*A manually
                //SparseMatrix<double,RowMajor> Jt=dampJ.transpose();
                //SparseMatrix<double> JtJ = Jt*dampJ;

                //solving to get the LM direction
                if (!LS->factorize(dampJ.transpose() * dampJ)) {
                    cout << "Solver Failed to factorize! " << endl;
                    return false;
                }

                LS->solve(rhs, direction);

                if (verbose) cout << "direction magnitude: " << direction.norm() << endl;

                if (direction.norm() < funcTolerance) {
                    x = prevx;
                    if (verbose) cout << "Stopping since direction magnitude small." << endl;
                    return true;
                }

                ST->objective_jacobian(prevx, EVec, J, false);
                double prevEnergy2 = EVec.squaredNorm();
                ST->objective_jacobian(prevx + direction, EVec, J, false);
                double newEnergy2 = EVec.squaredNorm();

                energy = newEnergy2;

                if (prevEnergy2 > newEnergy2) {
                    x = prevx + direction; {
                        if (std::abs(prevEnergy2 - newEnergy2) < funcTolerance) {
                            if (verbose) cout << "Stopping sincefunction didn't change above tolerance." << endl;
                            break;
                        }
                    }
                } else x = prevx;

                if (verbose) cout << "New energy: " << energy << endl;
                ST->objective_jacobian(x, EVec, J, true);
                energy = EVec.squaredNorm();

                DT->update(*ST, J, prevx, direction, verbose, dampJ);

                //The SolverTraits can order the optimization to stop by giving "true" of to continue by giving "false"
                if (ST->post_iteration(x)) {
                    if (verbose) cout << "ST->Post_iteration() gave a stop" << endl;
                    return true;
                }
                currIter++;
                prevx = x;
            } while (currIter <= maxIterations);

            return false;
        }
    };
}

#endif
