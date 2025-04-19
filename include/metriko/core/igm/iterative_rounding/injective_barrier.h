//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_INJECTIVE_BARRIER_H
#define METRIKO_INJECTIVE_BARRIER_H
#include "metriko/core/common/typedef.h"

namespace metriko {
    class InjectiveBarrier {
        const int N;
        const int nF;
        MatXi I;
        MatXi J;
        MatXd S;
        MatXd V; // field volume
    public:

        InjectiveBarrier(const int N_, const int nF_, const VecXd& F2) : N(N_), nF(nF_) {
            I.resize(N * nF, 4);
            J.resize(N * nF, 4);
            S.resize(N * nF, 4);
            V.resize(nF, N);
            for (int i = 0; i < nF; i++) {
            for (int j = 0; j < N; j++) {
                I.row(i * N + j) = VecXi::Constant(4, i * N + j);
                J.row(i * N + j) <<
                        2 * N * i + 2 * j,
                        2 * N * i + 2 * j + 1,
                        2 * N * i + 2 * ((j + 1) % N),
                        2 * N * i + 2 * ((j + 1) % N) + 1;
            }}

            for (int i = 0; i < nF; i++) {
            for (int j = 0; j < N; j++) {
                const int oft = 2 * N * i;
                Row2d v0 = F2.segment(oft + 2 * j, 2);
                Row2d v1 = F2.segment(oft + 2 * ((j + 1) % N), 2);
                V(i, j) = v0.norm() * v1.norm();
            }}

        }

        SprsD gen_image_filed(const VecXd &currField) const {
            SprsD m(N * nF, currField.size());
            std::vector<TripD> T;
            for (int i = 0; i < I.rows(); i++)
                for (int j = 0; j < I.cols(); j++)
                    T.emplace_back(I(i, j), J(i, j), S(i, j));
            m.setFromTriplets(T.begin(), T.end());
            return m;
        }

        void update(
            const VecXd &field, // todo: does not change??
            const bool updateJ,
            VecXd &fBarrier,
            SprsD &gBarrier
        ) {
            const double s = 0.5;

            fBarrier = VecXd::Zero(N * nF);
            VecXd barSpline = VecXd::Zero(N * nF);
            VecXd splineDerivative;
            if (updateJ) splineDerivative = VecXd::Zero(N * nF, 1);

            for (int i = 0; i < nF; i++) {
                for (int j = 0; j < N; j++) {
                    Row2d v0 = field.segment(2 * N * i + 2 * j, 2);
                    Row2d v1 = field.segment(2 * N * i + 2 * ((j + 1) % N), 2);
                    double v = V(i, j);
                    double r = (v0(0) * v1(1) - v0(1) * v1(0)) / v; // volume ratio
                    double bar = std::pow(r / s, 3) - 3 * std::pow(r / s, 2) + 3 * r / s;
                    double bar2 = 1. / bar - 1.;
                    if (r <= 0) bar2 = std::numeric_limits<double>::infinity();
                    if (r >= s) bar2 = 0.;
                    fBarrier(N * i + j) = bar2;
                    barSpline(N * i + j) = bar;

                    if (updateJ) {
                        double d = 3. * (r * r / (s * s * s)) - 6. * (r / (s * s)) + 3. / s;
                        if (r <= 0) d = std::numeric_limits<double>::infinity();
                        if (r >= s) d = 0;
                        splineDerivative(N * i + j) = d;
                        S.row(N * i + j) <<
                                v1(1) / v,
                               -v1(0) / v,
                               -v0(1) / v,
                                v0(0) / v;
                    }
                }
            }

            if (updateJ) {
                VecXd barDer = -splineDerivative.array() / (barSpline.array() * barSpline.array()).array();
                for (int i = 0; i < fBarrier.size(); i++)
                    if (std::abs(fBarrier(i)) < 10e-9) barDer(i) = 0.;
                    else if (fBarrier(i) == std::numeric_limits<double>::infinity())
                        barDer(i) = std::numeric_limits<double>::infinity();

                gBarrier.resize(barDer.size(), barDer.size());
                std::vector<TripD> T;
                for (int i = 0; i < barDer.size(); i++)
                    T.emplace_back(i, i, barDer(i));
                gBarrier.setFromTriplets(T.begin(), T.end());
            }

            //gBarrier = gBarrier * gen_image_filed(field);
        }
    };
}
#endif
