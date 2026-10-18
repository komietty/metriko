//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_INJECTIVE_BARRIER_H
#define METRIKO_INJECTIVE_BARRIER_H
#include "metriko/common/typedef.h"
#include "metriko/solver/matrix_ops.h"

namespace metriko {
    class InjectiveBarrier {
        const int N;
        const int nF;
        MatXi I;
        MatXi J;
        MatXd S;
        MatXd V; // field volume
        vec<vec<std::pair<int, int>>> readers; // for each field entry, the barrier rows (row, slot) reading it, rows ascending
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

            readers.resize(2 * N * nF);
            for (int i = 0; i < I.rows(); i++)
                for (int j = 0; j < I.cols(); j++)
                    readers[J(i, j)].emplace_back(I(i, j), j);
        }

        // diag(barDer) * (the image of the field) * G * w, summed in the order of the sparse products
        SprsD jacobian(const VecXd& barDer, const SprsD& G, double w) const {
            auto merge = [&](int c, vec<std::pair<int, double>>& acc) {
                acc.clear();
                for (SprsD::InnerIterator g(G, c); g; ++g)
                for (auto [r, k]: readers[g.row()]) {
                    double x = barDer(r) * S(r, k) * g.value();
                    auto it = rg::find(acc, r, &std::pair<int, double>::first);
                    if (it == acc.end()) acc.emplace_back(r, x); else it->second += x;
                }
                rg::sort(acc);
            };
            SprsD res;
            sparse_by_columns(N * nF, G.cols(), res,
                [&](int c) { thread_local vec<std::pair<int, double>> acc; merge(c, acc); return (int) acc.size(); },
                [&](int c, int* ri, double* rv) {
                    thread_local vec<std::pair<int, double>> acc;
                    merge(c, acc);
                    for (auto [r, x]: acc) { *ri++ = r; *rv++ = x * w; }
                });
            return res;
        }

        void update(
            const VecXd &field, // todo: does not change??
            const bool updateJ,
            VecXd &fBarrier,
            VecXd &barDer
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
                constexpr double inf = std::numeric_limits<double>::infinity();
                barDer = -splineDerivative.array() / barSpline.array().square();
                barDer = (fBarrier.array().abs() < 10e-9).select(0., (fBarrier.array() == inf).select(inf, barDer));
            }
        }
    };
}
#endif
