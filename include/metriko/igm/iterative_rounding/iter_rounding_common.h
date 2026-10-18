//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_ITER_ROUNDING_COMMON_H
#define METRIKO_ITER_ROUNDING_COMMON_H
#include "metriko/common/typedef.h"
#include "metriko/solver/matrix_ops.h"
#include "injective_barrier.h"

namespace metriko {
    class Integration {
    public:
        const SprsD &G2;
        const VecXd &F2;
        int xSize = 0;
        int ESize = 0;
        const int N;
        const int n;
        const double length; // grid scale
        const double wInteg;
        const double wConst; // a factor so that the fixed value do not move
        const double wClose;
        const double wBarrier;
        const bool localInjectivity;
        InjectiveBarrier *barrier;
        SprsD gInteg;
        SprsD gClose;
        SprsD gConst;
        SprsD X2F; // convert x to the field with parameter length

        Integration(
            const SprsD &intrinsicGrad_,
            const VecXd &intrinsicField_,
            const int N_,
            const int n_,
            const bool locInj_,
            const double length_,
            const double wInteg_,
            const double wConst_,
            const double wClose_,
            const double wBarrier_,
            InjectiveBarrier *barrier_
        ) : G2(intrinsicGrad_),
            F2(intrinsicField_),
            N(N_),
            n(n_),
            localInjectivity(locInj_),
            length(length_),
            wInteg(wInteg_),
            wConst(wConst_),
            wClose(wClose_),
            wBarrier(wBarrier_),
            barrier(barrier_)
        {
        }

        static void pre_iteration(const VecXd &prevx) { }

        void compute_jacobian(
            const VecXd &field,
            const VecXd& fInteg,
            const VecXd& fClose,
            const VecXd& fConst,
            const bool updateJ,
            VecXd &E,
            SprsD &J,
            const bool loop
        ) {
            if (localInjectivity) {
                VecXd fBarrier, barDer;
                barrier->update(field, updateJ, fBarrier, barDer);
                E.resize(fInteg.size() + fClose.size() + fConst.size() + fBarrier.size());
                E << fInteg * wInteg, fClose * wClose, fConst * wConst, fBarrier * wBarrier;
                if (updateJ) {
                    // the naive branch scales the barrier by wInteg here while the residual above uses
                    // wBarrier, so this J is not the exact derivative of E. it is load-bearing: the barrier
                    // then reaches the normal equations with weight wInteg * wBarrier = 1 in the gradient
                    // but wInteg^2 = 1e8 in the gauss-newton hessian, acting as curvature rather than as a
                    // driving force. using wBarrier here instead leaves it 1e8 times weaker and nefertiti
                    // at gridscale 0.02 hits the iteration cap with 10 inverted faces. the only consistent
                    // variant that still converges raises wBarrier to wInteg, at 24 instead of 13 iterations
                    SprsD gBarrier = loop?
                        barrier->jacobian(barDer, X2F, wBarrier) :
                        barrier->jacobian(barDer, gClose, wInteg); // maybe wBarrier??
                    sparse_vstack({&gInteg, &gClose, &gConst, &gBarrier}, J);
                }
                return;
            }

            E.resize(fInteg.size() + fClose.size() + fConst.size());
            E << fInteg * wInteg, fClose * wClose, fConst * wConst;
            if (updateJ) {
                sparse_vstack({&gInteg, &gClose, &gConst}, J);
            }
        }
    };
}
#endif
