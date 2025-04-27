//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_ITER_ROUNDING_COMMON_H
#define METRIKO_ITER_ROUNDING_COMMON_H
#include "metriko/core/common/typedef.h"
#include "metriko/core/solver/matrix_ops.h"
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
                SprsD gBarrier;
                VecXd fBarrier;
                barrier->update(field, updateJ, fBarrier, gBarrier);
                E.resize(fInteg.size() + fClose.size() + fConst.size() + fBarrier.size());
                E << fInteg * wInteg, fClose * wClose, fConst * wConst, fBarrier * wBarrier;
                if (updateJ) {
                    gBarrier = loop?
                        gBarrier * barrier->gen_image_filed(field) * X2F * wBarrier :
                        gBarrier * barrier->gen_image_filed(field) * gClose * wInteg; // maybe wBarrier??
                    MatXi I(4, 1);
                    I << 0, 1, 2, 3;
                    sparse_block(I, {&gInteg, &gClose, &gConst, &gBarrier}, J);
                }
                return;
            }

            E.resize(fInteg.size() + fClose.size() + fConst.size());
            E << fInteg * wInteg, fClose * wClose, fConst * wConst;
            if (updateJ) {
                MatXi I(3, 1);
                I << 0, 1, 2;
                sparse_block(I, {&gInteg, &gClose, &gConst}, J);
            }
        }
    };
}
#endif
