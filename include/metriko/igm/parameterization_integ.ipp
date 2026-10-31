//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#pragma once
namespace metriko {
    inline bool RosyParameterization::integ() {
        SprsD fullx2Nfunc = vtrans2cut * uncompress;
        SprsD Cfull = constraint * uncompress;
        reduce_to_linearly_independent(Cfull);

        //---- generating G: the matrix for generating a vector field from uv -----
        vec<TripD> eT;
        vec<TripD> iT;
        for (Face f: cut.faces) {
        for (Half h: f.adjHalfs()) {
            Row3d g = f.normal().cross(h.next().vec()) / (f.area() * 2);
            for (int k = 0; k < N; k++) {
                int c  = N * h.tail().id + k;
                int er = 3 * N * f.id + k * 3;
                int ir = 2 * N * f.id + k * 2;
                eT.emplace_back(er + 0, c, g.x());
                eT.emplace_back(er + 1, c, g.y());
                eT.emplace_back(er + 2, c, g.z());
                iT.emplace_back(ir + 0, c, g.dot(f.basisX()));
                iT.emplace_back(ir + 1, c, g.dot(f.basisY()));
            }
        }}
        SprsD G3(3 * N * cut.nF, N * cut.nV);
        SprsD G2(2 * N * cut.nF, N * cut.nV);
        G3.setFromTriplets(eT.begin(), eT.end());
        G2.setFromTriplets(iT.begin(), iT.end());

        MatXd rawF2(cut.nF, 2 * N);
        for (int i = 0; i < N; i++)
            rawF2.middleCols(2 * i, 2) << ext.middleCols(3 * i, 3).cwiseProduct(cut.faceBasisX).rowwise().sum(),
                                          ext.middleCols(3 * i, 3).cwiseProduct(cut.faceBasisY).rowwise().sum();

        auto locally_injectivity_check = [&](const VecXd &X) {
            if (!locInj) return false;
            const VecXd NF = fullx2Nfunc * X;
            const MatXd nfn_ = NF.reshaped<Eigen::RowMajor>(cut.nV, N);
            for (Face f: cut.faces) {
                Half h = f.half();
                Row2d a = nfn_.block(h.tail().id, 0, 1, 2);
                Row2d b = nfn_.block(h.next().tail().id, 0, 1, 2);
                Row2d c = nfn_.block(h.prev().tail().id, 0, 1, 2);
                Row2d d1 = b - a;
                Row2d d2 = c - a;
                if (d1.x() * d2.y() - d1.y() * d2.x() <= 0) return false;
            }
            return true;
        };

        SprsD W = SprsD(raw.faceArea.replicate(1, 2 * N).transpose().reshaped().asDiagonal());
        auto l = (cut.pos.colwise().maxCoeff() - cut.pos.colwise().minCoeff()).norm() * gridscale / avg_norm(ext, 3, N, cut.nF);
        auto s = iterative_rounding(fixedIdcs, fixedVals, W, singularIdcs, integerIdcs, l, Cfull, G2 * fullx2Nfunc, N, n, cut.nF, seamless, roundSeams, locInj, rawF2.reshaped<Eigen::RowMajor>().transpose(), locally_injectivity_check, fullx);

        VecXd NF = fullx2Nfunc * fullx;
        nfn = NF.reshaped<Eigen::RowMajor>(cut.nV, N);

        cfn.resize(raw.nF, N * 3);
        for (int j = 0; j < 3; j++) cfn.middleCols(N * j, N) = nfn(cut.idx.col(j), Eigen::indexing::all);

        return s;
    }
}
