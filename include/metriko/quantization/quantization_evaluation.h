//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QUANTIZATION_EVALUATION_H
#define METRIKO_QUANTIZATION_EVALUATION_H
#include "quantization_basisloop.h"
#include "metriko/solver/matrix_ops.h"

namespace metriko {
    template <typename Func>
    SprsD construct_generating_vectors(
        const Emesh& tm,
        const VecXi& th2side,
        const VecXd& R,
        Func compare
    ) {
        vec<const Ehalf*> canos;
        for (auto& th: tm.thalfs) if (th.cano) canos.push_back(&th);
        vec<vec<int>> loops(canos.size());

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < canos.size(); ++i)
            loops[i] = gen_basis_loop(tm, th2side, R, *canos[i], compare);

        vec<TripD> T;
        for (int i = 0; i < loops.size(); ++i)
        for (int thid: loops[i])
            T.emplace_back(i, tm.thalfs[thid].teid, 1);

        SprsD G(tm.tedges.size(), tm.tedges.size());
        G.setFromTriplets(T.begin(), T.end());

        for (const auto& [r, c, v]: nonzeros(G))
            METRIKO_CHECK(v <= 2, "generating loop {} passes tedge {} more than twice", r, c);

        reduce_to_linearly_independent(G);

        #if METRIKO_DEBUG
        std::cout << "The rank of the matrix is: " << G.rows() << std::endl;
        #endif

        return G.transpose();
    }
}
#endif
