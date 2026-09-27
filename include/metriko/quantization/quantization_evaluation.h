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
        const Tmesh& tmesh,
        const VecXd& R,
        Func compare
    ) {
        vec<const Thalf*> canos;
        for (auto& th: tmesh.thalfs) if (th.cano) canos.push_back(&th);
        vec<vec<int>> loops(canos.size());

        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < canos.size(); ++i)
            loops[i] = gen_basis_loop(tmesh.tquads, tmesh.thalfs, tmesh.th2quad, tmesh.th2side, R, *canos[i], compare);

        vec<TripD> T;
        for (int i = 0; i < loops.size(); ++i)
        for (int thid: loops[i])
            T.emplace_back(i, tmesh.thalfs[thid].teid, 1);

        SprsD G(tmesh.tedges.size(), tmesh.tedges.size());
        G.setFromTriplets(T.begin(), T.end());

        for (int k = 0; k < G.outerSize(); ++k)
        for (SprsD::InnerIterator it(G, k); it; ++it)
            METRIKO_CHECK(it.value() <= 2, "generating loop {} passes tedge {} more than twice", it.row(), it.col());

        reduce_to_linearly_independent(G);

        #if METRIKO_DEBUG
        std::cout << "The rank of the matrix is: " << G.rows() << std::endl;
        #endif

        return G.transpose();
    }
}
#endif
