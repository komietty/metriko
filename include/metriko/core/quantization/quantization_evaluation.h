//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_EVALUATION_H
#define METRIKO_QUANTIZATION_EVALUATION_H
#include "quantization_basisloop.h"

namespace metriko {

    inline VecXd construct_generating_vector(
        const std::vector<Tedge>& tedges,
        const std::vector<Thalf>& thalfs,
        const std::vector<int>& generating_loop
    ) {
        VecXd V = VecXd::Zero(tedges.size());
        for (int thid: generating_loop) {
            Thalf th = thalfs[thid];
            Tedge te = th.edge();
            V[te.id] += 1;
            assert(V[te.id] <= 2);
        }
        return V;
    }

    template <typename Func>
    MatXd construct_generating_vectors(
        const Tmesh& tmesh,
        const VecXd& R,
        Func compare
    ) {
        int rows = (int)tmesh.tedges.size();
        std::vector<std::pair<int, VecXd>> cache;

        for (const Thalf &th: tmesh.thalfs) {
            if (!th.cannonical) continue;
            auto l = gen_basis_loop(tmesh.tquads, tmesh.thalfs, tmesh.th2quad, tmesh.th2side, R, th, compare);
            auto g = construct_generating_vector(tmesh.tedges, tmesh.thalfs, l);
            cache.emplace_back((g.array() > 0).count(), g);
        }

        std::ranges::sort(cache, [](auto &a, auto &b) { return a.first < b.first; });

        int rank = 0;
        int cols = 0;
        MatXd G;

        for (const auto& g: cache | std::views::values) {
            G.conservativeResize(rows, cols + 1);
            G.col(cols) = g;
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(G);
            if (int r = qr.rank(); r > rank) { rank = r; cols++; }
            else  G.conservativeResize(rows, cols);
        }

        std::cout << "The rank of the matrix is: " << rank << std::endl;

        return G;
    }
}

#endif
