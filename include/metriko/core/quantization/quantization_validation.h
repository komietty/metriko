//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_VALIDATION_H
#define METRIKO_QUANTIZATION_VALIDATION_H
#include "metriko/core/tmesh/tmesh.h"

namespace metriko {
    inline bool compute_validation(
        const Tmesh &tm,
        const VecXd &X,
        const int thid
    ) {
        std::queue<int> q;
        std::vector<bool> visited;
        visited.assign(tm.nTE, false);
        q.emplace(tm.thalfs[thid].id);

        while (!q.empty()) {
            int thid = q.front();
            q.pop();

            const auto &th = tm.thalfs[thid];
            const auto &te = th.edge();
            if (!th.cannonical && te.seg_fr.id == 0) return false; /* is first seg */

            for (auto pair: th.adj_thalfs()) {
                int teid = pair.edge().id;
                if (!visited[teid] && X[teid] == 0) {
                    visited[teid] = true;
                    q.emplace(pair.id);
                }
            }
        }
        return true;
    }

    inline bool compute_validation(
        const Tmesh &tm,
        const VecXd &X
    ) {
        return (X.array() > 0).all() &&
               rg::all_of(
                   tm.thalfs | vw::filter([&](auto th) { return tm.th2sing[th.id] > -1; }),
                   [&](auto &th) { return compute_validation(tm, X, th.id); }
               );
    }
}

#endif
