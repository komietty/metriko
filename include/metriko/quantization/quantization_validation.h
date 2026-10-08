//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QUANTIZATION_VALIDATION_H
#define METRIKO_QUANTIZATION_VALIDATION_H
#include <queue>
#include "metriko/tmesh/emesh.h"

namespace metriko {
inline bool compute_validation(
    const Mgrph &mg,
    const Emesh &tm,
    const VecXd &X
) {
    if ((X.array() < 0).any()) return false;

    vec<vec<int>> adjs(mg.mnodes.size());
    for (const auto& [teid, nids] : tm.live_tedges())
        if (X[teid] == 0) {
            adjs[nids.front()].push_back(nids.back());
            adjs[nids.back()].push_back(nids.front());
        }

    vec visited(mg.mnodes.size(), false);

    for (int i = 0; i < mg.mnodes.size(); ++i) {
        if (mg.mnodes[i].jt != JunctionType::F || visited[i]) continue;

        int start_nid = i;
        std::queue<int> q;
        q.push(start_nid);
        visited[start_nid] = true;

        while (!q.empty()) {
            int curr_nid = q.front();
            q.pop();

            if (curr_nid != start_nid && mg.mnodes[curr_nid].jt == JunctionType::F) return false;

            for (int next_nid : adjs[curr_nid])
                if (!visited[next_nid]) {
                    visited[next_nid] = true;
                    q.push(next_nid);
                }
        }
    }

    return true;
}
}
#endif
