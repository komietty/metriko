//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_VALIDATION_H
#define METRIKO_QUANTIZATION_VALIDATION_H
#include <queue>
#include "metriko/core/tmesh/tmesh.h"

namespace metriko {
inline bool compute_validation(
    const mc::Mgrph &mg,
    const Tmesh &tm,
    const VecXd &X
) {
    if ((X.array() < 0).any()) return false;

    std::vector<std::vector<int>> adjs(mg.mnodes.size());
    for (auto& te : tm.tedges)
        if (X[te.id] == 0) {
            adjs[te.fr_nid].push_back(te.to_nid);
            adjs[te.to_nid].push_back(te.fr_nid);
        }

    std::vector visited(mg.mnodes.size(), false);

    for (int i = 0; i < mg.mnodes.size(); ++i) {
        if (mg.mnodes[i].jt != mc::JunctionType::F || visited[i]) continue;

        int start_nid = i;
        std::queue<int> q;
        q.push(start_nid);
        visited[start_nid] = true;

        while (!q.empty()) {
            int curr_nid = q.front();
            q.pop();

            if (curr_nid != start_nid && mg.mnodes[curr_nid].jt == mc::JunctionType::F) return false;

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
