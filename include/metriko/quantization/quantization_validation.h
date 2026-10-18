//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QUANTIZATION_VALIDATION_H
#define METRIKO_QUANTIZATION_VALIDATION_H
#include <numeric>
#include "metriko/tmesh/emesh.h"

namespace metriko {
// the tedges of length zero merge their end nodes: X is valid when no two fixed nodes end up merged
inline bool compute_validation(
    const Mgrph &mg,
    const Emesh &tm,
    const VecXd &X
) {
    if ((X.array() < 0).any()) return false;

    // 1: the nodes on a tedge of length zero
    vec<std::pair<int, int>> zeros;
    vec<int> nodes;
    for (const auto& [teid, nids] : tm.live_tedges()) {
        if (X[teid] != 0) continue;
        zeros.emplace_back(nids.front(), nids.back());
        nodes.push_back(nids.front());
        nodes.push_back(nids.back());
    }
    rg::sort(nodes);
    nodes.erase(rg::unique(nodes).begin(), nodes.end());
    auto at = [&](int nid) { return int(rg::lower_bound(nodes, nid) - nodes.begin()); };

    // 2: merge them
    vec<int> root(nodes.size());
    std::iota(root.begin(), root.end(), 0);
    auto find = [&](int i) { while (root[i] != i) i = root[i] = root[root[i]]; return i; };
    for (auto [a, b]: zeros) root[find(at(a))] = find(at(b));

    // 3: at most one fixed node per merged group
    vec<bool> fixed(nodes.size(), false);
    for (int i = 0; i < std::ssize(nodes); ++i) {
        if (mg.mnodes[nodes[i]].jt != JunctionType::F) continue;
        int r = find(i);
        if (fixed[r]) return false;
        fixed[r] = true;
    }
    return true;
}
}
#endif
