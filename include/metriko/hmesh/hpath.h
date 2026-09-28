//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_HPATH_H
#define METRIKO_HPATH_H
#include <format>
#include <unordered_set>
#include "hmloc.h"

namespace metriko {
inline vec<HmLoc> approx_shortest_path(
    const int n_div,
    const Hmesh& hm,
    const HmLoc& loc_bgn,
    const HmLoc& loc_end,
    const vec<Erng>& rngs
) {
    if (try_get_edge(hm, loc_bgn, loc_end)) return {loc_bgn, loc_end};
    if (try_get_face(hm, loc_bgn, loc_end)) return {loc_bgn, loc_end};

    struct Node { HmLoc loc; Row3d pos; };
    vec<Node> nodes;
    umap<int, vec<int>> e2n;
    std::unordered_set<int> corridor;

    for (auto& r: rngs) {
        corridor.insert(hm.edges[r.id].face0().id);
        corridor.insert(hm.edges[r.id].face1().id);
        int c = std::max(1, (int)std::lround(n_div * r.span()));
        for (int k = 0; k < c; ++k) {
            HmLocOnE l{.id=r.id, .r=r.lerp((k + 1.) / (c + 1.))};
            nodes.emplace_back(l, get_ptloc_pos(hm, l));
            e2n[r.id].push_back(nodes.size() - 1);
        }
    }

    // construct adjacency list
    vec<vec<std::pair<int, double>>> adj;
    adj.resize(nodes.size());
    auto link = [&](int a, int b) {
        double w = (nodes[a].pos - nodes[b].pos).norm();
        adj[a].emplace_back(b, w);
        adj[b].emplace_back(a, w);
    };

    // 同じ面に乗る「異なるエッジ」のノード同士のみ接続（同一エッジ接続は除外）
    for (int fid : corridor) {
        Face f = hm.faces[fid];
        const vec<int>* lists[3] = { nullptr, nullptr, nullptr };
        int nl = 0;
        for (Edge e : f.edges()) {
            if (auto it = e2n.find(e.id); it != e2n.end()) lists[nl++] = &it->second;
        }

        for (int i = 0; i < nl; ++i)
        for (int j = i + 1; j < nl; ++j)
        for (int a : *lists[i])
        for (int b : *lists[j])
            link(a, b);
    }

    // 始点・終点（面内任意点）をノードとして追加
    auto addPoint = [&](const HmLoc& l) -> int {
        int id = nodes.size();
        nodes.emplace_back(l, get_ptloc_pos(hm, l));
        adj.resize(nodes.size());
        std::unordered_set<int> linked;
        for (int fid : get_ptloc_faces(hm, l))
        for (Half h : hm.faces[fid].adjHalfs())
            if (auto it = e2n.find(h.edge().id); it != e2n.end())
                for (int n : it->second)
                    if (linked.insert(n).second) link(id, n);
        return id;
    };

    int s = addPoint(loc_bgn);
    int t = addPoint(loc_end);

    // Dijkstra
    constexpr double INF = std::numeric_limits<double>::infinity();
    vec dist(nodes.size(), INF);
    vec prev(nodes.size(), -1);
    std::priority_queue<std::pair<double, int>, vec<std::pair<double, int>>, std::greater<>> pq;
    dist[s] = 0;
    pq.emplace(0., s);
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        if (u == t) break;
        for (auto [v, w] : adj[u])
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.emplace(dist[v], v);
            }
    }

    if (dist[t] == INF) {
        std::cerr << std::format("no path within allowed regions: nodes {}, links(bgn) {}, links(end) {}\n", nodes.size() - 2, adj[s].size(), adj[t].size());
        return {};
    }

    vec<int> seq;
    for (int v = t; v != -1; v = prev[v]) seq.push_back(v);
    rg::reverse(seq);

    vec<HmLoc> path;
    path.reserve(seq.size());
    for (int v : seq) path.push_back(nodes[v].loc);
    return path;
}
}
#endif
