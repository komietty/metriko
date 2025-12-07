//
// Created by saki on 2025/12/07.
//

#ifndef TMESH_H_COMPUTE_DIJKSTRA_H
#define TMESH_H_COMPUTE_DIJKSTRA_H
#include <queue>

namespace metriko {
// need some refactor
std::vector<int> compute_dijkstra(
    const Hmesh& hm,
    const std::vector<bool>& visit,
    Vert v0,
    Vert v1
    ) {
    if (v0 == v1) return {};

    using len_half = std::tuple<double, int>;
    std::priority_queue<len_half, std::vector<len_half>, std::greater<>> pq;
    std::unordered_map<int, int> incoming_v2h;

    auto discovered = [&](Vert v) { return v == v0 || incoming_v2h.contains(v.id); };

    auto enqueue_iH = [&](Vert v, double d) {
        for (Half h: v.adjHalfs()) {
            Vec3d p1 = h.tail().pos();
            Vec3d p2 = h.head().pos();
            double mag = sqrt((p1 - p2).dot(p1 - p2));
            if (!discovered(h.head()) && !visit[h.id] && !visit[h.twin().id]) {
                pq.emplace(d + mag, h.id);
            }
        }
    };

    enqueue_iH(v0, 0.);

    while (!pq.empty()) {
        double dist = get<0>(pq.top());
        int iHc = get<1>(pq.top());
        int iVc = hm.head[iHc];
        pq.pop();
        if (discovered(hm.verts[iVc])) continue;
        incoming_v2h[iVc] = iHc;
        if (iVc == v1.id) {
            std::vector<int> path;
            int iV = iVc;
            while (iV != v0.id) {
                int iHp = incoming_v2h[iV];
                path.push_back(iHp);
                iV = hm.tail[iHp];
            }
            rg::reverse(begin(path), end(path));
            return path;
        }
        enqueue_iH(hm.verts[iVc], dist);
    }
    return {};
}
}

#endif
