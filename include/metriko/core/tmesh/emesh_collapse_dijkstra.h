#ifndef TMESH_H_COMPUTE_DIJKSTRA_H
#define TMESH_H_COMPUTE_DIJKSTRA_H
#include "./emesh.h"

namespace metriko {

constexpr double SAME_FACE_PENALTY = 10000.0;

// need some refactor
inline vec<int> compute_dijkstra(
    const Hmesh& hm,
    const vec<bool>& visit,
    const vec<bool>& allow_vert,
    Vert v0,
    Vert v1
) {
    if (v0 == v1) return {};

    using len_half = std::tuple<double, int>;
    std::priority_queue<len_half, vec<len_half>, std::greater<>> pq;
    std::unordered_map<int, int> incoming_v2h;

    auto discovered = [&](Vert v) { return v == v0 || incoming_v2h.contains(v.id); };

    auto enqueue_iH = [&](Vert v, double d, int incoming_iH) {
        int prev_fid = -1;
        if (incoming_iH != -1) { prev_fid = hm.face[incoming_iH]; }

        for (Half h: v.adjHalfs()) {
            Vec3d p1 = h.tail().pos();
            Vec3d p2 = h.head().pos();
            auto mag = sqrt((p1 - p2).dot(p1 - p2));
            auto pen = 0.;
            auto fid = h.face().id;

            // 変更点2: 同じFaceを連続して通る場合にペナルティを加算
            // (境界エッジなどでFace IDが -1 になる場合を考慮して -1 は除外)
            if (prev_fid != -1 && fid != -1 && prev_fid == fid) { pen = SAME_FACE_PENALTY; }
            if (!discovered(h.head()) && !visit[h.id] && !visit[h.twin().id] && allow_vert[h.head().id]) {
                pq.emplace(d + mag + pen, h.id);
            }
        }
    };

    enqueue_iH(v0, 0., -1);

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
        enqueue_iH(hm.verts[iVc], dist, iHc);
    }
    return {};
}

inline vec<Half> compute_dijkstra_for_tquad_temp(
    const Hmesh& hm,
    const vec<bool>& allow_vert,
    Vert v0,
    Vert v1
) {
    if (v0 == v1) return {};

    vec visit_edge = vec(hm.nH, false);

    using len_half = std::tuple<double, int>;
    std::priority_queue<len_half, vec<len_half>, std::greater<>> pq;
    std::unordered_map<int, int> incoming_v2h;

    auto discovered = [&](Vert v) { return v == v0 || incoming_v2h.contains(v.id); };

    auto enqueue_iH = [&](Vert v, double d, int incoming_iH) {

        int prev_fid = -1;
        if (incoming_iH != -1) { prev_fid = hm.face[incoming_iH]; }

        for (Half h: v.adjHalfs()) {
            // 条件: 境界エッジでなく、かつ、行き先の頂点が「許可リスト」に入っていること
            if (!discovered(h.head()) &&
                !visit_edge[h.id] &&
                !visit_edge[h.twin().id] &&
                allow_vert[h.head().id]) { // ← 変更: 許可されている頂点のみ進む

                Vec3d p1 = h.tail().pos();
                Vec3d p2 = h.head().pos();
                auto mag = sqrt((p1 - p2).dot(p1 - p2));
                auto pen = 0.;
                auto fid = h.face().id;
                if (prev_fid != -1 && fid != -1 && prev_fid == fid) { pen = SAME_FACE_PENALTY; }
                pq.emplace(d + mag + pen, h.id);
            }
        }
    };

    enqueue_iH(v0, 0., -1);

    while (!pq.empty()) {
        double dist = get<0>(pq.top());
        int iHc = get<1>(pq.top());
        int iVc = hm.head[iHc];
        pq.pop();

        if (discovered(hm.verts[iVc])) continue;
        incoming_v2h[iVc] = iHc;

        if (iVc == v1.id) {
            std::vector<Half> path;
            int iV = iVc;
            while (iV != v0.id) {
                int iHp = incoming_v2h[iV];
                path.push_back(hm.halfs[iHp]);
                iV = hm.tail[iHp];
            }
            rg::reverse(path);
            return path;
        }
        enqueue_iH(hm.verts[iVc], dist, iHc);
    }
    return {};
}

}

#endif
