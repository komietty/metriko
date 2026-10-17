//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_LOOP_H
#define METRIKO_LOOP_H
#include <queue>
#include "utilities.h"
#include "predicates.h"

// loop: a closed simple ccw polyline
namespace metriko {

inline int loop_next(const vec<complex>& loop, int i) { return (i + 1) % loop.size(); }
inline int loop_prev(const vec<complex>& loop, int i) { return (i + loop.size() - 1) % loop.size(); }

inline complex loop_inward(const vec<complex>& loop, int i) {
    auto o = loop[loop_next(loop, i)] - loop[i];
    auto e = loop[loop_prev(loop, i)] - loop[i];
    auto d = -(o / std::abs(o) + e / std::abs(e));
    if (std::abs(d) < 1e-6) d = o / std::abs(o) * complex(0, 1);
    return d / std::abs(d);
}

inline double loop_turn(const vec<complex>& loop, int i) {
    auto a = loop[i] - loop[loop_prev(loop, i)];
    auto b = loop[loop_next(loop, i)] - loop[i];
    return cross(a , b);
}

inline bool loop_into(const vec<complex>& loop, int i, complex d) {
    auto o = loop[loop_next(loop, i)] - loop[i];
    auto e = loop[loop_prev(loop, i)] - loop[i];
    return cross(o, e) > 0 ? cross(o, d) > EPS && cross(d, e) > EPS : cross(o, d) > EPS || cross(d, e) > EPS;
}

// the segment p -> q lies inside the loop. an end sitting on corner ia / ib (-1: none) must leave it inwards
inline bool loop_sees(const vec<complex>& loop, complex p, complex q, int ia, int ib) {
    if (ia >= 0 && !loop_into(loop, ia, q - p)) return false;
    if (ib >= 0 && !loop_into(loop, ib, p - q)) return false;
    for (int k = 0; k < loop.size(); ++k) {
        if (find_strict_intersection(p, q, loop[k], loop[loop_next(loop, k)])) return false;
        if (k != ia && k != ib && is_inside_segment(p, q, loop[k])) return false;
    }
    return true;
}

// the shortest path from corner a to corner b inside the loop (empty when there is no way through). it bends
// only at corners that are not convex, and there it is pushed inside by `hair` so that it touches no boundary
inline vec<complex> shortest_path_in_loop(const vec<complex>& loop, int a, int b, double hair) {
    vec idx = {a, b};
    vec pts = {loop[a], loop[b]};
    if (loop_sees(loop, pts[0], pts[1], a, b)) return pts;
    for (int i = 0; i < loop.size(); ++i)
        if (i != a && i != b && loop_turn(loop, i) < 1e-6) { idx.push_back(i); pts.push_back(loop[i] + loop_inward(loop, i) * hair); }

    const int m = idx.size();
    vec dist(m, std::numeric_limits<double>::infinity());
    vec from(m, -1);
    std::priority_queue<std::pair<double, int>, vec<std::pair<double, int>>, std::greater<>> pq;
    dist[0] = 0;
    pq.emplace(0., 0);
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        if (u == 1) break;
        for (int v = 0; v < m; ++v) {
            const double w = d + std::abs(pts[v] - pts[u]);
            if (w >= dist[v] || !loop_sees(loop, pts[u], pts[v], u < 2 ? idx[u] : -1, v < 2 ? idx[v] : -1)) continue;
            dist[v] = w;
            from[v] = u;
            pq.emplace(w, v);
        }
    }
    if (from[1] == -1) return {};

    vec<complex> path;
    for (int v = 1; v != -1; v = from[v]) path.push_back(pts[v]);
    return path | vw::reverse | rg::to<vec<complex>>();
}
}
#endif
