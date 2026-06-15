
#ifndef METRIKO_HPATH_H
#define METRIKO_HPATH_H
#include <unordered_set>
#include "hmloc.h"

namespace metriko {

inline vec<HmLoc> approx_shortest_path(
    const int n_div,
    const Hmesh& hm,
    const HmLoc& loc_bgn,
    const HmLoc& loc_end,
    const vec<std::tuple<int, double, double>>& allowed // (eid, r0, r1) allowed region ranges
) {

    auto incidentFaces = [&](const HmLoc& loc) -> vec<int> {
        return std::visit(overloaded{
            [&](const HmLocOnF& l) -> vec<int> { return { l.id }; },
            [&](const HmLocOnH& l) -> vec<int> {
                vec<int> fs; Half h = hm.halfs[l.id];
                if (h.face().id != -1)        fs.push_back(h.face().id);
                if (h.twin().face().id != -1) fs.push_back(h.twin().face().id);
                return fs;
            },
            [&](const HmLocOnE& l) -> vec<int> {
                vec<int> fs; Edge e = hm.edges[l.id];
                if (e.face0().id != -1) fs.push_back(e.face0().id);
                if (e.face1().id != -1) fs.push_back(e.face1().id);
                return fs;
            },
            [&](const HmLocOnV& l) -> vec<int> {
                vec<int> fs;
                for (Half h : hm.verts[l.id].adjHalfs())
                    if (h.face().id != -1) fs.push_back(h.face().id);
                return fs;
            },
            [&](const auto&) -> vec<int> { return {}; },
        }, loc);
    };
    auto incidentEdges = [&](const HmLoc& loc) -> vec<int> {
        return std::visit(overloaded{
            [&](const HmLocOnE& l) -> vec<int> { return { l.id }; },
            [&](const HmLocOnH& l) -> vec<int> { return { hm.halfs[l.id].edge().id }; },
            [&](const HmLocOnV& l) -> vec<int> {
                vec<int> es;
                for (Half h : hm.verts[l.id].adjHalfs()) es.push_back(h.edge().id);
                return es;
            },
            [&](const auto&) -> vec<int> { return {}; },
        }, loc);
    };
    auto shares = [](const vec<int>& a, const vec<int>& b) {
        for (int x : a) for (int y : b) if (x == y) return true;
        return false;
    };
    // same edge or shared face -> straight segment, skip the graph search.
    if (shares(incidentEdges(loc_bgn), incidentEdges(loc_end))) return { loc_bgn, loc_end };
    if (shares(incidentFaces(loc_bgn), incidentFaces(loc_end))) return { loc_bgn, loc_end };

    // if not on the same element, follows

    struct Node { HmLoc loc; Row3d pos; };
    vec<Node> nodes;
    umap<int, vec<int>> e2n;  // eid -> node

    for (auto& [eid, r0, r1] : allowed) {
        Half h = hm.edges[eid].half(); // h is canonical half
        int cnt = std::max(1, (int)std::lround(n_div * (r1 - r0)));
        for (int k = 0; k < cnt; ++k) {
            double r = r0 + (r1 - r0) * (k + 1.) / (cnt + 1.);
            HmLocOnH loc{h.id, r};
            nodes.push_back({HmLoc(loc), get_ptloc_pos(hm, HmLoc(loc))});
            e2n[eid].push_back(nodes.size() - 1);
        }
    }

    // ============================================================
    // 3. construct adjacency list
    // ============================================================
    vec<vec<std::pair<int, double>>> adj;
    adj.resize(nodes.size());
    auto link = [&](int a, int b) {
        double w = (nodes[a].pos - nodes[b].pos).norm();
        adj[a].emplace_back(b, w);
        adj[b].emplace_back(a, w);
    };

    // 同じ面に乗る「異なるエッジ」のノード同士のみ接続（同一エッジ接続は除外）
    for (Face f : hm.faces) {
        const vec<int>* lists[3] = { nullptr, nullptr, nullptr };
        int nl = 0;
        for (Half h : f.adjHalfs()) {
            auto it = e2n.find(h.edge().id);
            if (it != e2n.end()) lists[nl++] = &it->second;
        }
        for (int i = 0; i < nl; ++i)
        for (int j = i + 1; j < nl; ++j)
        for (int a : *lists[i])
        for (int b : *lists[j])
            link(a, b);
    }

    // ============================================================
    // 4. 始点・終点（面内任意点）をノードとして追加
    // ============================================================
    // 端点（vert / edge / face内点）をノード化し、属する面のエッジ候補へ接続
    auto addPoint = [&](const HmLoc& loc) -> int {
        int id = (int)nodes.size();
        nodes.push_back({ loc, get_ptloc_pos(hm, loc) });
        adj.resize(nodes.size());
        std::unordered_set<int> linked;
        for (int fid : incidentFaces(loc))
        for (Half h : hm.faces[fid].adjHalfs())
            if (auto it = e2n.find(h.edge().id); it != e2n.end())
                for (int n : it->second)
                    if (linked.insert(n).second) link(id, n);
        return id;
    };

    int s = addPoint(loc_bgn);
    int t = addPoint(loc_end);

    // ============================================================
    // 5. Dijkstra
    // ============================================================
    constexpr double INF = std::numeric_limits<double>::infinity();
    vec dist(nodes.size(), INF);
    vec prev(nodes.size(), -1);
    std::priority_queue< std::pair<double, int>, vec<std::pair<double, int>>, std::greater<>> pq;
    dist[s] = 0;
    pq.push({ 0., s });
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        if (u == t) break;
        for (auto [v, w] : adj[u])
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({ dist[v], v });
            }
    }

    if (dist[t] == INF) {
        std::cerr << "no path within allowed regions (disconnected by restriction?)" << std::endl;
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
