//
// Created by saki on 2026/04/09.
//

#ifndef TMESH_H_EMESH_COLLAPSE_EHALF_H
#define TMESH_H_EMESH_COLLAPSE_EHALF_H

#include "./emesh.h"
namespace metriko::tutte {

vec<Half> Emesh::collapse_half_find_path(int ehid) {
    auto& eh = ehalfs[ehid];
    auto& eq = equads[eh.eqid];
    auto  it = rg::find(eq.ehids, ehid);

    auto visit_edge = std::vector(hm.nH, false);

    // 1. 境界エッジを踏まないようにする (既存)
    for (int ehid_b: eq.ehids) {
        for (Half h: ehalfs[ehid_b].halfs) {
            visit_edge[h.id] = true;
            visit_edge[h.twin().id] = true; // 双対も念のため
        }
    }

    // ---------------------------------------------------------
    std::vector allow_vert(hm.nV, false);
    for (int vid: eq.verts_inside()) { allow_vert[vid] = true; }
    assert(it != eq.ehids.end());

    int idx = rg::distance(eq.ehids.begin(), it);
    int len = eq.sides.size();
    int idx_prev = (idx + len - 1) % len;
    int side_curr = eq.sides[idx];
    int side_prev = eq.sides[idx_prev];
    assert(side_curr != side_prev);

    Ehalf& eh_prev = ehalfs[eq.ehids[idx_prev]];
    Vert v0 = eh_prev.halfs.front().tail();
    Vert v1 = eh.halfs.back().head();

    // 始点と終点がAllowlistに入っていることを保証
    allow_vert[v0.id] = true;
    allow_vert[v1.id] = true;

    auto path = compute_dijkstra_for_tquad_temp(hm, visit_edge, allow_vert, v0, v1); // halfedge path

    // debug draw
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t count = 0;

        for (Half& h: path) {
            auto p1 = h.tail().pos();
            auto p2 = h.head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{count, count + 1});
            count += 2;
        }

        auto c = polyscope::registerCurveNetwork("test-"+ std::to_string(ehid), ns, es);
        c->resetTransform();
        c->setRadius(0.002);
    }

    return path;
}

bool Emesh::collapse_half(const int ehid) {
    auto& eh = ehalfs[ehid];
    auto& eq = equads[eh.eqid];
    auto  it = rg::find(eq.ehids, ehid);

    auto visit_edge = std::vector(hm.nH, false);

    // 1. 境界エッジを踏まないようにする (既存)
    for (int ehid_b: eq.ehids) {
        for (Half h: ehalfs[ehid_b].halfs) {
            visit_edge[h.id] = true;
            visit_edge[h.twin().id] = true; // 双対も念のため
        }
    }

    // ---------------------------------------------------------
    std::vector allow_vert(hm.nV, false);
    for (int vid: eq.verts_inside()) { allow_vert[vid] = true; }
    assert(it != eq.ehids.end());

    int idx = rg::distance(eq.ehids.begin(), it);
    int len = eq.sides.size();
    int idx_prev = (idx + len - 1) % len;
    int side_curr = eq.sides[idx];
    int side_prev = eq.sides[idx_prev];
    assert(side_curr != side_prev);

    Ehalf& eh_prev      = ehalfs[eq.ehids[idx_prev]];
    Ehalf& eh_prev_twin = ehalfs[ehalfs[eq.ehids[idx_prev]].twid];
    Vert v0 = eh_prev.halfs.front().tail();
    Vert v1 = eh.halfs.back().head();

    // guarantees bgn/end verts are in the allowed list
    allow_vert[v0.id] = true;
    allow_vert[v1.id] = true;

    vec<Half> path0 = compute_dijkstra_for_tquad_temp(hm, visit_edge, allow_vert, v0, v1); // halfedge path
    vec<Half> path1;
    for (Half& h: path0 | vw::reverse) { path1.emplace_back(h.twin()); }

    // todo 1: erase ehalf of this tquad
    eq.ehids.erase(eq.ehids.begin() + idx);
    eq.sides.erase(eq.sides.begin() + idx);

    // todo 2: Replace prev (and twin of prev) ehalf
    eh_prev.halfs = path0;
    eh_prev_twin.halfs = path1;

    // todo 3: extend ehalf of twin tquad
    int eh_prev_twin_id = eh_prev.twid;
    Equad& eq_twin = equads[eh_prev_twin.eqid];

    auto it_twin = rg::find(eq_twin.ehids, eh_prev_twin_id);
    if (it_twin != eq_twin.ehids.end()) {
        int idx_twin = rg::distance(eq_twin.ehids.begin(), it_twin);
        int insert_idx = (idx_twin + eq_twin.ehids.size() - 1) % eq_twin.ehids.size();
        Ehalf& eh_twin = ehalfs[eq_twin.ehids[insert_idx]];
        eh_twin.extend_next(eh);
    }

    // todo 4: twin of eh is consumed to its tquad
    int eh_twid = eh.twid;
    if (eh_twid >= 0) {
        Ehalf& eh_tw = ehalfs[eh_twid];
        Equad& eq_eh_tw = equads[eh_tw.eqid];

        auto it_eh_tw = rg::find(eq_eh_tw.ehids, eh_twid);
        if (it_eh_tw != eq_eh_tw.ehids.end()) {
            int idx_eh_tw = rg::distance(eq_eh_tw.ehids.begin(), it_eh_tw);

            // eh_tw (v1 -> v_mid) の次のエッジ (v_mid -> ...) を取得
            int next_idx = (idx_eh_tw + 1) % eq_eh_tw.ehids.size();
            int next_eh_id = eq_eh_tw.ehids[next_idx];
            Ehalf& next_eh = ehalfs[next_eh_id];

            // next_eh の先頭に eh_tw の halfs を挿入して前方に拡張 (extend_prev の役割)
            next_eh.halfs.insert(next_eh.halfs.begin(), eh_tw.halfs.begin(), eh_tw.halfs.end());

            // 吸収されて不要になった eh_tw を tquad のリストから削除
            eq_eh_tw.ehids.erase(eq_eh_tw.ehids.begin() + idx_eh_tw);
            eq_eh_tw.sides.erase(eq_eh_tw.sides.begin() + idx_eh_tw);
        }
    }


    // debug draw
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t count = 0;

        for (Half& h: path0) {
            auto p1 = h.tail().pos();
            auto p2 = h.head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{count, count + 1});
            count += 2;
        }

        auto c = polyscope::registerCurveNetwork("test-"+ std::to_string(ehid), ns, es);
        c->resetTransform();
        c->setRadius(0.002);
    }

}

}

#endif //TMESH_H_EMESH_COLLAPSE_EHALF_H
