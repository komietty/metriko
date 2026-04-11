//
// Created by saki on 2026/04/09.
//

#ifndef TMESH_H_EMESH_COLLAPSE_EHALF_H
#define TMESH_H_EMESH_COLLAPSE_EHALF_H

#include "./emesh.h"
namespace metriko::tutte {

constexpr auto circular_prev = [](auto& c, auto it) { return it == c.begin() ? std::prev(c.end()) : std::prev(it); };
constexpr auto circular_next = [](auto& c, auto it) { auto n = std::next(it); return n == c.end() ? c.begin() : n; };

inline vec<Half> Emesh::collapse_half_find_path(int ehid) {
    auto& eh = ehalfs[ehid];
    auto& eq = equads[eh.eqid];
    auto  it = rg::find(eq.edata, ehid, &Edata::ehid);
    assert(it != eq.edata.end());

    vec visit = std::vector(hm.nH, false);
    vec allow(hm.nV, false);

    for (auto [ehid_, _]: eq.edata) {
    for (Half h: ehalfs[ehid_].halfs) {
        visit[h.id] = true;
        visit[h.twin().id] = true; // 双対も念のため
    }}

    for (int vid: eq.verts_inside()) { allow[vid] = true; }
    auto it_prev = circular_prev(eq.edata, it);
    assert(it->side != it_prev->side);

    Ehalf& eh_prev = ehalfs[it_prev->ehid];
    Vert v0 = eh_prev.halfs.front().tail();
    Vert v1 = eh.halfs.back().head();

    // guarantees bgn/end vertex id is in the allowed list
    allow[v0.id] = true;
    allow[v1.id] = true;
    return compute_dijkstra_for_tquad_temp(hm, visit, allow, v0, v1);
}

inline bool Emesh::collapse_half(const int ehid) {
    auto& eh = ehalfs[ehid];
    auto& eq = equads[eh.eqid];
    auto it0 = rg::find(eq.edata, ehid, &Edata::ehid); // edata of ehid
    auto it1 = circular_prev(eq.edata, it0);           // edata of prev of ehid

    auto path0 = collapse_half_find_path(ehid);
    auto path1 = path0 | vw::reverse | vw::transform(&Half::twin) | rg::to<vec<Half>>();

    Ehalf& eh_prev_curr = ehalfs[it1->ehid];
    Ehalf& eh_prev_twin = ehalfs[eh_prev_curr.twid];

    // 1: erase ehalf of this tquad
    // 2: Replace prev (and twin of prev) ehalf
    eq.edata.erase(it0);
    eh_prev_curr.halfs = path0;
    eh_prev_twin.halfs = path1;

    // 3: extend ehalf of twin tquad
    Equad& eq_twin = equads[eh_prev_twin.eqid];
    auto it_twin = rg::find(eq_twin.edata, eh_prev_twin.id, &Edata::ehid);
    if (it_twin != eq_twin.edata.end()) {
        Ehalf& eh_twin = ehalfs[circular_prev(eq_twin.edata, it_twin)->ehid];
        eh_twin.extend_next(eh);
    }

    // 4: twin of eh is consumed to its tquad
    if (eh.twid >= 0) {
        Ehalf& eh_tw = ehalfs[eh.twid];
        Equad& eq_eh_tw = equads[eh_tw.eqid];
        auto it_eh_tw = rg::find(eq_eh_tw.edata, eh.twid, &Edata::ehid);
        if (it_eh_tw != eq_eh_tw.edata.end()) {
            Ehalf& next_eh = ehalfs[circular_next(eq_eh_tw.edata, it_eh_tw)->ehid];
            next_eh.extend_prev(eh_tw);
            eq_eh_tw.edata.erase(it_eh_tw);
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

        auto c = polyscope::registerCurveNetwork("collapse-half "+ std::to_string(ehid), ns, es);
        c->resetTransform();
        c->setRadius(0.002);
    }

    return true;
}

}

#endif //TMESH_H_EMESH_COLLAPSE_EHALF_H
