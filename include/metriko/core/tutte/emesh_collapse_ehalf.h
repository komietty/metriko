//
// Created by saki on 2026/04/09.
//

#ifndef TMESH_H_EMESH_COLLAPSE_EHALF_H
#define TMESH_H_EMESH_COLLAPSE_EHALF_H

#include "./emesh.h"
namespace metriko::tutte {

constexpr auto circular_prev = [](auto& c, auto it) { return it == c.begin() ? std::prev(c.end()) : std::prev(it); };
constexpr auto circular_next = [](auto& c, auto it) { auto n = std::next(it); return n == c.end() ? c.begin() : n; };

/*
inline std::optional<vec<Half>> Emesh::collapse_half_find_path(int ehid) {
    auto& eh = ehalfs[ehid];
    auto& eq = equads[eh.eqid]; // todo: inconsistent (maybe after collapse)!!!
    std::cout << "eqid: " << eq.id << std::endl;
    for (auto [ehid, side]: eq.edata) {
        std::cout << "ehid: " << ehid << ", side: " << side << std::endl;
    }

    auto  it = rg::find(eq.edata, ehid, &Edata::ehid);
    if(it == eq.edata.end()) {
        std::cout << "collapse_half_find_path... ehid: " << eh.id << " eh.eqid: " << eh.eqid << ", eq.id: " << eq.id << std::endl;
        eq.debug_draw();
        return std::nullopt;
    }

    vec allow(hm.nV, false);

    for (int vid: eq.verts_inside()) { allow[vid] = true; }
    auto it_prev = circular_prev(eq.edata, it);
    auto it_next = circular_next(eq.edata, it);
    Ehalf& eh_prev = ehalfs[it_prev->ehid];
    Ehalf& eh_next = ehalfs[it_next->ehid];

    assert(!(it->side != it_prev->side && it->side != it_next->side));
    assert(!(it->side == it_prev->side && it->side == it_next->side));

    if (it->side != it_prev->side) {
        Vert v0 = eh_prev.halfs.front().tail();
        Vert v1 = eh.halfs.back().head();
        allow[v0.id] = true;
        allow[v1.id] = true;
        return compute_dijkstra_for_tquad_temp(hm, allow, v0, v1);
    }

    if (it->side != it_next->side) {
        Vert v0 = eh.halfs.front().tail();
        Vert v1 = eh_next.halfs.back().head();
        allow[v0.id] = true;
        allow[v1.id] = true;
        return compute_dijkstra_for_tquad_temp(hm, allow, v0, v1);
    }

    return std::nullopt;
}
*/

inline bool Emesh::collapse_half(const int ehid) {
    auto& eh = ehalfs[ehid];
    auto& eq = equads[eh.eqid];
    auto it0 = rg::find(eq.edata, ehid, &Edata::ehid); // edata of ehid
    auto it_prev = circular_prev(eq.edata, it0);
    auto it_next = circular_next(eq.edata, it0);

    // try to find a path
    vec allow(hm.nV, false);

    for (int vid: eq.verts_inside()) { allow[vid] = true; }
    Ehalf& eh_prev = ehalfs[it_prev->ehid];
    Ehalf& eh_next = ehalfs[it_next->ehid];

    if(it0->side == it_prev->side && it0->side == it_next->side) {
        //std::cout << "next ehid: " << it_next->ehid << ", side: " << it_next->side << std::endl;
        //std::cout << "prev ehid: " << it_prev->ehid << ", side: " << it_prev->side << std::endl;
        eh.debug_draw();
        eq.debug_draw();
        //for (auto [ehid, side]: eq.edata) {
        //    std::cout << "ehid: " << ehid << ", side: " << side << std::endl;
        //}
        return false;
    }

    assert(!(it0->side != it_prev->side && it0->side != it_next->side));
    assert(!(it0->side == it_prev->side && it0->side == it_next->side));
    std::optional<vec<Half>> path_res = std::nullopt;
    bool merge_to_prev = it0->side != it_prev->side;

    if (merge_to_prev) {
        Vert v0 = eh_prev.halfs.front().tail();
        Vert v1 = eh.halfs.back().head();
        allow[v0.id] = true;
        allow[v1.id] = true;
        path_res = compute_dijkstra_for_tquad_temp(hm, allow, v0, v1);
    } else {
        Vert v0 = eh.halfs.front().tail();
        Vert v1 = eh_next.halfs.back().head();
        allow[v0.id] = true;
        allow[v1.id] = true;
        path_res = compute_dijkstra_for_tquad_temp(hm, allow, v0, v1);
    }

    if (!path_res.has_value()) {
        std::cout << "failed to find path_res" << std::endl;
        return false;
    }

    auto path0 = path_res.value();
    auto path1 = path0 | vw::reverse | vw::transform(&Half::twin) | rg::to<vec<Half>>();

    if (merge_to_prev) {
        Ehalf& eh_prev_curr = ehalfs[it_prev->ehid];
        Ehalf& eh_prev_twin = ehalfs[eh_prev_curr.twid];

        // 1: erase ehalf of this tquad
        // 2: Replace prev (and twin of prev) ehalf
        eq.edata.erase(it0);
        eh.eqid = -1;
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
                eh_tw.eqid = -1;
            }
        }
    } else {
        Ehalf& eh_next_curr = ehalfs[it_next->ehid];
        Ehalf& eh_next_twin = ehalfs[eh_next_curr.twid];

        // 1: erase ehalf of this tquad
        // 2: Replace next (and twin of next) ehalf
        eq.edata.erase(it0);
        eh.eqid = -1;
        eh_next_curr.halfs = path0;
        eh_next_twin.halfs = path1;

        // 3: extend ehalf of twin tquad
        Equad& eq_twin = equads[eh_next_twin.eqid];
        auto it_twin = rg::find(eq_twin.edata, eh_next_twin.id, &Edata::ehid);
        if (it_twin != eq_twin.edata.end()) {
            Ehalf& eh_twin = ehalfs[circular_next(eq_twin.edata, it_twin)->ehid];
            eh_twin.extend_prev(eh);
        }

        // 4: twin of eh is consumed to its tquad
        if (eh.twid >= 0) {
            Ehalf& eh_tw = ehalfs[eh.twid];
            Equad& eq_eh_tw = equads[eh_tw.eqid];
            auto it_eh_tw = rg::find(eq_eh_tw.edata, eh.twid, &Edata::ehid);

            if (it_eh_tw != eq_eh_tw.edata.end()) {
                Ehalf& prev_eh = ehalfs[circular_prev(eq_eh_tw.edata, it_eh_tw)->ehid];
                prev_eh.extend_next(eh_tw);
                eq_eh_tw.edata.erase(it_eh_tw);
                eh_tw.eqid = -1;
            }
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
