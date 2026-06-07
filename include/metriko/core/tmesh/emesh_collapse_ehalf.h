#ifndef EXAMPLE_EBD_CPP_EMESH_COLLAPSE_EHALF_H
#define EXAMPLE_EBD_CPP_EMESH_COLLAPSE_EHALF_H
#include "./emesh.h"
#include "./emesh_collapse_util.h"
#include "./emesh_collapse_dijkstra.h"

namespace metriko {

inline bool Emesh::collapse_ehalf(const int ehid) {
    auto& eh = ehalfs[ehid];
    auto& eq = equads[eh.eqid];
    auto it0 = rg::find(eq.data, ehid, &Edata::ehid); // edata of ehid
    auto it_prev = circular_prev(eq.data, it0);
    auto it_next = circular_next(eq.data, it0);

    // try to find a path
    vec allow(hm.nV, false);

    for (int vid: verts_inside(eq, *this, hm)) { allow[vid] = true; }
    Ehalf& eh_prev = ehalfs[it_prev->ehid];
    Ehalf& eh_next = ehalfs[it_next->ehid];

    assert(!(it0->side != it_prev->side && it0->side != it_next->side));
    assert(!(it0->side == it_prev->side && it0->side == it_next->side));
    std::optional<vec<Half>> path_res = std::nullopt;
    bool merge_to_prev = it0->side != it_prev->side;

    if (merge_to_prev) {
        Vert v0 = eh_prev.halfs.front().tail();
        Vert v1 = eh.halfs.back().head();

        for (Half& h : eh.halfs) {
            allow[h.tail().id] = false;
            allow[h.head().id] = false;
        }

        allow[v0.id] = true;
        allow[v1.id] = true;
        path_res = compute_dijkstra_for_tquad_temp(hm, allow, v0, v1);
    } else {
        Vert v0 = eh.halfs.front().tail();
        Vert v1 = eh_next.halfs.back().head();

        for (Half& h : eh.halfs) {
            allow[h.tail().id] = false;
            allow[h.head().id] = false;
        }

        allow[v0.id] = true;
        allow[v1.id] = true;
        path_res = compute_dijkstra_for_tquad_temp(hm, allow, v0, v1);
    }

    if (!path_res.has_value()) { std::cout << "failed to find path_res" << std::endl; return false; }

    auto path0 = path_res.value();
    auto path1 = path0 | vw::reverse | vw::transform(&Half::twin) | rg::to<vec<Half>>();

    if (merge_to_prev) {
        Ehalf& eh_prev_curr = ehalfs[it_prev->ehid];
        Ehalf& eh_prev_twin = ehalfs[eh_prev_curr.twid];

        // 1: erase ehalf of this tquad
        // 2: Replace prev (and twin of prev) ehalf
        eq.data.erase(it0);
        eh.eqid = -1;
        eh_prev_curr.halfs = path0;
        eh_prev_twin.halfs = path1;

        // 3: extend ehalf of twin tquad
        Equad& eq_twin = equads[eh_prev_twin.eqid];
        auto it_twin = rg::find(eq_twin.data, eh_prev_twin.id, &Edata::ehid);
        if (it_twin != eq_twin.data.end()) {
            Ehalf& eh_twin = ehalfs[circular_prev(eq_twin.data, it_twin)->ehid];
            eh_twin.extend_next(eh);
        }

        // 4: twin of eh is consumed to its tquad
        if (eh.twid >= 0) {
            Ehalf& eh_tw = ehalfs[eh.twid];
            Equad& eq_eh_tw = equads[eh_tw.eqid];
            auto it_eh_tw = rg::find(eq_eh_tw.data, eh.twid, &Edata::ehid);
            if (it_eh_tw != eq_eh_tw.data.end()) {
                Ehalf& next_eh = ehalfs[circular_next(eq_eh_tw.data, it_eh_tw)->ehid];
                next_eh.extend_prev(eh_tw);
                eq_eh_tw.data.erase(it_eh_tw);
                eh_tw.eqid = -1;
            }
        }
    } else {
        Ehalf& eh_next_curr = ehalfs[it_next->ehid];
        Ehalf& eh_next_twin = ehalfs[eh_next_curr.twid];

        // 1: erase ehalf of this tquad
        // 2: Replace next (and twin of next) ehalf
        eq.data.erase(it0);
        eh.eqid = -1;
        eh_next_curr.halfs = path0;
        eh_next_twin.halfs = path1;

        // 3: extend ehalf of twin tquad
        Equad& eq_twin = equads[eh_next_twin.eqid];
        auto it_twin = rg::find(eq_twin.data, eh_next_twin.id, &Edata::ehid);
        if (it_twin != eq_twin.data.end()) {
            Ehalf& eh_twin = ehalfs[circular_next(eq_twin.data, it_twin)->ehid];
            eh_twin.extend_prev(eh);
        }

        // 4: twin of eh is consumed to its tquad
        if (eh.twid >= 0) {
            auto& eh_tw = ehalfs[eh.twid];
            auto& eq_eh_tw = equads[eh_tw.eqid];
            auto  it_eh_tw = rg::find(eq_eh_tw.data, eh.twid, &Edata::ehid);

            if (it_eh_tw != eq_eh_tw.data.end()) {
                Ehalf& prev_eh = ehalfs[circular_prev(eq_eh_tw.data, it_eh_tw)->ehid];
                prev_eh.extend_next(eh_tw);
                eq_eh_tw.data.erase(it_eh_tw);
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
        c->setEnabled(false);
        c->resetTransform();
        c->setRadius(0.002);
    }

    return true;
}
}

#endif
