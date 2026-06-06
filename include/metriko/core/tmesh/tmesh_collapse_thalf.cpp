//
// Created by saki on 2026/06/06.
//

#include "./tmesh.h"

namespace metriko {

constexpr auto circular_prev = [](auto& c, auto it) { return it == c.begin() ? std::prev(c.end()) : std::prev(it); };
constexpr auto circular_next = [](auto& c, auto it) { auto n = std::next(it); return n == c.end() ? c.begin() : n; };


bool Tmesh::collapse_ehalf(int thid) {
    Thalf& th = thalfs[thid];
    Tquad& tq = tquads[th2quad[th.id]];
    auto it0 = rg::find(tq.data, thid, &Tdata::thid); // tdata of thid
    auto it_prev = circular_prev(tq.data, it0);
    auto it_next = circular_next(tq.data, it0);

    Thalf& th_prev = thalfs[it_prev->thid];
    Thalf& th_next = thalfs[it_next->thid];

    assert(!(it0->side != it_prev->side && it0->side != it_next->side));
    assert(!(it0->side == it_prev->side && it0->side == it_next->side));
    std::optional<vec<Half>> path_res = std::nullopt;
    bool merge_to_prev = it0->side != it_prev->side;

    Vert v0, v1;

    //if (merge_to_prev) {
    //    v0 = th_prev.halfs.front().tail();
    //    v1 = th.halfs.back().head();
    //} else {
    //    v0 = th.halfs.front().tail();
    //    v1 = th_next.halfs.back().head();
    //}
    //path_res = compute_dijkstra_for_tquad_temp(hm, allow, v0, v1);

    if (!path_res.has_value()) { std::cout << "failed to find path_res" << std::endl; return false; }

    auto path0 = path_res.value();
    auto path1 = path0 | vw::reverse | vw::transform(&Half::twin) | rg::to<vec<Half>>();

    /*
    if (merge_to_prev) {
        Thalf& th_prev_curr = thalfs[it_prev->thid];
        Thalf& th_prev_twin = thalfs[th_prev_curr.twid];

        // 1: erase ehalf of this tquad
        // 2: Replace prev (and twin of prev) ehalf
        tq.data.erase(it0);
        th2quad[th.id] = -1;
        th_prev_curr.halfs = path0;
        th_prev_twin.halfs = path1;

        // 3: extend ehalf of twin tquad
        Tquad& tq_twin = tquads[th_prev_twin.tqid];
        auto it_twin = rg::find(tq_twin.data, th_prev_twin.id, &Tdata::thid);
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
        Thalf& eh_next_curr = thalfs[it_next->thid];
        Thalf& eh_next_twin = thalfs[eh_next_curr.twid];

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
    */
    return true;
}
}