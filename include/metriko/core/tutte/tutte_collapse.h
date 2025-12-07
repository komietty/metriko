//
// Created by saki on 2025/12/05.
//

#ifndef TMESH_H_TUTTE_COLLAPSE_H
#define TMESH_H_TUTTE_COLLAPSE_H
#include "compute_dijkstra.h"

namespace metriko::tutte {

struct AuxDijkData {
    Vert vert;
    double value;

    bool operator<(const AuxDijkData& rhs) const noexcept {
        return value < rhs.value;
    }
};

inline void extract_polyline_from_tquad(
    Hmesh& hm,
    Tmesh& tm,
    Tquad& tq,
    set<HalfData>& data,
    VecXd& R,
    VecXd& X
) {

    int side = -1; // if 0 or 1, collapse
    for (int i = 0; i < 2; i++) {
        int sum = rg::fold_left(tq.thids_by_side(i), 0,
                                [&](int acc, int thid) { return acc + (int)X[tm.thalfs[thid].edge().id]; });
        if (sum == 1) side = i; // todo: right now not 0 but 1
    }

    if (side == -1) return;

    int bgn = -1;
    int end = -1;
    set<AuxDijkData> aux;

    auto thids_p = tq.thids_by_side(side == 0 ? 0 : 1); // collapse side p
    auto thids_q = tq.thids_by_side(side == 0 ? 2 : 3); // collapse side q
    auto thids_a = tq.thids_by_side(side == 0 ? 1 : 2); // remain side a
    auto thids_b = tq.thids_by_side(side == 0 ? 3 : 0); // remain side b

    auto hdata_tq = rg::equal_range(data, tq.id, {}, &HalfData::tqid);
    auto hdata_th_p = hdata_tq | vw::filter([&](const HalfData& hd) { return rg::contains(thids_p, hd.thid); });
    auto hdata_th_q = hdata_tq | vw::filter([&](const HalfData& hd) { return rg::contains(thids_q, hd.thid); });

    for (auto hd: hdata_th_p) { if (hd.first) { bgn = hd.half.tail().id; break; } }
    if (bgn == -1) { for (auto hd: hdata_th_p) { if (hd.crash) { bgn = hd.half.tail().id; break; } } }

    double sum = 0;

    for (int thid_: thids_a) {
        const Tedge& te = tm.tedges[tm.thalfs[thid_].teid];
        sum += X[te.id];

        auto hds = hdata_tq | vw::filter([&](const HalfData& hd) { return hd.thid == thid_; });
        auto last = rg::rbegin(hds);
        aux.emplace(AuxDijkData{ last->half.head(), sum });
    }

    for (auto hd: hdata_th_q) { if (hd.first) { end = hd.half.tail().id; break; } }
    if (end == -1) { for (auto hd: hdata_th_q) { if (hd.crash) { end = hd.half.tail().id; break; } } }

    aux.emplace(AuxDijkData{hm.verts[bgn], 0});
    aux.emplace(AuxDijkData{hm.verts[end], sum});

    for (int thid_: thids_b) {
        const Tedge& te = tm.tedges[tm.thalfs[thid_].teid];
        sum -= X[te.id];

        auto hds = hdata_tq | vw::filter([&](const HalfData& hd) { return hd.thid == thid_; });
        auto last = rg::rbegin(hds);
        aux.emplace(AuxDijkData{ last->half.head(), sum });
    }


    assert(sum == 0);
    std::vector<int> full_path;
    auto visit = std::vector(hm.nH, false);

    std::vector aux_sorted(aux.begin(), aux.end());

    {
        Vert v0 = aux_sorted[0].vert;
        Vert v1 = aux_sorted[1].vert;
        auto path = compute_dijkstra(hm, visit, v0, v1);
        for (int hid: path) { visit[hid] = true; }
        full_path.insert(full_path.end(), path.begin(), path.end());
    }
    {
        Vert v0 = aux_sorted[2].vert;
        Vert v1 = aux_sorted[3].vert;
        auto path = compute_dijkstra(hm, visit, v0, v1);
        for (int hid: path) { visit[hid] = true; }
        full_path.insert(full_path.end(), path.begin(), path.end());
    }
    {
        Vert v0 = aux_sorted[1].vert;
        Vert v1 = aux_sorted[2].vert;
        auto path = compute_dijkstra(hm, visit, v0, v1);
        for (int hid: path) { visit[hid] = true; }
        full_path.insert(full_path.end(), path.begin(), path.end());
    }


    // debug polyline
    {
        std::vector<glm::vec3> pts;
        std::vector<double> val;

        for (auto& ad: aux) {
            Vert v = ad.vert;
            pts.emplace_back(v.pos().x(), v.pos().y(), v.pos().z());
            val.emplace_back(ad.value);
        }

        auto p = polyscope::registerPointCloud("tquad polyline verts", pts);
        p->addScalarQuantity("value", val);
        p->setEnabled(true);
        p->setPointRadius(0.005);
        p->resetTransform();
    }
    {

        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> val1; //
        size_t count = 0;

        for (const int hid: full_path) {
            Half h = hm.halfs[hid];
            Row3d p1 = h.tail().pos();
            Row3d p2 = h.head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{count, count + 1});
            count += 2;
        }
        auto c = polyscope::registerCurveNetwork("dijkstra", ns, es);
        //auto v1 = c->addNodeScalarQuantity("val1", val1);
        c->setEnabled(false);
        c->resetTransform();
        c->setRadius(0.002);
    }
}
}

#endif //TMESH_H_TUTTE_COLLAPSE_H
