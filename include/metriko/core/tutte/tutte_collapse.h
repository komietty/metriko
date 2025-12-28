//
// Created by saki on 2025/12/05.
//

#ifndef TMESH_H_TUTTE_COLLAPSE_H
#define TMESH_H_TUTTE_COLLAPSE_H
#include "compute_dijkstra.h"

namespace metriko::tutte {

struct AuxDijkData {
    Vert vert;
    int side;
    double value;

    bool operator<(const AuxDijkData& rhs) const noexcept {
        return value < rhs.value;
    }
};

inline void extract_polyline_from_tquad(
    Hmesh& hm,
    tm::Tmesh& tm,
    tm::Tquad& tq,
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
    int bgn_ord = -1;
    int end_ord = -1;
    set<AuxDijkData> aux;

    int remain_side_a = side == 0 ? 1 : 2;
    int remain_side_b = side == 0 ? 3 : 0;
    auto thids_p = tq.thids_by_side(side == 0 ? 0 : 1); // collapse side p
    auto thids_q = tq.thids_by_side(side == 0 ? 2 : 3); // collapse side q
    auto thids_a = tq.thids_by_side(remain_side_a); // remain side a
    auto thids_b = tq.thids_by_side(remain_side_b); // remain side b

    auto hdata_tq = rg::equal_range(data, tq.id, {}, &HalfData::tqid);
    auto hdata_th_p = hdata_tq | vw::filter([&](const HalfData& hd) { return rg::contains(thids_p, hd.thid); });
    auto hdata_th_q = hdata_tq | vw::filter([&](const HalfData& hd) { return rg::contains(thids_q, hd.thid); });

    for (auto hd: hdata_th_p) {
        std::cout << "hid in p: " << hd.half.tail().pos() << std::endl;
        std::cout << "hid in p: " << hd.half.head().pos() << std::endl;
        if (hd.first) { bgn = hd.half.tail().id; bgn_ord = hd.order; break; }
        Half h = hd.half.twin();
        for (const HalfData& hd_ : data) {
            if (hd_.half == h && hd_.first) {
                bgn = hd_.half.tail().id;
                bgn_ord = 1; // todo: just wanted to represent not 0
                break;
            }
        }
    }
    if (bgn == -1) {
        for (auto hd: hdata_th_p) {
            if (hd.crash) {
                bgn = hd.half.tail().id; // todo: head??
                bgn_ord = hd.order;
                break;
            }
        }
    }

    // if not find
    if (bgn == -1) {
        for (auto hd: hdata_th_p) {
            Half h = hd.half.twin();
            for (const HalfData& hd_ : data) {
                if (hd_.half == h && hd_.crash) {
                    bgn = hd_.half.head().id;
                    bgn_ord = hd.order;
                    break;
                }
            }
        }
    }

    double sum = 0;

    for (int thid_: thids_a) {
        const auto& th = tm.thalfs[thid_];
        const auto& te = tm.tedges[th.teid];
        sum += X[te.id];

        auto hds = hdata_tq | vw::filter([&](const HalfData& hd) { return hd.thid == thid_; });
        auto last = rg::rbegin(hds);
        aux.emplace(AuxDijkData{ last->half.head(), tm.th2side(thid_), sum });
    }

    for (auto hd: hdata_th_q) {
        if (hd.first) { end = hd.half.tail().id; break; }

        Half h = hd.half.twin();
        for (const HalfData& hd_ : data) {
            if (hd_.half == h && hd_.first) {
                end = hd_.half.tail().id;
                end_ord = 1; // todo: just wanted to represent not 0
                break;
            }
        }
    }

    if (end == -1) {
        for (auto hd: hdata_th_q) {
            if (hd.crash) {
                end = hd.half.tail().id;
                break;
            }
        }
    }

    // if not find
    if (end == -1) {
        for (auto hd: hdata_th_q) {
            Half h = hd.half.twin();
            for (const HalfData& hd_ : data) {
                if (hd_.half == h && hd_.crash) {
                    end = hd_.half.head().id;
                    end_ord = hd.order;
                    break;
                }
            }
        }
    }

    int side_bgn = bgn_ord == 0 ? remain_side_b : remain_side_a;
    int side_end = end_ord == 0 ? remain_side_a : remain_side_b;
    aux.emplace(AuxDijkData{hm.verts[bgn], side_bgn, 0.});
    aux.emplace(AuxDijkData{hm.verts[end], side_end, sum});

    for (int thid_: thids_b) {
        const auto& th = tm.thalfs[thid_];
        const auto& te = tm.tedges[th.teid];
        sum -= X[te.id];

        auto hds = hdata_tq | vw::filter([&](const HalfData& hd) { return hd.thid == thid_; });
        auto last = rg::rbegin(hds);
        aux.emplace(AuxDijkData{ last->half.head(), tm.th2side(thid_), sum });
    }

    if (aux.size() == 0) return;


    assert(sum == 0);
    //assert(bgn != -1);
    //assert(end != -1);
    std::cout << "bgn: " << bgn << std::endl;
    std::cout << "end: " << end << std::endl;
    std::vector<int> full_path;
    auto visit = std::vector(hm.nH, false);

    std::vector aux_sorted(aux.begin(), aux.end());

    for (int i = 0; i < aux_sorted.size() - 1; i++) {
        auto& a0 = aux_sorted[i];
        auto& a1 = aux_sorted[i + 1];
        if (a0.side == a1.side) {
            auto path = compute_dijkstra(hm, visit, a0.vert, a1.vert);
            for (int hid: path) { visit[hid] = true; }
            full_path.insert(full_path.end(), path.begin(), path.end());
        }
    }

    for (int i = 0; i < aux_sorted.size() - 1; i++) {
        auto& a0 = aux_sorted[i];
        auto& a1 = aux_sorted[i + 1];
        if (a0.side != a1.side) {
            auto path = compute_dijkstra(hm, visit, a0.vert, a1.vert);
            for (int hid: path) { visit[hid] = true; }
            full_path.insert(full_path.end(), path.begin(), path.end());
        }
    }
    /*
    */

    // debug polyline
    {
        std::vector<glm::vec3> pts;
        std::vector<double> val;
        std::vector<double> side;

        std::cout << "aux size: " << aux.size() << std::endl;

        for (auto& ad: aux) {
            Vert v = ad.vert;
            pts.emplace_back(v.pos().x(), v.pos().y(), v.pos().z());
            val.emplace_back(ad.value);
            side.emplace_back(ad.side);
        }

        auto p = polyscope::registerPointCloud("tquad polyline verts " + std::to_string(tq.id), pts);
        p->addScalarQuantity("value", val);
        p->addScalarQuantity("side", side);
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
        auto c = polyscope::registerCurveNetwork("dijkstra " + std::to_string(tq.id), ns, es);
        //auto v1 = c->addNodeScalarQuantity("val1", val1);
        c->resetTransform();
        c->setRadius(0.002);
    }
    /*
     */
}
}

#endif //TMESH_H_TUTTE_COLLAPSE_H
