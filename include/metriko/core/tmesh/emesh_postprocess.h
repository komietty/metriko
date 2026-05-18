#ifndef EXAMPLE_EBD_CPP_EMESH_POSTPROCESS_H
#define EXAMPLE_EBD_CPP_EMESH_POSTPROCESS_H
#include "emesh.h"

namespace metriko::tutte {

struct HalfData {
    Half half;
    double v0;
    double v1;
    int thid;
    int tqid;
    int twin;
    int order; // order inside thalf
    bool first = false;
    bool crash = false;

    bool operator<(const HalfData& rhs) const noexcept {
        return std::tie(tqid, thid, order) < std::tie(rhs.tqid, rhs.thid, rhs.order);
    }
};

inline vec<HalfData> compute_half_data(const Emesh& em) {
    vec<HalfData> half_data;
    for (auto& eq: em.equads) {
        if (eq.id == -1) continue;
        for (int i = 0; i < 4; i++) {
            for (int ehid: eq.ehids(i)) {
                double sum = 0;
                int order = 0;
                auto& eh = em.ehalfs[ehid];
                double l = 0;
                for (Half h: eh.halfs) l += h.len();
                double d = 1. / l;
                for (Half h: eh.halfs) {
                    auto hl = h.len();
                    auto hd = HalfData { h, sum, sum + d * hl, ehid, eq.id, -1, order, false, false };
                    half_data.push_back(hd);
                    sum += d * hl;
                    order++;
                }
            }
        }
    }

    // --- ★ 追加: equal_range を機能させるために事前ソート ---
    // tqid -> thid -> order の優先順位で並び替える
    rg::sort(half_data.begin(), half_data.end(), [](const HalfData& a, const HalfData& b) {
        if (a.tqid != b.tqid) return a.tqid < b.tqid;
        if (a.thid != b.thid) return a.thid < b.thid;
        return a.order < b.order;
    });

    // --- 2パス目: twin のインデックスを解決する ---

    // ハーフエッジのIDをキーにして、half_data 配列内のインデックスを引けるマップを作成
    std::unordered_map<int, int> half_id_to_idx;
    for (int i = 0; i < half_data.size(); i++)
        half_id_to_idx[half_data[i].half.id] = i;

    // assign twin for each halfedge
    for (HalfData& d : half_data) {
        auto it = half_id_to_idx.find(d.half.twin().id);
        if (it != half_id_to_idx.end())  d.twin = it->second;
        else throw std::runtime_error("twin half does not exist");
    }

    return half_data;
}
}

#endif //EXAMPLE_EBD_CPP_EMESH_POSTPROCESS_H
