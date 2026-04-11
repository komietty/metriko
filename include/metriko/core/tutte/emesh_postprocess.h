//
// Created by saki on 2026/04/10.
//

#ifndef TMESH_H_EMESH_POSTPROCESS_H
#define TMESH_H_EMESH_POSTPROCESS_H
#include "./emesh.h"
#include <unordered_map>

namespace metriko::tutte {

vec<HalfData> compute_half_data(const Emesh& em) {
    vec<HalfData> half_data;
    for (auto& eq: em.equads) {
        for (int i = 0; i < 4; i++) {
            for (int ehid: eq.ehids_by_side(i)) {
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
    for (int i = 0; i < half_data.size(); i++) {
        half_id_to_idx[half_data[i].half.id] = i;
    }

    // 各 HalfData について、自分の twin となる Half の ID をマップで検索し、インデックスをセット
    for (int i = 0; i < half_data.size(); i++) {
        int twin_id = half_data[i].half.twin().id; // twin の ID を取得
        auto it = half_id_to_idx.find(twin_id);

        if (it != half_id_to_idx.end()) {
            half_data[i].twin = it->second; // 見つかったらそのインデックスを代入
        } else {
            half_data[i].twin = -1; // 境界エッジなどで twin がリストに存在しない場合は -1 のまま
        }
    }

    return half_data;
}
}

#endif
