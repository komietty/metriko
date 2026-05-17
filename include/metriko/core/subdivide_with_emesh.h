//
// Created by saki on 2026/04/15.
//

#ifndef TMESH_H_SUBDIVIDE_WITH_EMESH_H
#define TMESH_H_SUBDIVIDE_WITH_EMESH_H

#include <unordered_map>
#include <memory>
#include "./tutte_buk/emesh.h"

namespace metriko::tutte {

inline std::unique_ptr<Emesh> upgrade_emesh(
    const Emesh& old_em,
    const Hmesh& new_hm
) {
    auto new_em = std::make_unique<Emesh>(old_em, new_hm);

    // 2. 新しい Hmesh のハーフエッジを「tailとheadの頂点IDペア」で瞬時に引ける辞書を作る
    std::unordered_map<uint64_t, Half> v2h;
    v2h.reserve(new_hm.nH);
    for (Half h : new_hm.halfs) {
        uint64_t key = (uint64_t)h.tail().id << 32 | (uint32_t)h.head().id;
        v2h[key] = h;
    }

    // 3. 古いメッシュの頂点数（twelve_subdivideの中点IDの計算に使う）
    int orig_nV = old_em.hm.nV;

    // 4. Ehalf のパスを新しいHmeshのパス（長さ2倍）に置き換える
    for (auto& eh : new_em->ehalfs) {
        std::vector<Half> new_path;
        new_path.reserve(eh.halfs.size() * 2);

        for (Half old_h : eh.halfs) {
            int v0 = old_h.tail().id;
            int v1 = old_h.head().id;

            // twelve_subdivide における、外周エッジの中点IDの絶対法則
            int m = orig_nV + old_h.edge().id;

            // 新しい2つのハーフエッジ (v0 -> m) と (m -> v1) のキー
            uint64_t key1 = (uint64_t)v0 << 32 | (uint32_t)m;
            uint64_t key2 = (uint64_t)m << 32 | (uint32_t)v1;

            // 新しい経路に追加
            new_path.push_back(v2h.at(key1));
            new_path.push_back(v2h.at(key2));
        }

        eh.halfs = std::move(new_path);
    }

    return new_em;
}
}
#endif
