//
//--- Copyright (C) 2026 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
#ifndef METRIKO_SUBDIVIDE_TRACKER_H
#define METRIKO_SUBDIVIDE_TRACKER_H

#include "metriko/core/hmesh/hmesh.h"
#include <vector>
#include <complex>
#include <map>

#include "metriko/core/hmesh/utilities.h"

namespace metriko {

// 細分化されたメッシュと、その履歴（トラッキング情報）を保持する構造体
struct TrackedDenseMesh {
    vec<vec<int>> polygons; // 新しいメッシュを構築するためのFaceごとの頂点IDリスト (CCW順)
    vec<vec<complex>> uvs;  // 各Faceの各コーナーのUV座標（親FaceのローカルUV空間における座標）
    vec<int> face2parent;   // 新しいFace ID -> 元の(Level 0)親Face ID. Dijkstra探索時に「今元のメッシュのどこにいるか」を知るための最強の道標
    std::map<std::pair<int, int>, int> edge_to_old_id;
    int num_verts = -1;     // 新しいメッシュの総頂点数
};

inline TrackedDenseMesh compute_midpoint_subdivision(
    const Hmesh& base_hm,
    const VecXc& base_cf,
    const int levels
) {
    TrackedDenseMesh curr;
    curr.num_verts = base_hm.verts.size();

    // ==========================================
    // 1. Level 0 (元のメッシュ) の情報を抽出
    // ==========================================
    for (const auto& f : base_hm.faces) {
        vec<int> poly;
        vec<complex> f_uvs;
        for (auto h : f.adjHalfs()) {
            poly.push_back(h.crnr().vert().id);
            f_uvs.push_back(base_cf[h.crnr().id]);
        }
        curr.polygons.push_back(poly);
        curr.uvs.push_back(f_uvs);
        curr.face2parent.push_back(f.id); // 初期は自分自身が親
    }

    // エッジの重複分割を防ぐための無向エッジキー
    auto get_key = [](int v1, int v2) {
        return std::pair<int, int>(std::min(v1, v2), std::max(v1, v2));
    };

    for (int i = 0; i < base_hm.nE; ++i) {
        auto e = base_hm.edges[i];
        curr.edge_to_old_id[get_key(e.vert0().id, e.vert1().id)] = i;
    }

    // ==========================================
    // 2. 指定レベルだけ 1-to-4 の再帰的細分化
    // ==========================================
    for (int lvl = 0; lvl < levels; ++lvl) {
        TrackedDenseMesh next;
        next.num_verts = curr.num_verts;
        std::map<std::pair<int, int>, int> edge2mid;

        for (int i = 0; i < curr.polygons.size(); ++i) {
            const auto& poly = curr.polygons[i];
            const auto& f_uvs = curr.uvs[i];
            int parent_fid = curr.face2parent[i]; // 親の血統を引き継ぐ

            assert(poly.size() == 3 && "Only triangle meshes are supported for midpoint subdivision.");

            int     v0 = poly[0],  v1 = poly[1],  v2 = poly[2];
            complex u0 = f_uvs[0], u1 = f_uvs[1], u2 = f_uvs[2];

            // 各エッジの中点頂点IDを取得または新規作成
            int m01, m12, m20;
            auto k01 = get_key(v0, v1);
            if (edge2mid.contains(k01)) m01 = edge2mid[k01];
            else { m01 = next.num_verts++; edge2mid[k01] = m01; }

            auto k12 = get_key(v1, v2);
            if (edge2mid.contains(k12)) m12 = edge2mid[k12];
            else { m12 = next.num_verts++; edge2mid[k12] = m12; }

            auto k20 = get_key(v2, v0);
            if (edge2mid.contains(k20)) m20 = edge2mid[k20];
            else { m20 = next.num_verts++; edge2mid[k20] = m20; }

            // 中点のUVを親FaceのローカルUV空間で線形補間
            complex mu01 = (u0 + u1) * 0.5;
            complex mu12 = (u1 + u2) * 0.5;
            complex mu20 = (u2 + u0) * 0.5;

            // 4つの新しいFaceをCCW順で追加し、すべて同じ親FaceIDを割り当てる
            // Face 1: 角0
            next.polygons.push_back({v0, m01, m20});
            next.uvs.push_back({u0, mu01, mu20});
            next.face2parent.push_back(parent_fid);

            // Face 2: 角1
            next.polygons.push_back({v1, m12, m01});
            next.uvs.push_back({u1, mu12, mu01});
            next.face2parent.push_back(parent_fid);

            // Face 3: 角2
            next.polygons.push_back({v2, m20, m12});
            next.uvs.push_back({u2, mu20, mu12});
            next.face2parent.push_back(parent_fid);

            // Face 4: 中央
            next.polygons.push_back({m01, m12, m20});
            next.uvs.push_back({mu01, mu12, mu20});
            next.face2parent.push_back(parent_fid);

            // 【新設計】subdivide.h のように、親エッジが存在する場合のみ、分割後の2本のエッジにIDを引き継ぐ
            if (curr.edge_to_old_id.contains(k01)) {
                int old_id = curr.edge_to_old_id[k01];
                next.edge_to_old_id[get_key(v0, m01)] = old_id;
                next.edge_to_old_id[get_key(m01, v1)] = old_id;
            }
            if (curr.edge_to_old_id.contains(k12)) {
                int old_id = curr.edge_to_old_id[k12];
                next.edge_to_old_id[get_key(v1, m12)] = old_id;
                next.edge_to_old_id[get_key(m12, v2)] = old_id;
            }
            if (curr.edge_to_old_id.contains(k20)) {
                int old_id = curr.edge_to_old_id[k20];
                next.edge_to_old_id[get_key(v2, m20)] = old_id;
                next.edge_to_old_id[get_key(m20, v0)] = old_id;
            }

        }
        curr = std::move(next);
    }

    return curr;
}

inline std::vector<Row3d> reconstruct_3d_positions(
    const TrackedDenseMesh& dmesh,
    const Hmesh& base_hm,
    const VecXc& base_uv
) {
    std::vector<Row3d> dense_pos(dmesh.num_verts);
    std::vector<bool> visited(dmesh.num_verts, false);

    for (int i = 0; i < dmesh.polygons.size(); ++i) {
        const auto& poly = dmesh.polygons[i];
        const auto& f_uvs = dmesh.uvs[i];
        int parent_fid = dmesh.face2parent[i];

        for (int j = 0; j < poly.size(); ++j) {
            int vid = poly[j];

            if (!visited[vid]) {
                dense_pos[vid] = conversion_2d_3d(base_hm.faces[parent_fid], base_uv, f_uvs[j]);
                visited[vid] = true;
            }
        }
    }
    return dense_pos;
}

// subdivide.h と全く同じロジックで seam を写像する
inline std::vector<bool> compute_dense_seam(
    const std::vector<bool>& base_seam,
    const Hmesh& dense_hm,
    const TrackedDenseMesh& sdiv_data
) {
    std::vector<bool> dense_seam(dense_hm.nE, false);
    auto make_edge_key = [](int a, int b) {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };

    for (int i = 0; i < dense_hm.nE; ++i) {
        auto e = dense_hm.edges[i];
        auto key = make_edge_key(e.vert0().id, e.vert1().id);

        if (sdiv_data.edge_to_old_id.find(key) != sdiv_data.edge_to_old_id.end()) {
            int old_id = sdiv_data.edge_to_old_id.at(key);
            if (!base_seam.empty() && base_seam[old_id]) {
                dense_seam[i] = true;
            }
        }
    }
    return dense_seam;
}

}
#endif

