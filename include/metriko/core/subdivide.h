#ifndef TMESH_H_SUBDIVIDE_H
#define TMESH_H_SUBDIVIDE_H
#include "hmesh/hmesh.h"
#include <set>

namespace metriko {

struct Subdivide {
    std::unique_ptr<Hmesh> hm;
    VecXc uv;
    std::vector<bool> seam;
    VecXi matching; // 追加: エッジごとのperiod jump
    VecXi singular; // 追加: 頂点ごとの特異点チャージ
};

// 1-to-6 重心細分 (UV, Seam, Matching, Singular 対応版)
inline Subdivide barycentric_subdivide(
    const Hmesh& hm,
    const VecXc& uv,
    const std::vector<bool>& seam,
    const VecXi& matching_old = VecXi(), // 省略可能
    const VecXi& singular_old = VecXi(), // 省略可能
    int rosyN = 4 // クロスフィールド(4) または ラインフィールド(2)
) {
    int nV_new = hm.nV + hm.nE + hm.nF;
    int nF_new = hm.nF * 6;

    MatXd V_new(nV_new, 3);
    MatXi F_new(nF_new, 3);
    VecXc uv_new(nF_new * 3);

    int e_offset = hm.nV;
    int f_offset = hm.nV + hm.nE;

    // SeamとMatchingの引継ぎ用マップ: (min_v, max_v) -> old_edge_id
    std::map<std::pair<int, int>, int> edge_to_old_id;
    auto make_edge_key = [](int a, int b) {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };

    // 1. 既存の頂点のコピー
    for (int i = 0; i < hm.nV; ++i) {
        V_new.row(i) = hm.pos.row(i);
    }

    // 2. エッジ中点の座標計算と親子関係の記録
    for (int i = 0; i < hm.nE; ++i) {
        auto e = hm.edges[i];
        int v0 = e.vert0().id;
        int v1 = e.vert1().id;
        int m_id = e_offset + i;

        V_new.row(m_id) = (hm.pos.row(v0) + hm.pos.row(v1)) * 0.5;

        edge_to_old_id[make_edge_key(v0, m_id)] = i;
        edge_to_old_id[make_edge_key(m_id, v1)] = i;
    }

    // 3. 面の分割とコーナーUVの計算
    for (int i = 0; i < hm.nF; ++i) {
        auto f = hm.faces[i];
        auto h0 = f.half();
        auto h1 = h0.next();
        auto h2 = h1.next();

        int v0 = h0.tail().id;
        int v1 = h1.tail().id;
        int v2 = h2.tail().id;

        int m0 = e_offset + h0.edge().id;
        int m1 = e_offset + h1.edge().id;
        int m2 = e_offset + h2.edge().id;
        int c = f_offset + i;

        int c0 = i * 3 + 0;
        int c1 = i * 3 + 1;
        int c2 = i * 3 + 2;

        V_new.row(c) = (hm.pos.row(v0) + hm.pos.row(v1) + hm.pos.row(v2)) / 3.0;

        auto uv_v0 = uv(c0);
        auto uv_v1 = uv(c1);
        auto uv_v2 = uv(c2);
        auto uv_m0 = (uv_v0 + uv_v1) * 0.5;
        auto uv_m1 = (uv_v1 + uv_v2) * 0.5;
        auto uv_m2 = (uv_v2 + uv_v0) * 0.5;
        auto uv_c  = (uv_v0 + uv_v1 + uv_v2) / 3.0;

        int f_idx = i * 6;

        F_new.row(f_idx + 0) << v0, m0, c;
        uv_new((f_idx + 0) * 3 + 0) = uv_v0; uv_new((f_idx + 0) * 3 + 1) = uv_m0; uv_new((f_idx + 0) * 3 + 2) = uv_c;

        F_new.row(f_idx + 1) << m0, v1, c;
        uv_new((f_idx + 1) * 3 + 0) = uv_m0; uv_new((f_idx + 1) * 3 + 1) = uv_v1; uv_new((f_idx + 1) * 3 + 2) = uv_c;

        F_new.row(f_idx + 2) << v1, m1, c;
        uv_new((f_idx + 2) * 3 + 0) = uv_v1; uv_new((f_idx + 2) * 3 + 1) = uv_m1; uv_new((f_idx + 2) * 3 + 2) = uv_c;

        F_new.row(f_idx + 3) << m1, v2, c;
        uv_new((f_idx + 3) * 3 + 0) = uv_m1; uv_new((f_idx + 3) * 3 + 1) = uv_v2; uv_new((f_idx + 3) * 3 + 2) = uv_c;

        F_new.row(f_idx + 4) << v2, m2, c;
        uv_new((f_idx + 4) * 3 + 0) = uv_v2; uv_new((f_idx + 4) * 3 + 1) = uv_m2; uv_new((f_idx + 4) * 3 + 2) = uv_c;

        F_new.row(f_idx + 5) << m2, v0, c;
        uv_new((f_idx + 5) * 3 + 0) = uv_m2; uv_new((f_idx + 5) * 3 + 1) = uv_v0; uv_new((f_idx + 5) * 3 + 2) = uv_c;
    }

    auto hm_new = std::make_unique<Hmesh>(V_new, F_new);

    // --- Seam と Matching の再マッピング ---
    std::vector<bool> seam_new(hm_new->nE, false);
    VecXi matching_new = VecXi::Zero(hm_new->nE);
    bool has_matching = (matching_old.size() == hm.nE);

    for (int i = 0; i < hm_new->nE; ++i) {
        auto e_new = hm_new->edges[i];
        int ev0 = e_new.vert0().id;
        int ev1 = e_new.vert1().id;
        auto key = make_edge_key(ev0, ev1);

        if (edge_to_old_id.find(key) != edge_to_old_id.end()) {
            int old_id = edge_to_old_id[key];
            auto e_old = hm.edges[old_id];

            // Seam の復元
            if (!seam.empty() && seam[old_id]) {
                seam_new[i] = true;
            }

            // Matching の復元
            if (has_matching) {
                if (matching_old[old_id] == -1) {
                    matching_new[i] = -1; // 境界などの無効値はそのまま
                } else {
                    // 新しいエッジの face0 が、元の face0 由来かを判定
                    // 新しい面IDを6で割れば、親となる元の面IDがわかる
                    int p0 = (e_new.face0().id != -1) ? (e_new.face0().id / 6) : -1;

                    if (p0 == e_old.face0().id) {
                        // 進行方向が同じならそのまま
                        matching_new[i] = matching_old[old_id];
                    } else if (p0 == e_old.face1().id) {
                        // 進行方向が逆（face0とface1が反転して構築された）なら回転を反転
                        matching_new[i] = (-matching_old[old_id] + 1000 * rosyN) % rosyN;
                    } else {
                        matching_new[i] = matching_old[old_id];
                    }
                }
            }
        } else {
            // 元の面の内部に引かれた新しいエッジ
            // 同一面内のフィールドはフラットに接続するため period jump は無い (0)
            if (has_matching) matching_new[i] = 0;
        }
    }

    // --- Singular の再マッピング ---
    VecXi singular_new = VecXi::Zero(hm_new->nV);
    if (singular_old.size() == hm.nV) {
        // 特異点は元の頂点位置にのみ保存される（新しくできた中点や重心は常に0）
        singular_new.head(hm.nV) = singular_old;
    }

    return Subdivide{std::move(hm_new), uv_new, seam_new, matching_new, singular_new};
}

// 1-to-12 細分 (1-to-4 分割後、各サブ三角形を 1-to-3 重心細分)
// UV, Seam, Matching, Singular 対応版
inline Subdivide twelve_subdivide(
    const Hmesh& hm,
    const VecXc& uv,
    const std::vector<bool>& seam,
    const VecXi& matching_old = VecXi(),
    const VecXi& singular_old = VecXi(),
    int rosyN = 4
) {
    // 頂点数: 既存頂点 + エッジ中点 + 各面の4つの重心
    int nV_new = hm.nV + hm.nE + hm.nF * 4;
    // 面数: 1面につき12面
    int nF_new = hm.nF * 12;

    MatXd V_new(nV_new, 3);
    MatXi F_new(nF_new, 3);
    VecXc uv_new(nF_new * 3);

    int e_offset = hm.nV;
    int f_offset = hm.nV + hm.nE;

    // SeamとMatchingの引継ぎ用マップ: (min_v, max_v) -> old_edge_id
    std::map<std::pair<int, int>, int> edge_to_old_id;
    auto make_edge_key = [](int a, int b) {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };

    // 1. 既存の頂点のコピー
    for (int i = 0; i < hm.nV; ++i) {
        V_new.row(i) = hm.pos.row(i);
    }

    // 2. エッジ中点の座標計算と親子関係の記録
    for (int i = 0; i < hm.nE; ++i) {
        auto e = hm.edges[i];
        int v0 = e.vert0().id;
        int v1 = e.vert1().id;
        int m_id = e_offset + i;

        V_new.row(m_id) = (hm.pos.row(v0) + hm.pos.row(v1)) * 0.5;

        edge_to_old_id[make_edge_key(v0, m_id)] = i;
        edge_to_old_id[make_edge_key(m_id, v1)] = i;
    }

    // 3. 面の分割とコーナーUVの計算
    for (int i = 0; i < hm.nF; ++i) {
        auto f = hm.faces[i];
        auto h0 = f.half();
        auto h1 = h0.next();
        auto h2 = h1.next();

        int v0 = h0.tail().id;
        int v1 = h1.tail().id;
        int v2 = h2.tail().id;

        int m0 = e_offset + h0.edge().id; // 辺 v0-v1 の中点
        int m1 = e_offset + h1.edge().id; // 辺 v1-v2 の中点
        int m2 = e_offset + h2.edge().id; // 辺 v2-v0 の中点

        // 新しく追加される4つの重心のID
        int c0 = f_offset + i * 4 + 0;
        int c1 = f_offset + i * 4 + 1;
        int c2 = f_offset + i * 4 + 2;
        int c3 = f_offset + i * 4 + 3;

        auto p_v0 = hm.pos.row(v0);
        auto p_v1 = hm.pos.row(v1);
        auto p_v2 = hm.pos.row(v2);
        auto p_m0 = V_new.row(m0);
        auto p_m1 = V_new.row(m1);
        auto p_m2 = V_new.row(m2);

        // 4つの小三角形の重心位置
        V_new.row(c0) = (p_v0 + p_m0 + p_m2) / 3.0; // T0: v0側の角
        V_new.row(c1) = (p_m0 + p_v1 + p_m1) / 3.0; // T1: v1側の角
        V_new.row(c2) = (p_m2 + p_m1 + p_v2) / 3.0; // T2: v2側の角
        V_new.row(c3) = (p_m0 + p_m1 + p_m2) / 3.0; // T3: 中央の面

        // UVの計算
        auto uv_v0 = uv(i * 3 + 0);
        auto uv_v1 = uv(i * 3 + 1);
        auto uv_v2 = uv(i * 3 + 2);
        auto uv_m0 = (uv_v0 + uv_v1) * 0.5;
        auto uv_m1 = (uv_v1 + uv_v2) * 0.5;
        auto uv_m2 = (uv_v2 + uv_v0) * 0.5;

        auto uv_c0 = (uv_v0 + uv_m0 + uv_m2) / 3.0;
        auto uv_c1 = (uv_m0 + uv_v1 + uv_m1) / 3.0;
        auto uv_c2 = (uv_m2 + uv_m1 + uv_v2) / 3.0;
        auto uv_c3 = (uv_m0 + uv_m1 + uv_m2) / 3.0;

        int f_idx = i * 12;

        // --- T0 (v0側の三角形) を3分割 ---
        F_new.row(f_idx + 0) << v0, m0, c0;
        uv_new((f_idx + 0) * 3 + 0) = uv_v0; uv_new((f_idx + 0) * 3 + 1) = uv_m0; uv_new((f_idx + 0) * 3 + 2) = uv_c0;

        F_new.row(f_idx + 1) << m0, m2, c0;
        uv_new((f_idx + 1) * 3 + 0) = uv_m0; uv_new((f_idx + 1) * 3 + 1) = uv_m2; uv_new((f_idx + 1) * 3 + 2) = uv_c0;

        F_new.row(f_idx + 2) << m2, v0, c0;
        uv_new((f_idx + 2) * 3 + 0) = uv_m2; uv_new((f_idx + 2) * 3 + 1) = uv_v0; uv_new((f_idx + 2) * 3 + 2) = uv_c0;

        // --- T1 (v1側の三角形) を3分割 ---
        F_new.row(f_idx + 3) << m0, v1, c1;
        uv_new((f_idx + 3) * 3 + 0) = uv_m0; uv_new((f_idx + 3) * 3 + 1) = uv_v1; uv_new((f_idx + 3) * 3 + 2) = uv_c1;

        F_new.row(f_idx + 4) << v1, m1, c1;
        uv_new((f_idx + 4) * 3 + 0) = uv_v1; uv_new((f_idx + 4) * 3 + 1) = uv_m1; uv_new((f_idx + 4) * 3 + 2) = uv_c1;

        F_new.row(f_idx + 5) << m1, m0, c1;
        uv_new((f_idx + 5) * 3 + 0) = uv_m1; uv_new((f_idx + 5) * 3 + 1) = uv_m0; uv_new((f_idx + 5) * 3 + 2) = uv_c1;

        // --- T2 (v2側の三角形) を3分割 ---
        F_new.row(f_idx + 6) << m2, m1, c2;
        uv_new((f_idx + 6) * 3 + 0) = uv_m2; uv_new((f_idx + 6) * 3 + 1) = uv_m1; uv_new((f_idx + 6) * 3 + 2) = uv_c2;

        F_new.row(f_idx + 7) << m1, v2, c2;
        uv_new((f_idx + 7) * 3 + 0) = uv_m1; uv_new((f_idx + 7) * 3 + 1) = uv_v2; uv_new((f_idx + 7) * 3 + 2) = uv_c2;

        F_new.row(f_idx + 8) << v2, m2, c2;
        uv_new((f_idx + 8) * 3 + 0) = uv_v2; uv_new((f_idx + 8) * 3 + 1) = uv_m2; uv_new((f_idx + 8) * 3 + 2) = uv_c2;

        // --- T3 (中央の三角形) を3分割 ---
        F_new.row(f_idx + 9) << m0, m1, c3;
        uv_new((f_idx + 9) * 3 + 0) = uv_m0; uv_new((f_idx + 9) * 3 + 1) = uv_m1; uv_new((f_idx + 9) * 3 + 2) = uv_c3;

        F_new.row(f_idx + 10) << m1, m2, c3;
        uv_new((f_idx + 10) * 3 + 0) = uv_m1; uv_new((f_idx + 10) * 3 + 1) = uv_m2; uv_new((f_idx + 10) * 3 + 2) = uv_c3;

        F_new.row(f_idx + 11) << m2, m0, c3;
        uv_new((f_idx + 11) * 3 + 0) = uv_m2; uv_new((f_idx + 11) * 3 + 1) = uv_m0; uv_new((f_idx + 11) * 3 + 2) = uv_c3;
    }

    auto hm_new = std::make_unique<Hmesh>(V_new, F_new);

    // --- Seam と Matching の再マッピング ---
    std::vector<bool> seam_new(hm_new->nE, false);
    VecXi matching_new = VecXi::Zero(hm_new->nE);
    bool has_matching = (matching_old.size() == hm.nE);

    for (int i = 0; i < hm_new->nE; ++i) {
        auto e_new = hm_new->edges[i];
        int ev0 = e_new.vert0().id;
        int ev1 = e_new.vert1().id;
        auto key = make_edge_key(ev0, ev1);

        if (edge_to_old_id.find(key) != edge_to_old_id.end()) {
            int old_id = edge_to_old_id[key];
            auto e_old = hm.edges[old_id];

            // Seam の復元
            if (!seam.empty() && seam[old_id]) {
                seam_new[i] = true;
            }

            // Matching の復元
            if (has_matching) {
                if (matching_old[old_id] == -1) {
                    matching_new[i] = -1;
                } else {
                    // 面が12分割されたため、12で割って元の親面を特定
                    int p0 = (e_new.face0().id != -1) ? (e_new.face0().id / 12) : -1;

                    if (p0 == e_old.face0().id) {
                        matching_new[i] = matching_old[old_id];
                    } else if (p0 == e_old.face1().id) {
                        matching_new[i] = (-matching_old[old_id] + 1000 * rosyN) % rosyN;
                    } else {
                        matching_new[i] = matching_old[old_id];
                    }
                }
            }
        } else {
            // 元の面の内部エッジはジャンプなし
            if (has_matching) matching_new[i] = 0;
        }
    }

    // --- Singular の再マッピング ---
    VecXi singular_new = VecXi::Zero(hm_new->nV);
    if (singular_old.size() == hm.nV) {
        singular_new.head(hm.nV) = singular_old;
    }

    return Subdivide{std::move(hm_new), uv_new, seam_new, matching_new, singular_new};
}

// 1-to-12 細分 (4分割の真ん中を貫く中線分割)
// UV, Seam, Matching, Singular 対応版
inline Subdivide twelve_subdivide_2(
    const Hmesh& hm,
    const VecXc& uv,
    const std::vector<bool>& seam,
    const VecXi& matching_old = VecXi(),
    const VecXi& singular_old = VecXi(),
    int rosyN = 4
) {
    // 頂点数: 既存頂点 + 外周エッジ中点 + 各面の内部頂点4つ(重心1 + 内周エッジ中点3)
    int nV_new = hm.nV + hm.nE + hm.nF * 4;
    // 面数: 1面につき12面
    int nF_new = hm.nF * 12;

    MatXd V_new(nV_new, 3);
    MatXi F_new(nF_new, 3);
    VecXc uv_new(nF_new * 3);

    int e_offset = hm.nV;
    int f_offset = hm.nV + hm.nE;

    // SeamとMatchingの引継ぎ用マップ: (min_v, max_v) -> old_edge_id
    std::map<std::pair<int, int>, int> edge_to_old_id;
    auto make_edge_key = [](int a, int b) {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };

    // 1. 既存の頂点のコピー
    for (int i = 0; i < hm.nV; ++i) {
        V_new.row(i) = hm.pos.row(i);
    }

    // 2. 外周エッジ中点の座標計算と親子関係の記録
    for (int i = 0; i < hm.nE; ++i) {
        auto e = hm.edges[i];
        int v0 = e.vert0().id;
        int v1 = e.vert1().id;
        int m_id = e_offset + i;

        V_new.row(m_id) = (hm.pos.row(v0) + hm.pos.row(v1)) * 0.5;

        edge_to_old_id[make_edge_key(v0, m_id)] = i;
        edge_to_old_id[make_edge_key(m_id, v1)] = i;
    }

    // 3. 面の分割とコーナーUVの計算
    for (int i = 0; i < hm.nF; ++i) {
        auto f = hm.faces[i];
        auto h0 = f.half();
        auto h1 = h0.next();
        auto h2 = h1.next();

        int v0 = h0.tail().id;
        int v1 = h1.tail().id;
        int v2 = h2.tail().id;

        int m0 = e_offset + h0.edge().id; // 辺 v0-v1 の中点
        int m1 = e_offset + h1.edge().id; // 辺 v1-v2 の中点
        int m2 = e_offset + h2.edge().id; // 辺 v2-v0 の中点

        // 新しく追加される4つの内部頂点のID
        int m01 = f_offset + i * 4 + 0; // 内側エッジ m0-m1 の中点
        int m12 = f_offset + i * 4 + 1; // 内側エッジ m1-m2 の中点
        int m20 = f_offset + i * 4 + 2; // 内側エッジ m2-m0 の中点
        int c   = f_offset + i * 4 + 3; // 面全体の重心 (中線の交点)

        auto p_v0 = hm.pos.row(v0);
        auto p_v1 = hm.pos.row(v1);
        auto p_v2 = hm.pos.row(v2);
        auto p_m0 = V_new.row(m0);
        auto p_m1 = V_new.row(m1);
        auto p_m2 = V_new.row(m2);

        // 内部頂点の座標計算 (中線がまっすぐ貫くように配置)
        V_new.row(m01) = (p_m0 + p_m1) * 0.5;
        V_new.row(m12) = (p_m1 + p_m2) * 0.5;
        V_new.row(m20) = (p_m2 + p_m0) * 0.5;
        V_new.row(c)   = (p_v0 + p_v1 + p_v2) / 3.0;

        // UVの計算
        auto uv_v0 = uv(i * 3 + 0);
        auto uv_v1 = uv(i * 3 + 1);
        auto uv_v2 = uv(i * 3 + 2);
        auto uv_m0 = (uv_v0 + uv_v1) * 0.5;
        auto uv_m1 = (uv_v1 + uv_v2) * 0.5;
        auto uv_m2 = (uv_v2 + uv_v0) * 0.5;

        auto uv_m01 = (uv_m0 + uv_m1) * 0.5;
        auto uv_m12 = (uv_m1 + uv_m2) * 0.5;
        auto uv_m20 = (uv_m2 + uv_m0) * 0.5;
        auto uv_c   = (uv_v0 + uv_v1 + uv_v2) / 3.0;

        int f_idx = i * 12;

        // --- 外側の3つの小三角形をそれぞれ2分割 (中線が貫通する部分) ---
        // Corner 0 (v0周辺)
        F_new.row(f_idx + 0) << v0, m0, m20;
        uv_new((f_idx + 0) * 3 + 0) = uv_v0; uv_new((f_idx + 0) * 3 + 1) = uv_m0; uv_new((f_idx + 0) * 3 + 2) = uv_m20;

        F_new.row(f_idx + 1) << v0, m20, m2;
        uv_new((f_idx + 1) * 3 + 0) = uv_v0; uv_new((f_idx + 1) * 3 + 1) = uv_m20; uv_new((f_idx + 1) * 3 + 2) = uv_m2;

        // Corner 1 (v1周辺)
        F_new.row(f_idx + 2) << v1, m1, m01;
        uv_new((f_idx + 2) * 3 + 0) = uv_v1; uv_new((f_idx + 2) * 3 + 1) = uv_m1; uv_new((f_idx + 2) * 3 + 2) = uv_m01;

        F_new.row(f_idx + 3) << v1, m01, m0;
        uv_new((f_idx + 3) * 3 + 0) = uv_v1; uv_new((f_idx + 3) * 3 + 1) = uv_m01; uv_new((f_idx + 3) * 3 + 2) = uv_m0;

        // Corner 2 (v2周辺)
        F_new.row(f_idx + 4) << v2, m2, m12;
        uv_new((f_idx + 4) * 3 + 0) = uv_v2; uv_new((f_idx + 4) * 3 + 1) = uv_m2; uv_new((f_idx + 4) * 3 + 2) = uv_m12;

        F_new.row(f_idx + 5) << v2, m12, m1;
        uv_new((f_idx + 5) * 3 + 0) = uv_v2; uv_new((f_idx + 5) * 3 + 1) = uv_m12; uv_new((f_idx + 5) * 3 + 2) = uv_m1;

        // --- 内側の真ん中(4分割の中心)の小三角形を6分割 ---
        F_new.row(f_idx + 6) << c, m20, m0;
        uv_new((f_idx + 6) * 3 + 0) = uv_c; uv_new((f_idx + 6) * 3 + 1) = uv_m20; uv_new((f_idx + 6) * 3 + 2) = uv_m0;

        F_new.row(f_idx + 7) << c, m0, m01;
        uv_new((f_idx + 7) * 3 + 0) = uv_c; uv_new((f_idx + 7) * 3 + 1) = uv_m0; uv_new((f_idx + 7) * 3 + 2) = uv_m01;

        F_new.row(f_idx + 8) << c, m01, m1;
        uv_new((f_idx + 8) * 3 + 0) = uv_c; uv_new((f_idx + 8) * 3 + 1) = uv_m01; uv_new((f_idx + 8) * 3 + 2) = uv_m1;

        F_new.row(f_idx + 9) << c, m1, m12;
        uv_new((f_idx + 9) * 3 + 0) = uv_c; uv_new((f_idx + 9) * 3 + 1) = uv_m1; uv_new((f_idx + 9) * 3 + 2) = uv_m12;

        F_new.row(f_idx + 10) << c, m12, m2;
        uv_new((f_idx + 10) * 3 + 0) = uv_c; uv_new((f_idx + 10) * 3 + 1) = uv_m12; uv_new((f_idx + 10) * 3 + 2) = uv_m2;

        F_new.row(f_idx + 11) << c, m2, m20;
        uv_new((f_idx + 11) * 3 + 0) = uv_c; uv_new((f_idx + 11) * 3 + 1) = uv_m2; uv_new((f_idx + 11) * 3 + 2) = uv_m20;
    }

    auto hm_new = std::make_unique<Hmesh>(V_new, F_new);

    // --- Seam と Matching の再マッピング ---
    std::vector<bool> seam_new(hm_new->nE, false);
    VecXi matching_new = VecXi::Zero(hm_new->nE);
    bool has_matching = (matching_old.size() == hm.nE);

    for (int i = 0; i < hm_new->nE; ++i) {
        auto e_new = hm_new->edges[i];
        int ev0 = e_new.vert0().id;
        int ev1 = e_new.vert1().id;
        auto key = make_edge_key(ev0, ev1);

        // 外周のオリジナルエッジだった場合のみ転写
        if (edge_to_old_id.find(key) != edge_to_old_id.end()) {
            int old_id = edge_to_old_id[key];
            auto e_old = hm.edges[old_id];

            if (!seam.empty() && seam[old_id]) {
                seam_new[i] = true;
            }

            if (has_matching) {
                if (matching_old[old_id] == -1) {
                    matching_new[i] = -1;
                } else {
                    // 面が12分割されたため、12で割って元の親面を特定
                    int p0 = (e_new.face0().id != -1) ? (e_new.face0().id / 12) : -1;

                    if (p0 == e_old.face0().id) {
                        matching_new[i] = matching_old[old_id];
                    } else if (p0 == e_old.face1().id) {
                        matching_new[i] = (-matching_old[old_id] + 1000 * rosyN) % rosyN;
                    } else {
                        matching_new[i] = matching_old[old_id];
                    }
                }
            }
        } else {
            // 元の面の内部エッジはジャンプなし
            if (has_matching) matching_new[i] = 0;
        }
    }

    // --- Singular の再マッピング ---
    VecXi singular_new = VecXi::Zero(hm_new->nV);
    if (singular_old.size() == hm.nV) {
        singular_new.head(hm.nV) = singular_old; // 特異点は元の頂点位置にのみ引き継ぐ
    }

    return Subdivide{std::move(hm_new), uv_new, seam_new, matching_new, singular_new};
}

// 1-to-12 細分 (4分割の真ん中を貫く中線分割)
// UV, Seam, Matching, Singular 対応版
inline std::unique_ptr<Hmesh> twelve_subdivide_3(
    const Hmesh& hm,
    const std::vector<bool>& seam,
    std::vector<bool>& new_seam
) {
    // 頂点数: 既存頂点 + 外周エッジ中点 + 各面の内部頂点4つ(重心1 + 内周エッジ中点3)
    int nV_new = hm.nV + hm.nE + hm.nF * 4;
    // 面数: 1面につき12面
    int nF_new = hm.nF * 12;

    MatXd V_new(nV_new, 3);
    MatXi F_new(nF_new, 3);

    int e_offset = hm.nV;
    int f_offset = hm.nV + hm.nE;

    // SeamとMatchingの引継ぎ用マップ: (min_v, max_v) -> old_edge_id
    std::map<std::pair<int, int>, int> edge_to_old_id;
    auto make_edge_key = [](int a, int b) {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };

    // 1. 既存の頂点のコピー
    for (int i = 0; i < hm.nV; ++i) {
        V_new.row(i) = hm.pos.row(i);
    }

    // 2. 外周エッジ中点の座標計算と親子関係の記録
    for (int i = 0; i < hm.nE; ++i) {
        auto e = hm.edges[i];
        int v0 = e.vert0().id;
        int v1 = e.vert1().id;
        int m_id = e_offset + i;

        V_new.row(m_id) = (hm.pos.row(v0) + hm.pos.row(v1)) * 0.5;

        edge_to_old_id[make_edge_key(v0, m_id)] = i;
        edge_to_old_id[make_edge_key(m_id, v1)] = i;
    }

    // 3. 面の分割とコーナーUVの計算
    for (int i = 0; i < hm.nF; ++i) {
        auto f = hm.faces[i];
        auto h0 = f.half();
        auto h1 = h0.next();
        auto h2 = h1.next();

        int v0 = h0.tail().id;
        int v1 = h1.tail().id;
        int v2 = h2.tail().id;

        int m0 = e_offset + h0.edge().id; // 辺 v0-v1 の中点
        int m1 = e_offset + h1.edge().id; // 辺 v1-v2 の中点
        int m2 = e_offset + h2.edge().id; // 辺 v2-v0 の中点

        // 新しく追加される4つの内部頂点のID
        int m01 = f_offset + i * 4 + 0; // 内側エッジ m0-m1 の中点
        int m12 = f_offset + i * 4 + 1; // 内側エッジ m1-m2 の中点
        int m20 = f_offset + i * 4 + 2; // 内側エッジ m2-m0 の中点
        int c   = f_offset + i * 4 + 3; // 面全体の重心 (中線の交点)

        auto p_v0 = hm.pos.row(v0);
        auto p_v1 = hm.pos.row(v1);
        auto p_v2 = hm.pos.row(v2);
        auto p_m0 = V_new.row(m0);
        auto p_m1 = V_new.row(m1);
        auto p_m2 = V_new.row(m2);

        // 内部頂点の座標計算 (中線がまっすぐ貫くように配置)
        V_new.row(m01) = (p_m0 + p_m1) * 0.5;
        V_new.row(m12) = (p_m1 + p_m2) * 0.5;
        V_new.row(m20) = (p_m2 + p_m0) * 0.5;
        V_new.row(c)   = (p_v0 + p_v1 + p_v2) / 3.0;

        int f_idx = i * 12;

        // --- 外側の3つの小三角形をそれぞれ2分割 (中線が貫通する部分) ---
        // Corner 0 (v0周辺)
        F_new.row(f_idx + 0) << v0, m0, m20;
        F_new.row(f_idx + 1) << v0, m20, m2;

        // Corner 1 (v1周辺)
        F_new.row(f_idx + 2) << v1, m1, m01;
        F_new.row(f_idx + 3) << v1, m01, m0;

        // Corner 2 (v2周辺)
        F_new.row(f_idx + 4) << v2, m2, m12;
        F_new.row(f_idx + 5) << v2, m12, m1;

        // --- 内側の真ん中(4分割の中心)の小三角形を6分割 ---
        F_new.row(f_idx + 6) << c, m20, m0;
        F_new.row(f_idx + 7) << c, m0, m01;
        F_new.row(f_idx + 8) << c, m01, m1;
        F_new.row(f_idx + 9) << c, m1, m12;
        F_new.row(f_idx + 10) << c, m12, m2;
        F_new.row(f_idx + 11) << c, m2, m20;
    }

    auto hm_new = std::make_unique<Hmesh>(V_new, F_new);

    // --- Seam と Matching の再マッピング ---
    std::vector<bool> seam_new(hm_new->nE, false);
    VecXi matching_new = VecXi::Zero(hm_new->nE);

    for (int i = 0; i < hm_new->nE; ++i) {
        auto e_new = hm_new->edges[i];
        int ev0 = e_new.vert0().id;
        int ev1 = e_new.vert1().id;
        auto key = make_edge_key(ev0, ev1);

        // 外周のオリジナルエッジだった場合のみ転写
        if (edge_to_old_id.find(key) != edge_to_old_id.end()) {
            int old_id = edge_to_old_id[key];
            auto e_old = hm.edges[old_id];

            if (!seam.empty() && seam[old_id]) {
                seam_new[i] = true;
            }
        }
    }

    // --- Singular の再マッピング ---
    VecXc v1 = VecXc::Zero(hm_new->nC);
    VecXi v2 = VecXi::Zero(hm_new->nV);
    VecXi v3 = VecXi::Zero(hm_new->nV);
    new_seam = seam_new;

    return hm_new;
}

}


#endif