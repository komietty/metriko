#include "metriko/core/tmesh/emesh.h"

namespace metriko {
using namespace metriko::mc;

std::map<int, int> mnode2dense_v(
    const Tmesh& tm,
    const Mgrph& mg,
    const Hmesh& hm,
    const VecXc& uv2,
    const TrackedDenseMesh& sdiv_data
) {
    std::map<int, int> mnode2dense_v;
    vec occupied_v(sdiv_data.num_verts, false);

    auto map_node = [&](int nid, int fid) -> std::optional<std::pair<int, double>> {
        if (mnode2dense_v.contains(nid)) return std::nullopt;
        complex target_uv = get_face_uv(mg.mnodes[nid], fid, hm, uv2);

        int best_vid = -1;
        double min_dist = 1e9;
        for (size_t i = 0; i < sdiv_data.polygons.size(); ++i) {
            if (sdiv_data.face2parent[i] != fid) continue;
            for (int j = 0; j < 3; ++j) {
                double d = std::abs(sdiv_data.uvs[i][j] - target_uv);
                if (d < min_dist && !occupied_v[sdiv_data.polygons[i][j]]) {
                    min_dist = d;
                    best_vid = sdiv_data.polygons[i][j];
                }
            }
        }
        assert(best_vid != -1);
        return std::pair{best_vid, min_dist};
    };

    // for the singular mnodes
    for (const auto& te: tm.tedges) {
        auto id = te.fr_nid;
        auto mn = mg.mnodes[id];
        if (mg.mnodes[id].jt == JunctionType::F) {
            auto r = map_node(id, te.segs.front().face_id);
            if (r.has_value()) {
                int vid = r.value().first;
                mnode2dense_v[id] = vid;
                occupied_v[vid] = true;
            }
        }
    }

    vec<std::pair<int, int>> unassigned_mnodes;
    for (const auto& te : tm.tedges) {
        if (mg.mnodes[te.to_nid].jt == JunctionType::T)
            unassigned_mnodes.emplace_back(te.to_nid, te.segs.back().face_id);
    }

    rg::sort(unassigned_mnodes);
    unassigned_mnodes.erase(rg::unique(unassigned_mnodes).begin(), unassigned_mnodes.end());

    while (!unassigned_mnodes.empty()) {
        // vid -> 提案してきた mnode のリスト {mnode_id, face_id, distance}
        std::map<int, vec<std::tuple<int, int, double>>> bids;

        // 1. 全ての未割り当て mnode が一番近い vid に入札する
        for (const auto& [id, fid] : unassigned_mnodes) {
            auto r = map_node(id, fid);
            if (r.has_value()) {
                auto [vid, dist] = r.value();
                bids[vid].emplace_back(id, fid, dist);
            }
        }

        vec<std::pair<int, int>> next_unassigned;

        // 2. 入札の競合を解決する
        for (auto& [vid, proposers] : bids) {
            // 距離でソート（一番近い mnode が先頭に来るようにする）
            rg::sort(proposers, [](const auto& a, const auto& b) { return std::get<2>(a) < std::get<2>(b); });
            int winner_id = std::get<0>(proposers.front());
            mnode2dense_v[winner_id] = vid;
            occupied_v[vid] = true;

            for (size_t i = 1; i < proposers.size(); ++i)
                next_unassigned.emplace_back(std::get<0>(proposers[i]), std::get<1>(proposers[i]));
        }
        unassigned_mnodes = std::move(next_unassigned);
    }

    // Relaxation process.
    vec<vec<int>> dense_adj(sdiv_data.num_verts);
    for (const auto& poly : sdiv_data.polygons) {
        for (int i = 0; i < poly.size(); ++i) {
            int v0 = poly[i];
            int v1 = poly[(i + 1) % poly.size()];
            dense_adj[v0].push_back(v1);
            dense_adj[v1].push_back(v0);
        }
    }
    for (auto& adj : dense_adj) { rg::sort(adj); adj.erase(rg::unique(adj).begin(), adj.end()); }

    for (auto& [nid, vid] : mnode2dense_v) {
        auto& mn = mg.mnodes[nid];
        if (mn.jt == JunctionType::T && mn.adj.size() == 3) {

            // 3つの adj の中から、他の二つと直角なもの（branch）を選択
            int branch_idx = -1;
            double min_dot = 1e9;
            for (int i = 0; i < 3; ++i) {
                // i 番目以外の二つの間の角度（ドット積）を確認
                int j = (i + 1) % 3;
                int k = (i + 2) % 3;

                auto get_v = [&](int idx) {
                    auto& as = mn.adj[idx];
                    auto& sg = mg.mcurvs[as.curv_id].sgmts[as.sgmt_id];
                    int other = (sg.fr_nid == nid) ? sg.to_nid : sg.fr_nid;
                    complex d = get_face_uv(mg.mnodes[other], sg.face_id, hm, uv2) - get_face_uv(
                        mn, sg.face_id, hm, uv2);
                    return d / std::abs(d);
                };

                double dot_jk = (get_v(j) * std::conj(get_v(k))).real();
                if (dot_jk < min_dot) { // 最も反対方向（-1に近い）を向いているペアの「相方」が branch
                    min_dot = dot_jk;
                    branch_idx = i;
                }
            }

            int cid = mn.adj[branch_idx].curv_id;
            int sid = mn.adj[branch_idx].sgmt_id;
            auto& s = mg.mcurvs[cid].sgmts[sid];

            //int cid = mn.adj[1].curv_id;
            //int sid = mn.adj[1].sgmt_id;
            //auto& s = mg.mcurvs[cid].sgmts[sid];

            auto uv0 = get_face_uv(mg.mnodes[s.to_nid], s.face_id, hm, uv2); // center
            auto uv1 = get_face_uv(mg.mnodes[s.fr_nid], s.face_id, hm, uv2); // target
            auto dir = (uv1 - uv0) / std::abs(uv1 - uv0);
            auto best_vid = vid;
            auto get_dist = [&](int v) -> double {
                for (int i = 0; i < sdiv_data.polygons.size(); ++i) {
                    if (sdiv_data.face2parent[i] != s.face_id) continue;
                    for (int j = 0; j < 3; ++j) {
                        if (sdiv_data.polygons[i][j] == v) {
                            //auto uv = sdiv_data.uvs[i][j];
                            //auto d2 = (uv - uv0) / std::abs(uv - uv0);
                            //double dot = d2.real() * dir.real() + d2.imag() * dir.imag();
                            //if (dot < 0) continue;
                            //return dot;
                            return std::abs(sdiv_data.uvs[i][j] - uv1);
                        }
                    }
                }
                return 1e18;
            };

            double min_dist = get_dist(vid);
            for (int vid_ : dense_adj[vid]) {
                if (occupied_v[vid_]) continue;
                double d = get_dist(vid_);
                if (d < min_dist) { min_dist = d; best_vid = vid_; }
            }

            if (best_vid != vid) {
                occupied_v[vid] = false;
                vid = best_vid;
                occupied_v[vid] = true;
            }
        }
    }


    return mnode2dense_v;
}
}
