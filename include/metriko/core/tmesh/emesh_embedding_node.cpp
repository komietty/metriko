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
    // Move to another vert if a mnode is movable to other verts and if a vert is close to the bottom edge of the t-junction
    for (auto& [id, vid] : mnode2dense_v) {
        std::cout << "id: " << id << ", vid: " << vid << std::endl;
        auto mn = mg.mnodes[id];
        if (mn.jt == JunctionType::T && mn.adj.size() == 3) {

            int cid = mn.adj[1].curv_id;
            int sid = mn.adj[1].sgmt_id;
            auto& s = mg.mcurvs[cid].sgmts[sid];
            auto r = map_node(s.fr_nid, s.face_id);
            if (r.has_value()) {
                auto [vid_, dist] = r.value();
                std::cout << "candidate " << vid_ << std::endl;
                if (!occupied_v[vid_]) {
                    std::cout << "moved to " << vid_ << std::endl;
                    occupied_v[vid] = false;
                    occupied_v[vid_] = true;
                    mnode2dense_v[id] = vid_;
                }
            }
        }
    }


    return mnode2dense_v;
}
}
