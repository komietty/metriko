//
// Created by saki on 2026/05/22.
//

#include "./emesh.h"
using namespace metriko::mc;

namespace metriko {

static double compute_point_to_segment_distance(complex p, complex a, complex b) {
    complex ab = b - a;
    complex ap = p - a;
    if (std::abs(ab) < 1e-8) return std::abs(ap);
    auto t = (ap.real() * ab.real() + ap.imag() * ab.imag()) / std::norm(ab);
    auto proj = a + std::clamp(t, 0., 1.) * ab;
    return std::abs(p - proj);
}


void Emesh::assign_first_half(
    const Tmesh& tm,
    const Mgrph& mg,
    const std::map<int, int>& mnode2dense_v,
    const TrackedDenseMesh& data
) {
    for (int vid: sings) {
        auto tgts = eedges | vw::filter([&](const auto& ee) { return mnode2dense_v.at(ee.fr_nid) == vid; });
        auto v = hm.verts[vid];

        struct Data {
            Half h;   // current half
            Face f;   // original face
            int eeid; //
            double dist;
        };

        vec<Data> ee_data;

       // 1: assign closest half to each eedge with some overlaps
        for (auto& ee: tgts) {
            for (auto h: v.adjHalfs()) {
                auto s = &tm.tedges[ee.id].segs.front();
                Face f = h.face();
                if (data.face2parent[f.id] == s->face_id) {
                    auto uv = get_face_uv(mg.mnodes[s->to_nid], s->face_id, mg.hm, mg.cf);
                    auto cA = h.crnr();
                    auto cB = h.prev().crnr();
                    auto dA = abs(data.uvs[cA.face().id][cA.vert().id] - uv);
                    auto dB = abs(data.uvs[cB.face().id][cB.vert().id] - uv);
                    ee_data.push_back(dB < dA ? Data{h, f, ee.id, dB} : Data{h.prev().twin(), f, ee.id, dA});
                }
            }
        }

        int count1 = 0;
        int count2 = 0;
        for (auto _: v.adjHalfs()) count1++;
        for (auto _: tgts) count2++;
        assert(count1 >= count2);


        // 2: resolve overlaps
        while (true) {
            // if any ee in tgt overlaps, try switching to another side of half
            bool flag = true;

            for (auto h: v.adjHalfs()) {
                vec<Data*> overlaps;
                for (auto& d: ee_data) if (d.h == h) overlaps.push_back(&d);

                if (overlaps.size() > 1) {
                    assert(overlaps.size() == 2);
                    bool in_face_0 = overlaps[0]->f == h.face();
                    auto& d0 = in_face_0 ? overlaps[0] : overlaps[1]; // only ccw move is ok
                    auto& d1 = in_face_0 ? overlaps[1] : overlaps[0]; // only cw move is ok
                    if (d0->dist >= d1->dist) { d0->dist = -1; d0->h = h.prev().twin(); }
                    else                      { d1->dist = -1; d1->h = h.twin().next(); }
                    flag = false;
                }
            }

            if (flag) { std::cout << "singular vid: " << vid << " done!" << std::endl; break; }
        }

        for (const Data& d : ee_data) { eedges[d.eeid].halfs.push_back(d.h); }
    }
}

void Emesh::assign_last_half(
    const Tmesh& tm,
    const Mgrph& mg,
    const std::map<int, int>& mnode2dense_v
) {

}

static vec<Half> compute_dijkstra_snap(
    const Tedge& te,
    const Mgrph& mg,
    const Hmesh& hm,
    const TrackedDenseMesh& data,
    vec<bool>& occupied_verts,
    int bgn_vid,
    int end_vid
) {
    struct Node {
        int vid;
        double cost;
        bool operator>(const Node& rhs) const { return cost > rhs.cost; }
    };

    auto get_uv = [&](int dense_fid, int vid) -> complex {
        for (int i = 0; i < 3; ++i) { if (data.polygons[dense_fid][i] == vid) return data.uvs[dense_fid][i]; }
        return {0, 0};
    };

    std::priority_queue<Node, vec<Node>, std::greater<>> pq;
    vec min_cost(hm.verts.size(), std::numeric_limits<double>::infinity());
    vec came_from_hid(hm.verts.size(), -1);
    vec came_from_vid(hm.verts.size(), -1);

    pq.push({bgn_vid, 0.});
    min_cost[bgn_vid] = 0.;

    std::vector<int> valid_fids;
    for (const auto& sg : te.segs) {
        valid_fids.push_back(sg.face_id);
        for (auto h : mg.hm.faces[sg.face_id].adjHalfs()) valid_fids.push_back(h.twin().face().id);
    }
    rg::sort(valid_fids);
    valid_fids.erase(rg::unique(valid_fids).begin(), valid_fids.end());

    bool reached = false;

    while (!pq.empty()) {
        auto [curr_vid, curr_cost] = pq.top();
        pq.pop();

        if (curr_vid == end_vid) { reached = true; break; }
        if (curr_cost > min_cost[curr_vid]) continue;

        for (auto h : hm.verts[curr_vid].adjHalfs()) {
            if (h.head().id != end_vid && occupied_verts[h.head().id]) continue;
            int f_l = h.face().id;
            int f_r = h.twin().face().id;
            int p_l = data.face2parent[f_l];
            int p_r = data.face2parent[f_r];
            int next_vid = h.head().id;
            int this_fid = -1;
            int prnt_fid = -1;

            if      (rg::find(valid_fids, p_l) != valid_fids.end()) { this_fid = f_l; prnt_fid = p_l; }
            else if (rg::find(valid_fids, p_r) != valid_fids.end()) { this_fid = f_r; prnt_fid = p_r; }
            else { continue; }

            auto curr_uv = get_uv(this_fid, curr_vid);
            auto next_uv = get_uv(this_fid, next_vid);
            auto dist = std::abs(next_uv - curr_uv);

            auto  it = rg::find_if(te.segs, [&](const auto& s) { return s.face_id == prnt_fid; });
            auto* sg = it != te.segs.end() ? &(*it) : &te.segs.front();

            auto uv0 = get_face_uv(mg.mnodes[sg->fr_nid], sg->face_id, mg.hm, mg.cf);
            auto uv1 = get_face_uv(mg.mnodes[sg->to_nid], sg->face_id, mg.hm, mg.cf);
            auto d = compute_point_to_segment_distance(next_uv, uv0, uv1);
            auto c = curr_cost + dist * (1 + 100 * pow(d, 2));

            if (c < min_cost[next_vid]) {
                min_cost[next_vid] = c;
                came_from_hid[next_vid] = h.id;
                came_from_vid[next_vid] = curr_vid;
                pq.push({next_vid, c});
            }
        }
    }

    if (!reached) {
        std::cout << "Dijkstra failed!" << std::endl;
        return {};
    }

    vec<Half> path;
    int curr = end_vid;
    while (curr != bgn_vid) {
        int hid = came_from_hid[curr];
        Half h_ = hm.halfs[hid];
        path.push_back(h_);
        occupied_verts[h_.tail().id] = true;
        occupied_verts[h_.head().id] = true;
        curr = came_from_vid[curr];
    }
    rg::reverse(path);
    return path;
}

}
