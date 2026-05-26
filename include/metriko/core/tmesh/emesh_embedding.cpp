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

    pq.push({bgn_vid, 0.});
    min_cost[bgn_vid] = 0.;

    vec<int> valid_fids;
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

        if (reached) break;
        if (curr_vid == end_vid) { reached = true; break; }
        if (curr_cost > min_cost[curr_vid]) continue;

        for (auto h : hm.verts[curr_vid].adjHalfs()) {
            int next_vid = h.head().id;
            int this_fid = -1;
            int prnt_fid = -1;

            if (next_vid == end_vid) {
                came_from_hid[next_vid] = h.id;
                reached = true;
            }

            if (next_vid != end_vid && occupied_verts[next_vid]) continue;
            int fL = h.face().id;
            int fR = h.twin().face().id;
            int pL = data.face2parent[fL];
            int pR = data.face2parent[fR];

            if      (rg::find(valid_fids, pL) != valid_fids.end()) { this_fid = fL; prnt_fid = pL; }
            else if (rg::find(valid_fids, pR) != valid_fids.end()) { this_fid = fR; prnt_fid = pR; }
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
                pq.push({next_vid, c});
                came_from_hid[next_vid] = h.id;
            }
        }
    }

    if (!reached) { std::cout << "Dijkstra failed!" << std::endl; return {}; }

    vec<Half> path;
    int curr = end_vid;
    while (curr != bgn_vid) {
        int hid = came_from_hid[curr];
        Half h_ = hm.halfs[hid];
        path.push_back(h_);
        occupied_verts[h_.tail().id] = true;
        occupied_verts[h_.head().id] = true;
        curr = h_.tail().id;
    }
    rg::reverse(path);
    return path;
}


void Emesh::assign_first_half(
    const Tmesh& tm,
    const Mgrph& mg,
    const std::map<int, int>& mnode2dense_v,
    const TrackedDenseMesh& data
) {
    for (int vid: sings) {
        auto tgts = eedges | vw::filter([&](const auto& ee) {
            return mnode2dense_v.at(ee.fr_nid) == vid || mnode2dense_v.at(ee.to_nid) == vid;
        });
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
            auto s = &tm.tedges[ee.id].segs.front();
            auto uvFr = get_face_uv(mg.mnodes[s->fr_nid], s->face_id, mg.hm, mg.cf);
            auto uvTo = get_face_uv(mg.mnodes[s->to_nid], s->face_id, mg.hm, mg.cf);

            for (auto h: v.adjHalfs()) {
                Face f = h.face();
                assert(s->face_id >= 0);
                if (data.face2parent[f.id] == s->face_id) {
                    auto cL = h.crnr();
                    auto cR = h.prev().crnr();
                    auto uvL = data.uvs[f.id][cL.id % 3];
                    auto uvR = data.uvs[f.id][cR.id % 3];
                    auto dL = compute_point_to_segment_distance(uvL, uvFr, uvTo);
                    auto dR = compute_point_to_segment_distance(uvR, uvFr, uvTo);
                    if (is_points_into(uvFr, uvR, uvL, uvTo, 0)) {
                        ee_data.push_back(dR < dL ? Data{h, f, ee.id, dL} : Data{h.prev().twin(), f, ee.id, dR});
                    }
                }
            }
        }

        int count1 = 0;
        int count2 = 0;
        for (auto _: v.adjHalfs()) count1++;
        for (auto _: tgts) count2++;
        assert(count1 >= count2);


        /*
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

            if (flag) {
                //std::cout << "singular vid: " << vid << " done!" << std::endl;
                break;
            }
        }
        */
        for (const Data& d : ee_data) { eedges[d.eeid].halfs.push_back(d.h); }
    }
}

void Emesh::assign_inter_half(
    const Tmesh& tm,
    const Mgrph& mg,
    const std::map<int, int>& mnode2dense_v,
    const TrackedDenseMesh& data
) {
    vec<bool> occupied_verts = vec(hm.verts.size(), false);

    vec<const Tedge*> ptrs;
    ptrs.reserve(tm.tedges.size());
    for (const auto& te : tm.tedges) { ptrs.push_back(&te); }
    rg::sort(ptrs, {}, [&](const auto* e) {
        int priority = 2;
        if (mg.mnodes[e->fr_nid].jt == JunctionType::F) priority = 0;
        if (mg.mnodes[e->to_nid].jt == JunctionType::T) priority = 1;
        return std::pair{priority, e->len};
    });

    // eedges from singular points first
    for (const auto te: ptrs) {
        auto& ee = eedges[te->id];
        if (mg.mnodes[te->fr_nid].jt != JunctionType::F) continue;
        int i0 = ee.halfs.back().head().id;
        int i1 = mnode2dense_v.at(te->to_nid);
        vec<Half> res = compute_dijkstra_snap(*te, mg, hm, data, occupied_verts, i0, i1);
        ee.halfs.insert(ee.halfs.end(), res.begin(), res.end());
    }

    // other eedges
    for (const auto te: ptrs) {
        auto& ee = eedges[te->id];
        if (mg.mnodes[te->fr_nid].jt == JunctionType::F) continue;
        int i0 = mnode2dense_v.at(te->fr_nid);
        int i1 = mnode2dense_v.at(te->to_nid);
        vec<Half> res = compute_dijkstra_snap(*te, mg, hm, data, occupied_verts, i0, i1);
        ee.halfs.insert(ee.halfs.end(), res.begin(), res.end());
    }
}

void Emesh::assign_last_half(
    const Tmesh& tm,
    const Mgrph& mg,
    const std::map<int, int>& mnode2dense_v,
    const TrackedDenseMesh& data
) {
    for (int vid: crashes) {
        auto tgts = eedges | vw::filter([&](const auto& ee) {
            return mnode2dense_v.at(ee.fr_nid) == vid || mnode2dense_v.at(ee.to_nid) == vid;
        });
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
            auto s = &tm.tedges[ee.id].segs.back();
            auto uvFr = get_face_uv(mg.mnodes[s->to_nid], s->face_id, mg.hm, mg.cf); // todo: use to_nid temp
            auto uvTo = get_face_uv(mg.mnodes[s->fr_nid], s->face_id, mg.hm, mg.cf); // todo: use fr_nid temp

            for (auto h: v.adjHalfs()) {
                Face f = h.face();
                assert(s->face_id >= 0);
                if (data.face2parent[f.id] == s->face_id) {
                    auto cL = h.crnr();
                    auto cR = h.prev().crnr();
                    auto uvL = data.uvs[f.id][cL.id % 3];
                    auto uvR = data.uvs[f.id][cR.id % 3];
                    auto dL = compute_point_to_segment_distance(uvL, uvFr, uvTo);
                    auto dR = compute_point_to_segment_distance(uvR, uvFr, uvTo);
                    if (is_points_into(uvFr, uvR, uvL, uvTo, 0)) {
                        ee_data.push_back(dR < dL ? Data{h, f, ee.id, dL} : Data{h.prev().twin(), f, ee.id, dR});
                    }
                }
            }
        }

        for (const Data& d : ee_data) { eedges[d.eeid].halfs.push_back(d.h); }
    }
}


}
