//
//--- Copyright (C) 2026 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
#ifndef METRIKO_EMESH_H
#define METRIKO_EMESH_H

#include "tmesh.h"
#include "emesh_subdivide.h"
#include <queue>
#include <map>
#include <iostream>
#include <polyscope/point_cloud.h> // テストデバッグ可視化用

namespace metriko {
struct Emesh;

struct Eedge {
    int id = -1;
    int fr_nid = -1;
    int to_nid = -1;
    double len = -1;
    vec<Half> halfs;
};

struct Ehalf {
    const Emesh* em = nullptr;
    int id = -1;
    int twid = -1;
    int eeid = -1;
    int eqid = -1;
    int nxt_id = -1;
    int prv_id = -1;
    bool cano = false;
    bool bgn  = false;
    bool end  = false;
    double x = -1;
    vec<Half> halfs;

    const Ehalf& twin() const;
    const Ehalf& next() const;
    const Ehalf& prev() const;
    const Eedge& edge() const;
    int  nid_fr() const { return cano ? edge().fr_nid : edge().to_nid; }
    int  nid_to() const { return cano ? edge().to_nid : edge().fr_nid; }
    Vert tail() const;
    Vert head() const;
    void extend_prev(Ehalf prev) { halfs.insert(halfs.begin(), prev.halfs.begin(), prev.halfs.end()); }
    void extend_next(Ehalf next) { halfs.insert(halfs.end()  , next.halfs.begin(), next.halfs.end()); }
};

struct Edata {
    int ehid;
    int side;
};

struct Equad {
    int id = -1;
    std::vector<Edata> data;

    vec<int> ehids(int side) const {
        return data
            | vw::filter([&](const auto& d) { return d.side == side; })
            | vw::transform([](const auto& d) { return d.ehid; })
            | rg::to<vec<int>>();
    }

};

struct Emesh {
    const     Hmesh& hm;
    std::set<int> sings;
    vec<Equad> equads;
    vec<Ehalf> ehalfs;
    vec<Eedge> eedges;

    explicit Emesh(
        const     Tmesh& tm, //
        const mc::Mgrph& mg, //
        const     Hmesh& hm, // divided hm
        const TrackedDenseMesh& data,
        const std::map<int, int>& mnode2dense_v,
        const VecXd& X
    ) :hm(hm) {
        vec<int> occupied_halfs = {};

        for (const auto& te: tm.tedges) {
            auto i0 = mnode2dense_v.at(te.fr_nid);
            auto i1 = mnode2dense_v.at(te.to_nid);
            eedges.push_back({
                .id     = te.id,
                .fr_nid = te.fr_nid,
                .to_nid = te.to_nid,
                .len    = te.len,
                .halfs  = compute_dijkstra_snap(te, mg, hm, data, occupied_halfs, i0, i1)
            });
        }

        ehalfs.resize(tm.thalfs.size());
        equads.resize(tm.tquads.size());

        for (size_t i = 0; i < tm.thalfs.size(); ++i) {
            const auto& [_, id, twid, teid, nxt_id, prv_id, cano] = tm.thalfs[i];
            const auto& ee = eedges[teid];
            vec<Half> hs = {};

            if (cano) { for (Half h: ee.halfs)               { hs.push_back(h); } }
            else      { for (Half h: ee.halfs | vw::reverse) { hs.push_back(h.twin()); } }

            ehalfs[i] = {
                .em     = this,
                .id     = id,
                .twid   = twid,
                .eeid   = teid,
                .eqid   = tm.th2quad[id],
                .nxt_id = nxt_id,
                .prv_id = prv_id,
                .cano   = cano,
                .bgn    = cano && mg.mnodes[ee.fr_nid].jt == mc::JunctionType::F,
                .end    = cano && mg.mnodes[ee.to_nid].jt == mc::JunctionType::T,
                .x      = X[teid],
                .halfs  = hs
            };
        }

        for (size_t i = 0; i < tm.tquads.size(); ++i) {
            const auto& [id, thids, sides] = tm.tquads[i];
            vec<Edata> data_;
            for (int j = 0; j < thids.size(); ++j) { data_.push_back({thids[j], sides[j]}); }
            equads[i] = { .id = id, .data = data_ };
        }

        for (auto ee: eedges) {
            if (mg.mnodes[ee.fr_nid].jt == mc::JunctionType::F) sings.insert(ee.halfs.front().tail().id);
            if (mg.mnodes[ee.to_nid].jt == mc::JunctionType::F) sings.insert(ee.halfs.back().head().id);
        }
    }

    bool collapse_equad(int eqid);
    bool collapse_ehalf(int ehid);

    void replace_remain_side(
        int eqid,
        int side,
        const vec<int> &ehids_path,
        const vec<int> &ehids_prev,
        const vec<int> &ehids_next
    );
    void replace_ehalf(
        int eqid,
        int ehid,
        const vec<int>& reps, // ehids for replacing ehalf above
        const vec<int>& ext0, // ehids for extending next side of ehalfs
        const vec<int>& ext1  // ehids for extending prev side of ehalfs
    );

    static vec<Half> compute_dijkstra_snap(
        const     Tedge& te,
        const mc::Mgrph& mg,
        const     Hmesh& hm,
        const TrackedDenseMesh& data,
        vec<int>& occupied_halfs,
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
                if (rg::contains(occupied_halfs, h.id)) continue;
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

                auto uv0 = mc::get_face_uv(mg.mnodes[sg->fr_nid], sg->face_id, mg.hm, mg.cf);
                auto uv1 = mc::get_face_uv(mg.mnodes[sg->to_nid], sg->face_id, mg.hm, mg.cf);
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

        if (!reached) { std::cout << "Dijkstra failed!" << std::endl; return {}; }

        vec<Half> path;
        int curr = end_vid;
        while (curr != bgn_vid) {
            int hid = came_from_hid[curr];
            Half h_ = hm.halfs[hid];
            path.push_back(h_);
            occupied_halfs.push_back(hid);
            occupied_halfs.push_back(h_.twin().id);
            curr = came_from_vid[curr];
        }
        rg::reverse(path);
        return path;
    }

    static double compute_point_to_segment_distance(complex p, complex a, complex b) {
        complex ab = b - a;
        complex ap = p - a;
        if (std::abs(ab) < 1e-8) return std::abs(ap);
        auto t = (ap.real() * ab.real() + ap.imag() * ab.imag()) / std::norm(ab);
        auto proj = a + std::clamp(t, 0., 1.) * ab;
        return std::abs(p - proj);
    }
};


inline const Eedge& Ehalf::edge() const { return em->eedges[eeid]; }
inline const Ehalf& Ehalf::twin() const { return em->ehalfs[twid]; }
inline const Ehalf& Ehalf::next() const { return em->ehalfs[nxt_id]; }
inline const Ehalf& Ehalf::prev() const { return em->ehalfs[prv_id]; }
inline Vert Ehalf::tail() const { return halfs.front().tail(); }
inline Vert Ehalf::head() const { return halfs.back().head();  }
}
#endif