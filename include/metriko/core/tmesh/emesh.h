//
//--- Copyright (C) 2026 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
#ifndef METRIKO_EMESH_H
#define METRIKO_EMESH_H

#include "tmesh.h"
#include "emesh_subdivide.h"
#include <queue>
#include <map>
#include <set>
#include <iostream>

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

        // assign already known info about eedge
        for (auto& te: tm.tedges) {
            if (mg.mnodes[te.fr_nid].jt == mc::JunctionType::F) sings.insert(mnode2dense_v.at(te.fr_nid));
            if (mg.mnodes[te.to_nid].jt == mc::JunctionType::F) sings.insert(mnode2dense_v.at(te.to_nid));
            eedges.push_back({
                .id = te.id,
                .fr_nid = te.fr_nid,
                .to_nid = te.to_nid,
                .len = te.len,
                .halfs = {}
            });
        }

        assign_first_half(tm, mg, mnode2dense_v, data);
        assign_inter_half(tm, mg, mnode2dense_v, data);
        //assign_last_half(tm, mg, mnode2dense_v);

        //vec<bool> occupied_verts = vec(hm.verts.size(), false);
        //for (const auto& te: tm.tedges) {
        //    auto i0 = mnode2dense_v.at(te.fr_nid);
        //    auto i1 = mnode2dense_v.at(te.to_nid);
        //    eedges.push_back({
        //        .id     = te.id,
        //        .fr_nid = te.fr_nid,
        //        .to_nid = te.to_nid,
        //        .len    = te.len,
        //        .halfs  = compute_dijkstra_snap(te, mg, hm, data, occupied_verts, i0, i1)
        //    });
        //}

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

        //for (auto ee: eedges) {
        //    if (mg.mnodes[ee.fr_nid].jt == mc::JunctionType::F) sings.insert(ee.halfs.front().tail().id);
        //    if (mg.mnodes[ee.to_nid].jt == mc::JunctionType::F) sings.insert(ee.halfs.back().head().id);
        //}
    }

    bool collapse_equad(int eqid);
    bool collapse_ehalf(int ehid);

    void assign_first_half(
        const Tmesh& tm,
        const mc::Mgrph& mg,
        const std::map<int, int>& mnode2dense_v,
        const TrackedDenseMesh& data
    );

    void assign_inter_half(
        const Tmesh& tm,
        const mc::Mgrph& mg,
        const std::map<int, int>& mnode2dense_v,
        const TrackedDenseMesh& data
    );

    void assign_last_half (
        const Tmesh& tm,
        const mc::Mgrph& mg,
        const std::map<int, int>& mnode2dense_v
    );

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

};


inline const Eedge& Ehalf::edge() const { return em->eedges[eeid]; }
inline const Ehalf& Ehalf::twin() const { return em->ehalfs[twid]; }
inline const Ehalf& Ehalf::next() const { return em->ehalfs[nxt_id]; }
inline const Ehalf& Ehalf::prev() const { return em->ehalfs[prv_id]; }
inline Vert Ehalf::tail() const { return halfs.front().tail(); }
inline Vert Ehalf::head() const { return halfs.back().head();  }
}
#endif