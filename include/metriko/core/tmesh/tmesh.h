//
//--- Copyright (C) 2026 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
#ifndef METRIKO_TMESH_H
#define METRIKO_TMESH_H

#include "motorcycle.h"
#include <vector>
#include <cassert>
#include <complex>

namespace metriko {

struct Tmesh;

struct Tedge {
    int id;
    int fr_nid;
    int to_nid;
    int crv_id;
    bool isBgn;
    bool isEnd;
    double len = 0;
    vec<mc::Msgmt> segs;
};

struct Thalf {
    const Tmesh* tm = nullptr;
    int id     = -1;
    int twid   = -1;
    int teid   = -1;
    int nxt_id = -1;
    int prv_id = -1;
    bool cano  = false;

    const Thalf& twin() const;
    const Thalf& next() const;
    const Thalf& prev() const;
    const Tedge& edge() const;
    int nid_fr() const { return cano ? edge().fr_nid : edge().to_nid; }
    int nid_to() const { return cano ? edge().to_nid : edge().fr_nid; }
};

struct Tquad {
    int id = -1;
    vec<int> thids;
    vec<int> sides;

    vec<int> thids_by_side(int side) const {
        return vw::zip(thids, sides)
            | vw::filter([&](const auto& p) { return std::get<1>(p) == side; })
            | vw::elements<0>
            | rg::to<vec<int>>();
    }
};

struct Tmesh {
    vec<Tquad> tquads;
    vec<Thalf> thalfs;
    vec<Tedge> tedges;
    VecXi th2quad;
    VecXi th2side;
    VecXi th2iter;
    size_t nTQ;
    size_t nTE;
    size_t nTH;

    explicit Tmesh(const mc::Mgrph& mg) {
        //===== 1. Extract Thalfs and Tedges =====
        for (const auto& mc: mg.mcurvs) {
            vec<mc::Msgmt> sgs;
            int bgn_nid = mc.sgmts.front().fr_nid;

            for (auto& sg: mc.sgmts) {
                sgs.push_back(sg);

                if (mg.mnodes[sg.to_nid].jt != mc::JunctionType::None) {
                    int teid = tedges.size();
                    int thid = thalfs.size();
                    int end_nid = sg.to_nid;

                    double len = 0;
                    for (const auto& s: sgs) {
                        auto fr = mc::get_face_uv(mg.mnodes[s.fr_nid], s.face_id, mg.hm, mg.cf);
                        auto to = mc::get_face_uv(mg.mnodes[s.to_nid], s.face_id, mg.hm, mg.cf);
                        len += std::abs(to - fr);
                    }

                    thalfs.push_back({.tm = this, .id = thid,     .twid = thid + 1, .teid = teid, .cano = true });
                    thalfs.push_back({.tm = this, .id = thid + 1, .twid = thid,     .teid = teid, .cano = false});
                    tedges.push_back({
                        .id     = teid,
                        .fr_nid = bgn_nid,
                        .to_nid = end_nid,
                        .crv_id = mc.id,
                        .isBgn  = mg.mnodes[bgn_nid].jt == mc::JunctionType::F,
                        .isEnd  = mg.mnodes[end_nid].jt == mc::JunctionType::T,
                        .len    = len,
                        .segs   = std::move(sgs)
                    });

                    bgn_nid = end_nid;
                    sgs.clear();
                }
            }
        }

        // ===== 2. Determine Thalfs adjacency =====
        for (int nid = 0; nid < mg.mnodes.size(); ++nid) {
            const auto& mn = mg.mnodes[nid];
            if (mn.jt == mc::JunctionType::None) continue;

            vec<Thalf*> outgoing;
            for (auto& th: thalfs) if (th.nid_fr() == nid) outgoing.push_back(&th);
            if (outgoing.empty()) continue;

            auto get_rank = [&](const Thalf* th) {
                auto& sg = th->cano ? th->edge().segs.front() : th->edge().segs.back();
                auto  it = rg::find_if(mn.adj, [&](auto& as) { return as.curv_id == sg.curv_id && as.sgmt_id == sg.this_id; });
                return std::distance(mn.adj.begin(), it);
            };

            rg::sort(outgoing, [&](const Thalf* a, const Thalf* b) { return get_rank(a) < get_rank(b); });

            int n = outgoing.size();
            for (int i = 0; i < n; ++i) {
                Thalf* th_out  = outgoing[i];
                Thalf* th_in   = &thalfs[th_out->twid];
                Thalf* th_next = outgoing[(i - 1 + n) % n];
                th_in->nxt_id = th_next->id;
                th_next->prv_id = th_in->id;
            }
        }

        // ====== 3. Assign sides and thids =====
        vec visited(thalfs.size(), false);

        for (int i = 0; i < thalfs.size(); ++i) {
            if (visited[i]) continue;

            Tquad tq;
            tq.id = tquads.size();
            int curr_thid = i;
            int curr_side = 0;

            do {
                if (visited[curr_thid]) break;
                tq.thids.push_back(curr_thid);
                tq.sides.push_back(curr_side);
                visited[curr_thid] = true;
                auto& curr_th = thalfs[curr_thid];
                auto& next_th = thalfs[curr_th.nxt_id];
                if (curr_th.edge().crv_id != next_th.edge().crv_id) curr_side = (curr_side + 1) % 4;
                curr_thid = next_th.id;
            } while (curr_thid != i);

            // Rotate until sides data is sequential
            auto& fst_th = thalfs[tq.thids.front()];
            auto& lst_th = thalfs[tq.thids.back()];
            if (fst_th.edge().crv_id == lst_th.edge().crv_id) {
                int f = tq.sides.front();
                int n = rg::distance(tq.sides | vw::take_while([=](int x) { return x == f; }));
                rg::rotate(tq.sides, tq.sides.begin() + n);
                rg::rotate(tq.thids, tq.thids.begin() + n);

                int s = 0;
                int last_c = thalfs[tq.thids[0]].edge().crv_id;
                for (int j = 0; j < tq.thids.size(); ++j) {
                    int this_c = thalfs[tq.thids[j]].edge().crv_id;
                    if (this_c != last_c) s = (s + 1) % 4;
                    tq.sides[j] = s;
                    last_c = this_c;
                }
            }
            tquads.push_back(tq);
        }

        th2quad.resize(thalfs.size());
        th2side.resize(thalfs.size());
        th2iter.resize(thalfs.size());
        for (auto& [id, thids, sides] : tquads) {
        for (int j = 0; j < thids.size(); ++j) {
            th2quad[thids[j]] = id;
            th2side[thids[j]] = sides[j];
            th2iter[thids[j]] = j;
        }}

        nTE = tedges.size();
        nTH = thalfs.size();
        nTQ = tquads.size();
    }
};

inline const Tedge& Thalf::edge() const { return tm->tedges[teid]; }
inline const Thalf& Thalf::twin() const { return tm->thalfs[twid]; }
inline const Thalf& Thalf::next() const { return tm->thalfs[nxt_id]; }
inline const Thalf& Thalf::prev() const { return tm->thalfs[prv_id]; }
}
#endif