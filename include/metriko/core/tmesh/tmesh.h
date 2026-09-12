//
//--- Copyright (C) 2026 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
#ifndef METRIKO_TMESH_H
#define METRIKO_TMESH_H
#include "motorcycle.h"

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
    vec<Msgmt> segs;
};

struct Thalf {
    const Tmesh* tm = nullptr;
    int id   = -1;
    int twid = -1;
    int teid = -1;
    int nxid = -1;
    int pvid = -1;
    bool cano  = false;

    const Thalf& twin() const;
    const Thalf& next() const;
    const Thalf& prev() const;
    const Tedge& edge() const;
    int nid_fr() const { return cano ? edge().fr_nid : edge().to_nid; }
    int nid_to() const { return cano ? edge().to_nid : edge().fr_nid; }
};

struct Tdata {
    int thid;
    int side;
};

struct Tquad {
    int id = -1;
    vec<Tdata> data;

    vec<int> thids_by_side(int side) const {
        return data
            | vw::filter([&](const Tdata& d) { return d.side == side; })
            | vw::transform([](const Tdata& d) { return d.thid; })
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

    explicit Tmesh(const Mgrph& mg) {
        //===== 1. Extract Thalfs and Tedges =====
        for (const auto& mc: mg.mcurvs) {
            vec<Msgmt> sgs;
            int bgn_nid = mc.sgmts.front().fr_nid;

            for (auto& sg: mc.sgmts) {
                sgs.push_back(sg);

                if (mg.mnodes[sg.to_nid].jt != JunctionType::None) {
                    int teid = tedges.size();
                    int thid = thalfs.size();
                    int end_nid = sg.to_nid;

                    double len = 0;
                    for (const auto& s: sgs) {
                        auto fr = get_face_uv(mg.mnodes[s.fr_nid], s.face_id, mg.hm, mg.cf);
                        auto to = get_face_uv(mg.mnodes[s.to_nid], s.face_id, mg.hm, mg.cf);
                        len += std::abs(to - fr);
                    }

                    thalfs.push_back({.tm = this, .id = thid,     .twid = thid + 1, .teid = teid, .cano = true });
                    thalfs.push_back({.tm = this, .id = thid + 1, .twid = thid,     .teid = teid, .cano = false});
                    tedges.push_back({
                        .id     = teid,
                        .fr_nid = bgn_nid,
                        .to_nid = end_nid,
                        .crv_id = mc.id,
                        .isBgn  = mg.mnodes[bgn_nid].jt == JunctionType::F,
                        .isEnd  = mg.mnodes[end_nid].jt == JunctionType::T,
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
            if (mn.jt == JunctionType::None) continue;

            vec<Thalf*> outgoing;
            for (auto& th: thalfs) if (th.nid_fr() == nid) outgoing.push_back(&th);
            if (outgoing.empty()) continue;

            auto get_rank = [&](const Thalf* th) {
                auto& sg = th->cano ? th->edge().segs.front() : th->edge().segs.back();
                auto  it = rg::find_if(mn.adj, [&](const Row2i& ad) { return ad.x() == sg.curv_id && ad.y() == sg.this_id; });
                return std::distance(mn.adj.begin(), it);
            };

            rg::sort(outgoing, [&](const Thalf* a, const Thalf* b) { return get_rank(a) < get_rank(b); });

            int n = outgoing.size();
            for (int i = 0; i < n; ++i) {
                Thalf* th_out  = outgoing[i];
                Thalf* th_in   = &thalfs[th_out->twid];
                Thalf* th_next = outgoing[(i - 1 + n) % n];
                th_in->nxid = th_next->id;
                th_next->pvid = th_in->id;
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
                tq.data.push_back({curr_thid, curr_side});
                visited[curr_thid] = true;
                auto& curr_th = thalfs[curr_thid];
                auto& next_th = thalfs[curr_th.nxid];
                if (curr_th.edge().crv_id != next_th.edge().crv_id) curr_side = (curr_side + 1) % 4;
                curr_thid = next_th.id;
            } while (curr_thid != i);

            // Rotate until sides data is sequential
            auto& fst_th = thalfs[tq.data.front().thid];
            auto& lst_th = thalfs[tq.data.back().thid];
            if (fst_th.edge().crv_id == lst_th.edge().crv_id) {
                int f = tq.data.front().side;
                int n = rg::distance(tq.data | vw::take_while([=](const Tdata& d) { return d.side == f; }));
                rg::rotate(tq.data, tq.data.begin() + n);

                int s = 0;
                int last_c = thalfs[tq.data[0].thid].edge().crv_id;
                for (Tdata& td: tq.data) {
                    int this_c = thalfs[td.thid].edge().crv_id;
                    if (this_c != last_c) s = (s + 1) % 4;
                    td.side = s;
                    last_c = this_c;
                }
            }
            tquads.push_back(tq);
        }

        th2quad.resize(thalfs.size());
        th2side.resize(thalfs.size());
        th2iter.resize(thalfs.size());
        for (auto& tq : tquads) {
        for (int j = 0; j < tq.data.size(); ++j) {
            th2quad[tq.data[j].thid] = tq.id;
            th2side[tq.data[j].thid] = tq.data[j].side;
            th2iter[tq.data[j].thid] = j;
        }}

        nTE = tedges.size();
        nTH = thalfs.size();
        nTQ = tquads.size();
    }

    bool check_non_zero_tquad(const VecXd& X) const {
        for (const auto& [id, data] : tquads) {
            if (rg::all_of(data, [&](const Tdata& d) { return X[thalfs[d.thid].teid] == 0; })) return false;
        }
        return true;
    }
};

inline const Tedge& Thalf::edge() const { return tm->tedges[teid]; }
inline const Thalf& Thalf::twin() const { return tm->thalfs[twid]; }
inline const Thalf& Thalf::next() const { return tm->thalfs[nxid]; }
inline const Thalf& Thalf::prev() const { return tm->thalfs[pvid]; }
}
#endif