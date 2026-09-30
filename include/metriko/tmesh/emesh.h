//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_EMESH_H
#define METRIKO_EMESH_H
#include "motorcycle.h"
#include "metriko/hmesh/hmloc.h"
#include "metriko/hmesh/hpath.h"
#include "metriko/hmesh/utilities.h"

namespace metriko {
struct Emesh;

struct Tqpoint {
    int    nid = -1; // tnode of this point
    int    val = -1;
    double ord = -1; // geometric sort key: arc-length (th.r) position along the spine, per-tquad normalized
    bool   top = false;
};

struct Tqchain {
    vec<Tqpoint> pts = {};
    vec<int> tqids   = {};
    vec<int> bounds  = {};
    vec<int> thids_z = {}; // zero length thids from left to right
    vec<int> thids_t = {};
    vec<int> thids_b = {};
    int      thid_r  = -1;
    int      thid_l  = -1;
};

struct Eedge {
    int id = -1;
    vec<int> nids = {};
    void insert_locs(const vec<int>& locs);
};

struct Ehalf {
    const Emesh* tm = nullptr;
    int id   = -1;
    int twid = -1;
    int teid = -1;
    int tqid = -1;
    bool cano = false;
    bool bgn  = false;
    bool end  = false;
    double x  = -1;
    double r  = -1;
    int nid_fr() const;
    int nid_to() const;
    const HmLoc& loc_fr() const;
    const HmLoc& loc_to() const;
};

struct Edata {
    int thid = -1;
    int side = -1;
};

struct Equad {
    int id = -1;
    vec<Edata> data;

    auto find_of(this auto& self, int thid) { return rg::find(self.data, thid, &Edata::thid); }
    int  side_of(const Ehalf& th) const { return rg::find(data, th.id, &Edata::thid)->side; }

    // every side keeps at least one thalf: a tquad that lost a side is no longer a quad
    bool is_valid() const {
        for (int side = 0; side < 4; ++side) if (rg::none_of(data, [&](const Edata& d) { return d.side == side; })) return false;
        return true;
    }

    vec<int> thids(int side) const {
        return data | vw::filter([&](auto& d) { return d.side == side; })
                    | vw::transform([](auto& d) { return d.thid; })
                    | rg::to<vec<int>>();
    }
};

struct Emesh {
    const Hmesh& hm;
    vec<HmLoc> tnodes = {};
    vec<Eedge> tedges = {};
    vec<Ehalf> thalfs = {};
    vec<Equad> tquads = {};

    Emesh(const Emesh&) = delete;
    Emesh(Emesh&&)      = delete;

    explicit Emesh(const Hmesh& hm): hm(hm) {}

    // builds the t-mesh from the motorcycle graph: every curve is split into tedges at its junction nodes, the
    // thalfs leaving a junction are ordered by the node's adjacency, and every tquad is traced along the next
    // thalfs. r is the length of the tedge in the parameter domain (used by the quantization); x stays -1 until
    // set_x() is called with the quantization result
    explicit Emesh(const Mgrph& mg): hm(mg.hm) {
        tnodes.reserve(mg.mnodes.size());
        for (const Mnode& mn : mg.mnodes) {
            tnodes.push_back(std::visit(overloaded{
                [&](const auto&     _) -> HmLoc { METRIKO_FAIL("no impl"); },
                [&](const HmLocOnV& v) -> HmLoc { return HmLocOnV{v.id}; },
                [&](const HmLocOnE& e) -> HmLoc { return HmLocOnE{.id = e.id, .r = e.r}; },
                [&](const HmLocOnP& l) -> HmLoc { Face f = hm.faces[l.id]; return HmLocOnF{.id = l.id, .xy = f.to_local(f.uv2pos(l.uv))}; },
            }, mn.loc));
        }

        // 1: one tedge per run of segments between two junction nodes. a trailing run that ends at no junction
        //    does not make a tedge
        vec<int> crv;             // teid -> curve it lies on, delimits the tquad sides below
        vec<const Msgmt*> sg_bgn; // teid -> its first / last segment, ranks the thalfs at the junctions
        vec<const Msgmt*> sg_end;
        for (const Mcurv& mc: mg.mcurvs) {
            int    bgn_nid = mc.sgmts.front().fr_nid;
            size_t bgn_sg  = 0;
            double len     = 0;
            vec<int> nids  = {bgn_nid};
            for (size_t i = 0; i < mc.sgmts.size(); ++i) {
                const Msgmt& sg = mc.sgmts[i];
                len += std::abs(get_face_uv(mg.mnodes[sg.to_nid].loc, sg.face_id, hm) - get_face_uv(mg.mnodes[sg.fr_nid].loc, sg.face_id, hm));
                nids.push_back(sg.to_nid);
                if (mg.mnodes[sg.to_nid].jt == JunctionType::None) continue;

                const int  teid = tedges.size();
                const int  thid = thalfs.size();
                const bool bgn  = mg.mnodes[bgn_nid].jt == JunctionType::F;
                const bool end  = mg.mnodes[sg.to_nid].jt == JunctionType::T;
                tedges.push_back({.id = teid, .nids = std::move(nids)});
                thalfs.push_back({.tm = this, .id = thid,     .twid = thid + 1, .teid = teid, .cano = true,  .bgn = bgn, .end = end, .r = len});
                thalfs.push_back({.tm = this, .id = thid + 1, .twid = thid,     .teid = teid, .cano = false, .r = len});
                crv.push_back(mc.id);
                sg_bgn.push_back(&mc.sgmts[bgn_sg]);
                sg_end.push_back(&sg);

                bgn_nid = sg.to_nid;
                bgn_sg  = i + 1;
                len     = 0;
                nids    = {bgn_nid};
            }
        }

        // 2: at every junction, the thalf entering it continues with the next outgoing one clockwise
        vec<int>      nxid(thalfs.size(), -1);
        vec<vec<int>> outgoing(mg.mnodes.size());
        for (const Ehalf& th: thalfs) outgoing[th.nid_fr()].push_back(th.id);
        for (int nid = 0; nid < mg.mnodes.size(); ++nid) {
            const Mnode& mn  = mg.mnodes[nid];
            vec<int>&    out = outgoing[nid];
            if (mn.jt == JunctionType::None || out.empty()) continue;
            auto rank = [&](int thid) {
                const Ehalf& th = thalfs[thid];
                const Msgmt& sg = th.cano ? *sg_bgn[th.teid] : *sg_end[th.teid];
                return rg::distance(mn.adj.begin(), rg::find_if(mn.adj, [&](const Row2i& ad) { return ad.x() == sg.curv_id && ad.y() == sg.this_id; }));
            };
            rg::sort(out, {}, rank);
            const int n = out.size();
            for (int i = 0; i < n; ++i) nxid[thalfs[out[i]].twid] = out[(i - 1 + n) % n];
        }

        // 3: trace every tquad along nxid. a side is a run of thalfs on one curve
        vec visited(thalfs.size(), false);
        for (int i = 0; i < thalfs.size(); ++i) {
            if (visited[i]) continue;
            Equad tq;
            tq.id = tquads.size();
            int thid = i;
            int side = 0;
            do {
                if (visited[thid]) break;
                tq.data.emplace_back(thid, side);
                visited[thid] = true;
                const int next = nxid[thid];
                if (crv[thalfs[thid].teid] != crv[thalfs[next].teid]) side = (side + 1) % 4;
                thid = next;
            } while (thid != i);

            // the walk started inside a side: rotate so that every side is one contiguous run, then renumber
            if (crv[thalfs[tq.data.front().thid].teid] == crv[thalfs[tq.data.back().thid].teid]) {
                const int f = tq.data.front().side;
                const int n = rg::distance(tq.data | vw::take_while([=](const Edata& d) { return d.side == f; }));
                rg::rotate(tq.data, tq.data.begin() + n);
                int s = 0;
                int c = crv[thalfs[tq.data.front().thid].teid];
                for (Edata& d: tq.data) {
                    if (crv[thalfs[d.thid].teid] != c) s = (s + 1) % 4;
                    d.side = s;
                    c = crv[thalfs[d.thid].teid];
                }
            }
            for (const Edata& d: tq.data) thalfs[d.thid].tqid = tq.id;
            tquads.push_back(std::move(tq));
        }
    }

    void set_x(const VecXd& X) { for (Ehalf& th: thalfs) if (th.id != -1) th.x = X[th.teid]; }

    int step_next(int thid) const { auto& [_, d] = tquads[thalfs[thid].tqid]; auto it = rg::find(d, thid, &Edata::thid); METRIKO_CHECK(it != d.end(), "step next failed"); return circular_next(d, it)->thid; };
    int step_prev(int thid) const { auto& [_, d] = tquads[thalfs[thid].tqid]; auto it = rg::find(d, thid, &Edata::thid); METRIKO_CHECK(it != d.end(), "step prev failed"); return circular_prev(d, it)->thid; };
    int count_adj_tquads(int thid0) const {
        int count = 0, thid = thid0;
        do { ++count; thid = step_next(thalfs[thid].twid); }
        while (thid != thid0 && count <= thalfs.size());
        return count;
    }

    vec<Erng> allowed_range_thalfs(const vec<int>& thids) const;
    vec<Erng> allowed_range_tquads(const vec<int>& tqids) const;

    auto live_tedges() const { return tedges | vw::filter([](const Eedge& te) { return te.id != -1; }); }
    auto live_tedges()       { return tedges | vw::filter([](      Eedge& te) { return te.id != -1; }); }
    auto live_tquads() const { return tquads | vw::filter([](const Equad& tq) { return tq.id != -1; }); }
    auto live_tquads()       { return tquads | vw::filter([](      Equad& tq) { return tq.id != -1; }); }

    void collapse_tedge_snap(bool flag);
    void collapse_tedge_snap_dedup(int teid);
    bool reroute_tedge(int teid);

    void collapse_thalf(int thid);
    bool collapse_point_tquad(int tqid);
    bool collapse_tquad_chain_prepare(int tqid, Tqchain& chain) const;
    void collapse_tquad_chain_execute(Tqchain& chain);

    vec<int> add_new_path(const vec<HmLoc>& path, int nid0, int nid1) {
        vec<int> nids;
        for (int i = 0; i < path.size(); ++i) {
            if      (i == 0)                  nids.push_back(nid0);
            else if (i + 1 == path.size())    nids.push_back(nid1);
            else { tnodes.push_back(path[i]); nids.push_back(tnodes.size() - 1); }
        }
        return nids;
    }

    double path_length(const vec<int>& nids) const {
        double r = 0;
        for (size_t k = 0; k + 1 < nids.size(); ++k)
            r += (get_ptloc_pos(hm, tnodes[nids[k + 1]]) - get_ptloc_pos(hm, tnodes[nids[k]])).norm();
        return r;
    }
};

inline int Ehalf::nid_fr() const { const auto& nids = tm->tedges[teid].nids; return cano ? nids.front() : nids.back(); }
inline int Ehalf::nid_to() const { const auto& nids = tm->tedges[teid].nids; return cano ? nids.back() : nids.front(); }
inline const HmLoc& Ehalf::loc_fr() const { return tm->tnodes[nid_fr()]; }
inline const HmLoc& Ehalf::loc_to() const { return tm->tnodes[nid_to()]; }

inline void Eedge::insert_locs(const vec<int>& locs) {
    int f = locs.front();
    int b = locs.back();
    if      (nids.front() == b) { nids.insert(nids.begin(), locs.begin(), locs.end() - 1); } // prepend [f..b-1]
    else if (nids.back()  == f) { nids.insert(nids.end(),   locs.begin() + 1, locs.end()); } // append  [f+1..b]
    else if (nids.front() == f) { vec<int> s(locs.begin() + 1, locs.end()); rg::reverse(s); nids.insert(nids.begin(), s.begin(), s.end()); } // prepend reverse([f+1..b])
    else if (nids.back()  == b) { vec<int> s(locs.begin(), locs.end() - 1); rg::reverse(s); nids.insert(nids.end(),   s.begin(), s.end()); } // append  reverse([f..b-1])
    else METRIKO_FAIL("merge_into: the chains share no end node");
}
}
#endif
