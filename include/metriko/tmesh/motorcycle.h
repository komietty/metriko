//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_MOTORCYCLE_H
#define METRIKO_MOTORCYCLE_H
#include "../common/utilities.h"
#include "../common/predicates.h"
#include "../hmesh/hmesh.h"
#include "metriko/hmesh/hmloc.h"
#include "metriko/hmesh/hmray.h"
#include "metriko/hmesh/utilities.h"

namespace metriko {
constexpr double TOLERANCE_HALF = 1e-6;    //
constexpr double TOLERANCE_CRNR = 2e-6;    //
constexpr double TOLERANCE_EDGE_AB = 1e-6; // tolerance on two curvs crash close to an edge
constexpr double TOLERANCE_EDGE_CD = 1e-6; // tolerance on two curvs crash close to an edge

struct Mgrph;

struct Mport {
    complex uv;
    complex dr;
    int crnr_id = -1;
    int this_id = -1;
    int next_id = -1;
    int prev_id = -1;
};

enum class JunctionType { None, F, T, C };

struct Mnode {
    JunctionType jt = JunctionType::None;
    HmLoc      loc = {};
    vec<Row2i> adj = {}; // x: curv_id, y: sgmt_id
};

struct Msgmt {
    int fr_nid  = -1;
    int to_nid  = -1;
    int this_id = -1;
    int curv_id = -1;
    int face_id = -1;
};

struct Mbuff {
    complex uv = {0, 0};
    complex dr = {0, 0};
    int  cid = -1; // crnr_id; better using variant
    int  hid = -1; // half_id; better using variant
    double r = -1; // ratio.
    bool bgn = false;
    bool end = false;
};

struct Mcurv {
    Mgrph *mg = nullptr;
    int id = -1;
    vec<Msgmt> sgmts = {};
    Mbuff buff = {};

    bool operator==(const Mcurv &c) const { return id == c.id; }
    void add_segment(const Hmesh& hm, const VecXc& cf);
    int resolve_bgn_node(const Hmesh& hm, bool bgn, int cid) const;
};

inline void update_to_twin(const Hmesh& hm, const VecXc& cf, const VecXi& matching, Mbuff& buff) {
    auto get_m = [&](Half h) { return (h.isCanonical() ? -1 : 1) * matching[h.edge().id]; };

    if (buff.cid != -1) {
        auto c = hm.crnrs[buff.cid];
        auto v = c.vert();
        auto d = buff.dr;
        for (Half h: v.adjHalfs(c.half().next().twin())) {
            d *= std::polar(1., PI / 2 * get_m(h));
            auto c1  = h.next().crnr();
            auto uv0 = cf(c1.id);
            auto uv1 = cf(c1.half().crnr_t().id);
            auto uv2 = cf(c1.half().crnr_h().id);
            if (is_points_into(uv0, uv1, uv2, uv0 + d, 0) && c1 != c) { buff = Mbuff{.uv = uv0, .dr = d, .cid = c1.id}; return; }
        }
    }
    if (buff.hid != -1) {
        auto h = hm.halfs[buff.hid].twin();
        auto r = 1. - buff.r;
        auto uv = lerp(cf(h.next().crnr().id), cf(h.prev().crnr().id), r);
        auto dr = std::polar(1., PI / 2 * get_m(h)) * buff.dr;
        buff = Mbuff{.uv = uv, .dr = dr, .hid = h.id, .r = r};
        return;
    }

    METRIKO_FAIL("not implemented");
}

inline void update_to_oppo(const Hmesh& hm, const VecXc& cf, Mbuff& buff) {
    auto fr = buff.cid != -1 ? HmLoc{HmLocOnC{buff.cid}} : HmLoc{HmLocOnH{buff.hid, buff.r}};
    auto it = find_ray_intersection(hm, fr, cf, buff.dr);
    std::visit(overloaded{
        [&](const HmLocOnC& l) { buff = Mbuff{.uv = cf(l.id), .dr = buff.dr, .cid = l.id}; },
        [&](const HmLocOnH& l) { buff = Mbuff{.uv = lerp(cf(hm.halfs[l.id].crnr_t().id), cf(hm.halfs[l.id].crnr_h().id), l.r), .dr = buff.dr, .hid = l.id, .r = l.r}; },
        [&](const auto&) { METRIKO_FAIL("unexpected exit location"); },
    }, it);
}

struct  Mgrph {
    const Hmesh &hm;
    const VecXc &cf;
    vec<Mport> mports;
    vec<Mnode> mnodes;
    vec<Mcurv> mcurvs;

    Mgrph(
        const Hmesh &hm,
        const VecXc &cf,
        const VecXi &matching,
        const VecXi &singular
    ) : hm(hm), cf(cf) {
        gen_ports(singular);

        for (auto v: hm.verts | vw::filter([&](auto& v) { return singular[v.id]; }))
            mnodes.push_back({.jt = JunctionType::F, .loc = HmLocOnV{v.id}});

        // 1: Add the first segment for each curve
        mcurvs.reserve(mports.size());
        for (const auto& p : mports) {
            mcurvs.push_back({.mg = this, .id = p.this_id, .buff = {.uv = p.uv, .dr = p.dr, .cid = p.crnr_id, .bgn = true}});
            mcurvs.back().add_segment(hm, cf);
        }

        // 2: Add further segments until every curve crash to another curve
        while (rg::any_of(mcurvs, [](auto &c) { return !c.buff.end; })) {
            for (auto &mc: mcurvs) {
                if (mc.buff.end) continue;
                update_to_twin(hm, cf, matching, mc.buff);
                mc.add_segment(hm, cf);
            }
        }

        collect_node_adjacency();
        sort_node_adjacency();

        for (auto &c: mcurvs) {
        for (int i = 0; i < c.sgmts.size(); ++i) {
            c.sgmts[i].curv_id = c.id;
            c.sgmts[i].this_id  = i;
        }}
    }

    void gen_ports(const VecXi& singular);
    void collect_node_adjacency();
    void sort_node_adjacency();
};
}
#endif
