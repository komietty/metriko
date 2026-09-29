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
constexpr double TOLERANCE_EDGE_AB = 1e-6; // tolerance on two curvs crash close to an edge
constexpr double TOLERANCE_EDGE_CD = 1e-6; // tolerance on two curvs crash close to an edge

struct Mgrph;

struct Mport {
    complex dir;
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
    HmLoc   loc = {}; // curve front, chart-aware (OnC or OnH)
    complex dir = {}; // direction in the chart of get_chart_face(loc)
    bool end = false;
};

struct Mcurv {
    Mgrph *mg = nullptr;
    int id = -1;
    vec<Msgmt> sgmts = {};
    Mbuff buff = {};

    bool operator==(const Mcurv &c) const { return id == c.id; }
    void add_segment(const Hmesh& hm, const VecXc& cf);
    int resolve_fr_node(const Hmesh& hm, const HmLoc& fr) const;
};

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
            mcurvs.push_back({.mg = this, .id = p.this_id, .buff = {.loc = HmLocOnC{p.crnr_id}, .dir = p.dir}});
            mcurvs.back().add_segment(hm, cf);
        }

        // 2: Add further segments until every curve crash to another curve
        while (rg::any_of(mcurvs, [](auto &c) { return !c.buff.end; })) {
            for (auto &mc: mcurvs) {
                if (mc.buff.end) continue;
                std::tie(mc.buff.loc, mc.buff.dir) =
                    cross_to_twin(hm, mc.buff.loc, cf, matching, mc.buff.dir);
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
