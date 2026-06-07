//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_MOTORCYCLE_H
#define METRIKO_MOTORCYCLE_H
#include "../common/utilities.h"
#include "../common/predicates.h"
#include "../hmesh/hmesh.h"
#include "./common.h"
#include "metriko/core/hmesh/hmloc.h"

namespace metriko::mc {
constexpr double TOLERANCE_HALF = 1e-6;    //
constexpr double TOLERANCE_CRNR = 2e-6;    //
constexpr double TOLERANCE_EDGE_AB = 1e-4; // tolerance on two curvs crash close to an edge
constexpr double TOLERANCE_EDGE_CD = 1e-2; // tolerance on two curvs crash close to an edge


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

struct Asgmt {
    int curv_id;
    int sgmt_id;
};

struct Mnode {
    JunctionType jt = JunctionType::None;
    HmLoc       loc = {};
    vec<Asgmt>  adj = {};
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
            auto uv1 = cf(c1.half().next().crnr().id);
            auto uv2 = cf(c1.half().prev().crnr().id);
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

    throw std::runtime_error("Not implemented yet");
}

inline void update_to_oppo(const Hmesh& hm, const VecXc& cf, Mbuff& buff) {
    auto dir = buff.dr;
    auto uv0 = buff.uv;

    auto gen_buff_crnr = [&](Crnr c) { return Mbuff {.uv = cf(c.id), .dr = dir, .cid = c.id }; };
    auto gen_buff_half = [&](Half h) {
        double _, r;
        auto uv1 = cf(h.next().crnr().id);
        auto uv2 = cf(h.prev().crnr().id);
        find_extended_intersection(uv0, uv0 + dir, uv1, uv2, _, r);
        return Mbuff{.uv = lerp(uv1, uv2, r), .dr = dir, .hid = h.id, .r = r};
    };

    if (buff.cid != -1) {
        auto c  = hm.crnrs[buff.cid];
        auto h_ = try_get_opposite_half_from_crnr(c, cf, dir, TOLERANCE_HALF);
        auto c_ = try_get_opposite_crnr_from_crnr(c, cf, dir, TOLERANCE_CRNR);
        if (h_.has_value()) { buff = gen_buff_half(h_.value()); return; }
        if (c_.has_value()) { buff = gen_buff_crnr(c_.value()); return; }
    }
    if (buff.hid != -1) {
        auto h  = hm.halfs[buff.hid];
        auto h_ = try_get_opposite_half_from_half(h, cf, uv0, dir, TOLERANCE_HALF);
        auto c_ = try_get_opposite_crnr_from_half(h, cf, uv0, dir, TOLERANCE_CRNR);
        if (h_.has_value()) { buff = gen_buff_half(h_.value()); return; }
        if (c_.has_value()) { buff = gen_buff_crnr(c_.value()); return; }
    }
    throw std::runtime_error("Not implemented yet");
}

inline complex get_face_uv(const Mnode& mn, int fid, const Hmesh& hm, const VecXc& cf) {
    return std::visit(overloaded {
        [&](const HmLocOnP& f) -> complex { return f.uv; },
        [&](const HmLocOnV& v) -> complex { return cf[try_get_crnr(hm, v.id, fid).value().id]; },
        [&](const HmLocOnE& e) -> complex {
            auto h  = try_get_half(hm, e.id, fid).value();
            auto p0 = cf[h.next().crnr().id];
            auto p1 = cf[h.prev().crnr().id];
            return lerp(p0, p1, h.isCanonical() ? e.r : 1 - e.r);
        },
        [&](const auto& _) -> complex { throw std::runtime_error("invalid arguments"); },
    }, mn.loc);
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
