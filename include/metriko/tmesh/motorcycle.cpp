//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include "motorcycle.h"

using namespace metriko;

void Mgrph::gen_ports(const VecXi &singular) {
    for (Vert v: hm.verts) {
        if (singular[v.id] == 0) continue;
        vec<Mport> buff0 {}; // the outer scope buffer to assign next/prev
        vec<Mport> buff1 {}; // the inner scope buffer

        for (Half h: v.adjHalfs()) {
            buff1.clear();
            auto a = cf(h.next().crnr().id);
            auto b = cf(h.prev().crnr().id);
            auto c = cf(h.crnr().id);
            auto o = orientation(a, b, c);
            auto ab = b - a;
            auto ac = c - a;
            if (o < 0) throw std::invalid_argument("Orientation should be ccw order");

            int r;
            for (r = 0; r < 4; r++) {
                auto d = get_quater_rot(r);
                if (!is_points_into(a, b, c, a + d) && r > 0) break;
            }

            for (int i = 0; i < 4; i++) {
                auto d = get_quater_rot(r - i);
                if (is_points_into(a, b, c, a + d) ||
                    std::arg(d) == std::arg(ab) ||
                    std::arg(d) == std::arg(ac)
                ) buff1.push_back({.uv = a, .dr = d, .crnr_id = h.next().crnr().id});
            }

            rg::sort(buff1, [&](auto &p0, auto &p1) { return dot(p0.dr, ab) > dot(p1.dr, ab); });
            buff0.insert(buff0.end(), buff1.begin(), buff1.end());
        }

        for (int i = 0; i < buff0.size(); i++) { buff0[i].this_id = mports.size() + i; }
        for (int i = 0; i < buff0.size(); i++) {
            int s = buff0.size();
            buff0[i].prev_id = buff0[(i - 1 + s) % s].this_id;
            buff0[i].next_id = buff0[(i + 1 + s) % s].this_id;
        }
        mports.insert(mports.end(), buff0.begin(), buff0.end());
    }
}

void Mgrph::collect_node_adjacency() {
    for (auto& mn : mnodes) { mn.adj.clear(); }
    for (auto& mc : mcurvs) {
        for (int i = 0; i < mc.sgmts.size(); ++i) {
            auto& sg = mc.sgmts[i];
            Row2i ad = Row2i{mc.id, i};
            if (sg.fr_nid != -1) mnodes[sg.fr_nid].adj.push_back(ad);
            if (sg.to_nid != -1) mnodes[sg.to_nid].adj.push_back(ad);
        }
    }
}

void Mgrph::sort_node_adjacency() {
    for (int nid = 0; nid < mnodes.size(); ++nid) {
        auto& mn = mnodes[nid];
        if (mn.adj.size() < 3) continue;

        auto get_fid = [&](const Row2i& ad) { return mcurvs[ad.x()].sgmts[ad.y()].face_id; };
        auto get_dir = [&](const Row2i& ad) {
            auto& s  = mcurvs[ad.x()].sgmts[ad.y()];
            auto fr  = s.fr_nid == nid;
            auto uvA = get_face_uv(mnodes[s.fr_nid].loc, s.face_id, hm, cf);
            auto uvB = get_face_uv(mnodes[s.to_nid].loc, s.face_id, hm, cf);
            METRIKO_CHECK(abs(uvA - uvB) > 1e-8, "zero-length segment at node {}", nid);
            return fr ? uvB - uvA : uvA - uvB;
        };

        std::visit(overloaded {
            [&](const HmLocOnV& v) {
                umap<int, int> ccw_rank;
                int rank = 0;
                for (Half h : hm.verts[v.id].adjHalfs()) ccw_rank[h.face().id] = rank++;
                rg::sort(mn.adj, [&](auto& a, auto& b) {
                    int rA = ccw_rank.at(get_fid(a));
                    int rB = ccw_rank.at(get_fid(b));
                    if (rA != rB) return rA < rB;
                    return cross(get_dir(a), get_dir(b)) > 0;
                });
            },
            [&](const HmLocOnE& e) {
                Half h = hm.edges[e.id].half();
                int f0 = h.face().id;
                int f1 = h.twin().face().id;

                auto edge_rank = [&](int fid) {
                    if (fid == f0) return 0;
                    if (fid == f1) return 1;
                    return 2;
                };

                rg::sort(mn.adj, [&](auto& a, auto& b) {
                    int rA = edge_rank(get_fid(a));
                    int rB = edge_rank(get_fid(b));
                    if (rA != rB) return rA < rB;
                    return cross(get_dir(a), get_dir(b)) > 0;
                });
            },
            [&](const HmLocOnP& _) {
                rg::sort(mn.adj, [&](auto& a, auto& b) { return std::arg(get_dir(a)) < std::arg(get_dir(b)); });
            },
            [](const auto& _) { METRIKO_FAIL("unexpected location type"); }
        }, mn.loc);
    }
}

int Mcurv::resolve_bgn_node(const Hmesh& hm, const bool bgn, const HmLoc& fr) const {
    if (!bgn) return sgmts.back().to_nid;
    auto l = to_chart_free(hm, fr);
    auto i = rg::find(mg->mnodes, l, &Mnode::loc);
    METRIKO_CHECK(i != mg->mnodes.end(), "start node of the curve must exist");
    return std::distance(mg->mnodes.begin(), i);
}

void Mcurv::add_segment(const Hmesh &hm, const VecXc& cf) {
    auto fr  = buff.loc;
    auto bgn = buff.bgn;
    auto uv0 = get_chart_uv(hm, cf, fr);
    auto fid = get_chart_face(hm, fr).id;
    buff = {.loc = find_ray_intersection(hm, buff.loc, cf, buff.dir), .dir = buff.dir};
    auto uv3 = get_chart_uv(hm, cf, buff.loc);

    vec<std::tuple<double, double, Msgmt>> candidates;

    auto sgs = vw::all(mg->mcurvs) |
               vw::filter([&](const Mcurv &c) { return c.id != id; }) |
               vw::filter([&](const Mcurv &c) { if (bgn) return hm.crnrs[mg->mports[c.id].crnr_id].vert() != hm.crnrs[mg->mports[id].crnr_id].vert(); return true; }) |
               vw::transform([](const Mcurv &c) -> const auto& { return c.sgmts; }) |
               vw::join |
               vw::filter([&](const auto &s) { return s.face_id == fid; });

    for (auto &s: sgs) {
        auto ab  = 0.;
        auto cd  = 0.;
        auto uvA = get_face_uv(mg->mnodes[s.fr_nid].loc, fid, hm, cf);
        auto uvB = get_face_uv(mg->mnodes[s.to_nid].loc, fid, hm, cf);
        auto f1  = find_extended_intersection(uv0, uv3, uvA, uvB, ab, cd);
        auto f2  = ab >= 0. && cd >= 0. && ab <= 1. && cd <= 1.;
        if (f1 && f2) candidates.emplace_back(ab, cd, s);
    }

    Msgmt sg{
        .fr_nid  = resolve_bgn_node(hm, bgn, fr),
        .curv_id = id,
        .face_id = fid
    };

    // determine to_nid
    if (candidates.empty()) {
        auto ml = to_chart_free(hm, buff.loc); // OnC -> OnV, OnH -> OnE

        // if hit to other Mnode on a vert, return
        if (std::holds_alternative<HmLocOnV>(ml)) {
            for (int i = 0; i < mg->mnodes.size(); ++i) {
                auto& mn = mg->mnodes[i];
                if (mn.loc == ml) {
                    mn.jt = JunctionType::T;
                    sg.to_nid = i;
                    sgmts.push_back(sg);
                    buff.end = true;
                    return;
                }
            }
        }

        // otherwise
        mg->mnodes.push_back(Mnode{.loc = ml});
        sg.to_nid = mg->mnodes.size() - 1;
        sgmts.push_back(sg);
    }
    // intersection happens
    else {
        auto [ab, cd, sg_cd] = rg::min(candidates, [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });
        METRIKO_CHECK(ab >= TOLERANCE_EDGE_AB, "a curve passes too close to a vertex (not handled yet)");
        bool cd_bgn_snappable = cd < TOLERANCE_EDGE_CD;
        bool cd_end_snappable = cd > 1 - TOLERANCE_EDGE_CD;

        auto& cv = mg->mcurvs[sg_cd.curv_id];
        auto it = rg::find_if(cv.sgmts, [&](const Msgmt& s) { return s.fr_nid == sg_cd.fr_nid && s.to_nid == sg_cd.to_nid; });
        buff.end = true;

        auto snap = [&](int nid) {
            auto jt = mg->mnodes[nid].jt;
            if (jt == JunctionType::F || jt == JunctionType::C) return false;
            // if (jt == JunctionType::None) return false;

            for (const auto& c: mg->mcurvs) {
            for (const auto& s: c.sgmts) {
                if (s.to_nid == nid && s.face_id == fid) {
                    auto o  = get_face_uv(mg->mnodes[nid].loc, fid, hm, cf);
                    auto d0 = get_face_uv(mg->mnodes[s.fr_nid].loc , fid, hm, cf) - o;
                    auto d1 = get_face_uv(mg->mnodes[sg.fr_nid].loc, fid, hm, cf) - o;
                    auto l0 = abs(d0);
                    auto l1 = abs(d1);
                    if (l0 < EPS || l1 < EPS || dot(d0 / l0 , d1 / l1) > EPS) return false;
                };
            }}

            mg->mnodes[nid].jt = JunctionType::T;
            sg.to_nid = nid;
            sgmts.push_back(sg);
            return true;
        };

        if (cd_bgn_snappable && snap(it->fr_nid)) { std::cout << "cd bgn snappable, fid: " << fid << std::endl; return; }
        if (cd_end_snappable && snap(it->to_nid)) { std::cout << "cd end snappable, fid: " << fid << std::endl; return; }

        { // hit in middle case, or close to singular node
            auto mn = Mnode{.jt = JunctionType::T, .loc = HmLocOnP{fid, lerp(uv0, uv3, ab)}};
            mg->mnodes.push_back(mn);
            auto nid = mg->mnodes.size() - 1;

            sg.to_nid = nid;
            sgmts.push_back(sg);

            Msgmt s1   = *it;
            s1.fr_nid  = nid;
            it->to_nid = nid;

            cv.sgmts.insert(it + 1, s1);
        }
    }
}
