//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_MOTORCYCLE_H
#define METRIKO_MOTORCYCLE_H
#include "../common/utilities.h"
#include "../common/predicates.h"
#include "../hmesh/hmesh.h"
#include "./common.h"

namespace metriko::mc {

struct MotorcycleGraph;

constexpr double TOLERANCE_HALF = 1e-6; //
constexpr double TOLERANCE_CRNR = 2e-6; //
constexpr double TOLERANCE_EDGE = 0.;    // tolerance on two curvs crash close to an edge

struct Mport {
    complex uv;
    complex dr;
    int crnr_id = -1;
    int this_id = -1;
    int next_id = -1;
    int prev_id = -1;
};

struct OnVert { int vid;             bool operator==(const OnVert&) const = default; };
struct OnEdge { int eid; double r;   bool operator==(const OnEdge&) const = default; };
struct OnFace { int fid; complex uv; bool operator==(const OnFace&) const = default; };
using MnodeLoc = std::variant<std::monostate, OnVert, OnEdge, OnFace>;
enum class JunctionType { None, F, T };
enum class JunctionSide {
    None,
    B, // bottom
    L, // left
    R, // right
    T, // top
};

struct Asgmt {
    int curv_id;
    int sgmt_id;
};

struct Mnode {
    JunctionType jt = JunctionType::None;
    MnodeLoc    loc = {};
    vec<Asgmt>  adj = {};
};

struct Msgmt {
    JunctionSide fr_js = JunctionSide::None; // if fr side is crashing, assign side
    JunctionSide to_js = JunctionSide::None; // if to side is crashing, assign side
    int fr_nid  = -1;
    int to_nid  = -1;
    int curv_id = -1;
    int face_id = -1;
    int this_id = -1;
    int prev_id = -1;
    int next_id = -1;
    bool operator==(const Msgmt& s) const { return fr_nid == s.fr_nid && to_nid == s.to_nid; }
};

struct Mbuff {
    complex uv = {0, 0};
    complex dr = {0, 0};
    int  cid = -1; // crnr id
    int  hid = -1; // half id
    double r = -1; // ratio
    bool bgn = false;
    bool end = false;
};

struct Mcurv {
    MotorcycleGraph *mg = nullptr;
    int id = -1;
    vec<Msgmt> sgmts = {};
    Mbuff buff = {};

    bool operator==(const Mcurv &c) const { return id == c.id; }

    void post_process() {
        for (int i = 0; i < sgmts.size(); ++i) {
            sgmts[i].curv_id = id;
            sgmts[i].this_id  = i;
        }
        for (int i = 0; i < sgmts.size() - 1; ++i) {
            sgmts[i].next_id = i + 1;
            sgmts[i + 1].prev_id = i;
        }
    }

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

complex get_face_uv(const Mnode& mn, int fid, const Hmesh& hm, const VecXc& cf) {

    auto on_crnr_uv = [&](Face f, Vert v) {
        for (auto h: f.adjHalfs())
            if (h.crnr().vert() == v) return cf[h.crnr().id];
        throw std::runtime_error("invalid arguments");
    };

    auto on_edge_uv = [&](Face f, Edge e, double r) {
        for (auto h: f.adjHalfs())
            if (h.edge() == e) {
                auto uv0 = cf[h.next().crnr().id];
                auto uv1 = cf[h.prev().crnr().id];
                if ( h.isCanonical()) return lerp(uv0, uv1, r);
                if (!h.isCanonical()) return lerp(uv1, uv0, r);
            }
        throw std::runtime_error("invalid arguments");
    };

    return std::visit(overloaded {
        [&](const OnVert& v) -> complex { return on_crnr_uv(hm.faces[fid], hm.verts[v.vid]); },
        [&](const OnEdge& e) -> complex { return on_edge_uv(hm.faces[fid], hm.edges[e.eid], e.r); },
        [&](const OnFace& f) -> complex { return f.uv; },
        [](std::monostate) -> complex { throw std::runtime_error("Invalid Mnode"); }
    }, mn.loc);
}


struct  MotorcycleGraph {
    const Hmesh &hm;
    const VecXc &cf;
    vec<Mport> mports;
    vec<Mnode> mnodes;
    vec<Mcurv> mcurvs;

    MotorcycleGraph(
        const Hmesh &hm,
        const VecXc &cf,
        const VecXi &matching,
        const VecXi &singular
    ) : hm(hm), cf(cf) {
        gen_ports(singular);

        for (auto v: hm.verts | vw::filter([&](auto& v) { return singular[v.id]; }))
            mnodes.push_back({.jt = JunctionType::F, .loc = OnVert{v.id}});

        //1: Add the first segment for each curve
        mcurvs.reserve(mports.size());
        for (const auto& p : mports) {
            mcurvs.push_back({.mg = this, .id = p.this_id, .buff = {.uv = p.uv, .dr = p.dr, .cid = p.crnr_id, .bgn = true}});
            mcurvs.back().add_segment(hm, cf);
        }

        //2: Add further segments until every curve crash to another curve
        while (rg::any_of(mcurvs, [](auto &c) { return !c.buff.end; })) {
            for (auto &mc: mcurvs) {
                if (mc.buff.end) continue;
                update_to_twin(hm, cf, matching, mc.buff);
                mc.add_segment(hm, cf);
            }
        }

        collect_node_adjacency();
        sort_node_adjacency();

        for (auto &c: mcurvs) c.post_process();
    }

    void gen_ports(const VecXi& singular);
    void collect_node_adjacency();
    void sort_node_adjacency();
};
}

namespace metriko::mc {
inline void MotorcycleGraph::gen_ports(const VecXi &singular) {
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

inline void MotorcycleGraph::collect_node_adjacency() {
    for (auto& mn : mnodes) { mn.adj.clear(); }
    for (auto& mc : mcurvs) {
        for (int i = 0; i < mc.sgmts.size(); ++i) {
            auto& sg = mc.sgmts[i];
            Asgmt as = {mc.id, i};
            if (sg.fr_nid != -1) mnodes[sg.fr_nid].adj.push_back(as);
            if (sg.to_nid != -1) mnodes[sg.to_nid].adj.push_back(as);
        }
    }
}

inline void MotorcycleGraph::sort_node_adjacency() {
    for (int nid = 0; nid < mnodes.size(); ++nid) {
        auto& mn = mnodes[nid];
        if (mn.adj.size() < 3) continue;

        // 【新提案】T-Junction(T)の場合、JunctionSide を使ったトポロジーソートを試みる
        /*
        if (mn.jt == mc::JunctionType::T) {
            bool can_use_topology = true;

            // 接続しているすべてのセグメントが有効な JunctionSide を持っているかチェック
            for (const auto& as : mn.adj) {
                const auto& s = mcurvs[as.curv_id].sgmts[as.sgmt_id];
                mc::JunctionSide js = (s.fr_nid == nid) ? s.fr_js : s.to_js;

                if (js == mc::JunctionSide::None) {
                    can_use_topology = false; // 1つでも None があればトポロジーソートを断念
                    break;
                }
            }

            // 完全なトポロジー情報が揃っている場合のみ、JunctionSide で安全かつ高速にソート
            if (can_use_topology) {
                auto get_side_rank = [&](const Asgmt& as) {
                    const auto& s = mcurvs[as.curv_id].sgmts[as.sgmt_id];
                    mc::JunctionSide js = (s.fr_nid == nid) ? s.fr_js : s.to_js;
                    switch (js) {
                    case mc::JunctionSide::B: return 0;
                    case mc::JunctionSide::R: return 1;
                    case mc::JunctionSide::T: return 2;
                    case mc::JunctionSide::L: return 3;
                    default: return 4;
                    }
                };

                rg::sort(mn.adj, [&](auto& a, auto& b) {
                    return get_side_rank(a) < get_side_rank(b);
                });

                continue; // トポロジーソート完了、次のノードへ
            }

            // can_use_topology == false の場合は、下の「幾何ベースのソート」へフォールバックする
        }
        */

        auto get_fid = [&](const Asgmt& as) { return mcurvs[as.curv_id].sgmts[as.sgmt_id].face_id; };
        auto get_dir = [&](const Asgmt& as) {
            auto& s  = mcurvs[as.curv_id].sgmts[as.sgmt_id];
            auto fr  = s.fr_nid == nid;
            auto uvA = get_face_uv(mnodes[s.fr_nid], s.face_id, hm, cf);
            auto uvB = get_face_uv(mnodes[s.to_nid], s.face_id, hm, cf);
            return fr ? uvB - uvA : uvA - uvB;
        };

        std::visit(overloaded {
            [&](const OnVert& v) {
                std::unordered_map<int, int> ccw_rank;
                int rank = 0;
                for (Half h : hm.verts[v.vid].adjHalfs()) ccw_rank[h.face().id] = rank++;

                rg::sort(mn.adj, [&](auto& a, auto& b) {
                    int rA = ccw_rank.at(get_fid(a));
                    int rB = ccw_rank.at(get_fid(b));
                    if (rA != rB) return rA < rB;
                    return cross(get_dir(a), get_dir(b)) > 0;
                });
            },
            [&](const OnEdge& e) {
                Half h = hm.edges[e.eid].half();
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
            [&](const OnFace& f) {
                rg::sort(mn.adj, [&](auto& a, auto& b) {
                    return std::arg(get_dir(a)) < std::arg(get_dir(b));
                });
            },
            [](std::monostate) {}
        }, mn.loc);
    }
}

inline int Mcurv::resolve_bgn_node(const Hmesh& hm, const bool bgn, const int cid) const {
    if (!bgn) return sgmts.back().to_nid;
    auto l = MnodeLoc{OnVert{hm.crnrs[cid].vert().id}};
    auto i = rg::find(mg->mnodes, l, &Mnode::loc);
    assert(i != mg->mnodes.end() && "start node must exist!");
    return std::distance(mg->mnodes.begin(), i);
}

inline void Mcurv::add_segment(const Hmesh &hm, const VecXc& cf) {
    auto uv0 = buff.uv;
    auto cid = buff.cid;
    auto hid = buff.hid;
    auto bgn = buff.bgn;
    auto fid = cid != -1 ? hm.crnrs[cid].face().id : hm.halfs[hid].face().id;
    update_to_oppo(hm, cf, buff);
    auto uv3 = buff.uv;

    vec<std::tuple<double, double, Msgmt>> candidates;

    auto sgs = vw::all(mg->mcurvs) |
               vw::filter([&](const auto &e) { return e.id != id; }) |
               vw::transform([](const auto &e) -> const auto& { return e.sgmts; }) |
               vw::join |
               vw::filter([&](const auto &s) { return s.face_id == fid; });

    for (auto &s: sgs) {
        auto t   = TOLERANCE_EDGE;
        auto ab  = 0.;
        auto cd  = 0.;
        auto uvA = get_face_uv(mg->mnodes[s.fr_nid], fid, hm, cf);
        auto uvB = get_face_uv(mg->mnodes[s.to_nid], fid, hm, cf);
        auto f1  = find_extended_intersection(uv0, uv3, uvA, uvB, ab, cd);
        auto f2  = ab >= -t && cd >= -t && ab <= 1 + t && cd <= 1 + t;
        if (f1 && f2) candidates.emplace_back(ab, cd, s);
    }

    Msgmt sg{
        .fr_nid  = resolve_bgn_node(hm, bgn, cid),
        .curv_id = id,
        .face_id = fid
    };

    // determine to_nid
    if (candidates.empty()) {

        // if hit to crnr
        if (buff.cid != -1) {
            auto c1 = hm.crnrs[buff.cid];
            auto v1 = c1.vert();
            auto ml = MnodeLoc{OnVert{v1.id}};

            // if hit to other Mnode, return
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

            // otherwise
            auto n = Mnode{.loc = OnVert{v1.id}};
            mg->mnodes.push_back(n);
            sg.to_nid = mg->mnodes.size() - 1;
            sgmts.push_back(sg);
        }
        // otherwise it must hit half
        else {
            auto n = Mnode {};
            if (buff.hid != -1) {
                auto h = hm.halfs[buff.hid];
                n.loc = OnEdge{h.edge().id, h.isCanonical() ? buff.r : 1. - buff.r};
            }
            mg->mnodes.push_back(n);
            sg.to_nid = mg->mnodes.size() - 1;
            sgmts.push_back(sg);
        }
    }
    // intersection happens
    else {
        auto [ab, cd, s0] = rg::min(candidates, [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });

        if (cd >= TOLERANCE_EDGE && cd <= 1 - TOLERANCE_EDGE) {
            // inside face
            if (ab >= TOLERANCE_EDGE && ab <= 1 - TOLERANCE_EDGE) {
                auto mn = Mnode { .jt = JunctionType::T, .loc = OnFace{fid, lerp(uv0, uv3, ab)} };
                mg->mnodes.push_back(mn);
                auto nid = mg->mnodes.size() - 1;

                sg.to_nid = nid;
                sg.to_js  = JunctionSide::B;
                sgmts.push_back(sg);

                auto& cv = mg->mcurvs[s0.curv_id];
                auto  it = rg::find(cv.sgmts, s0);

                Msgmt s1   = *it;
                s1.fr_nid  = nid;
                s1.fr_js   = JunctionSide::R; // todo: need to check direction
                it->to_nid = nid;
                it->to_js  = JunctionSide::L; // same

                cv.sgmts.insert(it + 1, s1);
            }
            else throw std::runtime_error("Not implemented yet");
        }
        else if (cd < TOLERANCE_EDGE) {
            if (ab < TOLERANCE_EDGE) {
                auto uv4 = lerp(uv0, uv3, ab);
                auto uv5 = lerp(
                    get_face_uv(mg->mnodes[s0.fr_nid], fid, hm, cf),
                    get_face_uv(mg->mnodes[s0.to_nid], fid, hm, cf),
                    cd);

                for (Half h : hm.faces[fid].adjHalfs()) {
                    Crnr c = h.crnr();
                    auto uv6 = cf(c.id);
                    std::cout << "c.id: " << c.id << std::endl;
                    std::cout << "cid: " << cid << std::endl;
                    std::cout << "hid: " << hid << std::endl;
                    std::cout << "diff uv4 - uv6: " << abs(uv4 - uv6) << std::endl;
                    std::cout << "diff uv5 - uv6: " << abs(uv5 - uv6) << std::endl;
                }

                throw std::runtime_error("Not implemented yet");
            }
            if (ab > TOLERANCE_EDGE) throw std::runtime_error("Not implemented yet");
            auto& cv = mg->mcurvs[s0.curv_id];
            auto  it = rg::find(cv.sgmts, s0);
            sg.to_nid = it->fr_nid;
            sgmts.push_back(sg);
        }
        else if (cd > 1 - TOLERANCE_EDGE) {
            if (ab < TOLERANCE_EDGE) throw std::runtime_error("Not implemented yet");
            if (ab > TOLERANCE_EDGE) throw std::runtime_error("Not implemented yet");
            auto& cv = mg->mcurvs[s0.curv_id];
            auto  it = rg::find(cv.sgmts, s0);
            sg.to_nid = it->to_nid;
            sgmts.push_back(sg);
        }
        else {
            throw std::runtime_error("Not implemented yet");
        }

        buff.end = true;
    }
}
}

#endif
