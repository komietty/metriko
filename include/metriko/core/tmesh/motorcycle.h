//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_MOTORCYCLE_H
#define METRIKO_MOTORCYCLE_H
#include <utility>

#include "../common/utilities.h"
#include "../common/predicates.h"
#include "../hmesh/hmesh.h"
#include "../hmesh/utilities.h"

namespace metriko::mc {
template <class T> using vec = std::vector<T>;
template <class T> using set = std::set<T>;
template <class T> using opt = std::optional<T>;

class MotorcycleGraph;
class Mcurv;

class Mport {
public:
    int id;
    Vert vert;
    Face face;
    complex uv;
    complex dir;
    int next = -1;
    int prev = -1;

    Mport(
        const int id,
        const Vert vert,
        const Face face,
        const complex uv,
        const complex dir
    ): id(id), vert(vert), face(face), uv(uv), dir(dir) {
    }
};

class Melem {
public:
    const MotorcycleGraph *mg;
    explicit Melem(const MotorcycleGraph *mg) : mg(mg) { }
};

enum MvertType { None, First, HitL, HitR, HitB };

class Mvert : Melem {
public:
    complex uv;
    Mcurv *crash = nullptr; //
    MvertType type = None;
    opt<std::pair<Half, double>> cut; // part of the halfedge flag

    Mvert(
        const MotorcycleGraph *g,
        const complex uv,
        const opt<std::pair<Half, double>>& cut = std::nullopt,
        Mcurv *crash = nullptr,
        const MvertType side = None
    ) : Melem(g), uv(uv), crash(crash), type(side), cut(cut) { }
};

class Msgmt : Melem {
public:
    Mcurv *curv;
    Mvert fr;
    Mvert to;
    Face face;
    complex dir;
    int id = -1;
    int curv_id = -1;
    int prev_id = -1;
    int next_id = -1;

    Msgmt(
        const MotorcycleGraph *g,
        Mcurv *curv,
        Face face,
        Mvert fr,
        Mvert to,
        const complex dir
    ) : Melem(g),
        curv(curv),
        fr(std::move(fr)),
        to(std::move(to)),
        face(face),
        dir(dir) {}

    bool operator==(const Msgmt &rhs) const {
        return face  == rhs.face  &&
               fr.uv == rhs.fr.uv &&
               to.uv == rhs.to.uv;
    }

    complex diff() const { return to.uv - fr.uv; }
    const Msgmt &next() const;
    const Msgmt &prev() const;
};

class Mcurv : Melem {
public:
    const Mport &port;
    vec<Msgmt> sgmts;
    bool intersected = false;

    explicit Mcurv(const MotorcycleGraph *g, const Mport &p) : Melem(g), port(p) { }
    int id() const { return port.id; }
    bool operator==(const Mcurv &rhs) const { return port.id == rhs.port.id; }
    void add_segment(complex dir, complex uv0, double r0, Half h0, Half h1, bool first);

    void post_process() {
        for (int i = 0; i < sgmts.size(); ++i) {
            sgmts[i].curv_id = id();
            sgmts[i].id = i;
        }
        for (int i = 0; i < sgmts.size() - 1; ++i) {
            sgmts[i].next_id = i + 1;
            sgmts[i + 1].prev_id = i;
        }
    }
};

class MotorcycleGraph {
public:
    const Hmesh &hm;
    const VecXc &cf;
    vec<Mport> mports;
    vec<Mcurv> mcurvs;


    MotorcycleGraph(
        const Hmesh &hm,
        const VecXc &cf,
        const VecXi &matching,
        const VecXi &singular
    ) : hm(hm), cf(cf) {
        gen_ports(singular);
        for (auto &p: mports) { mcurvs.emplace_back(this, p); }

        // Add the first segment for each curve
        for (auto &c: mcurvs) {
            opt<Half> h0;
            opt<Half> h1;
            const auto &pr = c.port;
            const auto &vt = pr.vert;
            for (auto h: pr.face.adjHalfs()) {
                if      (h.tail() == vt) h0 = h;
                else if (h.head() != vt) h1 = h;
            }
            assert(h0.value().face() == h1.value().face());
            c.add_segment(pr.dir, pr.uv, 1, h0.value(), h1.value(), true); //todo ratio = 1 not 0 does not make sense. Need check.
        }

        // Add further segments until every curve crash to another curve
        // Need to consider parallel intersection (e.g., bumpy-cube case)
        while (rg::any_of(mcurvs, [](auto &e) { return !e.intersected; })) {
            for (auto &c: mcurvs) {
                if (c.intersected) continue;
                auto h_ = c.sgmts.back().to.cut.value().first;
                auto r_ = c.sgmts.back().to.cut.value().second;
                auto d_ = c.sgmts.back().dir;
                auto h0 = h_.twin();
                auto r0 = 1. - r_;
                auto uv = lerp(cf(h0.prev().crnr().id), cf(h0.next().crnr().id), r0);
                auto m  = (h0.isCanonical() ? -1 : 1) * matching[h0.edge().id];
                auto d  = std::polar(1., PI / 2 * m) * d_;
                auto h1 = get_opposite_half(h0, cf, uv, d);
                c.add_segment(d, uv, r0, h0, h1, false);
            }
        }

        // when done, assign index and next/prev info to each segment
        for (auto &c: mcurvs) c.post_process();
    }

    void gen_ports(const VecXi &singular);
};
}

namespace metriko::mc {
inline const Msgmt &Msgmt::next() const { return curv->sgmts[next_id]; }
inline const Msgmt &Msgmt::prev() const { return curv->sgmts[prev_id]; }

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
                )
                    buff1.emplace_back(-1, v, h.face(), a, d);
            }

            rg::sort(buff1, [&](auto &p0, auto &p1) { return dot(p0.dir, ab) > dot(p1.dir, ab); });
            buff0.insert(buff0.end(), buff1.begin(), buff1.end());
        }

        for (int i = 0; i < buff0.size(); i++) { buff0[i].id = mports.size() + i; }
        for (int i = 0; i < buff0.size(); i++) {
            int s = buff0.size();
            buff0[i].prev = buff0[(i - 1 + s) % s].id;
            buff0[i].next = buff0[(i + 1 + s) % s].id;
        }
        mports.insert(mports.end(), buff0.begin(), buff0.end());
    }
}

inline void split_segment(
    const MotorcycleGraph *mg,
    const Msgmt &s0, // segment to be split (might be reference of copy)
    const Msgmt &s1, // segment crashing
    const double r   // ratio for s0
) {
    Mcurv* c0 = s0.curv;
    Mcurv* c1 = s1.curv;
    auto it = rg::find(c0->sgmts, s0);
    const bool ccw = cross(s1.diff(), s0.diff()) > 0;
    const auto &fr = s0.fr;
    const auto &to = s0.to;
    const auto uv2 = lerp(fr.uv, to.uv, r);
    const auto vfr = Mvert(mg, uv2, std::nullopt, c1, ccw ? HitL : HitR);
    const auto vto = Mvert(mg, to.uv, to.cut, to.crash, to.type);
    it->to = Mvert(mg, uv2, std::nullopt, c1, ccw ? HitR : HitL);
    c0->sgmts.insert(it + 1, Msgmt(mg, c0, s0.face, vfr, vto, s0.dir));
}

inline void Mcurv::add_segment(
    const complex dir, // direction of the new segment in uv space
    const complex uv0, // origin in fr size
    const double r0,   // ratio of fr side
    const Half h0,     // halfedge fr in the face
    const Half h1,     // halfedge to in the face
    const bool first   // the first segment or not
) {
    complex uv1 = mg->cf(h1.prev().crnr().id);
    complex uv2 = mg->cf(h1.next().crnr().id);
    double r_, r1;
    bool r = find_extended_intersection(uv0, uv0 + dir, uv1, uv2, r_, r1);
    assert(r);
    r1 = std::clamp(r1, 0., 1.);
    auto uv3 = lerp(uv1, uv2, r1);
    auto f = h1.face();

    vec<std::tuple<double, double, Msgmt>> candidates;

    auto sgs = vw::all(mg->mcurvs) |
               vw::filter([&](auto &e) { return e != *this; }) |
               vw::transform([](auto &e) { return e.sgmts; }) |
               vw::join |
               vw::filter([&f](auto &sg) { return sg.face.id == f.id; });

    for (auto &sg: sgs) {
        double ab = 0;
        double cd = 0;
        if (find_strict_intersection(uv0, uv3, sg.fr.uv, sg.to.uv, ab, cd))
            candidates.emplace_back(ab, cd, sg);
    }

    if (candidates.empty()) {
        auto v1 = Mvert(mg, uv0, std::pair(h0, r0), this, first ? First : None);
        auto v2 = Mvert(mg, uv3, std::pair(h1, r1));
        sgmts.emplace_back(mg, this, f, v1, v2, dir);
    } else {
        auto [ab, cd, sg_copy] = rg::min(candidates, [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });
        auto v1 = Mvert(mg, uv0, std::pair(h0, r0), this, first ? First : None);
        auto v2 = Mvert(mg, lerp(uv0, uv3, ab), std::nullopt, sg_copy.curv, HitB);
        sgmts.emplace_back(mg, this, f, v1, v2, dir);
        split_segment(mg, sg_copy, sgmts.back(), cd);
        intersected = true;
    }
}
}

#endif
