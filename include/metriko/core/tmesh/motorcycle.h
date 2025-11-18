//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_MOTORCYCLE_H
#define METRIKO_MOTORCYCLE_H
#include "../common/utilities.h"
#include "../common/predicates.h"
#include "../hmesh/hmesh.h"
#include "../hmesh/utilities.h"

namespace metriko {
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

enum MvertType { None, HitL, HitR, HitB };

class Mvert : Melem {
public:
    complex uv;
    Mcurv *crash = nullptr;
    MvertType type = None;
    std::optional<Half> half; // part of the halfedge flag

    Mvert(
        const MotorcycleGraph *g,
        const complex uv,
        const std::optional<Half>& half = std::nullopt,
        Mcurv *crash = nullptr,
        const MvertType side = None
    ) : Melem(g), uv(uv), crash(crash), type(side), half(half)
    { }
};

class Msgmt : Melem {
public:
    Mcurv *curv;
    Face face;
    Mvert fr;
    Mvert to;
    int id = -1;
    int prev_id = -1;
    int next_id = -1;

    Msgmt(
        const MotorcycleGraph *g,
        Mcurv *curv,
        const Face &face,
        const Mvert &fr,
        const Mvert &to
    ) : Melem(g), curv(curv), face(face), fr(fr), to(to) {
    }

    bool operator==(const Msgmt &rhs) const {
        return face  == rhs.face &&
               fr.uv == rhs.fr.uv &&
               to.uv == rhs.to.uv;
    }

    complex diff() const { return to.uv - fr.uv; }
    const Msgmt &next() const;
    const Msgmt &prev() const;
};

class Cache {
public:
    Half half;
    double ratio;
    complex dir;
    bool intersected = false;
};

class Mcurv : Melem {
public:
    const Mport &port;
    std::vector<Msgmt> sgmts;
    Cache cache;

    explicit Mcurv(const MotorcycleGraph *g, const Mport &p) : Melem(g), port(p), cache() { }
    int id() const { return port.id; }
    bool operator==(const Mcurv &rhs) const { return port.id == rhs.port.id; }

    void add_segment(
        complex dir,
        complex uv0,
        Half h0,
        Half h1
    );

    void split_segment(
        const Msgmt &target,
        const Msgmt &crash,
        const double ratio
    ) {
        for (auto it = sgmts.begin(); it != sgmts.end(); ++it) {
            if (*it == target) {
                const bool ccw = cross(crash.diff(), it->diff()) > 0;
                const auto uv  = lerp(it->fr.uv, it->to.uv, ratio);
                const auto fr  = Mvert(mg, uv, std::nullopt, crash.curv, ccw ? HitL : HitR);
                const auto to  = Mvert(mg, it->to.uv, it->to.half, it->to.crash, it->to.type);
                it->to = Mvert(mg, uv, std::nullopt, crash.curv, ccw ? HitR : HitL);
                sgmts.insert(it + 1, Msgmt(mg, this, it->face, fr, to));
                return;
            }
        }
    }

    void post_process() {
        for (int i = 0; i < sgmts.size(); ++i) sgmts[i].id = i;
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
    std::vector<Mport> mports;
    std::vector<Mcurv> mcurvs;
    std::vector<Cache> caches;

    MotorcycleGraph(
        const Hmesh &hm,
        const VecXc &cf,
        const VecXi &matching,
        const VecXi &singular
    ) : hm(hm), cf(cf) {
        gen_ports(singular);
        for (auto &p: mports) mcurvs.emplace_back(this, p);

        // Add the first segment for each curve
        for (auto &c: mcurvs) {
            std::optional<Half> h0;
            std::optional<Half> h1;
            const auto &p = c.port;
            const auto &v = p.vert;
            for (auto h: p.face.adjHalfs()) {
                if      (h.tail() == v) h0 = h;
                else if (h.head() != v) h1 = h;
            }
            c.add_segment(p.dir, p.uv, h0.value(), h1.value());
        }

        // Add further segments until every curve crash to another curve
        // Need to consider parallel intersection (e.g., bumpy-cube case)
        while (rg::any_of(mcurvs, [](auto &e) { return !e.cache.intersected; })) {
            for (auto &c: mcurvs) {
                auto [ch, cr, cd, ci] = c.cache;
                if (ci) continue;
                auto h0 = ch.twin();
                auto uv = lerp(cf(h0.next().crnr().id), cf(h0.prev().crnr().id), cr);
                auto m  = (h0.isCanonical() ? -1 : 1) * matching[h0.edge().id];
                auto d  = std::polar(1., PI / 2 * m) * cd;
                auto h1 = get_opposite_half(h0, cf, uv, d);
                c.add_segment(d, uv, h0, h1);
            }
        }

        // when done, assign index and next/prev info to each segment
        for (auto &c: mcurvs) c.post_process();
    }

    void gen_ports(const VecXi &singular);
};
}

namespace metriko {
inline const Msgmt &Msgmt::next() const { return curv->sgmts[next_id]; }
inline const Msgmt &Msgmt::prev() const { return curv->sgmts[prev_id]; }

inline void MotorcycleGraph::gen_ports(const VecXi &singular) {
    for (Vert v: hm.verts) {
        if (singular[v.id] == 0) continue;
        std::vector<Mport> buff0 {}; // the outer scope buffer to assign next/prev
        std::vector<Mport> buff1 {}; // the inner scope buffer

        for (Half h: v.adjHalfs()) {
            buff1.clear();
            auto a  = cf(h.next().crnr().id);
            auto b  = cf(h.prev().crnr().id);
            auto c  = cf(h.crnr().id);
            auto o  = orientation(a, b, c);
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
                ) buff1.emplace_back(-1, v, h.face(), a, d);
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

inline void Mcurv::add_segment(
    const complex dir, // direction
    const complex uv0, // origin
    const Half h0,     // halfedge fr in the face
    const Half h1      // halfedge to in the face
) {
    complex uv1 = mg->cf(h1.prev().crnr().id);
    complex uv2 = mg->cf(h1.next().crnr().id);
    double r_ab, r_cd;
    bool r = find_extended_intersection(uv0, uv0 + dir, uv1, uv2, r_ab, r_cd);
    assert(r);
    r_cd = std::clamp(r_cd, 0., 1.);
    auto uv3 = lerp(uv1, uv2, r_cd);
    auto f = h1.face();

    std::vector<std::tuple<double, double, Msgmt>> candidates;

    auto sgs = vw::all(mg->mcurvs) |
               vw::filter([&](auto &e) { return e != *this; }) |
               vw::transform([](auto &e) { return e.sgmts; }) |
               vw::join |
               vw::filter([&f](auto &sg) { return sg.face.id == f.id; });

    for (auto &sg: sgs) {
        double ab = 0;
        double cd = 0;
        if (find_strict_intersection(uv0, uv3, sg.fr.uv, sg.to.uv, ab, cd, 0.))
            candidates.emplace_back(ab, cd, sg);
    }

    if (candidates.empty()) {
        auto v1 = Mvert(mg, uv0, h0);
        auto v2 = Mvert(mg, uv3, h1);
        sgmts.emplace_back(mg, this, f, v1, v2);
        cache = Cache(h1, r_cd, dir);
    } else {
        auto [ab, cd, sg] = rg::min(candidates, [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });
        auto v1 = Mvert(mg, uv0, h0);
        auto v2 = Mvert(mg, lerp(uv0, uv3, ab), std::nullopt, sg.curv, HitB);
        sgmts.emplace_back(mg, this, f, v1, v2);
        cache = Cache(h1, r_cd, dir, true);
        sg.curv->split_segment(sg, sgmts.back(), cd);
    }
}
}

#endif
