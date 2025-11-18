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
    ) : id(id), vert(vert), face(face), uv(uv), dir(dir) {
    }
};

class Melem {
public:
    const MotorcycleGraph *mg;

    explicit Melem(const MotorcycleGraph *mg) : mg(mg) {
    }
};

enum MvertType { None, HitL, HitR, HitB };

class Mvert : Melem {
public:
    complex uv;
    Mcurv *crash = nullptr;
    MvertType type = None;
    // std::optional<Half>  half // is part of the halfedge?
    // std::optional<Mport> port // is coming from the port?

    Mvert(
        const MotorcycleGraph *g,
        const complex uv,
        Mcurv *crash = nullptr,
        const MvertType side = None
    ) : Melem(g), uv(uv), crash(crash), type(side) {
    }
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

    explicit Mcurv(const MotorcycleGraph *g, const Mport &port) : Melem(g), port(port), cache() { }
    int id() const { return port.id; }
    bool operator==(const Mcurv &rhs) const { return port.id == rhs.port.id; }

    void add_segment(
        const VecXc &cfn,
        const std::vector<Mcurv> &curvs,
        const Mcurv &exception,
        complex dir,
        complex uv0,
        Half h
    );

    void split_segment(
        const Msgmt &target,
        const Msgmt &crash,
        const double ratio
    ) {
        for (auto it = sgmts.begin(); it != sgmts.end(); ++it) {
            if (*it == target) {
                const bool ccw = cross(crash.diff(), it->diff()) > 0;
                const auto uv = lerp(it->fr.uv, it->to.uv, ratio);
                const auto fr = Mvert(mg, uv, crash.curv, ccw ? HitL : HitR);
                const auto to = Mvert(mg, it->to.uv, it->to.crash, it->to.type);
                it->to = Mvert(mg, uv, crash.curv, ccw ? HitR : HitL);
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
    const VecXc &cfn;
    std::vector<Mport> mports;
    std::vector<Mcurv> mcurvs;
    std::vector<Cache> caches;

    MotorcycleGraph(
        const Hmesh &hm,
        const VecXc &cf,
        const VecXi &matching,
        const VecXi &singular
    ) : cfn(cf) {
        gen_ports(hm, cfn, singular);
        for (auto &p: mports) mcurvs.emplace_back(this, p);

        // Add the first segment for each curve
        for (auto &c: mcurvs) {
            std::optional<Half> found;
            const auto &p = c.port;
            const auto &v = p.vert;
            for (Half h: p.face.adjHalfs())
                if (h.head() != v && h.tail() != v) found = h;
            assert(found.has_value());
            c.add_segment(cfn, mcurvs, c, p.dir, p.uv, found.value());
        }

        // Add further segments until every curve crash to another curve
        while (rg::any_of(mcurvs, [](auto &e) { return !e.cache.intersected; })) {
            for (auto &c: mcurvs) {
                // need to consider parallel intersection (e.g., bumpy cube case)
                auto [ch, cr, cd, ci] = c.cache;
                if (ci) continue;
                auto h = ch.twin();
                auto m = (h.isCanonical() ? 1 : -1) * matching[h.edge().id];
                auto cf0 = cfn(h.next().crnr().id);
                auto cf1 = cfn(h.prev().crnr().id);
                auto uv  = lerp(cf0, cf1, cr);
                auto dir = std::polar(1., PI / 2 * m) * cd;
                c.add_segment(cfn, mcurvs, c, dir, uv, get_opposite_half(cfn, uv, dir, h));
            }
        }

        // debug
        std::vector<glm::vec3> vis_port;
        for (const auto &p: mports)
            vis_port.emplace_back(p.vert.pos().x(), p.vert.pos().y(), p.vert.pos().z());
        auto vq = polyscope::registerPointCloud("vis_port", vis_port);
        vq->setPointRadius(0.005);
        vq->resetTransform();

        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2> > es;
        size_t counter = 0;
        for (auto c: mcurvs) {
            for (auto &s: c.sgmts) {
                Row3d p1 = conversion_2d_3d(s.face, cf, s.fr.uv);
                Row3d p2 = conversion_2d_3d(s.face, cf, s.to.uv);
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array{counter, counter + 1});
                counter += 2;
            }
        }
        auto c = polyscope::registerCurveNetwork("segments", ns, es);
        c->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.0005);

        // when done, assign index and next/prev info to each segment
        for (auto &s: mcurvs) s.post_process();
    }

    void gen_ports(
        const Hmesh &hm,
        const VecXc &cf,
        const VecXi &singular
    );
};
}

namespace metriko {
inline const Msgmt &Msgmt::next() const { return curv->sgmts[next_id]; }
inline const Msgmt &Msgmt::prev() const { return curv->sgmts[prev_id]; }

inline void MotorcycleGraph::gen_ports(
    const Hmesh &hm,
    const VecXc &cf,
    const VecXi &singular
) {
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
            for (r = 0; r < 4; r++)
                if (!is_points_into(a, b, c, a + get_quater_rot(r)) && r > 0) break;

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

        for (int i = 0; i < buff0.size(); i++) {
            int s = buff0.size();
            buff0[i].prev = buff0[(i - 1 + s) % s].id;
            buff0[i].next = buff0[(i + 1 + s) % s].id;
        }
        mports.insert(mports.end(), buff0.begin(), buff0.end());
    }

    for (int i = 0; i < mports.size(); i++) mports[i].id = i; // assign index
}

inline void Mcurv::add_segment(
    const VecXc &cfn,
    const std::vector<Mcurv> &curvs,
    const Mcurv &exception,
    const complex dir,
    const complex uv0,
    const Half h
) {
    complex uv1 = cfn(h.prev().crnr().id);
    complex uv2 = cfn(h.next().crnr().id);
    double r_ab, r_cd;
    bool f = find_extended_intersection(uv0, uv0 + dir, uv1, uv2, r_ab, r_cd);
    assert(f);

    complex uv3 = lerp(uv1, uv2, r_cd);

    std::vector<std::tuple<double, double, Msgmt>> candidates;

    auto sgs = vw::all(curvs) |
               vw::filter([&](auto &e) { return e != exception; }) |
               vw::transform([](auto &e) { return e.sgmts; }) |
               vw::join |
               vw::filter([&h](auto &sg) { return sg.face.id == h.face().id; });

    for (auto &sg: sgs) {
        double ab = 0;
        double cd = 0;
        if (find_strict_intersection(uv0, uv3, sg.fr.uv, sg.to.uv, ab, cd))
            candidates.emplace_back(ab, cd, sg);
    }

    if (candidates.empty()) {
        auto v1 = Mvert(mg, uv0);
        auto v2 = Mvert(mg, uv3);
        sgmts.emplace_back(mg, this, h.face(), v1, v2);
        cache = Cache(h, r_cd, dir);
    } else {
        auto [ab, cd, sg] = rg::min(candidates, [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });
        auto v1 = Mvert(mg, uv0);
        auto v2 = Mvert(mg, lerp(uv0, uv3, ab), sg.curv, HitB);
        sgmts.emplace_back(mg, this, h.face(), v1, v2);
        cache = Cache(h, r_cd, dir, true);
        sg.curv->split_segment(sg, sgmts.back(), cd);
    }
}
}

#endif
