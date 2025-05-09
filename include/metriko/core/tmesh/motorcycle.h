//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_MOTORCYCLE_H
#define METRIKO_MOTORCYCLE_H
#include "../common/utilities.h"
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
        const MotorcycleGraph *graph;
        explicit Melem(const MotorcycleGraph *g): graph(g) {}
    };

    enum MvertType { None, HitL, HitR, HitB };

    class Mvert : Melem {
    public:
        complex uv;
        Mcurv *crash = nullptr;
        MvertType type = None;

        Mvert(
            const MotorcycleGraph *g,
            const complex uv,
            Mcurv *crash,
            const MvertType side
        ) : Melem(g), uv(uv), crash(crash), type(side) { }
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
        ): Melem(g), curv(curv), face(face), fr(fr), to(to) {
        }

        bool operator==(const Msgmt &rhs) const {
            return face == rhs.face &&
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

        explicit Mcurv(const MotorcycleGraph *g, const Mport &port): Melem(g), port(port), cache() { }

        bool operator==(const Mcurv &rhs) const { return port.id == rhs.port.id; }
        void add_segment_init(const VecXc &cfn, const std::vector<Mcurv> &curvs, const VecXi &matching);
        void add_segment_next(const VecXc &cfn, const std::vector<Mcurv> &curvs, const VecXi &matching);
        void add_segment(const VecXc &cfn, const std::vector<Mcurv> &curvs, const Mcurv &exception, complex dir, complex uv0, Half h);

        void split_segment(
            const Msgmt &target,
            const Msgmt &crash,
            const double ratio
        ) {
            for (auto it = sgmts.begin(); it != sgmts.end(); ++it) {
                if (*it == target) {
                    const bool ccw = cross(crash.diff(), it->diff()) > 0;
                    const auto uv = lerp(it->fr.uv, it->to.uv, ratio);
                    const auto fr = Mvert(graph, uv, crash.curv, ccw ? HitL : HitR);
                    const auto to = Mvert(graph, it->to.uv, it->to.crash, it->to.type);
                    it->to = Mvert(graph, uv, crash.curv, ccw ? HitR : HitL);
                    sgmts.insert(it + 1, Msgmt(graph, this, it->face, fr, to));
                    return;
                }
            }
        }

        void post_process() {
            assert(cache.intersected);
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

        MotorcycleGraph(
            const Hmesh &hmesh,
            const VecXc &cfn,
            const VecXi &matching,
            const VecXi &singular
        ): cfn(cfn) {
            gen_ports(hmesh, cfn, singular);
            for (auto &p: mports) mcurvs.emplace_back(this, p);
            for (auto &e: mcurvs) e.add_segment_init(cfn, mcurvs, matching);
            while (rg::any_of(mcurvs, [](auto &e) { return !e.cache.intersected; })) {
                for (auto &s: mcurvs) s.add_segment_next(cfn, mcurvs, matching);
            }
            for (auto &s: mcurvs) s.post_process();
        }

        void gen_ports(const Hmesh &mesh, const VecXc &cfn, const VecXi &singular);
    };
}

namespace metriko {
    inline const Msgmt &Msgmt::next() const { return curv->sgmts[next_id]; }
    inline const Msgmt &Msgmt::prev() const { return curv->sgmts[prev_id]; }

    inline void MotorcycleGraph::gen_ports(
        const Hmesh &mesh,
        const VecXc &cfn,
        const VecXi &singular
    ) {
        std::vector<Mport> joint;
        std::vector<Mport> ports;
        int counter = 0;
        for (auto v: mesh.verts | vw::filter([&](auto v_) { return singular[v_.id] != 0; })) {
            ports.clear();
            for (Half h: v.adjHalfs()) {
                auto a = cfn(h.next().crnr().id);
                auto b = cfn(h.prev().crnr().id);
                auto c = cfn(h.crnr().id);
                auto o = orientation(a, b, c);
                std::vector<Mport> temps;
                if (o < 0) throw std::invalid_argument("Orientation should be ccw order");
                int r;
                for (r = 0; r < 4; r++) { if (!points_into(get_quater_rot(r), a, b, c) && r > 0) break; }
                for (int i = 0; i < 4; i++) {
                    complex dir = get_quater_rot(r - i);
                    if (points_into(dir, a, b, c) ||
                        std::arg(dir) == std::arg(b - a) ||
                        std::arg(dir) == std::arg(c - a)
                    ) {
                        temps.emplace_back(-1, v, h.face(), a, dir);
                    }
                }

                rg::sort(temps, [&](const auto &p0, const auto &p1) {
                    complex dir = b - a;
                    double dotA = dot(p0.dir, dir);
                    double dotB = dot(p1.dir, dir);
                    return dotA > dotB;
                });

                for (auto &p: temps) {
                    p.id = counter;
                    counter++;
                }

                ports.insert(ports.end(), temps.begin(), temps.end());
            }

            for (int i = 0; i < ports.size(); i++) {
                int s = ports.size();
                ports[i].prev = ports[(i - 1 + s) % s].id;
                ports[i].next = ports[(i + 1 + s) % s].id;
            }

            joint.insert(joint.end(), ports.begin(), ports.end());
        }

        mports = joint;
    }

    inline void Mcurv::add_segment_init(
        const VecXc &cfn,
        const std::vector<Mcurv> &curvs,
        const VecXi &matching
    ) {
        Half h;
        for (Half h_: port.face.adjHalfs())
            if (h_.head().id != port.vert.id && h_.tail().id != port.vert.id) h = h_;
        complex uv0 = port.uv;
        complex dir = port.dir;
        add_segment(cfn, curvs, *this, dir, uv0, h);
    }


    inline void Mcurv::add_segment_next(
        const VecXc &cfn,
        const std::vector<Mcurv> &curvs,
        const VecXi &matching
    ) {
        // need to consider parallel intersection (e.g., bumpy cube case)
        if (cache.intersected) return;
        auto h = cache.half.twin();
        auto m = (h.isCanonical() ? 1 : -1) * matching[h.edge().id];
        auto uv0 = lerp(cfn(h.next().crnr().id), cfn(h.prev().crnr().id), cache.ratio);
        auto dir = std::polar(1., PI / 2 * m) * cache.dir;
        auto oh = get_oppsite_half(cfn, uv0, dir, h);
        add_segment(cfn, curvs, *this, dir, uv0, oh);
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
            auto v1 = Mvert(graph, uv0, nullptr, None);
            auto v2 = Mvert(graph, uv3, nullptr, None);
            sgmts.emplace_back(graph, this, h.face(), v1, v2);
            cache = Cache(h, r_cd, dir);
        } else {
            auto [ab, cd, sg] = rg::min(candidates, [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });
            auto v1 = Mvert(graph, uv0, nullptr, None);
            auto v2 = Mvert(graph, lerp(uv0, uv3, ab), sg.curv, HitB);
            sgmts.emplace_back(graph, this, h.face(), v1, v2);
            cache = Cache(h, r_cd, dir, true);
            sg.curv->split_segment(sg, sgmts.back(), cd);
        }
    }
}

#endif
