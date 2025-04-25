//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_MOTORCYCLE_H
#define METRIKO_MOTORCYCLE_H
#include "../common/utilities.h"
#include "../hmesh/hmesh.h"

namespace metriko {
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
        ): id(id), vert(vert), face(face), uv(uv), dir(dir) { }
    };


    enum CrashType { None, Left, Right, Crash };
    enum MvertType { None_, Sing, HitL, HitR, HitB };
    class Mcurv;
    class MotorcycleGraph;

    class Mvert {
    public:
        const MotorcycleGraph* graph;
        complex uv;
        Mcurv* crash = nullptr;
        CrashType side = None;
        std::optional<Half> half = std::nullopt;

        Mvert(
            const MotorcycleGraph* graph,
            const complex uv,
            Mcurv* crash,
            const CrashType side,
            const std::optional<Half> &half = std::nullopt
            ) : graph(graph), uv(uv), crash(crash), side(side), half(half) { }
    };

    class Msgmt {
    public:
        const MotorcycleGraph* graph;
        Mcurv* curv;
        Face face;
        Mvert fr;
        Mvert to;
        int id = -1;
        int prev_id = -1;
        int next_id = -1;

        Msgmt(
            const MotorcycleGraph* graph,
            Mcurv *curv,
            const Face &face,
            const Mvert& fr,
            const Mvert& to
        ): graph(graph), curv(curv), face(face), fr(fr), to(to) { }

        bool operator==(const Msgmt& rhs) const {
            return face.id == rhs.face.id &&
                   fr.uv == rhs.fr.uv &&
                   to.uv == rhs.to.uv;
        }


        complex diff() const { return to.uv - fr.uv; }
        const Msgmt& next() const;
        const Msgmt& prev() const;
    };

    class Cache {
    public:
        Half half;
        double ratio;
        complex dir;
        bool intersected = false;
    };


    class Mcurv {
    public:
        const MotorcycleGraph* graph;
        const Mport& port;
        std::vector<Msgmt> sgmts;
        Cache cache;

        explicit Mcurv(const MotorcycleGraph* graph, const Mport &port): graph(graph), port(port), cache() { }

        bool operator==(const Mcurv &rhs) const { return port.id == rhs.port.id; }

        void add_segment_init(const VecXc &cfn, const std::vector<Mcurv> &edges, const VecXi &matching);
        void add_segment_next(const VecXc &cfn, const std::vector<Mcurv> &edges, const VecXi &matching);
        void add_segment     (const VecXc &cfn, const std::vector<Mcurv> &edges, const Mcurv &exception, complex dir, complex uv0, Half h);

        void split_segment(
            const Msgmt& target,
            const Msgmt& crash,
            const double ratio
            ) {
            for (auto it = sgmts.begin(); it != sgmts.end(); ++it) {
                if (*it == target) {
                    const bool ccw = cross(crash.diff(), it->diff()) > 0;
                    const auto uv  = lerp(it->fr.uv, it->to.uv, ratio);
                    const auto fr  = Mvert(graph, uv, crash.curv, ccw ? Left : Right, std::nullopt);
                    const auto to  = Mvert(graph, it->to.uv, it->to.crash, it->to.side, it->to.half);
                    it->to = Mvert(graph, uv, crash.curv, ccw ? Right : Left, std::nullopt);
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

    ///
    /// The endpoint of the motorcycle graph system.
    ///
    class MotorcycleGraph {
    public:
        const VecXc& cfn;
        std::vector<Mport> mports;
        std::vector<Mcurv> medges;

        MotorcycleGraph(
            const Hmesh& hmesh,
            const VecXc &cfn,
            const VecXi &matching,
            const VecXi &singular
        ): cfn(cfn)
        {
            gen_ports(hmesh, cfn, singular);
            for (const auto &p: mports) { medges.emplace_back(this, p); }

            for (auto &e: medges) e.add_segment_init(cfn, medges, matching);

            while (rg::any_of(medges, [](auto& e){ return !e.cache.intersected; })) {
                for (auto &s: medges)  s.add_segment_next(cfn, medges, matching);
            }

            for (auto &s: medges) s.post_process();
        }

        void gen_ports(const Hmesh &mesh, const VecXc &cfn, const VecXi &singular);

        const Mcurv& get_medge(const int idx) const { return medges[idx]; }
    };

}

#include "motorcycle.ipp"
#endif
