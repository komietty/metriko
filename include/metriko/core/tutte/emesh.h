#ifndef TUTTE_H_EMESH_H
#define TUTTE_H_EMESH_H

#include "metriko/core/hmesh/hmesh.h"
#include "tutte_cutting.h"

namespace metriko::tutte {
struct Emesh;

struct Ehalf {
    const Emesh* em;
    vec<Half> halfs;
    int twid;
    double x;
    double r;

    Ehalf(
        const Emesh* em,
        const vec<Half>& halfs,
        const double x,
        const double r,
        const int twin
    ): em(em), halfs(halfs), x(x), r(r), twid(twin) {}

    const Ehalf& twin() const;
};
struct Equad {
    const Emesh* em;
    int id;
    vec<int> ehids;
    vec<int> sides;
    Equad(
        const Emesh* em,
        const int id,
        const vec<int>& ehids,
        const vec<int>& sides
    ): em(em), id(id), ehids(ehids), sides(sides) {}
};

struct Emesh {
    const Hmesh& hm;
    vec<Equad> equads;
    vec<Ehalf> ehalfs;

    explicit Emesh(
        const Hmesh& hm,
        const tm::Tmesh& tm,
        const set<HalfData>& hdata,
        const VecXd& X,
        const VecXd& R
    ): hm(hm) {
        equads.reserve(tm.tquads.size());
        ehalfs.reserve(tm.thalfs.size());
        for (const auto& tq: tm.tquads) {
            equads[tq.id] = Equad(this, tq.id, tq.thids, tq.sides);
            auto rg_tq = rg::equal_range(hdata, tq.id, {}, &HalfData::tqid);
            for (int thid: tq.thids) {
                auto ehs = rg_tq
                    | vw::filter([&](auto& hd) { return thid == hd.thid; })
                    | vw::transform([](auto& hd) { return hd.half; })
                    | rg::to<vec<Half>>();
                ehalfs[thid] = Ehalf(this, ehs, X[thid], R[thid], tm.thalfs[thid].twid);
            }
        }
    }
};

inline const Ehalf& Ehalf::twin() const { return em->ehalfs[twid]; }

}

#endif