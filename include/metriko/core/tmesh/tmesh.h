//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_TMESH_H
#define METRIKO_TMESH_H
#include "motorcycle.h"

namespace metriko::tm {
template <class T> using vec = std::vector<T>;
template <class T> using set = std::set<T>;
template <class T> using opt = std::optional<T>;

class Tmesh;

class Telem {
public:
    const Tmesh* tm;
    const int id;
    explicit Telem(const Tmesh* t, const int id): tm(t), id(id) {}
};


class Tsgmt : public Telem {
public:

    Tsgmt(const Tmesh* tm, const int id): Telem(tm, id) {}
};

class TempTedge : public Telem {};

class Tedge : public Telem {
public:
    const std::optional<mc::Msgmt> seg_fr;
    const std::optional<mc::Msgmt> seg_to;
    const vec<Tsgmt> segs;

    Tedge(
        const Tmesh* tm,
        const int id,
        const mc::Msgmt& fr,
        const mc::Msgmt& to
    ): Telem(tm, id), seg_fr(fr), seg_to(to) {}

    Tedge(
        const Tmesh* tm,
        const int id,
        const vec<Tsgmt>& segs
    ): Telem(tm, id), seg_fr(std::nullopt), seg_to(std::nullopt), segs(segs) {}

    [[nodiscard]] int n_segments() const { return seg_to.value().id - seg_fr.value().id + 1; }
    [[nodiscard]] complex uv_fr() const { return seg_fr.value().fr.uv; }
    [[nodiscard]] complex uv_to() const { return seg_to.value().to.uv; }
    [[nodiscard]] mc::MvertType type_fr() const { return seg_fr.value().fr.type; }
    [[nodiscard]] mc::MvertType type_to() const { return seg_to.value().to.type; }

    [[nodiscard]] auto segments() const {
        return seg_fr.value().curv->sgmts | vw::filter([&](auto& s) {
            return s.id >= seg_fr.value().id && s.id <= seg_to.value().id;
        });
    }
};

class Thalf : public Telem {
public:
    int teid;  //
    int twid;  //
    bool cano; // canonical flag

    Thalf(
        const Tmesh* tm,
        const int teid,
        const int id,
        const int twid,
        const bool cano
    ): Telem(tm, id), teid(teid), twid(twid), cano(cano) {}

    bool operator==(const Thalf& rhs) const {
        return id == rhs.id && twid == rhs.twid && teid == rhs.teid && cano == rhs.cano;
    }

    [[nodiscard]] const Thalf& twin() const;
    [[nodiscard]] const Thalf& next() const;
    [[nodiscard]] const Thalf& prev() const;
    [[nodiscard]] const Tedge& edge() const;
    [[nodiscard]] complex uv_fr() const;
    [[nodiscard]] complex uv_to() const;
    [[nodiscard]] complex dif_fr() const { return sg_fr().diff() * (cano ? 1. : -1.); }
    [[nodiscard]] complex dif_to() const { return sg_to().diff() * (cano ? 1. : -1.); }
    [[nodiscard]] mc::Msgmt sg_fr() const;
    [[nodiscard]] mc::Msgmt sg_to() const;
    [[nodiscard]] mc::MvertType type_fr() const;
    [[nodiscard]] mc::MvertType type_to() const;
    [[nodiscard]] vec<Thalf> adj_thalfs() const;
};

class Tquad {
public:
    int id;
    vec<int> thids;
    vec<int> sides;

    Tquad(
        const int id,                 //
        const vec<mc::Mcurv>& mcurvs, //
        const vec<Thalf>& thalfs,     //
        const int bgn_id              // beginning thalf index
    ): id(id) {
        int side_id = 0;
        int curr_id = bgn_id;

        do {
            thids.emplace_back(curr_id);
            sides.emplace_back(side_id);
            auto [th, f] = choose_next_thalf(mcurvs, thalfs, thalfs[curr_id]);
            if (f) side_id = (side_id + 1) % 4;
            curr_id = th.id;
        }
        while (bgn_id != curr_id);

        // sort thalfs as not to start from a middle of side
        if (sides.front() == sides.back()) {
            int i = sides.front();
            int n = (int)rg::distance(sides | vw::take_while([=](int x) { return x == i; }));
            rg::rotate(sides, sides.begin() + n);
            rg::rotate(thids, thids.begin() + n);
        }
        assert(sides.front() != sides.back());
    }

    int find_first_thid(int side) const {
        auto it = rg::find(sides, side);
        return it != sides.end() ? thids[rg::distance(sides.begin(), it)] : -1;
    }

    int find_side(const Thalf& th) const {
        for (int i = 0; i < thids.size(); i++)
            if (th.id == thids[i]) return sides[i];
        return -1;
    }

    std::vector<int> thids_by_side(int side, bool reverse = false) const {
        std::vector<int> result;
        for (int i = 0; i < thids.size(); i++) {
            int j = reverse ? (int)thids.size() - i - 1 : i;
            if (sides[j] == side) result.emplace_back(thids[j]);
        }
        return result;
    }

    static std::pair<Thalf, bool> choose_next_thalf(
        const std::vector<mc::Mcurv>& mcurvs,
        const std::vector<Thalf>& thalfs,
        const Thalf& curr
    );
};

class Tmesh {
public:
    vec<Tquad> tquads;
    vec<Thalf> thalfs;
    vec<Tedge> tedges;
    VecXi th2sing; // -1 if thalf is not from singular, otherwise vertex id
    VecXi th2quad;
    VecXi th2side;
    VecXi th2iter;
    size_t nTQ;
    size_t nTE;
    size_t nTH;

    explicit Tmesh() = default;

    explicit Tmesh(const vec<mc::Mcurv>& mcurvs) {

        for (auto& mc: mcurvs) {
            mc::Msgmt start = mc.sgmts.front();
            for (auto it = mc.sgmts.begin(); it != mc.sgmts.end(); ++it) {
                if (it->to.type != mc::None) {
                    int s = tedges.size();
                    tedges.emplace_back(this, s, start, *it);
                    thalfs.emplace_back(this, s, s * 2 + 0, s * 2 + 1, true);
                    thalfs.emplace_back(this, s, s * 2 + 1, s * 2 + 0, false);
                    if (it->next_id != -1) start = mc.sgmts[it->next_id];
                }
            }
        }

        std::vector visit(thalfs.size(), false);

        while (rg::any_of(visit, [](const bool f) { return !f; })) {
            auto it = rg::find(visit, false);
            auto id = std::distance(visit.begin(), it);
            auto tq = Tquad(tquads.size(), mcurvs, thalfs, id);
            tquads.emplace_back(tq);
            for (int thid: tq.thids) visit[thid] = true;
        }

        nTE = tedges.size();
        nTH = thalfs.size();
        nTQ = tquads.size();

        th2side.resize(nTH);
        th2quad.resize(nTH);
        th2iter.resize(nTH);
        th2sing.resize(nTH);

        for (int i = 0; i < nTQ; i++) {
            const Tquad& tq = tquads[i];
            for (int j = 0; j < tq.thids.size(); j++) {
                th2quad[tq.thids[j]] = i;
                th2side[tq.thids[j]] = tq.sides[j];
                th2iter[tq.thids[j]] = j;
            }
        }

        th2sing.setConstant(-1);
        for (auto& mc: mcurvs) {
            auto it = rg::find_if(thalfs, [&](auto& th) { return th.edge().seg_fr == mc.sgmts.front(); });
            assert(it != thalfs.end());
            th2sing[it->id] = mc.port.vert.id;
        }
    }

    int next_thid(const int curr) const {
        int iQ = th2quad[curr];
        int iT = th2iter[curr];
        auto& tq = tquads[iQ];
        return tq.thids[(iT + 1) % tq.thids.size()];
    }

    int prev_thid(const int curr) const {
        int iQ = th2quad[curr];
        int iT = th2iter[curr];
        auto& tq = tquads[iQ];
        return tq.thids[(iT - 1 + tq.thids.size()) % tq.thids.size()];
    }
};
}

// ipp
namespace metriko::tm {
inline const Tedge& Thalf::edge() const { return tm->tedges[teid]; }
inline const Thalf& Thalf::twin() const { return tm->thalfs[twid]; }
inline const Thalf& Thalf::next() const { return tm->thalfs[tm->next_thid(id)]; }
inline const Thalf& Thalf::prev() const { return tm->thalfs[tm->prev_thid(id)]; }

inline complex Thalf::uv_fr() const { return cano ? edge().uv_fr() : edge().uv_to(); }
inline complex Thalf::uv_to() const { return cano ? edge().uv_to() : edge().uv_fr(); }

inline mc::MvertType Thalf::type_fr() const { return cano ? edge().type_fr() : edge().type_to(); }
inline mc::MvertType Thalf::type_to() const { return cano ? edge().type_to() : edge().type_fr(); }

inline mc::Msgmt Thalf::sg_fr() const { return cano ? edge().seg_fr.value() : edge().seg_to.value(); }
inline mc::Msgmt Thalf::sg_to() const { return cano ? edge().seg_to.value() : edge().seg_fr.value(); }

inline std::vector<Thalf> Thalf::adj_thalfs() const {
    std::vector<Thalf> res;
    auto type = cano ? sg_to().to.type : sg_to().fr.type;
    auto& next = this->next();
    auto& twin = this->twin();
    switch (type) {
    case mc::HitR:
        res.emplace_back(next);
        res.emplace_back(next.twin().next());
        break;
    case mc::HitB:
    case mc::HitL:
        res.emplace_back(next);
        res.emplace_back(twin.prev().twin());
        break;
    default: ;
    }
    return res;
}

inline std::pair<Thalf, bool> Tquad::choose_next_thalf(
    const std::vector<mc::Mcurv>& mcurvs,
    const std::vector<Thalf>& thalfs,
    const Thalf& curr
) {
    auto find_th = [&thalfs](const mc::Msgmt& ms, bool cano) -> Thalf {
        auto it = rg::find_if(thalfs, [&](const Thalf &th) {
            const auto &te = th.edge();
            return th.cano == cano && (te.seg_fr == ms || te.seg_to == ms);
        });
        if (it == thalfs.end()) { throw std::runtime_error("thalf not found"); }
        return *it;
    };

    const Tedge& te = curr.edge();

    if (curr.cano) {
        switch (te.type_to()) {
        case mc::HitB: {
            for (const auto& ms: te.seg_to.value().to.crash->sgmts) {
                if (ms.face.id == te.seg_to.value().face.id) {
                    auto uv2 = te.seg_to.value().to.uv;
                    auto dif = te.seg_to.value().diff();
                    if (equal(ms.fr.uv, uv2) && cross(dif, ms.diff()) > 0)  return {find_th(ms, true) , true};
                    if (equal(ms.to.uv, uv2) && cross(dif, -ms.diff()) > 0) return {find_th(ms, false), true};
                }
            }
            throw std::runtime_error("thalf not found");
        }
        case mc::HitR: {
            const mc::Mcurv* c = te.seg_to.value().to.crash;
            return {find_th(c->sgmts.back(), false), true};
        }
        default: return {thalfs[curr.id + 2], false};
        }
    }

    //--- not cannonical --- //
    if (te.seg_fr == te.seg_fr.value().curv->sgmts.front()) {
        const mc::Mcurv& c = mcurvs[te.seg_fr.value().curv->port.prev];
        return {find_th(c.sgmts.front(), true), true};
    }
    if (te.type_fr() == mc::HitR) {
        const mc::Mcurv* c = te.seg_fr.value().fr.crash;
        return {find_th(c->sgmts.back(), false), true};
    }
    return {thalfs[curr.id - 2], false};
}
}

#endif
