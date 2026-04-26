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

struct Tmesh;

struct Telem {
    const Tmesh* tm;
    const int id;
    explicit Telem(const Tmesh* t, const int id): tm(t), id(id) {}
};

struct Tvert {
    complex uv;
    opt<std::pair<Half, double>> cut;
    explicit Tvert(const complex uv, opt<std::pair<Half, double>> cut): uv(uv), cut(std::move(cut)) {}
};

struct Tsgmt {
    Face  face;
    Tvert tvFr;
    Tvert tvTo;
    Tsgmt(const Face& f, Tvert fr, Tvert to): face(f), tvFr(std::move(fr)), tvTo(std::move(to)) {}
};

struct Ttemp {
    const int id;
    const mc::Msgmt fr;
    const mc::Msgmt to;
    Ttemp(const int id, const mc::Msgmt fr, const mc::Msgmt to): id(id), fr(fr), to(to) {}
};

struct Tedge : Telem {
    const vec<Tsgmt> segs;
    double len = 0;
    bool isBgn; // a flag to show that the bgn vert is from singular point
    bool isEnd; // a flag to show that the end vert is HitB
    // int sing_vid // todo: better replacing th2sing

    Tedge(
        const Tmesh* tm,
        const int teid,
        const vec<Tsgmt>& segs,
        const bool isBgn,
        const bool isEnd
    ): Telem(tm, teid), segs(segs), isBgn(isBgn), isEnd(isEnd) { for (auto& sg: segs) len += std::abs(sg.tvTo.uv - sg.tvFr.uv); }

    [[nodiscard]] complex uv_fr() const { return segs.front().tvFr.uv; }
    [[nodiscard]] complex uv_to() const { return segs.back().tvTo.uv;  }
};

struct Thalf : Telem {
    int      teid;
    int      twid;
    bool     cano;
    vec<int> adjs;

    Thalf(
        const Tmesh* tm,
        const int teid,
        const int thid,
        const int twid,
        const bool cano
    ): Telem(tm, thid), teid(teid), twid(twid), cano(cano) {}

    bool operator==(const Thalf& rhs) const { return id == rhs.id && twid == rhs.twid && teid == rhs.teid && cano == rhs.cano; }

    [[nodiscard]] const Thalf& twin() const;
    [[nodiscard]] const Thalf& next() const;
    [[nodiscard]] const Thalf& prev() const;
    [[nodiscard]] const Tedge& edge() const;
    [[nodiscard]] complex uv_fr() const;
    [[nodiscard]] complex uv_to() const;
};

struct Tquad {
    int id;
    vec<int> thids;
    vec<int> sides;

    Tquad(
        int id,
        const vec<mc::Mcurv>& mcurvs,
        const vec<Ttemp>& ttemps,
        const vec<Thalf>& thalfs,
        const int bgn_id
    ) : id(id) {

        auto find_th = [&](const mc::Msgmt& ms, bool cano) -> int {
            auto it = rg::find_if(thalfs, [&](auto& th) { return th.cano == cano && (ttemps[th.teid].fr == ms || ttemps[th.teid].to == ms); });
            assert(it != thalfs.end());
            return it->id;
        };

        auto get_next = [&](int thid) -> std::pair<int, bool> {
            const auto& th = thalfs[thid];
            const auto& tt = ttemps[th.teid];
            const auto& ct = tt.to.curv;
            assert(tt.fr.curv == tt.to.curv);

            if (th.cano) {
                auto& tov = tt.to.to;
                auto& sgs = tov.crash->sgmts;
                switch (tov.type) {
                case mc::HitR: { return {find_th(sgs.back(), false), true}; }
                case mc::HitB: {
                    for (auto& s: sgs) {
                        if (s.fr.crash == ct && s.fr.type == mc::HitL) return {find_th(s, true),  true};
                        if (s.to.crash == ct && s.to.type == mc::HitL) return {find_th(s, false), true};
                    }
                    throw std::runtime_error("thalf not found");
                }
                default: return {thid + 2, false};
                }
            }

            // not cannonical case
            if (tt.fr == ct->sgmts.front()) { return {find_th(mcurvs[ct->port.prev].sgmts.front(), true), true}; }
            if (tt.fr.fr.type == mc::HitR)  { return {find_th(tt.fr.fr.crash->sgmts.back()      , false), true}; }
            return {thid - 2, false};
        };

        int side_id = 0;
        int curr_id = bgn_id;

        do {
            thids.emplace_back(curr_id);
            sides.emplace_back(side_id);
            auto [thid, f] = get_next(curr_id);
            if (f) side_id = (side_id + 1) % 4;
            curr_id = thid;
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

    [[nodiscard]] vec<int> thids_by_side(int side) const {
        return vw::zip(thids, sides)
            | vw::filter([&](const auto& p) { return std::get<1>(p) == side; })
            | vw::elements<0>
            | rg::to<vec<int>>();
    }
};

struct Tmesh {
    vec<Tquad> tquads;
    vec<Thalf> thalfs;
    vec<Tedge> tedges;
    vec<Ttemp> ttemps;
    VecXi th2sing; // todo: remove. -1 if thalf is not from singular, otherwise vertex id
    VecXi th2quad;
    VecXi th2side;
    VecXi th2iter;
    size_t nTQ;
    size_t nTE;
    size_t nTH;

    explicit Tmesh(const vec<mc::Mcurv>& mcurvs) {
        for (auto& mc: mcurvs) {
            int bgnIdx = mc.sgmts.front().id;
            for (auto it = mc.sgmts.begin(); it != mc.sgmts.end(); ++it) {
                if (it->to.type != mc::None) {
                    int s = (int)ttemps.size();
                    ttemps.emplace_back(s, mc.sgmts[bgnIdx], *it);
                    thalfs.emplace_back(this, s, s * 2, s * 2 + 1, true);
                    thalfs.emplace_back(this, s, s * 2 + 1, s * 2, false);
                    if (it->next_id != -1) bgnIdx = mc.sgmts[it->next_id].id;
                }
            }
        }

        // assign tedges from ttemps
        for (auto& tt: ttemps) {
            auto sgs = tt.fr.curv->sgmts
                | vw::drop(tt.fr.id)
                | vw::take(tt.to.id - tt.fr.id + 1)
                | vw::transform([](auto& s) { return Tsgmt(s.face, Tvert(s.fr.uv, s.fr.cut), Tvert(s.to.uv, s.to.cut)); })
                | rg::to<vec<Tsgmt>>();
            bool isBgn = tt.fr.fr.type == mc::First;
            bool isEnd = tt.to.to.type == mc::HitB;
            tedges.emplace_back(this, tt.id, sgs, isBgn, isEnd);
        }


        th2side.resize((int)thalfs.size());
        th2quad.resize((int)thalfs.size());
        th2iter.resize((int)thalfs.size());
        th2sing.resize((int)thalfs.size());
        std::vector visit(thalfs.size(), false);

        while (rg::any_of(visit, [](const bool f) { return !f; })) {
            auto it = rg::find(visit, false);
            auto id = std::distance(visit.begin(), it);
            auto tq = Tquad((int)tquads.size(), mcurvs, ttemps, thalfs, id);
            tquads.emplace_back(tq);
            for (int i: tq.thids) visit[i] = true;
            for (int j = 0; j < tq.thids.size(); j++) {
                th2quad[tq.thids[j]] = tq.id;
                th2side[tq.thids[j]] = tq.sides[j];
                th2iter[tq.thids[j]] = j;
            }
        }

        nTE = tedges.size();
        nTH = thalfs.size();
        nTQ = tquads.size();

        th2sing.setConstant(-1);
        for (auto& mc: mcurvs) {
            auto it = rg::find_if(thalfs, [&](auto& th) { return ttemps[th.teid].fr == mc.sgmts.front(); });
            assert(it != thalfs.end());
            th2sing[it->id] = mc.port.vert.id;
        }

        // Assigns adjacent thalfs for each thalf
        for (auto& th: thalfs) {
            auto& tt = ttemps[th.teid];
            switch (th.cano ? tt.to.to.type : tt.fr.fr.type) {
                case mc::HitB:
                case mc::HitL: th.adjs = vec { th.next().id, th.twin().prev().twin().id }; break;
                case mc::HitR: th.adjs = vec { th.next().id, th.next().twin().next().id }; break;
                default: break;
            }
        }
    }

    [[nodiscard]] int next_thid(const int i) const { int iQ = th2quad[i]; int iT = th2iter[i]; auto& tq = tquads[iQ]; return tq.thids[(iT + 1) % tq.thids.size()]; }
    [[nodiscard]] int prev_thid(const int i) const { int iQ = th2quad[i]; int iT = th2iter[i]; auto& tq = tquads[iQ]; return tq.thids[(iT - 1 + tq.thids.size()) % tq.thids.size()]; }
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
}

#endif
