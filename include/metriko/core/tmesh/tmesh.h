//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_TMESH_H
#define METRIKO_TMESH_H
#include "motorcycle.h"

namespace metriko {
class Tmesh;

class Telem {
public:
    const Tmesh *tm;
    const int id;
    explicit Telem(const Tmesh *t, const int id): tm(t), id(id) { }
};

// todo: must remove the reference to the motorcycle graph
class Tedge : public Telem {
public:
    const Msgmt seg_fr;
    const Msgmt seg_to;

    Tedge(
        const Tmesh *tm,
        const int id,
        const Msgmt &fr,
        const Msgmt &to
    ): Telem(tm, id), seg_fr(fr), seg_to(to) {
    }

    [[nodiscard]] int n_segments() const { return seg_to.id - seg_fr.id + 1; }
    [[nodiscard]] complex uv_fr() const { return seg_fr.fr.uv; }
    [[nodiscard]] complex uv_to() const { return seg_to.to.uv; }
    [[nodiscard]] MvertType type_fr() const { return seg_fr.fr.type; }
    [[nodiscard]] MvertType type_to() const { return seg_to.to.type; }

    auto segments() const {
        return seg_fr.curv->sgmts | vw::filter([&](auto &s) {
            return s.id >= seg_fr.id && s.id <= seg_to.id;
        });
    }
};

class Thalf : public Telem {
public:
    int teid;  //
    int twid;  //
    bool cano; // canonical flag

    Thalf(
        const Tmesh *tm,
        const int teid,
        const int id,
        const int twid,
        const bool cano
    ): Telem(tm, id), teid(teid), twid(twid), cano(cano) {
    }

    bool operator==(const Thalf &rhs) const {
        return id == rhs.id && twid == rhs.twid && teid == rhs.teid && cano == rhs.cano;
    }

    const Thalf &twin() const;
    const Thalf &next() const;
    const Thalf &prev() const;
    const Tedge &edge() const;
    complex uv_fr() const;
    complex uv_to() const;
    complex dif_fr() const { return sg_fr().diff() * (cano ? 1. : -1.); }
    complex dif_to() const { return sg_to().diff() * (cano ? 1. : -1.); }
    Msgmt sg_fr() const;
    Msgmt sg_to() const;
    std::vector<Thalf> adj_thalfs() const;
};

class Tquad {
public:
    int id;
    std::vector<int> thids;
    std::vector<int> sides;

    Tquad(
        const int id,                     //
        const std::vector<Mcurv> &mcurvs, //
        const std::vector<Thalf> &thalfs, //
        const int bgn_id                  // beginning thalf index
    ) : id(id) {
        int side_id = 0;
        int curr_id = bgn_id;

        do {
            thids.emplace_back(curr_id);
            sides.emplace_back(side_id);
            auto [th, f] = choose_next_thalf(mcurvs, thalfs, thalfs[curr_id]);
            if (f) side_id = (side_id + 1) % 4;
            curr_id = th.id;
        } while (bgn_id != curr_id);

        // sort thalfs as not to start from a middle of side
        if (sides.front() == sides.back()) {
            int i = sides.front();
            int n = rg::distance(sides | vw::take_while([=](int x) { return x == i; }));
            rg::rotate(sides, sides.begin() + n);
            rg::rotate(thids, thids.begin() + n);
        }
        assert(sides.front() != sides.back());
    }

    int find_first_thid(int side) const {
        auto it = rg::find(sides, side);
        return it != sides.end() ? thids[rg::distance(sides.begin(), it)] : -1;
    }

    int find_side(const Thalf &th) const {
        for (int i = 0; i < thids.size(); i++)
            if (th.id == thids[i]) return sides[i];
        return -1;
    }

    std::vector<int> thids_by_side(int side, bool reverse = false) const {
        std::vector<int> result;
        for (int i = 0; i < thids.size(); i++) {
            int j = reverse ? (int) thids.size() - i - 1 : i;
            if (sides[j] == side) result.emplace_back(thids[j]);
        }
        return result;
    }

    static std::pair<Thalf, bool> choose_next_thalf(
        const std::vector<Mcurv> &mcurvs,
        const std::vector<Thalf> &thalfs,
        const Thalf &curr
    );
};

class Tmesh {
public:
    std::vector<Tquad> tquads;
    std::vector<Thalf> thalfs;
    std::vector<Tedge> tedges;
    VecXi th2sing; // -1 if thalf is not from singular, otherwise vertex id
    VecXi th2quad;
    VecXi th2side;
    VecXi th2iter;
    size_t nTQ;
    size_t nTE;
    size_t nTH;

    explicit Tmesh() = default;

    explicit Tmesh(const std::vector<Mcurv> &mcurvs) {
        for (auto &mc: mcurvs) {
            Msgmt start = mc.sgmts.front();
            for (auto it = mc.sgmts.begin(); it != mc.sgmts.end(); ++it) {
                if (it->to.type != None) {
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
            const Tquad &tq = tquads[i];
            for (int j = 0; j < tq.thids.size(); j++) {
                th2quad[tq.thids[j]] = i;
                th2side[tq.thids[j]] = tq.sides[j];
                th2iter[tq.thids[j]] = j;
            }
        }

        th2sing.setConstant(-1);
        for (auto &mc: mcurvs) {
            auto it = rg::find_if(thalfs, [&](auto &th) { return th.edge().seg_fr == mc.sgmts.front(); });
            assert(it != thalfs.end());
            th2sing[it->id] = mc.port.vert.id;
        }
    }

    int next_thid(const int curr) const {
        int iQ = th2quad[curr];
        int iT = th2iter[curr];
        auto &tq = tquads[iQ];
        return tq.thids[(iT + 1) % tq.thids.size()];
    }

    int prev_thid(const int curr) const {
        int iQ = th2quad[curr];
        int iT = th2iter[curr];
        auto &tq = tquads[iQ];
        return tq.thids[(iT - 1 + tq.thids.size()) % tq.thids.size()];
    }
};
}

// ipp
namespace metriko {
inline const Tedge &Thalf::edge() const { return tm->tedges[teid]; }
inline const Thalf &Thalf::twin() const { return tm->thalfs[twid]; }
inline const Thalf &Thalf::next() const { return tm->thalfs[tm->next_thid(id)]; }
inline const Thalf &Thalf::prev() const { return tm->thalfs[tm->prev_thid(id)]; }

inline complex Thalf::uv_fr() const { return cano ? edge().uv_fr() : edge().uv_to(); }
inline complex Thalf::uv_to() const { return cano ? edge().uv_to() : edge().uv_fr(); }

inline Msgmt Thalf::sg_fr() const { return cano ? edge().seg_fr : edge().seg_to; }
inline Msgmt Thalf::sg_to() const { return cano ? edge().seg_to : edge().seg_fr; }

inline std::vector<Thalf> Thalf::adj_thalfs() const {
    std::vector<Thalf> res;
    auto type = cano ? sg_to().to.type : sg_to().fr.type;
    auto &next = this->next();
    auto &twin = this->twin();
    switch (type) {
        case HitR:
            res.emplace_back(next);
            res.emplace_back(next.twin().next());
            break;
        case HitB:
        case HitL:
            res.emplace_back(next);
            res.emplace_back(twin.prev().twin());
            break;
        default: ;
    }
    return res;
}

inline std::pair<Thalf, bool> Tquad::choose_next_thalf(
    const std::vector<Mcurv> &mcurvs,
    const std::vector<Thalf> &thalfs,
    const Thalf &curr
) {
    auto find_th = [&thalfs](const Msgmt &ms, bool cano) -> Thalf {
        for (const Thalf &th: thalfs) {
            auto &te = th.edge();
            if (th.cano == cano && (te.seg_fr == ms || te.seg_to == ms)) return th;
        }
        throw std::runtime_error("thalf not found");
    };

    const Tedge &te = curr.edge();

    if (curr.cano) {
        if (te.type_to() == HitB) {
            for (const Msgmt &ms: te.seg_to.to.crash->sgmts) {
                if (ms.face.id == te.seg_to.face.id) {
                    auto uv2 = te.seg_to.to.uv;
                    auto dif = te.seg_to.diff();
                    if (equal(ms.fr.uv, uv2) && cross(dif,  ms.diff()) > 0) return {find_th(ms,  true), true};
                    if (equal(ms.to.uv, uv2) && cross(dif, -ms.diff()) > 0) return {find_th(ms, false), true};
                }
            }
            throw std::runtime_error("thalf not found");
        }
        if (te.type_to() == HitR) {
            const Mcurv *c = te.seg_to.to.crash;
            return {find_th(c->sgmts.back(), false), true};
        }
        return {thalfs[curr.id + 2], false};
    }

    //--- not cannonical --- //
    if (te.seg_fr == te.seg_fr.curv->sgmts.front()) {
        const Mcurv &c = mcurvs[te.seg_fr.curv->port.prev];
        return {find_th(c.sgmts.front(), true), true};
    }
    if (te.type_fr() == HitR) {
        const Mcurv *c = te.seg_fr.fr.crash;
        return {find_th(c->sgmts.back(), false), true};
    }
    return {thalfs[curr.id - 2], false};
}
}

#endif
