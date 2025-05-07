//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_TQUAD_H
#define METRIKO_TQUAD_H
#include "thalf.h"

namespace metriko {
    inline std::pair<Thalf, bool> choose_next_thalf(
        const std::vector<Mcurv> &mcurvs,
        const std::vector<Thalf> &thalfs,
        const Thalf &curr
    ) {
        auto find_th = [&thalfs]( const Msgmt &ms, bool cannonical) -> Thalf {
            for (const Thalf &th: thalfs) {
                auto &te = th.edge();
                if (th.cannonical == cannonical && (te.seg_fr == ms || te.seg_to == ms)) return th;
            }
            throw std::runtime_error("thalf not found");
        };

        const Tedge &te = curr.edge();

        if (curr.cannonical) {
            if (te.type_to() == HitB) {
                for (const Msgmt &ms: te.seg_to.to.crash->sgmts) {
                    if (ms.face.id == te.seg_to.face.id) {
                        auto uv2 = te.seg_to.to.uv;
                        auto dif = te.seg_to.diff();
                        if (equal(ms.fr.uv, uv2) && cross(dif,  ms.diff()) > 0) return {find_th(ms, true),  true};
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

    class Tquad {
    public:
        int id;
        std::vector<int> thids;
        std::vector<int> sides;

        Tquad(
            const int id,
            const std::vector<Mcurv> &mcurvs,
            const std::vector<Thalf> &thalfs,
            const Thalf &bgn
        ) : id(id) {
            int side_idx = 0;
            Thalf curr = bgn;

            do {
                thids.emplace_back(curr.id);
                sides.emplace_back(side_idx);
                auto [th, f] = choose_next_thalf(mcurvs, thalfs, curr);
                if (f) side_idx = (side_idx + 1) % 4;
                curr = th;
            } while (bgn != curr);

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
    };
}

#endif
