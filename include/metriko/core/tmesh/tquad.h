//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_TQUAD_H
#define METRIKO_TQUAD_H
#include "thalf.h"

namespace metriko {
    inline std::pair<Thalf, bool> choose_next_thalf(
        const std::vector<Mcurv> &splines,
        const std::vector<Thalf> &thalfs,
        const Thalf &curr
    ) {
        auto find_th = [&thalfs](const Msgmt &seg, bool cannonical) -> Thalf {
            for (const Thalf &th: thalfs) {
                auto &te = th.edge();
                if (th.cannonical == cannonical && (te.seg_fr == seg || te.seg_to == seg)) return th;
            }
            throw std::runtime_error("thalf not found");
        };

        const Tedge &te = curr.edge();

        if (curr.cannonical) {
            if (te.seg_to.to.side == Crash) {
                for (const Msgmt &seg: te.seg_to.to.crash->sgmts) {
                    //Task: fix bad implementation.
                    if (seg.face.id == te.seg_to.face.id) {
                        auto uv2 = te.seg_to.to.uv;
                        auto dif = te.seg_to.diff();
                        if (equal(seg.fr.uv, uv2) && cross(dif, seg.diff()) > 0) return {find_th(seg, true), true};
                        if (equal(seg.to.uv, uv2) && cross(dif, -seg.diff()) > 0) return {find_th(seg, false), true};
                    }
                }
                throw std::runtime_error("thalf not found");
            }
            if (te.seg_to.to.side == Right) {
                const Mcurv *spl = te.seg_to.to.crash;
                return {find_th(spl->sgmts.back(), false), true};
            }
            return {thalfs[curr.id + 2], false};
        }

        //--- not cannonical --- //
        if (te.seg_fr == te.seg_fr.curv->sgmts.front()) {
            const auto &spl = splines[te.seg_fr.curv->port.prev];
            return {find_th(spl.sgmts.front(), true), true};
        }
        if (te.seg_fr.fr.side == Right) {
            const Mcurv *spl = te.seg_fr.fr.crash;
            return {find_th(spl->sgmts.back(), false), true};
        }
        return {thalfs[curr.id - 2], false};
    }

    class Tquad {
    public:
        std::vector<int> thids;
        std::vector<int> sides;

        Tquad(
            const std::vector<Mcurv> &splines,
            const std::vector<Thalf> &thalfs,
            const Thalf &bgn
        ) {
            int side_idx = 0;
            Thalf curr = bgn;

            do {
                thids.emplace_back(curr.id);
                sides.emplace_back(side_idx);
                auto [thalf, flag] = choose_next_thalf(splines, thalfs, curr);
                if (flag) side_idx = (side_idx + 1) % 4;
                curr = thalf;
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
            //auto it = rg::find(sides, side);
            //return it != sides.end() ? thids[rg::distance(sides.begin(), it)] : -1;
            for (int i = 0; i < thids.size(); i++)
                if (sides[i] == side) return thids[i];
            return -1;
        }

        std::vector<int> thids_by_side(int side_id, bool reverse = false) const {
            std::vector<int> result;
            for (int i = 0; i < thids.size(); i++) {
                int j = reverse ? (int) thids.size() - i - 1 : i;
                if (sides[j] == side_id) result.emplace_back(thids[j]);
            }
            return result;
        }
    };
}

#endif
