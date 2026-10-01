//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include <format>
#include "emesh.h"
using namespace metriko;

bool Emesh::collapse_tquad_prepare(int tqid, Tqaux& aux) const {
    const auto& tq = tquads[tqid];
    if (tq.id == -1) return false;

    auto zero = [&](int thid) { return thalfs[thid].x == 0; };
    auto sumX = [&](const vec<int>& thids) { double s = 0; for (int t: thids) s += thalfs[t].x; return s; };
    auto sumR = [&](const vec<int>& thids) { double s = 0; for (int t: thids) s += thalfs[t].r; return s; };

    int side = -1;
    if (tq.thids(0).size() == 1 && tq.thids(2).size() == 1 && zero(tq.thids(0).front())) side = 0;
    if (tq.thids(1).size() == 1 && tq.thids(3).size() == 1 && zero(tq.thids(1).front())) side = 1;
    if (side == -1) return false;

    int s_r = side,
        s_t = side + 1,
        s_l = side + 2,
        s_b = (side + 3) % 4;

    if (sumX(tq.thids(s_t)) == 0 || sumX(tq.thids(s_b)) == 0) return false;

    aux.tqid    = tqid;
    aux.thid_r  = tq.thids(s_r).front();
    aux.thid_l  = tq.thids(s_l).front();
    aux.thids_t = tq.thids(s_t);
    aux.thids_b = tq.thids(s_b);
    const auto& th_l = thalfs[aux.thid_l];
    const auto& th_r = thalfs[aux.thid_r];

    auto find_terminal = [&](int thid, bool s0, bool s1) -> std::pair<int, bool> {
        const auto& th0 = thalfs[thid];
        const auto& th1 = thalfs[th0.twid];
        if (th0.bgn) return { th0.nid_fr(), s0 };
        if (th1.bgn) return { th1.nid_fr(), s1 };
        if (th0.end) return { th0.nid_to(), s1 };
        if (th1.end) return { th1.nid_to(), s0 };
        METRIKO_FAIL("zero thalf {} has no terminal node", thid);
    };

    auto is_junction = [&](int nid, int thid_at) {
        return !rg::contains(std::array{ th_l.nid_fr(), th_l.nid_to(), th_r.nid_fr(), th_r.nid_to() }, nid)
            && count_adj_tquads(thid_at) != 2;
    };
    int    oft = 0, btm = 0;
    double rt  = 0, rb  = 0;
    double span   = sumX(aux.thids_t),
           rt_sum = sumR(aux.thids_t),
           rb_sum = sumR(aux.thids_b);
    for (int thid: aux.thids_t | vw::reverse) { const auto& th = thalfs[thid]; oft += th.x; rt += th.r; if (is_junction(th.nid_fr(), th.id))   aux.pts.push_back({ .nid = th.nid_fr(), .val = oft, .ord = rt / rt_sum * span, .top = true  }); }
    for (int thid: aux.thids_b)               { const auto& th = thalfs[thid]; btm += th.x; rb += th.r; if (is_junction(th.nid_to(), th.twid)) aux.pts.push_back({ .nid = th.nid_to(), .val = btm, .ord = rb / rb_sum * span, .top = false }); }

    auto [nid_bgn, side_bgn] = find_terminal(aux.thid_l, true, false);
    auto [nid_end, side_end] = find_terminal(aux.thid_r, false, true);
    aux.pts.push_back({ .nid = nid_bgn, .val = 0,   .ord = 0.,          .top = side_bgn });
    aux.pts.push_back({ .nid = nid_end, .val = oft, .ord = (double)oft, .top = side_end });
    rg::stable_sort(aux.pts, {}, [](const Tqpoint& p) { return std::pair(p.val, p.ord); });

    for (const auto* thids: { &aux.thids_t, &aux.thids_b })
        for (int thid: *thids | vw::filter([&](int i){return zero(i); })) {
            auto p0 = rg::find(aux.pts, thalfs[thid].nid_fr(), &Tqpoint::nid);
            auto p1 = rg::find(aux.pts, thalfs[thid].nid_to(), &Tqpoint::nid);
            if (p0 == aux.pts.end() || p1 == aux.pts.end()) return false;
            if (std::abs(p0 - p1) != 1) return false;
            if (rg::any_of(aux.pts, [&](auto& p) { return p.val == p0->val && p.top != p0->top; })) return false;
        }
    return true;
}
