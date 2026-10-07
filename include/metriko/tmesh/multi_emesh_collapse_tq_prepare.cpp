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

bool Emesh::collapse_tquad_chain_prepare(int tqid, Tqchain& chain) const {
    auto& tq = tquads[tqid];
    if (tq.id == -1) return false;

    auto zero = [&](int thid) { return thalfs[thid].x == 0; };
    auto sumX = [&](const vec<int>& thids) { double s = 0; for (int t: thids) s += thalfs[t].x; return s; };
    auto sumR = [&](const vec<int>& thids) { double s = 0; for (int t: thids) s += thalfs[t].r; return s; };

    auto find_simple_chain = [&](vec<int>& seq) {
        while (true) {
            const auto& th_prev = thalfs[seq.back()];
            const auto& th_curr = thalfs[th_prev.twid];
            METRIKO_CHECK(th_curr.id != -1 && th_curr.tqid != -1, "thalf {}: twin {} is not on a live tquad", th_prev.id, th_prev.twid);
            const auto& tq_curr = tquads[th_curr.tqid];
            METRIKO_CHECK(tq_curr.id != -1 && tq_curr.find_of(th_curr.id) != tq_curr.data.end(), "thalf {} is not on its tquad {}", th_curr.id, th_curr.tqid);
            auto s0      = tq_curr.side_of(th_curr);
            auto thids_0 = tq_curr.thids(s0);
            auto thids_1 = tq_curr.thids((s0 + 1) % 4);
            auto thids_2 = tq_curr.thids((s0 + 2) % 4);
            auto thids_3 = tq_curr.thids((s0 + 3) % 4);
            if (!rg::all_of(thids_0, zero))               return true; // stop caz tq_curr is not collapsable
            if (thids_0.size() > 1 || thids_2.size() > 1) return true; // stop caz tq_curr has split chain
            if (n_adj_tqs(th_curr.id)   == 4)             return true; // stop caz cannot go father
            if (n_adj_tqs(th_curr.twid) == 4)             return true; // stop caz cannot go father
            if (sumX(thids_1) == 0 || sumX(thids_3) == 0) return true; // stop caz tq_curr's all sides are x == zero
            METRIKO_CHECK(!thids_2.empty(), "tquad {}: no thalf on the side opposite to thalf {} (side sizes {} {} {} {})", tq_curr.id, th_curr.id, thids_0.size(), thids_1.size(), thids_2.size(), thids_3.size());
            seq.push_back(thids_2.front());
        }
    };

    int side = -1;
    if (tq.thids(0).size() == 1 && tq.thids(2).size() == 1 && zero(tq.thids(0).front())) side = 0;
    if (tq.thids(1).size() == 1 && tq.thids(3).size() == 1 && zero(tq.thids(1).front())) side = 1;
    if (side == -1) return false;

    if (sumX(tq.thids(side == 0 ? 3 : 0)) == 0) return false;
    if (sumX(tq.thids(side == 0 ? 1 : 2)) == 0) return false;

    vec seq_l = {tq.thids(side == 0 ? 2 : 3)[0]};
    vec seq_r = {tq.thids(side == 0 ? 0 : 1)[0]};
    if (!find_simple_chain(seq_l)) return false;
    if (!find_simple_chain(seq_r)) return false;

    for (int thid: seq_l | vw::reverse) chain.thids_z.push_back(thalfs[thid].twid);
    for (int thid: seq_r)               chain.thids_z.push_back(thid);
    chain.thid_l = seq_l.back();
    chain.thid_r = seq_r.back();

    auto find_terminal = [&](int thid, bool s0, bool s1) -> std::pair<int, bool> {
        const auto& th0 = thalfs[thid];
        const auto& th1 = thalfs[th0.twid];
        if (th0.bgn) return { th0.nid_fr(), s0 };
        if (th1.bgn) return { th1.nid_fr(), s1 };
        if (th0.end) return { th0.nid_to(), s1 };
        if (th1.end) return { th1.nid_to(), s0 };
        METRIKO_FAIL("chain end thalf {} has no terminal node", thid);
    };

    // 1: push intermidiate pts
    auto oft = 0;
    for (size_t i = 1; i < chain.thids_z.size(); ++i) {
        const auto& th_l = thalfs[thalfs[chain.thids_z[i - 1]].twid];
        const auto& th_r = thalfs[chain.thids_z[i]];
        const auto& th_t = thalfs[th_r.twid];
        const auto& tq_  = tquads[th_r.tqid];
        auto n0 = th_r.nid_to();
        auto n1 = th_r.nid_fr();
        auto sr = tq_.side_of(th_r);
        auto thids_t = tq_.thids((sr + 1) % 4);
        auto thids_b = tq_.thids((sr + 3) % 4);

        chain.tqids.push_back(th_r.tqid);
        chain.thids_t.insert(chain.thids_t.end(), thids_t.begin(), thids_t.end());
        chain.thids_b.insert(chain.thids_b.end(), thids_b.begin(), thids_b.end());

        auto f = [&](int nid, int thid_at) {
            return !rg::contains(
                std::array{
                th_l.nid_fr(),
                th_l.nid_to(),
                th_r.nid_fr(),
                th_r.nid_to()
            }, nid) && n_adj_tqs(thid_at) != 2; // adjacency around the pushed node
        };

        auto   btm  = oft;
        double base = oft;
        double rt = 0;
        double rb = 0;
        double span   = sumX(thids_t);
        double rt_sum = sumR(thids_t);
        double rb_sum = sumR(thids_b);
        for (int thid: thids_t | vw::reverse) { auto& th = thalfs[thid]; oft += th.x; rt += th.r; if (f(th.nid_fr(), th.id))   chain.pts.push_back({ .nid = th.nid_fr(), .val = oft, .ord = base + rt / rt_sum * span, .top = true  }); }
        for (int thid: thids_b)               { auto& th = thalfs[thid]; btm += th.x; rb += th.r; if (f(th.nid_to(), th.twid)) chain.pts.push_back({ .nid = th.nid_to(), .val = btm, .ord = base + rb / rb_sum * span, .top = false }); }

        chain.bounds.push_back(oft);

        // push ladder thalf points (right side of tquad)
        if (i == chain.thids_z.size() - 1) break;
        if      (th_r.bgn) chain.pts.push_back({.nid = n1, .val = oft, .ord = (double)oft, .top = false });
        else if (th_t.bgn) chain.pts.push_back({.nid = n0, .val = oft, .ord = (double)oft, .top = true  });
        else if (th_r.end) chain.pts.push_back({.nid = n0, .val = oft, .ord = (double)oft, .top = true  });
        else if (th_t.end) chain.pts.push_back({.nid = n1, .val = oft, .ord = (double)oft, .top = false });
    }

    // the band folds onto itself around a tedge (a hairpin): both halves of that tedge are lateral thalfs of the
    // chain, and the merge would rewrite a tquad of the chain as if it were an outer one. not handled yet
    // for (const auto* thids: { &chain.thids_t, &chain.thids_b })
    // for (int thid: *thids)
    //     if (rg::contains(chain.thids_t, thalfs[thid].twid) ||
    //         rg::contains(chain.thids_b, thalfs[thid].twid)) return false;

    // 2: push the edge pts
    auto [nid_bgn, side_bgn] = find_terminal(chain.thid_l, true, false);
    auto [nid_end, side_end] = find_terminal(chain.thid_r, false, true);
    chain.pts.push_back({ .nid = nid_bgn, .val = 0,   .ord = 0.,          .top = side_bgn });
    chain.pts.push_back({ .nid = nid_end, .val = oft, .ord = (double)oft, .top = side_end });

    rg::stable_sort(chain.pts, {}, [](const Tqpoint& p) { return std::pair(p.val, p.ord); });

    // A zero thalf on a lateral side (the end of another band) stays as it is on the merged tedge, which needs
    // its two end points next to each other on one side. give up when an end point was not pushed, when something
    // falls between them, or when the other side has a point at the same position
    for (const auto* thids: { &chain.thids_t, &chain.thids_b }) {
        for (int thid: *thids) {
            if (!zero(thid)) continue;
            auto p0 = rg::find(chain.pts, thalfs[thid].nid_fr(), &Tqpoint::nid);
            auto p1 = rg::find(chain.pts, thalfs[thid].nid_to(), &Tqpoint::nid);
            if (p0 == chain.pts.end() || p1 == chain.pts.end()) return false;
            if (std::abs(p0 - p1) != 1) return false;
            if (rg::any_of(chain.pts, [&](const Tqpoint& p) { return p.val == p0->val && p.top != p0->top; })) return false;
        }
    }

    return true;
}
