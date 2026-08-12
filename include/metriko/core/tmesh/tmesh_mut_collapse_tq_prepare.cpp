#include "./tmesh_mut.h"
using namespace metriko;

bool find_simple_chain(const TmeshMut& tm, vec<int>& seq) {
    auto zero = [&](int thid) { return tm.thalfs[thid].x == 0; };
    while (true) {
        const auto& th_prev = tm.thalfs[seq.back()];
        const auto& th_curr = tm.thalfs[th_prev.twid];
        const auto& tq_curr = tm.tquads[th_curr.tqid];
        auto s0     = tq_curr.side_of(th_curr);
        auto thids_0 = tq_curr.thids(s0);
        auto thids_1 = tq_curr.thids((s0 + 1) % 4);
        auto thids_2 = tq_curr.thids((s0 + 2) % 4); // opposite
        auto thids_3 = tq_curr.thids((s0 + 3) % 4);
        if (thids_2.size() > 1 || thids_0.size() > 1) return true; // stop simple chain

        auto th_pair = tm.thalfs[thids_2.front()];
        auto a_curr_0 = tm.count_adj_tquads(th_curr.id); // top
        auto a_curr_1 = tm.count_adj_tquads(th_curr.twid); // top
        if (a_curr_0 == 4 || a_curr_1 == 4) return true;

        if (th_pair.x == 0 && rg::any_of(thids_1, zero)) { return false; }
        if (th_pair.x == 0 && rg::any_of(thids_3, zero)) { return false; }
        seq.push_back(th_pair.id);
    }
}

bool TmeshMut::collapse_tquad_chain_prepare(int tqid, Tqchain& chain) const {
    auto& tq = tquads[tqid];
    int side = -1;
    if (tq.id == -1) return false;
    if (tq.thids(0).size() == 1 && tq.thids(2).size() == 1 && thalfs[tq.thids(0).front()].x == 0) side = 0;
    if (tq.thids(1).size() == 1 && tq.thids(3).size() == 1 && thalfs[tq.thids(1).front()].x == 0) side = 1;
    if (side == -1) return false;

    if (tq.thids(side == 0 ? 2 : 3).size() != 1) return false;
    if (tq.thids(side == 0 ? 0 : 1).size() != 1) return false;
    for (int thid : tq.thids(side == 0 ? 3 : 0)) { if (thalfs[thid].x == 0) return false; }
    for (int thid : tq.thids(side == 0 ? 1 : 2)) { if (thalfs[thid].x == 0) return false; }

    vec seq_l = {tq.thids(side == 0 ? 2 : 3)[0]};
    vec seq_r = {tq.thids(side == 0 ? 0 : 1)[0]};
    if (!find_simple_chain(*this, seq_l)) { return false; }
    if (!find_simple_chain(*this, seq_r)) { return false; }

    for (int thid: seq_l | vw::reverse) chain.thids_z.push_back(thalfs[thid].twid);
    for (int thid: seq_r)               chain.thids_z.push_back(thid);
    chain.thid_l = seq_l.back();
    chain.thid_r = seq_r.back();

    for (size_t i = 1; i < chain.thids_z.size(); ++i) {
        const auto& th_ = thalfs[chain.thids_z[i]];
        const auto& tq_ = tquads[th_.tqid];
        auto side_0  = tq_.side_of(th_);
        auto thids_1 = tq_.thids((side_0 + 1) % 4);
        auto thids_3 = tq_.thids((side_0 + 3) % 4);
        chain.tqids.push_back(th_.tqid);
        chain.thids_t.insert(chain.thids_t.end(), thids_1.begin(), thids_1.end());
        chain.thids_b.insert(chain.thids_b.end(), thids_3.begin(), thids_3.end());
    }

    auto find_terminal = [&](int thid, bool s0, bool s1) -> std::pair<const HmLoc&, bool> {
        const auto& th0 = thalfs[thid];
        const auto& th1 = thalfs[th0.twid];
        if (th0.bgn) return { th0.loc_fr(), s0 };
        if (th1.bgn) return { th1.loc_fr(), s1 };
        if (th0.end) return { th0.loc_to(), s1 };
        if (th1.end) return { th1.loc_to(), s0 };
        throw std::runtime_error("error in find_terminal (chain)");
    };
    auto [loc_bgn, side_bgn] = find_terminal(chain.thid_l, true, false);
    auto [loc_end, side_end] = find_terminal(chain.thid_r, false, true);

    int s = 0; for (int t : chain.thids_t) s += thalfs[t].x;
    chain.pts.push_back({ .loc = loc_bgn, .val = 0, .ord = 0.,        .adj = count_adj_tquads(chain.thid_l), .top = side_bgn });
    chain.pts.push_back({ .loc = loc_end, .val = s, .ord = (double)s, .adj = count_adj_tquads(chain.thid_r), .top = side_end });

    auto oft = 0;

    for (size_t i = 1; i < chain.thids_z.size(); ++i) {
        const auto& th_l   = thalfs[thalfs[chain.thids_z[i - 1]].twid];
        const auto& th_r   = thalfs[chain.thids_z[i]];
        const auto& th_twn = thalfs[th_r.twid];
        const auto& tq_    = tquads[th_r.tqid];
        auto a0 = count_adj_tquads(th_twn.id); // top
        auto a1 = count_adj_tquads(th_r.id);   // btm
        auto l0 = th_r.loc_to();
        auto l1 = th_r.loc_fr();
        auto sr = tq_.side_of(th_r);
        auto thids_t = tq_.thids((sr + 1) % 4);
        auto thids_b = tq_.thids((sr + 3) % 4);

        auto btm = oft;
        auto f = [&](const HmLoc& l, int thid_at) {
            bool f1 = l != th_l.loc_fr() && l != th_l.loc_to() && l != th_r.loc_fr() && l != th_r.loc_to();
            bool f2 = count_adj_tquads(thid_at) != 2; // adjacency around the pushed node
            return f1 && f2;
        };

        double base = oft;
        double rt = 0;
        double rb = 0;
        int    span   = 0; for (int t: thids_t) span   += thalfs[t].x;
        double rt_sum = 0; for (int t: thids_t) rt_sum += thalfs[t].r;
        double rb_sum = 0; for (int t: thids_b) rb_sum += thalfs[t].r;
        for (int thid: thids_t | vw::reverse) { auto& th = thalfs[thid]; oft += th.x; rt += th.r; if (f(th.loc_fr(), th.id))   chain.pts.push_back({ .loc = th.loc_fr(), .val = oft, .ord = base + rt / rt_sum * span, .adj = 3, .top = true  }); }
        for (int thid: thids_b)               { auto& th = thalfs[thid]; btm += th.x; rb += th.r; if (f(th.loc_to(), th.twid)) chain.pts.push_back({ .loc = th.loc_to(), .val = btm, .ord = base + rb / rb_sum * span, .adj = 3, .top = false }); }

        chain.bounds.push_back(oft);

        // push ladder thalf points
        if (i == chain.thids_z.size() - 1) break;
        if      (th_r.bgn)   chain.pts.push_back({.loc = l1, .val = oft, .ord = (double)oft, .adj = a1, .top = false });
        else if (th_twn.bgn) chain.pts.push_back({.loc = l0, .val = oft, .ord = (double)oft, .adj = a0, .top = true  });
        else if (th_r.end)   chain.pts.push_back({.loc = l0, .val = oft, .ord = (double)oft, .adj = a0, .top = true  });
        else if (th_twn.end) chain.pts.push_back({.loc = l1, .val = oft, .ord = (double)oft, .adj = a1, .top = false });
    }

    rg::stable_sort(chain.pts, {}, [](const Tqpoint& p) { return std::pair(p.val, p.ord); });

    return true;
}
