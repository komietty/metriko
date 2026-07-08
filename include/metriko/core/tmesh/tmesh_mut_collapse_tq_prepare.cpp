#include "./tmesh_mut.h"
using namespace metriko;

bool TmeshMut::collapse_tquad_chain_prepare(int tqid, Tqchain& chain) const {
    auto& tq = tquads[tqid];
    int side = -1;
    if (tq.thids(0).size() == 1 && tq.thids(2).size() == 1 && thalfs[tq.thids(0).front()].x == 0) side = 0;
    if (tq.thids(1).size() == 1 && tq.thids(3).size() == 1 && thalfs[tq.thids(1).front()].x == 0) side = 1;
    if (side == -1) return false;

    auto find_simple_chain = [&](vec<int>& seq) {
        auto zero = [&](int thid) { return thalfs[thid].x == 0; };
        while (true) {
            const auto& th_curr = thalfs[seq.back()];
            const auto& th_twin = thalfs[th_curr.twid];
            const auto& tq_twin = tquads[th_twin.tqid];
            auto s0 = tq_twin.side_of(th_twin);
            auto h1 = tq_twin.thids((s0 + 1) % 4);
            auto h2 = tq_twin.thids((s0 + 2) % 4);
            auto h3 = tq_twin.thids((s0 + 3) % 4);
            double s = 0;
            for (int t : h2) s += thalfs[t].x;
            if (s > 0) break;
            if (h2.size() > 1 || rg::any_of(h1, zero) || rg::any_of(h3, zero)) return false;
            seq.push_back(h2.front());
        }
        return true;
    };

    vec seq_l = {tq.thids(side == 0 ? 2 : 3)[0]};
    vec seq_r = {tq.thids(side == 0 ? 0 : 1)[0]};
    if (!find_simple_chain(seq_l)) { return false; }
    if (!find_simple_chain(seq_r)) { return false; }

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

    auto find_terminal = [&](int thid, int s0, int s1) -> std::pair<const HmLoc&, int> {
        const auto& th0 = thalfs[thid];
        const auto& th1 = thalfs[th0.twid];
        if (th0.bgn) return { th0.loc_fr(), s0 };
        if (th1.bgn) return { th1.loc_fr(), s1 };
        if (th0.end) return { th0.loc_to(), s1 };
        if (th1.end) return { th1.loc_to(), s0 };
        throw std::runtime_error("error in find_terminal (chain)");
    };
    auto [loc_bgn, side_bgn] = find_terminal(chain.thid_l, 0, 1);
    auto [loc_end, side_end] = find_terminal(chain.thid_r, 1, 0);

    int s = 0; for (int t : chain.thids_t) s += thalfs[t].x;
    chain.pts.push_back({ .loc = loc_bgn, .val = 0, .adj = count_adj_tquads(chain.thid_l), .side = side_bgn });
    chain.pts.push_back({ .loc = loc_end, .val = s, .adj = count_adj_tquads(chain.thid_r), .side = side_end });

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

        // push inter tquad points
        auto btm = oft;
        auto f = [&](const auto& l) { return l != th_l.loc_fr() && l != th_l.loc_to() && l != th_r.loc_fr() && l != th_r.loc_to(); };
        for (int thid: thids_t) { auto& th = thalfs[thid]; auto& l = th.loc_to(); oft += th.x; if (f(l)) chain.pts.push_back({ .loc = l, .val = oft, .adj = 3, .side = 0}); }
        for (int thid: thids_b) { auto& th = thalfs[thid]; auto& l = th.loc_to(); btm += th.x; if (f(l)) chain.pts.push_back({ .loc = l, .val = btm, .adj = 3, .side = 1}); }
        chain.bounds.push_back(oft);

        // push ladder thalf points
        if (i == chain.thids_z.size() - 1) break;
        if      (th_r.bgn)   chain.pts.push_back({.loc = l1, .val = oft, .adj = a1, .side = 1 });
        else if (th_twn.bgn) chain.pts.push_back({.loc = l0, .val = oft, .adj = a0, .side = 0 });
        else if (th_r.end)   chain.pts.push_back({.loc = l0, .val = oft, .adj = a0, .side = 0 });
        else if (th_twn.end) chain.pts.push_back({.loc = l1, .val = oft, .adj = a1, .side = 1 });
    }

    rg::stable_sort(chain.pts, {}, &Tqpoint::val);
    return true;
}
