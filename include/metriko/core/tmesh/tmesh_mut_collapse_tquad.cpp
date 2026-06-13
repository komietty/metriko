#include "./tmesh_mut.h"
#include <cassert>
using namespace metriko;

bool TmeshMut::collapse_tquad_prepare(int tqid, Tqaux& tqaux) const {
    auto& tq = tquads[tqid];
    int side = -1;
    if (tq.thids(0).size() == 1 && tq.thids(2).size() == 1 && thalfs[tq.thids(0).front()].x == 0) side = 0;
    if (tq.thids(1).size() == 1 && tq.thids(3).size() == 1 && thalfs[tq.thids(1).front()].x == 0) side = 1;
    if (side == -1) return false;

    auto side_r  = side == 0 ? 0 : 1; // collapse
    auto side_t  = side == 0 ? 1 : 2; // remain
    auto side_l  = side == 0 ? 2 : 3; // collapse
    auto side_b  = side == 0 ? 3 : 0; // remain
    auto thid_l  = tq.thids(side_l)[0];
    auto thid_r  = tq.thids(side_r)[0];
    auto thids_t = tq.thids(side_t);
    auto thids_b = tq.thids(side_b);
    auto& th_l   = thalfs[thid_l];
    auto& th_r   = thalfs[thid_r];

    auto find_terminal = [&](int thid, int s0, int s1) -> std::pair<const HmLoc&, int> {
        auto& th0 = thalfs[thid];
        auto& th1 = thalfs[th0.twid];
        if (th0.bgn) return { th0.loc_fr(), s0 };
        if (th1.bgn) return { th1.loc_fr(), s1 };
        if (th0.end) return { th0.loc_to(), s1 };
        if (th1.end) return { th1.loc_to(), s0 };
        throw std::runtime_error("error in find_terminal");
    };

    auto [loc_bgn, side_bgn] = find_terminal(thid_r, side_b, side_t);
    auto [loc_end, side_end] = find_terminal(thid_l, side_t, side_b);
    auto s = rg::fold_left(thids_t | vw::transform([&](int t){ return thalfs[t].x; }), 0., std::plus{});

    vec<std::tuple<HmLoc, double, int>> aux;
    aux.emplace_back(loc_bgn, 0, side_bgn);
    aux.emplace_back(loc_end, s, side_end);

    auto v = 0.;
    auto f = [&](const HmLoc& l) {
        return l != th_l.loc_fr() &&
               l != th_l.loc_to() &&
               l != th_r.loc_fr() &&
               l != th_r.loc_to();
    };
    for (int i: thids_t) { auto& th = thalfs[i]; auto& l = th.loc_to(); v += th.x; if (f(l)) aux.emplace_back(l, v, side_t); }
    for (int i: thids_b) { auto& th = thalfs[i]; auto& l = th.loc_to(); v -= th.x; if (f(l)) aux.emplace_back(l, v, side_b); }

    // sort by val. on ties (exactly one per side), put the one matching the previous element's side first
    rg::stable_sort(aux, {}, [](const auto& t) { return std::get<1>(t); });
    for (size_t i = 1; i + 1 < aux.size(); ++i) {
        auto v1 = std::get<1>(aux[i]);
        auto v2 = std::get<1>(aux[i + 1]);
        auto s0 = std::get<2>(aux[i - 1]);
        auto s1 = std::get<2>(aux[i]);
        auto s2 = std::get<2>(aux[i + 1]);
        if (v1 == v2) { assert(s1 != s2); if(s1 != s0) std::swap(aux[i], aux[i + 1]); }
    }

    tqaux.checkpoints  = aux;
    tqaux.side_thid_l  = std::pair(side_l, thid_l);
    tqaux.side_thid_r  = std::pair(side_r, thid_r);
    tqaux.side_thids_t = std::pair(side_t, thids_t);
    tqaux.side_thids_b = std::pair(side_b, thids_b);
    tqaux.thid_r_merge_to_ahead = side_bgn == side_b;
    tqaux.thid_l_merge_to_ahead = side_end == side_t;
    return true;
}

void TmeshMut::collapse_tquad_execute(int tqid, Tqaux& tqaux) {
    auto& tq_crr = tquads[tqid];
    auto  region = allowed_range(tq_crr.id);

    struct Q {
        int thid;
        int side_fr;
        int side_to;
        double val_fr;
        double val_to;
    };

    vec<Q> qs;

    auto find_thid_in_tquad = [&](const HmLoc& fr, const HmLoc& to) -> std::optional<int> {
        for (auto& [thid, _]: tq_crr.data) {
            auto& th = thalfs[thid];
            bool  f1 = th.loc_fr() == fr && th.loc_to() == to;
            bool  f2 = th.loc_fr() == to && th.loc_to() == fr;
            if (f1 || f2) return thid; // use both dir in this case...
        }
        return std::nullopt;
    };

    for (int i = 0; i < tqaux.checkpoints.size() - 1; ++i) {
        auto& [loc0, val0, side0] = tqaux.checkpoints[i];
        auto& [loc1, val1, side1] = tqaux.checkpoints[i + 1];

        if (side0 == side1) {
            int thid = find_thid_in_tquad(loc0, loc1).value();
            qs.emplace_back(thid, side0, side1, val0, val1);
        }
        else {
            auto path = approx_shortest_path(20, hm, loc0, loc1, region);

            vec<int> nids;
            int fr_idx = rg::find(tnodes, loc0) - tnodes.begin();
            int to_idx = rg::find(tnodes, loc1) - tnodes.begin();
            for (int j = 0; j < path.size(); ++j) {
                if      (j == 0)                  nids.push_back(fr_idx);
                else if (j + 1 == path.size())    nids.push_back(to_idx);
                else { tnodes.push_back(path[j]); nids.push_back(tnodes.size() - 1); }
            }

            int teid  = tedges.size();
            int thid0 = thalfs.size();
            int thid1 = thalfs.size() + 1;
            double x  = std::abs(val1 - val0);
            tedges.push_back({ .nids = nids });
            thalfs.push_back({ .tm = this, .id = thid0, .twid = thid1, .teid = teid, .cano = true,  .x = x });
            thalfs.push_back({ .tm = this, .id = thid1, .twid = thid0, .teid = teid, .cano = false, .x = x });
            qs.emplace_back(thid0, side0, side1, val0, val1);
            qs.emplace_back(thid1, side1, side0, val1, val0);
        }
    }

    auto& [loc_bgn, val_bgn, side_bgn] = tqaux.checkpoints.front();
    auto& [loc_end, val_end, side_end] = tqaux.checkpoints.back();

    auto consume_pool = [&](const HmLoc& loc, int side_to_stop, bool invert) {
        vec<int> res;
        HmLoc cur = loc;
        while (true) {
            auto it = rg::find_if(qs, [&](auto& t) {
                const auto& th = thalfs[t.thid];
                if (invert) return th.loc_to() == cur;
                else        return th.loc_fr() == cur;
            });
            if (it == qs.end()) break;
            auto thid = it->thid;
            auto flag = it->side_to == side_to_stop;
            if (invert) flag = it->side_fr == side_to_stop;
            res.push_back(thid);
            if (invert) cur = thalfs[thid].loc_fr();
            else        cur = thalfs[thid].loc_to();
            qs.erase(it);
            if (flag) break;
        }
        if (invert) rg::reverse(res);
        return res;
    };

    vec<int> thids_bgn = consume_pool(loc_bgn, side_end, false);
    std::cout << "thids_bgn size (fwd): " << thids_bgn.size() << std::endl;
    if (thids_bgn.empty()) thids_bgn = consume_pool(loc_bgn, side_end, true);
    std::cout << "thids_bgn size (bwd): " << thids_bgn.size() << std::endl;

    vec<int> thids_end = consume_pool(loc_end, side_bgn, false);
    std::cout << "thids_end size (fwd): " << thids_end.size() << std::endl;
    if (thids_end.empty()) thids_end = consume_pool(loc_end, side_bgn, true);
    std::cout << "thids_end size (bwd): " << thids_end.size() << std::endl;

    vec<vec<int>> chains;
    while (!qs.empty()) {
        // from all remaining qs, get min and max val（qs 全体に対して1回）
        double min_val = qs.front().val_fr;
        double max_val = qs.front().val_fr;
        for (auto& q : qs) {
            min_val = std::min({ min_val, q.val_fr, q.val_to });
            max_val = std::max({ max_val, q.val_fr, q.val_to });
        }

        // 開始 q を1つ選ぶ（consume_pool が qs を erase するので、for で回しながら消さない）
        auto it = rg::find_if(qs, [&](auto& q) { return q.val_fr == min_val || q.val_fr == max_val; });
        if (it != qs.end()) {
            int thid = it->thid, side_fr = it->side_fr;
            auto res = consume_pool(thalfs[thid].loc_fr(), side_fr, false);
            assert(res.size() >= 2);
            chains.push_back(res);
        } else {
            auto it = rg::find_if(qs, [&](auto& q) { return q.val_to == min_val || q.val_to == max_val; });
            int thid = it->thid, side_to = it->side_to;
            auto res = consume_pool(thalfs[thid].loc_to(), side_to, true);
            assert(res.size() >= 2);
            chains.push_back(res);
        }
    }

    for (auto& chain: chains) {
        std::cout << "chain size: " << chain.size() << std::endl;
        for (int thid: chain) { std::cout << "thid: " << thid << std::endl; }
    }
    for (int thid: thids_bgn) { std::cout << "thid_bgn: " << thid << std::endl; }
    for (int thid: thids_end) { std::cout << "thid_end: " << thid << std::endl; }


    auto step_next = [&](int thid) { auto& q = tquads[thalfs[thid].tqid]; return circular_next(q.data, rg::find(q.data, thid, &TdataMut::thid))->thid; };
    auto step_prev = [&](int thid) { auto& q = tquads[thalfs[thid].tqid]; return circular_prev(q.data, rg::find(q.data, thid, &TdataMut::thid))->thid; };

    auto replace = [&](int thid_replace, const vec<int>& chain) {
        auto& th_replace = thalfs[thid_replace];
        auto& tq_replace = tquads[th_replace.tqid];
        int   side = tq_replace.side_of(th_replace);
        auto& data = tq_replace.data;
        int    pos = rg::find(data, side, &TdataMut::side) - data.begin();
        std::erase_if(data, [&](const TdataMut& d) { return d.side == side; });
        vec<TdataMut> repl;
        for (int thid : chain) { thalfs[thid].tqid = tq_replace.id; repl.push_back({ thid, side }); }
        data.insert(data.begin() + pos, repl.begin(), repl.end());
    };

    auto extend_and_replace = [&](const ThalfMut& th, bool ahd, const vec<int>& chain) {
        auto& th_nxt = thalfs[step_next(th.id)];
        auto& th_prv = thalfs[step_prev(th.id)];
        auto& te_ahd = tedges[thalfs[step_next(th_nxt.twid)].teid];
        auto& te_bhd = tedges[thalfs[step_prev(th_prv.twid)].teid];
        auto& [_, nids]  = tedges[th.teid];
        vec<int> wo_l(nids.begin(), nids.end() - 1);
        vec<int> wo_f(nids.begin() + 1, nids.end());
        if (ahd) { if (th.cano) te_ahd.insert_locs_front(wo_l); else te_ahd.insert_locs_after(wo_f); replace(th_nxt.twid, chain); }
        else     { if (th.cano) te_bhd.insert_locs_after(wo_f); else te_bhd.insert_locs_front(wo_l); replace(th_prv.twid, chain); }
    };

    const auto& th_r = thalfs[tqaux.side_thid_r.second];
    const auto& th_l = thalfs[tqaux.side_thid_l.second];
    extend_and_replace(th_r, tqaux.thid_r_merge_to_ahead, thids_bgn);
    extend_and_replace(th_l, tqaux.thid_l_merge_to_ahead, thids_end);

    for (auto& chain: chains) {
        const HmLoc& fr = thalfs[chain.front()].loc_fr();
        const HmLoc& to = thalfs[chain.back()].loc_to();
        auto o0 = find_thid_in_tquad(fr, to);
        auto o1 = find_thid_in_tquad(to, fr);
        if (o0.has_value()) replace(thalfs[o0.value()].twid, chain);
        if (o1.has_value()) replace(thalfs[o1.value()].twid, chain);
    }

    tedges[th_r.teid].id = -1;
    tedges[th_l.teid].id = -1;
    tq_crr.data.clear();
    tq_crr.id = -1;
}

