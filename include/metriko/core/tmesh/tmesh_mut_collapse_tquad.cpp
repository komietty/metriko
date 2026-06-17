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

    if (rg::any_of(thids_t, [&](int t) { return thalfs[t].x == 0; })) return false;
    if (rg::any_of(thids_b, [&](int t) { return thalfs[t].x == 0; })) return false;

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
        for (const auto& [thid, _]: tq_crr.data) {
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

        if (side0 == side1) { qs.emplace_back(find_thid_in_tquad(loc0, loc1).value(), side0, side1, val0, val1); }
        else {
            auto id0  = rg::find(tnodes, loc0) - tnodes.begin();
            auto id1  = rg::find(tnodes, loc1) - tnodes.begin();
            auto path = approx_shortest_path(20, hm, loc0, loc1, region);
            auto nids = add_new_path(path, id0, id1);
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
                return invert? (th.loc_to() == cur) : (th.loc_fr() == cur);
            });
            if (it == qs.end()) break;
            auto thid = it->thid;
            auto flag = invert ? (it->side_fr == side_to_stop) : (it->side_to == side_to_stop);
            res.push_back(thid);
            cur = invert ? thalfs[thid].loc_fr() : thalfs[thid].loc_to();
            qs.erase(it);
            if (flag) break;
        }
        if (invert) rg::reverse(res);
        return res;
    };

    vec<int> thids_bgn = consume_pool(loc_bgn, side_end, false); if (thids_bgn.empty()) thids_bgn = consume_pool(loc_bgn, side_end, true);
    vec<int> thids_end = consume_pool(loc_end, side_bgn, false); if (thids_end.empty()) thids_end = consume_pool(loc_end, side_bgn, true);

    vec<vec<int>> chains;
    while (!qs.empty()) {
        double lo = qs.front().val_fr;
        double hi = lo;
        for (auto& q : qs) {
            lo = std::min({ lo, q.val_fr, q.val_to });
            hi = std::max({ hi, q.val_fr, q.val_to });
        }
        auto ext = [&](double v) { return v == lo || v == hi; };
        if (auto it = rg::find_if(qs, [&](auto& q) { return ext(q.val_fr); }); it != qs.end()) { chains.push_back(consume_pool(thalfs[it->thid].loc_fr(), it->side_fr, false)); continue; }
        if (auto it = rg::find_if(qs, [&](auto& q) { return ext(q.val_to); }); it != qs.end()) { chains.push_back(consume_pool(thalfs[it->thid].loc_to(), it->side_to, true));  continue; }
        throw std::runtime_error("error in collapse_tquad_execute");
    }

    auto replace = [&](int thid_replace, const vec<int>& chain) {
        auto& [id, data] = tquads[thalfs[thid_replace].tqid];
        auto it = rg::find(data, thid_replace, &TdataMut::thid); assert(it != data.end());
        auto si = it->side;
        it = data.erase(it);
        vec<TdataMut> repl;
        for (int thid : chain) { thalfs[thid].tqid = id; repl.push_back({ thid, si }); }
        data.insert(it, repl.begin(), repl.end());
    };

    auto extend_and_replace = [&](const ThalfMut& th, bool ahd, const vec<int>& chain, int count_fr, int count_to) {
        auto& th_twn = thalfs[th.twid];
        auto& tq_twn = tquads[th_twn.tqid];
        auto& [_, nids] = tedges[th.teid];

        auto merge_into = [&](TedgeMut& te) {
            auto& tn = te.nids;
            int f = nids.front();
            int b = nids.back();
            if      (tn.front() == b) { tn.insert(tn.begin(), nids.begin(), nids.end() - 1); } // prepend [f..b-1]
            else if (tn.back()  == f) { tn.insert(tn.end(),   nids.begin() + 1, nids.end()); } // append  [f+1..b]
            else if (tn.front() == f) { vec<int> s(nids.begin() + 1, nids.end()); rg::reverse(s); tn.insert(tn.begin(), s.begin(), s.end()); } // prepend reverse([f+1..b])
            else if (tn.back()  == b) { vec<int> s(nids.begin(), nids.end() - 1); rg::reverse(s); tn.insert(tn.end(),   s.begin(), s.end()); } // append  reverse([f..b-1])
            else throw std::runtime_error("merge_into: no shared corner");
        };

        if (ahd) {
            auto& th_nxt = thalfs[step_next(th.id)];
            auto& th_ahd = thalfs[step_next(th_nxt.twid)];
            auto& tq_ahd = tquads[th_ahd.tqid];
            auto& te_ahd = tedges[th_ahd.teid];
            if (count_fr == 4 || count_to == 4) {
                auto it = rg::find(tq_ahd.data, th_ahd.id, &TdataMut::thid);
                auto si = it->side;
                thalfs[th.id].tqid = tq_ahd.id;
                tq_ahd.data.insert(it, TdataMut{ th.id, si });
            } else {
                std::erase_if(tq_twn.data, [&](const auto& d) { return d.thid == th_twn.id; });
                merge_into(te_ahd);
            }
            replace(th_nxt.twid, chain);
        } else {
            auto& th_prv = thalfs[step_prev(th.id)];
            auto& th_bhd = thalfs[step_prev(th_prv.twid)];
            auto& tq_bhd = tquads[th_bhd.tqid];
            auto& te_bhd = tedges[th_bhd.teid];
            if (count_fr == 4 || count_to == 4) {
                auto it = rg::find(tq_bhd.data, th_bhd.id, &TdataMut::thid);
                auto si = it->side;
                thalfs[th.id].tqid = tq_bhd.id;
                tq_bhd.data.insert(it + 1, TdataMut{ th.id, si });
            } else {
                std::erase_if(tq_twn.data, [&](const auto& d) { return d.thid == th_twn.id; });
                merge_into(te_bhd);
            }
            replace(th_prv.twid, chain);
        }
    };

    const auto& th_r = thalfs[tqaux.side_thid_r.second];
    const auto& th_l = thalfs[tqaux.side_thid_l.second];
    int c_r_fr = count_adj_tquads(th_r.id);
    int c_l_fr = count_adj_tquads(th_l.id);
    int c_r_to = count_adj_tquads(th_r.twid);
    int c_l_to = count_adj_tquads(th_l.twid);
    extend_and_replace(th_r, tqaux.thid_r_merge_to_ahead, thids_bgn, c_r_fr, c_r_to);
    extend_and_replace(th_l, tqaux.thid_l_merge_to_ahead, thids_end, c_l_fr, c_l_to);

    for (auto& chain: chains) {
        assert(chain.size() >= 2);
        const HmLoc& fr = thalfs[chain.front()].loc_fr();
        const HmLoc& to = thalfs[chain.back()].loc_to();
        if      (auto o = find_thid_in_tquad(fr, to); o.has_value()) replace(thalfs[o.value()].twid, chain);
        else if (auto o = find_thid_in_tquad(to, fr); o.has_value()) replace(thalfs[o.value()].twid, chain);
    }

    tedges[th_r.teid].id = -1;
    tedges[th_l.teid].id = -1;

    for (const auto& [thid, _]: tq_crr.data) {
        auto& th0 = thalfs[thid];
        auto& th1 = thalfs[th0.twid];
        if (th0.tqid == tq_crr.id) { th0 = {}; th1 = {}; }
    }
    tq_crr.data.clear();
    tq_crr.id = -1;
}