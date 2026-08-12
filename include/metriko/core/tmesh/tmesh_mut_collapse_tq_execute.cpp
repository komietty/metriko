#include "./tmesh_mut.h"
using namespace metriko;

void TmeshMut::collapse_tquad_chain_execute(Tqchain& chain) {

    struct Candidate {
        int  tqid;
        int  thid;
        bool t_fr;
        bool t_to;
        int  v_fr;
        int  v_to;
    };
    vec<Candidate> candidates;

    // create map of tqid and region
    std::map<int, vec<std::tuple<int, double, double>>> regions;
    for (int tqid : chain.tqids) regions[tqid] = allowed_range(tqid);

    auto tquad_of = [&](int v) {
        for (size_t k = 0; k < chain.bounds.size(); ++k)
            if (v < chain.bounds[k]) return chain.tqids[k]; // first tquad whose right edge exceeds v
        return chain.tqids.back();                          // v == s (last edge)
    };

    auto thid_of = [&](const HmLoc& fr, const HmLoc& to) -> std::optional<int> {
        for (auto* side: { &chain.thids_t, &chain.thids_b }) {
        for (auto thid: *side) {
            auto& th = thalfs[thid];
            bool  f1 = th.loc_fr() == fr && th.loc_to() == to;
            bool  f2 = th.loc_fr() == to && th.loc_to() == fr;
            if (f1 || f2) return thid;
        }}
        return std::nullopt;
    };

    auto thids_from_remainning_side = [&](const HmLoc& fr, const HmLoc& to) -> vec<int> {
        for (const HmLoc& a : { fr, to }) {
            const HmLoc& b = a == fr ? to : fr;

            for (const auto* side : { &chain.thids_t, &chain.thids_b }) {
                vec<int> res;
                HmLoc cur = a;
                while (cur != b) {
                    auto it = rg::find_if(*side, [&](int t) { return thalfs[t].loc_fr() == cur; });
                    if (it == side->end()) { break; } // cannot continue on this side
                    res.push_back(*it);
                    cur = thalfs[*it].loc_to(); // step loc_fr -> loc_to
                }
                if (cur == b) { if (a != fr) rg::reverse(res); return res; }  // normalize to fr -> to order
            }
        }
        throw std::runtime_error("could not find thids from remainning side");
    };

    for (int i = 0; i < chain.pts.size() - 1; ++i) {
        auto& [l0, v0, ord0, adj0, s0] = chain.pts[i];
        auto& [l1, v1, ord1, adj1, s1] = chain.pts[i + 1];
        auto tqid = tquad_of(std::min(v0, v1));

        if (s0 == s1) {
            candidates.push_back({ .tqid = tqid, .thid = thid_of(l0, l1).value(), .t_fr = s0, .t_to = s1, .v_fr = v0, .v_to = v1 });
        } else {
            auto id0  = rg::find(tnodes, l0) - tnodes.begin();
            auto id1  = rg::find(tnodes, l1) - tnodes.begin();
            auto path = approx_shortest_path(30, hm, l0, l1, regions[tqid]);
            auto nids = add_new_path(path, id0, id1);
            int teid  = tedges.size();
            int thid0 = thalfs.size();
            int thid1 = thalfs.size() + 1;
            double x  = std::abs(v1 - v0);
            double r = 0;   // geometric length of the traced path
            for (size_t k = 0; k + 1 < nids.size(); ++k)
                r += (get_ptloc_pos(hm, tnodes[nids[k + 1]]) - get_ptloc_pos(hm, tnodes[nids[k]])).norm();

            tedges.push_back({ .id = teid, .nids = nids });
            thalfs.push_back({ .tm = this, .id = thid0, .twid = thid1, .teid = teid, .cano = true,  .x = x, .r = r });
            thalfs.push_back({ .tm = this, .id = thid1, .twid = thid0, .teid = teid, .cano = false, .x = x, .r = r });
            candidates.push_back({ .tqid = tqid, .thid = thid0, .t_fr = s0, .t_to = s1, .v_fr = v0, .v_to = v1});
            candidates.push_back({ .tqid = tqid, .thid = thid1, .t_fr = s1, .t_to = s0, .v_fr = v1, .v_to = v0});
        }
    }

    auto consume_pool = [&](const HmLoc& loc, bool side_to_stop, bool invert) {
        vec<int> res;
        HmLoc cur = loc;
        while (true) {
            auto it = rg::find_if(candidates, [&](auto& t) {
                const auto& th = thalfs[t.thid];
                return invert? (th.loc_to() == cur) : (th.loc_fr() == cur);
            });
            if (it == candidates.end()) break;
            auto thid = it->thid;
            auto flag = invert ? (it->t_fr == side_to_stop) : (it->t_to == side_to_stop);
            res.push_back(thid);
            cur = invert ? thalfs[thid].loc_fr() : thalfs[thid].loc_to();
            candidates.erase(it);
            if (flag) break;
        }
        if (invert) rg::reverse(res);
        return res;
    };

    auto& [loc_bgn, v_bgn, o_bgn, a_bgn, is_top_bgn] = chain.pts.front(); // left
    auto& [loc_end, v_end, o_end, a_end, is_top_end] = chain.pts.back();  // right
    //std::println("is top bgn: {}", is_top_bgn);
    //std::println("is top end: {}", is_top_end);
    vec<int> thids_bgn = consume_pool(loc_bgn, !is_top_bgn,  is_top_bgn);
    vec<int> thids_end = consume_pool(loc_end, !is_top_end, !is_top_end);
    vec<vec<int>> chains;

    //std::println("l_bgn: {}, s_bgn: {}, a_bgn: {}, invert: {}, thdis_bgn: {}", loc_str(l_bgn), s_bgn, a_bgn, s_bgn == 1, thids_bgn);
    //std::println("l_end: {}, s_end: {}, a_end: {}, invert: {}, thdis_end: {}", loc_str(l_end), s_end, a_end, s_end == 0, thids_end);

    auto is_head = [&](const Candidate& q) { return rg::none_of(candidates, [&](const Candidate& r) { return thalfs[r.thid].loc_to() == thalfs[q.thid].loc_fr(); }); };
    auto is_tail = [&](const Candidate& q) { return rg::none_of(candidates, [&](const Candidate& r) { return thalfs[r.thid].loc_fr() == thalfs[q.thid].loc_to(); }); };

    while (!candidates.empty()) {
        if (auto it = rg::find_if(candidates, is_head); it != candidates.end()) { chains.push_back(consume_pool(thalfs[it->thid].loc_fr(), it->t_fr, false)); continue; }
        if (auto it = rg::find_if(candidates, is_tail); it != candidates.end()) { chains.push_back(consume_pool(thalfs[it->thid].loc_to(), it->t_to, true));  continue; }
        throw std::runtime_error("error in collapse_tquad_execute");
    }

    //std::cout << "thids_bgn size: " << thids_bgn.size() << std::endl;
    //std::cout << "thids_end size: " << thids_end.size() << std::endl;
    //for (auto& c: chains) { std::cout << "chain size: " << c.size() << std::endl; }

    auto replace = [&](const vec<int>& thids_replace, const vec<int>& chain) {
        size_t ci = 0;
        for (int old : thids_replace) {
            const HmLoc& e0 = thalfs[old].loc_fr();
            const HmLoc& e1 = thalfs[old].loc_to();
            vec<int> seg;
            while (ci < chain.size()) {
                seg.push_back(chain[ci]);
                auto& to = thalfs[chain[ci]].loc_to(); ++ci;
                if (ci >= chain.size() || to == e0 || to == e1) break;
            }

            // same body as single-tquad replace, per element
            auto& [id, data] = tquads[thalfs[old].tqid];
            auto it = rg::find(data, old, &TdataMut::thid); assert(it != data.end());
            auto si = it->side;
            it = data.erase(it);
            vec<TdataMut> repl;
            for (int t : seg) { thalfs[t].tqid = id; repl.push_back({ t, si }); }
            data.insert(it, repl.begin(), repl.end());
        }
    };

    auto extend = [&](const ThalfMut& th, bool ahd, bool not_consume) {
        auto step = [&](int thid) { return ahd ? step_next(thid) : step_prev(thid); };
        auto& th1 = thalfs[step(th.id)];    // nxt or prv
        auto& th2 = thalfs[step(th1.twid)]; // ahd or bhd
        auto& [id, data] = tquads[th2.tqid];
        if (not_consume) {
            auto it = rg::find(data, th2.id, &TdataMut::thid);
            auto si = it->side;
            thalfs[th.id].tqid = id;
            data.insert(ahd ? it : it + 1, TdataMut{ th.id, si });
        } else {
            auto& nids    = tedges[th.teid].nids;
            auto& th_twn  = thalfs[th.twid];
            auto& tq_twn  = tquads[th_twn.tqid];
            auto& th2_twn = thalfs[th2.twid];
            std::erase_if(tq_twn.data, [&](const auto& d) { return d.thid == th_twn.id; });
            tedges[th2.teid].insert_locs(nids);

            if      (th.bgn)     { if (th2.loc_fr() == th.loc_fr())     th2.bgn = true; else th2_twn.bgn = true; }
            else if (th_twn.bgn) { if (th2.loc_fr() == th_twn.loc_fr()) th2.bgn = true; else th2_twn.bgn = true; }
            else if (th.end)     { if (th2.loc_to() == th.loc_to())     th2.end = true; else th2_twn.end = true; }
            else if (th_twn.end) { if (th2.loc_to() == th_twn.loc_to()) th2.end = true; else th2_twn.end = true; }
        }
    };

    auto on_chain = [&](const HmLoc& l) { return rg::any_of(chain.pts, [&](const auto& c) { return c.loc == l; }); };
    const auto& th_l = thalfs[chain.thid_l];
    const auto& th_r = thalfs[chain.thid_r];
    bool cfr_r = count_adj_tquads(th_r.id) == 4   && !on_chain(th_r.loc_fr());
    bool cto_r = count_adj_tquads(th_r.twid) == 4 && !on_chain(th_r.loc_to());
    bool cfr_l = count_adj_tquads(th_l.id) == 4   && !on_chain(th_l.loc_fr());
    bool cto_l = count_adj_tquads(th_l.twid) == 4 && !on_chain(th_l.loc_to());

    if (!thids_bgn.empty() && thids_end.empty() && chain.pts.size() == 2) {
        assert(is_top_bgn == is_top_end);
        assert(thids_bgn.size() == 1);
        bool ahd_l = th_l.loc_fr() == loc_bgn;
        bool ahd_r = th_r.loc_fr() == loc_end;
        extend(th_l, ahd_l, cfr_l || cto_l);
        extend(th_r, ahd_r, cfr_r || cto_r);
        auto l0 = ahd_l ? th_l.loc_to() : th_l.loc_fr();
        auto l1 = ahd_r ? th_r.loc_to() : th_r.loc_fr();
        auto op = thids_from_remainning_side(l0, l1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>();
        replace(op, thids_bgn);
    } else {
        { // leftmost
            bool ahd = th_l.loc_fr() == loc_bgn;
            extend(th_l, ahd, cfr_l || cto_l);
            auto l0 = ahd ? th_l.loc_to(): th_l.loc_fr();
            auto la = thalfs[thids_bgn.back()].loc_to();
            auto lb = thalfs[thids_bgn.front()].loc_fr();
            auto l1 = la == th_l.loc_fr() || la == th_l.loc_to() ? lb : la;
            auto op = thids_from_remainning_side(l0, l1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>();
            replace(op, thids_bgn);
        }
        { // rightmost
            bool ahd = th_r.loc_fr() == loc_end;
            extend(th_r, ahd, cfr_r || cto_r);
            auto l0 = ahd ? th_r.loc_to() : th_r.loc_fr();
            auto la = thalfs[thids_end.back()].loc_to();
            auto lb = thalfs[thids_end.front()].loc_fr();
            auto l1 = la == th_r.loc_fr() || la == th_r.loc_to() ? lb : la;
            auto op = thids_from_remainning_side(l0, l1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>();
            replace(op, thids_end);
        }
        // in middle
        for (auto& c: chains) {
            assert(c.size() >= 2);
            auto l0 = thalfs[c.front()].loc_fr();
            auto l1 = thalfs[c.back()].loc_to();
            auto op = thids_from_remainning_side(l0, l1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>();
            replace(op, c);
        }
    }

    // clean up
    for (int tqid: chain.tqids) {
        auto& [id, data] = tquads[tqid];
        for (const auto& [thid, _]: data) {
            auto& th0 = thalfs[thid];
            auto& th1 = thalfs[th0.twid];
            auto& te  = tedges[th0.teid];
            if (th0.tqid == id) { th0 = {}; th1 = {}; te = {}; }
        }
        data.clear();
        id = -1;
    }
}