#include "./tmesh_mut.h"
using namespace metriko;

void TmeshMut::collapse_tquad_chain_execute(Tqchain& chain) {

    struct Candidate {
        int tqid;
        int thid;
        int s_fr;
        int s_to;
        int v_fr;
        int v_to;
    };
    vec<Candidate> candidates;

    // create map of tqid and region
    std::map<int, vec<std::tuple<int, double, double>>> regions;
    for (int tqid : chain.tqids) regions[tqid] = allowed_range(tqid);

    auto tquad_of = [&](int v) {
        for (size_t k = 0; k < chain.bounds.size(); ++k)
            if (v < chain.bounds[k]) return chain.tqids[k];   // first tquad whose right edge exceeds v
        return chain.tqids.back();                            // v == s (last edge)
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

    auto thids_of = [&](const HmLoc& fr, const HmLoc& to) -> std::optional<vec<int>> {
        for (const HmLoc& a : { fr, to }) {
            const HmLoc& b = a == fr ? to : fr;

            for (const auto* side : { &chain.thids_t, &chain.thids_b }) {
                vec<int> res;
                HmLoc cur = a;
                while (cur != b) {
                    //std::println("iter");
                    auto it = rg::find_if(*side, [&](int t) { return thalfs[t].loc_fr() == cur; });
                    if (it == side->end()) {
                        std::println("side_end");
                        break;
                    }        // cannot continue on this side
                    res.push_back(*it);
                    cur = thalfs[*it].loc_to();          // step loc_fr -> loc_to
                }
                if (cur == b) { if (a != fr) rg::reverse(res); return res; }  // normalize to fr -> to order
            }
        }
        return std::nullopt;
    };

    for (int i = 0; i < chain.pts.size() - 1; ++i) {
        auto& [l0, v0, adj0, s0] = chain.pts[i];
        auto& [l1, v1, adj1, s1] = chain.pts[i + 1];
        auto tqid = tquad_of(std::min(v0, v1));

        if (s0 == s1) {
            candidates.push_back({ .tqid = tqid, .thid = thid_of(l0, l1).value(), .s_fr = s0, .s_to = s1, .v_fr = v0, .v_to = v1 });
        } else {
            auto id0  = rg::find(tnodes, l0) - tnodes.begin();
            auto id1  = rg::find(tnodes, l1) - tnodes.begin();
            auto path = approx_shortest_path(20, hm, l0, l1, regions[tqid]);
            auto nids = add_new_path(path, id0, id1);
            int teid  = tedges.size();
            int thid0 = thalfs.size();
            int thid1 = thalfs.size() + 1;
            double x  = std::abs(v1 - v0);
            tedges.push_back({ .nids = nids });
            thalfs.push_back({ .tm = this, .id = thid0, .twid = thid1, .teid = teid, .cano = true,  .x = x });
            thalfs.push_back({ .tm = this, .id = thid1, .twid = thid0, .teid = teid, .cano = false, .x = x });
            candidates.push_back({ .tqid = tqid, .thid = thid0, .s_fr = s0, .s_to = s1, .v_fr = v0, .v_to = v1});
            candidates.push_back({ .tqid = tqid, .thid = thid1, .s_fr = s1, .s_to = s0, .v_fr = v1, .v_to = v0});
        }
    }

    auto& [l_bgn, v_bgn, a_bgn, s_bgn] = chain.pts.front(); // left
    auto& [l_end, v_end, a_end, s_end] = chain.pts.back();  // right

    auto consume_pool = [&](const HmLoc& loc, int side_to_stop, bool invert) {
        vec<int> res;
        HmLoc cur = loc;
        while (true) {
            auto it = rg::find_if(candidates, [&](auto& t) {
                const auto& th = thalfs[t.thid];
                return invert? (th.loc_to() == cur) : (th.loc_fr() == cur);
            });
            if (it == candidates.end()) break;
            auto thid = it->thid;
            auto flag = invert ? (it->s_fr == side_to_stop) : (it->s_to == side_to_stop);
            res.push_back(thid);
            cur = invert ? thalfs[thid].loc_fr() : thalfs[thid].loc_to();
            candidates.erase(it);
            if (flag) break;
        }
        if (invert) rg::reverse(res);
        return res;
    };

    vec<int> thids_bgn = consume_pool(l_bgn, s_end, s_bgn == 0);
    vec<int> thids_end = consume_pool(l_end, s_bgn, s_end == 1);
    vec<vec<int>> chains;

    std::println("l_bgn: {}, s_bgn: {}, a_bgn: {}, invert: {}, thdis_bgn: {}", loc_str(l_bgn), s_bgn, a_bgn, s_bgn == 1, thids_bgn);
    std::println("l_end: {}, s_end: {}, a_end: {}, invert: {}, thdis_end: {}", loc_str(l_end), s_end, a_end, s_end == 0, thids_end);


    while (!candidates.empty()) {
        int lo = candidates.front().v_fr;
        int hi = lo;
        for (auto& q : candidates) {
            lo = std::min({ lo, q.v_fr, q.v_to });
            hi = std::max({ hi, q.v_fr, q.v_to });
        }
        auto ext = [&](double v) { return v == lo || v == hi; };
        if (auto it = rg::find_if(candidates, [&](auto& q) { return ext(q.v_fr); }); it != candidates.end()) { chains.push_back(consume_pool(thalfs[it->thid].loc_fr(), it->s_fr, false)); continue; }
        if (auto it = rg::find_if(candidates, [&](auto& q) { return ext(q.v_to); }); it != candidates.end()) { chains.push_back(consume_pool(thalfs[it->thid].loc_to(), it->s_to, true));  continue; }
        throw std::runtime_error("error in collapse_tquad_execute");
    }

    std::cout << "thids_bgn size: " << thids_bgn.size() << std::endl;
    for (auto& c: chains) { std::cout << "chain size: " << c.size() << std::endl; }
    std::cout << "thids_end size: " << thids_end.size() << std::endl;

    auto replace = [&](const vec<int>& thids_replace, const vec<int>& chain) {
        size_t ci = 0;
        for (int old : thids_replace) {
            // take the chain slice covering old's span (aligned by the node shared with chain's current pos)
            const HmLoc& near = thalfs[chain[ci]].loc_fr();
            const HmLoc& b    = thalfs[old].loc_fr() == near ? thalfs[old].loc_to() : thalfs[old].loc_fr();
            vec<int> seg;
            do { seg.push_back(chain[ci]); }
            while (ci + 1 < chain.size() && thalfs[chain[ci++]].loc_to() != b);

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

    auto extend = [&](const ThalfMut& th, bool ahd, int count_fr, int count_to) {
        auto step = [&](int thid) { return ahd ? step_next(thid) : step_prev(thid); };
        auto& th1 = thalfs[step(th.id)];    // nxt or prv
        auto& th2 = thalfs[step(th1.twid)]; // ahd or bhd
        auto& [id, data] = tquads[th2.tqid];
        if (count_fr == 4 || count_to == 4) {
            auto it = rg::find(data, th2.id, &TdataMut::thid);
            auto si = it->side;
            thalfs[th.id].tqid = id;
            data.insert(ahd ? it : it + 1, TdataMut{ th.id, si });
        } else {
            auto& [nids]  = tedges[th.teid];
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

    { // leftmost
        const auto& th = thalfs[chain.thid_l];
        auto cfr = count_adj_tquads(th.id);
        auto cto = count_adj_tquads(th.twid);
        bool ahd = th.loc_fr() == l_bgn;
        extend(th, ahd, cfr, cto);

        auto l0 = ahd ? th.loc_to(): th.loc_fr();
        auto l1 = thalfs[thids_bgn.back()].loc_to();
        if (l1 == th.loc_fr() || l1 == th.loc_to()) { l1 = thalfs[thids_bgn.front()].loc_fr(); }
        auto op = thids_of(l0, l1).value()
            | vw::transform([&](int t) { return thalfs[t].twid; })
            | rg::to<vec<int>>();
        replace(op, thids_bgn);
    }
    { // rightmost
        const auto& th = thalfs[chain.thid_r];
        auto cfr = count_adj_tquads(th.id);
        auto cto = count_adj_tquads(th.twid);
        bool ahd = th.loc_fr() == l_end;
        extend(th, ahd, cfr, cto);

        auto l0 = ahd ? th.loc_to(): th.loc_fr();
        auto l1 = thalfs[thids_end.back()].loc_to();
        if (l1 == th.loc_fr() || l1 == th.loc_to()) { l1 = thalfs[thids_end.front()].loc_fr(); }
        auto op = thids_of(l0, l1).value()
            | vw::transform([&](int t) { return thalfs[t].twid; })
            | rg::to<vec<int>>();
        replace(op, thids_end);
    }

    for (auto& c: chains) {
        assert(c.size() >= 2);
        auto o = thids_of(thalfs[c.front()].loc_fr(), thalfs[c.back()].loc_to());
        auto o2 = o.value()                                        // twins (surviving side), kept in fr->to order
            | vw::transform([&](int t) { return thalfs[t].twid; })
            | rg::to<vec<int>>();
        std::println("chain");
        replace(o2, c);
    }

    for (int tqid: chain.tqids) {
        auto& [id, data] = tquads[tqid];
        for (const auto& [thid, _]: data) {
            auto& th0 = thalfs[thid];
            auto& th1 = thalfs[th0.twid];
            if (th0.tqid == id) { th0 = {}; th1 = {}; }
        }
        data.clear();
        id = -1;
    }
}
