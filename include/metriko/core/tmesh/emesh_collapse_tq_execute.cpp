#include "./emesh.h"
using namespace metriko;

void Emesh::collapse_tquad_chain_execute(Tqchain& chain) {

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
    for (int tqid : chain.tqids) regions[tqid] = allowed_range_tquads({tqid});

    auto tquad_of = [&](int v) {
        for (size_t k = 0; k < chain.bounds.size(); ++k)
            if (v < chain.bounds[k]) return chain.tqids[k]; // first tquad whose right edge exceeds v
        return chain.tqids.back();                          // v == s (last edge)
    };

    auto thid_of = [&](int fr, int to) -> std::optional<int> {
        for (auto* side: { &chain.thids_t, &chain.thids_b })
        for (auto thid: *side) {
            auto& th = thalfs[thid];
            if ((th.nid_fr() == fr && th.nid_to() == to) || (th.nid_fr() == to && th.nid_to() == fr)) return thid;
        }
        return std::nullopt;
    };

    auto thids_from_remainning_side = [&](int fr, int to) -> vec<int> {
        for (int a : { fr, to }) {
            int b = a == fr ? to : fr;

            for (const auto* side : { &chain.thids_t, &chain.thids_b }) {
                vec<int> res;
                int cur = a;
                while (cur != b) {
                    auto it = rg::find_if(*side, [&](int t) { return thalfs[t].nid_fr() == cur; });
                    if (it == side->end()) { break; } // cannot continue on this side
                    res.push_back(*it);
                    cur = thalfs[*it].nid_to(); // step nid_fr -> nid_to
                }
                if (cur == b) { if (a != fr) rg::reverse(res); return res; }  // normalize to fr -> to order
            }
        }
        throw std::runtime_error("could not find thids from remainning side");
    };

    for (int i = 0; i < chain.pts.size() - 1; ++i) {
        auto& [n0, v0, ord0, adj0, s0] = chain.pts[i];
        auto& [n1, v1, ord1, adj1, s1] = chain.pts[i + 1];
        auto tqid = tquad_of(std::min(v0, v1));

        if (s0 == s1) {
            candidates.push_back({ .tqid = tqid, .thid = thid_of(n0, n1).value(), .t_fr = s0, .t_to = s1, .v_fr = v0, .v_to = v1 });
        } else {
            auto path = approx_shortest_path(30, hm, tnodes[n0], tnodes[n1], regions[tqid]);

            // an empty path would silently create a tedge with empty nids, whose
            // loc_fr/loc_to dereference past a null buffer later — fail loudly here
            if (path.size() < 2) throw std::runtime_error(std::format( "[collapse tq] approx_shortest_path failed: tqid {}, {} -> {} (path size {}, allowed {})", tqid, loc_str(tnodes[n0]), loc_str(tnodes[n1]), path.size(), regions[tqid].size()));

            auto nids = add_new_path(path, n0, n1);
            int teid  = tedges.size();
            int thid0 = thalfs.size();
            int thid1 = thalfs.size() + 1;
            double x  = std::abs(v1 - v0);
            double r = path_length(nids);

            tedges.push_back({ .id = teid, .nids = nids });
            thalfs.push_back({ .tm = this, .id = thid0, .twid = thid1, .teid = teid, .cano = true,  .x = x, .r = r });
            thalfs.push_back({ .tm = this, .id = thid1, .twid = thid0, .teid = teid, .cano = false, .x = x, .r = r });
            candidates.push_back({ .tqid = tqid, .thid = thid0, .t_fr = s0, .t_to = s1, .v_fr = v0, .v_to = v1});
            candidates.push_back({ .tqid = tqid, .thid = thid1, .t_fr = s1, .t_to = s0, .v_fr = v1, .v_to = v0});
        }
    }

    auto consume_pool = [&](int nid, bool side_to_stop, bool invert) {
        vec<int> res;
        int cur = nid;
        while (true) {
            auto it = rg::find_if(candidates, [&](auto& t) {
                const auto& th = thalfs[t.thid];
                return (invert ? th.nid_to() : th.nid_fr()) == cur;
            });
            if (it == candidates.end()) break;
            auto thid = it->thid;
            auto flag = invert ? (it->t_fr == side_to_stop) : (it->t_to == side_to_stop);
            res.push_back(thid);
            cur = invert ? thalfs[thid].nid_fr() : thalfs[thid].nid_to();
            candidates.erase(it);
            if (flag) break;
        }
        if (invert) rg::reverse(res);
        return res;
    };

    auto& [nid_bgn, v_bgn, o_bgn, a_bgn, is_top_bgn] = chain.pts.front(); // left
    auto& [nid_end, v_end, o_end, a_end, is_top_end] = chain.pts.back();  // right
    vec<int> thids_bgn = consume_pool(nid_bgn, !is_top_bgn,  is_top_bgn);
    vec<int> thids_end = consume_pool(nid_end, !is_top_end, !is_top_end);
    vec<vec<int>> chains;

    auto is_head = [&](const Candidate& q) { return rg::none_of(candidates, [&](const Candidate& r) { return thalfs[r.thid].nid_to() == thalfs[q.thid].nid_fr(); }); };
    auto is_tail = [&](const Candidate& q) { return rg::none_of(candidates, [&](const Candidate& r) { return thalfs[r.thid].nid_fr() == thalfs[q.thid].nid_to(); }); };

    while (!candidates.empty()) {
        if (auto it = rg::find_if(candidates, is_head); it != candidates.end()) { chains.push_back(consume_pool(thalfs[it->thid].nid_fr(), it->t_fr, false)); continue; }
        if (auto it = rg::find_if(candidates, is_tail); it != candidates.end()) { chains.push_back(consume_pool(thalfs[it->thid].nid_to(), it->t_to, true));  continue; }
        throw std::runtime_error("error in collapse_tquad_execute");
    }

    auto replace = [&](const vec<int>& thids_replace, const vec<int>& chain) {
        size_t ci = 0;
        for (int old : thids_replace) {
            int e0 = thalfs[old].nid_fr();
            int e1 = thalfs[old].nid_to();
            vec<int> seg;
            while (ci < chain.size()) {
                seg.push_back(chain[ci]);
                int to = thalfs[chain[ci]].nid_to(); ++ci;
                if (ci >= chain.size() || to == e0 || to == e1) break;
            }

            // same body as single-tquad replace, per element
            auto& [id, data] = tquads[thalfs[old].tqid];
            auto it = rg::find(data, old, &Edata::thid); assert(it != data.end());
            auto si = it->side;
            it = data.erase(it);
            vec<Edata> repl;
            for (int t : seg) { thalfs[t].tqid = id; repl.push_back({ t, si }); }
            data.insert(it, repl.begin(), repl.end());
        }
    };

    auto extend = [&](const Ehalf& th, bool ahd, bool not_consume) {
        auto step = [&](int thid) { return ahd ? step_next(thid) : step_prev(thid); };
        auto& th1 = thalfs[step(th.id)];    // nxt or prv
        auto& th2 = thalfs[step(th1.twid)]; // ahd or bhd
        auto& [id, data] = tquads[th2.tqid];
        if (not_consume) {
            auto it = rg::find(data, th2.id, &Edata::thid);
            auto si = it->side;
            thalfs[th.id].tqid = id;
            data.insert(ahd ? it : it + 1, Edata{ th.id, si });
        } else {
            auto& nids    = tedges[th.teid].nids;
            auto& th_twn  = thalfs[th.twid];
            auto& tq_twn  = tquads[th_twn.tqid];
            auto& th2_twn = thalfs[th2.twid];
            std::erase_if(tq_twn.data, [&](const auto& d) { return d.thid == th_twn.id; });
            tedges[th2.teid].insert_locs(nids);

            if      (th.bgn)     { (th2.nid_fr() == th.nid_fr()     ? th2 : th2_twn).bgn = true; }
            else if (th_twn.bgn) { (th2.nid_fr() == th_twn.nid_fr() ? th2 : th2_twn).bgn = true; }
            else if (th.end)     { (th2.nid_to() == th.nid_to()     ? th2 : th2_twn).end = true; }
            else if (th_twn.end) { (th2.nid_to() == th_twn.nid_to() ? th2 : th2_twn).end = true; }
        }
    };

    auto on_chain = [&](int nid) { return rg::any_of(chain.pts, [&](const auto& c) { return c.nid == nid; }); };
    const auto& th_l = thalfs[chain.thid_l];
    const auto& th_r = thalfs[chain.thid_r];
    bool cfr_r = count_adj_tquads(th_r.id) == 4   && !on_chain(th_r.nid_fr());
    bool cto_r = count_adj_tquads(th_r.twid) == 4 && !on_chain(th_r.nid_to());
    bool cfr_l = count_adj_tquads(th_l.id) == 4   && !on_chain(th_l.nid_fr());
    bool cto_l = count_adj_tquads(th_l.twid) == 4 && !on_chain(th_l.nid_to());

    if (!thids_bgn.empty() && thids_end.empty()) {
        assert(is_top_bgn == is_top_end);
        assert(thids_bgn.size() == chain.pts.size() - 1);
        bool ahd_l = th_l.nid_fr() == nid_bgn;
        bool ahd_r = th_r.nid_fr() == nid_end;
        extend(th_l, ahd_l, cfr_l || cto_l);
        extend(th_r, ahd_r, cfr_r || cto_r);
        int n0 = ahd_l ? th_l.nid_to() : th_l.nid_fr();
        int n1 = ahd_r ? th_r.nid_to() : th_r.nid_fr();
        auto op = thids_from_remainning_side(n0, n1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>();
        replace(op, thids_bgn);
    } else {
        { // leftmost
            bool ahd = th_l.nid_fr() == nid_bgn;
            extend(th_l, ahd, cfr_l || cto_l);
            auto n0 = ahd ? th_l.nid_to() : th_l.nid_fr();
            auto na = thalfs[thids_bgn.back()].nid_to();
            auto nb = thalfs[thids_bgn.front()].nid_fr();
            auto n1 = na == th_l.nid_fr() || na == th_l.nid_to() ? nb : na;
            auto op = thids_from_remainning_side(n0, n1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>();
            replace(op, thids_bgn);
        }
        { // rightmost
            bool ahd = th_r.nid_fr() == nid_end;
            extend(th_r, ahd, cfr_r || cto_r);
            auto n0 = ahd ? th_r.nid_to() : th_r.nid_fr();
            auto na = thalfs[thids_end.back()].nid_to();
            auto nb = thalfs[thids_end.front()].nid_fr();
            auto n1 = na == th_r.nid_fr() || na == th_r.nid_to() ? nb : na;
            auto op = thids_from_remainning_side(n0, n1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>();
            replace(op, thids_end);
        }
        // in middle
        for (auto& c: chains) {
            assert(c.size() >= 2);
            auto n0 = thalfs[c.front()].nid_fr();
            auto n1 = thalfs[c.back()].nid_to();
            auto op = thids_from_remainning_side(n0, n1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>();
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