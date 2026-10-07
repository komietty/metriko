//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include <algorithm>
#include <format>
#include <map>
#include <set>
#include "emesh.h"
using namespace metriko;

void Emesh::collapse_tquad_chain_execute(Tqchain& chain) {
    struct Candidate {
        int  thid;
        bool t_fr;
        bool t_to;
    };
    vec<Candidate> candidates;

    // create map of tqid and region
    std::map<int, vec<Erng>> regions;
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
        METRIKO_FAIL("no thalf chain on the remaining side from node {} to node {}", fr, to);
    };

    // the merged tedge between the two chain points of a zero-width band, taken from the band itself: along n0's
    // side to its far end, then the zero side over to n1. a band that closes into a loop (around a tube) has both
    // points at the seam, and a shortest path between them would jump the seam instead of going around. empty when
    // n0 is not at the end of its side (an inner merge) or the walk does not reach n1
    auto along_band = [&](int n0, int n1) -> vec<HmLoc> {
        const auto& tz  = thalfs[rg::contains(std::array{thalfs[chain.thid_l].nid_fr(), thalfs[chain.thid_l].nid_to()}, n1) ? chain.thid_l : chain.thid_r];
        const int   end = tz.nid_fr() == n1 ? tz.nid_to() : tz.nid_fr();
        vec<int> thids;
        try { thids = thids_from_remainning_side(n0, end); } catch (const std::runtime_error&) { return {}; }
        thids.push_back(tz.id);
        vec<int> nids = { n0 };
        for (int t: thids) {
            auto ns = tedges[thalfs[t].teid].nids;
            if (ns.back() == nids.back()) rg::reverse(ns);
            if (ns.front() != nids.back()) return {};
            nids.insert(nids.end(), ns.begin() + 1, ns.end());
        }
        if (nids.back() != n1) return {};
        vec<HmLoc> path;
        for (int n: nids) path.push_back(tnodes[n]);
        return path;
    };

    auto euler = [&](const vec<Erng>& region) {
        std::set<int> fs, vs, es;

        for (auto& [eid, r0, r1]: region) {
            fs.insert(hm.edges[eid].face0().id);
            fs.insert(hm.edges[eid].face1().id);
        }

        for (int fid: fs)
        for (Half h: hm.faces[fid].adjHalfs()) {
            vs.insert(h.tail().id);
            es.insert(h.edge().id);
        }

        return (int)vs.size() - (int)es.size() + (int)fs.size();
    };

    for (int i = 0; i < chain.pts.size() - 1; ++i) {
        auto& [n0, v0, ord0, s0] = chain.pts[i];
        auto& [n1, v1, ord1, s1] = chain.pts[i + 1];
        auto tqid = tquad_of(std::min(v0, v1));

        if (s0 == s1) {
            candidates.push_back({.thid = thid_of(n0, n1).value(), .t_fr = s0, .t_to = s1 });
        } else {
            const int chi = euler(regions[tqid]);
            auto path = chi == 0 ? along_band(n0, n1) : vec<HmLoc>{};   // a band around a tube: follow the band, never a shortest path
            if (path.empty()) path = approx_shortest_path(30, hm, tnodes[n0], tnodes[n1], regions[tqid]);

            METRIKO_CHECK(path.size() >= 2, "no path within the allowed region: tqid {}", tqid);

            auto nids = add_new_path(path, n0, n1);
            int teid  = tedges.size();
            int thid0 = thalfs.size();
            int thid1 = thalfs.size() + 1;
            double x  = std::abs(v1 - v0);
            double r  = path_length(nids);

            tedges.push_back({ .id = teid, .nids = nids });
            thalfs.push_back({ .tm = this, .id = thid0, .twid = thid1, .teid = teid, .cano = true,  .x = x, .r = r });
            thalfs.push_back({ .tm = this, .id = thid1, .twid = thid0, .teid = teid, .cano = false, .x = x, .r = r });
            candidates.push_back({.thid = thid0, .t_fr = s0, .t_to = s1});
            candidates.push_back({.thid = thid1, .t_fr = s1, .t_to = s0});
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

    auto& [nid_bgn, v_bgn, o_bgn, is_top_bgn] = chain.pts.front(); // L
    auto& [nid_end, v_end, o_end, is_top_end] = chain.pts.back();  // R
    vec<int> thids_bgn = consume_pool(nid_bgn, !is_top_bgn,  is_top_bgn);
    vec<int> thids_end = consume_pool(nid_end, !is_top_end, !is_top_end);
    vec<vec<int>> chains;

    auto is_head = [&](const Candidate& q) { return rg::none_of(candidates, [&](const Candidate& r) { return thalfs[r.thid].nid_to() == thalfs[q.thid].nid_fr(); }); };
    auto is_tail = [&](const Candidate& q) { return rg::none_of(candidates, [&](const Candidate& r) { return thalfs[r.thid].nid_fr() == thalfs[q.thid].nid_to(); }); };

    while (!candidates.empty()) {
        if (auto it = rg::find_if(candidates, is_head); it != candidates.end()) { chains.push_back(consume_pool(thalfs[it->thid].nid_fr(), it->t_fr, false)); continue; }
        if (auto it = rg::find_if(candidates, is_tail); it != candidates.end()) { chains.push_back(consume_pool(thalfs[it->thid].nid_to(), it->t_to, true));  continue; }
        METRIKO_FAIL("chain with thalfs on one side only");
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
            auto it = rg::find(data, old, &Edata::thid); METRIKO_CHECK(it != data.end(), "thalf {} not on its tquad", old);
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

            if (th2.teid == th.teid) return; // when concave tquad attached, this case happens...
            tedges[th2.teid].insert_locs(nids);

            // the zero thalf closed th2 into a loop (a band around a tube): both attachments were valid, so make
            // sure the loop's endpoint is the node where the other tedges meet, not the node being absorbed
            if (auto& loop = tedges[th2.teid].nids; loop.front() == loop.back()) {
                auto deg = [&](int n) {
                    int c = 0;
                    for (auto& [id, ns]: live_tedges()) if (id != th2.teid && (ns.front() == n || ns.back() == n)) ++c;
                    return c;
                };
                const int crnr = deg(nids.front()) >= deg(nids.back()) ? nids.front() : nids.back();
                if (loop.front() != crnr) {
                    loop.pop_back();
                    rg::rotate(loop, rg::find(loop, crnr));
                    loop.push_back(crnr);
                }
            }

            if      (th.bgn)     { (th2.nid_fr() == th.nid_fr()     ? th2 : th2_twn).bgn = true; }
            else if (th_twn.bgn) { (th2.nid_fr() == th_twn.nid_fr() ? th2 : th2_twn).bgn = true; }
            else if (th.end)     { (th2.nid_to() == th.nid_to()     ? th2 : th2_twn).end = true; }
            else if (th_twn.end) { (th2.nid_to() == th_twn.nid_to() ? th2 : th2_twn).end = true; }
        }
    };

    auto on_chain = [&](int nid) { return rg::any_of(chain.pts, [&](const auto& c) { return c.nid == nid; }); };
    auto op = [&](int n0, int n1) { return thids_from_remainning_side(n0, n1) | vw::transform([&](int t) { return thalfs[t].twid; }) | rg::to<vec<int>>(); };
    const auto& th_l = thalfs[chain.thid_l];
    const auto& th_r = thalfs[chain.thid_r];

    auto consume_creates_invalid_tq = [&](const Ehalf& th) {
        Equad copy = tquads[thalfs[th.twid].tqid];
        std::erase_if(copy.data, [&](const Edata& d) { return d.thid == th_l.twid; });
        std::erase_if(copy.data, [&](const Edata& d) { return d.thid == th_r.twid; });
        return !copy.is_valid();
    };

    bool cfr_r  = n_adj_tqs(th_r.id) == 4   && !on_chain(th_r.nid_fr());
    bool cfr_l  = n_adj_tqs(th_l.id) == 4   && !on_chain(th_l.nid_fr());
    bool cto_r  = n_adj_tqs(th_r.twid) == 4 && !on_chain(th_r.nid_to());
    bool cto_l  = n_adj_tqs(th_l.twid) == 4 && !on_chain(th_l.nid_to());
    bool keep_l = consume_creates_invalid_tq(th_l);
    bool keep_r = consume_creates_invalid_tq(th_r);

    if (!thids_bgn.empty() && thids_end.empty()) {
        METRIKO_CHECK(is_top_bgn == is_top_end, "chain ends on diff sides");
        METRIKO_CHECK(thids_bgn.size() == chain.pts.size() - 1, "chain error");
        bool ahd_l = th_l.nid_fr() == nid_bgn;
        bool ahd_r = th_r.nid_fr() == nid_end;
        extend(th_l, ahd_l, cfr_l || cto_l || keep_l);
        extend(th_r, ahd_r, cfr_r || cto_r || keep_r);
        int n0 = ahd_l ? th_l.nid_to() : th_l.nid_fr();
        int n1 = ahd_r ? th_r.nid_to() : th_r.nid_fr();
        replace(op(n0, n1), thids_bgn);
    } else {
        { // leftmost
            bool ahd = th_l.nid_fr() == nid_bgn;
            extend(th_l, ahd, cfr_l || cto_l || keep_l);
            auto n0 = ahd ? th_l.nid_to() : th_l.nid_fr();
            auto na = thalfs[thids_bgn.back()].nid_to();
            auto nb = thalfs[thids_bgn.front()].nid_fr();
            auto n1 = na == th_l.nid_fr() || na == th_l.nid_to() ? nb : na;
            replace(op(n0, n1), thids_bgn);
        }
        { // rightmost
            bool ahd = th_r.nid_fr() == nid_end;
            extend(th_r, ahd, cfr_r || cto_r || keep_r);
            auto n0 = ahd ? th_r.nid_to() : th_r.nid_fr();
            auto na = thalfs[thids_end.back()].nid_to();
            auto nb = thalfs[thids_end.front()].nid_fr();
            auto n1 = na == th_r.nid_fr() || na == th_r.nid_to() ? nb : na;
            replace(op(n0, n1), thids_end);
        }
        // in middle
        for (auto& c: chains) {
            METRIKO_CHECK(c.size() >= 2, "chain error");
            auto n0 = thalfs[c.front()].nid_fr();
            auto n1 = thalfs[c.back()].nid_to();
            replace(op(n0, n1), c);
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
