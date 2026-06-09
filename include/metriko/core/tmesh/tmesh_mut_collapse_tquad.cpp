#include "./tmesh_mut.h"
#include <cassert>
using namespace metriko;

bool TmeshMut::collapse_tquad_prepare(int tqid, Tqaux& tqaux) const {
    auto& tq = tquads[tqid];
    int side = -1;
    if (tq.thids(0).size() == 1 && tq.thids(2).size() == 1 && thalfs[tq.thids(0).front()].x == 0) side = 0;
    if (tq.thids(1).size() == 1 && tq.thids(3).size() == 1 && thalfs[tq.thids(1).front()].x == 0) side = 1;
    if (side == -1) return false;

    auto side_p  = side == 0 ? 0 : 1; // collapse side p
    auto side_q  = side == 0 ? 2 : 3; // collapse side q
    auto side_a  = side == 0 ? 1 : 2; // remain side a
    auto side_b  = side == 0 ? 3 : 0; // remain side b
    auto thid_p  = tq.thids(side_p).front();
    auto thid_q  = tq.thids(side_q).front();
    auto thids_a = tq.thids(side_a);
    auto thids_b = tq.thids(side_b);

    auto find_terminal = [&](int thid, int s0, int s1) -> std::tuple<const HmLoc&, const HmLoc&, int> {
        auto& th0 = thalfs[thid];
        auto& th1 = thalfs[th0.twid];
        if (th0.bgn) return { th0.loc_fr(), th0.loc_to(), s0 };
        if (th1.bgn) return { th1.loc_fr(), th1.loc_to(), s1 };
        if (th0.end) return { th0.loc_to(), th0.loc_fr(), s1 };
        if (th1.end) return { th1.loc_to(), th1.loc_fr(), s0 };
        throw std::runtime_error("error in find_terminal");
    };

    auto [locB, locB_, sideB] = find_terminal(thid_p, side_b, side_a); // loc, another loc and side of the bgn
    auto [locE, locE_, sideE] = find_terminal(thid_q, side_a, side_b); // loc, another loc and side of the end
    auto l = rg::fold_left(thids_a | vw::transform([&](int t){ return thalfs[t].x; }), 0., std::plus{});

    vec<std::tuple<HmLoc, double, int>> aux;
    aux.emplace_back(locB, 0, sideB);
    aux.emplace_back(locE, l, sideE);

    auto v = 0.;
    auto f = [&](auto loc) { return loc != locB && loc != locE && loc != locB_ && loc != locE_; };
    for (int i: thids_a) { auto& th = thalfs[i]; auto lc = th.loc_to(); v += th.x; if (f(lc)) aux.emplace_back(lc, v, side_a); }
    for (int i: thids_b) { auto& th = thalfs[i]; auto lc = th.loc_to(); v -= th.x; if (f(lc)) aux.emplace_back(lc, v, side_b); }

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

    // todo: from find_terminal info, it need to pass data about which dirs(fr/to) of thids_p and q are to be merged into adjacent tquads
    tqaux.checkpoints  = aux;
    tqaux.side_thid_p  = std::pair(side_p, thid_p);
    tqaux.side_thid_q  = std::pair(side_q, thid_q);
    tqaux.side_thids_a = std::pair(side_a, thids_a);
    tqaux.side_thids_b = std::pair(side_b, thids_b);
    return true;
}

void TmeshMut::collapse_tquad_execute(int tqid, Tqaux& tqaux) {
    auto& tq_crr = tquads[tqid];
    auto  region = allowed_range(tq_crr.id);
    vec<std::tuple<int, int, int>> unused_thids; // thid, side of loc_fr, side of loc_to

    for (int i = 0; i < tqaux.checkpoints.size() - 1; ++i) {
        auto& [loc0, val0, side0] = tqaux.checkpoints[i];
        auto& [loc1, val1, side1] = tqaux.checkpoints[i + 1];

        if (side0 == side1) {
            auto it = rg::find_if(thalfs, [&](auto& th) { return th.loc_fr() == loc0 && th.loc_to() == loc1; });
            if (it != thalfs.end()) unused_thids.emplace_back(it->id, side0, side1);
            else throw std::runtime_error("cannot find the thalf");
        } else {
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
            unused_thids.emplace_back(thid0, side0, side1);
            unused_thids.emplace_back(thid1, side1, side0);
        }
    }

    auto& [loc_bgn, val_bgn, side_bgn] = tqaux.checkpoints.front();
    auto& [loc_end, val_end, side_end] = tqaux.checkpoints.back();

    // unused_thids はベクタ順≠接続順なので loc 連結（次は loc_fr==前の loc_to）を辿って seq を伸ばす。
    // break_same=false: side が prev と変わったら停止（同 side run）／true: 同じになったら停止（交互 run）。
    auto extend = [&](vec<int>& seq, HmLoc cur, int prev, bool fwd, bool break_same) {
        while (true) {
            auto it = rg::find_if(unused_thids, [&](auto& t) {
                auto& th = thalfs[std::get<0>(t)];
                return fwd ? th.loc_fr() == cur : th.loc_to() == cur;
            });
            if (it == unused_thids.end()) break;
            int side = fwd ? std::get<1>(*it) : std::get<2>(*it);
            int thid = std::get<0>(*it);
            if ((side == prev) == break_same) break;   // break_same ? s==prev : s!=prev
            prev = side;
            seq.push_back(thid);
            cur = fwd ? thalfs[thid].loc_to() : thalfs[thid].loc_fr();
            unused_thids.erase(it);
        }
    };

    // begin/end: 同 side の連鎖（side が変わったら break）
    // 残り(middle): 交互 side の連鎖に分割（side が同じになったら break）
    vec<int> thids_bgn, thids_end;
    vec<vec<int>> chains;

    extend(thids_bgn, loc_bgn, side_bgn, true,  false);
    extend(thids_end, loc_end, side_end, false, false);

    while (!unused_thids.empty()) {
        auto [thid, side_fr, _] = unused_thids.front();
        unused_thids.erase(unused_thids.begin());
        vec seq{ thid };
        extend(seq, thalfs[thid].loc_to(), side_fr, true, true);
        chains.push_back(seq);
    }

    // extend th and its twin on side_p

    // extend th and its twin on side_q
}

