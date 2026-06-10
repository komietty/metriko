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
    vec<std::tuple<int, int, int>> unused_thids; // thid, side of loc_fr, side of loc_to

    for (int i = 0; i < tqaux.checkpoints.size() - 1; ++i) {
        auto& [loc0, val0, side0] = tqaux.checkpoints[i];
        auto& [loc1, val1, side1] = tqaux.checkpoints[i + 1];

        if (side0 == side1) {
            for (auto& [thid, side]: tq_crr.data) {
                auto& th = thalfs[thid];
                if ((th.loc_fr() == loc0 && th.loc_to() == loc1) ||
                    (th.loc_fr() == loc1 && th.loc_to() == loc0))
                    unused_thids.emplace_back(thid, side0, side1);
            }
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
    auto extend = [&](vec<int>& seq, HmLoc cur, int prev, bool fwd, bool break_same, bool err_if_none) {
        while (true) {
            auto it = rg::find_if(unused_thids, [&](auto& t) {
                auto& th = thalfs[std::get<0>(t)];
                return fwd ? th.loc_fr() == cur : th.loc_to() == cur;
            });
            if (it == unused_thids.end()) {
                // bgn: side が変わるまで必ず次が見つかるはず。見つからなければ連結が壊れている。
                if (err_if_none) throw std::runtime_error("extend: no connecting thalf found before side change");
                break;
            }
            // 同 side run(break_same=false)は「進行先（次の cur）」の side で判定する：
            //   fwd → 次は loc_to(get<2>),  bwd → 次は loc_fr(get<1>)。
            // 交互 run(middle, break_same=true)は接続点 loc_fr(get<1>) の side で交互を判定。
            int side = break_same ? std::get<1>(*it)
                                  : (fwd ? std::get<2>(*it) : std::get<1>(*it));
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

    // bgn: loc_bgn から forward に辿り、to-side が初期 side_bgn の間は取り込み、
    //      side が変わる辺（コーナーに到達する cross）も含めて停止。
    //      接続先が尽きたら（既に他方が消費）そこで停止（throw しない）。
    {
        HmLoc cur = loc_bgn;
        while (true) {
            auto it = rg::find_if(unused_thids, [&](auto& t){ return thalfs[std::get<0>(t)].loc_fr() == cur; });
            if (it == unused_thids.end()) break;
            int  thid = std::get<0>(*it);
            bool flag = std::get<2>(*it) != side_bgn;   // to-side が初期 side から変わったか
            thids_bgn.push_back(thid);
            cur = thalfs[thid].loc_to();
            unused_thids.erase(it);
            if (flag) break;
        }
    }
    // end: 対称（backward, fr-side で判定）。
    {
        HmLoc cur = loc_end;
        while (true) {
            auto it = rg::find_if(unused_thids, [&](auto& t){ return thalfs[std::get<0>(t)].loc_fr() == cur; });
            if (it == unused_thids.end()) break;
            int  thid = std::get<0>(*it);
            bool flag = std::get<2>(*it) != side_end;   // fr-side が初期 side から変わったか
            thids_end.push_back(thid);
            cur = thalfs[thid].loc_to();
            unused_thids.erase(it);
            if (flag) break;
        }
    }

    while (!unused_thids.empty()) {
        auto [thid, side_fr, _] = unused_thids.front();
        unused_thids.erase(unused_thids.begin());
        vec seq{ thid };
        extend(seq, thalfs[thid].loc_to(), side_fr, true, true, false);
        chains.push_back(seq);
    }

    // tquad.data 順（CCW 境界順）を辿る next/prev。nxid/pvid フィールドの代わりに使う。
    auto step_next = [&](int thid) { auto& q = tquads[thalfs[thid].tqid]; return circular_next(q.data, rg::find(q.data, thid, &TdataMut::thid))->thid; };
    auto step_prev = [&](int thid) { auto& q = tquads[thalfs[thid].tqid]; return circular_prev(q.data, rg::find(q.data, thid, &TdataMut::thid))->thid; };

    ThalfMut& th_r = thalfs[tqaux.side_thid_r.second];
    ThalfMut& th_l = thalfs[tqaux.side_thid_l.second];
    TedgeMut& te_r = tedges[th_r.teid];
    TedgeMut& te_l = tedges[th_l.teid];

    // extend th and its twin on side_r
    {
        ThalfMut& th_nxt = thalfs[step_next(th_r.id)];
        ThalfMut& th_prv = thalfs[step_prev(th_r.id)];
        ThalfMut& th_ahd = thalfs[step_next(th_nxt.twid)];
        ThalfMut& th_bhd = thalfs[step_prev(th_prv.twid)];
        TedgeMut& te_ahd = tedges[th_ahd.teid];
        TedgeMut& te_bhd = tedges[th_bhd.teid];
        std::cout << "te_r     id: " << th_r.teid   << std::endl;
        std::cout << "te_r_ahd id: " << th_ahd.teid << std::endl;
        std::cout << "te_r_bhd id: " << th_bhd.teid << std::endl;
        std::cout << "merge to ahead: " << tqaux.thid_r_merge_to_ahead << std::endl;
        vec<int> r_wo_last (te_r.nids.begin(), te_r.nids.end() - 1);
        vec<int> r_wo_first(te_r.nids.begin() + 1, te_r.nids.end());
        bool cano = th_r.cano;
        if (tqaux.thid_r_merge_to_ahead) {
            if (cano) te_ahd.insert_locs_front(r_wo_last);
            else      te_ahd.insert_locs_after(r_wo_first);
        } else {
            if (cano) te_bhd.insert_locs_after(r_wo_first);
            else      te_bhd.insert_locs_front(r_wo_last);

            ThalfMut& th_replace = thalfs[th_prv.twid];
            TquadMut& tq_replace = tquads[th_replace.tqid];
            // find side of th_replace
            int side = tq_replace.side_of(th_replace);
            // replace all the data on `side` with thids_bgn（CCW で side は連続なので区間置換）
            auto& data = tq_replace.data;
            int pos = rg::find(data, side, &TdataMut::side) - data.begin();
            std::erase_if(data, [&](const TdataMut& d) { return d.side == side; });
            vec<TdataMut> repl;
            for (int thid : thids_bgn) { thalfs[thid].tqid = tq_replace.id; repl.push_back({ thid, side }); }
            data.insert(data.begin() + pos, repl.begin(), repl.end());
        }
    }

    {
        // extend th and its twin on side_l
        ThalfMut& th_nxt = thalfs[step_next(th_l.id)];
        ThalfMut& th_prv = thalfs[step_prev(th_l.id)];
        ThalfMut& th_ahd = thalfs[step_next(th_nxt.twid)];
        ThalfMut& th_bhd = thalfs[step_prev(th_prv.twid)];
        TedgeMut& te_ahd = tedges[th_ahd.teid];
        TedgeMut& te_bhd = tedges[th_bhd.teid];
        std::cout << "te_l     id: " << th_l.teid << std::endl;
        std::cout << "te_l_ahd id: " << th_ahd.teid << std::endl;
        std::cout << "te_l_bhd id: " << th_bhd.teid << std::endl;
        vec<int> l_wo_last(te_l.nids.begin(), te_l.nids.end() - 1);
        vec<int> l_wo_first(te_l.nids.begin() + 1, te_l.nids.end());
        bool cano = th_l.cano;

        if (tqaux.thid_l_merge_to_ahead) {
            if (cano) te_ahd.insert_locs_front(l_wo_last);
            else      te_ahd.insert_locs_after(l_wo_first);
        } else {
            if (cano) te_bhd.insert_locs_after(l_wo_first);
            else      te_bhd.insert_locs_front(l_wo_last);

            ThalfMut& th_replace = thalfs[th_prv.twid];
            TquadMut& tq_replace = tquads[th_replace.tqid];

            int side = tq_replace.side_of(th_replace);
            // replace all the data on `side` with thids_bgn（CCW で side は連続なので区間置換）
            auto& data = tq_replace.data;
            int pos = rg::find(data, side, &TdataMut::side) - data.begin();
            std::erase_if(data, [&](const TdataMut& d) { return d.side == side; });
            vec<TdataMut> repl;
            for (int thid : thids_end) { thalfs[thid].tqid = tq_replace.id; repl.push_back({ thid, side }); }
            data.insert(data.begin() + pos, repl.begin(), repl.end());

        }
    }

    te_r.id = -1;
    te_l.id = -1;
    tq_crr.id = -1;
}

