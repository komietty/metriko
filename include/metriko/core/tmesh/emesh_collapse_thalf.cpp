//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include "emesh.h"
using namespace metriko;

void Emesh::collapse_thalf(int thid) {
    auto& th_crr = thalfs[thid];
    auto& th_twn = thalfs[th_crr.twid];
    auto& tq_crr = tquads[th_crr.tqid];
    auto& tq_twn = tquads[th_twn.tqid];
    auto  it_crr = rg::find(tq_crr.data, thid     , &Edata::thid);
    auto  it_twn = rg::find(tq_twn.data, th_twn.id, &Edata::thid);
    auto  it_prv = circular_prev(tq_crr.data, it_crr);
    auto  it_nxt = circular_next(tq_crr.data, it_crr);
    auto& th_prv = thalfs[it_prv->thid];
    auto& th_nxt = thalfs[it_nxt->thid];
    auto& te_crr = tedges[th_crr.teid];
    auto& te_prv = tedges[th_prv.teid];
    auto& te_nxt = tedges[th_nxt.teid];

    if(it_crr->side != it_prv->side && it_crr->side != it_nxt->side) throw std::runtime_error("collapse thalf error");
    if(it_crr->side == it_prv->side && it_crr->side == it_nxt->side) throw std::runtime_error("collapse thalf error");

    int cout_fr = count_adj_tquads(th_crr.id);
    int cout_to = count_adj_tquads(th_twn.id);

    if      (cout_fr == 2) { te_prv.insert_locs(te_crr.nids); }
    else if (cout_to == 2) { te_nxt.insert_locs(te_crr.nids); }
    else {
        bool collapse_to_prev = it_crr->side != it_prv->side;

        auto [n_fr, n_to] = [&]{
            if ( collapse_to_prev &&  th_prv.cano) return std::pair{th_prv.nid_fr(), th_crr.nid_to()};
            if ( collapse_to_prev && !th_prv.cano) return std::pair{th_crr.nid_to(), th_prv.nid_fr()};
            if (!collapse_to_prev &&  th_nxt.cano) return std::pair{th_crr.nid_fr(), th_nxt.nid_to()};
            if (!collapse_to_prev && !th_nxt.cano) return std::pair{th_nxt.nid_to(), th_crr.nid_fr()};
            throw std::runtime_error("no impl");
        }();

        auto region = allowed_range_tquads({tq_crr.id});
        auto path   = approx_shortest_path(30, hm, tnodes[n_fr], tnodes[n_to], region);
        auto nids   = add_new_path(path, n_fr, n_to);

        { // Count duplication of verts for checking loc-inj
            vec<HmLoc> locs;
            for (const auto& [id, side]: tq_crr.data) {
                auto ns = tedges[thalfs[id].teid].nids;
                if (!thalfs[id].cano) rg::reverse(ns);
                for (size_t k = 0; k + 1 < ns.size(); ++k) locs.push_back(tnodes[ns[k]]);
            }

            int dups = 0;
            for (size_t i = 0; i < locs.size(); ++i)
            for (size_t j = 0; j < i; ++j)
                if (locs[i] == locs[j]) { dups++; break; }

            // todo: corner case! not loc-inj from tqad to hm verts. Dups == 1 is huge bottle neck for corse quad!
            // if dups >= 2, tquad is self intersected by one of its thalfs
            // if dups == 1, tquad is self intersected by one of its singular (now tempolary skip, and hope aother tquad is loc-inj)
            if (dups >= 2) throw std::runtime_error("collapse thalf error");
            if (dups == 1) {
                auto  r = 0.5;
                auto& te_tgt = collapse_to_prev ? te_prv : te_nxt;
                if (path_length(nids) < r * (path_length(te_crr.nids) + path_length(te_tgt.nids))) return;
            }
        }

        auto  it_twn_adj = collapse_to_prev ? circular_next(tq_twn.data, it_twn) : circular_prev(tq_twn.data, it_twn);
        auto& th_twn_adj = thalfs[it_twn_adj->thid];
        auto& te_twn_adj = tedges[th_twn_adj.teid];
        if(it_twn_adj->side != it_twn->side) throw std::runtime_error("collapse thalf error");

        // 1: collapse to the prv/nxt edge
        // 2: insert missing segments to the adjacent edge.
        if (collapse_to_prev) { te_prv.nids = nids; te_twn_adj.insert_locs(te_crr.nids); }
        else                  { te_nxt.nids = nids; te_twn_adj.insert_locs(te_crr.nids); }
    }

    // remove the data from tquads which have collapsed thalfs
    tq_crr.data.erase(it_crr);
    tq_twn.data.erase(it_twn);
    th_crr = {};
    th_twn = {};
    te_crr = {};
}
