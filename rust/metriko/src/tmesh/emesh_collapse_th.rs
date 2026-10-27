//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tmesh/emesh_collapse_th.ipp (see rust/PORTED_FROM).

use crate::common::utilities::{circular_next, circular_prev};
use crate::hmesh::hmloc::HmLoc;
use crate::tmesh::emesh::*;
use crate::tmesh::emesh_collapse_path::*;

impl<'a> Emesh<'a> {
    pub fn collapse_thalf(&mut self, thid: i32) {
        // the C++ binds references to the elements; their ids are fixed for the whole function
        let th_crr = self.thalfs[thid as usize];
        let th_twn = self.thalfs[th_crr.twid as usize];
        let tq_crr = th_crr.tqid as usize;
        let tq_twn = th_twn.tqid as usize;
        let it_crr = self.tquads[tq_crr].data.iter().position(|d| d.thid == thid).expect("rg::find");
        let it_twn = self.tquads[tq_twn].data.iter().position(|d| d.thid == th_twn.id).expect("rg::find");
        let n_crr = self.tquads[tq_crr].data.len();
        let it_prv = circular_prev(n_crr, it_crr);
        let it_nxt = circular_next(n_crr, it_crr);
        let d_crr = self.tquads[tq_crr].data[it_crr];
        let d_prv = self.tquads[tq_crr].data[it_prv];
        let d_nxt = self.tquads[tq_crr].data[it_nxt];
        let th_prv = self.thalfs[d_prv.thid as usize];
        let th_nxt = self.thalfs[d_nxt.thid as usize];
        let te_crr = th_crr.teid as usize;
        let te_prv = th_prv.teid as usize;
        let te_nxt = th_nxt.teid as usize;

        metriko_check!(d_crr.side == d_prv.side || d_crr.side == d_nxt.side, "collapse thalf error");
        metriko_check!(d_crr.side != d_prv.side || d_crr.side != d_nxt.side, "collapse thalf error");

        if self.n_adj_tqs(th_crr.id) == 2 && !self.is_fixed(th_crr.nid_fr(self)) {
            let ns = self.tedges[te_crr].nids.clone();
            self.tedges[te_prv].insert_locs(&ns);
        } else if self.n_adj_tqs(th_twn.id) == 2 && !self.is_fixed(th_crr.nid_to(self)) {
            let ns = self.tedges[te_crr].nids.clone();
            self.tedges[te_nxt].insert_locs(&ns);
        } else {
            let mut collapse_to_prev = d_crr.side != d_prv.side;
            if self.is_fixed(if collapse_to_prev { th_crr.nid_fr(self) } else { th_crr.nid_to(self) }) { collapse_to_prev = !collapse_to_prev; }

            // the merged tedge joins the far end of the neighbour (own) to the far end of te_crr (oth) by the shortest
            // line inside tq_crr, ordered like the neighbour
            let tq_thids: Vec<i32> = self.tquads[tq_crr].data.iter().map(|d| d.thid).collect();
            let chart = unfold(self, &tq_thids);
            let at = |it: usize, nid: i32| chart.at(it as i32, nid);
            let own = if collapse_to_prev { at(it_prv, th_prv.nid_fr(self)) } else { at(it_nxt, th_nxt.nid_to(self)) };
            let oth = if collapse_to_prev { at(it_crr, th_crr.nid_to(self)) } else { at(it_crr, th_crr.nid_fr(self)) };
            let (i_fr, i_to) = if collapse_to_prev == (if collapse_to_prev { th_prv } else { th_nxt }).cano { (own, oth) } else { (oth, own) };
            let path = retrace(self, &chart, i_fr, i_to);
            metriko_check!(!path.is_empty(), "no path through tquad {} between nodes {} and {}", self.tquads[tq_crr].id, chart.nids[i_fr as usize], chart.nids[i_to as usize]);
            let (n_fr, n_to) = (chart.nids[i_fr as usize], chart.nids[i_to as usize]);
            let nids = self.add_new_path(&path, n_fr, n_to);

            { // Count duplication of verts for checking loc-inj
                let mut locs: Vec<HmLoc> = vec![];
                for d in &self.tquads[tq_crr].data {
                    let id = d.thid as usize;
                    let mut ns = self.thalfs[id].tedge(self).nids.clone();
                    if !self.thalfs[id].cano { ns.reverse(); }
                    let mut k = 0;
                    while k + 1 < ns.len() { locs.push(self.tnodes[ns[k] as usize]); k += 1; }
                }

                let mut dups = 0;
                for i in 0..locs.len() {
                    for j in 0..i {
                        if locs[i] == locs[j] { dups += 1; break; }
                    }
                }

                // if dups >= 2, tquad is self intersected by one of its thalfs
                // if dups == 1, tquad is self intersected by one of its singular (temp skip, hope aother tquad is loc-inj)
                metriko_check!(dups < 2, "tquad {} is self-intersected by one of its thalfs", self.tquads[tq_crr].id);
                if dups == 1 {
                    let te_tgt = if collapse_to_prev { te_prv } else { te_nxt };
                    if self.path_length(&nids) < 0.5 * (self.path_length(&self.tedges[te_crr].nids) + self.path_length(&self.tedges[te_tgt].nids)) { return; }
                }
            }

            let n_twn = self.tquads[tq_twn].data.len();
            let it_twn_adj = if collapse_to_prev { circular_next(n_twn, it_twn) } else { circular_prev(n_twn, it_twn) };
            let th_twn_adj = self.thalfs[self.tquads[tq_twn].data[it_twn_adj].thid as usize];
            let te_twn_adj = th_twn_adj.teid as usize;
            metriko_check!(self.tquads[tq_twn].thids(self.tquads[tq_twn].data[it_twn].side).len() >= 2, "its twin is alone on a side of tquad");

            // 1: collapse to the prv/nxt edge
            // 2: insert missing segments to the adjacent edge.
            if collapse_to_prev {
                self.tedges[te_prv].nids = nids;
                let ns = self.tedges[te_crr].nids.clone();
                self.tedges[te_twn_adj].insert_locs(&ns);
            } else {
                self.tedges[te_nxt].nids = nids;
                let ns = self.tedges[te_crr].nids.clone();
                self.tedges[te_twn_adj].insert_locs(&ns);
            }
        }

        // remove the data from tquads which have collapsed thalfs
        // (vector::erase through the iterators taken at the start: positions, as the C++ pointers)
        self.tquads[tq_crr].data.remove(it_crr);
        self.tquads[tq_twn].data.remove(it_twn);
        self.thalfs[thid as usize] = Ehalf::default();
        self.thalfs[th_crr.twid as usize] = Ehalf::default();
        self.tedges[te_crr] = Eedge::default();
    }
}
