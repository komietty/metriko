//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tmesh/emesh_postprocess.ipp (see rust/PORTED_FROM).

use crate::common::stdcxx;
use crate::common::utilities::{circular_next};
use crate::hmesh::hmesh::*;
use crate::hmesh::hmloc::*;
use crate::tmesh::emesh::*;
use crate::tmesh::emesh_collapse_path::*;

#[derive(Clone, Copy, Debug)]
struct Cand<'h> {
    nid: i32,
    eid: i32,
    v: Vert<'h>,
    d: f64,
}

impl<'a> Emesh<'a> {
    fn snap_valid(&self, v: Vert<'_>) -> bool {
        for e in self.live_tedges() {
            if e.nids.iter().any(|&nid| matches!(self.tnodes[nid as usize], HmLoc::OnV { id } if id == v.id)) { return false; }
        }
        true
    }

    fn snap_joint(&mut self, teid: i32, v: Vert<'a>) {
        let hm = self.hm;
        let nid = *self.tedges[teid as usize].nids.last().unwrap();
        let HmLoc::OnF { id: lf, xy } = self.tnodes[nid as usize] else { return };
        if !self.snap_valid(v) { return; }
        let loc = HmLoc::OnF { id: lf, xy };

        let f = hm.face(lf);
        for e in self.live_tedges() {
            let nids = &e.nids;
            let mut k = 0;
            while k + 1 < nids.len() {
                if nids[k] == nid || nids[k + 1] == nid { k += 1; continue; }
                if find_strict_intersection_in_face(f, &loc, &HmLoc::OnV { id: v.id }, &self.tnodes[nids[k] as usize], &self.tnodes[nids[k + 1] as usize]) { return; }
                k += 1;
            }
        }

        let mut te_tails: Vec<(i32, i32)> = vec![]; // teid, the biggest  nid
        let mut te_heads: Vec<(i32, i32)> = vec![]; // teid, the smallest nid
        for e in self.live_tedges() {
            if nid == *e.nids.first().unwrap() { te_tails.push((e.id, 0)); }                   // initial value
            if nid == *e.nids.last().unwrap()  { te_heads.push((e.id, e.nids.len() as i32)); } // initial value
        }
        // only the consecutive run of in-star nodes at the joint end is trimmed: not leave and come back case...
        for (i, j) in te_tails.iter_mut() {
            let nids = &self.tedges[*i as usize].nids;
            *j = 0;
            while ((*j + 1) as usize) < nids.len() && is_in_star(v, &self.tnodes[nids[(*j + 1) as usize] as usize]) { *j += 1; }
        }
        for (i, j) in te_heads.iter_mut() {
            let nids = &self.tedges[*i as usize].nids;
            *j = nids.len() as i32 - 1;
            while *j > 0 && is_in_star(v, &self.tnodes[nids[(*j - 1) as usize] as usize]) { *j -= 1; }
        }

        // todo: apply this
        // the legs the snap leaves: from v to the first node each joint tedge keeps. they must not cross or touch any
        // other segment (the trimmed part of their own tedge aside), or the contact is one the repair cannot undo
        // (the C++ keeps this check commented out)

        for &(i, j) in &te_tails {
            let nids = &mut self.tedges[i as usize].nids;
            if j >= 2 { nids.drain(1..j as usize); }
        }
        for &(i, j) in &te_heads {
            let nids = &mut self.tedges[i as usize].nids;
            if ((j + 2) as usize) < nids.len() { let n = nids.len(); nids.drain((j + 1) as usize..n - 1); }
        }
        self.tnodes[nid as usize] = HmLoc::OnV { id: v.id };
    }

    fn snap_inter(&mut self, teid: i32, nid: i32, v: Vert<'a>) -> bool {
        let Some(it) = self.tedges[teid as usize].nids.iter().position(|&x| x == nid) else { return false };
        let i = it as i32;
        let mut lo = i;
        let mut hi = i;
        {
            let nids = &self.tedges[teid as usize].nids;
            while lo > 0                         && is_in_star(v, &self.tnodes[nids[(lo - 1) as usize] as usize]) { lo -= 1; }
            while (hi as usize) < nids.len() - 1 && is_in_star(v, &self.tnodes[nids[(hi + 1) as usize] as usize]) { hi += 1; }
        }
        let l_ = if lo < i { lo } else { i - 1 }; // neighbours after the trim: the two new legs are L -> v -> R
        let r_ = if hi > i { hi } else { i + 1 };

        let ends = {
            let nids = &self.tedges[teid as usize].nids;
            [self.tnodes[nids[l_ as usize] as usize], self.tnodes[nids[r_ as usize] as usize]]
        };
        for l in &ends {
            for f in v.adj_halfs().map(|h| h.face()).filter(|&f| is_in_face(f, l)) {
                for e in self.live_tedges() {
                    let ns = &e.nids;
                    let mut k = 0i32;
                    while ((k + 1) as usize) < ns.len() {
                        if e.id == teid && k >= l_ && k < r_ { k += 1; continue; }
                        if find_strict_intersection_in_face(f, l, &HmLoc::OnV { id: v.id }, &self.tnodes[ns[k as usize] as usize], &self.tnodes[ns[(k + 1) as usize] as usize]) { return false; }
                        k += 1;
                    }
                }
            }
        }

        let nids = &mut self.tedges[teid as usize].nids;
        if hi - i >= 2 { nids.drain((i + 1) as usize..hi as usize); }
        if i - lo >= 2 { nids.drain((lo + 1) as usize..i as usize); }
        self.tnodes[nid as usize] = HmLoc::OnV { id: v.id };
        true
    }

    pub fn collapse_tedge_snap(&mut self, flag: bool) {
        let hm = self.hm;

        // 1: snap joint tnodes. every (joint, face vertex) pair is a candidate,
        //    processed nearest first; the validation happens at apply time because
        //    earlier snaps change the occupancy and the geometry
        let mut cands: Vec<Cand<'a>> = vec![];
        for e in self.live_tedges() {
            let l = &self.tnodes[*e.nids.last().unwrap() as usize];
            let HmLoc::OnF { id: lf, .. } = *l else { continue };
            let p = get_ptloc_pos(hm, l);
            for v in hm.face(lf).verts() {
                // the designated initializer omits .nid: value-initialized to 0
                cands.push(Cand { nid: 0, eid: e.id, v, d: (v.pos() - p).length_squared() });
            }
        }
        stdcxx::sort_by(&mut cands, |a, b| a.d < b.d);
        for c in &cands { self.snap_joint(c.eid, c.v); }
        let mut candidates: Vec<Vec<Cand<'a>>> = vec![vec![]; self.tnodes.len()];

        // 2: snap inter tnodes
        for e in self.live_tedges() {
            let nids = &e.nids;
            let mut i = 1usize;
            while i < nids.len() - 1 {
                let nid = nids[i];
                let n = self.tnodes[nid as usize];
                let p = get_ptloc_pos(hm, &n);

                if matches!(n, HmLoc::OnV { .. }) { i += 1; continue; }
                for vid in get_ptloc_verts(hm, &n) {
                    let v = hm.vert(vid);
                    let d = (v.pos() - p).length_squared();
                    candidates[nid as usize].push(Cand { nid, eid: e.id, v, d });
                }
                i += 1;
            }
        }

        // 3: sort candidates in inner/outer order
        for c in candidates.iter_mut() { stdcxx::sort_by(c, |a, b| a.d < b.d); }
        let key = |c: &Vec<Cand<'_>>| if c.is_empty() { 1e9 } else { c[0].d };
        stdcxx::sort_by(&mut candidates, |a, b| key(a) < key(b));

        // 4: snap if it's valid
        if flag {
            for c in &candidates {
                for cd in c {
                    if self.snap_valid(cd.v) && self.snap_inter(cd.eid, cd.nid, cd.v) { break; }
                }
            }
        } else {
            for c in &candidates {
                if c.is_empty() { continue; }
                let cd = c[0];
                if self.snap_valid(cd.v) { self.snap_inter(cd.eid, cd.nid, cd.v); }
            }
        }
    }

    pub fn collapse_tedge_snap_dedup(&mut self, teid: i32) {
        let hm = self.hm;

        // nodes shared with other chains (junctions / crossings) must stay: dropping
        // one from this chain only would let the straightened chord cross the other
        // chain inside a face without a shared vertex
        let mut shared = vec![false; self.tnodes.len()];
        for e in self.live_tedges() {
            if e.id == teid { continue; }
            for &nid in &e.nids { shared[nid as usize] = true; }
        }

        let mut i = 1usize;
        while i < self.tedges[teid as usize].nids.len() - 1 {
            let nids = &self.tedges[teid as usize].nids;
            let nid_prev = nids[i - 1];
            let nid_curr = nids[i];
            let nid_next = nids[i + 1];

            if shared[nid_curr as usize] || nid_prev == nid_next { i += 1; continue; }

            // a face carrying all three nodes: the chord prev-next stays inside it
            let fids_curr = get_ptloc_faces(hm, &self.tnodes[nid_curr as usize]);
            let fids_prev = get_ptloc_faces(hm, &self.tnodes[nid_prev as usize]);
            let fids_next = get_ptloc_faces(hm, &self.tnodes[nid_next as usize]);

            let same_face = fids_curr.iter().any(|fid| fids_prev.contains(fid) && fids_next.contains(fid));
            if !same_face { i += 1; continue; }

            // drop the middle node; keep the removal only if the validation holds
            self.tedges[teid as usize].nids.remove(i);
        }
    }

    // retrace a tedge inside its two tquads together, to resolve a contact made by snapping
    pub fn reroute_tedge(&mut self, teid: i32) -> bool {
        let Some(th) = self.thalfs.iter().position(|th| th.id != -1 && th.teid == teid && th.cano) else { return false };
        let th = self.thalfs[th];
        if th.tqid < 0 || th.twin(self).tqid < 0 { return false; }
        let (front, back) = (*self.tedges[teid as usize].nids.first().unwrap(), *self.tedges[teid as usize].nids.last().unwrap());
        if front == back { return false; } // a closed tedge: no two ends to join

        // the boundary of the two tquads together, the tedge itself left out: from its head around to its tail
        let rest = |tq: &Equad, thid: i32| -> Vec<i32> {
            let mut res = vec![];
            let n = tq.data.len();
            let mut it = circular_next(n, tq.find_of(thid).unwrap_or(n));
            while tq.data[it].thid != thid { res.push(tq.data[it].thid); it = circular_next(n, it); }
            res
        };
        let mut lp = rest(th.tquad(self), th.id);
        let tail = lp.len() as i32;
        for t in rest(th.twin(self).tquad(self), th.twid) { lp.push(t); }
        let chart = unfold(self, &lp);
        let path = retrace(self, &chart, chart.at(tail, front), chart.at(0, back));
        if path.is_empty() { return false; }
        let nids = self.add_new_path(&path, front, back);
        self.tedges[teid as usize].nids = nids;

        let r = self.path_length(&self.tedges[teid as usize].nids);
        for t in self.thalfs.iter_mut() { if t.id != -1 && t.teid == teid { t.r = r; } }
        true
    }
}
