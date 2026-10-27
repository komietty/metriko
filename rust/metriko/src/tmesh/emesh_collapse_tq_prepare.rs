//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tmesh/emesh_collapse_tq_prepare.ipp (see rust/PORTED_FROM).

use crate::tmesh::emesh::*;

impl<'a> Emesh<'a> {
    pub fn collapse_tquad_prepare(&self, tqid: i32, aux: &mut Tqaux) -> bool {
        let tq = &self.tquads[tqid as usize];
        if tq.id == -1 { return false; }

        let zero = |thid: i32| self.thalfs[thid as usize].x == 0.;
        let sum_x = |thids: &[i32]| { let mut s = 0.; for &t in thids { s += self.thalfs[t as usize].x; } s };
        let sum_r = |thids: &[i32]| { let mut s = 0.; for &t in thids { s += self.thalfs[t as usize].r; } s };

        let mut side = -1;
        if tq.thids(0).len() == 1 && tq.thids(2).len() == 1 && zero(tq.thids(0)[0]) { side = 0; }
        if tq.thids(1).len() == 1 && tq.thids(3).len() == 1 && zero(tq.thids(1)[0]) { side = 1; }
        if side == -1 { return false; }

        let s_r = side;
        let s_t = side + 1;
        let s_l = side + 2;
        let s_b = (side + 3) % 4;

        let thids_r = tq.thids(s_r);
        let thids_l = tq.thids(s_l);
        let thids_t = tq.thids(s_t);
        let thids_b = tq.thids(s_b);
        // `int sum_x_* = sumX(...)`: the double is truncated to int
        let sum_x_r = sum_x(&thids_r) as i32;
        let sum_x_l = sum_x(&thids_l) as i32;
        let sum_x_t = sum_x(&thids_t) as i32;
        let sum_x_b = sum_x(&thids_b) as i32;
        let point = thids_r.len() == 1 && sum_x_r == 0 &&
                    thids_l.len() == 1 && sum_x_l == 0 &&
                    thids_t.len() == 1 && sum_x_t == 0 &&
                    thids_b.len() == 1 && sum_x_b == 0;
        if !point && (sum_x_t == 0 || sum_x_b == 0) { return false; }

        aux.tqid    = tqid;
        aux.thid_r  = thids_r[0];
        aux.thid_l  = thids_l[0];
        aux.thids_t = thids_t;
        aux.thids_b = thids_b;
        let th_l = self.thalfs[aux.thid_l as usize];
        let th_r = self.thalfs[aux.thid_r as usize];

        let find_terminal = |thid: i32, s0: bool, s1: bool| -> (i32, bool) {
            let th0 = &self.thalfs[thid as usize];
            let th1 = th0.twin(self);
            if self.is_fixed(th0.nid_fr(self)) { return (th0.nid_fr(self), s0); }
            if self.is_fixed(th0.nid_to(self)) { return (th0.nid_to(self), s1); }
            if th0.bgn { return (th0.nid_fr(self), s0); }
            if th1.bgn { return (th1.nid_fr(self), s1); }
            if th0.end { return (th0.nid_to(self), s1); }
            if th1.end { return (th1.nid_to(self), s0); }
            if self.n_adj_tqs(th0.id) >= self.n_adj_tqs(th0.twid) { (th0.nid_fr(self), s0) } else { (th0.nid_to(self), s1) }
        };

        let must_use_in_aux = |nid: i32, thid: i32| -> bool { // a crossing point or a singular point, not tq corners
            if [th_l.nid_fr(self), th_l.nid_to(self), th_r.nid_fr(self), th_r.nid_to(self)].contains(&nid) { return false; }
            self.n_adj_tqs(thid) != 2 || self.is_fixed(nid)
        };

        let mut oft = 0i32;
        let mut btm = 0i32;
        let mut rt = 0.;
        let mut rb = 0.;
        let span = sum_x(&aux.thids_t);
        let rt_sum = sum_r(&aux.thids_t);
        let rb_sum = sum_r(&aux.thids_b);
        // `oft += th.x` on an int: the sum is computed in double and truncated
        for &thid in aux.thids_t.iter().rev() {
            let th = &self.thalfs[thid as usize];
            oft = (oft as f64 + th.x) as i32;
            rt += th.r;
            if must_use_in_aux(th.nid_fr(self), th.id) { aux.pts.push(Tqpoint { nid: th.nid_fr(self), val: oft, ord: rt / rt_sum * span, top: true }); }
        }
        for &thid in aux.thids_b.iter() {
            let th = &self.thalfs[thid as usize];
            btm = (btm as f64 + th.x) as i32;
            rb += th.r;
            if must_use_in_aux(th.nid_to(self), th.twid) { aux.pts.push(Tqpoint { nid: th.nid_to(self), val: btm, ord: rb / rb_sum * span, top: false }); }
        }

        let (nid_bgn, side_bgn) = find_terminal(aux.thid_l, true, false);
        let (nid_end, side_end) = find_terminal(aux.thid_r, false, true);
        aux.pts.push(Tqpoint { nid: nid_bgn, val: 0, ord: 0., top: side_bgn });
        aux.pts.push(Tqpoint { nid: nid_end, val: oft, ord: oft as f64, top: side_end });
        // rg::stable_sort by std::pair(val, ord)
        aux.pts.sort_by(|a, b| {
            if a.val < b.val || (a.val == b.val && a.ord < b.ord) { std::cmp::Ordering::Less }
            else if b.val < a.val || (b.val == a.val && b.ord < a.ord) { std::cmp::Ordering::Greater }
            else { std::cmp::Ordering::Equal }
        });

        if !point {
            for thids in [&aux.thids_t, &aux.thids_b] {
                for &thid in thids.iter().filter(|&&i| zero(i)) {
                    let p0 = aux.pts.iter().position(|p| p.nid == self.thalfs[thid as usize].nid_fr(self));
                    let p1 = aux.pts.iter().position(|p| p.nid == self.thalfs[thid as usize].nid_to(self));
                    let (Some(p0), Some(p1)) = (p0, p1) else { return false; };
                    if (p0 as i64 - p1 as i64).abs() != 1 { return false; }
                    let (v0, t0) = (aux.pts[p0].val, aux.pts[p0].top);
                    if aux.pts.iter().any(|p| p.val == v0 && p.top != t0) { return false; }
                }
            }
        }

        true
    }
}
