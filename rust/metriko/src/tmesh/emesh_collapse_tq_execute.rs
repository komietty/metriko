//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tmesh/emesh_collapse_tq_execute.ipp (see rust/PORTED_FROM).

use crate::tmesh::emesh::*;
use crate::tmesh::emesh_collapse_path::*;

#[derive(Clone, Copy, Debug)]
struct Piece {
    thid: i32,
    t_fr: bool,
    t_to: bool,
}

// The C++ holds `th_l` / `th_r` as references into `thalfs` (and `extend` / `keep_of` take `const Ehalf&`): every read
// sees the current state of the element. Here they are ids, read from `self.thalfs` at each use.
// NOTE: in the C++, `thalfs.emplace_back` in the pieces loop may reallocate and leave those references dangling
// (undefined behavior); this port implements the intended reading, the live element.

impl<'a> Emesh<'a> {
    fn exe_lateral_run(&self, aux: &Tqaux, fr: i32, to: i32) -> Vec<i32> {
        for (a, b) in [(fr, to), (to, fr)] {
            for s in [&aux.thids_t, &aux.thids_b] {
                let mut res = vec![];
                let mut c = a;
                while c != b {
                    let Some(&t) = s.iter().find(|&&t| self.thalfs[t as usize].nid_fr(self) == c) else { break };
                    res.push(t);
                    c = self.thalfs[t as usize].nid_to(self);
                }
                if c == b {
                    if a != fr { res.reverse(); }
                    return res;
                }
            }
        }
        vec![]
    }

    fn exe_take_run(&self, pieces: &mut Vec<Piece>, nid: i32, side_to_stop: bool, backwards: bool) -> Vec<i32> {
        let mut res = vec![];
        let mut cur = nid;
        loop {
            let Some(it) = pieces.iter().position(|p| {
                (if backwards { self.thalfs[p.thid as usize].nid_to(self) } else { self.thalfs[p.thid as usize].nid_fr(self) }) == cur
            }) else { break };
            let p = pieces[it];
            let stop = if backwards { p.t_fr == side_to_stop } else { p.t_to == side_to_stop };
            res.push(p.thid);
            cur = if backwards { self.thalfs[p.thid as usize].nid_fr(self) } else { self.thalfs[p.thid as usize].nid_to(self) };
            pieces.remove(it);
            if stop { break; }
        }
        if backwards { res.reverse(); }
        res
    }

    fn exe_replace(&mut self, olds: &[i32], run: &[i32]) {
        let mut ri = 0usize;
        for &old in olds {
            let e0 = self.thalfs[old as usize].nid_fr(self);
            let e1 = self.thalfs[old as usize].nid_to(self);
            let mut seg = vec![];
            while ri < run.len() {
                seg.push(run[ri]);
                let to = self.thalfs[run[ri] as usize].nid_to(self);
                ri += 1;
                if ri >= run.len() || to == e0 || to == e1 { break; }
            }
            let tq = self.thalfs[old as usize].tqid as usize;
            let id = self.tquads[tq].id;
            let it = self.tquads[tq].data.iter().position(|d| d.thid == old);
            metriko_check!(it.is_some(), "thalf {} not on its tquad", old);
            let mut it = it.unwrap();
            let side = self.tquads[tq].data[it].side;
            self.tquads[tq].data.remove(it);
            for &t in &seg {
                self.thalfs[t as usize].tqid = id;
                self.tquads[tq].data.insert(it, Edata { thid: t, side });
                it += 1;
            }
        }
    }

    fn exe_outer(&self, aux: &Tqaux, n0: i32, n1: i32) -> Vec<i32> {
        let run = self.exe_lateral_run(aux, n0, n1);
        metriko_check!(!run.is_empty(), "no thalf run on a lateral side from node {} to node {}", n0, n1);
        run.iter().map(|&t| self.thalfs[t as usize].twid).collect()
    }

    /// `extend(const Ehalf& th, bool ahd, bool keep)`, `th` = `thalfs[thid]`
    fn exe_extend(&mut self, thid: i32, ahd: bool, keep: bool) {
        let step = |tm: &Emesh<'_>, t: i32| if ahd { tm.step_next(t) } else { tm.step_prev(t) };
        let th = thid as usize;
        let th1 = step(self, self.thalfs[th].id) as usize;
        let th2 = step(self, self.thalfs[th1].twid) as usize;
        let tq2 = self.thalfs[th2].tqid as usize;
        let id = self.tquads[tq2].id;
        if keep {
            let it = self.tquads[tq2].data.iter().position(|d| d.thid == self.thalfs[th2].id).expect("rg::find");
            let th_id_ = self.thalfs[th].id as usize;
            self.thalfs[th_id_].tqid = id;
            let side = self.tquads[tq2].data[it].side;
            let th_id = self.thalfs[th].id;
            self.tquads[tq2].data.insert(if ahd { it } else { it + 1 }, Edata { thid: th_id, side });
            return;
        }
        // nids: a reference to tedges[th.teid].nids; th_twn: thalfs[th.twid]; th2_twn: thalfs[th2.twid]
        let te = self.thalfs[th].teid as usize;
        let th_twn = self.thalfs[th].twid as usize;
        let th2_twn = self.thalfs[th2].twid as usize;
        {
            let (q, i) = (self.thalfs[th_twn].tqid as usize, self.thalfs[th_twn].id);
            self.tquads[q].data.retain(|d| d.thid != i);
        }

        if self.thalfs[th2].teid == self.thalfs[th].teid { return; }
        let te2 = self.thalfs[th2].teid as usize;
        let nids = self.tedges[te].nids.clone();
        self.tedges[te2].insert_locs(&nids);

        if self.tedges[te2].nids.first() == self.tedges[te2].nids.last() {
            let deg = |tm: &Emesh<'_>, n: i32| {
                let mut c = 0;
                for e in tm.live_tedges() {
                    if e.id != tm.thalfs[th2].teid && (*e.nids.first().unwrap() == n || *e.nids.last().unwrap() == n) { c += 1; }
                }
                c
            };
            let (nf, nb) = (*self.tedges[te].nids.first().unwrap(), *self.tedges[te].nids.last().unwrap());
            let crnr = if deg(self, nf) >= deg(self, nb) { nf } else { nb };
            let lp = &mut self.tedges[te2].nids;
            if *lp.first().unwrap() != crnr {
                lp.pop();
                let pos = lp.iter().position(|&x| x == crnr).unwrap_or(lp.len());
                lp.rotate_left(pos);
                lp.push(crnr);
            }
        }

        let (th_bgn, th_end) = (self.thalfs[th].bgn, self.thalfs[th].end);
        let (tw_bgn, tw_end) = (self.thalfs[th_twn].bgn, self.thalfs[th_twn].end);
        if th_bgn {
            let t = if self.thalfs[th2].nid_fr(self) == self.thalfs[th].nid_fr(self) { th2 } else { th2_twn };
            self.thalfs[t].bgn = true;
        } else if tw_bgn {
            let t = if self.thalfs[th2].nid_fr(self) == self.thalfs[th_twn].nid_fr(self) { th2 } else { th2_twn };
            self.thalfs[t].bgn = true;
        } else if th_end {
            let t = if self.thalfs[th2].nid_to(self) == self.thalfs[th].nid_to(self) { th2 } else { th2_twn };
            self.thalfs[t].end = true;
        } else if tw_end {
            let t = if self.thalfs[th2].nid_to(self) == self.thalfs[th_twn].nid_to(self) { th2 } else { th2_twn };
            self.thalfs[t].end = true;
        }
    }

    pub fn collapse_tquad_execute(&mut self, aux: &mut Tqaux) {
        let thid_l = aux.thid_l as usize;
        let thid_r = aux.thid_r as usize;

        // the band in its chart: a piece between two points on opposite sides is the shortest line between them
        let tq_thids: Vec<i32> = self.tquads[aux.tqid as usize].data.iter().map(|d| d.thid).collect();
        let chart = unfold(self, &tq_thids);
        let index_of = |tm: &Emesh<'_>, nid: i32, top: bool| -> i32 { // the position of a node of the top / bottom side
            let data = &tm.tquads[aux.tqid as usize].data;
            for k in 0..data.len() {
                if (if top { &aux.thids_t } else { &aux.thids_b }).contains(&data[k].thid)
                    && tm.thalfs[data[k].thid as usize].tedge(tm).nids.contains(&nid) { return chart.at(k as i32, nid); }
            }
            metriko_fail!("node {} is not on the {} side of tquad {}", nid, if top { "top" } else { "bottom" }, aux.tqid);
        };

        let mut pieces: Vec<Piece> = vec![];
        let mut i = 0usize;
        while i + 1 < aux.pts.len() {
            let Tqpoint { nid: n0, val: v0, top: s0, .. } = aux.pts[i];
            let Tqpoint { nid: n1, val: v1, top: s1, .. } = aux.pts[i + 1];
            i += 1;
            if s0 == s1 {
                let run = self.exe_lateral_run(aux, n0, n1);
                metriko_check!(!run.is_empty(), "no lateral thalfs between nodes {} and {}", n0, n1);
                for th in run { pieces.push(Piece { thid: th, t_fr: s0, t_to: s1 }); }
                continue;
            }

            let path = retrace(self, &chart, index_of(self, n0, s0), index_of(self, n1, s1));
            metriko_check!(!path.is_empty(), "no path through tquad {} between nodes {} and {}", aux.tqid, n0, n1);

            let nids  = self.add_new_path(&path, n0, n1);
            let teid  = self.tedges.len() as i32;
            let thid0 = self.thalfs.len() as i32;
            let thid1 = self.thalfs.len() as i32 + 1;
            let x     = (v1 - v0).abs() as f64;
            let r     = self.path_length(&nids);
            self.tedges.push(Eedge { id: teid, nids });
            self.thalfs.push(Ehalf { id: thid0, twid: thid1, teid, tqid: -1, cano: true,  bgn: false, end: false, x, r });
            self.thalfs.push(Ehalf { id: thid1, twid: thid0, teid, tqid: -1, cano: false, bgn: false, end: false, x, r });
            pieces.push(Piece { thid: thid0, t_fr: s0, t_to: s1 });
            pieces.push(Piece { thid: thid1, t_fr: s1, t_to: s0 });
        }

        let Tqpoint { nid: nid_bgn, top: top_bgn, .. } = *aux.pts.first().unwrap();
        let Tqpoint { nid: nid_end, top: top_end, .. } = *aux.pts.last().unwrap();
        let run_bgn = self.exe_take_run(&mut pieces, nid_bgn, !top_bgn, top_bgn);
        let run_end = self.exe_take_run(&mut pieces, nid_end, !top_end, !top_end);
        let mut runs: Vec<Vec<i32>> = vec![];
        let is_head = |tm: &Emesh<'_>, pieces: &[Piece], q: &Piece| !pieces.iter().any(|p| tm.thalfs[p.thid as usize].nid_to(tm) == tm.thalfs[q.thid as usize].nid_fr(tm));
        let is_tail = |tm: &Emesh<'_>, pieces: &[Piece], q: &Piece| !pieces.iter().any(|p| tm.thalfs[p.thid as usize].nid_fr(tm) == tm.thalfs[q.thid as usize].nid_to(tm));
        while !pieces.is_empty() {
            if let Some(it) = pieces.iter().position(|q| is_head(self, &pieces, q)) {
                let (n, s) = (self.thalfs[pieces[it].thid as usize].nid_fr(self), pieces[it].t_fr);
                let r = self.exe_take_run(&mut pieces, n, s, false);
                runs.push(r);
                continue;
            }
            if let Some(it) = pieces.iter().position(|q| is_tail(self, &pieces, q)) {
                let (n, s) = (self.thalfs[pieces[it].thid as usize].nid_to(self), pieces[it].t_to);
                let r = self.exe_take_run(&mut pieces, n, s, true);
                runs.push(r);
                continue;
            }
            metriko_fail!("pieces of the merged line on one side only");
        }

        let on_line = |nid: i32| aux.pts.iter().any(|p| p.nid == nid);
        let keep_of = |tm: &Emesh<'_>, th: usize| -> bool {
            let t = tm.thalfs[th];
            if tm.n_adj_tqs(t.id) == 4 && !on_line(t.nid_fr(tm)) { return true; }
            if tm.n_adj_tqs(t.twid) == 4 && !on_line(t.nid_to(tm)) { return true; }
            let mut across = t.twin(tm).tquad(tm).clone();
            let (l_tw, r_tw) = (tm.thalfs[thid_l].twid, tm.thalfs[thid_r].twid);
            across.data.retain(|d| !(d.thid == l_tw || d.thid == r_tw));
            !across.is_valid()
        };

        let keep_l = keep_of(self, thid_l);
        let keep_r = keep_of(self, thid_r);

        if !run_bgn.is_empty() && run_end.is_empty() {
            metriko_check!(top_bgn == top_end, "the merged line ends on different sides");
            metriko_check!(run_bgn.len() == aux.pts.len() - 1, "the merged line is incomplete");
            let ahd_l = self.thalfs[thid_l].nid_fr(self) == nid_bgn;
            let ahd_r = self.thalfs[thid_r].nid_fr(self) == nid_end;
            self.exe_extend(thid_l as i32, ahd_l, keep_l);
            self.exe_extend(thid_r as i32, ahd_r, keep_r);
            let n0 = if ahd_l { self.thalfs[thid_l].nid_to(self) } else { self.thalfs[thid_l].nid_fr(self) };
            let n1 = if ahd_r { self.thalfs[thid_r].nid_to(self) } else { self.thalfs[thid_r].nid_fr(self) };
            let olds = self.exe_outer(aux, n0, n1);
            self.exe_replace(&olds, &run_bgn);
        } else {
            let end_run = |tm: &mut Emesh<'a>, th: usize, nid_term: i32, keep: bool, run: &[i32]| {
                let ahd = tm.thalfs[th].nid_fr(tm) == nid_term;
                tm.exe_extend(th as i32, ahd, keep);
                let n0 = if ahd { tm.thalfs[th].nid_to(tm) } else { tm.thalfs[th].nid_fr(tm) };
                let na = tm.thalfs[*run.last().unwrap() as usize].nid_to(tm);
                let nb = tm.thalfs[*run.first().unwrap() as usize].nid_fr(tm);
                let n1 = if na == tm.thalfs[th].nid_fr(tm) || na == tm.thalfs[th].nid_to(tm) { nb } else { na };
                let olds = tm.exe_outer(aux, n0, n1);
                tm.exe_replace(&olds, run);
            };
            end_run(self, thid_l, nid_bgn, keep_l, &run_bgn);
            end_run(self, thid_r, nid_end, keep_r, &run_end);
            for run in &runs {
                metriko_check!(run.len() >= 2, "a middle run of the merged line has a single piece");
                let n0 = self.thalfs[*run.first().unwrap() as usize].nid_fr(self);
                let n1 = self.thalfs[*run.last().unwrap() as usize].nid_to(self);
                let olds = self.exe_outer(aux, n0, n1);
                self.exe_replace(&olds, run);
            }
        }

        let tq = aux.tqid as usize;
        let id = self.tquads[tq].id;
        let data: Vec<Edata> = self.tquads[tq].data.clone();
        for d in &data {
            let th0 = d.thid as usize;
            if self.thalfs[th0].tqid != id { continue; }
            let twid = self.thalfs[th0].twid as usize;
            self.thalfs[twid] = Ehalf::default();
            let teid = self.thalfs[th0].teid as usize;
            self.tedges[teid] = Eedge::default();
            self.thalfs[th0] = Ehalf::default();
        }
        self.tquads[tq].data.clear();
        self.tquads[tq].id = -1;
    }
}
