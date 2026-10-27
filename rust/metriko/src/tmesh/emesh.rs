//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tmesh/emesh.h (see rust/PORTED_FROM).

use crate::common::utilities::*;
use crate::common::stdcxx;
use crate::hmesh::hmesh::*;
use crate::hmesh::hmloc::*;
use crate::hmesh::utilities::get_face_uv;
use crate::tmesh::motorcycle::*;

#[derive(Clone, Copy, Debug)]
pub struct Tqpoint {
    pub nid: i32,  // tnode of this point
    pub val: i32,
    pub ord: f64,  // geometric sort key: arc-length (th.r) position along the spine, per-tquad normalized
    pub top: bool,
}

impl Default for Tqpoint {
    fn default() -> Self { Tqpoint { nid: -1, val: -1, ord: -1., top: false } }
}

#[derive(Clone, Debug)]
pub struct Tqchain {
    pub pts: Vec<Tqpoint>,
    pub tqids: Vec<i32>,
    pub bounds: Vec<i32>,
    pub thids_z: Vec<i32>, // zero length thids from left to right
    pub thids_t: Vec<i32>,
    pub thids_b: Vec<i32>,
    pub thid_r: i32,
    pub thid_l: i32,
}

impl Default for Tqchain {
    fn default() -> Self {
        Tqchain { pts: vec![], tqids: vec![], bounds: vec![], thids_z: vec![], thids_t: vec![], thids_b: vec![], thid_r: -1, thid_l: -1 }
    }
}

// a zero-width band about to be collapsed on its own: its two zero sides (l at val 0, r at val = length), its
// two lateral sides, and the points the merged line passes, sorted along it
#[derive(Clone, Debug)]
pub struct Tqaux {
    pub tqid: i32,
    pub thid_l: i32,
    pub thid_r: i32,
    pub thids_t: Vec<i32>,
    pub thids_b: Vec<i32>,
    pub pts: Vec<Tqpoint>,
}

impl Default for Tqaux {
    fn default() -> Self { Tqaux { tqid: -1, thid_l: -1, thid_r: -1, thids_t: vec![], thids_b: vec![], pts: vec![] } }
}

#[derive(Clone, Debug)]
pub struct Eedge {
    pub id: i32,
    pub nids: Vec<i32>,
}

impl Default for Eedge {
    fn default() -> Self { Eedge { id: -1, nids: vec![] } }
}

/// A thalf. The C++ keeps a back pointer `tm` for its navigation methods; here they take the `Emesh`.
#[derive(Clone, Copy, Debug)]
pub struct Ehalf {
    pub id: i32,
    pub twid: i32,
    pub teid: i32,
    pub tqid: i32,
    pub cano: bool,
    pub bgn: bool,
    pub end: bool,
    pub x: f64,
    pub r: f64,
}

impl Default for Ehalf {
    fn default() -> Self { Ehalf { id: -1, twid: -1, teid: -1, tqid: -1, cano: false, bgn: false, end: false, x: -1., r: -1. } }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Edata {
    pub thid: i32,
    pub side: i32,
}

impl Default for Edata {
    fn default() -> Self { Edata { thid: -1, side: -1 } }
}

#[derive(Clone, Debug)]
pub struct Equad {
    pub id: i32,
    pub data: Vec<Edata>,
}

impl Default for Equad {
    fn default() -> Self { Equad { id: -1, data: vec![] } }
}

impl Equad {
    /// `find_of(thid)`: the position of `thid` in `data`, `None` for `end()`
    pub fn find_of(&self, thid: i32) -> Option<usize> { self.data.iter().position(|d| d.thid == thid) }
    pub fn side_of(&self, th: &Ehalf) -> i32 {
        self.data[self.find_of(th.id).expect("side_of: thalf not in tquad")].side
    }

    // every side keeps at least one thalf: a tquad that lost a side is no longer a quad
    pub fn is_valid(&self) -> bool {
        for side in 0..4 { if !self.data.iter().any(|d| d.side == side) { return false; } }
        true
    }

    pub fn thids(&self, side: i32) -> Vec<i32> {
        self.data.iter().filter(|d| d.side == side).map(|d| d.thid).collect()
    }
}

pub struct Emesh<'a> {
    pub hm: &'a Hmesh,
    pub tnodes: Vec<HmLoc>,
    pub tedges: Vec<Eedge>,
    pub thalfs: Vec<Ehalf>,
    pub tquads: Vec<Equad>,
    pub fixed: Vec<bool>,
}

impl Ehalf {
    pub fn nid_fr(&self, tm: &Emesh<'_>) -> i32 {
        let nids = &tm.tedges[self.teid as usize].nids;
        if self.cano { *nids.first().unwrap() } else { *nids.last().unwrap() }
    }
    pub fn nid_to(&self, tm: &Emesh<'_>) -> i32 {
        let nids = &tm.tedges[self.teid as usize].nids;
        if self.cano { *nids.last().unwrap() } else { *nids.first().unwrap() }
    }
    pub fn side(&self, tm: &Emesh<'_>) -> i32 { self.tquad(tm).side_of(self) }
    pub fn twin<'t>(&self, tm: &'t Emesh<'_>) -> &'t Ehalf { &tm.thalfs[self.twid as usize] }
    pub fn tedge<'t>(&self, tm: &'t Emesh<'_>) -> &'t Eedge { &tm.tedges[self.teid as usize] }
    pub fn tquad<'t>(&self, tm: &'t Emesh<'_>) -> &'t Equad { &tm.tquads[self.tqid as usize] }
    pub fn loc_fr<'t>(&self, tm: &'t Emesh<'_>) -> &'t HmLoc { &tm.tnodes[self.nid_fr(tm) as usize] }
    pub fn loc_to<'t>(&self, tm: &'t Emesh<'_>) -> &'t HmLoc { &tm.tnodes[self.nid_to(tm) as usize] }
}

impl<'a> Emesh<'a> {
    /// `explicit Emesh(const Hmesh& hm)`
    pub fn empty(hm: &'a Hmesh) -> Emesh<'a> {
        Emesh { hm, tnodes: vec![], tedges: vec![], thalfs: vec![], tquads: vec![], fixed: vec![] }
    }

    /// `explicit Emesh(const Mgrph& mg)`
    pub fn new(mg: &Mgrph<'a>) -> Emesh<'a> {
        let hm = mg.hm;
        let mut tm = Emesh::empty(hm);
        tm.tnodes.reserve(mg.mnodes.len());
        for mn in &mg.mnodes {
            tm.tnodes.push(match mn.loc {
                HmLoc::OnV { id } => HmLoc::OnV { id },
                HmLoc::OnE { id, r } => HmLoc::OnE { id, r },
                HmLoc::OnP { id, uv } => { let f = hm.face(id); HmLoc::OnF { id, xy: f.to_local(f.uv2pos(uv)) } }
                _ => metriko_fail!("no impl"),
            });
            tm.fixed.push(mn.jt == JunctionType::F);
        }

        // 1: one tedge per run of segments between two junction nodes. a trailing run that ends at no junction
        //    does not make a tedge
        let mut crv: Vec<i32> = vec![];      // teid -> curve it lies on, delimits the tquad sides below
        let mut sg_bgn: Vec<Msgmt> = vec![]; // teid -> its first / last segment, ranks the thalfs at the junctions
        let mut sg_end: Vec<Msgmt> = vec![];
        for mc in &mg.mcurvs {
            let mut bgn_nid = mc.sgmts.first().unwrap().fr_nid;
            let mut bgn_sg = 0usize;
            let mut len = 0.;
            let mut nids = vec![bgn_nid];
            for i in 0..mc.sgmts.len() {
                let sg = &mc.sgmts[i];
                len += (get_face_uv(&mg.mnodes[sg.to_nid as usize].loc, sg.face_id, hm) - get_face_uv(&mg.mnodes[sg.fr_nid as usize].loc, sg.face_id, hm)).norm();
                nids.push(sg.to_nid);
                if mg.mnodes[sg.to_nid as usize].jt == JunctionType::None { continue; }

                let teid = tm.tedges.len() as i32;
                let thid = tm.thalfs.len() as i32;
                let bgn = mg.mnodes[bgn_nid as usize].jt == JunctionType::F;
                let end = mg.mnodes[sg.to_nid as usize].jt == JunctionType::T;
                tm.tedges.push(Eedge { id: teid, nids: std::mem::take(&mut nids) });
                tm.thalfs.push(Ehalf { id: thid,     twid: thid + 1, teid, cano: true,  bgn, end, r: len, ..Default::default() });
                tm.thalfs.push(Ehalf { id: thid + 1, twid: thid,     teid, cano: false, r: len, ..Default::default() });
                crv.push(mc.id);
                sg_bgn.push(mc.sgmts[bgn_sg]);
                sg_end.push(*sg);

                bgn_nid = sg.to_nid;
                bgn_sg = i + 1;
                len = 0.;
                nids = vec![bgn_nid];
            }
        }

        // 2: at every junction, the thalf entering it continues with the next outgoing one clockwise
        let mut nxid = vec![-1i32; tm.thalfs.len()];
        let mut outgoing: Vec<Vec<i32>> = vec![vec![]; mg.mnodes.len()];
        for th in &tm.thalfs { outgoing[th.nid_fr(&tm) as usize].push(th.id); }
        for nid in 0..mg.mnodes.len() {
            let mn = &mg.mnodes[nid];
            let out = &mut outgoing[nid];
            if mn.jt == JunctionType::None || out.is_empty() { continue; }
            let rank = |thid: i32| -> usize {
                let th = &tm.thalfs[thid as usize];
                let sg = if th.cano { &sg_bgn[th.teid as usize] } else { &sg_end[th.teid as usize] };
                mn.adj.iter().position(|ad| ad[0] == sg.curv_id && ad[1] == sg.this_id).unwrap_or(mn.adj.len())
            };
            stdcxx::sort_by(out, |a, b| rank(*a) < rank(*b));
            let n = out.len() as i32;
            for i in 0..n {
                nxid[tm.thalfs[out[i as usize] as usize].twid as usize] = out[((i - 1 + n) % n) as usize];
            }
        }

        // 3: trace every tquad along nxid. a side is a run of thalfs on one curve
        let mut visited = vec![false; tm.thalfs.len()];
        for i in 0..tm.thalfs.len() {
            if visited[i] { continue; }
            let mut tq = Equad { id: tm.tquads.len() as i32, data: vec![] };
            let mut thid = i as i32;
            let mut side = 0;
            loop {
                if visited[thid as usize] { break; }
                tq.data.push(Edata { thid, side });
                visited[thid as usize] = true;
                let next = nxid[thid as usize];
                if crv[tm.thalfs[thid as usize].teid as usize] != crv[tm.thalfs[next as usize].teid as usize] { side = (side + 1) % 4; }
                thid = next;
                if thid == i as i32 { break; }
            }

            // the walk started inside a side: rotate so that every side is one contiguous run, then renumber
            if crv[tm.thalfs[tq.data.first().unwrap().thid as usize].teid as usize] == crv[tm.thalfs[tq.data.last().unwrap().thid as usize].teid as usize] {
                let f = tq.data.first().unwrap().side;
                let n = tq.data.iter().take_while(|d| d.side == f).count();
                tq.data.rotate_left(n);
                let mut s = 0;
                let mut c = crv[tm.thalfs[tq.data.first().unwrap().thid as usize].teid as usize];
                for d in tq.data.iter_mut() {
                    if crv[tm.thalfs[d.thid as usize].teid as usize] != c { s = (s + 1) % 4; }
                    d.side = s;
                    c = crv[tm.thalfs[d.thid as usize].teid as usize];
                }
            }
            for d in &tq.data { tm.thalfs[d.thid as usize].tqid = tq.id; }
            tm.tquads.push(tq);
        }
        tm
    }

    pub fn set_x(&mut self, x: &[f64]) {
        for th in self.thalfs.iter_mut() { if th.id != -1 { th.x = x[th.teid as usize]; } }
    }

    pub fn is_fixed(&self, nid: i32) -> bool { (nid as usize) < self.fixed.len() && self.fixed[nid as usize] }

    pub fn n_adj_tqs(&self, thid: i32) -> i32 {
        let mut n = 0;
        let mut thid_ = thid;
        loop {
            n += 1;
            thid_ = self.step_next(self.thalfs[thid_ as usize].twid);
            if !(thid_ != thid && n as usize <= self.thalfs.len()) { break; }
        }
        n
    }

    pub fn step_next(&self, thid: i32) -> i32 {
        let d = &self.thalfs[thid as usize].tquad(self).data;
        let it = d.iter().position(|e| e.thid == thid);
        metriko_check!(it.is_some(), "step next failed");
        d[circular_next(d.len(), it.unwrap())].thid
    }

    pub fn step_prev(&self, thid: i32) -> i32 {
        let d = &self.thalfs[thid as usize].tquad(self).data;
        let it = d.iter().position(|e| e.thid == thid);
        metriko_check!(it.is_some(), "step prev failed");
        d[circular_prev(d.len(), it.unwrap())].thid
    }

    pub fn live_tedges(&self) -> impl Iterator<Item = &Eedge> + '_ { self.tedges.iter().filter(|te| te.id != -1) }
    pub fn live_tquads(&self) -> impl Iterator<Item = &Equad> + '_ { self.tquads.iter().filter(|tq| tq.id != -1) }
    pub fn live_thalfs(&self) -> impl Iterator<Item = &Ehalf> + '_ { self.thalfs.iter().filter(|th| th.id != -1) }
    pub fn cano_thalfs(&self) -> impl Iterator<Item = &Ehalf> + '_ { self.thalfs.iter().filter(|th| th.id != -1 && th.cano) }

    pub fn add_new_path(&mut self, path: &[HmLoc], nid0: i32, nid1: i32) -> Vec<i32> {
        let mut nids = vec![];
        for i in 0..path.len() {
            if i == 0 { nids.push(nid0); }
            else if i + 1 == path.len() { nids.push(nid1); }
            else { self.tnodes.push(path[i]); nids.push(self.tnodes.len() as i32 - 1); }
        }
        nids
    }

    pub fn path_length(&self, nids: &[i32]) -> f64 {
        let mut r = 0.;
        let mut k = 0;
        while k + 1 < nids.len() {
            r += (get_ptloc_pos(self.hm, &self.tnodes[nids[k + 1] as usize]) - get_ptloc_pos(self.hm, &self.tnodes[nids[k] as usize])).length();
            k += 1;
        }
        r
    }
}

impl Eedge {
    pub fn insert_locs(&mut self, locs: &[i32]) {
        let f = *locs.first().unwrap();
        let b = *locs.last().unwrap();
        if *self.nids.first().unwrap() == b {
            // prepend [f..b-1]
            let s: Vec<i32> = locs[..locs.len() - 1].to_vec();
            self.nids.splice(0..0, s);
        } else if *self.nids.last().unwrap() == f {
            // append  [f+1..b]
            self.nids.extend_from_slice(&locs[1..]);
        } else if *self.nids.first().unwrap() == f {
            // prepend reverse([f+1..b])
            let mut s: Vec<i32> = locs[1..].to_vec();
            s.reverse();
            self.nids.splice(0..0, s);
        } else if *self.nids.last().unwrap() == b {
            // append  reverse([f..b-1])
            let mut s: Vec<i32> = locs[..locs.len() - 1].to_vec();
            s.reverse();
            self.nids.extend_from_slice(&s);
        } else {
            metriko_fail!("merge_into: the chains share no end node");
        }
    }
}
