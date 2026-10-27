//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tmesh/emesh_collapse_path.h (see rust/PORTED_FROM).

use crate::common::typedef::*;
use crate::common::r#loop::shortest_path_in_loop;
use crate::hmesh::hmloc::*;
use crate::hmesh::hmflat::*;
use crate::tmesh::emesh::*;

pub const UNFOLD_HAIR: f64 = 1e-4;

#[derive(Clone, Debug, Default)]
pub struct Unfolded<'a> {
    pub nids: Vec<i32>,
    pub heads: Vec<i32>, // the k-th thalf of the boundary starts at nids[heads[k]]
    pub faces: Vec<HmFlat<'a>>,
}

impl<'a> Unfolded<'a> {
    pub fn size(&self) -> i32 { self.nids.len() as i32 }
    pub fn next(&self, i: i32) -> i32 { (i + 1) % self.size() }
    pub fn prev(&self, i: i32) -> i32 { (i + self.size() - 1) % self.size() }
    pub fn at(&self, k: i32, nid: i32) -> i32 {
        let mut i = self.heads[k as usize];
        loop {
            if self.nids[i as usize] == nid { return i; }
            i = self.next(i);
        }
    }
    pub fn uv(&self, tm: &Emesh<'_>, i: i32) -> Complex { flat_uv(&self.faces[i as usize], &tm.tnodes[self.nids[i as usize] as usize]) }
}

pub fn unfold<'a>(tm: &Emesh<'a>, thids: &[i32]) -> Unfolded<'a> {
    let hm = tm.hm;
    let mut c = Unfolded::default();

    for &thid in thids {
        let th = &tm.thalfs[thid as usize];
        let mut ns = th.tedge(tm).nids.clone();
        if !th.cano { ns.reverse(); }
        c.heads.push(c.nids.len() as i32);
        c.nids.extend_from_slice(&ns[..ns.len() - 1]);
    }

    let n = c.size();

    for k in 0..n {
        let a = &tm.tnodes[c.nids[k as usize] as usize];
        let b = &tm.tnodes[c.nids[c.next(k) as usize] as usize];
        let f = try_get_face(hm, a, b);
        metriko_check!(f.is_some(), "boundary segment {} -> {} lies in no face", loc_str(a), loc_str(b));
        c.faces.push(flat_own(f.unwrap()));
    }

    for k in 1..n as usize {
        let q = c.faces[k - 1];
        let p = c.faces[k];
        let l = tm.tnodes[c.nids[k] as usize];
        c.faces[k] = if p.f == q.f {
            q
        } else if let HmLoc::OnV { id } = l {
            *flat_around(&q, hm.vert(id), p.f, false).last().unwrap()
        } else {
            flat_across(&q, try_get_edge(hm, &l, &l).expect("bad optional access"))
        };
    }

    let gap = (flat_uv(&c.faces[n as usize - 1], &tm.tnodes[c.nids[0] as usize]) - c.uv(tm, 0)).norm();
    metriko_check!(gap < 1e-6, "the region does not unfold into one plane (closure error {})", gap);
    c
}

pub fn retrace(tm: &Emesh<'_>, c: &Unfolded<'_>, a: i32, b: i32) -> Vec<HmLoc> {
    let mut lp: Vec<Complex> = vec![];
    for i in 0..c.size() { lp.push(c.uv(tm, i)); }
    let pts = shortest_path_in_loop(&lp, a, b, UNFOLD_HAIR);
    if pts.is_empty() { return vec![]; }

    let mut path = vec![tm.tnodes[c.nids[a as usize] as usize]];
    let mut face = c.faces[a as usize];
    let mut s = 0;
    while s + 1 < pts.len() {
        let stop = if s == 0 { c.faces[c.prev(a) as usize].f } else { face.f };
        let from = *path.last().unwrap();
        for l in flat_trace(&mut face, &from, pts[s + 1], stop) { path.push(l); }
        s += 1;
    }
    *path.last_mut().unwrap() = tm.tnodes[c.nids[b as usize] as usize];
    path
}
