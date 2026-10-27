//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/hmesh/utilities.h (see rust/PORTED_FROM).

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use faer::sparse::Triplet;
use crate::common::typedef::*;
use crate::common::r#loop::DistId;
use crate::hmesh::hmesh::*;
use crate::hmesh::hmloc::*;
use crate::solver::matrix_ops::{SprsD, from_triplets};

pub fn try_get_crnr(hm: &Hmesh, vid: i32, fid: i32) -> Option<Crnr<'_>> {
    let f = hm.face(fid);
    let v = hm.vert(vid);
    for h in f.adj_halfs() {
        if h.crnr().vert() == v { return Some(h.crnr()); }
    }
    None
}

pub fn try_get_half(hm: &Hmesh, eid: i32, fid: i32) -> Option<Half<'_>> {
    let f = hm.face(fid);
    let e = hm.edge(eid);
    for h in f.adj_halfs() {
        if h.edge() == e { return Some(h); }
    }
    None
}

pub fn get_face_uv(loc: &HmLoc, fid: i32, hm: &Hmesh) -> Complex {
    match *loc {
        HmLoc::OnP { uv, .. } => uv,
        HmLoc::OnV { id } => try_get_crnr(hm, id, fid).expect("bad optional access").uv(),
        HmLoc::OnE { id, r } => {
            let h = try_get_half(hm, id, fid).expect("bad optional access");
            h.uv(if h.is_canonical() { r } else { 1. - r })
        }
        _ => metriko_fail!("no impl"),
    }
}

// per face weight of the poisson energy. around a cone the least squares integration admits a mode whose gradient
// grows like r^{-3/4} (z^{1/4} at index -1, conj(z)^{1/4} at +1): its dirichlet energy is finite, and once the mesh
// resolves it the cone is wound by -2pi. the weight (R / d)^beta, d the distance of the face to the nearest
// singular, makes the energy of that mode infinite for beta >= 1/2 and leaves the right modes (z^{5/4}, z^{3/4})
// finite for beta < 3/2. beta = 1; R = 8 grid cells is measured (4 still leaves wrong cones on bimba, rolling_stage).
// returns the diagonal matrix with the weight of each face repeated over its k rows (the 2N field components of the
// integration). gridscale: the grid spacing as a fraction of the bounding box diagonal
pub fn compute_poisson_weight_matrix(hm: &Hmesh, singular: &[i32], gridscale: f64, k: i32) -> SprsD {
    let beta = 1.;
    let cells = 8.;
    let mut mx = hm.pos[0];
    let mut mn = hm.pos[0];
    for p in &hm.pos[1..] {
        if p.x > mx.x { mx.x = p.x; }
        if p.y > mx.y { mx.y = p.y; }
        if p.z > mx.z { mx.z = p.z; }
        if p.x < mn.x { mn.x = p.x; }
        if p.y < mn.y { mn.y = p.y; }
        if p.z < mn.z { mn.z = p.z; }
    }
    let r_ = (mx - mn).length() * gridscale * cells;

    // edge path distance to the nearest inner singular, up to R
    let mut d = vec![f64::INFINITY; hm.n_v as usize];
    let mut pq: BinaryHeap<Reverse<DistId>> = BinaryHeap::new();
    for v in hm.verts() {
        if !v.is_boundary() && singular[v.id as usize] != 0 { d[v.id as usize] = 0.; pq.push(Reverse(DistId(0., v.id))); }
    }
    while let Some(Reverse(DistId(dv, v))) = pq.pop() {
        if dv > d[v as usize] || dv > r_ { continue; }
        for h in hm.vert(v).adj_halfs() {
            let w = h.head().id;
            let nd = dv + h.len();
            if nd < d[w as usize] { d[w as usize] = nd; pq.push(Reverse(DistId(nd, w))); }
        }
    }

    let mut w = vec![0.; hm.n_f as usize];
    for f in hm.faces() {
        let mut df = 0.;
        let mut el = 0.;
        for h in f.adj_halfs() { df += d[h.tail().id as usize].min(r_); el += h.len(); }
        let a = df / 3.;
        let b = el / 9.;
        // std::max(a, b): b when a < b, else a
        w[f.id as usize] = (r_ / if a < b { b } else { a }).powf(beta);
    }
    // w.replicate(1, k).transpose().reshaped().asDiagonal(): the weight of face f on rows f*k .. f*k+k-1
    let n = hm.n_f as usize * k as usize;
    let mut t = Vec::with_capacity(n);
    for f in 0..hm.n_f as usize {
        for j in 0..k as usize { t.push(Triplet::new(f * k as usize + j, f * k as usize + j, w[f])); }
    }
    from_triplets(n, n, &t)
}
