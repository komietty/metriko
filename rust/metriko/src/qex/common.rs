//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/qex/common.h (see rust/PORTED_FROM).

use crate::common::typedef::*;
use crate::common::utilities::get_quater_rot;
use crate::hmesh::hmesh::Hmesh;

#[derive(Clone, Copy, Debug)]
pub struct Qvert {
    pub uv: Complex,
    pub pos: Row3d,
    pub sid: i32, // simplex id (either vid, eid or fid)
}

impl Qvert {
    pub fn new(uv: Complex, pos: Row3d, sid: i32) -> Self { Qvert { uv, pos, sid } }
}

#[derive(Clone, Copy, Debug)]
pub struct Qport {
    pub idx: i32,
    pub vid: i32,
    pub eid: i32,
    pub fid: i32,
    pub uv: Complex,
    pub dir: Complex,
    pub pos: Row3d,
    pub is_connected: bool,
    pub next_id: i32,
    pub prev_id: i32,
}

impl Qport {
    pub fn new(idx: i32, vid: i32, eid: i32, fid: i32, uv: Complex, dir: Complex, pos: Row3d) -> Self {
        Qport { idx, vid, eid, fid, uv, dir, pos, is_connected: false, next_id: -1, prev_id: -1 }
    }
}

/// The C++ `Qedge` holds references to two elements of the `q_ports` vector; here their positions in it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Qedge {
    pub port1: usize,
    pub port2: usize,
}

/// The C++ `Qhalf` holds a reference to an element of the `q_edges` vector; here its position `qe` in it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Qhalf {
    pub qe: usize,
    pub idx: i32,
    pub cano: bool, // cannonical flag
}

impl Qhalf {
    /// the position in `q_ports` of `port1()`
    pub fn port1(&self, qes: &[Qedge]) -> usize { if self.cano { qes[self.qe].port1 } else { qes[self.qe].port2 } }
    /// the position in `q_ports` of `port2()`
    pub fn port2(&self, qes: &[Qedge]) -> usize { if self.cano { qes[self.qe].port2 } else { qes[self.qe].port1 } }
}

#[derive(Clone, Debug, Default)]
pub struct Qface {
    pub qhalfs: Vec<Qhalf>,
}

/// `compute_trs_matrix(hm, matching, rosyN, heR, heT)`: returns `(heR, heT)`
pub fn compute_trs_matrix(hm: &Hmesh, matching: &[i32], rosy_n: i32) -> (Vec<Complex>, Vec<Complex>) {
    let mut he_r = vec![cx(0., 0.); hm.n_h as usize];
    let mut he_t = vec![cx(0., 0.); hm.n_h as usize];
    for h in hm.halfs() {
        let c1 = h.next().crnr();
        let c2 = h.twin().prev().crnr();
        let uv1 = c1.uv();
        let uv2 = c2.uv();
        let mut m = (if h.is_canonical() { 1 } else { -1 }) * matching[h.edge().id as usize];
        m = if m < 0 { (rosy_n + m % rosy_n) % rosy_n } else { m % rosy_n };
        let r = get_quater_rot(m);
        let t = uv2 - r * uv1;
        he_r[h.id as usize] = r;
        he_t[h.id as usize] = cx(t.re.round(), t.im.round());
    }
    (he_r, he_t)
}

/// `nearby_grid(complex a)`
pub fn nearby_grid(a: Complex) -> Complex {
    cx(a.re.round(), a.im.round())
}

/// `nearby_grid(double x, double dir)`
pub fn nearby_grid_x(mut x: f64, dir: f64) -> f64 {
    if dir == 0. { return x; }
    let f = (x % 1.).abs();
    if f < EPS || 1. - f < EPS { x += dir * EPS; }
    x = if dir > 0. { x.ceil() } else { x.floor() };
    x
}

/// `nearby_grid(complex uv, complex dir)`
pub fn nearby_grid_dir(uv: Complex, dir: Complex) -> Complex {
    cx(nearby_grid_x(uv.re, dir.re), nearby_grid_x(uv.im, dir.im))
}
