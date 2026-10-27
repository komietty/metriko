//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/common/loop.h (see rust/PORTED_FROM).

// loop: a closed simple ccw polyline

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use crate::common::typedef::*;
use crate::common::utilities::*;
use crate::common::predicates::*;

pub fn loop_next(lp: &[Complex], i: usize) -> usize { (i + 1) % lp.len() }
pub fn loop_prev(lp: &[Complex], i: usize) -> usize { (i + lp.len() - 1) % lp.len() }

pub fn loop_inward(lp: &[Complex], i: usize) -> Complex {
    let o = lp[loop_next(lp, i)] - lp[i];
    let e = lp[loop_prev(lp, i)] - lp[i];
    let mut d = -(o / o.norm() + e / e.norm());
    if d.norm() < 1e-6 { d = o / o.norm() * cx(0., 1.); }
    d / d.norm()
}

pub fn loop_turn(lp: &[Complex], i: usize) -> f64 {
    let a = lp[i] - lp[loop_prev(lp, i)];
    let b = lp[loop_next(lp, i)] - lp[i];
    cross(a, b)
}

pub fn loop_into(lp: &[Complex], i: usize, d: Complex) -> bool {
    let o = lp[loop_next(lp, i)] - lp[i];
    let e = lp[loop_prev(lp, i)] - lp[i];
    if cross(o, e) > 0. { cross(o, d) > EPS && cross(d, e) > EPS } else { cross(o, d) > EPS || cross(d, e) > EPS }
}

/// the segment p -> q lies inside the loop. an end sitting on corner ia / ib (-1: none) must leave it inwards
pub fn loop_sees(lp: &[Complex], p: Complex, q: Complex, ia: i32, ib: i32) -> bool {
    if ia >= 0 && !loop_into(lp, ia as usize, q - p) { return false; }
    if ib >= 0 && !loop_into(lp, ib as usize, p - q) { return false; }
    for k in 0..lp.len() {
        if find_strict_intersection(p, q, lp[k], lp[loop_next(lp, k)], EPS) { return false; }
        if k as i32 != ia && k as i32 != ib && is_inside_segment(p, q, lp[k]) { return false; }
    }
    true
}

/// `(double, int)` ordered as `std::pair<double, int>`: the C++ min-heaps pop the smallest pair; equal pairs are
/// interchangeable, so a `BinaryHeap<Reverse<_>>` pops in the same order as libc++'s `priority_queue`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistId(pub f64, pub i32);
impl Eq for DistId {}
impl PartialOrd for DistId {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> { Some(self.cmp(o)) }
}
impl Ord for DistId {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        if self.0 < o.0 { return std::cmp::Ordering::Less; }
        if o.0 < self.0 { return std::cmp::Ordering::Greater; }
        self.1.cmp(&o.1)
    }
}

/// the shortest path from corner a to corner b inside the loop (empty when there is no way through). it bends
/// only at corners that are not convex, and there it is pushed inside by `hair` so that it touches no boundary
pub fn shortest_path_in_loop(lp: &[Complex], a: i32, b: i32, hair: f64) -> Vec<Complex> {
    let mut idx = vec![a, b];
    let mut pts = vec![lp[a as usize], lp[b as usize]];
    if loop_sees(lp, pts[0], pts[1], a, b) { return pts; }
    for i in 0..lp.len() {
        if i as i32 != a && i as i32 != b && loop_turn(lp, i) < 1e-6 {
            idx.push(i as i32);
            pts.push(lp[i] + loop_inward(lp, i) * hair);
        }
    }

    let m = idx.len();
    let mut dist = vec![f64::INFINITY; m];
    let mut from = vec![-1i32; m];
    let mut pq: BinaryHeap<Reverse<DistId>> = BinaryHeap::new();
    dist[0] = 0.;
    pq.push(Reverse(DistId(0., 0)));
    while let Some(Reverse(DistId(d, u))) = pq.pop() {
        let u = u as usize;
        if d > dist[u] { continue; }
        if u == 1 { break; }
        for v in 0..m {
            let w = d + (pts[v] - pts[u]).norm();
            if w >= dist[v] || !loop_sees(lp, pts[u], pts[v], if u < 2 { idx[u] } else { -1 }, if v < 2 { idx[v] } else { -1 }) { continue; }
            dist[v] = w;
            from[v] = u as i32;
            pq.push(Reverse(DistId(w, v as i32)));
        }
    }
    if from[1] == -1 { return vec![]; }

    let mut path = vec![];
    let mut v = 1i32;
    while v != -1 { path.push(pts[v as usize]); v = from[v as usize]; }
    path.reverse();
    path
}
