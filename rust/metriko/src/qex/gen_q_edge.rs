//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/qex/gen_q_edge.h (see rust/PORTED_FROM).

use std::collections::BTreeMap;
use crate::common::typedef::*;
use crate::common::predicates::orientation;
use crate::hmesh::hmesh::*;
use crate::qex::common::*;

// sign of orientation with the global tolerance
pub fn orient_sign(a: Complex, b: Complex, c: Complex) -> i32 {
    let o = orientation(a, b, c);
    if o > EPS { 1 } else if o < -EPS { -1 } else { 0 }
}

// closed point-in-triangle test (points on an edge or a vertex count as inside)
pub fn is_inside_face_closed(f: Face<'_>, p: Complex) -> bool {
    let cfn = &f.m.cfn;
    let (a, b, c) = (cfn[(f.id * 3) as usize], cfn[(f.id * 3 + 1) as usize], cfn[(f.id * 3 + 2) as usize]);
    let s = orient_sign(a, b, c);
    if s == 0 { return false; }
    orient_sign(a, b, p) * s >= 0 &&
    orient_sign(b, c, p) * s >= 0 &&
    orient_sign(c, a, p) * s >= 0
}

// does the closed segment [p,q] meet the half-open segment (a,b]?
pub fn meets_half_open(a: Complex, b: Complex, p: Complex, q: Complex) -> bool {
    let op = orient_sign(a, b, p);
    let oq = orient_sign(a, b, q);
    if op == 0 && oq == 0 { // collinear: overlap along the line
        let l2 = (b - a).norm_sqr();
        let t = |x: Complex| ((x - a) * (b - a).conj()).re / l2;
        let (tp, tq) = (t(p), t(q));
        let t0 = if tq < tp { tq } else { tp }; // std::min
        let t1 = if tp < tq { tq } else { tp }; // std::max
        return t1 > EPS && t0 <= 1. + EPS;
    }
    if op * oq > 0 { return false; } // p and q strictly on one side of the ray line
    let oa = orient_sign(p, q, a);
    let ob = orient_sign(p, q, b);
    if oa * ob > 0 { return false; }         // a and b strictly on one side of the edge line
    if oa == 0 && ob != 0 { return false; }  // the only contact is the excluded point a
    true
}

// p's outgoing direction is the reverse of the ray (a, d) traced in face f, compared in 3d so that
// ports living in a neighbouring chart can be matched without accumulating transitions
pub fn is_reverse_port(hm: &Hmesh, a: Complex, d: Complex, f: Face<'_>, p: &Qport) -> bool {
    let pb1 = f.uv2pos(a);
    let pb2 = f.uv2pos(a + d);
    let pa1 = hm.face(p.fid).uv2pos(p.uv);
    let pa2 = hm.face(p.fid).uv2pos(p.uv + p.dir);
    let da = normalized(pa2 - pa1);
    let db = normalized(pb2 - pb1);
    (1. + da.dot(db)).abs() < 1e-7
}

// PICK_NEXT_EDGE of Ebke et al. 2013 (Alg. 5): the edge of f other than the one entered through that
// meets (a,b]. if two edges meet it (the ray passes through a vertex or runs along an edge), take the one
// with fewer endpoints on the ray line, which steps around the vertex until the ray leaves properly
pub fn pick_next_half<'a>(a: Complex, b: Complex, f: Face<'a>, skip_hid: i32) -> Option<Half<'a>> {
    let mut pick = None;
    let mut best = 3;
    for h in f.adj_halfs() {
        if h.id == skip_hid { continue; }
        let p = h.cr_t().uv();
        let q = h.cr_h().uv();
        if !meets_half_open(a, b, p, q) { continue; }
        let on_line = (orient_sign(a, b, p) == 0) as i32 + (orient_sign(a, b, q) == 0) as i32;
        if on_line < best { best = on_line; pick = Some(h); }
    }
    pick
}

pub fn generate_q_edge(hm: &Hmesh, matching: &[i32], qports: &mut [Qport]) -> Vec<Qedge> {
    let mut qedges: Vec<Qedge> = vec![];
    let (he_r, he_t) = compute_trs_matrix(hm, matching, 4);

    // ports grouped by carrier for the arrival lookup (the C++ maps hold pointers into qports: positions here)
    let mut by_v: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    let mut by_e: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    let mut by_f: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (k, p) in qports.iter().enumerate() {
        if      p.vid >= 0 { by_v.entry(p.vid).or_default().push(k); }
        else if p.eid >= 0 { by_e.entry(p.eid).or_default().push(k); }
        else               { by_f.entry(p.fid).or_default().push(k); }
    }

    for pi in 0..qports.len() {
        if qports[pi].is_connected { continue; }
        let mut a = qports[pi].uv;
        let mut d = qports[pi].dir;
        let mut b = a + d;
        let mut fid = qports[pi].fid;
        let mut e_in = -1;
        if qports[pi].eid >= 0 { // an eqvert port starts on its own edge: never cross back over it
            let h0 = hm.edge(qports[pi].eid).half();
            e_in = (if h0.face().id == fid { h0 } else { h0.twin() }).id;
        }

        for _step in 0..1000 {
            let f = hm.face(fid);

            if is_inside_face_closed(f, b) {
                // the target grid point lies in the closure of f: find its q-vertex (vertex, edge, then face)
                let mut hit: Option<usize> = None;
                for h in f.adj_halfs() {
                    if (h.cr_t().uv() - b).norm() >= EPS || !by_v.contains_key(&h.tail().id) { continue; }
                    for &p in &by_v[&h.tail().id] {
                        if !qports[p].is_connected && p != pi && is_reverse_port(hm, a, d, f, &qports[p]) { hit = Some(p); break; }
                    }
                    if hit.is_some() { break; }
                }
                if hit.is_none() {
                    let pb = f.uv2pos(b);
                    for h in f.adj_halfs() {
                        if orient_sign(h.cr_t().uv(), h.cr_h().uv(), b) != 0 || !by_e.contains_key(&h.edge().id) { continue; }
                        for &p in &by_e[&h.edge().id] {
                            if !qports[p].is_connected && p != pi && (qports[p].pos - pb).length() < 1e-6 && is_reverse_port(hm, a, d, f, &qports[p]) { hit = Some(p); break; }
                        }
                        if hit.is_some() { break; }
                    }
                }
                if hit.is_none() && by_f.contains_key(&fid) {
                    for &p in &by_f[&fid] {
                        if !qports[p].is_connected && p != pi && (qports[p].uv - b).norm() < EPS && (qports[p].dir + d).norm() < EPS { hit = Some(p); break; }
                    }
                }
                if let Some(h) = hit {
                    qports[pi].is_connected = true;
                    qports[h].is_connected = true;
                    qedges.push(Qedge { port1: pi, port2: h });
                }
                break; // reached b: paired, or left dangling
            }

            let Some(nh) = pick_next_half(a, b, f, e_in) else { break }; // numerically inconsistent chart: leave dangling
            metriko_check!(!nh.twin().is_boundary(), "a q-edge reaches the mesh boundary (not handled yet)");
            let r = he_r[nh.id as usize];
            let t = he_t[nh.id as usize];
            a = r * a + t;
            b = r * b + t;
            d = r * d;
            fid = nh.twin().face().id;
            e_in = nh.twin().id;
        }
    }

    for p in qports.iter() {
        if !p.is_connected { println!("[qedge] unpaired port {} (vid {}, eid {}, fid {})", p.idx, p.vid, p.eid, p.fid); }
    }

    qedges
}
