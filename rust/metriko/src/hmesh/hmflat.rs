//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/hmesh/hmflat.h (see rust/PORTED_FROM).

use crate::common::typedef::*;
use crate::common::utilities::*;
use crate::common::stdcomplex::{std_div, std_polar};
use crate::common::predicates::side_distance;
use crate::hmesh::hmesh::*;
use crate::hmesh::hmloc::*;
use crate::hmesh::utilities::try_get_crnr;

#[derive(Clone, Copy, Debug)]
pub struct HmFlat<'a> {
    pub f: Face<'a>,
    pub uv: [Complex; 3],
}

pub fn flat_own(f: Face<'_>) -> HmFlat<'_> {
    let [h0, h1, h2] = f.halfs();
    HmFlat { f, uv: [h0.cr_t().uv(), h1.cr_t().uv(), h2.cr_t().uv()] }
}

pub fn flat_uv_vert(c: &HmFlat<'_>, v: Vert<'_>) -> Complex {
    let vs = c.f.verts();
    // rg::find(vs, v) - vs.begin()
    let i = vs.iter().position(|&x| x == v).unwrap_or(3);
    c.uv[i]
}

pub fn flat_uv(c: &HmFlat<'_>, loc: &HmLoc) -> Complex {
    let f = c.f;
    match *loc {
        HmLoc::OnV { id } => flat_uv_vert(c, f.m.vert(id)),
        HmLoc::OnE { id, r } => { let e = f.m.edge(id); lerp(flat_uv_vert(c, e.vert0()), flat_uv_vert(c, e.vert1()), r) }
        HmLoc::OnF { xy, .. } => {
            let [a, b, d] = f.verts();
            let w = calc_coefficient(f.to_local(a.pos()), f.to_local(b.pos()), f.to_local(d.pos()), xy);
            c.uv[0] + (c.uv[1] - c.uv[0]) * w.re + (c.uv[2] - c.uv[0]) * w.im
        }
        _ => metriko_fail!("{} has no uv", loc_str(loc)),
    }
}

pub fn flat_across<'a>(c: &HmFlat<'a>, e: Edge<'a>) -> HmFlat<'a> {
    let mut g = flat_own(if e.face0() == c.f { e.face1() } else { e.face0() });
    let d = std_div(flat_uv_vert(c, e.vert1()) - flat_uv_vert(c, e.vert0()), flat_uv_vert(&g, e.vert1()) - flat_uv_vert(&g, e.vert0()));
    let rot = std_polar(1., PI / 2. * (d.arg() / (PI / 2.)).round());
    let off = flat_uv_vert(c, e.vert0()) - rot * flat_uv_vert(&g, e.vert0());
    for uv in g.uv.iter_mut() { *uv = rot * *uv + off; }
    g
}

/// `flat_around(c, v, stop, ccw = true)`
pub fn flat_around<'a>(c: &HmFlat<'a>, v: Vert<'a>, stop: Face<'a>, ccw: bool) -> Vec<HmFlat<'a>> {
    let mut res = vec![*c];
    let h = try_get_crnr(c.f.m, v.id, c.f.id).expect("bad optional access").half().prev(); // the half of f leaving v
    let mut p = h;
    for g in v.adj_halfs_from(h, ccw) {
        if g.face() == c.f { continue; }
        let next = flat_across(res.last().unwrap(), (if ccw { g } else { p }).edge());
        res.push(next);
        p = g;
        if g.face() == stop { break; }
    }
    res
}

pub fn flat_loc(c: &HmFlat<'_>, q: Complex) -> HmLoc {
    let hs = c.f.halfs();
    for k in 0..3 {
        let a = c.uv[k];
        let b = c.uv[(k + 1) % 3];
        if side_distance(a, b, q).abs() > EPS { continue; }
        let r = (q - a).norm() / (b - a).norm();
        if r < EPS { return HmLoc::OnV { id: hs[k].tail().id }; }
        if r > 1. - EPS { return HmLoc::OnV { id: hs[k].head().id }; }
        return HmLoc::OnE { id: hs[k].edge().id, r: if hs[k].is_canonical() { r } else { 1. - r } };
    }
    let [a, b, d] = c.f.verts();
    HmLoc::OnF { id: c.f.id, xy: c.f.to_local(conversion_2d_3d(c.uv[0], c.uv[1], c.uv[2], a.pos(), b.pos(), d.pos(), q)) }
}

pub fn flat_trace<'a>(c: &mut HmFlat<'a>, from: &HmLoc, q: Complex, stop: Face<'a>) -> Vec<HmLoc> {
    let dist = |c: &HmFlat<'_>, q: Complex| -> f64 {
        // std::min({a, b, c})
        let a = side_distance(c.uv[0], c.uv[1], q);
        let b = side_distance(c.uv[1], c.uv[2], q);
        let d = side_distance(c.uv[2], c.uv[0], q);
        let mut m = a;
        if b < m { m = b; }
        if d < m { m = d; }
        m
    };

    let exit = |c: &HmFlat<'_>, u: Complex| -> (i32, f64) {
        let mut k = -1;
        let mut r = 0.;
        let mut t = 2.;
        for j in 0..3 {
            let Some((se, tp)) = find_extended_intersection(c.uv[j], c.uv[(j + 1) % 3], u, q) else { continue };
            if tp > EPS && se > -EPS && se < 1. + EPS && tp < t { t = tp; k = j as i32; r = se; }
        }
        (k, r)
    };

    let enter = |c: &HmFlat<'a>, loc: &HmLoc, u: Complex, stop: Face<'a>| -> HmFlat<'a> {
        let mut fan = vec![*c];
        if let HmLoc::OnE { id, .. } = *loc { fan.push(flat_across(c, c.f.m.edge(id))); }
        if let HmLoc::OnV { id } = *loc { fan = flat_around(c, c.f.m.vert(id), stop, true); }
        let step = u + (q - u) * 1e-6;
        if dist(&fan[0], step) > 0. { return fan[0]; }
        // rg::max_element(fan, {}, proj): the first maximum
        let mut best = 0;
        let mut best_d = dist(&fan[0], step);
        for i in 1..fan.len() {
            let di = dist(&fan[i], step);
            if best_d < di { best = i; best_d = di; }
        }
        fan[best]
    };

    let mut path = vec![];
    let mut u = flat_uv(c, from);
    *c = enter(c, from, u, stop);
    while dist(c, q) < -EPS {
        let (k, r) = exit(c, u);
        let k = k as usize;
        u = lerp(c.uv[k], c.uv[(k + 1) % 3], r);
        path.push(flat_loc(c, u));
        *c = if matches!(path.last().unwrap(), HmLoc::OnV { .. }) {
            enter(c, path.last().unwrap(), u, c.f)
        } else {
            flat_across(c, c.f.halfs()[k].edge())
        };
    }
    path.push(flat_loc(c, q));
    path
}
