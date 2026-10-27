//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/hmesh/hmray.h (see rust/PORTED_FROM).

use crate::common::typedef::*;
use crate::common::utilities::find_extended_intersection;
use crate::common::predicates::is_points_into;
use crate::common::stdcomplex::std_polar;
use crate::hmesh::hmesh::*;
use crate::hmesh::hmloc::*;

/// `find_ray_intersection(Half h, complex o, complex d, double tol)`
pub fn find_ray_intersection_half(h: Half<'_>, o: Complex, d: Complex, tol: f64) -> Option<HmLoc> {
    let a = h.cr_t().uv();
    let b = h.cr_h().uv();
    let (s, t) = find_extended_intersection(a, b, o, o + d)?;
    if t <= tol || s < -tol || s > 1. + tol { return None; }
    if s < tol { return Some(HmLoc::OnC { id: h.cr_t().id }); }
    if s > 1. - tol { return Some(HmLoc::OnC { id: h.cr_h().id }); }
    Some(HmLoc::OnH { id: h.id, r: s })
}

/// `find_ray_intersection(const Hmesh&, const HmLoc&, complex dir)`
pub fn find_ray_intersection(hm: &Hmesh, loc: &HmLoc, dir: Complex) -> HmLoc {
    let tol = EPS;
    match *loc {
        HmLoc::OnC { id } => {
            if let Some(x) = find_ray_intersection_half(hm.crnr(id).half(), hm.crnr(id).uv(), dir, tol) { return x; }
            metriko_fail!("no intersection");
        }
        HmLoc::OnH { id, r } => {
            let h = hm.half(id);
            if let Some(x) = find_ray_intersection_half(h.next(), h.uv(r), dir, tol) { return x; }
            if let Some(x) = find_ray_intersection_half(h.prev(), h.uv(r), dir, tol) { return x; }
            metriko_fail!("no intersection");
        }
        HmLoc::OnP { id, uv } => {
            let hs = hm.face(id).halfs();
            if let Some(x) = find_ray_intersection_half(hs[0], uv, dir, tol) { return x; }
            if let Some(x) = find_ray_intersection_half(hs[1], uv, dir, tol) { return x; }
            if let Some(x) = find_ray_intersection_half(hs[2], uv, dir, tol) { return x; }
            metriko_fail!("no intersection");
        }
        _ => metriko_fail!("no impl"),
    }
}

pub fn cross_to_twin(hm: &Hmesh, loc: &HmLoc, matching: &[i32], dir: Complex) -> (HmLoc, Complex) {
    let get_m = |h: Half<'_>| (if h.is_canonical() { -1 } else { 1 }) * matching[h.edge().id as usize];

    match *loc {
        HmLoc::OnC { id } => {
            let c = hm.crnr(id);
            let v = c.vert();
            let mut d = dir;
            for h in v.adj_halfs_from(c.half().next().twin(), true) {
                d *= std_polar(1., PI / 2. * get_m(h) as f64);
                let c1 = h.next().crnr();
                let uv0 = c1.uv();
                let uv1 = c1.next().uv();
                let uv2 = c1.prev().uv();
                if is_points_into(uv0, uv1, uv2, uv0 + d, 0.) && c1 != c { return (HmLoc::OnC { id: c1.id }, d); }
            }
            metriko_fail!("no face around vert {} admits the ray direction", v.id);
        }
        HmLoc::OnH { id, r } => {
            let h = hm.half(id).twin();
            (HmLoc::OnH { id: h.id, r: 1. - r }, std_polar(1., PI / 2. * get_m(h) as f64) * dir)
        }
        _ => metriko_fail!("no impl"),
    }
}
