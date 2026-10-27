//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/common/predicates.h (see rust/PORTED_FROM).

use crate::common::typedef::*;

// Todo: still using epsilon predicates... consider using strict predicates
pub fn orientation(pa: Complex, pb: Complex, pc: Complex) -> f64 {
    let acx = pa.re - pc.re;
    let bcx = pb.re - pc.re;
    let acy = pa.im - pc.im;
    let bcy = pb.im - pc.im;
    acx * bcy - acy * bcx
}

// Todo: consider using side_distance or side_ratio
/// `is_collinear(pa, pb, pc, eps = EPS)`
pub fn is_collinear(pa: Complex, pb: Complex, pc: Complex, eps: f64) -> bool {
    let o = orientation(pa, pb, pc);
    o.abs() < eps
}

// Todo: consider using side_distance or side_ratio
/// `is_points_into(p1, p2, p3, uv, eps = EPS)`
pub fn is_points_into(p1: Complex, p2: Complex, p3: Complex, uv: Complex, eps: f64) -> bool {
    orientation(p1, p2, uv) > eps && orientation(p1, p3, uv) < -eps
}

// Todo: consider using side_distance or side_ratio
pub fn is_inside_triangle(pa: Complex, pb: Complex, pc: Complex, uv: Complex) -> bool {
    orientation(pa, pb, uv) > EPS &&
    orientation(pb, pc, uv) > EPS &&
    orientation(pc, pa, uv) > EPS
}

pub fn side_distance(a: Complex, b: Complex, p: Complex) -> f64 {
    orientation(a, b, p) / (b - a).norm()
}

pub fn side_ratio(a: Complex, b: Complex, p: Complex) -> f64 {
    orientation(a, b, p) / (b - a).norm_sqr()
}

pub fn is_inside_segment(a: Complex, b: Complex, p: Complex) -> bool {
    let ab = b - a;
    let ap = p - a;
    let l = ab.norm();
    if l < EPS { return false; }
    let d = (ab.re * ap.re + ab.im * ap.im) / (l * l);
    side_ratio(a, b, p).abs() < EPS && d > EPS && d < 1. - EPS
}
