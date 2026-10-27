//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/common/utilities.h (see rust/PORTED_FROM).

use crate::common::typedef::*;

/// `circular_prev` on an index into a sequence of length `n`.
#[inline]
pub fn circular_prev(n: usize, i: usize) -> usize { if i == 0 { n - 1 } else { i - 1 } }

/// `circular_next` on an index into a sequence of length `n`.
#[inline]
pub fn circular_next(n: usize, i: usize) -> usize { if i + 1 == n { 0 } else { i + 1 } }

pub fn grid_unit(v: &[Row3d], scale: f64) -> f64 {
    // V.colwise().maxCoeff() - V.colwise().minCoeff(): the first extremum of each column, by strict comparison
    let mut mx = v[0];
    let mut mn = v[0];
    for p in &v[1..] {
        if p.x > mx.x { mx.x = p.x; }
        if p.y > mx.y { mx.y = p.y; }
        if p.z > mx.z { mx.z = p.z; }
        if p.x < mn.x { mn.x = p.x; }
        if p.y < mn.y { mn.y = p.y; }
        if p.z < mn.z { mn.z = p.z; }
    }
    (mx - mn).length() * scale
}

pub fn get_quater_rot(i: i32) -> Complex {
    match (i % 4 + 4) % 4 {
        0 => cx(1., 0.),
        1 => cx(0., 1.),
        2 => cx(-1., 0.),
        3 => cx(0., -1.),
        _ => panic!("Invalid index"),
    }
}

/// `equal(a, b, tolerance = 1e-10)`
pub fn equal(a: Complex, b: Complex, tolerance: f64) -> bool { (a - b).norm() < tolerance }

pub fn lerp(a: Complex, b: Complex, t: f64) -> Complex { a * (1. - t) + b * t }

pub fn normalize(a: Complex) -> Complex {
    let l = a.norm();
    metriko_check!(l > 0., "normalize of a zero vector");
    a / l
}

pub fn dot(a: Complex, b: Complex) -> f64 { a.re * b.re + a.im * b.im }

pub fn cross(a: Complex, b: Complex) -> f64 { a.re * b.im - a.im * b.re }

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MinmaxInt {
    pub min_x: i32,
    pub min_y: i32,
    pub max_x: i32,
    pub max_y: i32,
}

pub fn get_minmax_int(uvs: &[Complex]) -> MinmaxInt {
    // rg::min / rg::max over the projections: the extremum values
    let mut min_x = uvs[0].re;
    let mut min_y = uvs[0].im;
    let mut max_x = uvs[0].re;
    let mut max_y = uvs[0].im;
    for uv in &uvs[1..] {
        if uv.re < min_x { min_x = uv.re; }
        if uv.im < min_y { min_y = uv.im; }
        if max_x < uv.re { max_x = uv.re; }
        if max_y < uv.im { max_y = uv.im; }
    }
    MinmaxInt {
        min_x: min_x.floor() as i32,
        min_y: min_y.floor() as i32,
        max_x: max_x.ceil() as i32,
        max_y: max_y.ceil() as i32,
    }
}

pub fn calc_coefficient(ori: Complex, pt1: Complex, pt2: Complex, tgt: Complex) -> Complex {
    let v1 = pt1 - ori;
    let v2 = pt2 - ori;
    let dif = tgt - ori;
    let det = cross(v1, v2);
    metriko_check!(det != 0., "degenerate uv triangle");
    cx(cross(dif, v2) / det, cross(v1, dif) / det)
}

/// Find the intersection of a line passing through a and b and another line passing through c and d.
/// Returns `(ratio_a2b, ratio_c2d)` if the segments are not parallel. Check 0 <= ratio <= 1 if you want a
/// segment-segment intersection.
pub fn find_extended_intersection(a: Complex, b: Complex, c: Complex, d: Complex) -> Option<(f64, f64)> {
    let e1 = d - c;
    let e2 = b - a;
    let w = a - c;
    let det = cross(e1, e2);
    if det == 0. { return None; }
    let ratio_c2d = cross(w, e2) / det;
    let ratio_a2b = cross(w, e1) / det;
    Some((ratio_a2b, ratio_c2d))
}

/// `find_strict_intersection(a, b, c, d, eps = EPS)`: segments AB and CD cross strictly inside both.
pub fn find_strict_intersection(a: Complex, b: Complex, c: Complex, d: Complex, eps: f64) -> bool {
    find_strict_intersection_ratio(a, b, c, d, eps).is_some()
}

/// The overload of `find_strict_intersection` that also writes the ratios in AB and CD: `Some((ratio_a2b,
/// ratio_c2d))` when it returns true.
pub fn find_strict_intersection_ratio(a: Complex, b: Complex, c: Complex, d: Complex, eps: f64) -> Option<(f64, f64)> {
    let (ratio_a2b, ratio_c2d) = find_extended_intersection(a, b, c, d)?;
    if ratio_a2b >= 0. + eps && ratio_a2b <= 1. - eps && ratio_c2d >= 0. + eps && ratio_c2d <= 1. - eps {
        Some((ratio_a2b, ratio_c2d))
    } else {
        None
    }
}

pub fn conversion_2d_3d(
    o2: Complex, // origin of 2d
    a2: Complex,
    b2: Complex,
    o3: Row3d,   // origin of 3d
    a3: Row3d,
    b3: Row3d,
    uv: Complex, // target uv value
) -> Row3d {
    let c = calc_coefficient(o2, a2, b2, uv);
    o3 + (a3 - o3) * c.re + (b3 - o3) * c.im
}
