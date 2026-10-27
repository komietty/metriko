//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/common/typedef.h (see rust/PORTED_FROM).

pub type Complex = num_complex::Complex<f64>;
pub type Row3d = glam::DVec3;
pub type Row2d = glam::DVec2;

pub const PI: f64 = std::f64::consts::PI;
pub const TWO_PI: f64 = 2. * std::f64::consts::PI;
pub const EPS: f64 = 1e-9;

pub const fn cx(re: f64, im: f64) -> Complex { Complex::new(re, im) }

/// Eigen's `v.normalized()`: `v / v.norm()` component by component, `v` itself when its norm is 0.
#[inline]
pub fn normalized(v: Row3d) -> Row3d {
    let n = v.length_squared();
    if n > 0. { v / n.sqrt() } else { v }
}
