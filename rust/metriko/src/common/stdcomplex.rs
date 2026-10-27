//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- libc++ 20.1 / compiler-rt arithmetic of std::complex<double> as the reference build executes it.

use crate::common::typedef::*;

// The complex functions of libc++ 20.1 that this file calls, reproduced exactly (num_complex's `powf` / `powc` /
// `exp` / `ln` and its `/` use other formulas). A product of two complex numbers is clang's inline `_Complex`
// multiplication, the plain formula (no contraction under the reference `-ffp-contract=off`), as num_complex's `*`.
// A quotient of two complex numbers is a call to compiler-rt's `__divdc3`, see `std_div`.

/// `complex / complex`: libc++ 20.1 divides through the builtin `_Complex` division, which clang lowers to a call to
/// compiler-rt's `__divdc3`. That function is precompiled in the system library with floating point contraction, so
/// its `c*c + d*d`, `a*c + b*d` and `b*c - a*d` are fused multiply-adds whatever the flags of the reference build
/// (verified bitwise on the golden `02_rawf.connection` of spot and fandisk: 100% with the fused form below, 57% with
/// the plain formula).
pub fn std_div(z: Complex, w: Complex) -> Complex {
    let (mut a, mut b, mut c, mut d) = (z.re, z.im, w.re, w.im);
    let mut ilogbw = 0;
    let logbw = logb(c.abs().max(d.abs()));
    if logbw.is_finite() {
        ilogbw = logbw as i32;
        c = scalbn(c, -ilogbw);
        d = scalbn(d, -ilogbw);
    }
    let denom = c.mul_add(c, d * d);
    let mut x = scalbn(a.mul_add(c, b * d) / denom, -ilogbw);
    let mut y = scalbn(b.mul_add(c, -(a * d)) / denom, -ilogbw);
    if x.is_nan() && y.is_nan() {
        if denom == 0. && (!a.is_nan() || !b.is_nan()) {
            x = f64::INFINITY.copysign(c) * a;
            y = f64::INFINITY.copysign(c) * b;
        } else if (a.is_infinite() || b.is_infinite()) && c.is_finite() && d.is_finite() {
            a = (if a.is_infinite() { 1f64 } else { 0. }).copysign(a);
            b = (if b.is_infinite() { 1f64 } else { 0. }).copysign(b);
            x = f64::INFINITY * a.mul_add(c, b * d);
            y = f64::INFINITY * b.mul_add(c, -(a * d));
        } else if logbw.is_infinite() && logbw > 0. && a.is_finite() && b.is_finite() {
            c = (if c.is_infinite() { 1f64 } else { 0. }).copysign(c);
            d = (if d.is_infinite() { 1f64 } else { 0. }).copysign(d);
            x = 0. * a.mul_add(c, b * d);
            y = 0. * b.mul_add(c, -(a * d));
        }
    }
    cx(x, y)
}

/// C `logb`: the unbiased exponent of `x` as a double (-inf for 0, inf for inf, NaN for NaN)
pub fn logb(x: f64) -> f64 {
    if x.is_nan() { return x; }
    if x == 0. { return f64::NEG_INFINITY; }
    if x.is_infinite() { return f64::INFINITY; }
    let bits = x.to_bits();
    let e = ((bits >> 52) & 0x7ff) as i64;
    if e != 0 { return (e - 1023) as f64; }
    // subnormal: the position of the leading one of the mantissa
    let m = bits & ((1u64 << 52) - 1);
    (-1022 - (m.leading_zeros() as i64 - 11)) as f64
}

/// C `scalbn`: `x * 2^n`, exactly when representable, correctly rounded otherwise
pub fn scalbn(x: f64, mut n: i32) -> f64 {
    let x1p1023 = f64::from_bits(0x7fe0_0000_0000_0000);
    let x1p53 = f64::from_bits(0x4340_0000_0000_0000);
    let x1p_1022 = f64::from_bits(0x0010_0000_0000_0000);
    let mut y = x;
    if n > 1023 {
        y *= x1p1023;
        n -= 1023;
        if n > 1023 {
            y *= x1p1023;
            n -= 1023;
            if n > 1023 { n = 1023; }
        }
    } else if n < -1022 {
        // make sure the final result is rounded once, by scaling into the subnormal range in one step at the end
        y *= x1p_1022 * x1p53;
        n += 1022 - 53;
        if n < -1022 {
            y *= x1p_1022 * x1p53;
            n += 1022 - 53;
            if n < -1022 { n = -1022; }
        }
    }
    y * f64::from_bits(((0x3ff + n) as u64) << 52)
}

/// libc++ `std::log(complex)`: `(log(abs(x)), arg(x))`
pub fn std_log(x: Complex) -> Complex { cx(x.norm().ln(), x.arg()) }

/// libc++ `std::exp(complex)`
pub fn std_exp(x: Complex) -> Complex {
    let mut i = x.im;
    if i == 0. {
        return cx(x.re.exp(), 0f64.copysign(x.im));
    }
    if x.re.is_infinite() {
        if x.re < 0. {
            if !i.is_finite() { i = 1.; }
        } else if i == 0. || !i.is_finite() {
            if i.is_infinite() { i = f64::NAN; }
            return cx(x.re, i);
        }
    }
    let e = x.re.exp();
    cx(e * i.cos(), e * i.sin())
}

/// libc++ `std::pow(complex, complex)` = `exp(y * log(x))`; `pow(complex, int)` and `pow(complex, double)` promote
/// the exponent to `complex(y, 0)` and call it.
pub fn std_pow(x: Complex, y: Complex) -> Complex { std_exp(y * std_log(x)) }

/// libc++ `std::polar(rho, theta)`
pub fn std_polar(rho: f64, theta: f64) -> Complex {
    if rho.is_nan() || rho.is_sign_negative() { return cx(f64::NAN, f64::NAN); }
    if theta.is_nan() {
        if rho.is_infinite() { return cx(rho, theta); }
        return cx(theta, theta);
    }
    if theta.is_infinite() {
        if rho.is_infinite() { return cx(rho, f64::NAN); }
        return cx(f64::NAN, f64::NAN);
    }
    let mut x = rho * theta.cos();
    if x.is_nan() { x = 0.; }
    let mut y = rho * theta.sin();
    if y.is_nan() { y = 0.; }
    cx(x, y)
}

/// `std::round` to the nearest integer, halfway cases away from zero
#[inline]
pub fn std_round(x: f64) -> f64 { x.round() }
