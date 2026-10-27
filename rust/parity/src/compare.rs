//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Comparisons of a Rust result against a C++ dump. Discrete data must be equal; floating data either bitwise or
//! within `|rust - cpp| <= abs + rel * max(|rust|, |cpp|)`. Every check returns `Err(report)` with the number of
//! mismatches, the worst errors and the first few mismatching entries (as (row, col) when `cols` > 1); the
//! `assert_*` wrappers panic with that report.

use metriko::Complex;

#[derive(Clone, Copy, Debug)]
pub struct Tol {
    pub abs: f64,
    pub rel: f64,
}

impl Tol {
    pub const fn new(abs: f64, rel: f64) -> Tol { Tol { abs, rel } }
    pub fn ok(&self, a: f64, b: f64) -> bool {
        if a.to_bits() == b.to_bits() { return true; }
        if a.is_nan() || b.is_nan() { return false; }
        (a - b).abs() <= self.abs + self.rel * a.abs().max(b.abs())
    }
}

const SHOW: usize = 8;

fn at(i: usize, cols: usize) -> String {
    if cols > 1 { format!("[{}, {}]", i / cols, i % cols) } else { format!("[{i}]") }
}

fn check_len(what: &str, n_rust: usize, n_cpp: usize) -> Result<(), String> {
    if n_rust != n_cpp { Err(format!("{what}: {n_rust} entries in rust, {n_cpp} in c++")) } else { Ok(()) }
}

fn check_by<T: std::fmt::Debug>(what: &str, rust: &[T], cpp: &[T], cols: usize, same: impl Fn(&T, &T) -> bool) -> Result<(), String> {
    check_len(what, rust.len(), cpp.len())?;
    let bad: Vec<usize> = (0..rust.len()).filter(|&i| !same(&rust[i], &cpp[i])).collect();
    if bad.is_empty() { return Ok(()); }
    let mut s = format!("{what}: {}/{} entries differ", bad.len(), rust.len());
    for &i in bad.iter().take(SHOW) { s += &format!("\n  {} rust {:?} c++ {:?}", at(i, cols), rust[i], cpp[i]); }
    if bad.len() > SHOW { s += "\n  ..."; }
    Err(s)
}

pub fn check_eq_i32(what: &str, rust: &[i32], cpp: &[i32], cols: usize) -> Result<(), String> { check_by(what, rust, cpp, cols, |a, b| a == b) }
pub fn check_eq_bool(what: &str, rust: &[bool], cpp: &[bool]) -> Result<(), String> { check_by(what, rust, cpp, 1, |a, b| a == b) }

/// bitwise equality (NaNs with the same payload are equal, 0.0 and -0.0 are not)
pub fn check_bits_f64(what: &str, rust: &[f64], cpp: &[f64], cols: usize) -> Result<(), String> {
    check_close_f64(what, rust, cpp, cols, Tol::new(-1., 0.))
}

pub fn check_bits_c64(what: &str, rust: &[Complex], cpp: &[Complex], cols: usize) -> Result<(), String> {
    check_close_c64(what, rust, cpp, cols, Tol::new(-1., 0.))
}

/// within the tolerance (a negative `abs` with `rel` 0 means bitwise)
pub fn check_close_f64(what: &str, rust: &[f64], cpp: &[f64], cols: usize, tol: Tol) -> Result<(), String> {
    check_len(what, rust.len(), cpp.len())?;
    let bad: Vec<usize> = (0..rust.len()).filter(|&i| !tol.ok(rust[i], cpp[i])).collect();
    if bad.is_empty() { return Ok(()); }
    let (mut ia, mut ea, mut ir, mut er) = (0, 0f64, 0, 0f64);
    for &i in &bad {
        let e = (rust[i] - cpp[i]).abs();
        let r = e / rust[i].abs().max(cpp[i].abs());
        if e > ea || e.is_nan() { ea = e; ia = i; }
        if r > er || r.is_nan() { er = r; ir = i; }
    }
    let mut s = format!("{what}: {}/{} entries outside abs {:e} + rel {:e}; max abs err {:e} at {}, max rel err {:e} at {}",
        bad.len(), rust.len(), tol.abs.max(0.), tol.rel, ea, at(ia, cols), er, at(ir, cols));
    for &i in bad.iter().take(SHOW) { s += &format!("\n  {} rust {:e} c++ {:e} (diff {:e})", at(i, cols), rust[i], cpp[i], rust[i] - cpp[i]); }
    if bad.len() > SHOW { s += "\n  ..."; }
    Err(s)
}

/// real and imaginary parts within the tolerance
pub fn check_close_c64(what: &str, rust: &[Complex], cpp: &[Complex], cols: usize, tol: Tol) -> Result<(), String> {
    let re = |v: &[Complex]| -> Vec<f64> { v.iter().flat_map(|c| [c.re, c.im]).collect() };
    check_close_f64(&format!("{what} (re, im interleaved)"), &re(rust), &re(cpp), 2 * cols.max(1), tol)
}

pub fn assert_eq_i32(what: &str, rust: &[i32], cpp: &[i32], cols: usize) { if let Err(e) = check_eq_i32(what, rust, cpp, cols) { panic!("{e}"); } }
pub fn assert_eq_bool(what: &str, rust: &[bool], cpp: &[bool]) { if let Err(e) = check_eq_bool(what, rust, cpp) { panic!("{e}"); } }
pub fn assert_bits_f64(what: &str, rust: &[f64], cpp: &[f64], cols: usize) { if let Err(e) = check_bits_f64(what, rust, cpp, cols) { panic!("{e}"); } }
pub fn assert_bits_c64(what: &str, rust: &[Complex], cpp: &[Complex], cols: usize) { if let Err(e) = check_bits_c64(what, rust, cpp, cols) { panic!("{e}"); } }
pub fn assert_close_f64(what: &str, rust: &[f64], cpp: &[f64], cols: usize, tol: Tol) { if let Err(e) = check_close_f64(what, rust, cpp, cols, tol) { panic!("{e}"); } }
pub fn assert_close_c64(what: &str, rust: &[Complex], cpp: &[Complex], cols: usize, tol: Tol) { if let Err(e) = check_close_c64(what, rust, cpp, cols, tol) { panic!("{e}"); } }

/// flattens rows of N for the checks above
pub fn flat<T: Copy, const N: usize>(rows: &[[T; N]]) -> Vec<T> { rows.iter().flatten().copied().collect() }

/// flattens points (x, y, z per row)
pub fn flat3(p: &[glam::DVec3]) -> Vec<f64> { p.iter().flat_map(|v| v.to_array()).collect() }
