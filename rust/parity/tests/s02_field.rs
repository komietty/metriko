//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Stages 02, 03, 05, 06: the smoothest 4-rosy field, its seam, the combed field and the extrinsic field.
//!
//! The raw field comes from the smallest eigenvector of the connection Laplacian, solved by faer instead of Eigen:
//! `compressed` agrees to ~1e-14 and the field to ~1e-11 (the accepted solver deviation), so they are compared
//! within 1e-9 (absolute, values are unit complex numbers). Everything else is arithmetic on its inputs: the
//! connection bitwise, the discrete data (matching, singular, seam) exactly. Combing and the extrinsic field are
//! checked in isolation from the dumped raw field / seam / combed field, bitwise.

use faer::{Col, Mat};
use metriko::nvec::face_rosy_field::*;
use metriko::Complex;
use parity::*;
use parity::stage::*;

const SOLVER_TOL: Tol = Tol::new(1e-9, 0.);

fn mat_rows(m: &Mat<Complex>) -> Vec<Complex> {
    (0..m.nrows()).flat_map(|i| (0..m.ncols()).map(move |j| m[(i, j)])).collect()
}
fn col(c: &Col<Complex>) -> Vec<Complex> { (0..c.nrows()).map(|i| c[i]).collect() }

/// a FaceRosyField holding the dumped arrays of `stage`
fn dumped_field<'a>(hm: &'a metriko::hmesh::hmesh::Hmesh, g: &Golden, stage: &str) -> FaceRosyField<'a> {
    let mut f = FaceRosyField::new_empty(hm, 4);
    let fl = g.load(&format!("{stage}.field"));
    let v = fl.c64s();
    f.field = Mat::from_fn(fl.rows(), fl.cols(), |i, j| v[i * fl.cols() + j]);
    let c = g.load(&format!("{stage}.connection"));
    f.connection = Col::from_fn(c.len(), |i| c.c64s()[i]);
    let c = g.load(&format!("{stage}.compressed"));
    f.compressed = Col::from_fn(c.len(), |i| c.c64s()[i]);
    f.matching = g.load(&format!("{stage}.matching")).i32s().to_vec();
    f.singular = g.load(&format!("{stage}.singular")).i32s().to_vec();
    f
}

#[test]
fn raw_field_and_seam() {
    for (c, g) in golden_cases() {
        let n = |s: &str| format!("{} {s}", c.name());
        let hm = hmesh(&g, "01_hm", true);
        let f = FaceRosyField::new(&hm, 4, FieldType::Smoothest);
        assert_bits_c64(&n("02_rawf.connection"), &col(&f.connection), g.load("02_rawf.connection").c64s(), 1);
        assert_close_c64(&n("02_rawf.compressed"), &col(&f.compressed), g.load("02_rawf.compressed").c64s(), 1, SOLVER_TOL);
        assert_close_c64(&n("02_rawf.field"), &mat_rows(&f.field), g.load("02_rawf.field").c64s(), 4, SOLVER_TOL);
        assert_eq_i32(&n("02_rawf.matching"), &f.matching, g.load("02_rawf.matching").i32s(), 1);
        assert_eq_i32(&n("02_rawf.singular"), &f.singular, g.load("02_rawf.singular").i32s(), 1);
        // chained and isolated seam
        assert_eq_bool(&n("03_seam (chained)"), &compute_seam(&f), &g.load("03_seam.seam").bools());
        let rawf = dumped_field(&hm, &g, "02_rawf");
        assert_eq_bool(&n("03_seam"), &compute_seam(&rawf), &g.load("03_seam.seam").bools());
    }
}

#[test]
fn combed_and_extrinsic_field() {
    for (c, g) in golden_cases() {
        let n = |s: &str| format!("{} {s}", c.name());
        let hm = hmesh(&g, "01_hm", true);
        let rawf = dumped_field(&hm, &g, "02_rawf");
        let seam = g.load("03_seam.seam").bools();
        let cmbf = compute_combbed_field(&rawf, &seam);
        assert_bits_c64(&n("05_cmbf.field"), &mat_rows(&cmbf.field), g.load("05_cmbf.field").c64s(), 4);
        assert_bits_c64(&n("05_cmbf.connection"), &col(&cmbf.connection), g.load("05_cmbf.connection").c64s(), 1);
        assert_bits_c64(&n("05_cmbf.compressed"), &col(&cmbf.compressed), g.load("05_cmbf.compressed").c64s(), 1);
        assert_eq_i32(&n("05_cmbf.matching"), &cmbf.matching, g.load("05_cmbf.matching").i32s(), 1);
        assert_eq_i32(&n("05_cmbf.singular"), &cmbf.singular, g.load("05_cmbf.singular").i32s(), 1);

        let cmbf = dumped_field(&hm, &g, "05_cmbf");
        let extf = compute_extrinsic_field(&cmbf, 4);
        let e: Vec<f64> = (0..extf.nrows()).flat_map(|i| (0..extf.ncols()).map(move |j| (i, j))).map(|(i, j)| extf[(i, j)]).collect();
        assert_bits_f64(&n("06_extf"), &e, g.load("06_extf.extf").f64s(), 12);
    }
}
