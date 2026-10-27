//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Stages 04, 07, 08: the cut mesh, the setup of the integer grid map and its integration. Inputs from the dumps
//! (01_hm, 03_seam, 05_cmbf, 06_extf, 00_input.scale).
//!
//! The cut mesh and the setup are combinatorics plus plain arithmetic: exact / bitwise. The integration solves
//! sparse systems with faer instead of Eigen, so `fullx` / `nfn` / `cfn` agree up to rounding: observed max abs
//! difference 6.2e-10 for values up to ~180 (relative ~1e-12); the tolerance is abs 1e-8. Its return value and the
//! sign of every parameterized triangle must agree exactly.

use faer::Mat;
use metriko::igm::parameterization::{compute_cut_mesh, RosyParameterization};
use metriko::solver::matrix_ops::{nonzeros, SprsD};
use parity::*;
use parity::stage::*;

const SOLVER_TOL: Tol = Tol::new(1e-8, 0.);

fn sparse_eq(name: &str, g: &Golden, key: &str, m: &SprsD) {
    let s: Vec<i64> = g.meta(&format!("07_rp_setup.{key}_shape")).as_array().unwrap().iter().map(|v| v.as_i64().unwrap()).collect();
    assert_eq!((m.nrows() as i64, m.ncols() as i64), (s[0], s[1]), "{name} {key} shape");
    let mine: Vec<(usize, usize, f64)> = nonzeros(m.as_ref()).collect();
    let rc: Vec<i32> = mine.iter().flat_map(|&(r, c, _)| [r as i32, c as i32]).collect();
    assert_eq_i32(&format!("{name} {key} rc"), &rc, g.load(&format!("07_rp_setup.{key}_rc")).i32s(), 2);
    let v: Vec<f64> = mine.iter().map(|&(_, _, v)| v).collect();
    assert_bits_f64(&format!("{name} {key} val"), &v, g.load(&format!("07_rp_setup.{key}_val")).f64s(), 1);
}

fn rows(m: &Mat<f64>) -> Vec<f64> {
    (0..m.nrows()).flat_map(|i| (0..m.ncols()).map(move |j| (i, j))).map(|(i, j)| m[(i, j)]).collect()
}

fn signs(cfn: &[f64], nf: usize) -> Vec<bool> {
    (0..nf).map(|f| {
        let c = &cfn[f * 12..f * 12 + 12];
        (c[4] - c[0]) * (c[9] - c[1]) - (c[5] - c[1]) * (c[8] - c[0]) > 0.
    }).collect()
}

#[test]
fn cut_mesh_setup_and_integration() {
    for (c, g) in golden_cases() {
        let name = c.name();
        let n = |s: &str| format!("{name} {s}");
        let raw = hmesh(&g, "01_hm", true);
        let seam = g.load("03_seam.seam").bools();
        let cut = compute_cut_mesh(&raw, &seam);
        check_hmesh(&name, &cut, &g, "04_cutm", false).unwrap_or_else(|e| panic!("{e}"));

        let e = g.load("06_extf.extf");
        let ext = Mat::<f64>::from_fn(e.rows(), e.cols(), |i, j| e.f64s()[i * e.cols() + j]);
        let singular = g.load("05_cmbf.singular").i32s().to_vec();
        let matching = g.load("05_cmbf.matching").i32s().to_vec();
        let scale = g.load("00_input.scale").f64s()[0];

        let mut rp = RosyParameterization::new(&raw, &cut, &ext, &singular, &matching, &seam, 4, scale);
        rp.local_injectivity = true;
        rp.setup();
        assert_eq_i32(&n("he2matching"), &rp.he2matching, g.load("07_rp_setup.he2matching").i32s(), 1);
        assert_eq_i32(&n("he2transidx"), &rp.he2transidx, g.load("07_rp_setup.he2transidx").i32s(), 1);
        assert_eq_i32(&n("integerIdcs"), &rp.integer_idcs, g.load("07_rp_setup.integerIdcs").i32s(), 1);
        assert_eq_i32(&n("singularIdcs"), &rp.singular_idcs, g.load("07_rp_setup.singularIdcs").i32s(), 1);
        assert_eq_i32(&n("fixedIdcs"), &rp.fixed_idcs, g.load("07_rp_setup.fixedIdcs").i32s(), 1);
        assert_bits_f64(&n("fixedVals"), &rp.fixed_vals, g.load("07_rp_setup.fixedVals").f64s(), 1);
        let lr = g.load("07_rp_setup.lreductor");
        assert_eq_i32(&n("lreductor"), &rp.lreductor.data, lr.i32s(), lr.cols());
        sparse_eq(&name, &g, "vtrans2cut", &rp.vtrans2cut);
        sparse_eq(&name, &g, "constraint", &rp.constraint);
        sparse_eq(&name, &g, "uncompress", &rp.uncompress);

        let ret = rp.integ();
        assert_eq!(ret, g.meta("08_rp_integ.integ_returned").as_bool().unwrap(), "{name} integ returned");
        let fullx: Vec<f64> = (0..rp.fullx.nrows()).map(|i| rp.fullx[i]).collect();
        assert_close_f64(&n("08_rp_integ.fullx"), &fullx, g.load("08_rp_integ.fullx").f64s(), 1, SOLVER_TOL);
        assert_close_f64(&n("08_rp_integ.nfn"), &rows(&rp.nfn), g.load("08_rp_integ.nfn").f64s(), rp.nfn.ncols(), SOLVER_TOL);
        let cfn = rows(&rp.cfn);
        assert_close_f64(&n("08_rp_integ.cfn"), &cfn, g.load("08_rp_integ.cfn").f64s(), 12, SOLVER_TOL);
        let nf = raw.n_f as usize;
        assert_eq_bool(&n("uv triangle signs"), &signs(&cfn, nf), &signs(g.load("08_rp_integ.cfn").f64s(), nf));
    }
}
