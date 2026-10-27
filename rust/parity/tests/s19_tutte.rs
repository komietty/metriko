//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Stages 19–22: the mesh cut along the collapsed T-mesh, its cut mesh, the Tutte embedding of every patch and the
//! seamless SLIM relaxation, each from the dumped inputs of the previous stage (the Emesh from 18_tm_repaired).
//!
//! 19 and 20 are combinatorics and plain arithmetic: exact / bitwise. Tutte solves sparse LU systems with faer
//! instead of Eigen: observed max abs difference 1.7e-13, tolerance 1e-9. SLIM iterates LLT solves: observed max abs
//! 2.7e-11 (relative 6e-9), tolerance 1e-7 abs + rel; its iteration count is identical (checked against an
//! instrumented copy of the C++ by the porting worker).

use metriko::igm::parameterization::compute_cut_mesh;
use metriko::tutte::tutte::HalfData;
use metriko::tutte::tutte_cutting::compute_embedding_cut_hmesh;
use metriko::tutte::tutte_params::compute_tutte_parameterization;
use metriko::tutte::tutte_slim::seamless_slim;
use parity::*;
use parity::stage::*;

fn hdata(g: &Golden) -> Vec<HalfData> {
    let ids = g.load("19_emb.hdata_ids").rows_i32::<4>();
    let vals = g.load("19_emb.hdata_vals").rows_f64::<2>();
    ids.iter().zip(&vals).map(|(a, b)| HalfData { half: a[0], thid: a[1], tqid: a[2], order: a[3], v0: b[0], v1: b[1] }).collect()
}

#[test]
fn embedding_and_cut() {
    for (c, g) in golden_cases() {
        let n = c.name();
        let hm = hmesh(&g, "01_hm", true);
        let tm = emesh(&hm, &g, "18_tm_repaired");
        check_emesh(&n, &tm, &g, "18_tm_repaired", None).unwrap_or_else(|e| panic!("{e}"));
        let seam = g.load("03_seam.seam").bools();
        let matching = g.load("05_cmbf.matching").i32s().to_vec();
        let singular = g.load("05_cmbf.singular").i32s().to_vec();
        let (emb, seam1, matching1, singular1, hd) = compute_embedding_cut_hmesh(&hm, &tm, &seam, &matching, &singular);
        check_hmesh(&n, &emb, &g, "19_emb", true).unwrap_or_else(|e| panic!("{e}"));
        assert_eq_bool(&format!("{n} seam1"), &seam1, &g.load("19_emb.seam1").bools());
        assert_eq_i32(&format!("{n} matching1"), &matching1, g.load("19_emb.matching1").i32s(), 1);
        assert_eq_i32(&format!("{n} singular1"), &singular1, g.load("19_emb.singular1").i32s(), 1);
        let ids: Vec<i32> = hd.iter().flat_map(|d| [d.half, d.thid, d.tqid, d.order]).collect();
        let vals: Vec<f64> = hd.iter().flat_map(|d| [d.v0, d.v1]).collect();
        assert_eq_i32(&format!("{n} hdata_ids"), &ids, g.load("19_emb.hdata_ids").i32s(), 4);
        assert_bits_f64(&format!("{n} hdata_vals"), &vals, g.load("19_emb.hdata_vals").f64s(), 2);

        let hcut = compute_cut_mesh(&emb, &seam1);
        check_hmesh(&n, &hcut, &g, "20_hcut", false).unwrap_or_else(|e| panic!("{e}"));
    }
}

#[test]
fn tutte_and_slim() {
    for (c, g) in golden_cases() {
        let n = c.name();
        let hm = hmesh(&g, "01_hm", true);
        let tm = emesh(&hm, &g, "18_tm_repaired");
        let emb = hmesh(&g, "19_emb", true);
        let hcut = hmesh(&g, "20_hcut", false);
        let seam1 = g.load("19_emb.seam1").bools();
        let singular1 = g.load("19_emb.singular1").i32s().to_vec();

        let uv_ebd = compute_tutte_parameterization(&emb, &tm, &seam1, &hdata(&g));
        assert_close_f64(&format!("{n} 21_tutte.uv_ebd"), &flat(&uv_ebd), g.load("21_tutte.uv_ebd").f64s(), 2, Tol::new(1e-9, 0.));
        // lib.h: uv_cut(v) = uv_ebd(v.half().next().crnr()), bitwise on the dumped uv_ebd
        let cpp_ebd = g.load("21_tutte.uv_ebd").rows_f64::<2>();
        let gather: Vec<[f64; 2]> = hcut.verts().map(|v| cpp_ebd[v.half().next().crnr().id as usize]).collect();
        assert_bits_f64(&format!("{n} 21_tutte.uv_cut"), &flat(&gather), g.load("21_tutte.uv_cut").f64s(), 2);

        let uv_cut = g.load("21_tutte.uv_cut").rows_f64::<2>();
        let uv = seamless_slim(&emb, &hcut, &seam1, &singular1, &uv_cut, metriko::SLIM_MAX_ITER, metriko::SLIM_REL_TOL);
        assert_close_f64(&format!("{n} 22_slim.uv"), &flat(&uv), g.load("22_slim.uv").f64s(), 2, Tol::new(1e-7, 1e-7));
    }
}
