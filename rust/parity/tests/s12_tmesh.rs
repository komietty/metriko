//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Stages 12–18: quantization of the T-mesh, the collapse of zero lengths, the snaps, the dedup and the repair of
//! crossing tedges. All of it is combinatorics and plain arithmetic on its inputs, so everything is compared exactly
//! (structure) and bitwise (values).
//!
//! Inputs: the mesh (01_hm) with its corner uvs (09_hm_cfn) and the combed matching / singular (05_cmbf), from which
//! Mgrph and Emesh are rebuilt (they match stages 10 / 11, see s10_mg). An Emesh cannot be rebuilt from a dump
//! mid-collapse, so the stages are chained from 11 on, each one compared with its dump.

use metriko::quantization::quantization::{compute_quantization, validate_quantization};
use metriko::tmesh::emesh::Emesh;
use metriko::tmesh::emesh_validate::{repair_crossing_tedges, validate_collapse_done};
use metriko::tmesh::motorcycle::Mgrph;
use parity::*;
use parity::stage::*;

#[test]
fn quantization_collapse_snap_repair() {
    for (c, g) in golden_cases() {
        let name = c.name();
        let ck = |r: Result<(), String>| r.unwrap_or_else(|e| panic!("{e}"));
        let mut hm = hmesh(&g, "01_hm", true);
        hm.cfn = g.load("09_hm_cfn.cfn").c64s().to_vec();
        let matching = g.load("05_cmbf.matching").i32s().to_vec();
        let singular = g.load("05_cmbf.singular").i32s().to_vec();
        let mg = Mgrph::new(&hm, &matching, &singular);
        let mut tm = Emesh::new(&mg);

        let x = compute_quantization(&tm, &mg);
        assert_bits_f64(&format!("{name} 12_quantization.X"), &x, g.load("12_quantization.X").f64s(), 1);
        validate_quantization(&tm, &x);
        tm.set_x(&x);
        ck(check_emesh(&name, &tm, &g, "13_tm_setx", None));

        let rounds = metriko::collapse_zero_length(&mut tm);
        assert_eq!(rounds as i64, g.meta_i64("14_tm_collapsed.rounds"), "{name} collapse rounds");
        validate_collapse_done(&tm);
        ck(check_emesh(&name, &tm, &g, "14_tm_collapsed", None));

        tm.collapse_tedge_snap(false);
        ck(check_emesh(&name, &tm, &g, "15_tm_snap0", None));
        tm.collapse_tedge_snap(true);
        ck(check_emesh(&name, &tm, &g, "16_tm_snap1", None));
        for k in 0..tm.tedges.len() {
            if tm.tedges[k].id == -1 { continue; }
            let teid = tm.tedges[k].id;
            tm.collapse_tedge_snap_dedup(teid);
        }
        ck(check_emesh(&name, &tm, &g, "17_tm_dedup", None));
        repair_crossing_tedges(&mut tm);
        ck(check_emesh(&name, &tm, &g, "18_tm_repaired", None));
    }
}
