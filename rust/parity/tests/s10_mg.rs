//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Stages 10, 11: the motorcycle graph and the T-mesh built from it. Inputs from the dumps: the mesh (01_hm), its
//! corner uvs after the IGM (09_hm_cfn), and the combed field's matching / singular (05_cmbf). Everything is plain
//! arithmetic on those inputs, so the comparison is exact (structure) and bitwise (values).

use metriko::tmesh::emesh::Emesh;
use metriko::tmesh::motorcycle::Mgrph;
use parity::*;
use parity::stage::*;

#[test]
fn motorcycle_graph_and_tmesh() {
    for (c, g) in golden_cases() {
        let mut hm = hmesh(&g, "01_hm", true);
        hm.cfn = g.load("09_hm_cfn.cfn").c64s().to_vec();
        let matching = g.load("05_cmbf.matching").i32s().to_vec();
        let singular = g.load("05_cmbf.singular").i32s().to_vec();
        let mg = Mgrph::new(&hm, &matching, &singular);
        check_mgrph(&c.name(), &mg, &g, "10_mg").unwrap_or_else(|e| panic!("{e}"));
        let tm = Emesh::new(&mg);
        check_emesh(&c.name(), &tm, &g, "11_tm_init", None).unwrap_or_else(|e| panic!("{e}"));
    }
}
