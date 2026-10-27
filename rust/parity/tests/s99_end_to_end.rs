//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! The whole pipeline from the obj file, every stage chained on the Rust outputs (no dumped input), against the
//! C++ result: the quads and their patch ids exactly (hence the metriko_web cli hash), the quad vertex positions
//! within 1e-9 of the bounding box diagonal (they carry the last-bit differences of the faer solves). Also the
//! checks of test/test_remesh.cpp on the intermediate stages.

use metriko::compute_quadrangulation_impl;
use metriko::tmesh::emesh_validate::validate_no_crossing;
use parity::*;

#[test]
fn quadrangulation_matches_cpp() {
    for (c, g) in golden_cases() {
        let n = c.name();
        let m = load_mesh(&c.mesh_path()).unwrap_or_else(|e| panic!("{n}: {e}"));
        let (pos, idx, val) = compute_quadrangulation_impl(&m.v, &m.f, c.scale, false, true, |im| {
            // test/test_remesh.cpp
            assert_eq!(validate_no_crossing(im.emesh, "test"), 0, "{n}: tedge contacts");
            assert!(!im.q_faces.is_empty(), "{n}: no quad faces");
            assert!(im.q_ports.iter().all(|p| p.is_connected), "{n}: unconnected port");
            assert!(im.q_faces.iter().all(|f| f.qhalfs.len() == 4), "{n}: non-quad face");
            assert!(im.q_idx.len() == im.q_faces.len() && im.q_val.len() == im.q_idx.len(), "{n}: sizes");
            assert!(im.q_val.iter().all(|&v| v >= 0), "{n}: unlabelled quad");
            (im.q_pos.to_vec(), im.q_idx.to_vec(), im.q_val.to_vec())
        });

        assert_eq_i32(&format!("{n} q_idx"), &flat(&idx), g.load("29_quads.q_idx").i32s(), 4);
        assert_eq_i32(&format!("{n} q_val"), &val, g.load("29_quads.q_val").i32s(), 1);
        let diag = grid_unit(&m.v, 1.);
        assert_close_f64(&format!("{n} q_pos"), &flat3(&pos), g.load("29_quads.q_pos").f64s(), 3, Tol::new(1e-9 * diag, 0.));

        let s = Summary::new(m.v.len(), m.f.len(), &pos, &idx, &val);
        assert_eq!(format!("{:08x}", s.hash), g.summary["hash"].as_str().unwrap(), "{n}: cli hash");
        let expected: Vec<&str> = g.summary["lines"].as_array().unwrap().iter().map(|l| l.as_str().unwrap()).collect();
        assert_eq!(s.lines()[1], expected[1], "{n}: quad / patch / boundary counts");
    }
}
