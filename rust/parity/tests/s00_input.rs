//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Stage 00: the input as the C++ reads it (00_input.*, 01_hm.pos), and the cli summary of the C++ quads.

use parity::*;

/// the Rust obj loader, grid_unit and V / unit give the dumped V, F, unit and scaled positions bit for bit
#[test]
fn input_mesh() {
    for (c, g) in golden_cases() {
        let m = load_mesh(&c.mesh_path()).unwrap_or_else(|e| panic!("{}: {e}", c.name()));
        let name = |s: &str| format!("{} {s}", c.name());

        assert_bits_f64(&name("V"), &flat3(&m.v), g.load("00_input.V").f64s(), 3);
        assert_eq_i32(&name("F"), &flat(&m.f), g.load("00_input.F").i32s(), 3);
        assert_bits_f64(&name("scale"), &[c.scale], g.load("00_input.scale").f64s(), 1);

        let unit = grid_unit(&m.v, c.scale);
        assert_bits_f64(&name("unit"), &[unit], g.load("00_input.unit").f64s(), 1);
        assert_eq!(unit.to_bits(), g.meta_f64("00_input.unit").to_bits(), "{}: unit in manifest.json", c.name());

        // lib.h: Hmesh(MatXd(V / unit), F)
        let scaled: Vec<glam::DVec3> = m.v.iter().map(|p| *p / unit).collect();
        assert_bits_f64(&name("V / unit"), &flat3(&scaled), g.load("01_hm.pos").f64s(), 3);
        assert_eq_i32(&name("hm idx"), &flat(&m.f), g.load("01_hm.idx").i32s(), 3);

        assert_eq!(m.v.len() as u64, g.summary["n_input_vertices"].as_u64().unwrap(), "{}", c.name());
        assert_eq!(m.f.len() as u64, g.summary["n_input_triangles"].as_u64().unwrap(), "{}", c.name());
    }
}

/// parity::Summary on the dumped C++ quads prints the lines tools/dump (i.e. metriko_web's cli) printed
#[test]
fn cli_summary_of_cpp_quads() {
    for (c, g) in golden_cases() {
        let pos = g.load("29_quads.q_pos").dvec3s();
        let idx = g.load("29_quads.q_idx").rows_i32::<4>();
        let val = g.load("29_quads.q_val").i32s().to_vec();
        let s = Summary::new(
            g.summary["n_input_vertices"].as_u64().unwrap() as usize,
            g.summary["n_input_triangles"].as_u64().unwrap() as usize,
            &pos, &idx, &val,
        );
        let expected: Vec<&str> = g.summary["lines"].as_array().unwrap().iter().map(|l| l.as_str().unwrap()).collect();
        assert_eq!(s.lines().to_vec(), expected, "{}", c.name());
        assert_eq!(format!("{:08x}", s.hash), g.summary["hash"].as_str().unwrap(), "{}", c.name());
    }
}
