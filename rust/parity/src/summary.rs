//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of metriko_web/wasm/cli.cpp (summary lines and hash) and quadpatch.h (patch_boundaries, npatch).

use glam::DVec3;
use std::collections::{BTreeMap, BTreeSet};

/// what metriko_web's cli prints for a run (and tools/dump writes to summary.json)
#[derive(Clone, Debug, PartialEq)]
pub struct Summary {
    pub n_input_vertices: usize,
    pub n_input_triangles: usize,
    pub n_quad_vertices: usize,
    pub n_quads: usize,
    pub n_patches: usize,
    pub n_boundary_edges: usize,
    pub pos_sum: f64,
    pub hash: u32,
}

impl Summary {
    /// `pos`: the quad vertices (already multiplied by the unit), `idx`: the quads, `tqid`: patch id per quad
    pub fn new(n_input_vertices: usize, n_input_triangles: usize, pos: &[DVec3], idx: &[[i32; 4]], tqid: &[i32]) -> Summary {
        let bnd = patch_boundaries(idx, tqid);
        let ps: BTreeSet<i32> = tqid.iter().copied().filter(|&t| t >= 0).collect();
        let mut sum = 0.;
        for p in pos { sum += p.x; sum += p.y; sum += p.z; }
        Summary {
            n_input_vertices,
            n_input_triangles,
            n_quad_vertices: pos.len(),
            n_quads: idx.len(),
            n_patches: ps.len(),
            n_boundary_edges: bnd.len() / 2,
            pos_sum: sum,
            hash: quad_hash(idx, tqid),
        }
    }

    /// the four lines of cli.cpp before "time", byte for byte
    pub fn lines(&self) -> [String; 4] {
        [
            format!("input   {} vertices, {} triangles", self.n_input_vertices, self.n_input_triangles),
            format!("quads   {} vertices, {} quads, {} patches, {} boundary edges", self.n_quad_vertices, self.n_quads, self.n_patches, self.n_boundary_edges),
            format!("pos sum {}", fmt_e16(self.pos_sum)),
            format!("hash    {:08x}", self.hash),
        ]
    }
}

/// fnv-1a over the quad corners and then the patch ids, each as 4 little endian bytes (cli.cpp `mix`)
pub fn quad_hash(idx: &[[i32; 4]], tqid: &[i32]) -> u32 {
    let mut h: u32 = 2166136261;
    let mut mix = |x: u32| for b in 0..4 { h ^= (x >> (8 * b)) & 0xff; h = h.wrapping_mul(16777619); };
    for q in idx { for &i in q { mix(i as u32); } }
    for &t in tqid { mix(t as u32); }
    h
}

/// quadpatch.h `patch_boundaries`: quad edges between two patches, as vertex pairs (min, max) in edge order
pub fn patch_boundaries(qidx: &[[i32; 4]], tqid: &[i32]) -> Vec<i32> {
    let mut edge_quad: BTreeMap<(i32, i32), i32> = BTreeMap::new();
    for (q, quad) in qidx.iter().enumerate() {
        for j in 0..4 {
            let (a, b) = (quad[j], quad[(j + 1) % 4]);
            let key = (a.min(b), a.max(b));
            match edge_quad.get_mut(&key) {
                None => { edge_quad.insert(key, q as i32); }
                Some(v) => if *v >= 0 && tqid[q] == tqid[*v as usize] { *v = -1; },
            }
        }
    }
    edge_quad.iter().filter(|(_, q)| **q >= 0).flat_map(|(e, _)| [e.0, e.1]).collect()
}

/// std::format("{:.16e}", x): 17 significant digits, exponent with sign and at least two digits (1.5e+03)
pub fn fmt_e16(x: f64) -> String {
    if !x.is_finite() { return if x.is_nan() { "nan".into() } else if x > 0. { "inf".into() } else { "-inf".into() }; }
    let s = format!("{x:.16e}");
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{m}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
}
