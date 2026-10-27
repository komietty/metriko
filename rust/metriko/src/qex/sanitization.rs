//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/qex/sanitization.h (see rust/PORTED_FROM).

use crate::common::typedef::*;
use crate::hmesh::hmesh::Hmesh;
use crate::qex::common::*;

pub fn sanitization(mesh: &mut Hmesh, matching: &[i32], singular: &[i32], rosy_n: i32) {
    let (he_r, he_t) = compute_trs_matrix(mesh, matching, rosy_n);
    // the C++ edits `mesh.cfn` through a reference while walking the mesh; the walk reads topology only
    let mut cfn = std::mem::take(&mut mesh.cfn);
    {
        let mesh: &Hmesh = mesh;
        for v in mesh.verts() {
            let mut max = 0.;
            for cc in v.adj_crnrs() {
                let uv = cfn[cc.id as usize];
                // std::max(a, b): b when a < b
                if max < uv.re.abs() { max = uv.re.abs(); }
                if max < uv.im.abs() { max = uv.im.abs(); }
            }

            let mut init = false;
            for h in v.adj_halfs() {
                if !init {
                    if singular[v.id as usize] != 0 {
                        let cc = h.next().crnr();
                        let uv = cfn[cc.id as usize];
                        cfn[cc.id as usize] = cx(uv.re.round(), uv.im.round());
                    } else {
                        let cc = h.next().crnr();
                        let mut uv = cfn[cc.id as usize];
                        let delta = 2f64.powf(max.log2());
                        let sgn = |x: f64| ((0. < x) as i32 - (x < 0.) as i32) as f64;
                        let sign = cx(sgn(uv.re), sgn(uv.im));

                        uv = (uv + sign * delta) - sign * delta;
                        // every coordinate within tol of an integer is that integer: a vertex near a grid point
                        // becomes the point, otherwise the same point is found as an edge or face q-vertex in the
                        // neighbouring faces and never pairs. a vertex near a single isoline lands exactly on it,
                        // which the tracer passes by stepping around the vertex (pick_next_half)
                        let tol = 1e-4;
                        let g = nearby_grid(uv);
                        if (uv.re - g.re).abs() < tol { uv.re = g.re; }
                        if (uv.im - g.im).abs() < tol { uv.im = g.im; }
                        cfn[cc.id as usize] = uv;
                    }
                } else {
                    let cp = h.twin().prev().crnr();
                    let cc = h.next().crnr();
                    let r = he_r[h.twin().id as usize];
                    let t = he_t[h.twin().id as usize];
                    cfn[cc.id as usize] = r * cfn[cp.id as usize] + t;
                }
                init = true;
            }
        }
    }
    mesh.cfn = cfn;
}
