//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/qex/gen_q_vert.h (see rust/PORTED_FROM).

use crate::common::typedef::*;
use crate::common::utilities::get_minmax_int;
use crate::common::stdcomplex::std_div;
use crate::common::predicates::{is_collinear, is_inside_triangle};
use crate::hmesh::hmesh::Hmesh;
use crate::qex::common::*;

/// `generate_q_vert(mesh, vqvs, eqvs, fqvs)`: returns `(vqvs, eqvs, fqvs)`
pub fn generate_q_vert(mesh: &Hmesh) -> (Vec<Qvert>, Vec<Qvert>, Vec<Qvert>) {
    let mut vqvs = vec![];
    let mut eqvs = vec![];
    let mut fqvs = vec![];

    // vert_q_vert
    for v in mesh.verts() {
        let uv = v.half().next().crnr().uv(); // checks only one corner
        let x = uv.re.abs() % 1.;
        let y = uv.im.abs() % 1.;
        if (x < EPS || 1. - x < EPS) && (y < EPS || 1. - y < EPS) { vqvs.push(Qvert::new(cx(x, y), v.pos(), v.id)); }
    }

    // edge_q_vert
    for e in mesh.edges() {
        let uv1 = e.half().cr_t().uv();
        let uv2 = e.half().cr_h().uv();
        let mm = get_minmax_int(&[uv1, uv2]);
        for x in mm.min_x..=mm.max_x {
            for y in mm.min_y..=mm.max_y {
                let xy = cx(x as f64, y as f64);
                let a = std_div(xy - uv1, uv2 - uv1).re; // imag part is evaluated by is_collinear
                if is_collinear(uv1, uv2, xy, EPS) && a > EPS && a < 1. - EPS { eqvs.push(Qvert::new(xy, e.lerp(a), e.id)); }
            }
        }
    }

    // face_q_vert
    for f in mesh.faces() {
        let [c1, c2, c3] = f.crnrs();
        let mm = get_minmax_int(&[c1.uv(), c2.uv(), c3.uv()]);
        for x in mm.min_x..=mm.max_x {
            for y in mm.min_y..=mm.max_y {
                let c = cx(x as f64, y as f64);
                if is_inside_triangle(c1.uv(), c2.uv(), c3.uv(), c) { fqvs.push(Qvert::new(c, f.uv2pos(c), f.id)); }
            }
        }
    }
    (vqvs, eqvs, fqvs)
}
