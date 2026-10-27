//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/qex/gen_q_port.h (see rust/PORTED_FROM).

use crate::common::typedef::*;
use crate::common::utilities::{dot, get_quater_rot, lerp};
use crate::common::predicates::{is_collinear, is_points_into, orientation};
use crate::common::stdcxx;
use crate::common::stdcomplex::std_div;
use crate::hmesh::hmesh::Hmesh;
use crate::qex::common::*;

/// rg::sort(ps, {}, first), then the ring of next / prev ids, appended to `qports`
fn push_ring(mut ps: Vec<(f64, Qport)>, qports: &mut Vec<Qport>) {
    stdcxx::sort_by(&mut ps, |a, b| a.0 < b.0);

    let s = ps.len() as i32;
    for i in 0..s { ps[i as usize].1.idx = qports.len() as i32 + i; }
    for i in 0..s {
        ps[i as usize].1.prev_id = ps[((i - 1 + s) % s) as usize].1.idx;
        ps[i as usize].1.next_id = ps[((i + 1) % s) as usize].1.idx;
        qports.push(ps[i as usize].1);
    }
}

pub fn generate_eqvert_qport(mesh: &Hmesh, eqverts: &[Qvert], qports: &mut Vec<Qport>) {
    for qv in eqverts {
        let mut ps: Vec<(f64, Qport)> = vec![];
        let e = mesh.edge(qv.sid);
        for h in [e.half(), e.half().twin()] {
            let c1 = h.next().crnr();
            let c2 = h.prev().crnr();
            let uv1 = c1.uv();
            let uv2 = c2.uv();
            let p1 = c1.pos();
            let p2 = c2.pos();
            let uv = lerp(uv1, uv2, (qv.pos - p1).length() / (p2 - p1).length());

            for i in 0..4 {
                let d = get_quater_rot(i);

                if if is_collinear(uv1, uv2, uv + d, EPS) {
                    !h.is_canonical()                     // if collinear, always non cano face will be skipped
                } else {
                    orientation(uv1, uv2, uv + d) < 0.     // if not collinear dir must be inside face otherwise skipped
                } { continue; }

                let mut k = std_div(d, uv2 - uv1).arg(); // angle from the edge: [0, pi] in this chart
                if k < 0. { k += TWO_PI; }           // exact anti-parallel may yield -pi
                ps.push(((if h.is_canonical() { 0. } else { PI }) + k, Qport::new(-1, -1, e.id, h.face().id, nearby_grid(uv), d, qv.pos)));
            }
        }

        metriko_check!(ps.len() == 4, "edge {}: {} ports at the q-vert", e.id, ps.len());
        push_ring(ps, qports);
    }
}

pub fn generate_fqvert_qport(mesh: &Hmesh, fqverts: &[Qvert], qports: &mut Vec<Qport>) {
    for qv in fqverts {
        let f = mesh.face(qv.sid);
        for i in 0..4 {
            qports.push(Qport::new(qports.len() as i32, -1, -1, f.id, qv.uv, get_quater_rot(i), qv.pos));
        }
        for i in 0..4i32 {
            let l = qports.len() as i32;
            qports[(l - i - 1) as usize].next_id = qports[(l - (i - 1 + 4) % 4 - 1) as usize].idx;
            qports[(l - i - 1) as usize].prev_id = qports[(l - (i + 1 + 4) % 4 - 1) as usize].idx;
        }
    }
}

pub fn generate_vqvert_qport(mesh: &Hmesh, vqverts: &[Qvert], qports: &mut Vec<Qport>) {
    for qv in vqverts {
        let mut ps: Vec<(f64, Qport)> = vec![];
        let mut acc = 0.;
        let v = mesh.vert(qv.sid);
        for c in v.adj_crnrs() {
            let uv1 = c.uv();
            let uv2 = c.next().uv();
            let uv3 = c.prev().uv();
            for i in 0..4 {
                let d = get_quater_rot(i);
                let f1 = is_points_into(uv1, uv2, uv3, uv1 + d, EPS);
                let f2 = is_collinear(uv1, uv2, uv1 + d, EPS);
                let f3 = dot(uv2 - uv1, d) > 0.;
                if f1 || (f2 && f3) {
                    ps.push((acc + std_div(d, uv2 - uv1).arg(), Qport::new(-1, v.id, -1, c.face().id, uv1, d, qv.pos)));
                }
            }
            acc += std_div(uv3 - uv1, uv2 - uv1).arg();
        }

        push_ring(ps, qports);
    }
}
