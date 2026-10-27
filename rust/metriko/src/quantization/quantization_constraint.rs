//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/quantization/quantization_constraint.h (see rust/PORTED_FROM).

use faer::Mat;
use crate::tmesh::emesh::*;

pub fn compute_constraint(tm: &Emesh<'_>) -> Mat<f64> {
    let mut m = Mat::<f64>::zeros(tm.tquads.len() * 2, tm.tedges.len());
    for iq in 0..tm.tquads.len() {
        let tq = &tm.tquads[iq];
        for ih in 0..tq.data.len() {
            let side = tq.data[ih].side;
            let teid = tm.thalfs[tq.data[ih].thid as usize].teid as usize;
            match side {
                0 => { m[(iq * 2, teid)] += 1.; }
                2 => { m[(iq * 2, teid)] -= 1.; }
                1 => { m[(iq * 2 + 1, teid)] += 1.; }
                3 => { m[(iq * 2 + 1, teid)] -= 1.; }
                _ => panic!("side has to be 0 ~ 3"),
            }
        }
    }
    m
}
