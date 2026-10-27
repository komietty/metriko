//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/qex/gen_q_face.h (see rust/PORTED_FROM).

use crate::qex::common::*;

pub fn generate_q_faces(qps: &[Qport], qes: &[Qedge]) -> Vec<Qface> {
    let mut qfs: Vec<Qface> = vec![];
    let mut qhs: Vec<Qhalf> = vec![];
    for i in 0..qes.len() {
        qhs.push(Qhalf { qe: i, idx: (i * 2) as i32, cano: true });
        qhs.push(Qhalf { qe: i, idx: (i * 2 + 1) as i32, cano: false });
    }

    let mut po2qh = vec![-1i32; qps.len()];
    let mut visit = vec![false; qhs.len()];

    for qh in &qhs { po2qh[qps[qh.port1(qes)].idx as usize] = qh.idx; }

    for i in 0..qhs.len() {
        if visit[i] { continue; }
        let mut qfhs = vec![qhs[i]];
        visit[i] = true;
        let mut qp = qhs[i].port2(qes);
        for _ in 0..3 {
            let h = po2qh[qps[qps[qp].prev_id as usize].idx as usize];
            if h < 0 { break; }
            qfhs.push(qhs[h as usize]);
            visit[h as usize] = true;
            qp = qhs[h as usize].port2(qes);
        }
        if qfhs.len() == 4 {
            metriko_check!(qps[qps[qp].prev_id as usize].idx == qps[qhs[i].port1(qes)].idx, "q-face ring does not close");
            qfs.push(Qface { qhalfs: qfhs });
        }
    }
    qfs
}
