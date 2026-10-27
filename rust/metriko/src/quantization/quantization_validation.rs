//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/quantization/quantization_validation.h (see rust/PORTED_FROM).

use std::collections::VecDeque;
use crate::tmesh::emesh::*;
use crate::tmesh::motorcycle::*;

pub fn compute_validation(mg: &Mgrph<'_>, tm: &Emesh<'_>, x: &[f64]) -> bool {
    if x.iter().any(|&v| v < 0.) { return false; }

    let mut adjs: Vec<Vec<i32>> = vec![vec![]; mg.mnodes.len()];
    for e in tm.live_tedges() {
        if x[e.id as usize] == 0. {
            let (f, b) = (*e.nids.first().unwrap(), *e.nids.last().unwrap());
            adjs[f as usize].push(b);
            adjs[b as usize].push(f);
        }
    }

    let mut visited = vec![false; mg.mnodes.len()];

    for i in 0..mg.mnodes.len() {
        if mg.mnodes[i].jt != JunctionType::F || visited[i] { continue; }

        let start_nid = i as i32;
        let mut q: VecDeque<i32> = VecDeque::new();
        q.push_back(start_nid);
        visited[start_nid as usize] = true;

        while let Some(curr_nid) = q.pop_front() {
            if curr_nid != start_nid && mg.mnodes[curr_nid as usize].jt == JunctionType::F { return false; }

            for &next_nid in &adjs[curr_nid as usize] {
                if !visited[next_nid as usize] {
                    visited[next_nid as usize] = true;
                    q.push_back(next_nid);
                }
            }
        }
    }

    true
}
