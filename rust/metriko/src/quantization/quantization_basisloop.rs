//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/quantization/quantization_basisloop.h (see rust/PORTED_FROM).

use crate::common::stdcxx::PriorityQueue;
use crate::tmesh::emesh::*;

#[derive(Clone, Copy, Debug)]
pub struct Comparator {
    pub length: i32,
    pub weight: f64,
    pub thid: i32,
}

impl PartialEq for Comparator {
    fn eq(&self, rhs: &Self) -> bool { self.length == rhs.length && self.weight == rhs.weight }
}

impl Comparator {
    // flip weight comparator
    pub fn lt(&self, rhs: &Self) -> bool { self.length < rhs.length || (self.length == rhs.length && self.weight > rhs.weight) }
    // same
    pub fn gt(&self, rhs: &Self) -> bool { self.length > rhs.length || (self.length == rhs.length && self.weight < rhs.weight) }
    pub fn le(&self, rhs: &Self) -> bool { !self.gt(rhs) }
    pub fn ge(&self, rhs: &Self) -> bool { !self.lt(rhs) }
}

// side of every thalf within its tquad. emesh keeps it only inside the tquad data, so the quantization,
// which looks it up for every step of the search, builds the table once
pub fn thalf_sides(tm: &Emesh<'_>) -> Vec<i32> {
    let mut sides = vec![-1i32; tm.thalfs.len()];
    for tq in tm.live_tquads() {
        for d in &tq.data { sides[d.thid as usize] = d.side; }
    }
    sides
}

pub fn gen_basis_loop<F: Fn(&Comparator, &Comparator) -> bool>(
    tm: &Emesh<'_>,
    th2side: &[i32],
    r: &[f64],
    bgn: &Ehalf,
    compare: F,
) -> Vec<i32> {
    let thalfs = &tm.thalfs;
    let mut q = PriorityQueue::new(|a: &Comparator, b: &Comparator| compare(a, b));
    let mut parents = vec![-1i32; thalfs.len()]; // parent thalf in the search tree
    let mut visited = vec![false; thalfs.len()]; // pushed once already

    let mut counter = 0;
    q.push(Comparator { length: 0, weight: 0., thid: bgn.id });

    loop {
        let top = *q.top().unwrap();
        let len = top.length;
        let prev = thalfs[top.thid as usize];
        q.pop();

        for thid in prev.tquad(tm).thids((th2side[prev.id as usize] + 2) % 4) {
            let curr = *thalfs[thid as usize].twin(tm);

            if curr.id == bgn.id && counter > 0 {
                let mut idx = prev.id;
                let mut res = vec![idx];
                loop {
                    idx = parents[idx as usize];
                    res.push(idx);
                    if idx == bgn.id { break; }
                }
                return res;
            }

            // every ancestor of prev except bgn has been pushed, and bgn is
            // handled above, so the pushed flag covers the ancestor test too
            if !visited[curr.id as usize] {
                q.push(Comparator { length: len - 1, weight: r[curr.teid as usize], thid: curr.id });
                parents[curr.id as usize] = prev.id;
                visited[curr.id as usize] = true;
            }
        }
        counter += 1;
        if q.is_empty() { break; }
    }
    metriko_fail!("no loop through thalf {}", bgn.id);
}
