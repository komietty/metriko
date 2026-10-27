//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/quantization/quantization.h (see rust/PORTED_FROM).

use crate::common::stdcxx;
use crate::quantization::quantization_basisloop::*;
use crate::quantization::quantization_constraint::*;
use crate::quantization::quantization_evaluation::*;
use crate::quantization::quantization_validation::*;
use crate::quantization::quantization_weight::*;
use crate::solver::matrix_ops::*;
use crate::tmesh::emesh::*;
use crate::tmesh::motorcycle::*;

/// Eigen 5.0.1's `v.squaredNorm()` of a dynamic-size expression with element `i` = `coeff(i)`: the redux of
/// `cwiseAbs2()` in `LinearVectorizedTraversal` on 2-lane double packets (NEON, the reference build): two packet
/// accumulators over 4 elements per step, their sum, one more packet if any, the horizontal add, then the scalar
/// tail. The order of the additions is Eigen's, so the energies compared with `<` / `<=` below match bitwise.
fn eigen_squared_norm(n: usize, coeff: impl Fn(usize) -> f64) -> f64 {
    const PS: usize = 2;
    let sq = |i: usize| { let c = coeff(i); c * c };
    let aligned_size2 = (n / (2 * PS)) * (2 * PS);
    let aligned_size = (n / PS) * PS;
    if aligned_size == 0 {
        if n == 0 { return 0.; }
        let mut res = sq(0);
        for i in 1..n { res += sq(i); }
        return res;
    }
    let mut p0 = [sq(0), sq(1)];
    if aligned_size > PS {
        let mut p1 = [sq(2), sq(3)];
        let mut i = 2 * PS;
        while i < aligned_size2 {
            p0 = [p0[0] + sq(i), p0[1] + sq(i + 1)];
            p1 = [p1[0] + sq(i + 2), p1[1] + sq(i + 3)];
            i += 2 * PS;
        }
        p0 = [p0[0] + p1[0], p0[1] + p1[1]];
        if aligned_size > aligned_size2 { p0 = [p0[0] + sq(aligned_size2), p0[1] + sq(aligned_size2 + 1)]; }
    }
    let mut res = p0[0] + p0[1];
    for i in aligned_size..n { res += sq(i); }
    res
}

/// `(x.cwiseQuotient(R) - I).norm()`
fn energy_norm(x: &[f64], r: &[f64]) -> f64 {
    eigen_squared_norm(x.len(), |j| x[j] / r[j] - 1.).sqrt()
}

/// `(C * X).norm()`: integer entries, exact in any order
fn constraint_residual(c: &faer::Mat<f64>, x: &[f64]) -> f64 {
    let mut y = vec![0.; c.nrows()];
    for j in 0..c.ncols() {
        let xj = x[j];
        for i in 0..c.nrows() { y[i] += c[(i, j)] * xj; }
    }
    eigen_squared_norm(y.len(), |i| y[i]).sqrt()
}

/// `X += s * G.col(k)`
fn add_col(x: &mut [f64], g: &SprsD, k: usize, s: f64) {
    let rows = g.row_idx_of_col_raw(k);
    let vals = g.val_of_col(k);
    for (&i, &v) in rows.iter().zip(vals.iter()) { x[i] += s * v; }
}

pub fn validate_quantization(tm: &Emesh<'_>, x: &[f64]) {
    for tquad in tm.live_tquads() {
        for i in 0..2 {
            let mut sum_a = 0i32;
            let mut sum_b = 0i32;
            for thid in tquad.thids(i)     { sum_a += x[tm.thalfs[thid as usize].teid as usize] as i32; }
            for thid in tquad.thids(i + 2) { sum_b += x[tm.thalfs[thid as usize].teid as usize] as i32; }
            metriko_check!(sum_a == sum_b, "tquad {} side {}: {} vs {}", tquad.id, i, sum_a, sum_b);
        }
    }
}

// an arc of the strip graph: from thalf `fr` (entering its tquad) to thalf
// `to` = twin of a thalf on the opposite side. traversing the arc means the
// loop passes the tedge of `to`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StripArc { pub fr: i32, pub to: i32, pub teid: i32 }

pub fn build_strip_arcs(tm: &Emesh<'_>, th2side: &[i32]) -> Vec<StripArc> {
    let mut arcs = vec![];
    for th in &tm.thalfs {
        for d in &th.tquad(tm).data {
            if d.side == (th2side[th.id as usize] + 2) % 4 {
                let to = tm.thalfs[d.thid as usize].twin(tm);
                arcs.push(StripArc { fr: th.id, to: to.id, teid: to.teid });
            }
        }
    }
    arcs
}

// pricing: find one negative-cost cycle in the strip graph (bellman-ford
// with a virtual super source), or return empty when none exists. cycle
// cost = sum of w[teid] over its nodes, i.e. the exact energy change of
// applying +-1 along the loop (up to tedges hit twice, re-checked by the
// caller).
pub fn find_negative_cycle(n_thalfs: i32, arcs: &[StripArc], w: &[f64]) -> Vec<i32> {
    let n = n_thalfs as usize;
    let mut dist = vec![0.; n];
    let mut pred = vec![-1i32; n];
    let mut mark = vec![0i32; n];

    // any cycle of the predecessor graph is a negative cycle (every
    // relaxation strictly decreased a label), so a walk along pred that
    // returns to a node of the same walk yields one without waiting for
    // the n passes of bellman-ford to finish.
    let cycle_in_pred = |pred: &[i32], mark: &mut [i32]| -> Vec<i32> {
        mark.fill(0);
        for s in 0..n_thalfs {
            if mark[s as usize] != 0 { continue; }
            let mut x = s;
            while x != -1 && mark[x as usize] == 0 { mark[x as usize] = s + 1; x = pred[x as usize]; }
            if x == -1 || mark[x as usize] != s + 1 { continue; } // ended at a root or in an older walk
            let mut cyc = vec![x];
            let mut v = pred[x as usize];
            while v != x { cyc.push(v); v = pred[v as usize]; }
            return cyc;
        }
        vec![]
    };

    const CHECK_EVERY: i32 = 8;
    for it in 0..n_thalfs {
        let mut relaxed = false;
        for a in arcs {
            let (fr, to, teid) = (a.fr as usize, a.to as usize, a.teid as usize);
            if dist[fr] + w[teid] < dist[to] - 1e-12 {
                dist[to] = dist[fr] + w[teid];
                pred[to] = a.fr;
                relaxed = true;
            }
        }
        if !relaxed { return vec![]; } // converged: no negative cycle
        if it % CHECK_EVERY == CHECK_EVERY - 1 {
            let cyc = cycle_in_pred(&pred, &mut mark);
            if !cyc.is_empty() { return cyc; }
        }
    }
    cycle_in_pred(&pred, &mut mark) // still relaxing after n passes: the pred graph must contain a cycle
}

pub fn compute_quantization(tm: &Emesh<'_>, mg: &Mgrph<'_>) -> Vec<f64> {
    let nte = tm.tedges.len() as i32;
    let n = nte as usize;
    let mut r = vec![0.; n];
    for th in &tm.thalfs { r[th.teid as usize] = th.r; } // both thalfs of a tedge carry its (uv) length
    let th2side = thalf_sides(tm);

    let mut x = vec![0.; n];
    let c = compute_constraint(tm);
    let max1e9 = |w: f64| if w < 1e-9 { 1e-9 } else { w }; // std::max(w, 1e-9)
    let g = construct_generating_vectors(tm, &th2side, &r, |c1: &Comparator, c2: &Comparator| {
        c1.length as f64 / max1e9(c1.weight) < c2.length as f64 / max1e9(c2.weight)
    });
    let teid2clid = cols_by_row(&g); // teid -> columns of G whose loop passes it

    // ----- construct first step vector ----- //
    while x.contains(&0.) {
        let mut min = 1e+9;
        let mut thid = 0;
        for th in tm.thalfs.iter().filter(|th_| x[th_.teid as usize] == 0.) {
            let w = compute_weight(r[th.teid as usize], x[th.teid as usize], nte, false);
            if w < min { min = w; thid = th.id; }
        }

        let col_idcs = &teid2clid[tm.thalfs[thid as usize].teid as usize];
        if !col_idcs.is_empty() { add_col(&mut x, &g, col_idcs[0] as usize, 1.); }
        // todo: latter code might be better to reflect my intention...
        //const auto& ks = teid2clid[tm.thalfs[thid].teid];
        //auto k = rg::min_element(ks, {}, [&](int k) { return G.col(k).nonZeros(); });
        //if (k != ks.end()) X += G.col(*k);

        if compute_validation(mg, tm, &x) { break; }
    }

    metriko_check!(constraint_residual(&c, &x) < 1e-12, "first step breaks the quantization constraint");

    // ----- construct second step vector ----- //

    let mut e = energy_norm(&x, &r);

    // ----- second step: exact pricing (column generation) ----- //
    // apply negative-cost loops until none exists. the pricing searches ALL
    // strip loops, not a precomputed basis, so the result is locally optimal
    // w.r.t. every single-loop +-1 move.
    let arcs = build_strip_arcs(tm, &th2side);

    let mut improved = true;
    let mut guard = 0;
    while improved && guard < 100 * nte {
        improved = false;
        for sgn in [1i32, -1] {
            // marginal cost of changing tedge j by sgn:
            // ((X_j+sgn)/R_j - 1)^2 - (X_j/R_j - 1)^2
            let mut w = vec![0.; n];
            for j in 0..n {
                w[j] = (2. * sgn as f64 * (x[j] - r[j]) + 1.) / (r[j] * r[j]);
                if sgn < 0 && x[j] < 1. { w[j] = 1e18; } // keep X >= 0
            }
            let cyc = find_negative_cycle(tm.thalfs.len() as i32, &arcs, &w);
            if cyc.is_empty() { continue; }

            let mut gv = vec![0.; n];
            for &thid in &cyc { gv[tm.thalfs[thid as usize].teid as usize] += 1.; }
            let x1: Vec<f64> = (0..n).map(|j| x[j] + sgn as f64 * gv[j]).collect();
            let e1 = energy_norm(&x1, &r); // exact: a loop can hit one tedge twice
            if e1 < e && x1.iter().all(|&v| v >= 0.) && compute_validation(mg, tm, &x1) {
                x = x1;
                e = e1;
                improved = true;
            }
        }
        guard += 1;
    }

    let mut counter = 0;
    while counter < 30 {
        let prev_e = e;
        // `vec<std::tuple<int, int, bool>>`: the weight is truncated to int
        let mut es: Vec<(i32, i32, bool)> = vec![];
        let l = nte;
        for j in 0..l {
            es.push((compute_weight(r[j as usize], x[j as usize], l, false) as i32, j, false));
            es.push((compute_weight(r[j as usize], x[j as usize], l, true) as i32, j, true));
        }
        stdcxx::sort_by(&mut es, |a, b| a.0 < b.0);

        'exit_loops: for j in 0..(l * 2) as usize {
            let ei = es[j].1;
            for &k in &teid2clid[ei as usize] {
                let mut x1 = x.clone();
                let mut x2 = x.clone();
                add_col(&mut x1, &g, k as usize, 1.);
                add_col(&mut x2, &g, k as usize, -1.);
                let n1 = energy_norm(&x1, &r);
                let n2 = energy_norm(&x2, &r);
                if n1 <= e && compute_validation(mg, tm, &x1) { x = x1; e = n1; break 'exit_loops; }
                if n2 <= e && compute_validation(mg, tm, &x2) { x = x2; e = n2; break 'exit_loops; }
            }
        }
        counter = if prev_e == e { counter + 1 } else { 0 };
    }

    metriko_check!(constraint_residual(&c, &x) < 1e-12, "result breaks the quantization constraint");
    x
}

#[cfg(test)]
mod tests {
    use super::eigen_squared_norm;

    #[test]
    fn eigen_redux_order() {
        // n = 7: p0 = [x0, x1], p1 = [x2, x3], p0 + p1, + [x4, x5], horizontal add, + x6
        let v = [1e16, 1., -1e16, 1., 1., 1., 3.];
        let got = eigen_squared_norm(v.len(), |i| v[i]);
        let sq: Vec<f64> = v.iter().map(|x| x * x).collect();
        let p = [(sq[0] + sq[2]) + sq[4], (sq[1] + sq[3]) + sq[5]];
        assert_eq!(got.to_bits(), ((p[0] + p[1]) + sq[6]).to_bits());
    }
}
