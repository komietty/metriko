//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/igm/iterative_rounding/iter_rounding.h (see rust/PORTED_FROM).

use std::cell::RefCell;
use faer::Col;
use crate::solver::matrix_ops::*;
use crate::solver::sparse_solver::{SparseLlt, ComputationInfo};
use crate::solver::levenberg_marquardt::*;
use crate::igm::iterative_rounding::injective_barrier::InjectiveBarrier;
use crate::igm::iterative_rounding::iter_rounding_init::NaiveIntegration;
use crate::igm::iterative_rounding::iter_rounding_loop::SeamlessIntegration;

/// `CholeskyWrapper` without SuiteSparse (the reference configuration): the normal equations `J' J` by the sparse
/// LLT, the symbolic analysis redone only when the pattern of `J' J` changes.
#[derive(Default)]
pub struct CholeskyWrapper {
    outer: Vec<usize>,
    inner: Vec<usize>,
    llt: SparseLlt<f64>,
}

impl CholeskyWrapper {
    fn same_pattern(&self, a: &SprsD) -> bool {
        // a product result is compressed (faer: no per-column nnz); compare the column pointers and row indices
        a.symbolic().col_nnz().is_none()
            && self.outer.len() == a.ncols() + 1
            && self.inner.len() == non_zeros(a.as_ref())
            && self.outer.as_slice() == a.symbolic().col_ptr()
            && self.inner.as_slice() == a.symbolic().row_idx()
    }
}

impl LinearSolver for CholeskyWrapper {
    fn factorize(&mut self, j: &SprsD) -> bool {
        let a = matmul(transpose(j.as_ref()).as_ref(), j.as_ref());
        if !self.same_pattern(&a) {
            self.llt.analyze_pattern(a.as_ref());
            if self.llt.info() != ComputationInfo::Success { self.outer.clear(); self.inner.clear(); return false; }
            self.outer = a.symbolic().col_ptr().to_vec();
            self.inner = a.symbolic().row_idx().to_vec();
        }
        self.llt.factorize(a.as_ref());
        self.llt.info() == ComputationInfo::Success
    }

    fn solve(&mut self, rhs: &Col<f64>, x: &mut Col<f64>) -> bool { *x = self.llt.solve(rhs); true }
}

pub fn iterative_rounding(
    fixed_idcs: &[i32],
    fixed_vals: &[f64],
    weight_matrix: &SprsD,
    singular_idcs: &[i32],
    integer_idcs: &[i32],
    length: f64,
    cfull: &SprsD,
    g2: &SprsD,
    n_big: i32,
    n: i32,
    n_f: i32,
    seamless: bool,
    round_seams: bool,
    local_injectivity: bool,
    verbose: bool,
    intrinsic_field: &Col<f64>,
    iter_cb: &dyn Fn(&Col<f64>) -> bool,
    fullx: &mut Col<f64>,
) -> bool {
    let mut llt_ni = CholeskyWrapper::default();
    let mut llt_si = CholeskyWrapper::default();

    let mut dd_ni = DiagonalDamping::new(if local_injectivity { 0.01 } else { 0. });
    let mut is_lm = LMSolver::new();

    let barrier = RefCell::new(InjectiveBarrier::new(n_big, n_f, intrinsic_field));

    let mut ni = NaiveIntegration::new(
        g2,
        cfull,
        intrinsic_field,
        fixed_idcs,
        fixed_vals,
        weight_matrix,
        length,
        local_injectivity,
        n_big,
        n,
        iter_cb,
        &barrier,
    );

    //----- Naive Solution -----//
    if !seamless {
        // Either it doesn't need loc_inj or it is already loc_inj (because post_iteration has loc_inj check)
        if !local_injectivity || ni.post_iteration_impl(&ni.xf_small) { *fullx = ni.x0.clone(); return true; }
    }

    is_lm.init_default(&ni);
    is_lm.solve(&mut llt_ni, &mut ni, &mut dd_ni, verbose);

    // todo need to assign minix!
    let ux = spmv(&ni.uext, &is_lm.x);
    *fullx = Col::<f64>::from_fn(ni.ufull.nrows(), |i| ux[i]);
    if !seamless { return true; }

    //----- Seamless Solution -----//
    let x_ = Col::<f64>::from_fn(ni.ufull.ncols(), |i| is_lm.x[i]);
    let nn = (2 * n_big * n_f) as usize;
    let f2_ = Col::<f64>::from_fn(nn, |i| ux[ux.nrows() - nn + i]); // x is changed from F2, but how is this important??
    let mut si = SeamlessIntegration::new(&ni, &x_, &f2_, integer_idcs, singular_idcs, round_seams);

    let mut success = true;
    let mut rounded = false;
    let mut dd_si = DiagonalDamping::new(if local_injectivity { 0.01 } else { 0. });
    let mut ir_lm = LMSolver::new();

    while !si.left_idcs.is_empty() {
        if !si.init_fixed_indices() { continue; }
        rounded = true;
        dd_si.curr_lambda = if local_injectivity { 0.01 } else { 0. };
        ir_lm.init(&si, 100, 1e-7, 1e-7);
        ir_lm.solve(&mut llt_si, &mut si, &mut dd_si, verbose);
        if !si.post_checking(&ir_lm.x) {
            success = false;
            break;
        }
    }

    *fullx = if rounded { spmv(&ni.ufull, &ir_lm.x) } else { ni.x0.clone() };
    success
}
