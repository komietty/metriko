//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/igm/iterative_rounding/iter_rounding_loop.h (see rust/PORTED_FROM).

use faer::Col;
use faer::sparse::Triplet;
use crate::solver::matrix_ops::*;
use crate::solver::levenberg_marquardt::{SolverTraits, spmv, vsub};
use crate::igl::speye::speye;
use crate::igm::iterative_rounding::iter_rounding_common::*;
use crate::igm::iterative_rounding::iter_rounding_init::NaiveIntegration;

/// `SeamlessIntegration<T>`, built from the naive integration it continues (`T = NaiveIntegration`).
pub struct SeamlessIntegration<'a> {
    pub base: Integration<'a>,
    pub ufull: &'a SprsD,
    pub frac: f64,
    pub x_init: Col<f64>,
    pub x_curr: Col<f64>,
    pub x_prev: Col<f64>,
    pub left_idcs: Vec<i32>,
    pub fixed_idcs: Vec<i32>,
    pub fixed_vals: Vec<i32>,
    pub integer_idcs: Vec<i32>,  // used when rounding uv on seams
    pub singular_idcs: Vec<i32>, // used when rounding uv on singularities
    pub rounded_singulars: bool,
    pub round_seams: bool,
}

impl<'a> SeamlessIntegration<'a> {
    pub fn init_fixed_indices(&mut self) -> bool {
        self.x_prev = self.x_curr.clone();
        let x = spmv(self.ufull, &self.x_curr);
        let mut min_dif = 3276700.0;
        let mut min_idx: i32 = -1;
        for i in 0..self.left_idcs.len() {
            let v = self.frac * x[self.left_idcs[i] as usize];
            let d = (v - v.round()).abs();
            if d < min_dif { min_idx = i as i32; min_dif = d; }
        }

        let mi = min_idx as usize;
        self.fixed_idcs.push(self.left_idcs[mi]);
        self.fixed_vals.push(x[self.left_idcs[mi] as usize].round() as i32);
        self.left_idcs.remove(mi);

        // completed rounding singularities; starting to round rest of seams
        if self.left_idcs.is_empty() && !self.round_seams && !self.rounded_singulars {
            self.left_idcs = self.integer_idcs.clone();
            self.rounded_singulars = true;
        }

        //Poisson error
        self.base.g_integ = scale(&self.base.x2f, self.base.w_integ);

        //Closeness
        self.base.g_close = speye(self.x_curr.nrows());
        self.base.g_close = scale(&self.base.g_close, self.base.w_close);

        //fixedIndices constantness
        let mut gc_t: Vec<TripD> = vec![];
        for i in 0..self.fixed_idcs.len() {
            gc_t.push(Triplet::new(i, self.fixed_idcs[i] as usize, 1.));
        }

        self.base.g_const = from_triplets(self.fixed_idcs.len(), x.nrows(), &gc_t);
        self.base.g_const = scale(&matmul(self.base.g_const.as_ref(), self.ufull.as_ref()), self.base.w_const);

        if self.base.e_size == 0 { // why needed this inside code
            let mut j = empty_sprs();
            let mut e = Col::<f64>::zeros(0);
            // the C++ passes VecXd::Random(UFull.cols()): only the size of E is read, which does not depend on the
            // values (and the barrier update without the jacobian has no side effect), so zeros give the same ESize
            self.objective_jacobian_impl(&Col::<f64>::zeros(self.ufull.ncols()), &mut e, &mut j, false);
            self.base.e_size = e.nrows();
        }

        min_dif > 10e-7
    }

    // memo: J * E is rhs of LMSolver
    pub fn objective_jacobian_impl(&self, x_curr_small: &Col<f64>, e: &mut Col<f64>, j: &mut SprsD, update_j: bool) {
        let b = &self.base;
        let x = spmv(self.ufull, x_curr_small); // here UFull
        let field = spmv(&b.x2f, x_curr_small);
        let f_integ = vsub(&field, b.f2);
        let f_close = vsub(x_curr_small, &self.x_prev);
        let mut f_const = Col::<f64>::zeros(self.fixed_idcs.len());
        for i in 0..self.fixed_idcs.len() {
            f_const[i] = x[self.fixed_idcs[i] as usize] - self.fixed_vals[i] as f64 * self.frac;
        }
        b.compute_jacobian(&field, &f_integ, &f_close, &f_const, update_j, e, j, true);
    }

    // assure that all values of round diff by fixed indices are less than 1e-7
    pub fn post_checking(&mut self, x: &Col<f64>) -> bool {
        let f = spmv(self.ufull, x);
        let v: Vec<f64> = self.fixed_idcs.iter().map(|&i| self.frac * f[i as usize]).collect();
        self.x_curr = x.clone();
        self.x_init = x.clone();
        self.x_prev = x.clone();
        if self.fixed_idcs.is_empty() { return true; }
        // (v - v.round()).abs().maxCoeff()
        let mut m = (v[0] - v[0].round()).abs();
        for &vi in &v[1..] { let d = (vi - vi.round()).abs(); if d > m { m = d; } }
        m <= 10e-7
    }

    /// `SeamlessIntegration(ni, xInit_, F2_, iIdcs, sIdcs, roundSeams_)`
    pub fn new(
        ni: &'a NaiveIntegration<'a>,
        x_init: &Col<f64>,
        f2: &'a Col<f64>,
        i_idcs: &[i32],
        s_idcs: &[i32],
        round_seams: bool,
    ) -> Self {
        let nb = &ni.base;
        let base = Integration::new(nb.g2, f2, nb.n_big, nb.n, nb.local_injectivity, nb.length, 1., 10e4, 0.01, nb.w_barrier, nb.barrier);
        let integer_idcs = i_idcs.to_vec();
        let singular_idcs = s_idcs.to_vec();
        let left_idcs = if round_seams { integer_idcs.clone() } else { singular_idcs.clone() };
        let mut si = SeamlessIntegration {
            base,
            ufull: &ni.ufull,
            frac: 1.,
            x_init: x_init.clone(),
            x_curr: x_init.clone(),
            x_prev: x_init.clone(),
            left_idcs,
            fixed_idcs: vec![],
            fixed_vals: vec![],
            integer_idcs,
            singular_idcs,
            rounded_singulars: false,
            round_seams,
        };
        si.base.x_size = si.ufull.ncols();
        si.base.x2f = pruned(scale(&matmul(si.base.g2.as_ref(), si.ufull.as_ref()), si.base.length).as_ref());
        si
    }
}

impl SolverTraits for SeamlessIntegration<'_> {
    fn x_size(&self) -> usize { self.base.x_size }
    fn e_size(&self) -> usize { self.base.e_size }
    fn initial_solution(&self, x0_: &mut Col<f64>) { *x0_ = self.x_init.clone(); }
    fn pre_iteration(&mut self, prevx: &Col<f64>) { Integration::pre_iteration(prevx); }
    fn objective_jacobian(&mut self, x: &Col<f64>, e: &mut Col<f64>, j: &mut SprsD, update_j: bool) {
        self.objective_jacobian_impl(x, e, j, update_j);
    }
    fn post_iteration(&mut self, _x: &Col<f64>) -> bool { false }
}
