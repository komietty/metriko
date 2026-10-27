//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/igm/iterative_rounding/iter_rounding_init.h (see rust/PORTED_FROM).

use std::cell::RefCell;
use faer::Col;
use faer::sparse::Triplet;
use crate::solver::matrix_ops::*;
use crate::solver::sparse_solver::{SparseLlt, SparseLu, ComputationInfo};
use crate::solver::levenberg_marquardt::{SolverTraits, spmv, spmv_scaled, tspmv, vnorm, vadd, vsub};
use crate::igl::slice::slice;
use crate::igm::iterative_rounding::injective_barrier::InjectiveBarrier;
use crate::igm::iterative_rounding::iter_rounding_common::*;

pub struct NaiveIntegration<'a> {
    pub base: Integration<'a>,
    pub fixed_idcs: &'a [i32],
    pub fixed_vals: &'a [f64],
    pub ufull: SprsD,      // parmMat * URaw.
    pub uext: SprsD,       // direct sum of UFull and rawF2. Just pass through rawF2.
    pub xf_small: Col<f64>,
    pub x0: Col<f64>,
    pub iter_cb: &'a dyn Fn(&Col<f64>) -> bool,
}

fn head(v: &Col<f64>, n: usize) -> Col<f64> { Col::<f64>::from_fn(n, |i| v[i]) }
fn tail(v: &Col<f64>, n: usize) -> Col<f64> { let o = v.nrows() - n; Col::<f64>::from_fn(n, |i| v[o + i]) }

impl<'a> NaiveIntegration<'a> {
    pub fn post_iteration_impl(&self, xf_: &Col<f64>) -> bool { (self.iter_cb)(&head(&spmv(&self.uext, xf_), self.x0.nrows())) }

    pub fn prepare_jacobian_component(&mut self) {
        let b = &mut self.base;
        let l1 = b.f2.nrows();
        let l2 = self.uext.nrows();
        let l3 = self.fixed_idcs.len();

        let mut t: Vec<TripD> = vec![];

        for (r, c, v) in nonzeros(b.g2.as_ref()) { t.push(Triplet::new(r, c, -b.length * v)); }
        for i in 0..l1 { t.push(Triplet::new(i, b.g2.ncols() + i, 1.)); }
        b.g_integ = from_triplets(b.g2.nrows(), b.g2.ncols() + l1, &t);

        t.clear();
        for i in 0..l1 { t.push(Triplet::new(i, self.x0.nrows() + i, 1.)); }
        b.g_close = from_triplets(l1, l2, &t);

        t.clear();
        for i in 0..l3 { t.push(Triplet::new(i, self.fixed_idcs[i] as usize, 1.)); }
        b.g_const = from_triplets(l3, l2, &t);

        b.g_integ = scale(&matmul(b.g_integ.as_ref(), self.uext.as_ref()), b.w_integ);
        b.g_close = scale(&matmul(b.g_close.as_ref(), self.uext.as_ref()), b.w_close);
        b.g_const = scale(&matmul(b.g_const.as_ref(), self.uext.as_ref()), b.w_const);
    }

    pub fn objective_jacobian_impl(&self, xf_: &Col<f64>, e: &mut Col<f64>, j: &mut SprsD, update_j: bool) {
        let b = &self.base;
        let xf = spmv(&self.uext, xf_); // here UExt...
        let x = head(&xf, self.x0.nrows());
        let f = tail(&xf, b.f2.nrows());
        let f_integ = vsub(&f, &spmv_scaled(b.g2, b.length, &x));
        let f_close = vsub(&f, b.f2);
        //const VecXd fInteg = F2 - length * G2 * xCurr;
        //const VecXd fClose = VecXd::Zero(F2.size());
        let mut f_const = Col::<f64>::zeros(self.fixed_idcs.len());
        for i in 0..self.fixed_idcs.len() {
            f_const[i] = x[self.fixed_idcs[i] as usize] - self.fixed_vals[i];
        }
        b.compute_jacobian(&f, &f_integ, &f_close, &f_const, update_j, e, j, false);
    }

    pub fn new(
        g2: &'a SprsD,
        cfull: &SprsD,
        f2: &'a Col<f64>,
        fixed_idcs: &'a [i32],
        fixed_vals: &'a [f64],
        weight_matrix: &SprsD, // per face weight of the poisson energy, repeated over the 2N rows of the face
        length: f64,
        locinj: bool,
        n_big: i32,
        n: i32,
        iter_cb: &'a dyn Fn(&Col<f64>) -> bool,
        barrier: &'a RefCell<InjectiveBarrier>,
    ) -> Self {
        let mut base = Integration::new(g2, f2, n_big, n, locinj, length, 10e3, 10e3, 1., 1e-4, barrier);

        //----- Generating naive poisson solution -----
        // Compute poisson eq. so that the gradient of enegy function equal to zero.
        // Conceptually it computes uv to follow the given nvec with the constraint.
        // See eq. (6) in the report by Bommes(2012)
        // if isometricity is not important, then the solver below might have room for optimization
        // e.g. consider only the conformality... use conformal optimization

        let ufull = sparse_null_space(cfull);
        base.x2f = pruned(matmul(g2.as_ref(), ufull.as_ref()).as_ref());
        let x2f_t_w = matmul(transpose(base.x2f.as_ref()).as_ref(), weight_matrix.as_ref());
        let e_mat = scale(&matmul(x2f_t_w.as_ref(), base.x2f.as_ref()), length);
        let f = spmv(&x2f_t_w, f2);

        let const_mat = slice(&ufull, fixed_idcs, 1);

        let mut bm_t: Vec<TripD> = vec![];

        for (r, c, v) in nonzeros(e_mat.as_ref()) { bm_t.push(Triplet::new(r, c, v)); }
        for (r, c, v) in nonzeros(const_mat.as_ref()) {
            bm_t.push(Triplet::new(r + e_mat.nrows(), c, v));
            bm_t.push(Triplet::new(c, r + e_mat.nrows(), v));
        }

        // the fixed values are constraints C x = v on a few variables only: write x = x0 + Z y, with Z the identity
        // on the untouched variables and the (dense, small) kernel of C on the touched ones. E restricted to y is
        // symmetric positive definite and is solved by cholesky; the LU of the full KKT system is the fallback
        let fixed_vals_col = Col::<f64>::from_fn(fixed_vals.len(), |i| fixed_vals[i]);
        let (z, xp) = sparse_null_space_rhs(&const_mat, &fixed_vals_col);
        let ey = matmul(matmul(transpose(z.as_ref()).as_ref(), e_mat.as_ref()).as_ref(), z.as_ref());
        let rhs = tspmv(&z, &vsub(&f, &spmv(&e_mat, &xp)));
        let llt = SparseLlt::<f64>::from_matrix(ey.as_ref());
        let y = llt.solve(&rhs);
        // a (nearly) singular E passes the factorization but not the solve
        let all_finite = (0..y.nrows()).all(|i| y[i].is_finite());
        let rhs_norm = vnorm(&rhs);
        let xsmall = if llt.info() == ComputationInfo::Success && all_finite
            && vnorm(&vsub(&spmv(&ey, &y), &rhs)) <= 1e-8 * if rhs_norm < 1. { 1. } else { rhs_norm }
        {
            vadd(&xp, &spmv(&z, &y))
        } else {
            let nb = e_mat.nrows() + const_mat.nrows();
            let big_mat = from_triplets(nb, nb, &bm_t);
            let mut big_rhs = Col::<f64>::zeros(f.nrows() + fixed_vals.len());
            for i in 0..f.nrows() { big_rhs[i] = f[i]; }
            for i in 0..fixed_vals.len() { big_rhs[f.nrows() + i] = fixed_vals[i]; }
            let solver = SparseLu::from_matrix(big_mat.as_ref());
            if solver.info() != ComputationInfo::Success { metriko_fail!("initial Poisson solve (SparseLU) failed"); }
            let xsmall_full = solver.solve(&big_rhs);
            head(&xsmall_full, ufull.ncols())
        };

        let x0 = spmv(&ufull, &xsmall);
        let xf_small = vconcat(&[&xsmall, f2]);

        let mut ue_t: Vec<TripD> = vec![];
        for (r, c, v) in nonzeros(ufull.as_ref()) { ue_t.push(Triplet::new(r, c, v)); }

        for k in 0..f2.nrows() { ue_t.push(Triplet::new(ufull.nrows() + k, ufull.ncols() + k, 1.)); }

        let uext = from_triplets(ufull.nrows() + f2.nrows(), ufull.ncols() + f2.nrows(), &ue_t);

        let mut ni = NaiveIntegration { base, fixed_idcs, fixed_vals, ufull, uext, xf_small, x0, iter_cb };

        let mut e_ = Col::<f64>::zeros(0);
        let mut j_ = empty_sprs();
        let xf = ni.xf_small.clone();
        ni.objective_jacobian_impl(&xf, &mut e_, &mut j_, false);
        ni.base.e_size = e_.nrows();
        ni.base.x_size = ni.uext.ncols();

        ni.prepare_jacobian_component();
        ni
    }
}

impl SolverTraits for NaiveIntegration<'_> {
    fn x_size(&self) -> usize { self.base.x_size }
    fn e_size(&self) -> usize { self.base.e_size }
    fn initial_solution(&self, xf_: &mut Col<f64>) { *xf_ = self.xf_small.clone(); }
    fn pre_iteration(&mut self, prevx: &Col<f64>) { Integration::pre_iteration(prevx); }
    fn objective_jacobian(&mut self, xf_: &Col<f64>, e: &mut Col<f64>, j: &mut SprsD, update_j: bool) {
        self.objective_jacobian_impl(xf_, e, j, update_j);
    }
    fn post_iteration(&mut self, xf_: &Col<f64>) -> bool { self.post_iteration_impl(xf_) }
}
