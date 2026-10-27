//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/igm/iterative_rounding/iter_rounding_common.h (see rust/PORTED_FROM).

use std::cell::RefCell;
use faer::Col;
use faer::sparse::Triplet;
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets, nonzeros, matmul, sparse_block};
use crate::igm::iterative_rounding::injective_barrier::InjectiveBarrier;

/// Eigen's `A * s` for a sparse `A`: every stored value times `s`, the pattern unchanged.
pub fn scale(a: &SprsD, s: f64) -> SprsD {
    let t: Vec<TripD> = nonzeros(a.as_ref()).map(|(i, j, v)| Triplet::new(i, j, v * s)).collect();
    from_triplets(a.nrows(), a.ncols(), &t)
}

/// `[a; b; c; ...]` of dense vectors
pub fn vconcat(parts: &[&Col<f64>]) -> Col<f64> {
    let n: usize = parts.iter().map(|p| p.nrows()).sum();
    let mut out = Col::<f64>::zeros(n);
    let mut k = 0;
    for p in parts { for i in 0..p.nrows() { out[k] = p[i]; k += 1; } }
    out
}

/// `v * s` of a dense vector
pub fn vscale(v: &Col<f64>, s: f64) -> Col<f64> { Col::<f64>::from_fn(v.nrows(), |i| v[i] * s) }

/// `class Integration`, the common part of the naive and the seamless integration. `N` is `n_big`.
pub struct Integration<'a> {
    pub g2: &'a SprsD,
    pub f2: &'a Col<f64>,
    pub x_size: usize,
    pub e_size: usize,
    pub n_big: i32,
    pub n: i32,
    pub length: f64, // grid scale
    pub w_integ: f64,
    pub w_const: f64, // a factor so that the fixed value do not move
    pub w_close: f64,
    pub w_barrier: f64,
    pub local_injectivity: bool,
    pub barrier: &'a RefCell<InjectiveBarrier>,
    pub g_integ: SprsD,
    pub g_close: SprsD,
    pub g_const: SprsD,
    pub x2f: SprsD, // convert x to the field with parameter length
}

pub(crate) fn empty_sprs() -> SprsD { from_triplets(0, 0, &[]) }

impl<'a> Integration<'a> {
    pub fn new(
        intrinsic_grad: &'a SprsD,
        intrinsic_field: &'a Col<f64>,
        n_big: i32,
        n: i32,
        loc_inj: bool,
        length: f64,
        w_integ: f64,
        w_const: f64,
        w_close: f64,
        w_barrier: f64,
        barrier: &'a RefCell<InjectiveBarrier>,
    ) -> Self {
        Integration {
            g2: intrinsic_grad,
            f2: intrinsic_field,
            x_size: 0,
            e_size: 0,
            n_big,
            n,
            length,
            w_integ,
            w_const,
            w_close,
            w_barrier,
            local_injectivity: loc_inj,
            barrier,
            g_integ: empty_sprs(),
            g_close: empty_sprs(),
            g_const: empty_sprs(),
            x2f: empty_sprs(),
        }
    }

    pub fn pre_iteration(_prevx: &Col<f64>) {}

    pub fn compute_jacobian(
        &self,
        field: &Col<f64>,
        f_integ: &Col<f64>,
        f_close: &Col<f64>,
        f_const: &Col<f64>,
        update_j: bool,
        e: &mut Col<f64>,
        j: &mut SprsD,
        lp: bool,
    ) {
        if self.local_injectivity {
            let mut g_barrier = empty_sprs();
            let mut f_barrier = Col::<f64>::zeros(0);
            self.barrier.borrow_mut().update(field, update_j, &mut f_barrier, &mut g_barrier);
            *e = vconcat(&[
                &vscale(f_integ, self.w_integ),
                &vscale(f_close, self.w_close),
                &vscale(f_const, self.w_const),
                &vscale(&f_barrier, self.w_barrier),
            ]);
            if update_j {
                // the naive branch scales the barrier by wInteg here while the residual above uses
                // wBarrier, so this J is not the exact derivative of E. it is load-bearing: the barrier
                // then reaches the normal equations with weight wInteg * wBarrier = 1 in the gradient
                // but wInteg^2 = 1e8 in the gauss-newton hessian, acting as curvature rather than as a
                // driving force. using wBarrier here instead leaves it 1e8 times weaker and nefertiti
                // at gridscale 0.02 hits the iteration cap with 10 inverted faces. the only consistent
                // variant that still converges raises wBarrier to wInteg, at 24 instead of 13 iterations
                let img = self.barrier.borrow().gen_image_filed(field);
                let gb_img = matmul(g_barrier.as_ref(), img.as_ref());
                g_barrier = if lp {
                    scale(&matmul(gb_img.as_ref(), self.x2f.as_ref()), self.w_barrier)
                } else {
                    scale(&matmul(gb_img.as_ref(), self.g_close.as_ref()), self.w_integ) // maybe wBarrier??
                };
                let idx = vec![vec![0usize], vec![1], vec![2], vec![3]];
                *j = sparse_block(&idx, &[&self.g_integ, &self.g_close, &self.g_const, &g_barrier]);
            }
            return;
        }

        *e = vconcat(&[
            &vscale(f_integ, self.w_integ),
            &vscale(f_close, self.w_close),
            &vscale(f_const, self.w_const),
        ]);
        if update_j {
            let idx = vec![vec![0usize], vec![1], vec![2]];
            *j = sparse_block(&idx, &[&self.g_integ, &self.g_close, &self.g_const]);
        }
    }
}
