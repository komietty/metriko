//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/igm/parameterization_integ.ipp (see rust/PORTED_FROM).

use faer::{Col, Mat};
use faer::sparse::Triplet;
use crate::common::typedef::*;
use crate::common::utilities::grid_unit;
use crate::hmesh::utilities::compute_poisson_weight_matrix;
use crate::igm::parameterization::*;
use crate::igm::iterative_rounding::iter_rounding::iterative_rounding;
use crate::solver::matrix_ops::*;
use crate::solver::levenberg_marquardt::spmv;

impl<'a> RosyParameterization<'a> {
    pub fn integ(&mut self) -> bool {
        let raw = self.raw;
        let cut = self.cut;
        let ext = self.ext;
        let nb = self.n_big as usize;

        let fullx2_nfunc = matmul(self.vtrans2cut.as_ref(), self.uncompress.as_ref());
        let mut cfull = matmul(self.constraint.as_ref(), self.uncompress.as_ref());
        reduce_to_linearly_independent(&mut cfull);

        //---- generating G: the matrix for generating a vector field from uv -----
        // (the C++ also assembles G3, the 3D gradient, which only the METRIKO_DEBUG evaluation below reads)
        let mut i_t: Vec<TripD> = vec![];
        for f in cut.faces() {
            for h in f.adj_halfs() {
                let g = f.normal().cross(h.next().vec()) / (f.area() * 2.);
                for k in 0..nb {
                    let c = nb * h.tail().id as usize + k;
                    let ir = 2 * nb * f.id as usize + k * 2;
                    i_t.push(Triplet::new(ir, c, g.dot(f.basis_x())));
                    i_t.push(Triplet::new(ir + 1, c, g.dot(f.basis_y())));
                }
            }
        }
        let g2 = from_triplets(2 * nb * cut.n_f as usize, nb * cut.n_v as usize, &i_t);

        // rawF2 (cut.nF x 2N), flattened row-major: the intrinsic field of face f, vector i, at 2N f + 2i + {0, 1}
        let n_f = cut.n_f as usize;
        let mut raw_f2 = Col::<f64>::zeros(n_f * 2 * nb);
        for i in 0..nb {
            for f in 0..n_f {
                let e = Row3d::new(ext[(f, 3 * i)], ext[(f, 3 * i + 1)], ext[(f, 3 * i + 2)]);
                let bx = cut.face_basis_x[f];
                let by = cut.face_basis_y[f];
                // cwiseProduct(basis).rowwise().sum()
                raw_f2[f * 2 * nb + 2 * i] = e.x * bx.x + e.y * bx.y + e.z * bx.z;
                raw_f2[f * 2 * nb + 2 * i + 1] = e.x * by.x + e.y * by.y + e.z * by.z;
            }
        }

        //------ generating locally_injectivity_check -----
        let local_injectivity = self.local_injectivity;
        let verbose = self.verbose;
        let locally_injectivity_check = |x: &Col<f64>| -> bool {
            if !local_injectivity { return false; }
            let nf = spmv(&fullx2_nfunc, x);
            // nfn_ = NF.reshaped<RowMajor>(cut.nV, N): vertex v, function k at N v + k
            for f in cut.faces() {
                let h = f.half();
                let a = Row2d::new(nf[nb * h.tail().id as usize], nf[nb * h.tail().id as usize + 1]);
                let b = Row2d::new(nf[nb * h.next().tail().id as usize], nf[nb * h.next().tail().id as usize + 1]);
                let c = Row2d::new(nf[nb * h.prev().tail().id as usize], nf[nb * h.prev().tail().id as usize + 1]);
                let d1 = b - a;
                let d2 = c - a;
                if d1.x * d2.y - d1.y * d2.x <= 0. { return false; }
            }
            if verbose { println!("now the map is locally injective"); }
            true
        };

        let length = grid_unit(&cut.pos, self.gridscale) / RosyParameterization::avg_norm(ext, 3, nb, n_f);
        let weight = compute_poisson_weight_matrix(raw, self.singlars, self.gridscale, 2 * self.n_big);
        let g2_full = matmul(g2.as_ref(), fullx2_nfunc.as_ref());
        let mut fullx = Col::<f64>::zeros(0);

        let success = iterative_rounding(
            &self.fixed_idcs,
            &self.fixed_vals,
            &weight,
            &self.singular_idcs,
            &self.integer_idcs,
            length,
            &cfull,
            &g2_full,
            self.n_big,
            self.n,
            cut.n_f,
            self.seamless,
            self.round_seams,
            self.local_injectivity,
            self.verbose,
            &raw_f2,
            &locally_injectivity_check,
            &mut fullx,
        );
        self.fullx = fullx;

        if !success && self.verbose { println!("Rounding has failed!"); }

        let nf = spmv(&fullx2_nfunc, &self.fullx);
        self.nfn = Mat::<f64>::from_fn(cut.n_v as usize, nb, |v, k| nf[nb * v + k]);

        self.cfn = Mat::<f64>::zeros(raw.n_f as usize, nb * 3);
        for j in 0..3 {
            for f in 0..raw.n_f as usize {
                let v = cut.idx[f][j] as usize;
                for k in 0..nb { self.cfn[(f, nb * j + k)] = self.nfn[(v, k)]; }
            }
        }

        // (the METRIKO_DEBUG block, which reports the deviation of the recovered field from the original one, is not
        // ported: it only prints)

        success
    }
}
