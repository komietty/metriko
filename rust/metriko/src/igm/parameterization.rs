//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/igm/parameterization.h (see rust/PORTED_FROM).

use faer::{Col, Mat};
use faer::sparse::Triplet;
use crate::common::typedef::*;
use crate::hmesh::hmesh::*;
use crate::solver::matrix_ops::{SprsD, from_triplets};

/// A small dense integer matrix (`MatXi` in the C++), row-major.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MatI {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<i32>,
}

impl MatI {
    pub fn zero(rows: usize, cols: usize) -> Self { MatI { rows, cols, data: vec![0; rows * cols] } }
    pub fn identity(rows: usize, cols: usize) -> Self {
        let mut m = Self::zero(rows, cols);
        for i in 0..rows.min(cols) { m.data[i * cols + i] = 1; }
        m
    }
    #[inline] pub fn get(&self, i: usize, j: usize) -> i32 { self.data[i * self.cols + j] }
    #[inline] pub fn set(&mut self, i: usize, j: usize, v: i32) { self.data[i * self.cols + j] = v; }
    /// `a * b`
    pub fn mul(&self, b: &MatI) -> MatI {
        let mut m = MatI::zero(self.rows, b.cols);
        for i in 0..self.rows {
            for j in 0..b.cols {
                let mut s = 0;
                for k in 0..self.cols { s += self.get(i, k) * b.get(k, j); }
                m.set(i, j, s);
            }
        }
        m
    }
    /// `-a`
    pub fn neg(&self) -> MatI { MatI { rows: self.rows, cols: self.cols, data: self.data.iter().map(|v| -v).collect() } }
    /// `a += b`
    pub fn add_assign(&mut self, b: &MatI) { for (x, y) in self.data.iter_mut().zip(b.data.iter()) { *x += *y; } }
    /// `a -= b`
    pub fn sub_assign(&mut self, b: &MatI) { for (x, y) in self.data.iter_mut().zip(b.data.iter()) { *x -= *y; } }
    /// `a.cwiseAbs().maxCoeff()`
    pub fn abs_max(&self) -> i32 { self.data.iter().map(|v| v.abs()).max().unwrap_or(0) }
}

pub struct RosyParameterization<'a> {
    pub(crate) raw: &'a Hmesh,
    pub(crate) cut: &'a Hmesh,
    pub(crate) ext: &'a Mat<f64>,
    pub(crate) singlars: &'a [i32],
    pub(crate) matching: &'a [i32],
    pub(crate) seam: &'a [bool],
    pub(crate) gridscale: f64, // Global scaling of the grid
    pub(crate) n_t: i32,       // Number of transitions
    pub(crate) n_r: i32,       // Number of vertices that are inside the raw mesh plus nT
    pub(crate) n_s: i32,       // Number of vertices that are singular but not on  boundaries
    pub n_big: i32,             // `N`: Uncompressed parametric functions
    pub n: i32,                 // Independent parameteric functions
    pub vtrans2cut: SprsD,      // Map between the whole mesh (vertex + translational jump) representation to the vertex-based representation on the cut mesh
    pub constraint: SprsD,      // Linear constraints (resulting from non-singular nodes)
    pub uncompress: SprsD,      // Global uncompression of n->N
    pub he2matching: Vec<i32>,  // Table of halfedge index to matching index
    pub he2transidx: Vec<i32>,  // Table of halfedge index to transiton index. Range from -nT to nT except 0. Sub 1 (or Add 1) when using as idx
    pub integer_idcs: Vec<i32>, // Indices to be integer (n * transition)
    pub singular_idcs: Vec<i32>,// Indices of singulars (used if singular has to be integer)
    pub fixed_idcs: Vec<i32>,   // Translation fixing indices
    pub fixed_vals: Vec<f64>,   // Translation fixed values at fixed indices
    pub lreductor: MatI,        // Linear Reduction tying the n dofs to the full N
    pub fullx: Col<f64>,        // Result: computed result
    pub nfn: Mat<f64>,          // Result: vertex-uv value on cut mesh (cut.nV x N)
    pub cfn: Mat<f64>,          // Result: corner-uv value on raw mesh (raw.nF x 3N)
    pub verbose: bool,          // Output the integration log.
    pub seamless: bool,         // Whether to do full translational seamless.
    pub round_seams: bool,      // Whether to round seams or round singularities
    pub local_injectivity: bool,// Enforce local injectivity; might result in failure!
}

impl<'a> RosyParameterization<'a> {
    /// `RosyParameterization(raw, cut, ext, singlars, matching, seam, degree, ratio = 0.02)`.
    /// `ext` is the extrinsic field, `raw.nF` x `3 degree` (the C++ `MatXd`).
    pub fn new(
        raw: &'a Hmesh,
        cut: &'a Hmesh,
        ext: &'a Mat<f64>,
        singlars: &'a [i32],
        matching: &'a [i32],
        seam: &'a [bool],
        degree: i32,
        ratio: f64,
    ) -> Self {
        let mut rp = RosyParameterization {
            raw,
            cut,
            ext,
            singlars,
            matching,
            seam,
            gridscale: ratio,
            n_t: 0,
            n_r: 0,
            n_s: 0,
            n_big: 0,
            n: 0,
            vtrans2cut: from_triplets(0, 0, &[]),
            constraint: from_triplets(0, 0, &[]),
            uncompress: from_triplets(0, 0, &[]),
            he2matching: vec![],
            he2transidx: vec![],
            integer_idcs: vec![],
            singular_idcs: vec![],
            fixed_idcs: vec![],
            fixed_vals: vec![],
            lreductor: MatI::zero(0, 0),
            fullx: Col::zeros(0),
            nfn: Mat::zeros(0, 0),
            cfn: Mat::zeros(0, 0),
            verbose: false,
            seamless: false,
            round_seams: true,
            local_injectivity: false,
        };
        rp.n_big = degree;
        if rp.n_big % 2 == 0 { rp.set_sign_symmetry(rp.n_big); }
        else { rp.lreductor = MatI::identity(rp.n_big as usize, rp.n_big as usize); rp.n = rp.n_big; }
        rp.fixed_vals = vec![0.; rp.n as usize];
        rp.compute_he2transidx();
        rp.compute_he2matching();
        rp.n_r = raw.n_v + rp.n_t;
        rp.n_s = raw.verts().filter(|&v| rp.is_inside_singular(v)).count() as i32;
        rp
    }

    pub(crate) fn is_inside_singular(&self, v: Vert<'_>) -> bool { !v.is_boundary() && self.singlars[v.id as usize] != 0 }

    pub(crate) fn find_first_seam_he<'m>(&self, v: Vert<'m>) -> Half<'m> {
        let mut s = v.half();
        for h in v.adj_halfs_dir(false) { if self.seam[h.edge().id as usize] { s = h; break; } }
        s
    }

    pub(crate) fn find_first_bndr_he<'m>(&self, v: Vert<'m>) -> Half<'m> {
        let mut b = v.half();
        for h in v.adj_halfs_dir(false) { if h.edge().is_boundary() { b = h; break; } }
        b
    }

    fn set_sign_symmetry(&mut self, n_big: i32) {
        if n_big % 2 != 0 { panic!("sign symmetry needs an even N"); }
        self.n = n_big / 2;
        let (rows, cols) = (n_big as usize, (n_big / 2) as usize);
        let mut l = MatI::zero(rows, cols);
        for i in 0..cols { l.set(i, i, 1); l.set(cols + i, i, -1); }
        self.lreductor = l;
    }

    #[allow(dead_code)]
    fn set_tris_symmetry(&mut self, n_big: i32) {
        // the entire first N/3 lines are symmetric w.r.t. to the next two (N/3) packets,
        if n_big % 3 != 0 { panic!("triangle symmetry needs N divisible by 3"); }
        let nb = n_big as usize;
        if n_big % 2 == 0 {
            self.n = n_big / 3;
            let c = nb / 3;
            let mut l = MatI::zero(nb, c);
            // block(0, 0, N/2, N/3) << I(N/3), -I(N/6), I(N/6): the comma initializer fills row bands, so the
            // last two N/6 blocks sit side by side under the first
            for i in 0..c { l.set(i, i, 1); }
            for i in 0..nb / 6 { l.set(c + i, i, -1); l.set(c + i, nb / 6 + i, 1); }
            // block(N/2, 0, N/2, N/3) = -block(0, 0, N/2, N/3)
            for i in 0..nb / 2 { for j in 0..c { let v = l.get(i, j); l.set(nb / 2 + i, j, -v); } }
            self.lreductor = l;
        } else {
            self.n = 2 * n_big / 3;
            let c = 2 * nb / 3;
            let mut l = MatI::zero(nb, c);
            for i in 0..c { l.set(i, i, 1); }
            for i in 0..nb / 3 { l.set(c + i, i, -1); }
            for i in 0..nb / 3 { l.set(c + i, nb / 3 + i, -1); }
            self.lreductor = l;
        }
    }

    /// `field.reshaped<RowMajor>(nF * N, dim).rowwise().norm().mean()`
    pub(crate) fn avg_norm(field: &Mat<f64>, dim: usize, n_big: usize, n_f: usize) -> f64 {
        let mut sum = 0.;
        for f in 0..n_f {
            for k in 0..n_big {
                let mut sq = 0.;
                for d in 0..dim { let x = field[(f, k * dim + d)]; sq += x * x; }
                sum += sq.sqrt();
            }
        }
        sum / (n_f * n_big) as f64
    }
}

//----- utility function to cut mesh by seam -----//
pub fn compute_cut_mesh(m: &Hmesh, seam: &[bool]) -> Hmesh {
    let mut cut_f = vec![[0i32; 3]; m.n_f as usize];
    let mut cut2pos: Vec<Row3d> = vec![];

    for v in m.verts() {
        let mut bgn = v.half();
        if v.is_boundary() { for h in v.adj_halfs_dir(false) { if h.edge().is_boundary() { bgn = h; break; } } }
        else               { for h in v.adj_halfs_dir(false) { if seam[h.edge().id as usize] { bgn = h; break; } } }

        for h in v.adj_halfs_from(bgn, true) {
            if h.is_boundary() { continue; }
            if seam[h.edge().id as usize] || h.id == bgn.id { cut2pos.push(v.pos()); }
            for j in 0..3 {
                let i_f = h.face().id as usize;
                if m.idx[i_f][j] == v.id { cut_f[i_f][j] = cut2pos.len() as i32 - 1; }
            }
        }
    }

    Hmesh::new(&cut2pos, &cut_f, false)
}

/// `assign_block(triplets, mat, row_bgn, col_bgn)`: the nonzeros of `mat`, row by row, shifted.
pub fn assign_block(triplets: &mut Vec<Triplet<usize, usize, f64>>, mat: &MatI, row_bgn: i32, col_bgn: i32) {
    for i in 0..mat.rows {
        for j in 0..mat.cols {
            let val = mat.get(i, j);
            if val != 0 { // need to check it is possible to cast T1 to T2
                triplets.push(Triplet::new((row_bgn + i as i32) as usize, (col_bgn + j as i32) as usize, val as f64));
            }
        }
    }
}
