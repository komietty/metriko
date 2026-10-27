//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/solver/sparse_solver.h (see rust/PORTED_FROM).

// The C++ uses Eigen's SimplicialLDLT / SimplicialLLT / SparseLU (the reference configuration, without SuiteSparse).
// Here they are faer 0.24.4 behind the same interface: analyze_pattern / factorize / compute / info / solve.
// DEVIATION (accepted, see PORTING.md): faer's fill reducing orderings and its simplicial/supernodal kernels differ
// from Eigen's, so solutions agree with the C++ up to rounding, not bitwise. The failure semantics are Eigen's:
// LDLT fails on an exactly zero pivot, LLT on a non positive one, LU on a structurally or numerically zero pivot.

use faer::{Col, Conj, Mat, MatMut, Par, Side};
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::perm::PermRef;
use faer::sparse::SparseColMatRef;
use faer::sparse::linalg::cholesky::{
    factorize_symbolic_cholesky, CholeskySymbolicParams, LdltRef, LltRef, SymbolicCholesky, SymmetricOrdering,
};
use faer::sparse::linalg::colamd;
use faer::sparse::linalg::lu::simplicial::{
    factorize_simplicial_numeric_lu, factorize_simplicial_numeric_lu_scratch, SimplicialLu,
};
use faer::linalg::cholesky::ldlt::factor::LdltRegularization;
use faer::linalg::cholesky::llt::factor::LltRegularization;
use faer::traits::ComplexField;
use crate::common::typedef::*;
use crate::solver::matrix_ops::SprsC;

/// Eigen's `ComputationInfo`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComputationInfo {
    Success,
    NumericalIssue,
    InvalidInput,
}

/// Eigen's `SimplicialLLT` (`LDLT = false`) / `SimplicialLDLT` (`LDLT = true`) with `UpLo = Lower` and AMD ordering:
/// only the lower triangle of the matrix is read; for complex scalars the factorization is `L D L^H`.
pub struct SparseCholesky<T: ComplexField, const LDLT: bool> {
    symbolic: Option<SymbolicCholesky<usize>>,
    values: Vec<T>,
    info: ComputationInfo,
    factorized: bool,
}

pub type SparseLdlt<T> = SparseCholesky<T, true>;
pub type SparseLlt<T> = SparseCholesky<T, false>;

impl<T: ComplexField, const LDLT: bool> Default for SparseCholesky<T, LDLT> {
    fn default() -> Self { Self::new() }
}

impl<T: ComplexField, const LDLT: bool> SparseCholesky<T, LDLT> {
    pub fn new() -> Self {
        SparseCholesky { symbolic: None, values: vec![], info: ComputationInfo::Success, factorized: false }
    }

    /// `SparseLDLT<M> ldlt(A)` / `ldlt.compute(A)`
    pub fn from_matrix(a: SparseColMatRef<'_, usize, T>) -> Self {
        let mut s = Self::new();
        s.compute(a);
        s
    }

    pub fn compute(&mut self, a: SparseColMatRef<'_, usize, T>) {
        self.analyze_pattern(a);
        self.factorize(a);
    }

    pub fn analyze_pattern(&mut self, a: SparseColMatRef<'_, usize, T>) {
        self.factorized = false;
        match factorize_symbolic_cholesky(a.symbolic(), Side::Lower, SymmetricOrdering::Amd, CholeskySymbolicParams::default()) {
            Ok(s) => { self.symbolic = Some(s); self.info = ComputationInfo::Success; }
            Err(_) => { self.symbolic = None; self.info = ComputationInfo::InvalidInput; }
        }
    }

    pub fn factorize(&mut self, a: SparseColMatRef<'_, usize, T>) {
        self.factorized = false;
        let Some(sym) = self.symbolic.as_ref() else { self.info = ComputationInfo::InvalidInput; return; };
        self.values.clear();
        self.values.resize(sym.len_val(), T::zero_impl());
        let ok = if LDLT {
            let mut mem = MemBuffer::new(sym.factorize_numeric_ldlt_scratch::<T>(Par::Seq, Default::default()));
            sym.factorize_numeric_ldlt(&mut self.values, a, Side::Lower, LdltRegularization::default(), Par::Seq, MemStack::new(&mut mem), Default::default()).is_ok()
        } else {
            let mut mem = MemBuffer::new(sym.factorize_numeric_llt_scratch::<T>(Par::Seq, Default::default()));
            sym.factorize_numeric_llt(&mut self.values, a, Side::Lower, LltRegularization::default(), Par::Seq, MemStack::new(&mut mem), Default::default()).is_ok()
        };
        self.info = if ok { ComputationInfo::Success } else { ComputationInfo::NumericalIssue };
        self.factorized = ok;
    }

    pub fn info(&self) -> ComputationInfo { self.info }

    /// `solve(B)` for a matrix of right hand sides. After a failed factorization Eigen returns values from the
    /// partial factor; callers only use them after checking `info()`, here they are NaN.
    pub fn solve_mat(&self, b: &Mat<T>) -> Mat<T> {
        let mut x = b.clone();
        self.solve_in_place(x.as_mut());
        x
    }

    pub fn solve(&self, b: &Col<T>) -> Col<T> {
        let mut x = Mat::<T>::from_fn(b.nrows(), 1, |i, _| b[i].clone());
        self.solve_in_place(x.as_mut());
        Col::<T>::from_fn(b.nrows(), |i| x[(i, 0)].clone())
    }

    fn solve_in_place(&self, mut x: MatMut<'_, T>) {
        let sym = match (&self.symbolic, self.factorized) {
            (Some(s), true) => s,
            _ => {
                let nan = T::from_f64_impl(f64::NAN);
                for j in 0..x.ncols() { for i in 0..x.nrows() { x[(i, j)] = nan.clone(); } }
                return;
            }
        };
        let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T>(x.ncols(), Par::Seq));
        if LDLT {
            LdltRef::new(sym, &self.values).solve_in_place_with_conj(Conj::No, x, Par::Seq, MemStack::new(&mut mem));
        } else {
            LltRef::new(sym, &self.values).solve_in_place_with_conj(Conj::No, x, Par::Seq, MemStack::new(&mut mem));
        }
    }
}

/// Eigen's `SparseLU<SparseMatrix<double>>`: a fill reducing column ordering (COLAMD) and LU with partial pivoting.
/// `info()` is `NumericalIssue` when a column has no nonzero pivot candidate, structurally or numerically, as in
/// Eigen's `pivotL` (`pivmax <= 0`).
pub struct SparseLu {
    n: usize,
    col_perm_fwd: Vec<usize>,
    col_perm_inv: Vec<usize>,
    row_perm_fwd: Vec<usize>,
    row_perm_inv: Vec<usize>,
    lu: SimplicialLu<usize, f64>,
    info: ComputationInfo,
}

impl SparseLu {
    /// `Eigen::SparseLU solver(A)` / `lu.compute(A)`
    pub fn from_matrix(a: SparseColMatRef<'_, usize, f64>) -> Self {
        let n = a.nrows();
        let mut s = SparseLu {
            n,
            col_perm_fwd: vec![0; n],
            col_perm_inv: vec![0; n],
            row_perm_fwd: vec![0; n],
            row_perm_inv: vec![0; n],
            lu: SimplicialLu::new(),
            info: ComputationInfo::Success,
        };
        if a.nrows() != a.ncols() { s.info = ComputationInfo::InvalidInput; return s; }
        {
            let nnz = a.compute_nnz();
            let mut mem = MemBuffer::new(colamd::order_scratch::<usize>(n, n, nnz));
            if colamd::order(&mut s.col_perm_fwd, &mut s.col_perm_inv, a.symbolic(), Default::default(), MemStack::new(&mut mem)).is_err() {
                s.info = ComputationInfo::InvalidInput;
                return s;
            }
        }
        let col_perm = PermRef::new_checked(&s.col_perm_fwd, &s.col_perm_inv, n);
        let mut mem = MemBuffer::new(factorize_simplicial_numeric_lu_scratch::<usize, f64>(n, n));
        let res = factorize_simplicial_numeric_lu(&mut s.row_perm_fwd, &mut s.row_perm_inv, &mut s.lu, a, col_perm, MemStack::new(&mut mem));
        if res.is_err() { s.info = ComputationInfo::NumericalIssue; return s; }
        // faer stores the row perm inverse in `row_perm_fwd` during the factorization and fixes it at the end; the
        // pivot of column j is the last entry of U's column j
        let u = s.lu.u_factor_unsorted();
        for j in 0..n {
            let vals = u.val_of_col(j);
            let rows = u.row_idx_of_col_raw(j);
            match rows.iter().position(|&r| r == j) {
                Some(p) if vals[p] != 0. => {}
                _ => { s.info = ComputationInfo::NumericalIssue; break; }
            }
        }
        s
    }

    pub fn info(&self) -> ComputationInfo { self.info }

    pub fn solve_mat(&self, b: &Mat<f64>) -> Mat<f64> {
        let mut x = b.clone();
        let row_perm = PermRef::new_checked(&self.row_perm_fwd, &self.row_perm_inv, self.n);
        let col_perm = PermRef::new_checked(&self.col_perm_fwd, &self.col_perm_inv, self.n);
        let mut mem = MemBuffer::new(faer::perm::permute_rows_in_place_scratch::<usize, f64>(self.n, x.ncols()));
        self.lu.solve_in_place_with_conj(row_perm, col_perm, Conj::No, x.as_mut(), Par::Seq, MemStack::new(&mut mem));
        x
    }

    pub fn solve(&self, b: &Col<f64>) -> Col<f64> {
        let m = Mat::<f64>::from_fn(b.nrows(), 1, |i, _| b[i]);
        let x = self.solve_mat(&m);
        Col::<f64>::from_fn(b.nrows(), |i| x[(i, 0)])
    }
}

// A is hermitian; only its lower triangle is read
pub fn solve_hermitian(a: &SprsC, b: &Col<Complex>) -> Col<Complex> {
    let ldlt = SparseLdlt::<Complex>::from_matrix(a.as_ref());
    if ldlt.info() != ComputationInfo::Success {
        metriko_fail!("failed to factorize a hermitian system (info {}): the mesh is likely degenerate, too coarsely tessellated, or contains sharp features / bad triangles", ldlt.info() as i32);
    }
    ldlt.solve(b)
}

/// `std::mt19937`
pub struct Mt19937 {
    mt: [u32; 624],
    idx: usize,
}

impl Mt19937 {
    pub fn new(seed: u32) -> Self {
        let mut mt = [0u32; 624];
        mt[0] = seed;
        for i in 1..624 {
            mt[i] = 1812433253u32.wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_add(i as u32);
        }
        Mt19937 { mt, idx: 624 }
    }

    pub fn next_u32(&mut self) -> u32 {
        if self.idx >= 624 {
            for i in 0..624 {
                let y = (self.mt[i] & 0x8000_0000) | (self.mt[(i + 1) % 624] & 0x7fff_ffff);
                let mut v = self.mt[(i + 397) % 624] ^ (y >> 1);
                if y & 1 != 0 { v ^= 0x9908_b0df; }
                self.mt[i] = v;
            }
            self.idx = 0;
        }
        let mut y = self.mt[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }
}

/// libc++ `std::uniform_real_distribution<double>(a, b)(g)` for `g = std::mt19937`:
/// `(b - a) * generate_canonical<double, 53>(g) + a`, where generate_canonical draws k = 2 words and returns
/// `(w0 + w1 * 2^32) / 2^64`.
pub fn uniform_real(g: &mut Mt19937, a: f64, b: f64) -> f64 {
    let rp = 4294967296.0f64; // _URNG::max() - _URNG::min() + 1
    let mut base = rp;
    let mut sp = g.next_u32() as f64;
    sp += g.next_u32() as f64 * base;
    base *= rp;
    (b - a) * (sp / base) + a
}

// inverse power iteration for the smallest generalized eigenvector of L x = s M x. iterates until the rayleigh
// quotient s = x* L x / x* M x stops changing (relative change below tol), or maxIter is reached
/// `solve_smallest_eig(L, M, maxIter = 500, tol = 1e-8)`
pub fn solve_smallest_eig(l: &SprsC, m: &SprsC, max_iter: i32, tol: f64) -> Col<Complex> {
    let ldlt = SparseLdlt::<Complex>::from_matrix(l.as_ref());
    if ldlt.info() != ComputationInfo::Success {
        metriko_fail!("failed to factorize the connection Laplacian (info {}): the mesh is likely degenerate, too coarsely tessellated, or contains sharp features / bad triangles", ldlt.info() as i32);
    }

    // fixed-seed start vector: the result must not depend on what else consumed std::rand in the process.
    // real and imaginary parts are drawn in separate statements to keep the draw order defined
    let mut rng = Mt19937::new(0);
    let n = l.nrows();
    let mut x = Col::<Complex>::zeros(n);
    for i in 0..n {
        let re = uniform_real(&mut rng, -1., 1.);
        let im = uniform_real(&mut rng, -1., 1.);
        x[i] = cx(re, im);
    }

    // Eigen's complex dot: sum conj(a_i) * b_i
    let dotc = |a: &Col<Complex>, b: &Col<Complex>| -> Complex {
        let mut s = cx(0., 0.);
        for i in 0..a.nrows() { s += a[i].conj() * b[i]; }
        s
    };

    let mut s_prev = f64::INFINITY;
    for _ in 0..max_iter {
        let mx: Col<Complex> = m * &x;
        x = ldlt.solve(&mx);
        let mx: Col<Complex> = m * &x;
        let nrm = dotc(&x, &mx).norm().sqrt();
        for i in 0..n { x[i] /= nrm; } // x* M x = 1
        let lx: Col<Complex> = l * &x;
        let s = dotc(&x, &lx).norm(); // rayleigh quotient
        if (s - s_prev).abs() <= tol * s.abs() { break; }
        s_prev = s;
    }
    x
}
