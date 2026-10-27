//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of Eigen 5.0.1 Eigen/src/SparseQR/SparseQR.h (MPL-2.0): SparseQR<SparseMatrix<double>,
//--- COLAMDOrdering<int>>, analyzePattern + factorize with the default (rank revealing) threshold.
//
// metriko only reads rank() and colsPermutation(). The factors Q (Householder vectors), R and hcoeffs are still
// built as Eigen builds them, since every column's elimination reads the previous Householder vectors in their
// storage order (Q columns are filled unordered), and R is kept for inspection.

use super::colamd::colamd_ordering;
use super::coletree::{coletree, ColPattern};
use faer::sparse::SparseColMatRef;

/// `ComputationInfo` of the factorization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QrInfo {
    Success,
    InvalidInput,
}

/// A column of a sparse factor filled the way SparseMatrix's `insertBack*` fills it: rows in insertion order.
#[derive(Clone, Default, Debug)]
pub struct SparseCol {
    pub rows: Vec<usize>,
    pub vals: Vec<f64>,
}

/// Sparse left-looking QR factorization with numerical column pivoting.
pub struct SparseQr {
    nrows: usize,
    ncols: usize,
    info: QrInfo,
    last_error: String,
    r: Vec<SparseCol>,           // the triangular factor matrix, by column (before the pivot permutation)
    q: Vec<SparseCol>,           // the Householder vectors, by column
    hcoeffs: Vec<f64>,           // the Householder coefficients
    perm_c: Vec<usize>,          // fill-reducing column permutation (perm_c[old] = new)
    pivotperm: Vec<usize>,       // the permutation for rank revealing
    output_perm_c: Vec<usize>,   // the final column permutation
    nonzeropivots: usize,        // number of non zero pivots
    etree: Vec<i32>,             // column elimination tree
    first_row_elt: Vec<i32>,     // first element in each row
}

impl SparseQr {
    /// `compute(mat)`: `analyzePattern(mat)` then `factorize(mat)`. `mat` must be column-major and compressed.
    pub fn compute(mat: SparseColMatRef<'_, usize, f64>) -> Self {
        let mut qr = Self {
            nrows: mat.nrows(),
            ncols: mat.ncols(),
            info: QrInfo::Success,
            last_error: String::new(),
            r: vec![],
            q: vec![],
            hcoeffs: vec![],
            perm_c: vec![],
            pivotperm: vec![],
            output_perm_c: vec![],
            nonzeropivots: 0,
            etree: vec![],
            first_row_elt: vec![],
        };
        qr.analyze_pattern(mat);
        qr.factorize(mat);
        qr
    }

    /// Preprocessing step of a QR factorization: the fill-reducing permutation is computed and applied to the
    /// columns of A and the column elimination tree is computed as well. Only the sparsity pattern of `mat` is
    /// exploited. It is assumed that there is no empty row in the matrix `mat`.
    fn analyze_pattern(&mut self, mat: SparseColMatRef<'_, usize, f64>) {
        let n = mat.ncols();
        let m = mat.nrows();
        let col_ptr = mat.col_ptr();
        let row_idx = mat.row_idx();
        assert!(mat.col_nnz().is_none(), "SparseQR requires a sparse matrix in compressed mode");

        // Compute the column fill reducing ordering
        self.perm_c = colamd_ordering(m, n, col_ptr, row_idx);
        if self.perm_c.is_empty() { self.perm_c = (0..n).collect(); }

        // Compute the column elimination tree of the permuted matrix
        self.output_perm_c = inverse(&self.perm_c);
        let col_len: Vec<usize> = (0..n).map(|j| col_ptr[j + 1] - col_ptr[j]).collect();
        let pattern = ColPattern { nrows: m, ncols: n, col_beg: &col_ptr[..n], col_len: &col_len, row_idx };
        let perm: Vec<i32> = self.output_perm_c.iter().map(|&x| x as i32).collect();
        coletree(&pattern, &mut self.etree, &mut self.first_row_elt, Some(&perm));
    }

    /// Performs the numerical QR factorization of the input matrix (same sparsity pattern as in analyze_pattern).
    fn factorize(&mut self, mat: SparseColMatRef<'_, usize, f64>) {
        let m = mat.nrows();
        let n = mat.ncols();
        let diag_size = m.min(n);
        // Record the visited nodes. Eigen sizes it max(m, n); the etree root sentinel n can be marked when m <= n,
        // which Eigen writes one past the end: the extra entry holds that write.
        let mut mark = vec![-1i64; m.max(n) + 1];
        // Store temporarily the row indexes for the current column of R and Q. Eigen sizes Qidx m, but writes one
        // past it when the column count exceeds the row count (see DEVIATION below); m + 1 entries hold every write.
        let mut ridx = vec![0usize; n];
        let mut qidx = vec![0usize; m + 1];
        let mut tval = vec![0f64; m]; // The dense vector used to compute the current column

        self.r = vec![SparseCol::default(); n];
        self.q = vec![SparseCol::default(); diag_size];
        self.hcoeffs = vec![0.; diag_size];

        // m_pmat = mat, with the fill-in reducing permutation applied lazily: the outer index of column p is the one
        // of the original column i with perm_c[i] = p
        let col_ptr = mat.col_ptr();
        let row_idx = mat.row_idx();
        let values = mat.val();
        let mut pcol_beg = vec![0usize; n];
        let mut pcol_len = vec![0usize; n];
        for i in 0..n {
            let p = self.perm_c[i];
            pcol_beg[p] = col_ptr[i];
            pcol_len[p] = col_ptr[i + 1] - col_ptr[i];
        }
        let pmat = ColPattern { nrows: m, ncols: n, col_beg: &pcol_beg, col_len: &pcol_len, row_idx };

        // Compute the default threshold as in MatLab, see: Tim Davis, "Algorithm 915, SuiteSparseQR: Multifrontal
        // Multithreaded Rank-Revealing Sparse QR Factorization, ACM Trans. on Math. Soft. 38(1), 2011, Page 8:3
        let pivot_threshold = {
            let mut max2_norm = 0f64;
            for j in 0..n {
                // m_pmat.col(j).norm(): sqrt of the sequential sum of the squares
                let mut s = 0.;
                for p in pcol_beg[j]..pcol_beg[j] + pcol_len[j] { s += values[p] * values[p]; }
                max2_norm = max2_norm.max(s.sqrt());
            }
            if max2_norm == 0. { max2_norm = 1.; }
            (20 * (m as i32 + n as i32)) as f64 * max2_norm * f64::EPSILON
        };

        // Initialize the numerical permutation
        self.pivotperm = (0..n).collect();

        let mut nonzero_col = 0usize; // Record the number of valid pivots

        // Left looking rank-revealing QR factorization: compute a column of R and Q at a time
        for col in 0..n {
            for x in mark.iter_mut() { *x = -1; }
            mark[nonzero_col] = col as i64;
            qidx[0] = nonzero_col;
            let mut nzcol_r = 0usize;
            let mut nzcol_q = 1usize;
            let mut found_diag = nonzero_col >= m;
            for x in tval.iter_mut() { *x = 0.; }

            // Symbolic factorization: find the nonzero locations of the column k of the factors R and Q, i.e., all
            // the nodes (with indexes lower than rank) reachable through the column elimination tree (etree) rooted
            // at node k. Note: if the diagonal entry does not exist, then its contribution must be explicitly added,
            // thus the trick with found_diag that permits to do one more iteration on the diagonal element if this
            // one has not been found.
            let beg = pcol_beg[col];
            let len = pcol_len[col];
            let mut itp = 0usize;
            while itp < len || !found_diag {
                let mut cur_idx = nonzero_col;
                if itp < len { cur_idx = row_idx[beg + itp]; }
                if cur_idx == nonzero_col { found_diag = true; }

                // Get the nonzeros indexes of the current column of R
                let mut st = self.first_row_elt[cur_idx]; // The traversal of the etree starts here
                if st < 0 {
                    self.last_error = "Empty row found during numerical factorization".to_string();
                    self.info = QrInfo::InvalidInput;
                    return;
                }

                // Traverse the etree
                let bi = nzcol_r;
                while mark[st as usize] != col as i64 {
                    ridx[nzcol_r] = st as usize; // Add this row to the list,
                    mark[st as usize] = col as i64; // and mark this row as visited
                    nzcol_r += 1;
                    st = self.etree[st as usize];
                }

                // Reverse the list to get the topological ordering
                let nt = nzcol_r - bi;
                for i in 0..nt / 2 { ridx.swap(bi + i, nzcol_r - i - 1); }

                // Copy the current (curIdx,pcol) value of the input matrix
                if itp < len { tval[cur_idx] = values[beg + itp]; } else { tval[cur_idx] = 0.; }

                // Compute the pattern of Q(:,k)
                if cur_idx > nonzero_col && mark[cur_idx] != col as i64 {
                    qidx[nzcol_q] = cur_idx; // Add this row to the pattern of Q,
                    mark[cur_idx] = col as i64; // and mark it as visited
                    nzcol_q += 1;
                }
                itp += 1;
            }

            // Browse all the indexes of R(:,col) in reverse order
            for i in (0..nzcol_r).rev() {
                let cur_idx = ridx[i];

                // Apply the curIdx-th householder vector to the current column (temporarily stored into tval)
                // First compute q' * tval (SparseMatrixBase::dot: two interleaved accumulators, x * y + z unfused)
                // a Q column not started yet is empty (Eigen's iterator over it is empty too)
                let empty = SparseCol::default();
                let qc = self.q.get(cur_idx).unwrap_or(&empty);
                let mut res1 = 0.;
                let mut res2 = 0.;
                let mut it = 0;
                while it < qc.rows.len() {
                    res1 = qc.vals[it] * tval[qc.rows[it]] + res1;
                    it += 1;
                    if it < qc.rows.len() { res2 = qc.vals[it] * tval[qc.rows[it]] + res2; }
                    it += 1;
                }
                let mut tdot = res1 + res2;
                tdot *= self.hcoeffs.get(cur_idx).copied().unwrap_or(0.);

                // Then update tval = tval - q * tau
                for (k, &row) in qc.rows.iter().enumerate() { tval[row] -= tdot * qc.vals[k]; }

                // Detect fill-in for the current column of Q
                if self.etree[ridx[i]] == nonzero_col as i32 {
                    for &iq in &qc.rows {
                        if mark[iq] != col as i64 {
                            qidx[nzcol_q] = iq; // Add this row to the pattern of Q,
                            nzcol_q += 1;
                            mark[iq] = col as i64; // and mark it as visited
                        }
                    }
                }
            } // End update current column

            let mut tau = 0.;
            let mut beta = 0.;

            if nonzero_col < diag_size {
                // Compute the Householder reflection that eliminate the current column
                let c0 = if nzcol_q > 0 { tval[qidx[0]] } else { 0. };

                // First, the squared norm of Q((col+1):m, col)
                let mut sqr_norm = 0.;
                for itq in 1..nzcol_q { sqr_norm += tval[qidx[itq]] * tval[qidx[itq]]; }
                if sqr_norm == 0. {
                    beta = c0;
                    tval[qidx[0]] = 1.;
                } else {
                    beta = (c0 * c0 + sqr_norm).sqrt();
                    if c0 >= 0. { beta = -beta; }
                    tval[qidx[0]] = 1.;
                    for itq in 1..nzcol_q { tval[qidx[itq]] /= c0 - beta; }
                    tau = (beta - c0) / beta;
                }
            }

            // Insert values in R
            for i in (0..nzcol_r).rev() {
                let cur_idx = ridx[i];
                if cur_idx < nonzero_col {
                    self.r[col].rows.push(cur_idx);
                    self.r[col].vals.push(tval[cur_idx]);
                    tval[cur_idx] = 0.;
                }
            }

            if nonzero_col < diag_size && beta.abs() >= pivot_threshold {
                self.r[col].rows.push(nonzero_col);
                self.r[col].vals.push(beta);
                // The householder coefficient
                self.hcoeffs[nonzero_col] = tau;
                // Record the householder reflections
                for itq in 0..nzcol_q {
                    let iq = qidx[itq];
                    self.q[nonzero_col].rows.push(iq);
                    self.q[nonzero_col].vals.push(tval[iq]);
                    tval[iq] = 0.;
                }
                nonzero_col += 1;
            } else {
                // Zero pivot found: move implicitly this column to the end
                for j in nonzero_col..n - 1 { self.pivotperm.swap(j, j + 1); }

                // Recompute the column elimination tree
                let perm: Vec<i32> = self.pivotperm.iter().map(|&x| x as i32).collect();
                coletree(&pmat, &mut self.etree, &mut self.first_row_elt, Some(&perm));
            }
        }

        for h in self.hcoeffs.iter_mut().skip(nonzero_col) { *h = 0.; }

        self.nonzeropivots = nonzero_col;

        if nonzero_col < n {
            // Permute the triangular factor to put the 'dead' columns to the end: column j of R becomes column
            // pivotperm[j]
            // (permute_outer of SparsePermutation.h copies each column as it is stored)
            let old = std::mem::take(&mut self.r);
            self.r = (0..n).map(|j| old[self.pivotperm[j]].clone()).collect();

            // Update the column permutation: m_outputPerm_c * m_pivotperm
            self.output_perm_c = (0..n).map(|i| self.output_perm_c[self.pivotperm[i]]).collect();
        }

        self.info = QrInfo::Success;
    }

    pub fn rows(&self) -> usize { self.nrows }
    pub fn cols(&self) -> usize { self.ncols }

    /// The number of non linearly dependent columns as determined by the pivoting threshold.
    pub fn rank(&self) -> usize { self.nonzeropivots }

    /// `colsPermutation().indices()`: the combined column permutation (fill-reducing, then rank revealing).
    pub fn cols_permutation(&self) -> &[usize] { &self.output_perm_c }

    /// `Success` if the factorization went fine, `InvalidInput` on an empty row.
    pub fn info(&self) -> QrInfo { self.info }

    /// The reason of an `InvalidInput`.
    pub fn last_error_message(&self) -> &str { &self.last_error }

    /// The columns of R, after the permutation that moves the dead columns to the end; rows in storage order.
    pub fn matrix_r(&self) -> &[SparseCol] { &self.r }

    /// The Householder vectors (columns of Q's compact representation), rows in storage order.
    pub fn householder_vectors(&self) -> &[SparseCol] { &self.q }

    /// The Householder coefficients.
    pub fn householder_coeffs(&self) -> &[f64] { &self.hcoeffs }
}

/// `PermutationMatrix::inverse().indices()`.
fn inverse(perm: &[usize]) -> Vec<usize> {
    let mut inv = vec![0; perm.len()];
    for (i, &p) in perm.iter().enumerate() { inv[p] = i; }
    inv
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::sparse::{SparseColMat, SymbolicSparseColMat};

    // Constraint matrices C (rows x cols triplets) of the comparison against Eigen 5.0.1: rank() and
    // colsPermutation() of SparseQR<SparseMatrix<double>, COLAMDOrdering<int>>::compute(C.transpose()), as
    // matrix_ops.h's reduce_to_linearly_independent calls it.
    const QR0: (usize, usize, &[(usize, usize, f64)], usize, &str) = (9, 7, &[(0, 4, -1.0), (0, 2, -0.21970433928072453),
        (1, 0, 1.0), (1, 5, -0.3930897144600749), (2, 0, -1.0), (2, 4, 1.0), (2, 4, -1.0), (2, 5, 0.9559773681685328),
        (3, 0, -1.0), (3, 4, 1.0), (3, 4, -1.0), (3, 5, 0.9559773681685328), (3, 4, -1.0), (3, 2, -0.21970433928072453),
        (4, 0, -1.0), (4, 4, 1.0), (4, 4, -1.0), (4, 5, 0.9559773681685328), (4, 4, -1.0), (4, 2, -0.21970433928072453),
        (4, 0, -1.0), (4, 4, 1.0), (4, 4, -1.0), (4, 5, 0.9559773681685328), (4, 4, -1.0), (4, 2, -0.21970433928072453),
        (5, 5, 1.0), (5, 4, -1.0), (5, 3, -1.0), (5, 4, -1.0), (5, 4, -0.47402756102383137), (6, 0, -1.0), (6, 0, -1.0),
        (6, 4, -0.9989459738135338), (7, 0, -1.0), (7, 4, 1.0), (7, 4, -1.0), (7, 5, 0.9559773681685328), (7, 4, -1.0),
        (7, 2, -0.21970433928072453), (7, 0, -1.0), (7, 0, -1.0), (7, 4, -0.9989459738135338), (8, 2, -1.0), (8, 2, 1.0),
        (8, 4, -1.0), (8, 2, -1.0), (8, 4, -0.5946209393441677)], 5, "0 8 5 1 2 6 3 4 7");
    const QR1: (usize, usize, &[(usize, usize, f64)], usize, &str) = (6, 8, &[(0, 2, -1.0), (0, 7, -1.0), (0, 1, -1.0),
        (0, 0, -1.0), (1, 4, -1.0), (1, 1, -1.0), (1, 4, 1.0), (2, 2, -1.0), (2, 7, -1.0), (2, 1, -1.0), (2, 0, -1.0),
        (2, 4, -1.0), (2, 1, -1.0), (2, 4, 1.0), (3, 0, -1.0), (3, 7, 1.0), (3, 2, -1.0), (3, 6, -1.0), (4, 2, -1.0), (4,
        7, -1.0), (4, 1, -1.0), (4, 0, -1.0), (4, 4, -1.0), (4, 1, -1.0), (4, 4, 1.0), (5, 1, 1.0), (5, 1, 1.0)], 3,
        "1 0 3 5 2 4");

    /// `SparseMatrix<double>(C.transpose())` for C = setFromTriplets(trips): duplicates summed in insertion order,
    /// explicit zeros kept, rows of every column sorted.
    fn transposed(rows: usize, cols: usize, trips: &[(usize, usize, f64)]) -> SparseColMat<usize, f64> {
        let mut c: Vec<std::collections::BTreeMap<usize, f64>> = vec![Default::default(); cols];
        for &(i, j, v) in trips { *c[j].entry(i).or_insert(0.) += v; }
        let mut t: Vec<Vec<(usize, f64)>> = vec![vec![]; rows];
        for (j, col) in c.iter().enumerate() { for (&i, &v) in col { t[i].push((j, v)); } }
        let mut cp = vec![0usize];
        let mut ri = vec![];
        let mut vals = vec![];
        for col in &t { for &(r, v) in col { ri.push(r); vals.push(v); } cp.push(ri.len()); }
        SparseColMat::new(SymbolicSparseColMat::new_checked(cols, rows, cp, None, ri), vals)
    }

    #[test]
    fn matches_eigen_0() {
        let (rows, cols, trips, rank, perm) = QR0;
        let ct = transposed(rows, cols, trips);
        let qr = SparseQr::compute(ct.as_ref());
        assert_eq!(qr.info(), QrInfo::Success);
        assert_eq!(qr.rank(), rank);
        assert_eq!(qr.cols_permutation().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "), perm);
    }

    #[test]
    fn matches_eigen_1() {
        let (rows, cols, trips, rank, perm) = QR1;
        let ct = transposed(rows, cols, trips);
        let qr = SparseQr::compute(ct.as_ref());
        assert_eq!(qr.info(), QrInfo::Success);
        assert_eq!(qr.rank(), rank);
        assert_eq!(qr.cols_permutation().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "), perm);
    }
}
