//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/solver/matrix_ops.h (see rust/PORTED_FROM).

use faer::{Col, Mat};
use faer::sparse::{SparseColMat, SparseColMatRef, Triplet};
use faer::traits::ComplexField;
use crate::common::typedef::*;

pub type SprsD = SparseColMat<usize, f64>;
pub type SprsC = SparseColMat<usize, Complex>;
pub type TripD = Triplet<usize, usize, f64>;
pub type TripC = Triplet<usize, usize, Complex>;

/// `SparseMatrix(rows, cols)` + `setFromTriplets`: duplicates are summed, row indices sorted in every column.
pub fn from_triplets<T: ComplexField>(nrows: usize, ncols: usize, t: &[Triplet<usize, usize, T>]) -> SparseColMat<usize, T> {
    SparseColMat::<usize, T>::try_new_from_triplets(nrows, ncols, t).expect("setFromTriplets")
}

/// the nonzeros of a column-major sparse matrix as `(row, col, value)`, column by column, rows in storage order
pub fn nonzeros<'a, T: Clone>(s: SparseColMatRef<'a, usize, T>) -> impl Iterator<Item = (usize, usize, T)> + 'a {
    (0..s.ncols()).flat_map(move |j| {
        let rows = s.row_idx_of_col_raw(j);
        let vals = s.val_of_col(j);
        rows.iter().zip(vals.iter()).map(move |(&i, v)| (i, j, v.clone()))
    })
}

/// an explicit column-major copy of `A^T`, rows sorted in every column
pub fn transpose<T: ComplexField>(a: SparseColMatRef<'_, usize, T>) -> SparseColMat<usize, T> {
    let t: Vec<Triplet<usize, usize, T>> = nonzeros(a).map(|(i, j, v)| Triplet::new(j, i, v)).collect();
    from_triplets(a.ncols(), a.nrows(), &t)
}

/// the sparse product `A * B` with rows sorted in every column. The pattern is the symbolic product (explicit zeros
/// from cancellation are kept, as in Eigen's product).
pub fn matmul<T: ComplexField>(a: SparseColMatRef<'_, usize, T>, b: SparseColMatRef<'_, usize, T>) -> SparseColMat<usize, T> {
    let p: SparseColMat<usize, T> = a * b;
    sorted(p.as_ref())
}

/// a copy with the row indices of every column sorted ascending
pub fn sorted<T: ComplexField>(a: SparseColMatRef<'_, usize, T>) -> SparseColMat<usize, T> {
    let t: Vec<Triplet<usize, usize, T>> = nonzeros(a).map(|(i, j, v)| Triplet::new(i, j, v)).collect();
    from_triplets(a.nrows(), a.ncols(), &t)
}

/// `A * x`
pub fn mul_vec<T: ComplexField>(a: SparseColMatRef<'_, usize, T>, x: &Col<T>) -> Col<T> { a * x }

/// Eigen's `pruned()` with the default reference 0: drops the exact zeros
pub fn pruned<T: ComplexField + PartialEq>(a: SparseColMatRef<'_, usize, T>) -> SparseColMat<usize, T> {
    let zero = T::zero_impl();
    let t: Vec<Triplet<usize, usize, T>> = nonzeros(a).filter(|(_, _, v)| *v != zero).map(|(i, j, v)| Triplet::new(i, j, v)).collect();
    from_triplets(a.nrows(), a.ncols(), &t)
}

/// the diagonal matrix of a vector (`v.asDiagonal()` converted to sparse: every diagonal entry is stored)
pub fn diagonal(v: &[f64]) -> SprsD {
    let t: Vec<TripD> = v.iter().enumerate().map(|(i, &x)| Triplet::new(i, i, x)).collect();
    from_triplets(v.len(), v.len(), &t)
}

/// `SparseMatrix::nonZeros()`
pub fn non_zeros<T>(a: SparseColMatRef<'_, usize, T>) -> usize { a.compute_nnz() }

// the solutions of C x = b for a sparse C whose nonzeros sit in a few columns, as x = x0 + Z y. the columns C does not
// touch stay free: identity columns of Z, in their order. the touched ones are confined to the kernel of C restricted to
// them, computed densely by FullPivLU, whose vectors follow in Z. x0 solves C x0 = b and is zero off the touched columns
/// `sparse_null_space(C, b, x0)`: returns `(Z, x0)`
pub fn sparse_null_space_rhs(c: &SprsD, b: &Col<f64>) -> (SprsD, Col<f64>) {
    use crate::solver::eigen::full_piv_lu::FullPivLu;
    let nx = c.ncols();
    let mut cols: Vec<usize> = vec![]; // the touched columns

    for (_, col, _) in nonzeros(c.as_ref()) { cols.push(col); }

    cols.sort_unstable();
    cols.dedup();
    let mut loc = vec![-1i64; nx];
    for (j, &col) in cols.iter().enumerate() { loc[col] = j as i64; }

    let mut t: Vec<TripD> = vec![];
    let mut ny = 0usize;
    for i in 0..nx { if loc[i] == -1 { t.push(Triplet::new(i, ny, 1.)); ny += 1; } }
    let mut x0 = Col::<f64>::zeros(nx);

    if !cols.is_empty() {
        let mut cd = Mat::<f64>::zeros(c.nrows(), cols.len());
        for (r, col, v) in nonzeros(c.as_ref()) { cd[(r, loc[col] as usize)] = v; }

        let lu = FullPivLu::new(&cd);
        let xc = lu.solve(b);
        for j in 0..cols.len() { x0[cols[j]] = xc[j]; }

        if lu.dimension_of_kernel() > 0 {
            let k = lu.kernel();
            for kk in 0..k.ncols() {
                for j in 0..k.nrows() {
                    if k[(j, kk)] != 0. { t.push(Triplet::new(cols[j], ny, k[(j, kk)])); }
                }
                ny += 1;
            }
        }
    }

    (from_triplets(nx, ny, &t), x0)
}

/// the homogeneous case, C x = 0: x = Z y
pub fn sparse_null_space(c: &SprsD) -> SprsD {
    sparse_null_space_rhs(c, &Col::<f64>::zeros(c.nrows())).0
}

pub fn reduce_to_linearly_independent(mat: &mut SprsD) {
    use crate::solver::eigen::sparse_qr::SparseQr;
    if mat.nrows() == 0 { return; }
    let qr = SparseQr::compute(transpose(mat.as_ref()).as_ref());
    let rank = qr.rank();
    let idcs = qr.cols_permutation(); // the remaining row idcs of original mtx

    let mut rid_map = vec![-1i64; mat.nrows()];
    for j in 0..rank { rid_map[idcs[j]] = j as i64; } // original mat row idx -> QR decomped mat row idx

    let mut t: Vec<TripD> = Vec::with_capacity(non_zeros(mat.as_ref()));
    for (r, c, v) in nonzeros(mat.as_ref()) {
        let rid = rid_map[r];
        if rid != -1 { t.push(Triplet::new(rid as usize, c, v)); }
    }

    *mat = from_triplets(rank, mat.ncols(), &t);
}

/// `sparse_block(idcs, mats, result)`: the block matrix whose block (i, j) is `mats[idcs(i, j)]`
pub fn sparse_block(idcs: &[Vec<usize>], mats: &[&SprsD]) -> SprsD {
    // assessing dimensions
    let row_oft = idcs.len();
    let col_oft = idcs[0].len();
    let mut row_offsets = vec![0usize; row_oft];
    let mut col_offsets = vec![0usize; col_oft];
    for i in 1..row_oft { row_offsets[i] = row_offsets[i - 1] + mats[idcs[i - 1][0]].nrows(); }
    for i in 1..col_oft { col_offsets[i] = col_offsets[i - 1] + mats[idcs[0][i - 1]].ncols(); }

    let nrows = row_offsets[row_oft - 1] + mats[idcs[row_oft - 1][0]].nrows();
    let ncols = col_offsets[col_oft - 1] + mats[idcs[0][col_oft - 1]].ncols();

    // NOTE: as in the C++, block (i, j) takes the entries of mats[i] (not mats[idcs(i, j)])
    let mut t: Vec<TripD> = vec![];
    for i in 0..row_offsets.len() {
        for j in 0..col_offsets.len() {
            for (r, c, v) in nonzeros(mats[i].as_ref()) {
                t.push(Triplet::new(row_offsets[i] + r, col_offsets[j] + c, v));
            }
        }
    }
    from_triplets(nrows, ncols, &t)
}

/// for each row, the columns holding a nonzero on it, in column order
pub fn cols_by_row(s: &SprsD) -> Vec<Vec<i32>> {
    let mut res = vec![vec![]; s.nrows()];
    for (r, c, _) in nonzeros(s.as_ref()) { res[r].push(c as i32); }
    res
}
