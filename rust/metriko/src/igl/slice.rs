//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/slice.cpp (MPL-2.0).

use faer::sparse::Triplet;
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets, nonzeros};

/// `igl::slice(X, R, C, Y)` for a sparse X: `Y = X(R, C)`, rows and columns may repeat.
pub fn slice_rc(x: &SprsD, r: &[i32], c: &[i32]) -> SprsD {
    let xm = x.nrows();
    let xn = x.ncols();
    let ym = r.len();
    let yn = c.len();
    // special case when R or C is empty
    if ym == 0 || yn == 0 { return from_triplets(ym, yn, &[]); }
    // Build reindexing maps for columns and rows
    let mut ri: Vec<Vec<usize>> = vec![vec![]; xm];
    for (i, &rr) in r.iter().enumerate() { ri[rr as usize].push(i); }
    let mut ci: Vec<Vec<usize>> = vec![vec![]; xn];
    for (i, &cc) in c.iter().enumerate() { ci[cc as usize].push(i); }
    let mut entries: Vec<TripD> = vec![];
    // Iterate over outside, then inside
    for (row, col, v) in nonzeros(x.as_ref()) {
        for &rit in &ri[row] {
            for &cit in &ci[col] { entries.push(Triplet::new(rit, cit, v)); }
        }
    }
    from_triplets(ym, yn, &entries)
}

/// `igl::slice(X, R, dim, Y)`: dim 1 selects rows `R`, dim 2 selects columns `R`.
pub fn slice(x: &SprsD, r: &[i32], dim: i32) -> SprsD {
    match dim {
        1 => {
            // boring base case
            if x.ncols() == 0 { return from_triplets(r.len(), 0, &[]); }
            let c: Vec<i32> = (0..x.ncols() as i32).collect(); // igl::colon(0, X.cols() - 1, C)
            slice_rc(x, r, &c)
        }
        2 => {
            // boring base case
            if x.nrows() == 0 { return from_triplets(0, r.len(), &[]); }
            let c: Vec<i32> = (0..x.nrows() as i32).collect();
            slice_rc(x, &c, r)
        }
        _ => panic!("Unsupported dimension"),
    }
}
