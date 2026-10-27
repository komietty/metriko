//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of Eigen 5.0.1 Eigen/src/SparseCore/SparseColEtree.h (MPL-2.0), itself a modified version of
//--- sp_coletree.c of SuperLU.
//
// Only `etree_find` and `coletree` are ported: SparseQR does not use `nr_etdfs` / `treePostorder`.
// Indices are `i32` as Eigen's `StorageIndex` of `SparseMatrix<double>`.

/// A column-major sparse pattern seen through Eigen's `InnerIterator`: column `j` holds the rows
/// `row_idx[col_beg[j] .. col_beg[j] + col_len[j]]`, in storage order (this covers both compressed matrices and the
/// uncompressed, lazily column-permuted `m_pmat` of SparseQR).
pub struct ColPattern<'a> {
    pub nrows: usize,
    pub ncols: usize,
    pub col_beg: &'a [usize],
    pub col_len: &'a [usize],
    pub row_idx: &'a [usize],
}

impl ColPattern<'_> {
    #[inline]
    pub fn col(&self, j: usize) -> &[usize] {
        let b = self.col_beg[j];
        &self.row_idx[b..b + self.col_len[j]]
    }
}

/// Find the root of the tree/set containing the vertex i: use path halving.
fn etree_find(mut i: i32, pp: &mut [i32]) -> i32 {
    let mut p = pp[i as usize]; // parent
    let mut gp = pp[p as usize]; // grand parent
    while gp != p {
        pp[i as usize] = gp; // parent pointer on find path is changed to former grand parent
        i = gp;
        p = pp[i as usize];
        gp = pp[p as usize];
    }
    p
}

/// Compute the column elimination tree of a sparse matrix.
/// - `mat`: the matrix in column-major format
/// - `parent`: the elimination tree
/// - `first_row_elt`: the column index of the first element in each row
/// - `perm`: the permutation to apply to the column of `mat` (column `col` of the permuted matrix is column
///   `perm[col]` of `mat`)
pub fn coletree(mat: &ColPattern, parent: &mut Vec<i32>, first_row_elt: &mut Vec<i32>, perm: Option<&[i32]>) -> i32 {
    let nc = mat.ncols as i32; // number of columns
    let m = mat.nrows as i32;
    let diag_size = nc.min(m);
    let mut root = vec![0i32; nc as usize]; // root of subtree of etree
    let mut pp = vec![0i32; nc as usize]; // disjoint sets

    parent.clear();
    parent.resize(nc as usize, 0);
    // Compute first nonzero column in each row
    first_row_elt.clear();
    first_row_elt.resize(m as usize, nc);
    for i in 0..diag_size { first_row_elt[i as usize] = i; }
    for col in 0..nc {
        let pcol = match perm { Some(p) => p[col as usize], None => col };
        for &row in mat.col(pcol as usize) {
            first_row_elt[row] = first_row_elt[row].min(col);
        }
    }
    // Compute etree by Liu's algorithm for symmetric matrices, except use (firstRowElt[r],c) in place of an edge
    // (r,c) of A. Thus each row clique in A'*A is replaced by a star centered at its first vertex, which has the same
    // fill.
    for col in 0..nc {
        let mut found_diag = col >= m;
        pp[col as usize] = col;
        let mut cset = col;
        root[cset as usize] = col;
        parent[col as usize] = nc;
        // The diagonal element is treated here even if it does not exist in the matrix hence the loop is executed
        // once more
        let pcol = match perm { Some(p) => p[col as usize], None => col };
        let rows = mat.col(pcol as usize);
        let mut it = 0usize;
        while it < rows.len() || !found_diag {
            // A sequence of interleaved find and union is performed
            let i = if it < rows.len() { rows[it] as i32 } else { col };
            it += 1;
            if i == col { found_diag = true; }

            let row = first_row_elt[i as usize];
            if row >= col { continue; }
            let rset = etree_find(row, &mut pp); // Find the name of the set containing row
            let rroot = root[rset as usize];
            if rroot != col {
                parent[rroot as usize] = col;
                pp[cset as usize] = rset;
                cset = rset;
                root[cset as usize] = col;
            }
        }
    }
    0
}
