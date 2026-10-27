//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of Eigen 5.0.1 Eigen/src/OrderingMethods/Eigen_Colamd.h and the COLAMDOrdering functor of
//--- Eigen/src/OrderingMethods/Ordering.h (MPL-2.0). Eigen_Colamd.h is itself modified from COLAMD 2.3 by
//--- Stefan I. Larimore and Timothy A. Davis, Copyright (c) 1998-2003 by the University of Florida (see the
//--- notice reproduced at the end of this file).
//
// The index type is `i32`, Eigen's `StorageIndex` of `SparseMatrix<double>`, so that the arithmetic (dense
// thresholds, `INT_MAX - n_col`, the unsigned hash) is the one of the C++. The C code carves the `Col` and `Row`
// arrays out of the tail of `A`; here they are separate arrays and `A` keeps only the part the C code uses for the
// matrix (`Alen - Col_size - Row_size` entries), which is the only part it ever indexes.

const N_KNOBS: usize = 20;
const N_STATS: usize = 20;

// indices into knobs and stats
const DENSE_ROW: usize = 0;
const DENSE_COL: usize = 1;
const DEFRAG_COUNT: usize = 2;
const STATUS: usize = 3;
const INFO1: usize = 4;
const INFO2: usize = 5;
const INFO3: usize = 6;

// error codes returned in stats[3]
const OK: i32 = 0;
const OK_BUT_JUMBLED: i32 = 1;
const ERROR_NROW_NEGATIVE: i32 = -3;
const ERROR_NCOL_NEGATIVE: i32 = -4;
const ERROR_NNZ_NEGATIVE: i32 = -5;
const ERROR_P0_NONZERO: i32 = -6;
const ERROR_A_TOO_SMALL: i32 = -7;
const ERROR_COL_LENGTH_NEGATIVE: i32 = -8;
const ERROR_ROW_INDEX_OUT_OF_BOUNDS: i32 = -9;

#[inline]
fn ones_complement(r: i32) -> i32 { -r - 1 }

const EMPTY: i32 = -1;

// row and column status
const ALIVE: i32 = 0;
const DEAD: i32 = -1;

// column status
const DEAD_PRINCIPAL: i32 = -1;
const DEAD_NON_PRINCIPAL: i32 = -2;

/// `ColStructure`. The C unions are single fields:
/// shared1 = thickness (alive) | parent (dead), shared2 = score (alive) | order (dead),
/// shared3 = headhash | hash | prev, shared4 = degree_next | hash_next.
#[derive(Clone, Copy, Default)]
struct ColStructure {
    start: i32,  // index for A of first row in this column, or Dead if column is dead
    length: i32, // number of rows in this column
    shared1: i32,
    shared2: i32,
    shared3: i32,
    shared4: i32,
}

impl ColStructure {
    #[inline] fn is_dead(&self) -> bool { self.start < ALIVE }
    #[inline] fn is_alive(&self) -> bool { self.start >= ALIVE }
    #[inline] fn is_dead_principal(&self) -> bool { self.start == DEAD_PRINCIPAL }
    #[inline] fn kill_principal(&mut self) { self.start = DEAD_PRINCIPAL; }
    #[inline] fn kill_non_principal(&mut self) { self.start = DEAD_NON_PRINCIPAL; }
}

/// `RowStructure`. shared1 = degree | p, shared2 = mark | first_column.
#[derive(Clone, Copy, Default)]
struct RowStructure {
    start: i32,  // index for A of first col in this row
    length: i32, // number of principal columns in this row
    shared1: i32,
    shared2: i32,
}

impl RowStructure {
    #[inline] fn is_dead(&self) -> bool { self.shared2 < ALIVE }
    #[inline] fn is_alive(&self) -> bool { self.shared2 >= ALIVE }
    #[inline] fn kill(&mut self) { self.shared2 = DEAD; }
}

/// `colamd_c`: `(n_col + 1) * sizeof(ColStructure) / sizeof(IndexType)`, six indices per column.
#[inline]
fn colamd_c(n_col: i32) -> i32 { (n_col + 1) * 6 }

/// `colamd_r`: `(n_row + 1) * sizeof(RowStructure) / sizeof(IndexType)`, four indices per row.
#[inline]
fn colamd_r(n_row: i32) -> i32 { (n_row + 1) * 4 }

/// Returns the recommended value of Alen, -1 if any argument is negative.
fn recommended(nnz: i32, n_row: i32, n_col: i32) -> i32 {
    if nnz < 0 || n_row < 0 || n_col < 0 { -1 } else { 2 * nnz + colamd_c(n_col) + colamd_r(n_row) + n_col + nnz / 5 }
}

/// Sets default parameters: rows / columns over 50% dense are ignored.
fn set_defaults(knobs: &mut [f64; N_KNOBS]) {
    for k in knobs.iter_mut() { *k = 0.; }
    knobs[DENSE_ROW] = 0.5;
    knobs[DENSE_COL] = 0.5;
}

/// Computes a column ordering (Q) of A such that P(AQ)=LU or (AQ)'AQ=LL' have less fill-in. On return
/// `p[0..n_col]` is the permutation: `p[k]` is the column placed k-th.
fn compute_ordering(n_row: i32, n_col: i32, mut alen: i32, a: &mut Vec<i32>, p: &mut [i32], knobs: &[f64; N_KNOBS], stats: &mut [i32; N_STATS]) -> bool {
    for s in stats.iter_mut() { *s = 0; }
    stats[STATUS] = OK;
    stats[INFO1] = -1;
    stats[INFO2] = -1;

    if n_row < 0 { stats[STATUS] = ERROR_NROW_NEGATIVE; stats[INFO1] = n_row; return false; }
    if n_col < 0 { stats[STATUS] = ERROR_NCOL_NEGATIVE; stats[INFO1] = n_col; return false; }

    let nnz = p[n_col as usize];
    if nnz < 0 { stats[STATUS] = ERROR_NNZ_NEGATIVE; stats[INFO1] = nnz; return false; }
    if p[0] != 0 { stats[STATUS] = ERROR_P0_NONZERO; stats[INFO1] = p[0]; return false; }

    // === Allocate the Row and Col arrays from array A ===
    let col_size = colamd_c(n_col);
    let row_size = colamd_r(n_row);
    let need = 2 * nnz + n_col + col_size + row_size;
    if need > alen {
        stats[STATUS] = ERROR_A_TOO_SMALL;
        stats[INFO1] = need;
        stats[INFO2] = alen;
        return false;
    }
    alen -= col_size + row_size;
    a.truncate(alen as usize);
    let mut col = vec![ColStructure::default(); (n_col + 1) as usize];
    let mut row = vec![RowStructure::default(); (n_row + 1) as usize];

    // === Construct the row and column data structures ===
    if !init_rows_cols(n_row, n_col, &mut row, &mut col, a, p, stats) { return false; }

    // === Initialize scores, kill dense rows/columns ===
    let (n_row2, n_col2, max_deg) = init_scoring(n_row, n_col, &mut row, &mut col, a, p, knobs);

    // === Order the supercolumns ===
    let ngarbage = find_ordering(n_row, n_col, alen, &mut row, &mut col, a, p, n_col2, max_deg, 2 * nnz);

    // === Order the non-principal columns ===
    order_children(n_col, &mut col, p);

    // === Return statistics in stats ===
    stats[DENSE_ROW] = n_row - n_row2;
    stats[DENSE_COL] = n_col - n_col2;
    stats[DEFRAG_COUNT] = ngarbage;
    true
}

/// Takes the column form of the matrix in A and creates the row form of the matrix. Also, row and column
/// attributes are stored in the Col and Row structs. If the columns are un-sorted or contain duplicate row indices,
/// this routine will also sort and remove duplicate row indices from the column form of the matrix.
fn init_rows_cols(n_row: i32, n_col: i32, row: &mut [RowStructure], col: &mut [ColStructure], a: &mut [i32], p: &mut [i32], stats: &mut [i32; N_STATS]) -> bool {
    // === Initialize columns, and check column pointers ===
    for c in 0..n_col as usize {
        col[c].start = p[c];
        col[c].length = p[c + 1] - p[c];
        if col[c].length < 0 {
            stats[STATUS] = ERROR_COL_LENGTH_NEGATIVE;
            stats[INFO1] = c as i32;
            stats[INFO2] = col[c].length;
            return false;
        }
        col[c].shared1 = 1;     // thickness
        col[c].shared2 = 0;     // score
        col[c].shared3 = EMPTY; // prev
        col[c].shared4 = EMPTY; // degree_next
    }
    // p [0..n_col] no longer needed, used as "head" in subsequent routines

    // === Scan columns, compute row degrees, and check row indices ===
    stats[INFO3] = 0; // number of duplicate or unsorted row indices
    for r in 0..n_row as usize {
        row[r].length = 0;
        row[r].shared2 = -1; // mark
    }

    for c in 0..n_col {
        let mut last_row = -1;
        let mut cp = p[c as usize] as usize;
        let cp_end = p[c as usize + 1] as usize;
        while cp < cp_end {
            let r = a[cp];
            cp += 1;
            // make sure row indices within range
            if r < 0 || r >= n_row {
                stats[STATUS] = ERROR_ROW_INDEX_OUT_OF_BOUNDS;
                stats[INFO1] = c;
                stats[INFO2] = r;
                stats[INFO3] = n_row;
                return false;
            }
            if r <= last_row || row[r as usize].shared2 == c {
                // row index are unsorted or repeated (or both), thus col is jumbled. This is a notice, not an
                // error condition.
                stats[STATUS] = OK_BUT_JUMBLED;
                stats[INFO1] = c;
                stats[INFO2] = r;
                stats[INFO3] += 1;
            }
            if row[r as usize].shared2 != c {
                row[r as usize].length += 1;
            } else {
                // this is a repeated entry in the column, it will be removed
                col[c as usize].length -= 1;
            }
            // mark the row as having been seen in this column
            row[r as usize].shared2 = c;
            last_row = r;
        }
    }

    // === Compute row pointers ===
    // row form of the matrix starts directly after the column form of matrix in A
    row[0].start = p[n_col as usize];
    row[0].shared1 = row[0].start; // p
    row[0].shared2 = -1;           // mark
    for r in 1..n_row as usize {
        row[r].start = row[r - 1].start + row[r - 1].length;
        row[r].shared1 = row[r].start;
        row[r].shared2 = -1;
    }

    // === Create row form ===
    if stats[STATUS] == OK_BUT_JUMBLED {
        // if cols jumbled, watch for repeated row indices
        for c in 0..n_col {
            let mut cp = p[c as usize] as usize;
            let cp_end = p[c as usize + 1] as usize;
            while cp < cp_end {
                let r = a[cp] as usize;
                cp += 1;
                if row[r].shared2 != c {
                    a[row[r].shared1 as usize] = c;
                    row[r].shared1 += 1;
                    row[r].shared2 = c;
                }
            }
        }
    } else {
        // if cols not jumbled, we don't need the mark (this is faster)
        for c in 0..n_col {
            let mut cp = p[c as usize] as usize;
            let cp_end = p[c as usize + 1] as usize;
            while cp < cp_end {
                let r = a[cp] as usize;
                cp += 1;
                a[row[r].shared1 as usize] = c;
                row[r].shared1 += 1;
            }
        }
    }

    // === Clear the row marks and set row degrees ===
    for r in 0..n_row as usize {
        row[r].shared2 = 0;             // mark
        row[r].shared1 = row[r].length; // degree
    }

    // === See if we need to re-create columns ===
    if stats[STATUS] == OK_BUT_JUMBLED {
        // === Compute col pointers ===
        // col form of the matrix starts at A [0]. Note, we may have a gap between the col form and the row form if
        // there were duplicate entries, if so, it will be removed upon the first garbage collection
        col[0].start = 0;
        p[0] = col[0].start;
        for c in 1..n_col as usize {
            // note that the lengths here are for pruned columns, i.e. no duplicate row indices will exist for
            // these columns
            col[c].start = col[c - 1].start + col[c - 1].length;
            p[c] = col[c].start;
        }

        // === Re-create col form ===
        for r in 0..n_row {
            let mut rp = row[r as usize].start as usize;
            let rp_end = rp + row[r as usize].length as usize;
            while rp < rp_end {
                let c = a[rp] as usize;
                rp += 1;
                a[p[c] as usize] = r;
                p[c] += 1;
            }
        }
    }

    // === Done. Matrix is not (or no longer) jumbled ===
    true
}

/// Kills dense or empty columns and rows, calculates an initial score for each column, and places all columns in
/// the degree lists. Returns (n_row2, n_col2, max_deg).
fn init_scoring(n_row: i32, n_col: i32, row: &mut [RowStructure], col: &mut [ColStructure], a: &mut [i32], head: &mut [i32], knobs: &[f64; N_KNOBS]) -> (i32, i32, i32) {
    // === Extract knobs ===
    let dense_row_count = 0.max(((knobs[DENSE_ROW] * n_col as f64) as i32).min(n_col));
    let dense_col_count = 0.max(((knobs[DENSE_COL] * n_row as f64) as i32).min(n_row));
    let mut max_deg = 0;
    let mut n_col2 = n_col; // number of non-dense, non-empty columns
    let mut n_row2 = n_row; // number of non-dense, non-empty rows

    // === Kill empty columns ===
    // Put the empty columns at the end in their natural order, so that LU factorization can proceed as far as
    // possible.
    for c in (0..n_col as usize).rev() {
        let deg = col[c].length;
        if deg == 0 {
            // this is a empty column, kill and order it last
            n_col2 -= 1;
            col[c].shared2 = n_col2; // order
            col[c].kill_principal();
        }
    }

    // === Kill dense columns ===
    // Put the dense columns at the end, in their natural order
    for c in (0..n_col as usize).rev() {
        // skip any dead columns
        if col[c].is_dead() { continue; }
        let deg = col[c].length;
        if deg > dense_col_count {
            // this is a dense column, kill and order it last
            n_col2 -= 1;
            col[c].shared2 = n_col2;
            // decrement the row degrees
            let mut cp = col[c].start as usize;
            let cp_end = cp + col[c].length as usize;
            while cp < cp_end {
                row[a[cp] as usize].shared1 -= 1;
                cp += 1;
            }
            col[c].kill_principal();
        }
    }

    // === Kill dense and empty rows ===
    for r in 0..n_row as usize {
        let deg = row[r].shared1;
        if deg > dense_row_count || deg == 0 {
            // kill a dense or empty row
            row[r].kill();
            n_row2 -= 1;
        } else {
            // keep track of max degree of remaining rows
            max_deg = max_deg.max(deg);
        }
    }

    // === Compute initial column scores ===
    // At this point the row degrees are accurate. They reflect the number of "live" (non-dense) columns in each
    // row. No empty rows exist. Some "live" columns may contain only dead rows, however. These are pruned in the
    // code below.

    // now find the initial matlab score for each column
    for c in (0..n_col as usize).rev() {
        // skip dead column
        if col[c].is_dead() { continue; }
        let mut score = 0;
        let mut cp = col[c].start as usize;
        let mut new_cp = cp;
        let cp_end = cp + col[c].length as usize;
        while cp < cp_end {
            // get a row
            let r = a[cp];
            cp += 1;
            // skip if dead
            if row[r as usize].is_dead() { continue; }
            // compact the column
            a[new_cp] = r;
            new_cp += 1;
            // add row's external degree
            score += row[r as usize].shared1 - 1;
            // guard against integer overflow
            score = score.min(n_col);
        }
        // determine pruned column length
        let col_length = (new_cp - col[c].start as usize) as i32;
        if col_length == 0 {
            // a newly-made null column (all rows in this col are "dense" and have already been killed)
            n_col2 -= 1;
            col[c].shared2 = n_col2;
            col[c].kill_principal();
        } else {
            // set column length and set score
            col[c].length = col_length;
            col[c].shared2 = score;
        }
    }

    // At this point, all empty rows and columns are dead. All live columns are "clean" (containing no dead rows)
    // and simplicial (no supercolumns yet). Rows may contain dead columns, but all live rows contain at least one
    // live column.

    // === Initialize degree lists ===
    // clear the hash buckets
    for c in 0..=n_col as usize { head[c] = EMPTY; }
    let mut min_score = n_col;
    // place in reverse order, so low column indices are at the front of the lists. This is to encourage natural
    // tie-breaking
    for c in (0..n_col).rev() {
        // only add principal columns to degree lists
        if col[c as usize].is_alive() {
            // === Add columns score to DList ===
            let score = col[c as usize].shared2;
            // now add this column to dList at proper score location
            let next_col = head[score as usize];
            col[c as usize].shared3 = EMPTY;    // prev
            col[c as usize].shared4 = next_col; // degree_next
            // if there already was a column with the same score, set its previous pointer to this new column
            if next_col != EMPTY { col[next_col as usize].shared3 = c; }
            head[score as usize] = c;
            // see if this score is less than current min
            min_score = min_score.min(score);
        }
    }
    let _ = min_score;

    // === Return number of remaining columns, and max row degree ===
    (n_row2, n_col2, max_deg)
}

/// Order the principal columns of the supercolumn form of the matrix (no supercolumns on input). Uses a minimum
/// approximate column minimum degree ordering method. Returns the number of garbage collections.
fn find_ordering(n_row: i32, n_col: i32, alen: i32, row: &mut [RowStructure], col: &mut [ColStructure], a: &mut [i32], head: &mut [i32], n_col2: i32, mut max_deg: i32, mut pfree: i32) -> i32 {
    // === Initialization and clear mark ===
    let max_mark = i32::MAX - n_col; // INT_MAX defined in <limits.h>
    let mut tag_mark = clear_mark(n_row, row);
    let mut min_score = 0;
    let mut ngarbage = 0;

    // === Order the columns ===
    let mut k = 0;
    while k < n_col2 {
        // === Select pivot column, and order it ===
        // get pivot column from head of minimum degree list
        while min_score < n_col && head[min_score as usize] == EMPTY { min_score += 1; }
        let pivot_col = head[min_score as usize];
        let mut next_col = col[pivot_col as usize].shared4;
        head[min_score as usize] = next_col;
        if next_col != EMPTY { col[next_col as usize].shared3 = EMPTY; }

        // remember score for defrag check
        let pivot_col_score = col[pivot_col as usize].shared2;
        // the pivot column is the kth column in the pivot order
        col[pivot_col as usize].shared2 = k;
        // increment order count by column thickness
        let pivot_col_thickness = col[pivot_col as usize].shared1;
        k += pivot_col_thickness;

        // === Garbage_collection, if necessary ===
        let needed_memory = pivot_col_score.min(n_col - k);
        if pfree + needed_memory >= alen {
            pfree = garbage_collection(n_row, n_col, row, col, a, pfree);
            ngarbage += 1;
            // after garbage collection we will have enough. garbage collection has wiped out the Row[].shared2.mark
            // array
            tag_mark = clear_mark(n_row, row);
        }

        // === Compute pivot row pattern ===
        // get starting location for this new merged row
        let pivot_row_start = pfree;
        // initialize new row counts to zero
        let mut pivot_row_degree = 0;
        // tag pivot column as having been visited so it isn't included in merged pivot row
        col[pivot_col as usize].shared1 = -pivot_col_thickness;

        // pivot row is the union of all rows in the pivot column pattern
        let mut cp = col[pivot_col as usize].start as usize;
        let mut cp_end = cp + col[pivot_col as usize].length as usize;
        while cp < cp_end {
            // get a row
            let r = a[cp] as usize;
            cp += 1;
            // skip if row is dead
            if row[r].is_dead() { continue; }
            let mut rp = row[r].start as usize;
            let rp_end = rp + row[r].length as usize;
            while rp < rp_end {
                // get a column
                let c = a[rp];
                rp += 1;
                // add the column, if alive and untagged
                let col_thickness = col[c as usize].shared1;
                if col_thickness > 0 && col[c as usize].is_alive() {
                    // tag column in pivot row
                    col[c as usize].shared1 = -col_thickness;
                    // place column in pivot row
                    a[pfree as usize] = c;
                    pfree += 1;
                    pivot_row_degree += col_thickness;
                }
            }
        }

        // clear tag on pivot column
        col[pivot_col as usize].shared1 = pivot_col_thickness;
        max_deg = max_deg.max(pivot_row_degree);

        // === Kill all rows used to construct pivot row ===
        // also kill pivot row, temporarily
        cp = col[pivot_col as usize].start as usize;
        cp_end = cp + col[pivot_col as usize].length as usize;
        while cp < cp_end {
            // may be killing an already dead row
            let r = a[cp] as usize;
            cp += 1;
            row[r].kill();
        }

        // === Select a row index to use as the new pivot row ===
        let pivot_row_length = pfree - pivot_row_start;
        let pivot_row = if pivot_row_length > 0 {
            // pick the "pivot" row arbitrarily (first row in col)
            a[col[pivot_col as usize].start as usize]
        } else {
            // there is no pivot row, since it is of zero length
            EMPTY
        };

        // === Approximate degree computation ===
        // Here begins the computation of the approximate degree. The column score is the sum of the pivot row
        // "length", plus the size of the set differences of each row in the column minus the pattern of the pivot
        // row itself. The column ("thickness") itself is also excluded from the column score (we thus use an
        // approximate external degree).

        // === Compute set differences ===
        // pivot row is currently dead - it will be revived later.

        // for each column in pivot row
        let mut rp = pivot_row_start as usize;
        let mut rp_end = rp + pivot_row_length as usize;
        while rp < rp_end {
            let c = a[rp] as usize;
            rp += 1;
            // clear tags used to construct pivot row pattern
            let col_thickness = -col[c].shared1;
            col[c].shared1 = col_thickness;

            // === Remove column from degree list ===
            let cur_score = col[c].shared2;
            let prev_col = col[c].shared3;
            next_col = col[c].shared4;
            if prev_col == EMPTY {
                head[cur_score as usize] = next_col;
            } else {
                col[prev_col as usize].shared4 = next_col;
            }
            if next_col != EMPTY { col[next_col as usize].shared3 = prev_col; }

            // === Scan the column ===
            let mut cp = col[c].start as usize;
            let cp_end = cp + col[c].length as usize;
            while cp < cp_end {
                // get a row
                let r = a[cp] as usize;
                cp += 1;
                // skip if dead
                if row[r].is_dead() { continue; }
                let row_mark = row[r].shared2;
                let mut set_difference = row_mark - tag_mark;
                // check if the row has been seen yet
                if set_difference < 0 { set_difference = row[r].shared1; }
                // subtract column thickness from this row's set difference
                set_difference -= col_thickness;
                // absorb this row if the set difference becomes zero
                if set_difference == 0 {
                    row[r].kill();
                } else {
                    // save the new mark
                    row[r].shared2 = set_difference + tag_mark;
                }
            }
        }

        // === Add up set differences for each column ===
        // for each column in pivot row
        rp = pivot_row_start as usize;
        rp_end = rp + pivot_row_length as usize;
        while rp < rp_end {
            // get a column
            let c = a[rp] as usize;
            rp += 1;
            let mut hash: u32 = 0;
            let mut cur_score = 0;
            let mut cp = col[c].start as usize;
            // compact the column
            let mut new_cp = cp;
            let cp_end = cp + col[c].length as usize;
            while cp < cp_end {
                // get a row
                let r = a[cp];
                cp += 1;
                // skip if dead
                if row[r as usize].is_dead() { continue; }
                let row_mark = row[r as usize].shared2;
                // compact the column
                a[new_cp] = r;
                new_cp += 1;
                // compute hash function
                hash = hash.wrapping_add(r as u32);
                // add set difference
                cur_score += row_mark - tag_mark;
                // integer overflow...
                cur_score = cur_score.min(n_col);
            }

            // recompute the column's length
            col[c].length = (new_cp - col[c].start as usize) as i32;

            // === Further mass elimination ===
            if col[c].length == 0 {
                // nothing left but the pivot row in this column
                col[c].kill_principal();
                pivot_row_degree -= col[c].shared1;
                // order it
                col[c].shared2 = k;
                // increment order count by column thickness
                k += col[c].shared1;
            } else {
                // === Prepare for supercolumn detection ===
                // save score so far
                col[c].shared2 = cur_score;
                // add column to hash table, for supercolumn detection
                hash %= (n_col + 1) as u32;
                let head_column = head[hash as usize];
                let first_col;
                if head_column > EMPTY {
                    // degree list "hash" is non-empty, use prev (shared3) of first column in degree list as head
                    // of hash bucket
                    first_col = col[head_column as usize].shared3;
                    col[head_column as usize].shared3 = c as i32;
                } else {
                    // degree list "hash" is empty, use head as hash bucket
                    first_col = -(head_column + 2);
                    head[hash as usize] = -(c as i32 + 2);
                }
                col[c].shared4 = first_col;
                // save hash function in Col [col].shared3.hash
                col[c].shared3 = hash as i32;
            }
        }

        // The approximate external column degree is now computed.

        // === Supercolumn detection ===
        detect_super_cols(col, a, head, pivot_row_start, pivot_row_length);

        // === Kill the pivotal column ===
        col[pivot_col as usize].kill_principal();

        // === Clear mark ===
        tag_mark += max_deg + 1;
        if tag_mark >= max_mark { tag_mark = clear_mark(n_row, row); }

        // === Finalize the new pivot row, and column scores ===
        // for each column in pivot row
        rp = pivot_row_start as usize;
        // compact the pivot row
        let mut new_rp = rp;
        rp_end = rp + pivot_row_length as usize;
        while rp < rp_end {
            let c = a[rp];
            rp += 1;
            // skip dead columns
            if col[c as usize].is_dead() { continue; }
            a[new_rp] = c;
            new_rp += 1;
            // add new pivot row to column
            let at = (col[c as usize].start + col[c as usize].length) as usize;
            col[c as usize].length += 1;
            a[at] = pivot_row;

            // retrieve score so far and add on pivot row's degree. (we wait until here for this in case the pivot
            // row's degree was reduced due to mass elimination).
            let mut cur_score = col[c as usize].shared2 + pivot_row_degree;
            // calculate the max possible score as the number of external columns minus the 'k' value minus the
            // columns thickness
            let max_score = n_col - k - col[c as usize].shared1;
            // make the score the external degree of the union-of-rows
            cur_score -= col[c as usize].shared1;
            // make sure score is less or equal than the max score
            cur_score = cur_score.min(max_score);
            // store updated score
            col[c as usize].shared2 = cur_score;

            // === Place column back in degree list ===
            next_col = head[cur_score as usize];
            col[c as usize].shared4 = next_col;
            col[c as usize].shared3 = EMPTY;
            if next_col != EMPTY { col[next_col as usize].shared3 = c; }
            head[cur_score as usize] = c;

            // see if this score is less than current min
            min_score = min_score.min(cur_score);
        }

        // === Resurrect the new pivot row ===
        if pivot_row_degree > 0 {
            // update pivot row length to reflect any cols that were killed during super-col detection and mass
            // elimination
            let pr = pivot_row as usize;
            row[pr].start = pivot_row_start;
            row[pr].length = (new_rp - pivot_row_start as usize) as i32;
            row[pr].shared1 = pivot_row_degree;
            row[pr].shared2 = 0;
            // pivot row is no longer dead
        }
    }

    // === All principal columns have now been ordered ===
    ngarbage
}

/// Orders the non-principal columns by walking up the parent tree (a column is a child of the column which
/// absorbed it). The final permutation vector is then placed in p [0 ... n_col-1], with p [0] being the first
/// column, and p [n_col-1] being the last.
fn order_children(n_col: i32, col: &mut [ColStructure], p: &mut [i32]) {
    // === Order each non-principal column ===
    for i in 0..n_col as usize {
        // find an un-ordered non-principal column
        if !col[i].is_dead_principal() && col[i].shared2 == EMPTY {
            let mut parent = i as i32;
            // once found, find its principal parent
            loop {
                parent = col[parent as usize].shared1;
                if col[parent as usize].is_dead_principal() { break; }
            }

            // now, order all un-ordered non-principal columns along path to this parent. collapse tree at the same
            // time
            let mut c = i as i32;
            // get order of parent
            let mut order = col[parent as usize].shared2;

            loop {
                // order this column
                col[c as usize].shared2 = order;
                order += 1;
                // collapse tree
                col[c as usize].shared1 = parent;
                // get immediate parent of this column
                c = col[c as usize].shared1;
                // continue until we hit an ordered column. There are guaranteed not to be anymore unordered columns
                // above an ordered column
                if col[c as usize].shared2 != EMPTY { break; }
            }

            // re-order the super_col parent to largest order for this group
            col[parent as usize].shared2 = order;
        }
    }

    // === Generate the permutation ===
    for c in 0..n_col {
        p[col[c as usize].shared2 as usize] = c;
    }
}

/// Detects supercolumns by finding matches between columns in the hash buckets. Check amongst columns in the set
/// A [row_start ... row_start + row_length-1]. The columns under consideration are currently *not* in the degree
/// lists, and have already been placed in the hash buckets.
///
/// The hash bucket for columns whose hash function is equal to h is stored as follows: if head [h] is >= 0, then
/// head [h] contains a degree list, so head [h] is the first column in degree bucket h, and
/// Col [head [h]].headhash gives the first column in hash bucket h. Otherwise, the degree list is empty, and
/// -(head [h] + 2) is the first column in hash bucket h. For a column c in a hash bucket, Col [c].shared3.prev is
/// NOT a "previous column" pointer. Col [c].shared3.hash is used instead as the hash number for that column. The
/// value of Col [c].shared4.hash_next is the next column in the same hash bucket.
fn detect_super_cols(col: &mut [ColStructure], a: &[i32], head: &mut [i32], row_start: i32, row_length: i32) {
    // === Consider each column in the row ===
    let mut rp = row_start as usize;
    let rp_end = rp + row_length as usize;
    while rp < rp_end {
        let c0 = a[rp] as usize;
        rp += 1;
        if col[c0].is_dead() { continue; }

        // get hash number for this column
        let hash = col[c0].shared3;

        // === Get the first column in this hash bucket ===
        let head_column = head[hash as usize];
        let first_col = if head_column > EMPTY { col[head_column as usize].shared3 } else { -(head_column + 2) };

        // === Consider each column in the hash bucket ===
        let mut super_c = first_col;
        while super_c != EMPTY {
            let su = super_c as usize;
            let length = col[su].length;
            // prev_c is the column preceding column c in the hash bucket
            let mut prev_c = super_c;

            // === Compare super_c with all columns after it ===
            let mut c = col[su].shared4;
            while c != EMPTY {
                let cu = c as usize;
                // not identical if lengths or scores are different
                if col[cu].length != length || col[cu].shared2 != col[su].shared2 {
                    prev_c = c;
                    c = col[cu].shared4;
                    continue;
                }

                // compare the two columns
                let mut cp1 = col[su].start as usize;
                let mut cp2 = col[cu].start as usize;
                let mut i = 0;
                while i < length {
                    // the columns are "clean" (no dead rows)
                    // row indices will same order for both supercols, no gather scatter necessary
                    let differ = a[cp1] != a[cp2];
                    cp1 += 1;
                    cp2 += 1;
                    if differ { break; }
                    i += 1;
                }

                // the two columns are different if the for-loop "broke"
                if i != length {
                    prev_c = c;
                    c = col[cu].shared4;
                    continue;
                }

                // === Got it! two columns are identical ===
                col[su].shared1 += col[cu].shared1;
                col[cu].shared1 = super_c; // parent
                col[cu].kill_non_principal();
                // order c later, in order_children()
                col[cu].shared2 = EMPTY;
                // remove c from hash bucket
                col[prev_c as usize].shared4 = col[cu].shared4;
                c = col[cu].shared4;
            }
            super_c = col[su].shared4;
        }

        // === Empty this hash bucket ===
        if head_column > EMPTY {
            // corresponding degree list "hash" is not empty
            col[head_column as usize].shared3 = EMPTY;
        } else {
            // corresponding degree list "hash" is empty
            head[hash as usize] = EMPTY;
        }
    }
}

/// Defragments and compacts columns and rows in the workspace A. Used when all available memory has been used
/// while performing row merging. Returns the index of the first free position in A, after garbage collection.
fn garbage_collection(n_row: i32, n_col: i32, row: &mut [RowStructure], col: &mut [ColStructure], a: &mut [i32], pfree: i32) -> i32 {
    // === Defragment the columns ===
    let mut pdest = 0usize;
    for c in 0..n_col as usize {
        if col[c].is_alive() {
            let mut psrc = col[c].start as usize;
            // move and compact the column
            col[c].start = pdest as i32;
            let length = col[c].length;
            for _ in 0..length {
                let r = a[psrc];
                psrc += 1;
                if row[r as usize].is_alive() {
                    a[pdest] = r;
                    pdest += 1;
                }
            }
            col[c].length = (pdest - col[c].start as usize) as i32;
        }
    }

    // === Prepare to defragment the rows ===
    for r in 0..n_row {
        let ru = r as usize;
        if row[ru].is_alive() {
            if row[ru].length == 0 {
                // this row is of zero length. cannot compact it, so kill it
                row[ru].kill();
            } else {
                // save first column index in Row [r].shared2.first_column
                let psrc = row[ru].start as usize;
                row[ru].shared2 = a[psrc];
                // flag the start of the row with the one's complement of row
                a[psrc] = ones_complement(r);
            }
        }
    }

    // === Defragment the rows ===
    let mut psrc = pdest;
    while psrc < pfree as usize {
        // find a negative number ... the start of a row
        let v = a[psrc];
        psrc += 1;
        if v < 0 {
            psrc -= 1;
            // get the row index
            let r = ones_complement(a[psrc]) as usize;
            // restore first column index
            a[psrc] = row[r].shared2;
            // move and compact the row
            row[r].start = pdest as i32;
            let length = row[r].length;
            for _ in 0..length {
                let c = a[psrc];
                psrc += 1;
                if col[c as usize].is_alive() {
                    a[pdest] = c;
                    pdest += 1;
                }
            }
            row[r].length = (pdest - row[r].start as usize) as i32;
        }
    }

    // === Return the new value of pfree ===
    pdest as i32
}

/// Clears the Row [].shared2.mark array, and returns the new tag_mark.
fn clear_mark(n_row: i32, row: &mut [RowStructure]) -> i32 {
    for r in 0..n_row as usize {
        if row[r].is_alive() { row[r].shared2 = 0; }
    }
    1
}

/// `COLAMDOrdering<int>::operator()`: the column approximate minimum degree ordering of a column-major,
/// compressed sparse pattern. Returns Eigen's `perm.indices()`: `perm[old_col] = new_pos`.
pub fn colamd_ordering(nrows: usize, ncols: usize, col_ptr: &[usize], row_idx: &[usize]) -> Vec<usize> {
    let m = nrows as i32;
    let n = ncols as i32;
    let nnz = col_ptr[ncols] as i32;
    // Get the recommended value of Alen to be used by colamd
    let alen = recommended(nnz, m, n);
    // Set the default parameters
    let mut knobs = [0.; N_KNOBS];
    let mut stats = [0; N_STATS];
    set_defaults(&mut knobs);

    let mut p = vec![0i32; ncols + 1];
    let mut a = vec![0i32; alen as usize];
    for i in 0..=ncols { p[i] = col_ptr[i] as i32; }
    for i in 0..nnz as usize { a[i] = row_idx[i] as i32; }
    // Call Colamd routine to compute the ordering
    let info = compute_ordering(m, n, alen, &mut a, &mut p, &knobs, &mut stats);
    assert!(info, "COLAMD failed ");

    let mut perm = vec![0usize; ncols];
    for i in 0..ncols { perm[p[i] as usize] = i; }
    perm
}

// Notice of the original COLAMD (as reproduced in Eigen_Colamd.h):
//
//   The authors of the code itself are Stefan I. Larimore and Timothy A. Davis (davis@cise.ufl.edu), University of
//   Florida. The algorithm was developed in collaboration with John Gilbert, Xerox PARC, and Esmond Ng, Oak Ridge
//   National Laboratory.
//
//   Date: September 8, 2003. Version 2.3.
//
//   Acknowledgements: This work was supported by the National Science Foundation, under grants DMS-9504974 and
//   DMS-9803599.
//
//   Notice: Copyright (c) 1998-2003 by the University of Florida. All Rights Reserved.
//
//   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
//
//   Permission is hereby granted to use, copy, modify, and/or distribute this program, provided that the Copyright,
//   this License, and the Availability of the original version is retained on all copies and made accessible to the
//   end-user of any code or package that includes COLAMD or any modified version of COLAMD.
//
//   Availability: The colamd/symamd library is available at http://www.suitesparse.com

#[cfg(test)]
mod tests {
    use super::*;

    // Orderings of the comparison against Eigen 5.0.1's COLAMDOrdering<int> (every one of 444 random patterns was
    // identical).
    const CO0: (usize, usize, &[usize], &[usize], &str) = (6, 10, &[0, 3, 5, 6, 7, 9, 10, 11, 14, 15, 16], &[0, 2, 5, 1,
        5, 1, 4, 0, 5, 0, 4, 3, 4, 5, 0, 0], "5 4 0 1 6 7 2 3 8 9");
    const CO1: (usize, usize, &[usize], &[usize], &str) = (10, 10, &[0, 2, 3, 4, 5, 8, 12, 14, 14, 16, 18], &[4, 7, 0, 8,
        3, 3, 4, 6, 1, 2, 3, 9, 2, 5, 5, 9, 3, 9], "2 0 1 7 8 5 3 9 4 6");

    #[test]
    fn matches_eigen_0() {
        let (m, n, cp, ri, want) = CO0;
        let perm = colamd_ordering(m, n, cp, ri);
        assert_eq!(perm.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "), want);
    }

    #[test]
    fn matches_eigen_1() {
        let (m, n, cp, ri, want) = CO1;
        let perm = colamd_ordering(m, n, cp, ri);
        assert_eq!(perm.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "), want);
    }
}
