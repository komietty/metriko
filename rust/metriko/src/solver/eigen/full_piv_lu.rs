//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of Eigen 5.0.1 Eigen/src/LU/FullPivLU.h (MPL-2.0): the decomposition, rank, kernel and solve that
//--- include/metriko/solver/matrix_ops.h uses. The dense kernels these call are ported with the arithmetic Eigen
//--- runs on the reference platform (Apple arm64, NEON with FMA, -ffp-contract=off):
//---   Eigen/src/Core/products/TriangularSolverVector.h  (triangular_solve_vector, col-major, OnTheLeft)
//---   Eigen/src/Core/products/GeneralMatrixVector.h     (general_matrix_vector_product, col-major)
//---   Eigen/src/Core/products/TriangularSolverMatrix.h  (triangular_solve_matrix, OnTheLeft, col-major)
//---   Eigen/src/Core/products/GeneralBlockPanelKernel.h (gebp_kernel: the per element accumulation of each path,
//---                                                      and evaluateProductBlockingSizesHeuristic)
//
// Arithmetic of the reference: Eigen's packet `pmadd` is fused (vfmaq_f64), its scalar `madd` is not
// (EIGEN_SCALAR_MADD_USE_FMA is 0 there), and plain expressions are never contracted. `f64::mul_add` reproduces the
// fused operations, plain `a * b + c` the others.
//
// DEVIATION: Eigen sizes the blocks of its level-3 kernels from the cache sizes it queries at run time. They are
// fixed here to the values Eigen reports on the reference platform (l1 = 16 KiB, l2 = l3 = 512 KiB), so the kernel
// of matrices of rank >= 48 is computed with the blocking of the reference machine on every machine.

use faer::mat::AsMatRef;
use faer::{Col, Mat, MatRef};

/// Cache sizes `manage_caching_sizes(GetAction, ...)` reports on the reference platform (see DEVIATION above).
const L1: i64 = 16 * 1024;
const L2: i64 = 512 * 1024;
const L3: i64 = 512 * 1024;
/// gebp_traits<double, double> on NEON.
const GEBP_MR: i64 = 6;
const GEBP_NR: i64 = 8;
/// EIGEN_TUNE_TRIANGULAR_PANEL_WIDTH.
const TRIANGULAR_PANEL_WIDTH: usize = 8;

/// LU decomposition of a matrix with complete pivoting, `P A Q = L U`.
pub struct FullPivLu {
    lu: Mat<f64>,
    p: Vec<usize>,
    q: Vec<usize>,
    rows_transpositions: Vec<usize>,
    cols_transpositions: Vec<usize>,
    nonzero_pivots: usize,
    l1_norm: f64,
    maxpivot: f64,
    det_pq: i32,
}

impl FullPivLu {
    /// Computes the decomposition of `matrix` (`FullPivLU(const EigenBase&)`), given as a `Mat`, a `&Mat` or a
    /// `MatRef`.
    pub fn new(matrix: impl AsMatRef<T = f64, Rows = usize, Cols = usize>) -> Self {
        let matrix = matrix.as_mat_ref();
        let rows = matrix.nrows();
        let cols = matrix.ncols();
        let mut this = Self {
            lu: matrix.to_owned(),
            p: vec![0; rows],
            q: vec![0; cols],
            rows_transpositions: vec![0; rows],
            cols_transpositions: vec![0; cols],
            nonzero_pivots: 0,
            l1_norm: 0.,
            maxpivot: 0.,
            det_pq: 1,
        };
        this.compute_in_place();
        this
    }

    fn compute_in_place(&mut self) {
        let rows = self.lu.nrows();
        let cols = self.lu.ncols();
        let size = rows.min(cols);
        let lu = &mut self.lu;

        // m_l1_norm = m_lu.cwiseAbs().colwise().sum().maxCoeff(); only used by rcond(), which metriko does not call
        self.l1_norm = (0..cols).map(|j| (0..rows).fold(0., |s, i| s + lu[(i, j)].abs())).fold(0., f64::max);

        // will store the transpositions, before we accumulate them at the end.
        // can't accumulate on-the-fly because that will be done in reverse order for the rows.
        let mut number_of_transpositions = 0; // number of NONTRIVIAL transpositions, i.e. m_rowsTranspositions[i]!=i

        self.nonzero_pivots = size; // the generic case is that in which all pivots are nonzero (invertible case)
        self.maxpivot = 0.;

        for k in 0..size {
            // First, we need to find the pivot.

            // biggest coefficient in the remaining bottom-right corner (starting at row k, col k): the first maximum
            // of |a| in column-major order (FindCoeff.h replaces the incumbent only on a strictly larger candidate,
            // its vectorized loop gives the same position)
            let mut row_of_biggest_in_corner = k;
            let mut col_of_biggest_in_corner = k;
            let mut biggest_in_corner = lu[(k, k)].abs();
            for j in k..cols {
                for i in k..rows {
                    let s = lu[(i, j)].abs();
                    if s > biggest_in_corner {
                        biggest_in_corner = s;
                        row_of_biggest_in_corner = i;
                        col_of_biggest_in_corner = j;
                    }
                }
            }

            if biggest_in_corner == 0. {
                // before exiting, make sure to initialize the still uninitialized transpositions in a sane state
                // without destroying what we already have.
                self.nonzero_pivots = k;
                for i in k..size {
                    self.rows_transpositions[i] = i;
                    self.cols_transpositions[i] = i;
                }
                break;
            }

            let abs_pivot = biggest_in_corner;
            if abs_pivot > self.maxpivot { self.maxpivot = abs_pivot; }

            // Now that we've found the pivot, we need to apply the row/col swaps to bring it to the location (k,k).
            self.rows_transpositions[k] = row_of_biggest_in_corner;
            self.cols_transpositions[k] = col_of_biggest_in_corner;
            if k != row_of_biggest_in_corner {
                for j in 0..cols {
                    let t = lu[(k, j)];
                    lu[(k, j)] = lu[(row_of_biggest_in_corner, j)];
                    lu[(row_of_biggest_in_corner, j)] = t;
                }
                number_of_transpositions += 1;
            }
            if k != col_of_biggest_in_corner {
                for i in 0..rows {
                    let t = lu[(i, k)];
                    lu[(i, k)] = lu[(i, col_of_biggest_in_corner)];
                    lu[(i, col_of_biggest_in_corner)] = t;
                }
                number_of_transpositions += 1;
            }

            // Now that the pivot is at the right location, we update the remaining bottom-right corner by Gaussian
            // elimination.
            if k < rows - 1 {
                let d = lu[(k, k)];
                for i in k + 1..rows { lu[(i, k)] /= d; }
            }
            if k < size - 1 {
                // the outer product col(k).tail * row(k).tail (outer_product_selector_run): column j of the block
                // gets `a - u_j * l_i`, product then difference
                for j in k + 1..cols {
                    let u = lu[(k, j)];
                    for i in k + 1..rows {
                        let l = lu[(i, k)];
                        lu[(i, j)] -= u * l;
                    }
                }
            }
        }

        // the main loop is over, we still have to accumulate the transpositions to find the permutations P and Q
        for i in 0..rows { self.p[i] = i; }
        for k in (0..size).rev() { self.p.swap(k, self.rows_transpositions[k]); }

        for i in 0..cols { self.q[i] = i; }
        for k in 0..size { self.q.swap(k, self.cols_transpositions[k]); }

        self.det_pq = if number_of_transpositions % 2 == 1 { -1 } else { 1 };
    }

    /// The LU decomposition matrix: the upper-triangular part is U, the unit-lower-triangular part is L (at least
    /// for square matrices; in the non-square case, special care is needed, see the Eigen documentation).
    pub fn matrix_lu(&self) -> MatRef<'_, f64> { self.lu.as_ref() }

    /// The number of nonzero pivots in the LU decomposition, nonzero in the exact sense.
    pub fn nonzero_pivots(&self) -> usize { self.nonzero_pivots }

    /// The absolute value of the biggest pivot, i.e. the biggest diagonal coefficient of U.
    pub fn max_pivot(&self) -> f64 { self.maxpivot }

    /// `permutationP().indices()`.
    pub fn permutation_p(&self) -> &[usize] { &self.p }

    /// `permutationQ().indices()`.
    pub fn permutation_q(&self) -> &[usize] { &self.q }

    /// The determinant of the permutations P and Q (`m_det_pq`).
    pub fn det_pq(&self) -> i32 { self.det_pq }

    /// The L1 norm of the decomposed matrix (`m_l1_norm`).
    pub fn l1_norm(&self) -> f64 { self.l1_norm }

    /// The threshold used by rank() and kernel(): `epsilon * diagonalSize`.
    pub fn threshold(&self) -> f64 {
        // this formula comes from experimenting (see "LU precision tuning" thread on the list) and turns out to be
        // identical to Higham's formula used already in LDLt.
        f64::EPSILON * self.lu.nrows().min(self.lu.ncols()) as f64
    }

    /// The rank of the decomposed matrix: the pivots larger than `|maxpivot| * threshold`.
    pub fn rank(&self) -> usize {
        let premultiplied_threshold = self.maxpivot.abs() * self.threshold();
        let mut result = 0;
        for i in 0..self.nonzero_pivots {
            result += (self.lu[(i, i)].abs() > premultiplied_threshold) as usize;
        }
        result
    }

    /// The dimension of the kernel of the decomposed matrix.
    pub fn dimension_of_kernel(&self) -> usize { self.lu.ncols() - self.rank() }

    /// A matrix whose columns form a basis of the kernel (`kernel_retval::evalTo`). When the kernel is {0}, a
    /// single zero column.
    pub fn kernel(&self) -> Mat<f64> {
        let cols = self.lu.ncols();
        let rank = self.rank();
        let dimker = cols - rank;
        if dimker == 0 {
            // The Kernel is just {0}, so it doesn't have a basis properly speaking, but let's avoid
            // crashing/asserting as that depends on floating point calculations. Let's just return a single column
            // vector filled with zeros.
            return Mat::zeros(cols, 1);
        }

        // Lemma: If the matrix A has the LU decomposition PAQ = LU, then Ker A = Q(Ker U).
        // Thus, all we need to do is to compute Ker U, and then apply Q.
        //
        // U is upper triangular, with eigenvalues sorted so that any zeros appear at the end. Thus, the diagonal of
        // U ends with exactly dimKer zero's. Let us use that to construct dimKer linearly independent vectors in
        // Ker U.

        let premultiplied_threshold = self.maxpivot * self.threshold();
        let mut pivots = Vec::with_capacity(rank);
        for i in 0..self.nonzero_pivots {
            if self.lu[(i, i)].abs() > premultiplied_threshold { pivots.push(i); }
        }
        assert_eq!(pivots.len(), rank);

        // we construct a temporaty trapezoid matrix m, by taking the U matrix and permuting the rows and cols to
        // bring the nonnegligible pivots to the top of the main diagonal. We need that to be able to apply our
        // triangular solvers.
        let mut m = Mat::<f64>::from_fn(rank, cols, |i, j| self.lu[(i, j)]);
        for i in 0..rank {
            for j in 0..i { m[(i, j)] = 0.; }
            for j in i..cols { m[(i, j)] = self.lu[(pivots[i], j)]; }
        }
        for j in 0..rank {
            for i in j + 1..rank { m[(i, j)] = 0.; }
        }
        for i in 0..rank { swap_cols(&mut m, i, pivots[i]); }

        // ok, we have our trapezoid matrix, we can apply the triangular solver. notice that the math behind this
        // suggests that we should apply this to the negative of the RHS, but for performance we just put the
        // negative sign elsewhere, see below.
        {
            let (data, stride) = col_major_data_mut(&mut m);
            let (tri, other) = data.split_at_mut(rank * stride);
            trsm_left_upper(rank, dimker, tri, stride, other, stride);
        }

        // now we must undo the column permutation that we had applied!
        for i in (0..rank).rev() { swap_cols(&mut m, i, pivots[i]); }

        // see the negative sign in the next line, that's what we were talking about above.
        let mut dst = Mat::<f64>::zeros(cols, dimker);
        for i in 0..rank {
            for k in 0..dimker { dst[(self.q[i], k)] = -m[(i, rank + k)]; }
        }
        for i in rank..cols {
            for k in 0..dimker { dst[(self.q[i], k)] = 0.; }
        }
        for k in 0..dimker { dst[(self.q[rank + k], k)] = 1.; }
        dst
    }

    /// A solution x of `A x = rhs` (`_solve_impl`), the free variables set to zero. When no solution exists, a
    /// least-squares-like answer as Eigen computes it.
    pub fn solve(&self, rhs: &Col<f64>) -> Col<f64> {
        // The decomposition PAQ = LU can be rewritten as A = P^{-1} L U Q^{-1}. So we proceed as follows:
        // Step 1: compute c = P * rhs.
        // Step 2: replace c by the solution x to Lx = c. Exists because L is invertible.
        // Step 3: replace c by the solution x to Ux = c. May or may not exist.
        // Step 4: result = Q * c;
        let rows = self.lu.nrows();
        let cols = self.lu.ncols();
        let nonzero_pivots = self.rank();
        let smalldim = rows.min(cols);
        assert_eq!(rhs.nrows(), rows);

        let mut dst = Col::<f64>::zeros(cols);
        if nonzero_pivots == 0 { return dst; }

        // Step 1
        let mut c = vec![0.; rows];
        for i in 0..rows { c[self.p[i]] = rhs[i]; }

        // Step 2
        let lu = self.lu.as_ref();
        trsv_col_major(smalldim, |i, j| lu[(i, j)], &mut c[..smalldim], true, true);
        if rows > cols {
            let (top, bottom) = c.split_at_mut(cols);
            if rows - cols == 1 {
                // runtime 1 x n times n x 1: generic_product_impl falls back to the inner product
                bottom[0] += -1. * inner_product(cols, |k| lu[(cols, k)], |k| top[k]);
            } else {
                gemv_col_major(rows - cols, cols, rows, |i, j| lu[(cols + i, j)], |j| top[j], bottom, -1.);
            }
        }

        // Step 3
        trsv_col_major(nonzero_pivots, |i, j| lu[(i, j)], &mut c[..nonzero_pivots], false, false);

        // Step 4
        for i in 0..nonzero_pivots { dst[self.q[i]] = c[i]; }
        for i in nonzero_pivots..cols { dst[self.q[i]] = 0.; }
        dst
    }
}

fn swap_cols(m: &mut Mat<f64>, a: usize, b: usize) {
    if a == b { return; }
    for i in 0..m.nrows() {
        let t = m[(i, a)];
        m[(i, a)] = m[(i, b)];
        m[(i, b)] = t;
    }
}

/// The column-major storage of `m`, with its column stride.
fn col_major_data_mut(m: &mut Mat<f64>) -> (&mut [f64], usize) {
    let (nrows, ncols) = (m.nrows(), m.ncols());
    let stride = m.col_stride() as usize;
    assert!(m.row_stride() == 1 && (ncols <= 1 || stride >= nrows));
    let len = if ncols == 0 { 0 } else { stride * (ncols - 1) + nrows };
    // SAFETY: a faer Mat stores its columns at `ptr + j * col_stride` with unit row stride; the slice covers exactly
    // that storage and is borrowed mutably from `m`.
    let data = unsafe { std::slice::from_raw_parts_mut(m.as_ptr_mut(), len) };
    (data, stride)
}

/// The inner product `lhs.row(0).dot(rhs.col(0))` of InnerProduct.h, non-vectorized: `a0 * b0`, then
/// `(a_k * b_k) + result` in order.
fn inner_product(n: usize, a: impl Fn(usize) -> f64, b: impl Fn(usize) -> f64) -> f64 {
    if n == 0 { return 0.; }
    let mut result = a(0) * b(0);
    for k in 1..n { result = a(k) * b(k) + result; }
    result
}

/// `general_matrix_vector_product<..., ColMajor, ...>::run`: `res += alpha * lhs * rhs`. Rows covered by packets of
/// two accumulate with fused multiply-adds, a trailing odd row with plain ones; the columns are taken in blocks.
fn gemv_col_major(rows: usize, cols: usize, lhs_stride: usize, lhs: impl Fn(usize, usize) -> f64, rhs: impl Fn(usize) -> f64, res: &mut [f64], alpha: f64) {
    const RES_PACKET_SIZE: usize = 2;
    let block_cols = if cols < 128 { cols } else if lhs_stride * 8 < 32000 { 16 } else { 4 };
    let packet_rows = rows - rows % RES_PACKET_SIZE;
    let mut j2 = 0;
    while j2 < cols {
        let jend = (j2 + block_cols).min(cols);
        for i in 0..rows {
            if i < packet_rows {
                let mut c = 0.;
                for j in j2..jend { c = lhs(i, j).mul_add(rhs(j), c); }
                res[i] = c.mul_add(alpha, res[i]);
            } else {
                let mut c = 0.;
                for j in j2..jend { c += lhs(i, j) * rhs(j); }
                res[i] += alpha * c;
            }
        }
        j2 += block_cols;
    }
}

/// `triangular_solve_vector<..., OnTheLeft, Mode, false, ColMajor>::run` on the leading `size` x `size` triangle of
/// `lhs`: `lower` selects Lower or Upper, `unit_diag` the UnitDiag flag.
fn trsv_col_major(size: usize, lhs: impl Fn(usize, usize) -> f64 + Copy, rhs: &mut [f64], lower: bool, unit_diag: bool) {
    let pw = TRIANGULAR_PANEL_WIDTH;
    if lower {
        let mut pi = 0;
        while pi < size {
            let actual_panel_width = (size - pi).min(pw);
            let start_block = pi;
            let end_block = pi + actual_panel_width;
            for k in 0..actual_panel_width {
                let i = pi + k;
                if rhs[i] != 0. {
                    if !unit_diag { rhs[i] /= lhs(i, i); }
                    let r = actual_panel_width - k - 1; // remaining size
                    let s = i + 1;
                    let ri = rhs[i];
                    for t in 0..r { rhs[s + t] -= ri * lhs(s + t, i); }
                }
            }
            let r = size - end_block; // remaining size
            if r > 0 {
                let (head, tail) = rhs.split_at_mut(end_block);
                gemv_col_major(r, actual_panel_width, 0, |a, b| lhs(end_block + a, start_block + b), |b| head[start_block + b], tail, -1.);
            }
            pi += pw;
        }
    } else {
        let mut pi = size;
        while pi > 0 {
            let actual_panel_width = pi.min(pw);
            let start_block = pi - actual_panel_width;
            for k in 0..actual_panel_width {
                let i = pi - k - 1;
                if rhs[i] != 0. {
                    if !unit_diag { rhs[i] /= lhs(i, i); }
                    let r = actual_panel_width - k - 1; // remaining size
                    let s = i - r;
                    let ri = rhs[i];
                    for t in 0..r { rhs[s + t] -= ri * lhs(s + t, i); }
                }
            }
            let r = start_block; // remaining size
            if r > 0 {
                let (head, tail) = rhs.split_at_mut(start_block);
                gemv_col_major(r, actual_panel_width, 0, |a, b| lhs(a, start_block + b), |b| tail[b], head, -1.);
            }
            pi = pi.saturating_sub(pw);
        }
    }
}

/// `evaluateProductBlockingSizesHeuristic<double, double, KcFactor = 4>` with one thread: the (kc, mc, nc) block
/// sizes of a k-deep product of m x k by k x n.
fn product_blocking_sizes(mut k: i64, mut m: i64, mut n: i64) -> (i64, i64, i64) {
    if k.max(m.max(n)) < 48 { return (k, m, n); }
    const K_PEELING: i64 = 8;
    const KC_FACTOR: i64 = 4;
    let k_div = KC_FACTOR * (GEBP_MR * 8 + GEBP_NR * 8);
    let k_sub = GEBP_MR * (GEBP_NR * 8);
    let max_kc = (((L1 - k_sub) / k_div) & !(K_PEELING - 1)).max(1);
    let old_k = k;
    if k > max_kc {
        k = if k % max_kc == 0 { max_kc } else { max_kc - K_PEELING * ((max_kc - 1 - (k % max_kc)) / (K_PEELING * (k / max_kc + 1))) };
    }
    let actual_l2 = 1572864; // == 1.5 MB
    let lhs_bytes = m * k * 8;
    let remaining_l1 = L1 - k_sub - lhs_bytes;
    let max_nc = if remaining_l1 >= GEBP_NR * 8 * k { remaining_l1 / (k * 8) } else { (3 * actual_l2) / (2 * 2 * max_kc * 8) };
    let nc = (actual_l2 / (2 * k * 8)).min(max_nc) & !(GEBP_NR - 1);
    if n > nc {
        n = if n % nc == 0 { nc } else { nc - GEBP_NR * ((nc - (n % nc)) / (GEBP_NR * (n / nc + 1))) };
    } else if old_k == k {
        let problem_size = k * n * 8;
        let mut actual_lm = actual_l2;
        let mut max_mc = m;
        if problem_size <= 1024 {
            actual_lm = L1;
        } else if L3 != 0 && problem_size <= 32768 {
            actual_lm = L2;
            max_mc = 576.min(max_mc);
        }
        let mut mc = (actual_lm / (3 * k * 8)).min(max_mc);
        if mc > GEBP_MR {
            mc -= mc % GEBP_MR;
        } else if mc == 0 {
            return (k, m, n);
        }
        m = if m % mc == 0 { mc } else { mc - GEBP_MR * ((mc - (m % mc)) / (GEBP_MR * (m / mc + 1))) };
    }
    (k, m, n)
}

/// The per element arithmetic of `gebp_kernel` (mr = 6, nr = 8, double, NEON): `C += alpha * A B` for a `rows` x
/// `depth` block `A` and a `depth` x `cols` block `B` (`a[i + k * rows]`, `b[k + j * depth]`, the packed blocks),
/// `C` the `rows` x `cols` block of `c` at (`ci`, `cj`) with column stride `cs`. Rows covered by packets of two (the
/// 3, 2 and 1 packet paths) accumulate the depth with fused multiply-adds, in order; the 1-packet path on a 4-column
/// block splits the peeled depth between two accumulators (even and odd k) and adds them before the remaining
/// depth. A trailing odd row accumulates with plain multiply-adds. Every result is `R + C * alpha` (fused in the
/// packet paths).
fn gebp(rows: usize, depth: usize, cols: usize, a: &[f64], b: &[f64], c: &mut [f64], ci: usize, cj: usize, cs: usize, alpha: f64) {
    const LHS_PROGRESS: usize = 2;
    let packet_cols4 = (cols / 4) * 4;
    let packet_cols8 = (cols / 8) * 8;
    let peeled_mc3 = (rows / (3 * LHS_PROGRESS)) * (3 * LHS_PROGRESS);
    let peeled_mc2 = peeled_mc3 + ((rows - peeled_mc3) / (2 * LHS_PROGRESS)) * (2 * LHS_PROGRESS);
    let peeled_mc1 = peeled_mc2 + ((rows - peeled_mc2) / LHS_PROGRESS) * LHS_PROGRESS;
    const PK: usize = 8;
    let peeled_kc = depth & !(PK - 1);
    for i in 0..rows {
        for j in 0..cols {
            let at = |k: usize| a[i + k * rows];
            let bt = |k: usize| b[k + j * depth];
            let r = &mut c[(ci + i) + (cj + j) * cs];
            if i < peeled_mc1 {
                let split = i >= peeled_mc2 && j >= packet_cols8 && j < packet_cols4;
                let mut acc = 0.;
                if split {
                    let mut d = 0.;
                    for k in 0..peeled_kc {
                        if k % 2 == 0 { acc = at(k).mul_add(bt(k), acc); } else { d = at(k).mul_add(bt(k), d); }
                    }
                    acc = acc + d;
                    for k in peeled_kc..depth { acc = at(k).mul_add(bt(k), acc); }
                } else {
                    for k in 0..depth { acc = at(k).mul_add(bt(k), acc); }
                }
                *r = acc.mul_add(alpha, *r);
            } else {
                let mut acc = 0.;
                for k in 0..depth { acc = at(k) * bt(k) + acc; }
                *r += alpha * acc;
            }
        }
    }
}

/// `triangular_solve_matrix<double, Index, OnTheLeft, Upper, false, ColMajor, ColMajor, 1>::run`: solves
/// `T X = B` in place, `T` the leading `size` x `size` upper triangle stored at `tri` (stride `tri_stride`), `B` the
/// `size` x `cols` block stored at `other` (stride `other_stride`).
fn trsm_left_upper(size: usize, cols: usize, tri: &[f64], tri_stride: usize, other: &mut [f64], other_stride: usize) {
    if size == 0 || cols == 0 { return; } // Nothing to solve.
    let t = |i: usize, j: usize| tri[i + j * tri_stride];
    let os = other_stride;

    // gemm_blocking_space(rows = size, cols, depth = size, 1 thread, no l3 blocking)
    let (kc, mc, _) = product_blocking_sizes(size as i64, size as i64, cols as i64);
    let kc = kc as usize;
    let mc = (size as i64).min(mc) as usize;

    let small_panel_width = GEBP_MR.max(GEBP_NR) as usize;

    // the goal here is to subdivise the Rhs panels such that we keep some cache coherence when accessing the rhs
    // elements
    let mut subcols = if cols > 0 { L2 / (4 * 8 * (os.max(size) as i64)) } else { 0 };
    subcols = ((subcols / GEBP_NR) * GEBP_NR).max(GEBP_NR);
    let subcols = subcols as usize;

    // pack a rows x depth block of the triangle (gemm_pack_lhs) and a depth x cols block of other (gemm_pack_rhs)
    let pack_a = |i0: usize, k0: usize, rows: usize, depth: usize| -> Vec<f64> {
        let mut v = vec![0.; rows * depth];
        for k in 0..depth { for i in 0..rows { v[i + k * rows] = t(i0 + i, k0 + k); } }
        v
    };
    let pack_b = |o: &[f64], k0: usize, j0: usize, depth: usize, cols: usize| -> Vec<f64> {
        let mut v = vec![0.; depth * cols];
        for j in 0..cols { for k in 0..depth { v[k + j * depth] = o[(k0 + k) + (j0 + j) * os]; } }
        v
    };

    let mut k2 = size;
    while k2 > 0 {
        let actual_kc = k2.min(kc);

        // We have selected and packed a big horizontal panel R1 of rhs. Let B be the packed copy of this panel, and R2
        // the remaining part of rhs. The corresponding vertical panel of lhs is split into A11 (the triangular part)
        // and A21 the remaining rectangular part. Then the high level algorithm is:
        //  - B = R1                    => general block copy (done during the next step)
        //  - R1 = A11^-1 B             => tricky part
        //  - update B from the new R1  => actually this has to be performed continuously during the above step
        //  - R2 -= A21 * B             => GEPP

        // The tricky part: compute R1 = A11^-1 B while updating B from R1. The idea is to split A11 into multiple
        // small vertical panels. Each panel can be split into a small triangular part T1k which is processed without
        // optimization, and the remaining small part T2k which is processed using gebp with appropriate block
        // strides
        let mut j2 = 0;
        while j2 < cols {
            let actual_cols = (cols - j2).min(subcols);
            // for each small vertical panels [T1k^T, T2k^T]^T of lhs
            let mut k1 = 0;
            while k1 < actual_kc {
                let actual_panel_width = (actual_kc - k1).min(small_panel_width);
                // tr solve (trsmKernelL, Upper, col-major triangle)
                {
                    let top = k2 - k1; // one past the last row of the panel
                    for k in 0..actual_panel_width {
                        let i = top - k - 1;
                        let rs = actual_panel_width - k - 1; // remaining size
                        let s = i - rs;
                        let a = 1. / t(i, i);
                        for j in j2..j2 + actual_cols {
                            other[i + j * os] *= a;
                            let b = other[i + j * os];
                            for i3 in 0..rs { other[s + i3 + j * os] -= b * t(s + i3, i); }
                        }
                    }
                }

                let length_target = actual_kc - k1 - actual_panel_width;
                let start_block = k2 - k1 - actual_panel_width;

                // GEBP
                if length_target > 0 {
                    let start_target = k2 - actual_kc;
                    let block_a = pack_a(start_target, start_block, length_target, actual_panel_width);
                    let block_b = pack_b(other, start_block, j2, actual_panel_width, actual_cols);
                    gebp(length_target, actual_panel_width, actual_cols, &block_a, &block_b, other, start_target, j2, os, -1.);
                }
                k1 += small_panel_width;
            }
            j2 += subcols;
        }

        // R2 -= A21 * B => GEPP
        {
            let end = k2 as i64 - kc as i64;
            let mut i2 = 0i64;
            if end > 0 {
                let kb = k2 - kc;
                let block_b = pack_b(other, kb, 0, actual_kc, cols);
                while i2 < end {
                    let actual_mc = (mc as i64).min(end - i2) as usize;
                    if actual_mc > 0 {
                        let block_a = pack_a(i2 as usize, kb, actual_mc, actual_kc);
                        gebp(actual_mc, actual_kc, cols, &block_a, &block_b, other, i2 as usize, 0, os, -1.);
                    }
                    i2 += mc as i64;
                }
            }
        }
        k2 = k2.saturating_sub(kc);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Cases of the comparison against Eigen 5.0.1 (clang++ -std=c++20 -O2 -ffp-contract=off, Apple arm64): the
    // permutations, the LU matrix, the kernel and the solution are bit for bit Eigen's.
    const LU0_DIM: (usize, usize) = (4, 4);
    const LU0_A: [u64; 16] = [0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xbff0000000000000,
        0xbff0000000000000, 0x4000000000000000, 0xbff0000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000];
    const LU0_RHS: [u64; 4] = [0x4000000000000000, 0x3ff0000000000000, 0x0000000000000000, 0x4000000000000000];
    const LU0_P: &str = "P 2 1 0 3";
    const LU0_Q: &str = "Q 2 1 0 3";
    const LU0_LU: &str = "LU 4000000000000000 bfe0000000000000 bfe0000000000000 bfe0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000";
    const LU0_N: &str = "N 1 4000000000000000 1";
    const LU0_K: &str = "K 4 3 0000000000000000 3ff0000000000000 8000000000000000 0000000000000000 3ff0000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 3ff0000000000000";
    const LU0_X: &str = "X 0000000000000000 0000000000000000 0000000000000000 0000000000000000";
    const LU1_DIM: (usize, usize) = (8, 17);
    const LU1_A: [u64; 136] = [0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x8000000000000000, 0xbff0000000000000, 0xbff0000000000000, 0x8000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x3ff0000000000000, 0xbff0000000000000, 0xbff0000000000000,
        0x3ff0000000000000, 0x3ff0000000000000, 0xbff0000000000000, 0x3ff0000000000000, 0x3ff0000000000000,
        0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x3ff0000000000000, 0x0000000000000000,
        0xbff0000000000000, 0x0000000000000000, 0x0000000000000000, 0xbff0000000000000, 0x3ff0000000000000,
        0x3ff0000000000000, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0xbff0000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x3ff0000000000000, 0x8000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x8000000000000000, 0x0000000000000000, 0x0000000000000000, 0xbff0000000000000,
        0x0000000000000000, 0x3ff0000000000000, 0x0000000000000000, 0x0000000000000000, 0x3ff0000000000000,
        0xbff0000000000000, 0xbff0000000000000];
    const LU1_RHS: [u64; 8] = [0x4008000000000000, 0xc008000000000000, 0x3ff0000000000000, 0xc000000000000000,
        0xc000000000000000, 0xc008000000000000, 0xc000000000000000, 0x4000000000000000];
    const LU1_P: &str = "P 3 1 2 0 4 5 6 7";
    const LU1_Q: &str = "Q 1 7 9 3 4 5 6 0 8 2 10 11 12 13 14 15 16";
    const LU1_LU: &str = "LU bff0000000000000 8000000000000000 0000000000000000 8000000000000000 3ff0000000000000 0000000000000000 8000000000000000 8000000000000000 3ff0000000000000 bff0000000000000 3ff0000000000000 bff0000000000000 8000000000000000 3ff0000000000000 bff0000000000000 bff0000000000000 0000000000000000 0000000000000000 bff0000000000000 bff0000000000000 8000000000000000 3ff0000000000000 bff0000000000000 bff0000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 bff0000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 bff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000";
    const LU1_N: &str = "N 3 3ff0000000000000 3";
    const LU1_K: &str = "K 17 14 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 bff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 bff0000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 8000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 bff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000";
    const LU1_X: &str = "X 0000000000000000 4014000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 4008000000000000 0000000000000000 c010000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000";
    const LU2_DIM: (usize, usize) = (7, 8);
    const LU2_A: [u64; 56] = [0x3fdf11db55000000, 0x3fe57dfc18000000, 0x3fe27590ec000000, 0xbfb5fde24c000000,
        0xbfb4075a70000000, 0xbfc3134568000000, 0xbfe13c8357000000, 0xbfeddc3ffd800000, 0xbfdff42d6a000000,
        0x3fecd2a518800000, 0xbfec145e1e800000, 0x3faa8b1b68000000, 0x3feee0e2fa800000, 0xbfe9817500800000,
        0xbfedefab0d800000, 0x3fd35c9c2c000000, 0xbfcea44256000000, 0x3fd4e9317d000000, 0xbfeb6f9746000000,
        0xbfdb5c2424000000, 0x3fd37d5ee6000000, 0x3fed353a9b800000, 0xbfe2e91e87000000, 0x3fdf7f0028000000,
        0x3fd8985458000000, 0x3fbf12db24000000, 0x3fc082b310000000, 0xbfdeecd2a9000000, 0x3fd30de502000000,
        0x3fe997adc0800000, 0x3fcff5c9f2000000, 0xbfc431fad8000000, 0xbfeff5523b000000, 0xbfdbc5a0e2000000,
        0x3fd9725723000000, 0x3fe03e247e000000, 0xbfc0e570ee000000, 0x3fb73951b8000000, 0xbfee9e8630800000,
        0xbfdfd3e28b000000, 0xbfd98a7e13000000, 0xbfd5f60b98000000, 0xbfd367e5b7000000, 0xbfda40df66000000,
        0xbfeb4bc779800000, 0xbfe9beca85000000, 0xbfda2aa862000000, 0x3fdaa9e4b2000000, 0x3fa6a8e758000000,
        0x3fe996de30000000, 0xbfe6421d51000000, 0x3fe2f4eef7000000, 0x3fe250f2a7000000, 0xbfedf4c28b000000,
        0xbfe34d0f3b000000, 0x3fe8cd3702800000];
    const LU2_RHS: [u64; 7] = [0xbfd953c2c1000000, 0x3fe6cb359b800000, 0xbf7ec13dc0000000, 0x3fd9e373f6000000,
        0x3fe2d6fc2f800000, 0xbfd24224a7000000, 0x3fe0624e7a800000];
    const LU2_P: &str = "P 2 1 4 3 0 5 6";
    const LU2_Q: &str = "Q 4 7 2 1 6 5 3 0";
    const LU2_LU: &str = "LU bfeff5523b000000 bfe9a03afaefb2a1 bfd31442f1ea4130 3fc438ba616b4b9b bfd0003bef804fe1 3fdbcee88c76e42b bfd97ad7e3727a79 bfedf4c28b000000 bff71faa286c178a bfd70dc946e38ea9 bfdfe56f23706f05 bfcfbd067dd53dce 3fc16524ddadb06e bfd1d15b97ef7570 bfeb6f9746000000 bfd894813ceafe96 bff545683e9c633e bfca11f6d632dbbb 3fda6d85b9beb5dd 3f61174242da534b 3fbbad9e8b7bc441 3faa8b1b68000000 bfdd4bfaae2e1748 bff15239cdebd27e bff55952286eb10e bfede938b5de9a18 bfe825d8487fda60 3fe2dea20700ce86 bfda2aa862000000 bfe79a9c3b7e3ae4 bfe61b11fca37c67 bff3f8e25eabbb71 c0002836fd9d3b1c 3fbf2ef38df7f8e2 bfcedb984bf80328 bfdfd3e28b000000 bfe0f7cfc4f7aefe 3fc58bc2f9fd719b bff1bb76ac1bb22f bff4540ec56949ee bfe95efba20cbb51 3fdd09b6afb9d39f 3fbf12db24000000 bfdf99a85d94af1b 3fe8aca4fd4d6e0d 3fd1ab803a0441a2 3fd5bdcee354f25a 3fd3c1d179cb6906 bfec0dcaad81492a bfb4075a70000000 3fe37cbb0d3d9378 3fe5cee2bf78551a 3fd799c1b9a03d34 3fe8b09ea47dc3dd bf8e8ff378201f58 bfdfe9104ed22868";
    const LU2_N: &str = "N 7 40002836fd9d3b1c 7";
    const LU2_K: &str = "K 8 1 3ff0000000000000 bfaa23420187fa73 bfa11e0de137636d bfe2331270f82817 bfe4f406114aa36b bfcec10d760f4534 3fdc003acff06f86 3fe02db73284f9d3";
    const LU2_X: &str = "X 0000000000000000 bfc177f052bd804e bfd350c1acc35e32 bfef4c8f97d0d001 bfc0b8fafccd224e bfca1a0164da5155 bfe4a5b4a87d194e 3fa5d36f9deab55e";
    const LU3_DIM: (usize, usize) = (8, 6);
    const LU3_A: [u64; 48] = [0x0000000000000000, 0x0000000000000000, 0xbff0000000000000, 0xbff0000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x3ff0000000000000,
        0xbff0000000000000, 0x0000000000000000, 0x3ff0000000000000, 0xbff0000000000000, 0x0000000000000000,
        0xbff0000000000000, 0xbff0000000000000, 0xbff0000000000000, 0x0000000000000000, 0x0000000000000000,
        0xbff0000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0xbff0000000000000, 0xbff0000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
        0x0000000000000000, 0xbff0000000000000, 0x3ff0000000000000, 0x3ff0000000000000, 0x0000000000000000,
        0x0000000000000000, 0xbff0000000000000, 0x0000000000000000, 0x0000000000000000];
    const LU3_RHS: [u64; 8] = [0x4000000000000000, 0xc008000000000000, 0x4008000000000000, 0xbff0000000000000,
        0xbff0000000000000, 0xc000000000000000, 0x3ff0000000000000, 0xc008000000000000];
    const LU3_P: &str = "P 2 1 0 4 3 5 6 7";
    const LU3_Q: &str = "Q 0 1 2 5 4 3";
    const LU3_LU: &str = "LU bff0000000000000 8000000000000000 8000000000000000 8000000000000000 3ff0000000000000 8000000000000000 8000000000000000 8000000000000000 0000000000000000 bff0000000000000 bff0000000000000 3ff0000000000000 bff0000000000000 8000000000000000 3ff0000000000000 3ff0000000000000 0000000000000000 0000000000000000 bff0000000000000 8000000000000000 3ff0000000000000 8000000000000000 8000000000000000 8000000000000000 3ff0000000000000 3ff0000000000000 0000000000000000 bff0000000000000 8000000000000000 3ff0000000000000 3ff0000000000000 3ff0000000000000 bff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000";
    const LU3_N: &str = "N 4 3ff0000000000000 4";
    const LU3_K: &str = "K 6 2 bff0000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 3ff0000000000000 0000000000000000 0000000000000000";
    const LU3_X: &str = "X c014000000000000 3ff0000000000000 3ff0000000000000 0000000000000000 0000000000000000 c000000000000000";

    fn hx(x: f64) -> String { format!("{:016x}", x.to_bits()) }

    fn run(dim: (usize, usize), a: &[u64], rhs: &[u64]) -> [String; 6] {
        let (rows, cols) = dim;
        let a = Mat::<f64>::from_fn(rows, cols, |i, j| f64::from_bits(a[i + j * rows]));
        let lu = FullPivLu::new(a.as_ref());
        let m = lu.matrix_lu();
        let k = if lu.dimension_of_kernel() > 0 {
            let k = lu.kernel();
            let v: Vec<String> = (0..k.ncols()).flat_map(|j| (0..k.nrows()).map(move |i| (i, j))).map(|(i, j)| hx(k[(i, j)])).collect();
            format!("K {} {} {}", k.nrows(), k.ncols(), v.join(" "))
        } else {
            "K none".to_string()
        };
        let x = lu.solve(&Col::<f64>::from_fn(rows, |i| f64::from_bits(rhs[i])));
        [
            format!("P {}", lu.permutation_p().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ")),
            format!("Q {}", lu.permutation_q().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ")),
            format!("LU {}", (0..cols).flat_map(|j| (0..rows).map(move |i| (i, j))).map(|(i, j)| hx(m[(i, j)])).collect::<Vec<_>>().join(" ")),
            format!("N {} {} {}", lu.nonzero_pivots(), hx(lu.max_pivot()), lu.rank()),
            k,
            format!("X {}", (0..cols).map(|i| hx(x[i])).collect::<Vec<_>>().join(" ")),
        ]
    }

    #[test]
    fn matches_eigen_0() {
        // case 401: 4x4
        let r = run(LU0_DIM, &LU0_A, &LU0_RHS);
        for (got, want) in r.iter().zip([LU0_P, LU0_Q, LU0_LU, LU0_N, LU0_K, LU0_X]) { assert_eq!(got, want); }
    }

    #[test]
    fn matches_eigen_1() {
        // case 550: 8x17
        let r = run(LU1_DIM, &LU1_A, &LU1_RHS);
        for (got, want) in r.iter().zip([LU1_P, LU1_Q, LU1_LU, LU1_N, LU1_K, LU1_X]) { assert_eq!(got, want); }
    }

    #[test]
    fn matches_eigen_2() {
        // case 741: 7x8
        let r = run(LU2_DIM, &LU2_A, &LU2_RHS);
        for (got, want) in r.iter().zip([LU2_P, LU2_Q, LU2_LU, LU2_N, LU2_K, LU2_X]) { assert_eq!(got, want); }
    }

    #[test]
    fn matches_eigen_3() {
        // case 654: 8x6
        let r = run(LU3_DIM, &LU3_A, &LU3_RHS);
        for (got, want) in r.iter().zip([LU3_P, LU3_Q, LU3_LU, LU3_N, LU3_K, LU3_X]) { assert_eq!(got, want); }
    }
}
