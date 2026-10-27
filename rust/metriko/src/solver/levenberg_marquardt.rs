//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/solver/levenberg_marquardt.h (see rust/PORTED_FROM).

// The C++ templates `DiagonalDamping<SolverTraits>` and `LMSolver<LinearSolver, SolverTraits, DampingTraits>` become
// the traits below. `LMSolver::init` stored pointers to the three collaborators; here `init` takes the sizes from the
// traits object and `solve` receives the collaborators.

use faer::Col;
use faer::sparse::Triplet;
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets, nonzeros};

/// Eigen's `A * x` for a column-major sparse `A` (`sparse_time_dense_product`, ColMajor): every column `j` scatters
/// `a_ij * x_j` into the result, so each row sums its terms in increasing column order, starting from zero.
pub fn spmv(a: &SprsD, x: &Col<f64>) -> Col<f64> {
    let mut res = Col::<f64>::zeros(a.nrows());
    for (i, j, v) in nonzeros(a.as_ref()) { res[i] += v * x[j]; }
    res
}

/// `A * x` with every coefficient of `A` scaled by `s` first: Eigen's `s * A * x` evaluates `(s * a_ij) * x_j`.
pub fn spmv_scaled(a: &SprsD, s: f64, x: &Col<f64>) -> Col<f64> {
    let mut res = Col::<f64>::zeros(a.nrows());
    for (i, j, v) in nonzeros(a.as_ref()) { res[i] += (s * v) * x[j]; }
    res
}

/// Eigen's `A.transpose() * x` for a column-major sparse `A` (a row-major product): entry `j` is the sum over the
/// stored entries of column `j`, in storage order, starting from zero.
pub fn tspmv(a: &SprsD, x: &Col<f64>) -> Col<f64> {
    let mut res = Col::<f64>::zeros(a.ncols());
    for j in 0..a.ncols() {
        let rows = a.row_idx_of_col_raw(j);
        let vals = a.val_of_col(j);
        let mut tmp = 0.;
        for (&i, &v) in rows.iter().zip(vals.iter()) { tmp += v * x[i]; }
        res[j] += tmp;
    }
    res
}

/// `v.squaredNorm()`
pub fn squared_norm(v: &Col<f64>) -> f64 {
    let mut s = 0.;
    for i in 0..v.nrows() { s += v[i] * v[i]; }
    s
}

/// `v.norm()`
pub fn vnorm(v: &Col<f64>) -> f64 { squared_norm(v).sqrt() }

/// `v.lpNorm<Infinity>()`
pub fn inf_norm(v: &Col<f64>) -> f64 {
    let mut m = 0.;
    for i in 0..v.nrows() { let a = v[i].abs(); if a > m { m = a; } }
    m
}

/// `a + b`
pub fn vadd(a: &Col<f64>, b: &Col<f64>) -> Col<f64> { Col::<f64>::from_fn(a.nrows(), |i| a[i] + b[i]) }

/// `a - b`
pub fn vsub(a: &Col<f64>, b: &Col<f64>) -> Col<f64> { Col::<f64>::from_fn(a.nrows(), |i| a[i] - b[i]) }

/// the linear solver of an LM step: factorizes the damped jacobian and solves the normal equations
pub trait LinearSolver {
    fn factorize(&mut self, j: &SprsD) -> bool;
    fn solve(&mut self, rhs: &Col<f64>, x: &mut Col<f64>) -> bool;
}

/// the problem: `xSize`, `ESize`, `initial_solution`, `pre_iteration`, `objective_jacobian`, `post_iteration`
pub trait SolverTraits {
    fn x_size(&self) -> usize;
    fn e_size(&self) -> usize;
    fn initial_solution(&self, x0: &mut Col<f64>);
    fn pre_iteration(&mut self, prevx: &Col<f64>);
    fn objective_jacobian(&mut self, x: &Col<f64>, e: &mut Col<f64>, j: &mut SprsD, update_j: bool);
    fn post_iteration(&mut self, x: &Col<f64>) -> bool;
}

pub trait DampingTraits {
    fn init(&self, j: &SprsD, damp_j: &mut SprsD);
    fn update(&mut self, j: &SprsD, prv_energy2: f64, new_energy2: f64, damp_j: &mut SprsD);
}

/// `DiagonalDamping<SolverTraits>` (the template parameter is unused in the C++)
pub struct DiagonalDamping {
    pub curr_lambda: f64,
}

impl DiagonalDamping {
    /// `DiagonalDamping(double _currLambda = 0.01)`
    pub fn new(curr_lambda: f64) -> Self { DiagonalDamping { curr_lambda } }

    pub fn build_damped_matrix(&self, j: &SprsD, damp_j: &mut SprsD) {
        let mut damp_vec = vec![0.; j.ncols()];
        let mut damp_j_tris: Vec<TripD> = vec![];
        for (r, c, v) in nonzeros(j.as_ref()) {
            damp_vec[c] += self.curr_lambda * v * v;
            damp_j_tris.push(Triplet::new(r, c, v));
        }
        for i in 0..damp_vec.len() {
            damp_j_tris.push(Triplet::new(j.nrows() + i, i, damp_vec[i].sqrt()));
        }
        *damp_j = from_triplets(j.nrows() + j.ncols(), j.ncols(), &damp_j_tris);
    }
}

impl DampingTraits for DiagonalDamping {
    fn init(&self, j: &SprsD, damp_j: &mut SprsD) { self.build_damped_matrix(j, damp_j); }

    fn update(&mut self, j: &SprsD, prv_energy2: f64, new_energy2: f64, damp_j: &mut SprsD) {
        let f = prv_energy2 > new_energy2 && new_energy2 != f64::INFINITY;
        if f { self.curr_lambda /= 10.; } else { self.curr_lambda *= 10.; }
        self.build_damped_matrix(j, damp_j);
    }
}

/// `LMSolver<LinearSolver, SolverTraits, DampingTraits>`
pub struct LMSolver {
    pub x: Col<f64>,              // current solution; always updated
    pub prevx: Col<f64>,          // the solution of the previous iteration
    pub x0: Col<f64>,             // the initial solution to the system
    pub d: Col<f64>,              // the direction taken.
    pub curr_objective: Col<f64>, // the current value of the energy
    pub prev_objective: Col<f64>, // the previous value of the energy

    pub func_tolerance: f64,
    pub foo_tolerance: f64,

    // always updated to the current iteration
    pub energy: f64,
    pub foo_optimality: f64,
    pub max_iter: i32,
    pub cur_iter: i32,
}

impl Default for LMSolver {
    fn default() -> Self { Self::new() }
}

/// Eigen's `resize`: a no-op when the size does not change, otherwise new storage.
/// DEVIATION: the new storage of a resized Eigen vector is uninitialized; here it is zero (only observable when the
/// C++ reads `x` before any iteration wrote it, i.e. when the very first factorization fails).
fn resize(v: &mut Col<f64>, n: usize) { if v.nrows() != n { *v = Col::<f64>::zeros(n); } }

impl LMSolver {
    pub fn new() -> Self {
        LMSolver {
            x: Col::zeros(0),
            prevx: Col::zeros(0),
            x0: Col::zeros(0),
            d: Col::zeros(0),
            curr_objective: Col::zeros(0),
            prev_objective: Col::zeros(0),
            func_tolerance: 0.,
            foo_tolerance: 0.,
            energy: 0.,
            foo_optimality: 0.,
            max_iter: 0,
            cur_iter: 0,
        }
    }

    /// `init(LS, ST, DT, _maxIterations = 100, _funcTolerance = 10e-10, _fooTolerance = 10e-10)`
    pub fn init<S: SolverTraits>(&mut self, st: &S, max_iterations: i32, func_tolerance: f64, foo_tolerance: f64) {
        self.max_iter = max_iterations;
        self.func_tolerance = func_tolerance;
        self.foo_tolerance = foo_tolerance;

        resize(&mut self.d, st.x_size());
        resize(&mut self.x, st.x_size());
        resize(&mut self.x0, st.x_size());
        resize(&mut self.prevx, st.x_size());
        resize(&mut self.curr_objective, st.e_size());
        resize(&mut self.curr_objective, st.e_size());
    }

    /// `init(LS, ST, DT)` with the default iteration count and tolerances
    pub fn init_default<S: SolverTraits>(&mut self, st: &S) { self.init(st, 100, 10e-10, 10e-10); }

    pub fn solve<L: LinearSolver, S: SolverTraits, D: DampingTraits>(&mut self, ls: &mut L, st: &mut S, dt: &mut D, verbose: bool) -> bool {
        st.initial_solution(&mut self.x0);
        self.prevx = self.x0.clone();

        let mut rhs: Col<f64>;
        let mut dir = Col::<f64>::zeros(0);
        if verbose { println!("******Beginning Optimization******"); }

        // estimating initial miu
        let mut damp_j = SprsD::try_new_from_triplets(0, 0, &[]).unwrap();
        let mut e_cur = Col::<f64>::zeros(0);
        let mut e_prv: Col<f64>;
        let mut j = SprsD::try_new_from_triplets(0, 0, &[]).unwrap();

        self.cur_iter = 0;
        st.objective_jacobian(&self.prevx, &mut e_cur, &mut j, true);
        dt.init(&j, &mut damp_j);

        loop {
            st.pre_iteration(&self.prevx);
            let prv_energy2 = squared_norm(&e_cur);
            rhs = tspmv(&j, &e_cur);
            for i in 0..rhs.nrows() { rhs[i] = -rhs[i]; }

            self.foo_optimality = inf_norm(&rhs);
            if self.foo_optimality < self.foo_tolerance { self.x = self.prevx.clone(); break; }

            // solving to get the LM direction
            if !ls.factorize(&damp_j) { println!("Solver Failed to factorize! "); return false; }

            ls.solve(&rhs, &mut dir);

            if vnorm(&dir) < self.func_tolerance { self.x = self.prevx.clone(); return true; }

            e_prv = e_cur.clone();

            st.objective_jacobian(&vadd(&self.prevx, &dir), &mut e_cur, &mut j, false);
            let new_energy2 = squared_norm(&e_cur);
            self.energy = new_energy2;
            let accepted = prv_energy2 > new_energy2;
            self.x = if accepted { vadd(&self.prevx, &dir) } else { self.prevx.clone() };

            if (prv_energy2 - new_energy2).abs() < self.func_tolerance { break; }
            if accepted { st.objective_jacobian(&self.x, &mut e_cur, &mut j, true); } else { e_cur = e_prv; }

            self.energy = squared_norm(&e_cur);

            dt.update(&j, prv_energy2, new_energy2, &mut damp_j);

            if st.post_iteration(&self.x) { return true; }

            self.cur_iter += 1;
            self.prevx = self.x.clone();
            if self.cur_iter > self.max_iter { break; } // do { ... } while (curIter <= maxIter)
        }

        false
    }
}
