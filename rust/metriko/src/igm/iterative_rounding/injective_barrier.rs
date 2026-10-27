//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/igm/iterative_rounding/injective_barrier.h (see rust/PORTED_FROM).

use faer::Col;
use faer::sparse::Triplet;
use crate::common::typedef::Row2d;
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets, diagonal};

pub struct InjectiveBarrier {
    n_big: i32,
    n_f: i32,
    i: Vec<[i32; 4]>,
    j: Vec<[i32; 4]>,
    s: Vec<[f64; 4]>,
    v: Vec<Vec<f64>>, // field volume, nF x N
}

impl InjectiveBarrier {
    /// `InjectiveBarrier(N, nF, F2)`
    pub fn new(n_big: i32, n_f: i32, f2: &Col<f64>) -> Self {
        let nn = (n_big * n_f) as usize;
        // DEVIATION: the C++ leaves S uninitialized until the first `update(.., updateJ = true)`, which always
        // precedes its use; here it starts at zero.
        let mut b = InjectiveBarrier { n_big, n_f, i: vec![[0; 4]; nn], j: vec![[0; 4]; nn], s: vec![[0.; 4]; nn], v: vec![vec![0.; n_big as usize]; n_f as usize] };
        let nb = n_big;
        for i in 0..n_f {
            for j in 0..nb {
                let r = (i * nb + j) as usize;
                b.i[r] = [i * nb + j; 4];
                b.j[r] = [
                    2 * nb * i + 2 * j,
                    2 * nb * i + 2 * j + 1,
                    2 * nb * i + 2 * ((j + 1) % nb),
                    2 * nb * i + 2 * ((j + 1) % nb) + 1,
                ];
            }
        }

        for i in 0..n_f {
            for j in 0..nb {
                let oft = (2 * nb * i) as usize;
                let v0 = Row2d::new(f2[oft + 2 * j as usize], f2[oft + 2 * j as usize + 1]);
                let k = 2 * ((j + 1) % nb) as usize;
                let v1 = Row2d::new(f2[oft + k], f2[oft + k + 1]);
                b.v[i as usize][j as usize] = v0.length() * v1.length();
            }
        }
        b
    }

    pub fn gen_image_filed(&self, curr_field: &Col<f64>) -> SprsD {
        let mut t: Vec<TripD> = vec![];
        for i in 0..self.i.len() {
            for j in 0..4 {
                t.push(Triplet::new(self.i[i][j] as usize, self.j[i][j] as usize, self.s[i][j]));
            }
        }
        from_triplets((self.n_big * self.n_f) as usize, curr_field.nrows(), &t)
    }

    pub fn update(
        &mut self,
        field: &Col<f64>, // todo: does not change??
        update_j: bool,
        f_barrier: &mut Col<f64>,
        g_barrier: &mut SprsD,
    ) {
        let s = 0.5;
        let nb = self.n_big as usize;
        let nn = (self.n_big * self.n_f) as usize;

        *f_barrier = Col::<f64>::zeros(nn);
        let mut bar_spline = Col::<f64>::zeros(nn);
        let mut spline_derivative = Col::<f64>::zeros(0);
        if update_j { spline_derivative = Col::<f64>::zeros(nn); }

        for i in 0..self.n_f as usize {
            for j in 0..nb {
                let a = 2 * nb * i + 2 * j;
                let b = 2 * nb * i + 2 * ((j + 1) % nb);
                let v0 = [field[a], field[a + 1]];
                let v1 = [field[b], field[b + 1]];
                let v = self.v[i][j];
                let r = (v0[0] * v1[1] - v0[1] * v1[0]) / v; // volume ratio
                let bar = (r / s).powf(3.) - 3. * (r / s).powf(2.) + 3. * r / s;
                let mut bar2 = 1. / bar - 1.;
                if r <= 0. { bar2 = f64::INFINITY; }
                if r >= s { bar2 = 0.; }
                f_barrier[nb * i + j] = bar2;
                bar_spline[nb * i + j] = bar;

                if update_j {
                    let mut d = 3. * (r * r / (s * s * s)) - 6. * (r / (s * s)) + 3. / s;
                    if r <= 0. { d = f64::INFINITY; }
                    if r >= s { d = 0.; }
                    spline_derivative[nb * i + j] = d;
                    self.s[nb * i + j] = [
                         v1[1] / v,
                        -v1[0] / v,
                        -v0[1] / v,
                         v0[0] / v,
                    ];
                }
            }
        }

        if update_j {
            let inf = f64::INFINITY;
            let mut bar_der = vec![0.; nn];
            for k in 0..nn {
                let bd = -spline_derivative[k] / (bar_spline[k] * bar_spline[k]);
                bar_der[k] = if f_barrier[k].abs() < 10e-9 { 0. } else if f_barrier[k] == inf { inf } else { bd };
            }
            *g_barrier = diagonal(&bar_der);
        }

        //gBarrier = gBarrier * gen_image_filed(field);
    }
}
