//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/polar_svd.cpp (MPL-2.0), for the fixed 2x2 instantiation metriko uses, with
//--- the Eigen 5.0.1 code it runs: Eigen/src/SVD/JacobiSVD.h (2x2 real, ComputeFullU | ComputeFullV),
//--- Eigen/src/misc/RealSvd2x2.h and Eigen/src/Jacobi/Jacobi.h (MPL-2.0).

/// A 2x2 matrix, `m[row][col]`.
pub type Mat2 = [[f64; 2]; 2];

pub const MAT2_ZERO: Mat2 = [[0., 0.], [0., 0.]];
pub const MAT2_IDENTITY: Mat2 = [[1., 0.], [0., 1.]];

/// Eigen's fixed size 2x2 product `a * b` as the reference build evaluates it: the lazy coefficient based product,
/// vectorized along the columns of the result with `pmadd` (fused on NEON): `fma(a(i,1), b(1,j), a(i,0) * b(0,j))`.
pub fn mat2_mul(a: &Mat2, b: &Mat2) -> Mat2 {
    let mut r = MAT2_ZERO;
    for i in 0..2 {
        for j in 0..2 {
            r[i][j] = a[i][1].mul_add(b[1][j], a[i][0] * b[0][j]);
        }
    }
    r
}

pub fn mat2_transpose(a: &Mat2) -> Mat2 { [[a[0][0], a[1][0]], [a[0][1], a[1][1]]] }

/// Eigen's `Matrix2d::determinant()`: `m(0,0) * m(1,1) - m(1,0) * m(0,1)`
pub fn mat2_det(a: &Mat2) -> f64 { a[0][0] * a[1][1] - a[1][0] * a[0][1] }

/// Eigen's `JacobiRotation<double>` (real case)
#[derive(Clone, Copy, Debug)]
struct JacobiRotation {
    c: f64,
    s: f64,
}

impl JacobiRotation {
    /// `operator*`: concatenates two planar rotations
    fn mul(&self, other: &JacobiRotation) -> JacobiRotation {
        JacobiRotation {
            c: self.c * other.c - self.s * other.s,
            s: self.c * other.s + self.s * other.c,
        }
    }

    /// `transpose()`
    fn transpose(&self) -> JacobiRotation { JacobiRotation { c: self.c, s: -self.s } }

    /// `makeJacobi(x, y, z)`
    fn make_jacobi(x: f64, y: f64, z: f64) -> JacobiRotation {
        let deno = 2. * y.abs();
        if deno < f64::MIN_POSITIVE {
            JacobiRotation { c: 1., s: 0. }
        } else {
            let tau = (x - z) / deno;
            let w = (tau * tau + 1.).sqrt();
            let t = if tau > 0. { 1. / (tau + w) } else { 1. / (tau - w) };
            let sign_t = if t > 0. { 1. } else { -1. };
            let n = 1. / (t * t + 1.).sqrt();
            JacobiRotation { c: n, s: -sign_t * (y / y.abs()) * t.abs() * n }
        }
    }
}

/// `internal::apply_rotation_in_the_plane` on two vectors of size 2 (the scalar loop; the vectorized paths compute
/// the same values: `c*x + s*y` and `c*y - s*x == -s*x + c*y`)
fn apply_rotation_in_the_plane(x: &mut [f64; 2], y: &mut [f64; 2], j: JacobiRotation) {
    let (c, s) = (j.c, j.s);
    if c == 1. && s == 0. { return; }
    for i in 0..2 {
        let xi = x[i];
        let yi = y[i];
        x[i] = c * xi + s * yi;
        y[i] = -s * xi + c * yi;
    }
}

/// `m.applyOnTheLeft(p, q, j)`: rotates rows p and q
fn apply_on_the_left(m: &mut Mat2, p: usize, q: usize, j: JacobiRotation) {
    let mut x = m[p];
    let mut y = m[q];
    apply_rotation_in_the_plane(&mut x, &mut y, j);
    m[p] = x;
    m[q] = y;
}

/// `m.applyOnTheRight(p, q, j)`: rotates columns p and q by `j.transpose()`
fn apply_on_the_right(m: &mut Mat2, p: usize, q: usize, j: JacobiRotation) {
    let mut x = [m[0][p], m[1][p]];
    let mut y = [m[0][q], m[1][q]];
    apply_rotation_in_the_plane(&mut x, &mut y, j.transpose());
    m[0][p] = x[0];
    m[1][p] = x[1];
    m[0][q] = y[0];
    m[1][q] = y[1];
}

/// `internal::real_2x2_jacobi_svd(matrix, p, q, &j_left, &j_right)`
fn real_2x2_jacobi_svd(matrix: &Mat2, p: usize, q: usize) -> (JacobiRotation, JacobiRotation) {
    let mut m: Mat2 = [[matrix[p][p], matrix[p][q]], [matrix[q][p], matrix[q][q]]];
    let t = m[0][0] + m[1][1];
    let d = m[1][0] - m[0][1];
    let rot1 = if d.abs() < f64::MIN_POSITIVE {
        JacobiRotation { c: 1., s: 0. }
    } else {
        // If d!=0, then t/d cannot overflow because the magnitude of the
        // entries forming d are not too small compared to the ones forming t.
        let u = t / d;
        let tmp = (1. + u * u).sqrt();
        JacobiRotation { c: u / tmp, s: 1. / tmp }
    };
    apply_on_the_left(&mut m, 0, 1, rot1);
    let j_right = JacobiRotation::make_jacobi(m[0][0], m[0][1], m[1][1]);
    let j_left = rot1.mul(&j_right.transpose());
    (j_left, j_right)
}

/// The output of `JacobiSVD<Matrix2d>(A, ComputeFullU | ComputeFullV)`: `A = U diag(s) V^T`, `s` descending.
/// `valid` is false for a non finite input (Eigen's `InvalidInput`).
pub struct JacobiSvd2 {
    pub u: Mat2,
    pub s: [f64; 2],
    pub v: Mat2,
    pub valid: bool,
}

/// `JacobiSVD<Matrix2d>::compute_impl` (square case: no QR preconditioner, real: no 2x2 complex preconditioning)
pub fn jacobi_svd2(matrix: &Mat2) -> JacobiSvd2 {
    // currently we stop when we reach precision 2*epsilon as the last bit of precision can require an unreasonable
    // number of iterations, only worsening the precision of U and V as we accumulate more rotations
    let precision = 2. * f64::EPSILON;

    // limit for denormal numbers to be considered zero in order to avoid infinite loops (see bug 286)
    let consider_as_zero = f64::MIN_POSITIVE;

    // Scaling factor to reduce over/under-flows: cwiseAbs().maxCoeff<PropagateNaN>()
    let mut scale = 0f64;
    let mut has_nan = false;
    for j in 0..2 {
        for i in 0..2 {
            let a = matrix[i][j].abs();
            if a.is_nan() { has_nan = true; }
            if a > scale { scale = a; }
        }
    }
    if has_nan || !scale.is_finite() {
        // DEVIATION: Eigen returns with info() == InvalidInput and leaves U, V and the singular values unset (igl's
        // polar_svd then reads them uninitialized); here they are NaN.
        return JacobiSvd2 { u: [[f64::NAN; 2]; 2], s: [f64::NAN; 2], v: [[f64::NAN; 2]; 2], valid: false };
    }
    if scale == 0. { scale = 1.; }

    /*** step 1. The R-SVD step: square, nothing to precondition ***/
    let mut w: Mat2 = [[matrix[0][0] / scale, matrix[0][1] / scale], [matrix[1][0] / scale, matrix[1][1] / scale]];
    let mut u = MAT2_IDENTITY;
    let mut v = MAT2_IDENTITY;

    /*** step 2. The main Jacobi SVD iteration. ***/
    // cwiseAbs().diagonal().maxCoeff()
    let mut max_diag_entry = { let a = w[0][0].abs(); let b = w[1][1].abs(); if b > a { b } else { a } };
    let maxi = |a: f64, b: f64| if a < b { b } else { a };

    let mut finished = false;
    while !finished {
        finished = true;

        // do a sweep: for all index pairs (p,q), perform SVD of the corresponding 2x2 sub-matrix
        let (p, q) = (1usize, 0usize);
        // if this 2x2 sub-matrix is not diagonal already...
        let threshold = maxi(consider_as_zero, precision * max_diag_entry);
        if w[p][q].abs() > threshold || w[q][p].abs() > threshold {
            finished = false;
            // perform SVD decomposition of 2x2 sub-matrix corresponding to indices p,q to make it diagonal
            let (j_left, j_right) = real_2x2_jacobi_svd(&w, p, q);

            // accumulate resulting Jacobi rotations
            apply_on_the_left(&mut w, p, q, j_left);
            apply_on_the_right(&mut u, p, q, j_left.transpose());

            apply_on_the_right(&mut w, p, q, j_right);
            apply_on_the_right(&mut v, p, q, j_right);

            // keep track of the largest diagonal coefficient
            max_diag_entry = maxi(max_diag_entry, maxi(w[p][p].abs(), w[q][q].abs()));
        }
    }

    /*** step 3. The work matrix is now diagonal, so ensure it's positive so its diagonal entries are the singular
     * values ***/
    let mut s = [0f64; 2];
    for i in 0..2 {
        let a = w[i][i];
        s[i] = a.abs();
        if a < 0. { u[0][i] = -u[0][i]; u[1][i] = -u[1][i]; }
    }

    s[0] *= scale;
    s[1] *= scale;

    /*** step 4. Sort singular values in descending order ***/
    for i in 0..2 {
        // tail(2 - i).maxCoeff(&pos): the first maximum
        let mut pos = 0usize;
        let mut mx = s[i];
        for k in i + 1..2 { if s[k] > mx { mx = s[k]; pos = k - i; } }
        if mx == 0. { break; }
        if pos != 0 {
            let pos = pos + i;
            s.swap(i, pos);
            for row in u.iter_mut() { row.swap(pos, i); }
            for row in v.iter_mut() { row.swap(pos, i); }
        }
    }

    JacobiSvd2 { u, s, v, valid: true }
}

/// The outputs of `igl::polar_svd(A, R, T, U, S, V)`: `A = R T` with R a rotation (no reflection), and the SVD.
pub struct PolarSvd2 {
    pub r: Mat2,
    pub t: Mat2,
    pub u: Mat2,
    pub s: [f64; 2],
    pub v: Mat2,
}

/// `igl::polar_svd(A, R, T, U, S, V)` = `polar_svd(A, includeReflections = false, ...)` for a 2x2 matrix
pub fn polar_svd(a: &Mat2) -> PolarSvd2 {
    let svd = jacobi_svd2(a);
    let u = svd.u;
    let v = svd.v;
    let s = svd.s;
    let mut r = mat2_mul(&u, &mat2_transpose(&v));
    // SVT = S.asDiagonal() * V.adjoint()
    let svt: Mat2 = [[s[0] * v[0][0], s[0] * v[1][0]], [s[1] * v[0][1], s[1] * v[1][1]]];
    let t;
    // Check for reflection
    if mat2_det(&r) < 0. {
        let mut w = v;
        w[0][1] *= -1.;
        w[1][1] *= -1.;
        r = mat2_mul(&u, &mat2_transpose(&w));
        t = mat2_mul(&w, &svt);
    } else {
        t = mat2_mul(&v, &svt);
    }
    PolarSvd2 { r, t, u, s, v }
}
