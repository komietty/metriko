//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of Eigen 5.0.1 Eigen/src/Geometry/Umeyama.h (MPL-2.0) for a 3 x 4 point set, with the parts of
//--- Eigen/src/SVD/JacobiSVD.h (square 3 x 3, full U and V), Eigen/src/misc/RealSvd2x2.h, Eigen/src/Jacobi/Jacobi.h
//--- and Eigen/src/LU/Determinant.h it runs.

/// a 3 x 3 matrix, `m[row][col]`
pub type Mat3 = [[f64; 3]; 3];
/// a 3 x 4 matrix, `m[row][col]`
pub type Mat34 = [[f64; 4]; 3];
/// a 4 x 4 matrix, `m[row][col]`
pub type Mat4 = [[f64; 4]; 4];

/// `JacobiRotation<double>` (c, s)
#[derive(Clone, Copy, Debug)]
struct Jrot {
    c: f64,
    s: f64,
}

impl Jrot {
    fn transpose(self) -> Jrot { Jrot { c: self.c, s: -self.s } }

    /// `operator*`: (c1 c2 - s1 s2, c1 s2 + s1 c2)
    fn mul(self, o: Jrot) -> Jrot { Jrot { c: self.c * o.c - self.s * o.s, s: self.c * o.s + self.s * o.c } }

    /// `makeJacobi(x, y, z)`
    fn make_jacobi(x: f64, y: f64, z: f64) -> Jrot {
        let deno = 2. * y.abs();
        if deno < f64::MIN_POSITIVE {
            Jrot { c: 1., s: 0. }
        } else {
            let tau = (x - z) / deno;
            let w = (tau * tau + 1.).sqrt();
            let t = if tau > 0. { 1. / (tau + w) } else { 1. / (tau - w) };
            let sign_t = if t > 0. { 1. } else { -1. };
            let n = 1. / (t * t + 1.).sqrt();
            Jrot { c: n, s: -sign_t * (y / y.abs()) * t.abs() * n }
        }
    }
}

/// `internal::apply_rotation_in_the_plane(x, y, j)` on two 3-vectors given by element accessors
#[inline]
fn rotate(xs: [&mut f64; 3], ys: [&mut f64; 3], j: Jrot) {
    let (c, s) = (j.c, j.s);
    if c == 1. && s == 0. { return; }
    for (x, y) in xs.into_iter().zip(ys) {
        let xi = *x;
        let yi = *y;
        *x = c * xi + s * yi;
        *y = -s * xi + c * yi;
    }
}

/// `m.applyOnTheLeft(p, q, j)`: rows p and q
fn apply_on_the_left(m: &mut Mat3, p: usize, q: usize, j: Jrot) {
    let (rp, rq) = if p < q { let (a, b) = m.split_at_mut(q); (&mut a[p], &mut b[0]) } else { let (a, b) = m.split_at_mut(p); (&mut b[0], &mut a[q]) };
    let [x0, x1, x2] = rp;
    let [y0, y1, y2] = rq;
    rotate([x0, x1, x2], [y0, y1, y2], j);
}

/// `m.applyOnTheRight(p, q, j)`: columns p and q, rotated by `j.transpose()`
fn apply_on_the_right(m: &mut Mat3, p: usize, q: usize, j: Jrot) {
    let j = j.transpose();
    let mut xs = [m[0][p], m[1][p], m[2][p]];
    let mut ys = [m[0][q], m[1][q], m[2][q]];
    {
        let [x0, x1, x2] = &mut xs;
        let [y0, y1, y2] = &mut ys;
        rotate([x0, x1, x2], [y0, y1, y2], j);
    }
    for r in 0..3 { m[r][p] = xs[r]; m[r][q] = ys[r]; }
}

/// `internal::real_2x2_jacobi_svd(matrix, p, q, &j_left, &j_right)`
fn real_2x2_jacobi_svd(matrix: &Mat3, p: usize, q: usize) -> (Jrot, Jrot) {
    let mut m = [[matrix[p][p], matrix[p][q]], [matrix[q][p], matrix[q][q]]];
    let t = m[0][0] + m[1][1];
    let d = m[1][0] - m[0][1];
    let rot1 = if d.abs() < f64::MIN_POSITIVE {
        Jrot { c: 1., s: 0. }
    } else {
        let u = t / d;
        let tmp = (1. + u * u).sqrt();
        Jrot { c: u / tmp, s: 1. / tmp }
    };
    // m.applyOnTheLeft(0, 1, rot1)
    if !(rot1.c == 1. && rot1.s == 0.) {
        for k in 0..2 {
            let xi = m[0][k];
            let yi = m[1][k];
            m[0][k] = rot1.c * xi + rot1.s * yi;
            m[1][k] = -rot1.s * xi + rot1.c * yi;
        }
    }
    let j_right = Jrot::make_jacobi(m[0][0], m[0][1], m[1][1]);
    let j_left = rot1.mul(j_right.transpose());
    (j_left, j_right)
}

/// `JacobiSVD<Matrix3d, ComputeFullU | ComputeFullV>`: `(U, singular values, V)`
pub fn jacobi_svd3(matrix: &Mat3) -> (Mat3, [f64; 3], Mat3) {
    let precision = 2. * f64::EPSILON;
    let consider_as_zero = f64::MIN_POSITIVE;
    let maxi = |a: f64, b: f64| if a < b { b } else { a }; // numext::maxi = std::max

    let mut u: Mat3 = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let mut v: Mat3 = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let mut sv = [0.; 3];

    // scale = matrix.cwiseAbs().maxCoeff<PropagateNaN>()
    let mut scale = matrix[0][0].abs();
    for c in 0..3 { for r in 0..3 { let a = matrix[r][c].abs(); if a > scale || a.is_nan() { scale = a; } } }
    if !scale.is_finite() { return ([[f64::NAN; 3]; 3], [f64::NAN; 3], [[f64::NAN; 3]; 3]); }
    if scale == 0. { scale = 1.; }

    /*** step 1. square: the work matrix is the matrix divided by scale ***/
    let mut w = *matrix;
    for r in 0..3 { for c in 0..3 { w[r][c] /= scale; } }

    /*** step 2. The main Jacobi SVD iteration. ***/
    let mut max_diag_entry = w[0][0].abs();
    for i in 1..3 { if w[i][i].abs() > max_diag_entry { max_diag_entry = w[i][i].abs(); } }

    let mut finished = false;
    while !finished {
        finished = true;
        for p in 1..3 {
            for q in 0..p {
                let threshold = maxi(consider_as_zero, precision * max_diag_entry);
                if w[p][q].abs() > threshold || w[q][p].abs() > threshold {
                    finished = false;
                    let (j_left, j_right) = real_2x2_jacobi_svd(&w, p, q);
                    apply_on_the_left(&mut w, p, q, j_left);
                    apply_on_the_right(&mut u, p, q, j_left.transpose());
                    apply_on_the_right(&mut w, p, q, j_right);
                    apply_on_the_right(&mut v, p, q, j_right);
                    max_diag_entry = maxi(max_diag_entry, maxi(w[p][p].abs(), w[q][q].abs()));
                }
            }
        }
    }

    /*** step 3. The work matrix is now diagonal, so ensure it's positive so its diagonal entries are the singular values ***/
    for i in 0..3 {
        let a = w[i][i];
        sv[i] = a.abs();
        if a < 0. { for r in 0..3 { u[r][i] = -u[r][i]; } }
    }
    for x in sv.iter_mut() { *x *= scale; }

    /*** step 4. Sort singular values in descending order ***/
    for i in 0..3 {
        // m_singularValues.tail(3 - i).maxCoeff(&pos): the first maximum
        let mut pos = 0;
        let mut max_rem = sv[i];
        for k in (i + 1)..3 { if sv[k] > max_rem { max_rem = sv[k]; pos = k - i; } }
        if max_rem == 0. { break; }
        if pos != 0 {
            let pos = pos + i;
            sv.swap(i, pos);
            for r in 0..3 { u[r].swap(i, pos); v[r].swap(i, pos); }
        }
    }
    (u, sv, v)
}

/// `MatrixBase::determinant()` of a 3 x 3 matrix
pub fn determinant3(m: &Mat3) -> f64 {
    let h = |a: usize, b: usize, c: usize| m[0][a] * (m[1][b] * m[2][c] - m[1][c] * m[2][b]);
    h(0, 1, 2) - h(1, 0, 2) + h(2, 0, 1)
}

/// `Eigen::umeyama(src, dst, with_scaling)` for 3 x 4 point sets: the 4 x 4 similarity mapping src onto dst
pub fn umeyama(src: &Mat34, dst: &Mat34, with_scaling: bool) -> Mat4 {
    let m = 3;
    let n = 4;
    let one_over_n = 1. / n as f64;

    // src.rowwise().sum() * one_over_n: the unrolled reduction of a strided row halves, (x0 + x1) + (x2 + x3)
    let row_sum = |a: &Mat34, r: usize| (a[r][0] + a[r][1]) + (a[r][2] + a[r][3]);
    let src_mean = [row_sum(src, 0) * one_over_n, row_sum(src, 1) * one_over_n, row_sum(src, 2) * one_over_n];
    let dst_mean = [row_sum(dst, 0) * one_over_n, row_sum(dst, 1) * one_over_n, row_sum(dst, 2) * one_over_n];

    let mut src_demean = [[0.; 4]; 3];
    let mut dst_demean = [[0.; 4]; 3];
    for r in 0..3 {
        for c in 0..4 {
            src_demean[r][c] = src[r][c] - src_mean[r];
            dst_demean[r][c] = dst[r][c] - dst_mean[r];
        }
    }

    // sigma = one_over_n * dst_demean * src_demean.transpose(): each coefficient is the vectorized reduction of the 4
    // products of two contiguous rows, [t0, t1] + [t2, t3] then the two lanes: (t0 + t2) + (t1 + t3)
    let mut sigma: Mat3 = [[0.; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let t = |k: usize| (one_over_n * dst_demean[i][k]) * src_demean[j][k];
            sigma[i][j] = (t(0) + t(2)) + (t(1) + t(3));
        }
    }

    let (u, sv, v) = jacobi_svd3(&sigma);

    let mut rt: Mat4 = [[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]];
    let mut s = [1.; 3];
    if determinant3(&u) * determinant3(&v) < 0. { s[m - 1] = -1.; }

    // Rt.block(0, 0, m, m) = U * S.asDiagonal() * V.transpose(): a lazy product assigned to a 3 x 3 block by slices.
    // Eigen evaluates rows 0, 1 as one packet with pmul / pmadd, which on arm64 NEON (the reference build) is a fused
    // multiply-add whatever -ffp-contract says, and row 2 by the scalar unroller, which halves: t0 + (t1 + t2).
    // Checked bitwise against Eigen 5.0.1 (Apple clang 17, -O2 -ffp-contract=off) on 2000 random point sets.
    for i in 0..3 {
        let l = |k: usize| u[i][k] * s[k];
        for j in 0..3 {
            rt[i][j] = if i < 2 {
                l(2).mul_add(v[j][2], l(1).mul_add(v[j][1], l(0) * v[j][0]))
            } else {
                l(0) * v[j][0] + (l(1) * v[j][1] + l(2) * v[j][2])
            };
        }
    }

    if with_scaling {
        // src_demean.rowwise().squaredNorm().sum() * one_over_n
        let sq = |r: usize| { let a = &src_demean[r]; ((a[0] * a[0] + a[1] * a[1]) + a[2] * a[2]) + a[3] * a[3] };
        let src_var = ((sq(0) + sq(1)) + sq(2)) * one_over_n;
        let c = 1. / src_var * ((sv[0] * s[0] + sv[1] * s[1]) + sv[2] * s[2]);
        for i in 0..3 {
            let t = |k: usize| (c * rt[i][k]) * src_mean[k];
            rt[i][3] = dst_mean[i] - ((t(0) + t(1)) + t(2));
        }
        for i in 0..3 { for j in 0..3 { rt[i][j] *= c; } }
    } else {
        for i in 0..3 {
            let t = |k: usize| rt[i][k] * src_mean[k];
            rt[i][3] = dst_mean[i] - ((t(0) + t(1)) + t(2));
        }
    }
    rt
}
