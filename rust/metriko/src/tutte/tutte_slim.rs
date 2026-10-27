//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tutte/tutte_slim.h (see rust/PORTED_FROM).

use faer::Col;
use faer::sparse::Triplet;
use crate::common::typedef::*;
use crate::common::utilities::get_quater_rot;
use crate::common::stdcomplex::std_div;
use crate::hmesh::hmesh::*;
use crate::igl::doublearea::doublearea;
use crate::igl::flip_avoiding_line_search::flip_avoiding_line_search;
use crate::igl::grad::grad;
use crate::igl::local_basis::local_basis;
use crate::igl::polar_svd::{polar_svd, Mat2, mat2_mul, mat2_transpose, mat2_det, MAT2_ZERO};
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets, transpose};
use crate::solver::sparse_solver::{ComputationInfo, SparseLdlt, SparseLlt};

// Arithmetic of the Eigen expressions below, as the reference build evaluates them (checked bitwise against Eigen
// 5.0.1, clang -O2 -ffp-contract=off, arm64): fixed size products with a plain left hand side are lazy coefficient
// based products vectorized with pmadd, i.e. a fused chain `fma(l(i,k), r(k,j), ...)` started by `l(i,0) * r(0,j)`;
// products whose left hand side is a transpose (`A.transpose() * A`) reduce each coefficient as
// `(l0 r0 + l2 r2) + (l1 r1 + l3 r3)` without fusing; `Matrix2d::squaredNorm()` is `(a00^2 + a01^2) + (a10^2 + a11^2)`;
// outer products and the determinant are plain.

/// Eigen's `Matrix2d::squaredNorm()`: two packets (the columns) added, then reduced
fn squared_norm2(j: &Mat2) -> f64 {
    (j[0][0] * j[0][0] + j[0][1] * j[0][1]) + (j[1][0] * j[1][0] + j[1][1] * j[1][1])
}

struct Link { to: usize, r: Complex, t: Complex }

// slim (rabinovich et al. 2017) with the symmetric dirichlet energy on the mesh cut along its seam, kept seamless:
// the copies of a vertex on the seam move together under the transitions between their charts (a rotation by a
// multiple of 90 degrees and an integer translation), so the result stays an integer-grid map while the seam is
// free to move. only the cones are fixed
/// Returns the uv per vertex of `hc`.
pub fn seamless_slim(
    hm: &Hmesh,          // the mesh with its seam
    hc: &Hmesh,          // hm cut along the seam: the same faces, corner by corner
    seam: &[bool],       // per edge of hm
    singular: &[i32],    // per vertex of hm
    uv0: &[[f64; 2]],    // per vertex of hc: a seamless map to start from
    max_iter: i32,
    rel_tol: f64,
) -> Vec<[f64; 2]> {
    let n = hc.n_v as usize;
    let nf = hc.n_f as usize;

    // 1: the transitions across the seam: z_b = r z_a + t from a copy a of a vertex to the copy b on the other side
    let mut links: Vec<Vec<Link>> = (0..n).map(|_| vec![]).collect();
    let mut org = vec![0i32; n];
    for f in 0..nf { for j in 0..3 { org[hc.idx[f][j] as usize] = hm.idx[f][j]; } }

    let copy_at = |f: Face<'_>, v: Vert<'_>| -> usize {
        for j in 0..3 { if hm.idx[f.id as usize][j] == v.id { return hc.idx[f.id as usize][j] as usize; } }
        metriko_fail!("vertex {} is not on face {}", v.id, f.id);
    };
    let uv_at = |k: usize| cx(uv0[k][0], uv0[k][1]);

    let mut gap = 0f64;
    for e in hm.edges() {
        if !seam[e.id as usize] || e.is_boundary() { continue; }
        let f0 = e.face0();
        let f1 = e.face1();
        let p = e.vert0();
        let q = e.vert1();
        let d0 = uv_at(copy_at(f0, q)) - uv_at(copy_at(f0, p));
        let d1 = uv_at(copy_at(f1, q)) - uv_at(copy_at(f1, p));
        let r = get_quater_rot((std_div(d1, d0).arg() / (PI / 2.)).round() as i32);
        let t0 = uv_at(copy_at(f1, p)) - r * uv_at(copy_at(f0, p));
        let t = cx(t0.re.round(), t0.im.round());
        let g = (t0 - t).norm();
        if gap < g { gap = g; }
        for v in [p, q] {
            let a = copy_at(f0, v);
            let b = copy_at(f1, v);
            links[a].push(Link { to: b, r, t });
            links[b].push(Link { to: a, r: r.conj(), t: -r.conj() * t });
        }
    }
    metriko_check!(gap < 1e-6, "the seam transitions are not integer (off by {})", gap);

    // 2: one free point per vertex of hm: its first copy, the others follow it. the cones stay where they are, and so
    //    does any vertex whose transitions do not close around it
    let mut var = vec![-1i32; n];            // the free point a copy follows, -1 when fixed
    let mut rk = vec![cx(1., 0.); n];        // the copy is rk z + tk of its free point z
    let mut tk = vec![cx(0., 0.); n];
    let mut first: Vec<usize> = vec![];      // the first copy of each free point
    let mut copies: Vec<usize> = vec![];     // the number of copies of each free point
    let mut seen = vec![false; n];
    for k in 0..n {
        if seen[k] { continue; }
        let mut grp = vec![k];
        seen[k] = true;
        let mut ok = singular[org[k] as usize] == 0;
        let mut i = 0;
        while i < grp.len() {
            let a = grp[i];
            for l in &links[a] {
                let b = l.to;
                let rb = l.r * rk[a];
                let tb = l.r * tk[a] + l.t;
                if !seen[b] { seen[b] = true; rk[b] = rb; tk[b] = tb; grp.push(b); }
                else if (rk[b] - rb).norm() > 1e-9 || (tk[b] - tb).norm() > 1e-6 { ok = false; }
            }
            i += 1;
        }
        if !ok { continue; }
        for &a in &grp { var[a] = first.len() as i32; }
        first.push(k);
        copies.push(grp.len());
    }
    let m = first.len();

    // the copy k at z: rk z + tk, or where it is when fixed
    let expand = |z: &Col<f64>| -> Vec<[f64; 2]> {
        let mut uv = vec![[0f64; 2]; n];
        for k in 0..n {
            let p = if var[k] < 0 { uv_at(k) } else { rk[k] * cx(z[2 * var[k] as usize], z[2 * var[k] as usize + 1]) + tk[k] };
            uv[k] = [p.re, p.im];
        }
        uv
    };
    let shrink = |uv: &[[f64; 2]]| -> Col<f64> {
        let mut z = Col::<f64>::zeros(2 * m);
        for g in 0..m { z[2 * g] = uv[first[g]][0]; z[2 * g + 1] = uv[first[g]][1]; }
        z
    };

    // 3: per face, the gradient of the three corner values in a frame of the face, and the face area
    let (f1, f2, _f3) = local_basis(&hc.pos, &hc.idx);
    let g = grad(&hc.pos, &hc.idx);
    let gr = transpose(g.as_ref()); // Eigen::SparseMatrix<double, RowMajor> Gr = G: column r of gr is row r of G
    let mut dx = vec![[0f64; 3]; nf];
    let mut dy = vec![[0f64; 3]; nf];
    for f in 0..nf {
        let b1 = f1[f].to_array();
        let b2 = f2[f].to_array();
        for d in 0..3 {
            let row = f + d * nf;
            let cols = gr.as_ref().row_idx_of_col_raw(row);
            let vals = gr.as_ref().val_of_col(row);
            for (&c, &val) in cols.iter().zip(vals.iter()) {
                let c = c as i32;
                let j = if c == hc.idx[f][0] { 0 } else if c == hc.idx[f][1] { 1 } else { 2 };
                dx[f][j] += b1[d] * val;
                dy[f][j] += b2[d] * val;
            }
        }
    }
    let area = doublearea(&hc.pos, &hc.idx);

    // 4: the system over the free points. its pattern is fixed: where each face adds its 6 x 6 block is found once
    let mut k_mat: SprsD = {
        let mut t: Vec<TripD> = vec![];
        for g in 0..2 * m { t.push(Triplet::new(g, g, 0.)); }
        for f in 0..nf {
            for a in 0..3 {
                for b in 0..3 {
                    let ga = var[hc.idx[f][a] as usize];
                    let gb = var[hc.idx[f][b] as usize];
                    if ga < 0 || gb < 0 { continue; }
                    for i in 0..2 { for j in 0..2 { t.push(Triplet::new(2 * ga as usize + i, 2 * gb as usize + j, 0.)); } }
                }
            }
        }
        from_triplets(2 * m, 2 * m, &t)
    };
    let slot = |k_mat: &SprsD, r: usize, c: usize| -> i32 {
        let col_ptr = k_mat.symbolic().col_ptr();
        let rows = &k_mat.symbolic().row_idx()[col_ptr[c]..col_ptr[c + 1]];
        (col_ptr[c] + rows.partition_point(|&x| x < r)) as i32
    };
    let mut slots: Vec<[i32; 36]> = vec![[0; 36]; nf];
    for f in 0..nf {
        for a in 0..6 {
            for b in 0..6 {
                let ga = var[hc.idx[f][a / 2] as usize];
                let gb = var[hc.idx[f][b / 2] as usize];
                slots[f][a * 6 + b] = if ga < 0 || gb < 0 { -1 } else { slot(&k_mat, 2 * ga as usize + a % 2, 2 * gb as usize + b % 2) };
            }
        }
    }
    let mut diag = vec![0i32; 2 * m];
    for g in 0..2 * m { diag[g] = slot(&k_mat, g, g); }

    // the jacobian of a face, rows u, v and columns along the frame
    let jacobian = |uv: &[[f64; 2]], f: usize| -> Mat2 {
        let mut j_: Mat2 = MAT2_ZERO;
        for j in 0..3 {
            let u = uv[hc.idx[f][j] as usize];
            let d = [dx[f][j], dy[f][j]];
            for r in 0..2 { for c in 0..2 { j_[r][c] += d[c] * u[r]; } }
        }
        j_
    };
    let energy = |uv: &[[f64; 2]]| -> f64 {
        let mut e = 0.;
        for f in 0..nf {
            let j_ = jacobian(uv, f);
            let det = mat2_det(&j_);
            e += if det > 0. { area[f] * squared_norm2(&j_) * (1. + 1. / (det * det)) } else { f64::INFINITY };
        }
        e
    };

    let proximal = 1e-4; // as in libigl
    let mut llt = SparseLlt::<f64>::new();   // fast, but a nearly degenerate face can make it fail
    let mut ldlt = SparseLdlt::<f64>::new(); // then this one
    let mut z = shrink(uv0);
    let mut uv = expand(&z);
    let mut e = energy(&uv);
    let mut small = 0;

    for it in 0..max_iter {
        // 5: per face, the closest rotation R and the weights W of the quadratic proxy (the local step), and its
        //    block of the weighted least squares || W (J - R) ||^2 over the free points of the face (the global step)
        for x in k_mat.val_mut().iter_mut() { *x = 0.; }
        let mut rhs = Col::<f64>::zeros(2 * m);
        for f in 0..nf {
            let mut c_: Mat2 = MAT2_ZERO;            // the part of J that does not move
            let mut b_ = [[0f64; 6]; 4];             // J00, J01, J10, J11 by the free points
            for j in 0..3 {
                let k = hc.idx[f][j] as usize;
                let d = [dx[f][j], dy[f][j]];
                if var[k] < 0 {
                    for r in 0..2 { for c in 0..2 { c_[r][c] += d[c] * uv[k][r]; } }
                    continue;
                }
                let tv = [tk[k].re, tk[k].im];
                for r in 0..2 { for c in 0..2 { c_[r][c] += d[c] * tv[r]; } }
                let rk_m: Mat2 = [[rk[k].re, -rk[k].im], [rk[k].im, rk[k].re]];
                for c in 0..2 {
                    b_[0][2 * j + c] = d[0] * rk_m[0][c];
                    b_[1][2 * j + c] = d[1] * rk_m[0][c];
                    b_[2][2 * j + c] = d[0] * rk_m[1][c];
                    b_[3][2 * j + c] = d[1] * rk_m[1][c];
                }
            }

            let ps = polar_svd(&jacobian(&uv, f));
            let (r_, u_, s) = (ps.r, ps.u, ps.s);
            let mut w = [0f64; 2];
            for i in 0..2 { w[i] = if (s[i] - 1.).abs() < 1e-8 { 2. } else { ((s[i] - s[i].powf(-3.)) / (s[i] - 1.)).sqrt() }; }
            // W = U * w.asDiagonal() * U.transpose()
            let ud: Mat2 = [[u_[0][0] * w[0], u_[0][1] * w[1]], [u_[1][0] * w[0], u_[1][1] * w[1]]];
            let w_m = mat2_mul(&ud, &mat2_transpose(&u_));

            // the residual W (J - R) = Wb B z - vec(W (R - C)), row by row
            let mut wb = [[0f64; 4]; 4];
            wb[0][0] = w_m[0][0]; wb[0][2] = w_m[0][1]; wb[1][1] = w_m[0][0]; wb[1][3] = w_m[0][1];
            wb[2][0] = w_m[1][0]; wb[2][2] = w_m[1][1]; wb[3][1] = w_m[1][0]; wb[3][3] = w_m[1][1];
            let rc: Mat2 = [[r_[0][0] - c_[0][0], r_[0][1] - c_[0][1]], [r_[1][0] - c_[1][0], r_[1][1] - c_[1][1]]];
            let d_m = mat2_mul(&w_m, &rc);
            let bv = [d_m[0][0], d_m[0][1], d_m[1][0], d_m[1][1]];
            // A = Wb * B: fused chain over k
            let mut a_ = [[0f64; 6]; 4];
            for i in 0..4 {
                for j in 0..6 {
                    let mut acc = wb[i][0] * b_[0][j];
                    for kk in 1..4 { acc = wb[i][kk].mul_add(b_[kk][j], acc); }
                    a_[i][j] = acc;
                }
            }
            // Kf = area(f) * A.transpose() * A, rf = area(f) * A.transpose() * b
            let ar = area[f];
            let mut kf = [[0f64; 6]; 6];
            let mut rf = [0f64; 6];
            for i in 0..6 {
                let l = [ar * a_[0][i], ar * a_[1][i], ar * a_[2][i], ar * a_[3][i]];
                for j in 0..6 {
                    kf[i][j] = (l[0] * a_[0][j] + l[2] * a_[2][j]) + (l[1] * a_[1][j] + l[3] * a_[3][j]);
                }
                rf[i] = (l[0] * bv[0] + l[2] * bv[2]) + (l[1] * bv[1] + l[3] * bv[3]);
            }

            let kv = k_mat.val_mut();
            for a in 0..6 {
                let ga = var[hc.idx[f][a / 2] as usize];
                if ga < 0 { continue; }
                rhs[2 * ga as usize + a % 2] += rf[a];
                for c in 0..6 {
                    let p = slots[f][a * 6 + c];
                    if p < 0 { continue; }
                    kv[p as usize] += kf[a][c];
                }
            }
        }
        {
            let kv = k_mat.val_mut();
            for g in 0..m {
                for i in 0..2 {
                    kv[diag[2 * g + i] as usize] += proximal * copies[g] as f64;
                    rhs[2 * g + i] += proximal * copies[g] as f64 * z[2 * g + i];
                }
            }
        }
        if it == 0 { llt.analyze_pattern(k_mat.as_ref()); }
        llt.factorize(k_mat.as_ref());
        let zd = if llt.info() == ComputationInfo::Success { llt.solve(&rhs) } else {
            ldlt.compute(k_mat.as_ref());
            metriko_check!(ldlt.info() == ComputationInfo::Success, "the slim system is singular");
            ldlt.solve(&rhs)
        };
        let dst = expand(&zd);

        // 6: the step toward it that flips no face, and the energy there
        let e_new = flip_avoiding_line_search(&hc.idx, &mut uv, &dst, &energy, e);
        z = shrink(&uv);
        small = if (e - e_new).abs() < rel_tol * e_new.abs() { small + 1 } else { 0 }; // a single short step is the flip barrier, not convergence
        e = e_new;
        if small == 2 { break; }
    }
    uv
}
