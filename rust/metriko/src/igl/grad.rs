//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/grad.cpp (MPL-2.0).

use faer::sparse::Triplet;
use crate::common::typedef::*;
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets};

/// `igl::grad(V, F, G, uniform = false)` for a triangle mesh in 3D (`grad_tri`): the `3m x n` operator whose rows
/// `f + d*m` give the d-th component of the gradient on face f of a function given per vertex.
pub fn grad(v: &[Row3d], f: &[[i32; 3]]) -> SprsD {
    // Number of faces
    let m = f.len();
    // Number of vertices
    let nv = v.len();
    // Number of dimensions
    let dims = 3;
    let mut eperp21 = vec![Row3d::ZERO; m];
    let mut eperp13 = vec![Row3d::ZERO; m];

    for i in 0..m {
        // renaming indices of vertices of triangles for convenience
        let i1 = f[i][0] as usize;
        let i2 = f[i][1] as usize;
        let i3 = f[i][2] as usize;

        // #F x 3 matrices of triangle edge vectors, named after opposite vertices
        let v32 = v[i3] - v[i2];
        let v13 = v[i1] - v[i3];
        let v21 = v[i2] - v[i1];
        let n = v32.cross(v13);
        // area of parallelogram is twice area of triangle
        let dbl_a = n.dot(n).sqrt();
        // now normalize normals to get unit normals
        let u = n / dbl_a;

        // rotate each vector 90 degrees around normal
        let norm21 = v21.dot(v21).sqrt();
        let norm13 = v13.dot(v13).sqrt();
        let mut e21 = u.cross(v21);
        e21 = e21 / e21.dot(e21).sqrt();
        e21 *= norm21 / dbl_a;
        let mut e13 = u.cross(v13);
        e13 = e13 / e13.dot(e13).sqrt();
        e13 *= norm13 / dbl_a;
        eperp21[i] = e21;
        eperp13[i] = e13;
    }

    // create sparse gradient operator matrix
    let mut gijv: Vec<TripD> = Vec::with_capacity(4 * dims * m);
    for fi in 0..m {
        let e13 = eperp13[fi].to_array();
        let e21 = eperp21[fi].to_array();
        for d in 0..dims {
            gijv.push(Triplet::new(fi + d * m, f[fi][1] as usize, e13[d]));
            gijv.push(Triplet::new(fi + d * m, f[fi][0] as usize, -e13[d]));
            gijv.push(Triplet::new(fi + d * m, f[fi][2] as usize, e21[d]));
            gijv.push(Triplet::new(fi + d * m, f[fi][0] as usize, -e21[d]));
        }
    }
    // setFromTriplets: the two entries on (f + d*m, F(f,0)) are summed
    from_triplets(dims * m, nv, &gijv)
}
