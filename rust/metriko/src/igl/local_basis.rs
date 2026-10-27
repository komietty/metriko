//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/local_basis.cpp (MPL-2.0).

use crate::common::typedef::*;

/// `igl::local_basis(V, F, B1, B2, B3)`: an orthonormal frame per face, returns `(B1, B2, B3)`.
/// B1 along the first edge, B3 the normal, B2 completes it.
pub fn local_basis(v: &[Row3d], f: &[[i32; 3]]) -> (Vec<Row3d>, Vec<Row3d>, Vec<Row3d>) {
    let mut b1 = vec![Row3d::ZERO; f.len()];
    let mut b2 = vec![Row3d::ZERO; f.len()];
    let mut b3 = vec![Row3d::ZERO; f.len()];

    for i in 0..f.len() {
        let v1 = normalized(v[f[i][1] as usize] - v[f[i][0] as usize]);
        let t = v[f[i][2] as usize] - v[f[i][0] as usize];
        let v3 = normalized(v1.cross(t));
        let v2 = normalized(v1.cross(v3));

        b1[i] = v1;
        b2[i] = -v2;
        b3[i] = v3;
    }
    (b1, b2, b3)
}
