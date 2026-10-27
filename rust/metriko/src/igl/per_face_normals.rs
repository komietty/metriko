//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/per_face_normals.cpp (MPL-2.0).

use crate::common::typedef::Row3d;

/// `igl::per_face_normals(V, F, N)`: the unit normal of every triangle, `Z = (0, 0, 0)` for a degenerate one.
pub fn per_face_normals(v: &[Row3d], f: &[[i32; 3]]) -> Vec<Row3d> {
    per_face_normals_z(v, f, Row3d::ZERO)
}

/// `igl::per_face_normals(V, F, Z, N)`
pub fn per_face_normals_z(v: &[Row3d], f: &[[i32; 3]], z: Row3d) -> Vec<Row3d> {
    f.iter().map(|t| {
        let v1 = v[t[1] as usize] - v[t[0] as usize];
        let v2 = v[t[2] as usize] - v[t[0] as usize];
        let n = v1.cross(v2);
        let r = n.length();
        if r == 0. { z } else { n / r }
    }).collect()
}
