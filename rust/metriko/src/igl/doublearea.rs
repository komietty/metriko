//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/doublearea.cpp (MPL-2.0).

use crate::common::typedef::*;

/// `igl::doublearea(V, F, dblA)` for a triangle mesh in 3D (the `dim == 3` case): twice the area of every face, as
/// the norm of its three projected double areas.
pub fn doublearea(v: &[Row3d], f: &[[i32; 3]]) -> Vec<f64> {
    // Projected area helper
    let proj_doublearea = |x: usize, y: usize, fi: usize| -> f64 {
        let p0 = v[f[fi][0] as usize].to_array();
        let p1 = v[f[fi][1] as usize].to_array();
        let p2 = v[f[fi][2] as usize].to_array();
        let rx = p0[x] - p2[x];
        let sx = p1[x] - p2[x];
        let ry = p0[y] - p2[y];
        let sy = p1[y] - p2[y];
        rx * sy - ry * sx
    };

    let m = f.len();
    let mut dbl_a = vec![0.; m];
    for fi in 0..m {
        for d in 0..3 {
            let dbl_ad = proj_doublearea(d, (d + 1) % 3, fi);
            dbl_a[fi] += dbl_ad * dbl_ad;
        }
    }
    for a in dbl_a.iter_mut() { *a = a.sqrt(); }
    dbl_a
}
