//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/triangle_triangle_adjacency.cpp (MPL-2.0).

use crate::igl::vertex_triangle_adjacency::vertex_triangle_adjacency_csr;

/// `igl::triangle_triangle_adjacency(F, TT)`: `TT(f, k)` is the face across the edge `F(f,k) -> F(f,(k+1)%3)`, `-1`
/// on a boundary. The neighbor is the first face, in the vertex-face order of `F(f,k)`, that also has `F(f,(k+1)%3)`.
pub fn triangle_triangle_adjacency_tt(f: &[[i32; 3]]) -> Vec<[i32; 3]> {
    let n = f.iter().flat_map(|t| t.iter()).copied().max().unwrap_or(-1) + 1;
    let (vf, ni) = vertex_triangle_adjacency_csr(f, n as usize);
    let mut tt = vec![[-1i32; 3]; f.len()];
    // Loop over faces (igl::parallel_for: every face writes its own row only)
    for fi in 0..f.len() {
        // Loop over corners
        for k in 0..3 {
            let vi = f[fi][k];
            let vin = f[fi][(k + 1) % 3];
            // Loop over face neighbors incident on this corner
            for j in ni[vi as usize]..ni[vi as usize + 1] {
                let fne = vf[j as usize];
                // Not this face
                if fne != fi as i32 {
                    // Face neighbor also has [vi,vin] edge
                    let g = &f[fne as usize];
                    if g[0] == vin || g[1] == vin || g[2] == vin {
                        tt[fi][k] = fne;
                        break;
                    }
                }
            }
        }
    }
    tt
}

/// `igl::triangle_triangle_adjacency(F, TT, TTi)`: returns `(TT, TTi)`; `TTi(f, k)` is the edge index of the shared
/// edge inside face `TT(f, k)`.
pub fn triangle_triangle_adjacency(f: &[[i32; 3]]) -> (Vec<[i32; 3]>, Vec<[i32; 3]>) {
    let tt = triangle_triangle_adjacency_tt(f);
    let mut tti = vec![[-1i32; 3]; tt.len()];
    for fi in 0..f.len() {
        for k in 0..3 {
            let vi = f[fi][k];
            let vj = f[fi][(k + 1) % 3];
            let fne = tt[fi][k];
            if fne >= 0 {
                for kn in 0..3 {
                    let vin = f[fne as usize][kn];
                    let vjn = f[fne as usize][(kn + 1) % 3];
                    if vi == vjn && vin == vj {
                        tti[fi][k] = kn as i32;
                        break;
                    }
                }
            }
        }
    }
    (tt, tti)
}
