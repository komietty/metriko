//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/vertex_triangle_adjacency.cpp (MPL-2.0).

/// `igl::vertex_triangle_adjacency(n, F, VF, VFi)` (also the `(V, F, VF, VFi)` overload with `n = V.rows()`):
/// returns `(VF, VFi)`, for every vertex the faces around it and the corner of the vertex in each, in face order.
pub fn vertex_triangle_adjacency(n: usize, f: &[[i32; 3]]) -> (Vec<Vec<i32>>, Vec<Vec<i32>>) {
    let mut vf: Vec<Vec<i32>> = vec![vec![]; n];
    let mut vfi: Vec<Vec<i32>> = vec![vec![]; n];
    for (fi, tri) in f.iter().enumerate() {
        for i in 0..3 {
            vf[tri[i] as usize].push(fi as i32);
            vfi[tri[i] as usize].push(i as i32);
        }
    }
    (vf, vfi)
}

/// `igl::vertex_triangle_adjacency(F, n, VF, NI)`: the compressed form, returns `(VF, NI)`. The faces around vertex
/// `i` are `VF[NI[i] .. NI[i + 1]]`, in face order.
pub fn vertex_triangle_adjacency_csr(f: &[[i32; 3]], n: usize) -> (Vec<i32>, Vec<i32>) {
    // vfd  #V list so that vfd(i) contains the vertex-face degree (number of faces incident on vertex i)
    let mut vfd = vec![0i32; n];
    for tri in f {
        for j in 0..3 { vfd[tri[j] as usize] += 1; }
    }
    // igl::cumsum(vfd, 1, NI), then prepend a zero
    let mut ni = vec![0i32; n + 1];
    let mut sum = 0i32;
    for i in 0..n {
        sum += vfd[i];
        ni[i + 1] = sum;
    }
    // vfd now acts as a counter
    let mut vfd = ni.clone();

    let mut vf = vec![0i32; 3 * f.len()];
    for (i, tri) in f.iter().enumerate() {
        for j in 0..3 {
            vf[vfd[tri[j] as usize] as usize] = i as i32;
            vfd[tri[j] as usize] += 1;
        }
    }
    (vf, ni)
}
