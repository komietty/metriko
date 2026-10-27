//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/edge_topology.cpp (MPL-2.0).

/// `igl::edge_topology(V, F, EV, FE, EF)` for a triangle mesh: returns `(EV, FE, EF)`.
/// EV: the two vertices of every edge (smaller id first), FE: the edge of every face corner `(f, i)` = edge
/// `F(f,i) -> F(f,(i+1)%3)`, EF: the faces on both sides of every edge (`-1` on a boundary), the first being the
/// face in which the edge runs from EV(e,0) to EV(e,1).
pub fn edge_topology(n_v: usize, f: &[[i32; 3]]) -> (Vec<[i32; 2]>, Vec<[i32; 3]>, Vec<[i32; 2]>) {
    if n_v == 0 || f.is_empty() {
        return (vec![], vec![], vec![]);
    }
    let mut ett: Vec<[i32; 4]> = Vec::with_capacity(f.len() * 3);
    for (fi, tri) in f.iter().enumerate() {
        for i in 0..3 {
            let mut v1 = tri[i];
            let mut v2 = tri[(i + 1) % 3];
            if v1 > v2 { std::mem::swap(&mut v1, &mut v2); }
            ett.push([v1, v2, fi as i32, i as i32]);
        }
    }
    // std::sort of vector<vector<int>>: lexicographic, all rows distinct
    ett.sort_unstable();

    let mut en = 1usize; // the last is always counted
    for i in 0..ett.len() - 1 {
        if !(ett[i][0] == ett[i + 1][0] && ett[i][1] == ett[i + 1][1]) { en += 1; }
    }

    let mut ev = vec![[-1i32; 2]; en];
    let mut fe = vec![[-1i32; 3]; f.len()];
    let mut ef = vec![[-1i32; 2]; en];
    let mut en = 0usize;

    let mut i = 0usize;
    while i < ett.len() {
        if i == ett.len() - 1 || !(ett[i][0] == ett[i + 1][0] && ett[i][1] == ett[i + 1][1]) {
            let r1 = ett[i];
            ev[en] = [r1[0], r1[1]];
            ef[en][0] = r1[2];
            fe[r1[2] as usize][r1[3] as usize] = en as i32;
        } else {
            let r1 = ett[i];
            let r2 = ett[i + 1];
            ev[en] = [r1[0], r1[1]];
            ef[en][0] = r1[2];
            ef[en][1] = r2[2];
            fe[r1[2] as usize][r1[3] as usize] = en as i32;
            fe[r2[2] as usize][r2[3] as usize] = en as i32;
            i += 1; // skip the next one
        }
        en += 1;
        i += 1;
    }

    for i in 0..ef.len() {
        let fid = ef[i][0];
        let mut flip = true;
        for j in 0..3 {
            if f[fid as usize][j] == ev[i][0] && f[fid as usize][(j + 1) % 3] == ev[i][1] { flip = false; }
        }
        if flip { ef[i].swap(0, 1); }
    }
    (ev, fe, ef)
}
