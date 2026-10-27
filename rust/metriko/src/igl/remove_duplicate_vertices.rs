//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/remove_duplicate_vertices.cpp (MPL-2.0),
//--- with the parts of round.cpp, sortrows.cpp and unique_rows.cpp it calls.

use crate::common::typedef::Row3d;
use crate::common::stdcxx;

/// `igl::round(r)`: half away from zero, as `floor(r + 0.5)` / `ceil(r - 0.5)`
pub fn round(r: f64) -> f64 {
    if r > 0.0 { (r + 0.5).floor() } else { (r - 0.5).ceil() }
}

/// `igl::sortrows(X, ascending = true, Y, IX)` on rows of 3: returns `IX`
fn sortrows_ascending(x: &[[f64; 3]]) -> Vec<i32> {
    let mut ix: Vec<i32> = (0..x.len() as i32).collect();
    // std::sort with index_less_than: lexicographic on the rows
    stdcxx::sort_by(&mut ix, |&i, &j| {
        for c in 0..3 {
            if x[i as usize][c] < x[j as usize][c] { return true; }
            else if x[j as usize][c] < x[i as usize][c] { return false; }
        }
        false
    });
    ix
}

/// `igl::unique_rows(A, C, IA, IC)` on rows of 3: returns `(IA, IC)`
fn unique_rows(a: &[[f64; 3]]) -> (Vec<i32>, Vec<i32>) {
    let im = sortrows_ascending(a);
    let sort_a: Vec<[f64; 3]> = im.iter().map(|&i| a[i as usize]).collect();
    let num_rows = sort_a.len();
    let mut v_ia: Vec<usize> = (0..num_rows).collect();
    // std::unique with index_equal: keeps the first of every run of equal rows
    v_ia.dedup_by(|j, i| sort_a[*i] == sort_a[*j]);

    let mut ic = vec![0i32; a.len()];
    {
        let mut j = 0usize;
        for i in 0..num_rows {
            if sort_a[v_ia[j]] != sort_a[i] { j += 1; }
            ic[im[i] as usize] = j as i32;
        }
    }
    let ia: Vec<i32> = v_ia.iter().map(|&k| im[k]).collect();
    (ia, ic)
}

/// `igl::remove_duplicate_vertices(V, epsilon, SV, SVI, SVJ)`: returns `(SV, SVI, SVJ)`
pub fn remove_duplicate_vertices(v: &[Row3d], epsilon: f64) -> (Vec<Row3d>, Vec<i32>, Vec<i32>) {
    let rows: Vec<[f64; 3]> = if epsilon > 0. {
        v.iter().map(|p| [round(p.x / epsilon), round(p.y / epsilon), round(p.z / epsilon)]).collect()
    } else {
        v.iter().map(|p| p.to_array()).collect()
    };
    let (svi, svj) = unique_rows(&rows);
    let sv = svi.iter().map(|&i| v[i as usize]).collect();
    (sv, svi, svj)
}
