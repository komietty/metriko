//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/adjacency_list.cpp (MPL-2.0).

/// `igl::adjacency_list(F, A, sorted = false)` for a triangle mesh: the vertices adjacent to every vertex, ascending
/// and without duplicates; with `sorted`, reordered counter-clockwise (input assumed to be manifold).
pub fn adjacency_list(f: &[[i32; 3]], sorted: bool) -> Vec<Vec<i32>> {
    let max = f.iter().flat_map(|t| t.iter()).copied().max().unwrap_or(-1);
    let mut a: Vec<Vec<i32>> = vec![vec![]; (max + 1) as usize];

    // Loop over faces
    for tri in f {
        // Loop over this face
        for j in 0..3 {
            // Get indices of edge: s --> d
            let s = tri[j];
            let d = tri[(j + 1) % 3];
            a[s as usize].push(d);
            a[d as usize].push(s);
        }
    }

    // Remove duplicates
    for ai in a.iter_mut() {
        ai.sort_unstable();
        ai.dedup();
    }

    // If needed, sort every VV
    if sorted {
        // for every vertex v store a set of ordered edges not incident to v that belongs to triangle incident on v.
        let mut sr: Vec<Vec<[i32; 2]>> = vec![vec![]; a.len()];
        for tri in f {
            for j in 0..3 {
                // Get indices of edge: s --> d
                let s = tri[j];
                let d = tri[(j + 1) % 3];
                // Get index of opposing vertex v
                let v = tri[(j + 2) % 3];
                sr[s as usize].push([d, v]);
            }
        }

        for v in 0..sr.len() {
            let vv = &mut a[v];
            if vv.is_empty() { continue; }
            let sr = &sr[v];

            let mut pn = sr.clone();

            // Compute previous/next for every element in sr
            for i in 0..sr.len() {
                let a_ = sr[i][0];
                let b_ = sr[i][1];

                // search for previous
                let mut p = -1;
                for j in 0..sr.len() { if sr[j][1] == a_ { p = j as i32; } }
                pn[i][0] = p;

                // search for next
                let mut n = -1;
                for j in 0..sr.len() { if sr[j][0] == b_ { n = j as i32; } }
                pn[i][1] = n;
            }

            // assume manifoldness (look for beginning of a single chain)
            let mut c = 0usize;
            for _ in 0..=sr.len() {
                if pn[c][0] != -1 { c = pn[c][0] as usize; }
            }

            if pn[c][0] == -1 {
                // border case: finally produce the new vv relation
                for j in 0..sr.len() {
                    vv[j] = sr[c][0];
                    if pn[c][1] != -1 { c = pn[c][1] as usize; }
                }
                *vv.last_mut().unwrap() = sr[c][1];
            } else {
                // finally produce the new vv relation
                for j in 0..sr.len() {
                    vv[j] = sr[c][0];
                    c = pn[c][1] as usize;
                }
            }
        }
    }
    a
}
