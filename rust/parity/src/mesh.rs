//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of metriko_web/wasm/quadpatch.h load_mesh, with libigl v2.6.4 read_triangle_mesh (obj path: readOBJ,
//--- list_to_matrix, polygon_corners, polygons_to_triangles), is_edge_manifold, is_vertex_manifold and
//--- boundary_facets (MPL-2.0).

//! The input mesh as the C++ cli reads it. Only .obj is ported (the models are obj); off / ply / stl return an error.

use glam::DVec3;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::Path;

#[derive(Clone, Debug, Default)]
pub struct InputMesh {
    pub v: Vec<DVec3>,
    pub f: Vec<[i32; 3]>,
}

/// quadpatch.h `load_mesh`: read, then reject meshes metriko cannot take, with the same messages in the same order
pub fn load_mesh(path: &Path) -> Result<InputMesh, String> {
    let ext = path.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
    if ext != "obj" {
        // DEVIATION: igl::read_triangle_mesh also reads off / ply / stl / wrl / mesh / msh; only obj is ported
        return Err(format!("failed to read the mesh: .{ext} is not supported by the Rust loader (obj only)"));
    }
    let bytes = std::fs::read(path).map_err(|_| "failed to read the mesh".to_string())?;
    let (vv, vf) = read_obj(&bytes).ok_or_else(|| "failed to read the mesh".to_string())?;

    // read_triangle_mesh: list_to_matrix(vV, V) (rows truncated to 3, at least 3 each), fan triangulation of vF
    let mut m = InputMesh::default();
    if !vv.is_empty() {
        m.v = vv.iter().map(|v| DVec3::new(v[0], v[1], v[2])).collect();
        for p in &vf {
            for c in 1..p.len().saturating_sub(1) { m.f.push([p[0], p[c], p[c + 1]]); }
        }
    }

    if m.v.is_empty() || m.f.is_empty() { return Err("the mesh is empty".into()); }

    let nv = m.v.len();
    if m.f.iter().flatten().any(|&i| i < 0 || i as usize >= nv) {
        // DEVIATION: the C++ indexes out of bounds here (undefined behavior); report it instead
        return Err("the mesh has face indices out of range".into());
    }
    let mut used = vec![false; nv];
    for &i in m.f.iter().flatten() { used[i as usize] = true; }
    let n = used.iter().filter(|u| !**u).count();
    if n > 0 { return Err(format!("the mesh has {n} unreferenced vertices")); }

    if !is_edge_manifold(&m.f) { return Err("the mesh is not edge-manifold".into()); }
    if !is_vertex_manifold(&m.f) { return Err("the mesh is not vertex-manifold".into()); }
    let nb = boundary_facets(&m.f);
    if nb > 0 { return Err(format!("the mesh has {nb} boundary edges: metriko needs a closed mesh")); }
    Ok(m)
}

/// fgets(line, 2048, fp): lines of at most 2047 bytes, the newline included
fn fgets_lines(b: &[u8]) -> Vec<&[u8]> {
    const MAX: usize = 2047;
    let mut res = vec![];
    let mut i = 0;
    while i < b.len() {
        let mut j = i;
        while j < b.len() && j - i < MAX {
            j += 1;
            if b[j - 1] == b'\n' { break; }
        }
        res.push(&b[i..j]);
        i = j;
    }
    res
}

fn is_space(c: u8) -> bool { matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c) }

/// sscanf("%s"): the first whitespace-delimited word
fn first_word(s: &[u8]) -> Option<&[u8]> {
    let i = s.iter().position(|&c| !is_space(c))?;
    let j = s[i..].iter().position(|&c| is_space(c)).map_or(s.len(), |k| i + k);
    Some(&s[i..j])
}

/// the longest prefix of `s` that is a decimal floating point number: [+-]?(d+[.d*]|.d+)([eE][+-]?d+)?
fn float_prefix(s: &[u8]) -> usize {
    let mut i = 0;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') { i += 1; }
    let d0 = i;
    while i < s.len() && s[i].is_ascii_digit() { i += 1; }
    let mut nd = i - d0;
    if i < s.len() && s[i] == b'.' {
        i += 1;
        let d1 = i;
        while i < s.len() && s[i].is_ascii_digit() { i += 1; }
        nd += i - d1;
    }
    if nd == 0 { return 0; }
    if i < s.len() && (s[i] == b'e' || s[i] == b'E') {
        let mut k = i + 1;
        if k < s.len() && (s[k] == b'+' || s[k] == b'-') { k += 1; }
        let e0 = k;
        while k < s.len() && s[k].is_ascii_digit() { k += 1; }
        if k > e0 { i = k; }
    }
    i
}

/// successive `>> double` extractions (std::istream_iterator<double>) until one fails. libc++'s num_get is emulated
/// for decimal numbers; DEVIATION: hexadecimal and inf / nan spellings are not recognized
fn stream_doubles(s: &[u8]) -> Vec<f64> {
    let mut res = vec![];
    let mut i = 0;
    loop {
        while i < s.len() && is_space(s[i]) { i += 1; }
        let n = float_prefix(&s[i..]);
        if n == 0 { break; }
        match std::str::from_utf8(&s[i..i + n]).ok().and_then(|t| t.parse::<f64>().ok()) {
            Some(x) => res.push(x),
            None => break,
        }
        i += n;
    }
    res
}

/// sscanf("%lf %lf %lf"): the number of conversions
fn scanf_doubles(s: &[u8], max: usize) -> usize {
    let mut i = 0;
    let mut k = 0;
    while k < max {
        while i < s.len() && is_space(s[i]) { i += 1; }
        let n = float_prefix(&s[i..]);
        if n == 0 { break; }
        i += n;
        k += 1;
    }
    k
}

/// sscanf "%ld" at the start of `s` (leading whitespace skipped): the value and the rest
fn scan_long(s: &[u8]) -> Option<(i64, &[u8])> {
    let mut i = 0;
    while i < s.len() && is_space(s[i]) { i += 1; }
    let b = i;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') { i += 1; }
    let d = i;
    while i < s.len() && s[i].is_ascii_digit() { i += 1; }
    if i == d { return None; }
    let v = std::str::from_utf8(&s[b..i]).ok()?.parse::<i64>().ok()?;
    Some((v, &s[i..]))
}

/// sscanf(word, pattern) for patterns made of %ld separated by literal text: the number of conversions
fn scan_pattern<'a>(word: &'a [u8], seps: &[&[u8]], out: &mut [i64; 3]) -> usize {
    let mut s = word;
    let mut n = 0;
    for k in 0..=seps.len() {
        match scan_long(s) {
            Some((v, rest)) => { out[n] = v; n += 1; s = rest; }
            None => return n,
        }
        if k < seps.len() {
            if !s.starts_with(seps[k]) { return n; }
            s = &s[seps[k].len()..];
        }
    }
    n
}

/// libigl readOBJ (FILE* overload): vertices (all their coordinates) and polygons (0-based); None when it fails
fn read_obj(bytes: &[u8]) -> Option<(Vec<Vec<f64>>, Vec<Vec<i32>>)> {
    let mut v: Vec<Vec<f64>> = vec![];
    let mut f: Vec<Vec<i32>> = vec![];
    let (mut ntc, mut nn) = (0usize, 0usize);
    for line in fgets_lines(bytes) {
        let Some(ty) = first_word(line) else { continue };
        // the rest of the line right after the type, counted from the start of the line as libigl does
        let l = &line[ty.len().min(line.len())..];
        match ty {
            b"v" => {
                let mut x = stream_doubles(&line[1..]);
                if x.len() < 3 { return None; }
                x.truncate(3); // read_triangle_mesh resizes every vertex to at most 3
                v.push(x);
            }
            b"vn" => { if scanf_doubles(l, 3) != 3 { return None; } nn += 1; }
            b"vt" => { let c = scanf_doubles(l, 3); if c != 2 && c != 3 { return None; } ntc += 1; }
            b"f" => {
                let shift = |i: i64, n: usize| -> i32 { let i = i as i32; if i < 0 { i + n as i32 } else { i - 1 } };
                let (mut face, mut ftc, mut fnn) = (vec![], 0usize, 0usize);
                let mut rest = l;
                while let Some(word) = first_word(rest) {
                    let off = word.as_ptr() as usize - rest.as_ptr() as usize + word.len();
                    rest = &rest[off..];
                    let mut x = [0i64; 3];
                    if scan_pattern(word, &[b"/", b"/"], &mut x) == 3 {
                        face.push(shift(x[0], v.len())); let _ = (shift(x[1], ntc), shift(x[2], nn)); ftc += 1; fnn += 1;
                    } else if scan_pattern(word, &[b"/"], &mut x) == 2 {
                        face.push(shift(x[0], v.len())); ftc += 1;
                    } else if scan_pattern(word, &[b"//"], &mut x) == 2 {
                        face.push(shift(x[0], v.len())); fnn += 1;
                    } else if scan_pattern(word, &[], &mut x) == 1 {
                        face.push(shift(x[0], v.len()));
                    } else {
                        return None;
                    }
                }
                let n = face.len();
                if n > 0 && (fnn == 0 || fnn == n) && (ftc == 0 || ftc == n) { f.push(face); } else { return None; }
            }
            _ => {} // usemtl, comments, groups, smoothing, mtllib and unknown lines are ignored
        }
    }
    Some((v, f))
}

/// igl::is_edge_manifold: no undirected edge is shared by more than two faces
pub fn is_edge_manifold(f: &[[i32; 3]]) -> bool {
    let mut count: BTreeMap<(i32, i32), usize> = BTreeMap::new();
    for t in f { for c in 0..3 { let (a, b) = (t[c], t[(c + 1) % 3]); *count.entry((a.min(b), a.max(b))).or_default() += 1; } }
    count.values().all(|&c| c <= 2)
}

/// igl::is_vertex_manifold: the faces around every vertex 0..=max index form one fan, connected through shared edges
pub fn is_vertex_manifold(f: &[[i32; 3]]) -> bool {
    let n = f.iter().flatten().copied().max().map_or(0, |m| m + 1) as usize;
    // triangle_triangle_adjacency: the other faces sharing each edge of a face
    let mut by_edge: BTreeMap<(i32, i32), Vec<usize>> = BTreeMap::new();
    for (i, t) in f.iter().enumerate() { for c in 0..3 { let (a, b) = (t[c], t[(c + 1) % 3]); by_edge.entry((a.min(b), a.max(b))).or_default().push(i); } }
    let mut v2f: Vec<Vec<usize>> = vec![vec![]; n];
    for (i, t) in f.iter().enumerate() { for &v in t { v2f[v as usize].push(i); } }
    let mut all = true;
    for v in 0..n {
        let ring: BTreeSet<usize> = v2f[v].iter().copied().collect();
        if ring.is_empty() { all = false; continue; }
        let mut seen = BTreeSet::new();
        let mut q = VecDeque::from([*ring.iter().next().unwrap()]);
        while let Some(g) = q.pop_front() {
            if !seen.insert(g) { continue; }
            let t = f[g];
            for c in 0..3 {
                let (a, b) = (t[c], t[(c + 1) % 3]);
                for &h in &by_edge[&(a.min(b), a.max(b))] {
                    if h != g && !seen.contains(&h) && f[h].contains(&(v as i32)) { q.push_back(h); }
                }
            }
        }
        all &= ring.len() == seen.len();
    }
    all
}

/// igl::boundary_facets on triangles: the number of undirected edges used by exactly one face
pub fn boundary_facets(f: &[[i32; 3]]) -> usize {
    let mut count: BTreeMap<(i32, i32), usize> = BTreeMap::new();
    for t in f { for c in 0..3 { let (a, b) = (t[c], t[(c + 1) % 3]); *count.entry((a.min(b), a.max(b))).or_default() += 1; } }
    count.values().filter(|&&c| c == 1).count()
}

/// lib.h `grid_unit`: the bounding box diagonal (Eigen: column max - column min, sqrt((x*x + y*y) + z*z)) times scale
pub fn grid_unit(v: &[DVec3], scale: f64) -> f64 {
    let (mut lo, mut hi) = (v[0], v[0]);
    for p in &v[1..] { lo = lo.min(*p); hi = hi.max(*p); }
    let d = hi - lo;
    (d.x * d.x + d.y * d.y + d.z * d.z).sqrt() * scale
}
