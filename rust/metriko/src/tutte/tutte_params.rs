//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tutte/tutte_params.h (see rust/PORTED_FROM).

use std::collections::{HashMap, HashSet, VecDeque};
use faer::Mat;
use faer::sparse::Triplet;
use crate::common::typedef::*;
use crate::common::utilities::get_quater_rot;
use crate::hmesh::hmesh::*;
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets};
use crate::solver::sparse_solver::SparseLu;
use crate::tmesh::emesh::*;
use crate::tutte::tutte::HalfData;

pub fn boundary_snap_laplacian(mesh: &Hmesh) -> SprsD {
    let mut t: Vec<TripD> = vec![];

    for v in mesh.verts() {
        if v.is_boundary() { t.push(Triplet::new(v.id as usize, v.id as usize, 1.)); }
        else {
            let mut sum = 0.;
            for h in v.adj_halfs() {
                let l = (v.pos() - h.head().pos()).length();
                let w = 1. / (l + 1e-12);
                sum += w;
                t.push(Triplet::new(v.id as usize, h.head().id as usize, -w));
            }
            t.push(Triplet::new(v.id as usize, v.id as usize, sum));
        }
    }
    from_triplets(mesh.n_v as usize, mesh.n_v as usize, &t)
}

/// `rg::equal_range(data, key, {}, proj)` on a range sorted by `proj`
fn equal_range<F: Fn(&HalfData) -> i32>(data: &[HalfData], key: i32, proj: F) -> &[HalfData] {
    let lo = data.partition_point(|d| proj(d) < key);
    let hi = data.partition_point(|d| key >= proj(d)); // !(key < proj(d))
    &data[lo..hi]
}

/// `m.coeff(i, c)` of a sparse matrix: 0 when the entry is not stored
fn coeff(m: &SprsD, i: usize, c: usize) -> f64 { m.as_ref().get(i, c).copied().unwrap_or(0.) }

/// The tutte embedding of one tquad: the faces enclosed by its rim (the halfedges annotated with it), their
/// corners welded across interior edges, the rim fixed on the boundary of the integer rectangle given by the
/// quantized lengths of its sides, everything else by the tutte (inverse edge length) laplacian.
/// Returns the uv of every corner of the tquad's faces as an `nC x 2` sparse matrix (explicit zeros kept).
pub fn embedding_tutte_for_tquad(
    tqid: i32,
    data: &[HalfData],
    hm: &Hmesh,      // the sub hm
    tm: &Emesh<'_>,  // the org tm
) -> SprsD {
    let tq_rg = equal_range(data, tqid, |d| d.tqid);

    let mut rim: HashSet<i32> = HashSet::new(); // the tquad's own boundary, as halfedge ids
    for it in tq_rg { rim.insert(it.half); }

    // find all faces
    let mut queue: VecDeque<i32> = VecDeque::new();
    let mut visit = vec![false; hm.n_f as usize];

    for it in tq_rg {
        let fid = hm.half(it.half).face().id;
        queue.push_back(fid);
        visit[fid as usize] = true;
    }

    while let Some(front) = queue.pop_front() {
        let f0 = hm.face(front);
        for h0 in f0.adj_halfs() {
            let f1 = h0.twin().face();
            if visit[f1.id as usize] || rim.contains(&h0.id) { continue; }
            queue.push_back(f1.id);
            visit[f1.id as usize] = true;
        }
    }

    // a map from sum hm cid -> rep cid, bacause in latter use reps is mutable, fan is for temp container
    let mut fan: Vec<i32> = (0..hm.n_c).collect();
    let fan_of = |fan: &[i32], mut c: i32| { while fan[c as usize] != c { c = fan[c as usize]; } c };

    for f in hm.faces() {
        if !visit[f.id as usize] { continue; }
        for h in f.adj_halfs() {
            if h.twin().is_boundary() || !visit[h.twin().face().id as usize] { continue; }
            if rim.contains(&h.id) || rim.contains(&h.twin().id) { continue; }
            let r = fan_of(&fan, h.twin().prev().crnr().id);
            let l = fan_of(&fan, h.next().crnr().id);
            fan[l as usize] = r; // corners at h.tail()
            let r = fan_of(&fan, h.twin().next().crnr().id);
            let l = fan_of(&fan, h.prev().crnr().id);
            fan[l as usize] = r; // corners at h.head()
        }
    }

    let mut fids: Vec<i32> = vec![];
    for f in hm.faces() {
        if visit[f.id as usize] { fids.push(f.id); }
    }

    // number the fans, then let every corner point straight at its number
    let mut reps = vec![-1i32; hm.n_c as usize]; // sub hm cid -> sub hm vid
    let mut gids: Vec<i32> = vec![];               // sub hm vid -> org hm vid
    for &fid in &fids {
        for j in 0..3 {
            let cid = fan_of(&fan, fid * 3 + j);
            if reps[cid as usize] < 0 { reps[cid as usize] = gids.len() as i32; gids.push(hm.idx[fid as usize][j as usize]); }
            reps[(fid * 3 + j) as usize] = reps[cid as usize];
        }
    }

    let v: Vec<Row3d> = gids.iter().map(|&g| hm.pos[g as usize]).collect();
    let mut f: Vec<[i32; 3]> = vec![[0; 3]; fids.len()];
    let mut uv_b = Mat::<f64>::zeros(gids.len(), 2);
    for i in 0..f.len() {
        for j in 0..3 {
            f[i][j] = reps[(fids[i] * 3 + j as i32) as usize];
        }
    }

    let mut dir = cx(1., 0.);
    let mut sum = cx(0., 0.);
    for i in 0..4 {
        for thid in tm.tquads[tqid as usize].thids(i) {
            let x = tm.thalfs[thid as usize].x;
            for it in equal_range(tq_rg, thid, |d| d.thid) {
                let val = x * it.v0;
                let row = reps[hm.half(it.half).next().crnr().id as usize] as usize;
                uv_b[(row, 0)] = val * dir.re + sum.re;
                uv_b[(row, 1)] = val * dir.im + sum.im;
            }
            sum += dir * x;
        }
        dir *= cx(0., 1.);
    }

    let lu = SparseLu::from_matrix(boundary_snap_laplacian(&Hmesh::new(&v, &f, true)).as_ref());
    let uv = lu.solve_mat(&uv_b);

    let mut t: Vec<TripD> = vec![];
    for f in hm.faces() {
        if !visit[f.id as usize] { continue; }
        for h in f.adj_halfs() {
            let r = reps[h.crnr().id as usize] as usize;
            t.push(Triplet::new(h.crnr().id as usize, 0, uv[(r, 0)]));
            t.push(Triplet::new(h.crnr().id as usize, 1, uv[(r, 1)]));
        }
    }
    from_triplets(hm.n_c as usize, 2, &t)
}

// take the tutte result as the input, embed it until seam intersection.
// computes halfedges to search with in the next loop at the same time.
pub fn sequential_mapping(
    hm: &Hmesh,
    uv_in: &SprsD,
    half_in: i32,
    boun: &[bool],          // flag if a halfedge is the boundary of the tquad
    seam: &[bool],          // need to be altered for new cut hmesh
    flag: &mut [bool],      // the flag to check a face is already marked
    uv_all: &mut [[f64; 2]],
) -> Vec<i32> {
    let mut next_h: Vec<i32> = vec![]; // the halfedges to the other tquad
    let mut visit = vec![false; hm.n_e as usize];
    let mut q: VecDeque<i32> = VecDeque::new();
    q.push_back(half_in);
    if flag[hm.half(half_in).face().id as usize] { return next_h; }

    while let Some(&front) = q.front() {
        let f = hm.half(front).face();
        let [c0, c1, c2] = f.crnrs();
        for c in [c0, c1, c2] {
            uv_all[c.id as usize] = [coeff(uv_in, c.id as usize, 0), coeff(uv_in, c.id as usize, 1)];
        }
        flag[f.id as usize] = true;
        q.pop_front();

        for h in f.adj_halfs() {
            let e = h.edge();
            if seam[e.id as usize] { continue; }                                   // 1: if hit the seam, just stops
            if visit[e.id as usize] { continue; }                                  // 2: if hit the visited edge, just stops
            if boun[h.id as usize] { next_h.push(h.twin().id); continue; }         // 3: if hit boundary, puts it as a bridge to the next tquad
            visit[e.id as usize] = true;                                           // 4: inside tquad. add it to the queue
            q.push_back(h.twin().id);
        }
    }
    next_h
}

// try to multiply rotation until halfedge coner values corresponds
// need to consider: is there any possibility of flip?
pub fn apply_transition(
    h: Half<'_>,        // the halfedge of unfixed side
    m0: &[[f64; 2]],    // the fixed uv information
    m1: &mut SprsD,     // the unfixed adjacent uv information
) {
    let at = |m: &SprsD, i: usize| cx(coeff(m, i, 0), coeff(m, i, 1));
    let c0 = h.twin().next().crnr();
    let c1 = h.twin().prev().crnr();
    let uv0 = cx(m0[c0.id as usize][0], m0[c0.id as usize][1]);
    let uv0a = cx(m0[c1.id as usize][0], m0[c1.id as usize][1]);
    let uv1 = at(m1, h.prev().crnr().id as usize);
    let uv1a = at(m1, h.next().crnr().id as usize);
    let dir0 = uv0a - uv0;
    let dir1 = uv1a - uv1;

    for i in 0..4 {
        let rot = get_quater_rot(i);
        if (rot * dir1 - dir0).norm() > 1e-5 { continue; }

        let rows: Vec<usize> = m1.as_ref().row_idx_of_col_raw(0).to_vec();
        for r in rows {
            let p = rot * (at(m1, r) - uv1) + uv0;
            *m1.get_mut(r, 0).expect("coeffRef") = p.re;
            *m1.get_mut(r, 1).expect("coeffRef") = p.im;
        }
        return;
    }
    metriko_fail!("failed to apply_transition");
}

/// The tutte parameterization of the cut mesh: every tquad embedded on its own (`embedding_tutte_for_tquad`), then
/// the tquads are laid out one after the other across the tedges, each rotated / translated to match the corner
/// values already fixed. Returns the uv of every crnr of `hm` (`nC` rows).
pub fn compute_tutte_parameterization(
    hm: &Hmesh,          // hmesh after tutte cutting
    tm: &Emesh<'_>,      // tmesh original
    seam: &[bool],       // seam adapted to tutte cutting
    data: &[HalfData],   //
) -> Vec<[f64; 2]> {
    let mut uv = vec![[0f64; 2]; hm.n_c as usize];

    // compute uv per tquad first...
    let empty = || from_triplets::<f64>(0, 0, &[]);
    #[cfg(not(feature = "rayon"))]
    let mut uv_tq: Vec<SprsD> = (0..tm.tquads.len())
        .map(|i| if tm.tquads[i].id == -1 { empty() } else { embedding_tutte_for_tquad(i as i32, data, hm, tm) })
        .collect();
    // the tquads are independent: the parallel loop gives the same result
    #[cfg(feature = "rayon")]
    let mut uv_tq: Vec<SprsD> = {
        use rayon::prelude::*;
        (0..tm.tquads.len()).into_par_iter()
            .map(|i| if tm.tquads[i].id == -1 { empty() } else { embedding_tutte_for_tquad(i as i32, data, hm, tm) })
            .collect()
    };

    let mut flag = vec![false; hm.n_f as usize];
    let mut stack: Vec<i32> = vec![];

    // unordered_map<Half, const HalfData*>: lookups only; emplace keeps the first entry of a half
    let mut data_by_half: HashMap<i32, usize> = HashMap::with_capacity(data.len());
    for (k, d) in data.iter().enumerate() { data_by_half.entry(d.half).or_insert(k); }

    { // 1: process the first tquad
        let h = data[0].half;
        let i = data[0].tqid;
        let mut b = vec![false; hm.n_h as usize];
        for d in data { if d.tqid == i { b[d.half as usize] = true; } }
        for nh in sequential_mapping(hm, &uv_tq[i as usize], h, &b, seam, &mut flag, &mut uv) { stack.push(nh); }
    }

    // 2: other tquads
    while let Some(top) = stack.pop() {
        let h = hm.half(top);
        if flag[h.face().id as usize] { continue; }

        let it = data_by_half.get(&h.id);
        metriko_check!(it.is_some(), "no half data for this half");

        let curr = &data[*it.unwrap()];
        apply_transition(h, &uv, &mut uv_tq[curr.tqid as usize]);

        let mut b = vec![false; hm.n_h as usize];
        for d in data { if d.tqid == curr.tqid { b[d.half as usize] = true; } }
        for nh in sequential_mapping(hm, &uv_tq[curr.tqid as usize], h.id, &b, seam, &mut flag, &mut uv) { stack.push(nh); }
    }
    uv
}
