//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tutte/tutte_cutting.h (see rust/PORTED_FROM).

use std::collections::{BTreeMap, BTreeSet, HashMap};
use crate::common::typedef::*;
use crate::common::stdcxx;
use crate::common::stdcomplex::std_div;
use crate::hmesh::hmesh::*;
use crate::hmesh::hmloc::*;
use crate::tmesh::emesh::*;
use crate::tutte::tutte::HalfData;

/// An element of `std::set<std::pair<double, int>>`, ordered as the pair (by `<` on the double, then the int).
#[derive(Clone, Copy, Debug)]
pub struct Split(pub f64, pub i32);

impl PartialEq for Split {
    fn eq(&self, o: &Self) -> bool { self.cmp(o) == std::cmp::Ordering::Equal }
}
impl Eq for Split {}
impl PartialOrd for Split {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> { Some(self.cmp(o)) }
}
impl Ord for Split {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        if self.0 < o.0 { return std::cmp::Ordering::Less; }
        if o.0 < self.0 { return std::cmp::Ordering::Greater; }
        self.1.cmp(&o.1)
    }
}

// split points on an original halfedge: (ratio, vid)
pub type EdgeSplits = BTreeSet<Split>;

#[derive(Clone, Copy, Debug)]
pub struct SegRec {
    pub i0: i32,
    pub i1: i32,
    pub d0: HalfData,
    pub d1: HalfData,
}

/// `std::minmax(a, b)`
#[inline]
fn minmax(a: i32, b: i32) -> (i32, i32) { if b < a { (b, a) } else { (a, b) } }

pub fn is_on_side(sides: &[Vec<i32>], vid: i32, sid: i32) -> bool {
    (vid as usize) < sides.len() && sides[vid as usize].contains(&sid)
}

pub fn common_side(sides: &[Vec<i32>], vid0: i32, vid1: i32) -> i32 {
    if vid0 as usize >= sides.len() { return -1; }
    for &l in &sides[vid0 as usize] {
        if is_on_side(sides, vid1, l) { return l; }
    }
    -1
}

/// the overload of `common_side` with three vertices
pub fn common_side3(sides: &[Vec<i32>], vid0: i32, vid1: i32, vid2: i32) -> i32 {
    if vid0 as usize >= sides.len() { return -1; }
    for &l in &sides[vid0 as usize] {
        if is_on_side(sides, vid1, l) &&
           is_on_side(sides, vid2, l) { return l; }
    }
    -1
}

// vertex chain along an original halfedge, tail -> head, split points included
pub fn chain_of(splits: &[EdgeSplits], h: Half<'_>) -> Vec<i32> {
    let mut chain = vec![h.tail().id];
    for s in splits[h.id as usize].iter().rev() { chain.push(s.1); } // descending r = tail -> head
    chain.push(h.head().id);
    chain
}

// triangulate one face against the segments crossing it. the resulting pieces may
// be CONCAVE (traced tedges bend inside a face), so the trace allows reflex turns
// and the pieces are ear-clipped rather than fanned.
pub fn face_cutting(
    f: Face<'_>,                // face to be cut
    sgs: &[[i32; 2]],           // segment endpoints (cut-mesh vertex indices)
    splits: &[EdgeSplits],      // split points per original halfedge
    sides: &[Vec<i32>],         // vertex -> tquad * side
    vpos: &[Row3d],             // vertex position
    tris: &mut Vec<[i32; 3]>,   // output triangles (appended, CCW wrt f)
) {
    // 1: collect the directed edges bounding the pieces:
    //    segments (both directions) + subdivided boundary halfedges (one direction)
    let mut halfs: Vec<[i32; 2]> = vec![];
    for s in sgs {
        halfs.push([s[0], s[1]]);
        halfs.push([s[1], s[0]]);
    }

    for h in f.adj_halfs() {
        let ch = chain_of(splits, h);
        let mut k = 0;
        while k + 1 < ch.len() { halfs.push([ch[k], ch[k + 1]]); k += 1; }
    }

    let cr = |u: Complex, v: Complex| (u.conj() * v).im;

    let penalty = 100.; // dominates min_angle in (-pi, pi]

    // 2: trace the pieces by the leftmost-turn rule (max signed CCW angle).
    //    reflex turns are legal, and every chain must close.
    while let Some(h) = halfs.pop() {
        let mut poly = vec![h];
        let sta = h[0];
        let mut cur = h[1];

        while cur != sta {
            let prev_v = poly.last().unwrap()[0];
            let curr_v = poly.last().unwrap()[1];
            let p_curr = f.to_local(vpos[curr_v as usize]);
            let d_prev = p_curr - f.to_local(vpos[prev_v as usize]);

            let mut best_it = halfs.len();
            let mut best = -10.; // signed turn angle in (-pi, pi], in the face plane like the ear clipping below
            for (it, e) in halfs.iter().enumerate() {
                if e[0] != curr_v { continue; }
                if e[1] == prev_v { continue; } // no immediate u-turn
                let ang = std_div(f.to_local(vpos[e[1] as usize]) - p_curr, d_prev).arg();
                if ang > best { best = ang; best_it = it; }
            }

            metriko_check!(best_it != halfs.len(), "open chain in face {}", f.id);

            poly.push(halfs[best_it]);
            cur = halfs[best_it][1];
            halfs.remove(best_it);
        }

        // 3: ear-clip the piece in the face plane (concave-capable, keeps CCW
        //    orientation so no flipped faces can be produced)
        let mut cyc: Vec<i32> = poly.iter().map(|p| p[0]).collect();

        // interior angles of a CCW triangle in the face plane (all positive);
        // negative for a CW triangle, so it also ranks invalid choices last
        let min_angle = |va: i32, vb: i32, vc: i32| -> f64 {
            let a = f.to_local(vpos[va as usize]);
            let b = f.to_local(vpos[vb as usize]);
            let c = f.to_local(vpos[vc as usize]);
            // std::min({x, y, z})
            let x = std_div(c - a, b - a).arg();
            let y = std_div(a - b, c - b).arg();
            let z = std_div(b - c, a - c).arg();
            let mut m = x;
            if y < m { m = y; }
            if z < m { m = z; }
            m
        };

        metriko_check!(cyc.len() >= 3, "polygon with {} vertices", cyc.len());

        while cyc.len() > 3 {
            let n = cyc.len();
            let mut best_k = n;
            let mut best_q = f64::MIN;
            for k in 0..n {
                let va = cyc[(k + n - 1) % n];
                let vb = cyc[k];
                let vc = cyc[(k + 1) % n];
                let a = f.to_local(vpos[va as usize]);
                let b = f.to_local(vpos[vb as usize]);
                let c = f.to_local(vpos[vc as usize]);
                if cr(b - a, c - b) <= EPS * (b - a).norm() * (c - b).norm() { continue; } // reflex or flat corner (sine below EPS): not an ear

                let mut empty = true; // no other cycle vertex inside the ear
                let mut m = 0;
                while m < n && empty {
                    if m == k || m == (k + 1) % n || m == (k + n - 1) % n { m += 1; continue; }
                    let q = f.to_local(vpos[cyc[m] as usize]);
                    empty = cr(b - a, q - a) < -EPS * (b - a).norm() * (q - a).norm()
                         || cr(c - b, q - b) < -EPS * (c - b).norm() * (q - b).norm()
                         || cr(a - c, q - c) < -EPS * (a - c).norm() * (q - c).norm();
                    m += 1;
                }
                if !empty { continue; }

                // quality of this clip: the worst interior angle it commits to.
                // a quad also fixes the remaining triangle, so score the pair —
                // this picks the better diagonal and avoids slivers whenever the
                // other choice is available
                let mut q = min_angle(va, vb, vc);
                if n == 4 { let q2 = min_angle(vc, cyc[(k + 2) % n], va); if q2 < q { q = q2; } }
                if common_side3(sides, va, vb, vc) >= 0                         { q -= penalty; }
                if n == 4 && common_side3(sides, vc, cyc[(k + 2) % n], va) >= 0 { q -= penalty; }
                if q > best_q { best_q = q; best_k = k; }
            }
            metriko_check!(best_k != n, "ear clipping failed in face {} (polygon of {})", f.id, cyc.len());
            tris.push([cyc[(best_k + n - 1) % n], cyc[best_k], cyc[(best_k + 1) % n]]);
            cyc.remove(best_k);
        }

        tris.push([cyc[0], cyc[1], cyc[2]]);
    }
}

/// directed edge -> tri on its left
fn index_tri(tri_of: &mut BTreeMap<(i32, i32), i32>, tris: &[[i32; 3]], i: usize) {
    for j in 0..3 { tri_of.insert((tris[i][j], tris[i][(j + 1) % 3]), i as i32); }
}

fn unindex_tri(tri_of: &mut BTreeMap<(i32, i32), i32>, tris: &[[i32; 3]], i: usize) {
    for j in 0..3 { tri_of.remove(&(tris[i][j], tris[i][(j + 1) % 3])); }
}

pub fn split_degenerate_faces(
    hm: &Hmesh,
    seam0: &[bool],
    sgs: &[SegRec],
    h_aux: &mut [EdgeSplits],
    sides: &[Vec<i32>],
    vpos: &mut Vec<Row3d>,
    tris: &mut Vec<[i32; 3]>,
) {
    let mut sg_set: BTreeSet<(i32, i32)> = BTreeSet::new();
    for s in sgs { sg_set.insert(minmax(s.i0, s.i1)); }

    // seam sub-edge -> (halfedge id, ratio of key.first, ratio of key.second):
    // lets a midpoint inserted on it be threaded back into h_aux so the
    // seam1 / matching1 propagation still finds every sub-halfedge
    let mut sm_info: BTreeMap<(i32, i32), (i32, f64, f64)> = BTreeMap::new();

    for e in hm.edges() {
        if !seam0[e.id as usize] { continue; }
        let h = e.half();
        let ch = chain_of(h_aux, h);
        let mut rs = vec![1.]; // r = 1 at tail, 0 at head along h
        for s in h_aux[h.id as usize].iter().rev() { rs.push(s.0); }
        rs.push(0.);
        let mut k = 0;
        while k + 1 < ch.len() {
            let key = minmax(ch[k], ch[k + 1]);
            let v = if key.0 == ch[k] { (h.id, rs[k], rs[k + 1]) } else { (h.id, rs[k + 1], rs[k]) };
            sm_info.insert(key, v);
            k += 1;
        }
    }

    let mut tri_of: BTreeMap<(i32, i32), i32> = BTreeMap::new(); // directed edge -> tri on its left, kept up to date by split_at
    for i in 0..tris.len() { index_tri(&mut tri_of, tris, i); }

    // split edge j of triangle i (and the neighbor across it) at the midpoint.
    // s: the side the degeneracy is on; rejected when the neighbor's opposite
    // vertex is also on it (the split would not resolve anything)
    let mut split_at = |tris: &mut Vec<[i32; 3]>, vpos: &mut Vec<Row3d>, tri_of: &mut BTreeMap<(i32, i32), i32>, h_aux: &mut [EdgeSplits], i: usize, j: usize, s: i32| -> bool {
        let a = tris[i][j];
        let b = tris[i][(j + 1) % 3];
        let c = tris[i][(j + 2) % 3];
        if sg_set.contains(&minmax(a, b)) { return false; }
        let k = *tri_of.get(&(b, a)).expect("map::at") as usize;
        let d = tris[k][0] + tris[k][1] + tris[k][2] - a - b;
        if is_on_side(sides, d, s) { return false; }
        let m = vpos.len() as i32;

        if let Some(&(hid, r_lo, r_hi)) = sm_info.get(&minmax(a, b)) {
            // the edge lies on a seam: thread the midpoint into h_aux
            let rm = (r_lo + r_hi) / 2.;
            h_aux[hid as usize].insert(Split(rm, m));
            h_aux[hm.half(hid).twin().id as usize].insert(Split(1. - rm, m));
            let (lo, hi) = minmax(a, b);
            sm_info.remove(&minmax(a, b));
            sm_info.insert((lo, m), (hid, r_lo, rm)); // m is always the largest id
            sm_info.insert((hi, m), (hid, r_hi, rm));
        }

        let mid = (vpos[a as usize] + vpos[b as usize]) / 2.;
        unindex_tri(tri_of, tris, i);
        unindex_tri(tri_of, tris, k);
        tris[i] = [a, m, c];
        tris[k] = [b, m, d];
        tris.push([m, b, c]);
        tris.push([m, a, d]);
        let n = tris.len();
        for t in [i, k, n - 2, n - 1] { index_tri(tri_of, tris, t); }
        vpos.push(mid);
        true
    };

    let mut again = true;
    while std::mem::replace(&mut again, false) {
        // fully collinear triangles
        let mut i = 0;
        while i < tris.len() && !again {
            let l = common_side3(sides, tris[i][0], tris[i][1], tris[i][2]);
            if l >= 0 {
                let mut j = 0;
                while j < 3 && !again { again = split_at(tris, vpos, &mut tri_of, h_aux, i, j, l); j += 1; }
            }
            i += 1;
        }

        // interior chords on a side
        let mut i = 0;
        while i < tris.len() && !again {
            let mut j = 0;
            while j < 3 && !again {
                let a = tris[i][j];
                let b = tris[i][(j + 1) % 3];
                let c = tris[i][(j + 2) % 3];
                let s = common_side(sides, a, b);
                if s >= 0 && !is_on_side(sides, c, s) { again = split_at(tris, vpos, &mut tri_of, h_aux, i, j, s); }
                j += 1;
            }
            i += 1;
        }
    }
}

/// `vid_of(nid)` of `compute_embedding_cut_hmesh`: nid -> cut-mesh vertex index. tnode identity replaces
/// epsilon-based position matching. creates the vertex and registers the edge split exactly once.
fn vid_of(
    nid: i32,
    hm: &Hmesh,
    tm: &Emesh<'_>,
    vid_of_nid: &mut HashMap<i32, i32>,
    vpos: &mut Vec<Row3d>,
    auxs: &mut [EdgeSplits],
) -> i32 {
    if let Some(&v) = vid_of_nid.get(&nid) { return v; }
    {
        let l = tm.tnodes[nid as usize];
        let v = match l {
            HmLoc::OnV { id } => id,
            _ => {
                vpos.push(get_ptloc_pos(hm, &l));
                let i = vpos.len() as i32 - 1;
                if let Some((h0, r)) = try_get_half_ratio(hm, &l) {
                    auxs[h0.id as usize].insert(Split(r, i));
                    auxs[h0.twin().id as usize].insert(Split(1. - r, i));
                }
                i
            }
        };
        vid_of_nid.insert(nid, v);
        v
    }
}

/// The mesh `hm` cut along the traced tedges of `tm` (every tedge becomes a chain of halfedges), with the seam,
/// matching and singularities carried over, and the halfedges lying on tedges annotated (`data`, sorted by
/// `(tqid, thid, order)`, both directions of a cut halfedge adjacent). Returns
/// `(hm_cut, seam1, matching1, singular1, data)`; `HalfData::half` are half ids of `hm_cut`.
pub fn compute_embedding_cut_hmesh(
    hm: &Hmesh,
    tm: &Emesh<'_>,
    seam0: &[bool],
    matching0: &[i32],
    singular0: &[i32],
) -> (Hmesh, Vec<bool>, Vec<i32>, Vec<i32>, Vec<HalfData>) {
    let mut cuts: BTreeMap<i32, Vec<[i32; 2]>> = BTreeMap::new(); // face id -> segment endpoints

    let mut sgms: Vec<SegRec> = vec![]; // flat annotation list; matched to cut halfedges at the end
    let mut auxs: Vec<EdgeSplits> = vec![EdgeSplits::new(); hm.n_h as usize];

    let mut vpos: Vec<Row3d> = hm.pos.clone();

    let mut vid_of_nid: HashMap<i32, i32> = HashMap::new();
    let mut sides: Vec<Vec<i32>> = vec![]; // cut vertex -> tquad * side

    for th in tm.cano_thalfs() {
        let nids = &th.tedge(tm).nids;
        let nsgs = nids.len() as i32 - 1;
        let tw = th.twin(tm);
        let sa = th.tqid * 4 + th.side(tm);
        let sb = tw.tqid * 4 + tw.side(tm);

        let mut total = 0.;
        let mut len = vec![0.; nsgs as usize];
        for k in 0..nsgs as usize {
            len[k] = (get_ptloc_pos(hm, &tm.tnodes[nids[k + 1] as usize]) - get_ptloc_pos(hm, &tm.tnodes[nids[k] as usize])).length();
            total += len[k];
        }
        metriko_check!(total >= EPS, "zero-length tedge");

        let mut sum = 0.;
        for k in 0..nsgs as usize {
            let fr = &tm.tnodes[nids[k] as usize];
            let to = &tm.tnodes[nids[k + 1] as usize];
            let v0 = sum / total;
            sum += len[k];
            let v1 = sum / total;
            let i0 = vid_of(nids[k], hm, tm, &mut vid_of_nid, &mut vpos, &mut auxs);
            let i1 = vid_of(nids[k + 1], hm, tm, &mut vid_of_nid, &mut vpos, &mut auxs);
            for i in [i0, i1] {
                if i as usize >= sides.len() { sides.resize(vpos.len(), vec![]); }
                let ss = &mut sides[i as usize];
                if !ss.contains(&sa) { ss.push(sa); }
                if !ss.contains(&sb) { ss.push(sb); }
            }
            let eo = try_get_edge(hm, fr, to);
            let fo = try_get_face(hm, fr, to);
            let t = nsgs - k as i32 - 1;
            if eo.is_none() { cuts.entry(fo.expect("bad optional access").id).or_default().push([i0, i1]); }
            sgms.push(SegRec {
                i0,
                i1,
                d0: HalfData { v0,          v1,          thid: th.id,   tqid: th.tqid, order: k as i32, ..Default::default() },
                d1: HalfData { v0: 1. - v1, v1: 1. - v0, thid: th.twid, tqid: tw.tqid, order: t,        ..Default::default() },
            });
        }
    }

    let mut tris: Vec<[i32; 3]> = vec![];
    for f in hm.faces() {
        let it = cuts.get(&f.id);
        let cut = it.is_some();
        let split = f.halfs().iter().any(|h| !auxs[h.id as usize].is_empty());
        if !cut && !split { let [a, b, c] = f.verts(); tris.push([a.id, b.id, c.id]); }
        else { face_cutting(f, it.map(|v| v.as_slice()).unwrap_or(&[]), &auxs, &sides, &vpos, &mut tris); }
    }

    split_degenerate_faces(hm, seam0, &sgms, &mut auxs, &sides, &mut vpos, &mut tris);

    let hm_cut = Hmesh::new(&vpos, &tris, false);
    let mut seam1 = vec![false; hm_cut.n_e as usize];
    let mut matching1 = vec![0i32; hm_cut.n_e as usize];
    let mut singular1 = vec![0i32; hm_cut.n_v as usize];
    singular1[..hm.n_v as usize].copy_from_slice(&singular0[..hm.n_v as usize]);

    let mut half_by_verts: BTreeMap<(i32, i32), i32> = BTreeMap::new();
    for h in hm_cut.halfs() {
        half_by_verts.entry((h.tail().id, h.head().id)).or_insert(h.id);
    }

    for e in hm.edges() {
        if !seam0[e.id as usize] { continue; }
        let ch = chain_of(&auxs, e.half());
        let mut k = 0;
        while k + 1 < ch.len() {
            let h = hm_cut.half(*half_by_verts.get(&(ch[k], ch[k + 1])).expect("map::at")); // same direction as h
            seam1[h.edge().id as usize] = true;
            matching1[h.edge().id as usize] = if h.is_canonical() { matching0[e.id as usize] } else { -matching0[e.id as usize] };
            k += 1;
        }
    }

    // annotate the cut halfedges: both directions pushed adjacently, twin = index ^ 1
    let mut data: Vec<HalfData> = vec![];
    for sg in sgms.iter_mut() {
        let it = half_by_verts.get(&(sg.i0, sg.i1));
        metriko_check!(it.is_some(), "no half between vertices");
        let h = hm_cut.half(*it.unwrap());
        sg.d0.half = h.id;
        sg.d1.half = h.twin().id;
        data.push(sg.d0);
        data.push(sg.d1);
    }

    stdcxx::sort_by(&mut data, |a, b| a.lt(b));
    (hm_cut, seam1, matching1, singular1, data)
}
