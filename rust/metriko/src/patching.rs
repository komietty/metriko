//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/patching.h (see rust/PORTED_FROM).

// quad patching: build the quad mesh from the extracted q-faces and replay the
// t-mesh on it, so that every quad knows which tquad (patch) it belongs to.

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use crate::common::typedef::*;
use crate::common::utilities::circular_prev;
use crate::common::stdcxx::LibcxxIntMap;
use crate::hmesh::hmesh::Hmesh;
use crate::hmesh::hmloc::HmLoc;
use crate::igl::aabb::Aabb;
use crate::igl::per_face_normals::per_face_normals;
use crate::igl::remove_duplicate_vertices::remove_duplicate_vertices;
use crate::igl::umeyama::{umeyama, Mat34};
use crate::qex::common::*;
use crate::tmesh::emesh::Emesh;

//------------------------------------------------------------------------------
// quad mesh from the q-faces
//------------------------------------------------------------------------------

/// quad soup of the extracted q-faces: 4 corners per quad, corner (i, j) is qfaces[i].qhalfs[j].port1()
/// (`q_ports` / `q_edges`: the vectors the C++ q-halfs reference)
pub fn quad_soup(qports: &[Qport], qedges: &[Qedge], qfaces: &[Qface]) -> (Vec<Row3d>, Vec<[i32; 4]>) {
    let l = qfaces.len();
    let mut pos = vec![Row3d::ZERO; l * 4];
    let mut idx = vec![[0i32; 4]; l];
    for i in 0..l {
        for j in 0..4 {
            pos[i * 4 + j] = qports[qfaces[i].qhalfs[j].port1(qedges)].pos;
            idx[i][j] = (i * 4 + j) as i32;
        }
    }
    (pos, idx)
}

/// weld the per-quad duplicated corners into a connected quad mesh
pub fn weld_quad_soup(pos: &[Row3d], idx: &[[i32; 4]]) -> (Vec<Row3d>, Vec<[i32; 4]>) {
    // (pos.colwise().maxCoeff() - pos.colwise().minCoeff()).norm()
    let mut mx = pos[0];
    let mut mn = pos[0];
    for p in &pos[1..] {
        if p.x > mx.x { mx.x = p.x; }
        if p.y > mx.y { mx.y = p.y; }
        if p.z > mx.z { mx.z = p.z; }
        if p.x < mn.x { mn.x = p.x; }
        if p.y < mn.y { mn.y = p.y; }
        if p.z < mn.z { mn.z = p.z; }
    }
    let eps = 1e-7 * (mx - mn).length();
    let (pos_w, _svi, svj) = remove_duplicate_vertices(pos, eps);
    let idx_w = idx.iter().map(|q| q.map(|i| svj[i as usize])).collect();
    (pos_w, idx_w)
}

// conformal relaxation of a quad mesh on the input surface: shape-up (bouaziz et al. 2012) with the set of squares as
// the shape constraint. local step: every quad is projected onto its closest square, a similarity of the unit square
// fitted to its four corners (umeyama 1991). global step: every vertex takes the least-squares compromise of the
// corners its quads ask for, which for this constraint set alone is their mean. the surface closeness is enforced as
// a projection onto the input surface after every global step (the limit of an infinite weight on that constraint),
// so the vertices never leave it. a square has right angles and equal sides, and its size is free per quad: the
// iteration drives every corner towards 90 degrees while the quad size varies smoothly over the mesh.
// not part of shape-up: a step must not fold a quad. a quad convex with respect to the surface normal (summed over the
// faces closest to its corners) before a step keeps its corners where they were if the step would make it non-convex;
// reverting corners can fold a neighbour in turn, so this repeats until nothing changes. quads that came out of the
// quad extraction folded move freely, and are guarded once they come out convex
/// `conformalize_on_surface(pos, idx, V, F, max_iter = 200, tol = 1e-6)`
pub fn conformalize_on_surface(pos: &mut Vec<Row3d>, idx: &[[i32; 4]], v: &[Row3d], f: &[[i32; 3]], max_iter: i32, tol: f64) {
    // the unit square, corners in the cyclic order of a quad
    let square: Mat34 = [
        [-1.,  1., 1., -1.],
        [-1., -1., 1.,  1.],
        [ 0.,  0., 0.,  0.],
    ];
    let mut n_quads = vec![0i32; pos.len()];
    for q in idx { for j in 0..4 { n_quads[q[j] as usize] += 1; } }

    let n = per_face_normals(v, f);
    let tree = Aabb::init(v, f);
    let (_sqr_d, mut ii, _c) = tree.squared_distance(v, f, pos);

    let convex = |p: &[Row3d], face_of: &[i32], q: usize| -> bool {
        let mut nn = Row3d::ZERO;
        for j in 0..4 { nn += n[face_of[idx[q][j] as usize] as usize]; }
        for j in 0..4 {
            let c = p[idx[q][j] as usize];
            let a = p[idx[q][(j + 3) % 4] as usize];
            let b = p[idx[q][(j + 1) % 4] as usize];
            if (b - c).cross(a - c).dot(nn) <= 0. { return false; }
        }
        true
    };
    let mut guarded: Vec<bool> = (0..idx.len()).map(|q| convex(pos, &ii, q)).collect();

    // fit[i] = (T.topLeftCorner<3, 3>() * square).colwise() + T.topRightCorner<3, 1>(), corners as columns
    let fit_of = |pos: &[Row3d], q: &[i32; 4]| -> [Row3d; 4] {
        let mut p: Mat34 = [[0.; 4]; 3];
        for j in 0..4 { let c = pos[q[j] as usize]; p[0][j] = c.x; p[1][j] = c.y; p[2][j] = c.z; }
        let t = umeyama(&square, &p, true);
        let mut fit = [Row3d::ZERO; 4];
        for j in 0..4 {
            // the products with the square's 0 / +-1 entries are exact, so the order of the lazy product's sums
            // (fused for rows 0, 1, halved for row 2, see igl::umeyama) gives these values
            let e = |i: usize| (t[i][0] * square[0][j] + t[i][1] * square[1][j]) + t[i][2] * square[2][j];
            fit[j] = Row3d::new(e(0) + t[0][3], e(1) + t[1][3], e(2) + t[2][3]);
        }
        fit
    };

    for _it in 0..max_iter {
        // #pragma omp parallel for: the quads are independent
        #[cfg(feature = "rayon")]
        let fit: Vec<[Row3d; 4]> = { use rayon::prelude::*; idx.par_iter().map(|q| fit_of(pos, q)).collect() };
        #[cfg(not(feature = "rayon"))]
        let fit: Vec<[Row3d; 4]> = idx.iter().map(|q| fit_of(pos, q)).collect();

        let mut next = vec![Row3d::ZERO; pos.len()];
        for i in 0..idx.len() { for j in 0..4 { next[idx[i][j] as usize] += fit[i][j]; } }
        for vv in 0..pos.len() { next[vv] = if n_quads[vv] > 0 { next[vv] / n_quads[vv] as f64 } else { pos[vv] }; }

        let (_sqr_d, mut i_next, mut c) = tree.squared_distance(v, f, &next);
        let mut changed = true;
        while changed {
            changed = false;
            for q in 0..idx.len() {
                if !guarded[q] || convex(&c, &i_next, q) { continue; }
                for j in 0..4 {
                    let vv = idx[q][j] as usize;
                    if c[vv] == pos[vv] { continue; }
                    c[vv] = pos[vv];
                    i_next[vv] = ii[vv];
                    changed = true;
                }
            }
        }
        for q in 0..idx.len() { if !guarded[q] { guarded[q] = convex(&c, &i_next, q); } }

        // (C - pos).rowwise().norm().maxCoeff(), in grid units: the pipeline runs on V / grid_unit
        let mut mv = (c[0] - pos[0]).length();
        for vv in 1..pos.len() { let d = (c[vv] - pos[vv]).length(); if d > mv { mv = d; } }
        *pos = c;
        ii = i_next;
        if mv < tol { break; }
    }
}

/// welded (and optionally smoothed) quad mesh of the q-faces. face (i, j) still corresponds to qfaces[i].qhalfs[j].port1()
/// `extract_quad_mesh(hm, qfaces, refine = true)`; the C++ q-faces reach their ports by reference, here through
/// `q_ports` / `q_edges`.
pub fn extract_quad_mesh(hm: &Hmesh, qports: &[Qport], qedges: &[Qedge], qfaces: &[Qface], refine: bool) -> (Vec<Row3d>, Vec<[i32; 4]>) {
    let (pos, idx) = quad_soup(qports, qedges, qfaces);
    let (mut pos_w, idx_w) = weld_quad_soup(&pos, &idx);
    if refine { conformalize_on_surface(&mut pos_w, &idx_w, &hm.pos, &hm.idx, 200, 1e-6); }
    (pos_w, idx_w)
}

//------------------------------------------------------------------------------
// patch labeling: replay the t-mesh on the quad mesh
//------------------------------------------------------------------------------

/// the quad mesh as directed edges
pub struct QuadGraph<'a> {
    idx: &'a [[i32; 4]],
    corner: BTreeMap<(i32, i32), (i32, i32)>, // directed edge a -> b: (quad on its left, corner a)
    nbrs: Vec<Vec<i32>>,                      // quad vertex -> neighbours
    qv_of_vid: HashMap<i32, i32>,             // mesh vertex -> quad vertex (q-vertices on a vertex)
}

impl<'a> QuadGraph<'a> {
    pub fn new(qports: &[Qport], qedges: &[Qedge], qfaces: &[Qface], idx: &'a [[i32; 4]]) -> Self {
        let n = if !idx.is_empty() { idx.iter().flatten().copied().max().unwrap() + 1 } else { 0 };
        let mut g = QuadGraph { idx, corner: BTreeMap::new(), nbrs: vec![vec![]; n as usize], qv_of_vid: HashMap::new() };
        for i in 0..idx.len() {
            for j in 0..4 {
                let a = idx[i][j];
                let b = idx[i][(j + 1) % 4];
                g.corner.insert((a, b), (i as i32, j as i32));
                g.nbrs[a as usize].push(b);
                let v = qports[qfaces[i].qhalfs[j].port1(qedges)].vid;
                if v >= 0 { g.qv_of_vid.insert(v, a); }
            }
        }
        g
    }
    /// quad on the left of a -> b
    fn quad(&self, a: i32, b: i32) -> i32 { self.corner.get(&(a, b)).expect("map::at").0 }
    fn valence(&self, v: i32) -> i32 { self.nbrs[v as usize].len() as i32 }
    /// neighbour of v after w, CCW
    fn rot(&self, v: i32, w: i32) -> i32 { let (q, j) = *self.corner.get(&(v, w)).expect("map::at"); self.idx[q as usize][((j + 3) % 4) as usize] }
    /// arrived t -> v: keep going
    fn straight(&self, t: i32, v: i32) -> i32 { self.rot(v, self.rot(v, t)) }
    /// neighbours of v, CCW from `first`
    fn ring(&self, v: i32, first: i32) -> Vec<i32> {
        let mut r = vec![first];
        while (r.len() as i32) < self.valence(v) { let b = *r.last().unwrap(); r.push(self.rot(v, b)); }
        r
    }
}

/// canonical thalf, quantized length, tquads along nids
#[derive(Clone, Copy, Debug)]
struct TvEdge { cano: i32, steps: i32, left: i32, right: i32 }
impl Default for TvEdge { fn default() -> Self { TvEdge { cano: -1, steps: 0, left: -1, right: -1 } } }

/// a tedge leaving a node, forward along its nids or not, and its slot offset around the node
#[derive(Clone, Copy, Debug)]
struct Branch { teid: i32, fwd: bool, slot: i32 }

/// the t-mesh as seen from the quad mesh, built once and never changed
struct TmeshView<'a, 'b> {
    tm: &'a Emesh<'b>,
    singular: &'a [i32],
    g: &'a QuadGraph<'a>,
    te: Vec<TvEdge>,
    // node -> (teid, leaves the node forward). `label_quad_patches` iterates this unordered_map, so it keeps libc++'s
    // order; `branch_count` answers `branches.at(nid).size()`
    branches: LibcxxIntMap<Vec<(i32, bool)>>,
    branch_count: HashMap<i32, usize>,
}

impl<'a, 'b> TmeshView<'a, 'b> {
    fn new(tm: &'a Emesh<'b>, singular: &'a [i32], g: &'a QuadGraph<'a>) -> Self {
        let mut te = vec![TvEdge::default(); tm.tedges.len()];
        for th in tm.cano_thalfs() { te[th.teid as usize] = TvEdge { cano: th.id, steps: th.x.round() as i32, left: th.tqid, right: th.twin(tm).tqid }; }
        let mut branches: LibcxxIntMap<Vec<(i32, bool)>> = LibcxxIntMap::new();
        for e in tm.live_tedges() {
            branches.entry(*e.nids.first().unwrap()).push((e.id, true));
            branches.entry(*e.nids.last().unwrap()).push((e.id, false));
        }
        let branch_count = branches.iter().map(|(k, v)| (k, v.len())).collect();
        TmeshView { tm, singular, g, te, branches, branch_count }
    }
    fn far_node(&self, teid: i32, fwd: bool) -> i32 { let nids = &self.tm.tedges[teid as usize].nids; if fwd { *nids.last().unwrap() } else { *nids.first().unwrap() } }
    /// the singular vertex a node sits on, or None
    fn singular_vert(&self, nid: i32) -> Option<i32> {
        match self.tm.tnodes[nid as usize] {
            HmLoc::OnV { id } if self.singular[id as usize] != 0 => Some(id),
            _ => None,
        }
    }
    /// quad vertex of a singular node, -1 otherwise
    fn singular_qv(&self, nid: i32) -> i32 {
        match self.singular_vert(nid) {
            Some(v) => self.g.qv_of_vid.get(&v).copied().unwrap_or(-1),
            None => -1,
        }
    }
    // all tedges at a node, CCW from (t, fwd), each with its slot offset in quad edges. rotating CCW around the node
    // from an outgoing tedge sweeps through the tquad P on its left and reaches P's previous boundary thalf, which
    // ends at the node. this follows P's boundary order, so a tquad touching itself at the node is handled; the slot
    // advances by 2 when the node lies inside a side of P
    fn chain(&self, nid: i32, mut t: i32, mut fwd: bool) -> Vec<Branch> {
        let tm = self.tm;
        let mut ch = vec![Branch { teid: t, fwd, slot: 0 }];
        let mut slot = 0;
        for _k in 1..*self.branch_count.get(&nid).expect("unordered_map::at") {
            let c = self.te[t as usize].cano;
            let th = if fwd { &tm.thalfs[c as usize] } else { tm.thalfs[c as usize].twin(tm) }; // P on its left, pointing away from the node
            let tq = th.tquad(tm);
            let cur = tq.data.iter().position(|d| d.thid == th.id).expect("find");
            let prv = circular_prev(tq.data.len(), cur);
            let th2 = &tm.thalfs[tq.data[prv].thid as usize]; // ends at the node
            slot += if tq.data[prv].side == tq.data[cur].side { 2 } else { 1 };
            t = th2.teid;
            fwd = !th2.cano; // outgoing from the node: backwards along nids if cano
            ch.push(Branch { teid: t, fwd, slot });
        }
        ch
    }
}

/// walk teid starting with the quad edge from -> to
#[derive(Clone, Copy, Debug)]
struct Job { teid: i32, fwd: bool, from: i32, to: i32 }

// walk every tedge on the quad mesh, labeling the quads on its two sides. the state of the walk only: a replay is
// copied to try the rotations at a singular
#[derive(Clone)]
struct Replay<'v, 'a, 'b> {
    tv: &'v TmeshView<'a, 'b>,
    tqid_of_quad: Vec<i32>,     // the labels so far, -1 when none
    track: BTreeSet<(i32, i32)>, // quad edges (min, max) lying on a tedge
    node_qv: HashMap<i32, i32>, // node -> quad vertex
    done: Vec<bool>,            // per tedge
    jobs: VecDeque<Job>,
    ok: bool,
}

impl<'v, 'a, 'b> Replay<'v, 'a, 'b> {
    fn new(tv: &'v TmeshView<'a, 'b>) -> Self {
        Replay {
            tv,
            tqid_of_quad: vec![-1; tv.g.idx.len()],
            track: BTreeSet::new(),
            node_qv: HashMap::new(),
            done: vec![false; tv.tm.tedges.len()],
            jobs: VecDeque::new(),
            ok: true,
        }
    }

    // anchor a node at a quad vertex and queue the tedges leaving it, the first at slot `rot`
    fn anchor(&mut self, nid: i32, qv: i32, branches: &[Branch], ring: &[i32], rot: i32) {
        self.node_qv.insert(nid, qv);
        for &Branch { teid: t, fwd, slot } in branches {
            if !self.done[t as usize] { self.jobs.push_back(Job { teid: t, fwd, from: qv, to: ring[((slot + rot) as usize) % ring.len()] }); }
        }
    }

    // walk a tedge straight along the quad edges; false at the first contradiction
    fn walk(&mut self, jb: &Job) -> bool {
        let tv = self.tv;
        let label = |tq: &mut Vec<i32>, q: i32, tqid: i32| -> bool {
            let l = &mut tq[q as usize];
            if *l >= 0 && *l != tqid { return false; }
            *l = tqid;
            true
        };
        let e = tv.te[jb.teid as usize];
        let lq = if jb.fwd { e.left } else { e.right };
        let rq = if jb.fwd { e.right } else { e.left };
        let mut a = jb.from;
        let mut b = jb.to;
        for k in 0..e.steps {
            if !label(&mut self.tqid_of_quad, tv.g.quad(a, b), lq) || !label(&mut self.tqid_of_quad, tv.g.quad(b, a), rq) { return false; }
            self.track.insert((a.min(b), a.max(b)));
            if k + 1 == e.steps { break; }
            if tv.g.valence(b) != 4 { return false; } // irregular vertex before the far node
            let nb = tv.g.straight(a, b);
            a = b;
            b = nb;
        }
        let nid = tv.far_node(jb.teid, jb.fwd);
        let s = tv.singular_qv(nid);
        if s >= 0 && s != b { return false; }                              // landed away from the singular
        if let Some(&q) = self.node_qv.get(&nid) { return q == b; }        // anchored elsewhere
        self.anchor(nid, b, &tv.chain(nid, jb.teid, !jb.fwd), &tv.g.ring(b, a), 0); // slots CCW from the reverse of arrival
        true
    }

    fn run(&mut self) {
        while self.ok && !self.jobs.is_empty() {
            let jb = self.jobs.pop_front().unwrap();
            if self.done[jb.teid as usize] { continue; }
            self.done[jb.teid as usize] = true;
            self.ok = self.walk(&jb);
        }
    }
}

// replay the t-mesh on the quad mesh and label every quad by its tquad (patch): the tquad of quad i is the i-th entry,
// -1 when none was found, with no geometry involved.
//
// singular vertices are q-vertices, so their quad vertex is exact. every tedge walks exactly x quad edges straight
// ahead, and its landing vertex IS its far node, which anchors the tedges there in CCW order. the only unknown is
// the rotation of the chain at the first singular of each component: every rotation is replayed and the one
// without contradictions is kept. the labels are then flooded into the patch interiors without crossing a track
/// `label_quad_patches(tm, singular, qfaces, qidx)`; `qidx`: welded quad mesh faces, corner (i, j) <->
/// qfaces[i].qhalfs[j].port1() (reached through `q_ports` / `q_edges`)
pub fn label_quad_patches(tm: &Emesh<'_>, singular: &[i32], qports: &[Qport], qedges: &[Qedge], qfaces: &[Qface], qidx: &[[i32; 4]]) -> Vec<i32> {
    let g = QuadGraph::new(qports, qedges, qfaces, qidx);
    let tv = TmeshView::new(tm, singular, &g);
    let mut rp = Replay::new(&tv);
    let mut no_rotation = 0; // singular nodes with no consistent rotation
    let mut several = 0;     // with more than one
    for (nid, brs) in tv.branches.iter() {
        let s = tv.singular_qv(nid);
        if s < 0 || rp.node_qv.contains_key(&nid) { continue; }
        let ring = g.ring(s, *g.nbrs[s as usize].first().unwrap());
        let chain = tv.chain(nid, brs[0].0, brs[0].1);
        let mut good: Vec<Replay> = vec![];
        for rot in 0..ring.len() {
            let mut t = rp.clone();
            t.anchor(nid, s, &chain, &ring, rot as i32);
            t.run();
            if t.ok { good.push(t); }
        }
        if good.is_empty() { no_rotation += 1; continue; }
        several += (good.len() > 1) as i32;
        rp = good.swap_remove(0);
    }
    let mut unreached = 0; // tedges never walked
    for e in tm.live_tedges() { if tv.te[e.id as usize].steps > 0 && !rp.done[e.id as usize] { unreached += 1; } }

    let mut tqid = std::mem::take(&mut rp.tqid_of_quad);
    let mut que: VecDeque<i32> = VecDeque::new(); // flood
    for i in 0..tqid.len() { if tqid[i] >= 0 { que.push_back(i as i32); } }
    while let Some(q) = que.pop_front() {
        for j in 0..4 {
            let a = qidx[q as usize][j];
            let b = qidx[q as usize][(j + 1) % 4];
            if rp.track.contains(&(a.min(b), a.max(b))) { continue; }
            let nb = g.quad(b, a);
            if tqid[nb as usize] < 0 { tqid[nb as usize] = tqid[q as usize]; que.push_back(nb); }
        }
    }
    let unlabeled = tqid.iter().filter(|&&x| x < 0).count();
    if no_rotation != 0 || several != 0 || unreached != 0 || unlabeled != 0 {
        eprintln!("[quad patch] singulars with no consistent rotation {} | with several {} | unreached tedges {} | unlabeled quads {}", no_rotation, several, unreached, unlabeled);
    }
    tqid
}
