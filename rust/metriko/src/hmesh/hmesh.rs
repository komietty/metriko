//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/hmesh/hmesh.h (see rust/PORTED_FROM).

use std::hash::{Hash, Hasher};
use crate::common::typedef::*;
use crate::common::utilities::*;
use crate::igl::edge_topology::edge_topology;

// The C++ elements are `{int id; Hmesh* m;}` compared by id only. Here they are Copy handles borrowing the mesh.
macro_rules! elem {
    ($name:ident) => {
        #[derive(Clone, Copy)]
        pub struct $name<'a> {
            pub id: i32,
            pub m: &'a Hmesh,
        }
        impl PartialEq for $name<'_> {
            fn eq(&self, rhs: &Self) -> bool { self.id == rhs.id }
        }
        impl Eq for $name<'_> {}
        impl Hash for $name<'_> {
            fn hash<H: Hasher>(&self, state: &mut H) { self.id.hash(state); }
        }
        impl std::fmt::Debug for $name<'_> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "{}({})", stringify!($name), self.id)
            }
        }
    };
}
elem!(Half);
elem!(Vert);
elem!(Edge);
elem!(Face);
elem!(Crnr);
elem!(Loop);

pub struct Hmesh {
    pub n_v: i32,
    pub n_f: i32,
    pub n_e: i32,
    pub n_h: i32,
    pub n_l: i32,
    pub n_c: i32,
    pub n_eular_chars: i32,

    pub vert2half: Vec<i32>,
    pub edge2half: Vec<i32>,
    pub face2half: Vec<i32>,
    pub loop2half: Vec<i32>,
    pub crnr2half: Vec<i32>,
    pub next: Vec<i32>,
    pub prev: Vec<i32>,
    pub twin: Vec<i32>,
    pub tail: Vec<i32>,
    pub head: Vec<i32>,
    pub edge: Vec<i32>,
    pub face: Vec<i32>,
    pub crnr: Vec<i32>,
    pub is_bv: Vec<bool>,
    pub pos: Vec<Row3d>,
    pub idx: Vec<[i32; 3]>,
    pub cfn: Vec<Complex>,
    pub edge2vert: Vec<[i32; 2]>,
    pub edge2face: Vec<[i32; 2]>,
    pub face2edge: Vec<[i32; 3]>,
    pub vert_basis_x: Vec<Row3d>,
    pub vert_basis_y: Vec<Row3d>,
    pub vert_normal: Vec<Row3d>,
    pub face_basis_x: Vec<Row3d>,
    pub face_basis_y: Vec<Row3d>,
    pub face_normal: Vec<Row3d>,
    pub face_area: Vec<f64>,
    pub bary_dual_area: Vec<f64>,
    pub circ_dual_area: Vec<f64>,
    pub half_cotan: Vec<f64>,
    pub edge_cotan: Vec<f64>,
    pub bary_center: Vec<Row3d>,
    pub angle_defect: Vec<f64>,
    pub dihedral_arg: Vec<f64>,
    pub he_arg_on_vert: Vec<f64>,
    pub he_arg_on_face: Vec<f64>,
}

/// Adjacency walk: `AdjVH` (around a vertex), `AdjFH` (around a face), `AdjLH` (along a boundary loop) and `AdjVC`
/// (the crnrs of the halfs around a vertex). The walk yields its start first and stops when it comes back to it,
/// exactly as the C++ begin/end pair.
#[derive(Clone, Copy, PartialEq, Eq)]
enum AdjKind { VH, FH }

#[derive(Clone, Copy)]
pub struct AdjIter<'a> {
    m: &'a Hmesh,
    start: i32,
    h: i32,
    bgn: bool,
    ccw: bool,
    kind: AdjKind,
}

impl<'a> Iterator for AdjIter<'a> {
    type Item = Half<'a>;
    fn next(&mut self) -> Option<Half<'a>> {
        if self.bgn && self.h == self.start { return None; }
        let h = Half { id: self.h, m: self.m };
        self.h = match (self.kind, self.ccw) {
            (AdjKind::VH, true)  => h.prev().twin().id,
            (AdjKind::VH, false) => h.twin().next().id,
            (AdjKind::FH, true)  => h.next().id,
            (AdjKind::FH, false) => h.prev().id,
        };
        self.bgn = true;
        Some(h)
    }
}

impl<'a> AdjIter<'a> {
    fn new(m: &'a Hmesh, hid: i32, ccw: bool, kind: AdjKind) -> Self {
        AdjIter { m, start: hid, h: hid, bgn: false, ccw, kind }
    }
}

impl<'a> Half<'a> {
    #[inline] pub fn next(self) -> Half<'a> { Half { id: self.m.next[self.id as usize], m: self.m } }
    #[inline] pub fn prev(self) -> Half<'a> { Half { id: self.m.prev[self.id as usize], m: self.m } }
    #[inline] pub fn twin(self) -> Half<'a> { Half { id: self.m.twin[self.id as usize], m: self.m } }
    #[inline] pub fn tail(self) -> Vert<'a> { Vert { id: self.m.tail[self.id as usize], m: self.m } }
    #[inline] pub fn head(self) -> Vert<'a> { Vert { id: self.m.head[self.id as usize], m: self.m } }
    #[inline] pub fn edge(self) -> Edge<'a> { Edge { id: self.m.edge[self.id as usize], m: self.m } }
    #[inline] pub fn face(self) -> Face<'a> { Face { id: self.m.face[self.id as usize], m: self.m } }
    #[inline] pub fn crnr(self) -> Crnr<'a> { Crnr { id: self.m.crnr[self.id as usize], m: self.m } }
    #[inline] pub fn cr_t(self) -> Crnr<'a> { self.next().crnr() }
    #[inline] pub fn cr_h(self) -> Crnr<'a> { self.prev().crnr() }
    #[inline] pub fn is_canonical(self) -> bool { self.edge().half().id == self.id }
    #[inline] pub fn is_boundary(self) -> bool { self.face().id == -1 }
    #[inline] pub fn len(self) -> f64 { self.vec().length() }
    #[inline] pub fn cot(self) -> f64 { self.m.half_cotan[self.id as usize] }
    #[inline] pub fn varg(self) -> f64 { self.m.he_arg_on_vert[self.id as usize] }
    #[inline] pub fn farg(self) -> f64 { self.m.he_arg_on_face[self.id as usize] }
    #[inline] pub fn darg(self) -> f64 { self.m.dihedral_arg[self.id as usize] }
    #[inline] pub fn vec(self) -> Row3d { self.head().pos() - self.tail().pos() }
    #[inline] pub fn nml(self) -> Row3d {
        let d = self.face().normal() + self.twin().face().normal();
        if d.length() > 0. { normalized(d) } else { Row3d::ZERO }
    }
    #[inline] pub fn lerp(self, r: f64) -> Row3d { self.tail().pos() * (1. - r) + self.head().pos() * r }
    #[inline] pub fn uv(self, r: f64) -> Complex { lerp(self.cr_t().uv(), self.cr_h().uv(), r) }
}

impl<'a> Vert<'a> {
    #[inline] pub fn half(self) -> Half<'a> { Half { id: self.m.vert2half[self.id as usize], m: self.m } }
    #[inline] pub fn pos(self) -> Row3d { self.m.pos[self.id as usize] }
    #[inline] pub fn basis_x(self) -> Row3d { self.m.vert_basis_x[self.id as usize] }
    #[inline] pub fn basis_y(self) -> Row3d { self.m.vert_basis_y[self.id as usize] }
    #[inline] pub fn normal(self) -> Row3d { self.m.vert_normal[self.id as usize] }
    #[inline] pub fn is_boundary(self) -> bool { self.m.is_bv[self.id as usize] }
    #[inline] pub fn is_interior(self) -> bool { !self.m.is_bv[self.id as usize] }
    #[inline] pub fn bary_area(self) -> f64 { self.m.bary_dual_area[self.id as usize] }
    #[inline] pub fn circ_area(self) -> f64 { self.m.circ_dual_area[self.id as usize] }
    /// `adjHalfs(ccw = true)`
    #[inline] pub fn adj_halfs(self) -> AdjIter<'a> { self.adj_halfs_dir(true) }
    #[inline] pub fn adj_halfs_dir(self, ccw: bool) -> AdjIter<'a> {
        AdjIter::new(self.m, self.m.vert2half[self.id as usize], ccw, AdjKind::VH)
    }
    /// `adjHalfs(Half h, ccw = true)`
    #[inline] pub fn adj_halfs_from(self, h: Half<'a>, ccw: bool) -> AdjIter<'a> {
        metriko_check!(h.tail().id == self.id, "half {} is not on vert {}", h.id, self.id);
        AdjIter::new(self.m, h.id, ccw, AdjKind::VH)
    }
    /// `adjCrnrs(ccw = true)`: the crnrs of the adjHalfs, one per half
    #[inline] pub fn adj_crnrs(self) -> impl Iterator<Item = Crnr<'a>> + 'a { self.adj_crnrs_dir(true) }
    #[inline] pub fn adj_crnrs_dir(self, ccw: bool) -> impl Iterator<Item = Crnr<'a>> + 'a {
        self.adj_halfs_dir(ccw).map(|h| h.cr_t())
    }
}

impl<'a> Edge<'a> {
    #[inline] pub fn half(self) -> Half<'a> { Half { id: self.m.edge2half[self.id as usize], m: self.m } }
    #[inline] pub fn vert0(self) -> Vert<'a> { Vert { id: self.m.edge2vert[self.id as usize][0], m: self.m } }
    #[inline] pub fn vert1(self) -> Vert<'a> { Vert { id: self.m.edge2vert[self.id as usize][1], m: self.m } }
    #[inline] pub fn face0(self) -> Face<'a> { Face { id: self.m.edge2face[self.id as usize][0], m: self.m } }
    #[inline] pub fn face1(self) -> Face<'a> { Face { id: self.m.edge2face[self.id as usize][1], m: self.m } }
    #[inline] pub fn len(self) -> f64 { Half { id: self.m.edge2half[self.id as usize], m: self.m }.len() }
    #[inline] pub fn cot(self) -> f64 { self.m.edge_cotan[self.id as usize] }
    #[inline] pub fn is_boundary(self) -> bool { self.half().is_boundary() || self.half().twin().is_boundary() }
    #[inline] pub fn vec(self) -> Row3d { self.vert1().pos() - self.vert0().pos() }
    #[inline] pub fn nml(self) -> Row3d {
        let d = self.face0().normal() + self.face1().normal();
        if d.length() > 0. { normalized(d) } else { Row3d::ZERO }
    }
    #[inline] pub fn lerp(self, r: f64) -> Row3d { self.vert0().pos() * (1. - r) + self.vert1().pos() * r }
}

impl<'a> Face<'a> {
    #[inline] pub fn half(self) -> Half<'a> { Half { id: self.m.face2half[self.id as usize], m: self.m } }
    #[inline] pub fn basis_x(self) -> Row3d { self.m.face_basis_x[self.id as usize] }
    #[inline] pub fn basis_y(self) -> Row3d { self.m.face_basis_y[self.id as usize] }
    #[inline] pub fn normal(self) -> Row3d { self.m.face_normal[self.id as usize] }
    #[inline] pub fn center(self) -> Row3d { self.m.bary_center[self.id as usize] }
    #[inline] pub fn area(self) -> f64 { self.m.face_area[self.id as usize] }
    pub fn to_local(self, v: Row3d) -> Complex {
        let d = v - self.half().tail().pos();
        cx(d.dot(self.basis_x()), d.dot(self.basis_y()))
    }
    pub fn to_world(self, v: Complex) -> Row3d {
        self.half().tail().pos() + self.basis_x() * v.re + self.basis_y() * v.im
    }
    pub fn uv2pos(self, uv: Complex) -> Row3d {
        let [a, b, c] = self.crnrs();
        conversion_2d_3d(a.uv(), b.uv(), c.uv(), a.pos(), b.pos(), c.pos(), uv)
    }
    /// `adjHalfs(ccw = true)`
    #[inline] pub fn adj_halfs(self) -> AdjIter<'a> { self.adj_halfs_dir(true) }
    #[inline] pub fn adj_halfs_dir(self, ccw: bool) -> AdjIter<'a> {
        AdjIter::new(self.m, self.m.face2half[self.id as usize], ccw, AdjKind::FH)
    }
    /// `adjHalfs(Half h, ccw = true)`
    #[inline] pub fn adj_halfs_from(self, h: Half<'a>, ccw: bool) -> AdjIter<'a> {
        metriko_check!(h.face().id == self.id, "half {} is not on face {}", h.id, self.id);
        AdjIter::new(self.m, h.id, ccw, AdjKind::FH)
    }
    /// curr, next, prev
    #[inline] pub fn halfs(self) -> [Half<'a>; 3] { let h = self.half(); [h, h.next(), h.prev()] }
    /// tails of halfs above
    #[inline] pub fn verts(self) -> [Vert<'a>; 3] { let [h0, h1, h2] = self.halfs(); [h0.tail(), h1.tail(), h2.tail()] }
    /// edges of halfs above
    #[inline] pub fn edges(self) -> [Edge<'a>; 3] { let [h0, h1, h2] = self.halfs(); [h0.edge(), h1.edge(), h2.edge()] }
    /// crnrs of halfs above
    #[inline] pub fn crnrs(self) -> [Crnr<'a>; 3] { let [h0, h1, h2] = self.halfs(); [h0.crnr(), h1.crnr(), h2.crnr()] }
}

impl<'a> Crnr<'a> {
    #[inline] pub fn half(self) -> Half<'a> { Half { id: self.m.crnr2half[self.id as usize], m: self.m } }
    #[inline] pub fn vert(self) -> Vert<'a> { self.half().prev().tail() }
    #[inline] pub fn face(self) -> Face<'a> { self.half().face() }
    /// next crnr of the same face in ccw
    #[inline] pub fn next(self) -> Crnr<'a> { self.half().cr_t() }
    /// prev crnr of the same face in ccw
    #[inline] pub fn prev(self) -> Crnr<'a> { self.half().cr_h() }
    #[inline] pub fn pos(self) -> Row3d { self.vert().pos() }
    #[inline] pub fn uv(self) -> Complex { self.m.cfn[self.id as usize] }
}

impl<'a> Loop<'a> {
    #[inline] pub fn half(self) -> Half<'a> { Half { id: self.m.loop2half[self.id as usize], m: self.m } }
    #[inline] pub fn adj_halfs(self) -> AdjIter<'a> { self.adj_halfs_dir(true) }
    #[inline] pub fn adj_halfs_dir(self, ccw: bool) -> AdjIter<'a> {
        AdjIter::new(self.m, self.m.loop2half[self.id as usize], ccw, AdjKind::FH)
    }
}

impl Hmesh {
    #[inline] pub fn half(&self, id: i32) -> Half<'_> { Half { id, m: self } }
    #[inline] pub fn vert(&self, id: i32) -> Vert<'_> { Vert { id, m: self } }
    #[inline] pub fn edge(&self, id: i32) -> Edge<'_> { Edge { id, m: self } }
    #[inline] pub fn face(&self, id: i32) -> Face<'_> { Face { id, m: self } }
    #[inline] pub fn crnr(&self, id: i32) -> Crnr<'_> { Crnr { id, m: self } }
    #[inline] pub fn lp(&self, id: i32) -> Loop<'_> { Loop { id, m: self } }
    pub fn halfs(&self) -> impl Iterator<Item = Half<'_>> + '_ { (0..self.n_h).map(move |id| Half { id, m: self }) }
    pub fn verts(&self) -> impl Iterator<Item = Vert<'_>> + '_ { (0..self.n_v).map(move |id| Vert { id, m: self }) }
    pub fn edges(&self) -> impl Iterator<Item = Edge<'_>> + '_ { (0..self.n_e).map(move |id| Edge { id, m: self }) }
    pub fn faces(&self) -> impl Iterator<Item = Face<'_>> + '_ { (0..self.n_f).map(move |id| Face { id, m: self }) }
    pub fn crnrs(&self) -> impl Iterator<Item = Crnr<'_>> + '_ { (0..self.n_c).map(move |id| Crnr { id, m: self }) }
    pub fn loops(&self) -> impl Iterator<Item = Loop<'_>> + '_ { (0..self.n_l).map(move |id| Loop { id, m: self }) }
}

/// Output of `dcel`: VH, EH, FH, HV, HE, HF, nextH, prevH, twinH.
struct Dcel {
    vh: Vec<i32>,
    eh: Vec<[i32; 2]>,
    fh: Vec<[i32; 3]>,
    hv: Vec<i32>,
    he: Vec<i32>,
    hf: Vec<i32>,
    next_h: Vec<i32>,
    prev_h: Vec<i32>,
    twin_h: Vec<i32>,
}

fn dcel(
    d: &[i32],
    f: &[[i32; 3]],
    ev: &[[i32; 2]],
    ef: &[[i32; 2]],
    efi: &[[i32; 2]],
) -> Dcel {
    let mut eh = vec![[-1i32; 2]; ev.len()];
    let mut num_h = 0i32;

    for i in 0..ef.len() {
        if ef[i][0] != -1 { eh[i][0] = num_h; num_h += 1; }
        if ef[i][1] != -1 { eh[i][1] = num_h; num_h += 1; }
    }

    let mut he = vec![-1i32; num_h as usize];
    for i in 0..eh.len() {
        if eh[i][0] != -1 { he[eh[i][0] as usize] = i as i32; }
        if eh[i][1] != -1 { he[eh[i][1] as usize] = i as i32; }
    }

    // DEVIATION: the C++ `conservativeResize`s leave entries never written uninitialized (only reachable for vertices
    // in no edge, which a valid mesh does not have); here they are -1.
    let max_v = ev.iter().flat_map(|e| e.iter()).copied().max().unwrap_or(-1);
    let mut hv = vec![-1i32; num_h as usize];
    let mut vh = vec![-1i32; (max_v + 1) as usize];
    for i in 0..ev.len() {
        if eh[i][0] != -1 { hv[eh[i][0] as usize] = ev[i][0]; vh[ev[i][0] as usize] = eh[i][0]; }
        if eh[i][1] != -1 { hv[eh[i][1] as usize] = ev[i][1]; vh[ev[i][1] as usize] = eh[i][1]; }
    }

    let mut twin_h = vec![-1i32; num_h as usize];
    for i in 0..eh.len() {
        if eh[i][0] != -1 && eh[i][1] != -1 {
            twin_h[eh[i][0] as usize] = eh[i][1];
            twin_h[eh[i][1] as usize] = eh[i][0];
        }
    }

    let mut fh = vec![[-1i32; 3]; f.len()];
    let mut hf = vec![-1i32; num_h as usize];
    for i in 0..ef.len() {
        if ef[i][0] != -1 { fh[ef[i][0] as usize][efi[i][0] as usize] = eh[i][0]; hf[eh[i][0] as usize] = ef[i][0]; }
        if ef[i][1] != -1 { fh[ef[i][1] as usize][efi[i][1] as usize] = eh[i][1]; hf[eh[i][1] as usize] = ef[i][1]; }
    }

    let mut next_h = vec![-1i32; he.len()];
    let mut prev_h = vec![-1i32; he.len()];
    for i in 0..d.len() {
        let di = d[i] as usize;
        for j in 0..di {
            next_h[fh[i][j] as usize] = fh[i][(j + 1) % di];
            prev_h[fh[i][(j + 1) % di] as usize] = fh[i][j];
        }
    }
    Dcel { vh, eh, fh, hv, he, hf, next_h, prev_h, twin_h }
}

impl Hmesh {
    /// `Hmesh(V, F, only_topology = false)`
    pub fn new(v: &[Row3d], f: &[[i32; 3]], only_topology: bool) -> Hmesh {
        let (edge2vert, face2edge, edge2face) = edge_topology(v.len(), f);
        let pos = v.to_vec();
        let idx = f.to_vec();
        let n_v = pos.len() as i32;
        let n_f = idx.len() as i32;
        let n_e = edge2vert.len() as i32;
        let n_h = (edge2vert.len() * 2) as i32;
        let n_c = n_f * 3;
        let n_eular_chars = n_v - n_e + n_f;
        let mut m = Hmesh {
            n_v, n_f, n_e, n_h, n_l: 0, n_c, n_eular_chars,
            vert2half: vec![-1; n_v as usize],
            edge2half: vec![-1; n_e as usize],
            face2half: vec![-1; n_f as usize],
            loop2half: vec![],
            crnr2half: vec![-1; n_c as usize],
            next: vec![-1; n_h as usize],
            prev: vec![-1; n_h as usize],
            twin: vec![-1; n_h as usize],
            tail: vec![-1; n_h as usize],
            head: vec![-1; n_h as usize],
            edge: vec![-1; n_h as usize],
            face: vec![-1; n_h as usize],
            crnr: vec![-1; n_h as usize],
            is_bv: vec![],
            pos,
            idx,
            cfn: vec![],
            edge2vert,
            edge2face,
            face2edge,
            vert_basis_x: vec![],
            vert_basis_y: vec![],
            vert_normal: vec![],
            face_basis_x: vec![],
            face_basis_y: vec![],
            face_normal: vec![],
            face_area: vec![],
            bary_dual_area: vec![],
            circ_dual_area: vec![],
            half_cotan: vec![],
            edge_cotan: vec![],
            bary_center: vec![],
            angle_defect: vec![],
            dihedral_arg: vec![],
            he_arg_on_vert: vec![],
            he_arg_on_face: vec![],
        };

        let n_p = 3usize;

        // set up topology
        let d = vec![3i32; m.idx.len()];
        let mut efi = vec![[-1i32; 2]; m.edge2face.len()]; // number of an edge inside the face
        for i in 0..m.edge2face.len() {
            for k in 0..2 {
                if m.edge2face[i][k] == -1 { continue; }
                for j in 0..3 { if m.face2edge[m.edge2face[i][k] as usize][j] == i as i32 { efi[i][k] = j as i32; } }
            }
        }
        let dc = dcel(&d, &m.idx, &m.edge2vert, &m.edge2face, &efi);

        for iv in 0..n_v as usize { m.vert2half[iv] = dc.vh[iv]; }
        for ie in 0..n_e as usize { m.edge2half[ie] = dc.eh[ie][0]; }
        for i_f in 0..n_f as usize { m.face2half[i_f] = dc.fh[i_f][0]; }
        for ih in 0..dc.hv.len() {
            m.head[ih] = dc.hv[dc.next_h[ih] as usize];
            m.tail[ih] = dc.hv[ih];
            m.edge[ih] = dc.he[ih];
            m.face[ih] = dc.hf[ih];
            m.next[ih] = dc.next_h[ih];
            m.prev[ih] = dc.prev_h[ih];
            if dc.twin_h[ih] != -1 { m.twin[ih] = dc.twin_h[ih]; }
        }

        // setup boundaries
        let n_fp = n_f as usize * n_p;
        let mut ib = n_fp;
        for ih in 0..n_fp {
            if m.twin[ih] != -1 { continue; }
            let mut jh = ih;
            let jb = ib;
            'walk: loop {
                m.tail[ib] = m.head[jh];
                m.head[ib] = m.tail[jh];
                m.edge[ib] = m.edge[jh];
                m.prev[ib] = ib as i32 - 1;
                m.next[ib] = ib as i32 + 1;
                m.twin[jh] = ib as i32;
                m.twin[ib] = jh as i32;

                jh = m.prev[jh] as usize;
                while m.twin[jh] != -1 {
                    if jh == ih { break 'walk; }
                    jh = m.prev[m.twin[jh] as usize] as usize;
                }
                ib += 1;
            }
            // loop_done:
            m.prev[jb] = ib as i32;
            m.next[ib] = jb as i32;
            m.loop2half.push(jb as i32);
            ib += 1;
        }
        m.n_l = m.loop2half.len() as i32;

        // setup edge->half on boundary
        for ie in 0..n_e as usize {
            if m.edge2half[ie] == -1 { m.edge2half[ie] = m.twin[dc.eh[ie][1] as usize]; }
        }

        // setup crnr_half & half_crnr
        for i_f in 0..n_f as usize {
            for it in 0..3 {
                let vid = f[i_f][it];
                let cid = (i_f * 3 + it) as i32;
                let hs = m.face(i_f as i32).halfs().map(|h| h.id);
                for hid in hs {
                    let (t, h) = (m.tail[hid as usize], m.head[hid as usize]);
                    if t != vid && h != vid {
                        m.crnr2half[cid as usize] = hid;
                        m.crnr[hid as usize] = cid;
                    }
                }
            }
        }

        // set up a filter to check whether a vertex is on boundary or not
        let mut is_bv = vec![false; n_v as usize];
        for e in m.edges() {
            if e.is_boundary() {
                is_bv[e.half().tail().id as usize] = true;
                is_bv[e.half().head().id as usize] = true;
            }
        }
        m.is_bv = is_bv;

        if only_topology { return m; }

        let mut vert_normal = vec![Row3d::ZERO; n_v as usize];
        let mut vert_basis_x = vec![Row3d::ZERO; n_v as usize];
        let mut vert_basis_y = vec![Row3d::ZERO; n_v as usize];
        let mut face_basis_x = vec![Row3d::ZERO; n_f as usize];
        let mut face_basis_y = vec![Row3d::ZERO; n_f as usize];
        let mut face_normal = vec![Row3d::ZERO; n_f as usize];
        let mut face_area = vec![0.; n_f as usize];
        let mut half_cotan = vec![0.; n_h as usize];
        let mut edge_cotan = vec![0.; n_e as usize];
        let mut angle_defect = vec![0.; n_v as usize];
        let mut dihedral_arg = vec![0.; n_h as usize];
        let mut he_arg_on_vert = vec![0.; n_h as usize];
        let mut he_arg_on_face = vec![0.; n_h as usize];
        let mut bary_center = vec![Row3d::ZERO; n_f as usize];
        let mut bary_dual_area = vec![0.; n_v as usize];
        let mut circ_dual_area = vec![0.; n_v as usize];

        for f in m.faces() {
            let x = f.half().vec();
            let t = f.half().prev().vec() * -1.;
            let n = x.cross(t);
            let mut p = Row3d::ZERO;

            // set up  basis of face
            face_basis_x[f.id as usize] = normalized(x);
            face_basis_y[f.id as usize] = -normalized(x.cross(n));
            face_normal[f.id as usize] = normalized(n);

            // set up barycenter of face
            for h in f.adj_halfs() { p += h.tail().pos(); }
            bary_center[f.id as usize] = p / n_p as f64;

            // set up area of face
            face_area[f.id as usize] = n.length() * 0.5;
        }
        m.face_basis_x = face_basis_x;
        m.face_basis_y = face_basis_y;
        m.face_normal = face_normal;
        m.bary_center = bary_center;
        m.face_area = face_area;

        // set up halfedge angle on face coordinate
        for f in m.faces() {
            let mut sum = 0.;
            for h in f.adj_halfs() {
                let v1 = normalized(h.vec());
                let v2 = normalized(h.prev().vec());
                if (v1 - f.basis_x()).length() > 1e-10 { sum += v1.dot(v2).acos(); }
                he_arg_on_face[h.id as usize] = sum;
            }
        }
        m.he_arg_on_face = he_arg_on_face;

        // set up vertex orthogonal coordinate
        for i_f in 0..n_f as usize {
            for i_p in 0..n_p {
                let i_curr = m.idx[i_f][i_p] as usize;
                let i_next = m.idx[i_f][(i_p + n_p + 1) % n_p] as usize;
                let i_prev = m.idx[i_f][(i_p + n_p - 1) % n_p] as usize;
                let d0 = normalized(m.pos[i_next] - m.pos[i_curr]);
                let d1 = normalized(m.pos[i_prev] - m.pos[i_curr]);
                let d = d1.dot(d0);
                let phi = if d >= 1. { 0. } else if d <= -1. { PI } else { d.acos() };
                vert_normal[m.idx[i_f][i_p] as usize] += m.face_normal[i_f] * phi;
            }
        }
        for n in vert_normal.iter_mut() { *n = normalized(*n); }
        m.vert_normal = vert_normal;

        for v in m.verts() {
            let v_ = m.pos[v.half().head().id as usize] - m.pos[v.id as usize];
            let n = m.vert_normal[v.id as usize];
            let x = normalized(v_ - v_.dot(n) * n);
            vert_basis_x[v.id as usize] = x;
            vert_basis_y[v.id as usize] = n.cross(x);
        }
        m.vert_basis_x = vert_basis_x;
        m.vert_basis_y = vert_basis_y;

        // set up a dihedral angle for each halfedge
        for h in m.halfs() {
            if h.edge().is_boundary() { dihedral_arg[h.id as usize] = 0.; continue; }
            let n1 = m.face_normal[h.face().id as usize];
            let n2 = m.face_normal[h.twin().face().id as usize];
            let v = h.vec() / h.len();
            let c = n1.cross(n2);
            dihedral_arg[h.id as usize] = v.dot(c).atan2(n1.dot(n2));
        }
        m.dihedral_arg = dihedral_arg;

        // set up halfedge cotan and edge cotan
        for h in m.halfs() {
            if h.is_boundary() { half_cotan[h.id as usize] = 0.; continue; }
            let vn = h.next().vec();
            let vp = h.prev().vec() * -1.;
            half_cotan[h.id as usize] = vp.dot(vn) / vp.cross(vn).length();
        }
        m.half_cotan = half_cotan;
        for e in m.edges() {
            edge_cotan[e.id as usize] = (e.half().cot() + e.half().twin().cot()) * 0.5;
        }
        m.edge_cotan = edge_cotan;

        // set up dual area of vertex
        for v in m.verts() {
            let mut a1 = 0.;
            let mut a2 = 0.;
            for h in v.adj_halfs() {
                if h.face().id != -1 { a1 += m.face_area[h.face().id as usize]; }
                a2 += h.cot() * h.vec().length_squared();
                a2 += h.prev().cot() * h.prev().vec().length_squared();
            }
            bary_dual_area[v.id as usize] = a1 / 3.;
            circ_dual_area[v.id as usize] = a2 * 0.125;
        }
        m.bary_dual_area = bary_dual_area;
        m.circ_dual_area = circ_dual_area;

        // set up halfedge angle on vertex coordinate, and angle defect on vertex
        for v in m.verts() {
            let mut sum = 0.;
            for h in v.adj_halfs() {
                he_arg_on_vert[h.id as usize] = sum;
                let v1 = normalized(h.vec());
                let v2 = normalized(h.prev().twin().vec());
                sum += v1.dot(v2).acos();
            }
            angle_defect[v.id as usize] = if v.is_boundary() { PI - sum } else { TWO_PI - sum };
            for h in v.adj_halfs() { he_arg_on_vert[h.id as usize] *= TWO_PI / sum; }
        }
        m.he_arg_on_vert = he_arg_on_vert;
        m.angle_defect = angle_defect;
        m
    }
}
