//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/nvec/face_rosy_field.h (see rust/PORTED_FROM).

use std::collections::{BTreeSet, VecDeque};
use faer::{Col, Mat};
use faer::sparse::Triplet;
use crate::common::typedef::*;
use crate::hmesh::hmesh::*;
use crate::solver::matrix_ops::{SprsC, TripC, from_triplets, nonzeros};
use crate::solver::sparse_solver::{solve_hermitian, solve_smallest_eig};
use crate::igl::adjacency_list::adjacency_list;
use crate::igl::cut_mesh_from_singularities::cut_mesh_from_singularities;
use crate::igl::dijkstra::{dijkstra, dijkstra_path};
use crate::igl::triangle_triangle_adjacency::triangle_triangle_adjacency;
use crate::igl::vertex_triangle_adjacency::vertex_triangle_adjacency;

use crate::common::stdcomplex::{std_div, std_polar, std_pow, std_round};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldType { UnSpecified, Smoothest, CurvatureAligned }

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatchingType { Principal, Curl }

pub struct FaceRosyField<'a> {
    pub mesh: &'a Hmesh,
    pub rosy_n: i32,
    /// `MatXc field`: nF x rosyN
    pub field: Mat<Complex>,
    pub connection: Col<Complex>,
    pub compressed: Col<Complex>,
    pub matching: Vec<i32>,
    pub singular: Vec<i32>,
}

impl<'a> FaceRosyField<'a> {
    /// `FaceRosyField(const Hmesh &m, const int nRosy)`
    pub fn new_empty(m: &'a Hmesh, n_rosy: i32) -> Self {
        FaceRosyField {
            mesh: m,
            rosy_n: n_rosy,
            field: Mat::<Complex>::zeros(0, 0),
            connection: Col::<Complex>::zeros(0),
            compressed: Col::<Complex>::zeros(0),
            matching: vec![],
            singular: vec![],
        }
    }

    /// `FaceRosyField(const Hmesh &m, const int nRosy, FieldType type, double guidanceSmooth = 4.)`
    pub fn new(m: &'a Hmesh, n_rosy: i32, ty: FieldType) -> Self { Self::new_with_smooth(m, n_rosy, ty, 4.) }

    // guidanceSmooth: how far the per-face curvature estimate is diffused before it is used as
    // alignment guidance, in mean face areas (the radius is roughly its square root in mean
    // edge lengths). 0 keeps the raw per-face estimate
    pub fn new_with_smooth(m: &'a Hmesh, n_rosy: i32, ty: FieldType, guidance_smooth: f64) -> Self {
        let mut s = Self::new_empty(m, n_rosy);
        let mesh = s.mesh;
        let rosy_n = s.rosy_n;
        let l = s.connection_laplacian();
        // SprsC(mesh.faceArea.cast<complex>().asDiagonal())
        let mt: Vec<TripC> = (0..mesh.n_f as usize).map(|i| Triplet::new(i, i, cx(mesh.face_area[i], 0.))).collect();
        let m_ = from_triplets(mesh.n_f as usize, mesh.n_f as usize, &mt);
        match ty {
            FieldType::Smoothest => {
                s.compressed = solve_smallest_eig(&l, &m_, 500, 1e-8);
                //int unitSmoothIter = 32
                // the eigenproblem minimises the dirichlet energy relative to the L2 norm, so
                // |psi| is free and sags to zero over whole regions; the field that is actually
                // used throws |psi| away and is therefore not a critical point of the energy it
                // is meant to minimise. the constrained problem, min psi* L psi with |psi_f| = 1,
                // is the harmonic cross field of viertel and osting 2019, and its diffusion
                // generated form is just diffuse, renormalise, repeat. sweeping the diffusion
                // time from coarse to fine plays the role of their ginzburg-landau annealing and
                // lets close singularity pairs annihilate on the way down. measured: dirichlet
                // -23..-38%, curl -28..-43%, nefertiti 126 -> 108 singularities, and its poisson
                // solution becomes locally injective so the whole LM stage drops out
                // (the disabled C++ block that implements it is not ported)
            }
            FieldType::CurvatureAligned => {
                if rosy_n != 2 && rosy_n != 4 { panic!("curvature alignment needs a 2- or 4-rosy field"); }
                let lambda = 0.;
                let mut d = s.principal_curvature_dir();
                // D.array().square()
                if rosy_n == 4 { for i in 0..d.nrows() { d[i] = d[i] * d[i]; } }
                // the estimate above only sees one face's three dihedral angles, so on an
                // irregular triangulation it is dominated by meshing noise. diffusing it first
                // amounts to estimating the curvature over a neighbourhood: on nefertiti it
                // takes the field from 205 to 138 singularities and cuts the curl by 21% for a
                // 3% loss of alignment, and leaves fandisk and spot essentially unchanged
                if guidance_smooth > 0. {
                    // mesh.faceArea.mean(): a plain sum (Eigen's vectorized reduction only differs in the last bits)
                    let mut sum = 0.;
                    for &a in &mesh.face_area { sum += a; }
                    let mean = sum / mesh.n_f as f64;
                    let lhs_ = sparse_add(&m_, &sparse_scale(guidance_smooth * mean, &l));
                    let rhs_: Col<Complex> = &m_ * &d;
                    d = solve_hermitian(&lhs_, &rhs_);
                }
                // M * D / sqrt(abs((D.adjoint() * M * D)[0]))
                let md: Col<Complex> = &m_ * &d;
                let mut dmd = cx(0., 0.);
                {
                    // D.adjoint() * M: the row vector conj(D)^T M, then its product with D
                    let mut r = vec![cx(0., 0.); d.nrows()];
                    for (i, j, v) in nonzeros(m_.as_ref()) { r[j] += d[i].conj() * v; }
                    for j in 0..d.nrows() { dmd += r[j] * d[j]; }
                }
                let den = dmd.norm().sqrt();
                let rhs = Col::<Complex>::from_fn(md.nrows(), |i| md[i] / den);
                let lhs = sparse_sub(&l, &sparse_scale(lambda, &m_));
                s.compressed = solve_hermitian(&lhs, &rhs);
            }
            FieldType::UnSpecified => {}
        }
        let r = std_polar(1., TWO_PI / rosy_n as f64);
        let mut field = Mat::<Complex>::zeros(mesh.n_f as usize, rosy_n as usize);
        for f in mesh.faces() {
            for i in 0..rosy_n {
                field[(f.id as usize, i as usize)] =
                    std_pow(std_polar(1., s.compressed[f.id as usize].arg()), cx(1. / rosy_n as f64, 0.)) * std_pow(r, cx(i as f64, 0.));
            }
        }
        s.field = field;
        s.compute_matching(MatchingType::Principal);
        s
    }

    pub fn connection_laplacian(&self) -> SprsC {
        let mesh = self.mesh;
        let mut transport = vec![cx(0., 0.); mesh.n_h as usize];
        let mut t: Vec<TripC> = vec![];

        for e in mesh.edges() {
            let h1 = e.half();
            let h2 = e.half().twin();
            let v1 = std_polar(1., h1.farg());
            let v2 = std_polar(1., h2.farg());
            let r = std_div(-v2, v1);
            transport[h1.id as usize] = r;
            transport[h2.id as usize] = std_div(cx(1., 0.), r);
        }

        for f in mesh.faces() {
            for h in f.adj_halfs() {
                if h.twin().is_boundary() { continue; }
                let w = 1.;
                let r = transport[h.twin().id as usize];
                let twin_f = h.twin().face();
                t.push(Triplet::new(f.id as usize, f.id as usize, cx(w, 0.)));
                t.push(Triplet::new(f.id as usize, twin_f.id as usize, std_pow(r, cx(self.rosy_n as f64, 0.)) * -w));
            }
        }

        let n = mesh.n_f as usize;
        let s = from_triplets(n, n, &t);
        // S + 1e-9 * I: the identity scaled by the real 1e-9 is (1e-9, 0) on the diagonal; the sum adds the two
        // operands entry by entry over the union of their patterns
        let mut i_ = Vec::with_capacity(n);
        for k in 0..n { i_.push(Triplet::new(k, k, cx(1e-9 * 1., 1e-9 * 0.))); }
        let i_ = from_triplets(n, n, &i_);
        sparse_add(&s, &i_)
    }

    pub fn principal_curvature_dir(&self) -> Col<Complex> {
        let mesh = self.mesh;
        let mut d = Col::<Complex>::zeros(mesh.n_f as usize);
        for f in mesh.faces() {
            let mut dir = cx(0., 0.);
            for h in f.adj_halfs() {
                let l = h.len();
                let c = std_polar(1., h.farg()) * l;
                dir += -c * c / l * h.darg();
            }
            d[f.id as usize] = dir * 0.25;
        }
        d
    }

    pub fn compute_singular(&mut self, effort: &[f64]) {
        let mesh = self.mesh;
        self.singular = vec![1; mesh.n_v as usize];
        for v in mesh.verts().filter(|v| !v.is_boundary()) {
            let mut sum = mesh.angle_defect[v.id as usize] * self.rosy_n as f64;
            for h in v.adj_halfs() { sum += (if h.is_canonical() { -1. } else { 1. }) * effort[h.edge().id as usize]; }
            self.singular[v.id as usize] = std_round(sum / TWO_PI) as i32;
        }
    }

    fn compute_connection(&mut self) {
        let mesh = self.mesh;
        self.connection = Col::<Complex>::zeros(mesh.n_e as usize);
        for e in mesh.edges().filter(|e| !e.is_boundary()) {
            let v = normalized(e.half().vec());
            let ef = cx(v.dot(e.face0().basis_x()), v.dot(e.face0().basis_y()));
            let eg = cx(v.dot(e.face1().basis_x()), v.dot(e.face1().basis_y()));
            self.connection[e.id as usize] = std_div(eg, ef);
        }
    }

    fn compute_matching(&mut self, ty: MatchingType) {
        self.compute_connection();
        self.matching = vec![-1; self.mesh.n_e as usize];
        let mut effort = vec![0.; self.mesh.n_e as usize];
        match ty {
            MatchingType::Principal => self.p_matching(&mut effort),
            MatchingType::Curl => self.c_matching(&mut effort),
        }
        self.compute_singular(&effort);
    }

    fn p_matching(&mut self, effort: &mut [f64]) {
        let mesh = self.mesh;
        let n = self.rosy_n;
        for e in mesh.edges().filter(|e| !e.is_boundary()) {
            let f0 = e.face0().id as usize;
            let f1 = e.face1().id as usize;
            let conn = self.connection[e.id as usize];
            let mut min_rot = 1000.;
            let mut offset = 0;
            let mut coef = cx(1., 0.);
            let transport0 = conn * self.field[(f0, 0)];
            for i in 0..n {
                let jf = self.field[(f0, i as usize)];
                let jg = self.field[(f1, i as usize)];
                coef *= std_div(jg, conn * jf);
                let r = std_div(jg, transport0).arg();
                if r.abs() < f64::abs(min_rot) {
                    offset = i;
                    min_rot = r;
                }
            }
            effort[e.id as usize] = coef.arg();

            let mut residuals = 0.;
            for i in 0..n {
                let v0 = self.field[(f0, i as usize)];
                let v1 = self.field[(f1, ((i + offset) % n) as usize)];
                residuals += std_div(v1, conn * v0).arg();
            }
            // int = int - double: the double is truncated toward zero on assignment
            self.matching[e.id as usize] = (offset as f64 - std_round((residuals - effort[e.id as usize]) / TWO_PI)) as i32;
            //std::cout << residuals - effort[e.id] << std::endl;
            //matching[e.id] = offset;
        }
    }

    fn c_matching(&mut self, effort: &mut [f64]) {
        let mesh = self.mesh;
        let n = self.rosy_n;
        for e in mesh.edges().filter(|e| !e.is_boundary()) {
            let mut i_m = 0;
            let mut min = 32767000.;
            let v = normalized(e.half().vec());
            let f0 = e.face0();
            let f1 = e.face1();

            for i in 0..n {
                let mut curl = 0.;
                for j in 0..n {
                    let c0 = self.field[(f0.id as usize, j as usize)];
                    let c1 = self.field[(f1.id as usize, ((i + j) % n) as usize)];
                    let v0 = c0.re * f0.basis_x() + c0.im * f0.basis_y();
                    let v1 = c1.re * f1.basis_x() + c1.im * f1.basis_y();
                    curl += v.dot(v1 - v0).powf(2.);
                }
                if curl < min { i_m = i; min = curl; }
            }

            self.matching[e.id as usize] = i_m;

            let mut coef = cx(1., 0.);
            for i in 0..n {
                let v0 = self.field[(e.face0().id as usize, i as usize)];
                let v1 = self.field[(e.face1().id as usize, ((i + i_m) % n) as usize)];
                coef *= std_div(v1, v0 * self.connection[e.id as usize]);
            }
            effort[e.id as usize] = coef.arg();
        }
    }
}

/// Eigen's sparse `A + B`: entry by entry over the union of the patterns, `a + b` where both are stored.
fn sparse_add(a: &SprsC, b: &SprsC) -> SprsC {
    let mut t: Vec<TripC> = nonzeros(a.as_ref()).map(|(i, j, v)| Triplet::new(i, j, v)).collect();
    t.extend(nonzeros(b.as_ref()).map(|(i, j, v)| Triplet::new(i, j, v)));
    from_triplets(a.nrows(), a.ncols(), &t)
}

/// Eigen's sparse `A - B`: `a - b` where both are stored, `-b` where only `b` is.
fn sparse_sub(a: &SprsC, b: &SprsC) -> SprsC {
    let mut t: Vec<TripC> = vec![];
    let n = a.ncols();
    for j in 0..n {
        let ra = a.as_ref().row_idx_of_col_raw(j);
        let va = a.as_ref().val_of_col(j);
        let rb = b.as_ref().row_idx_of_col_raw(j);
        let vb = b.as_ref().val_of_col(j);
        let (mut p, mut q) = (0, 0);
        while p < ra.len() || q < rb.len() {
            if q == rb.len() || (p < ra.len() && ra[p] < rb[q]) { t.push(Triplet::new(ra[p], j, va[p])); p += 1; }
            else if p == ra.len() || rb[q] < ra[p] { t.push(Triplet::new(rb[q], j, -vb[q])); q += 1; }
            else { t.push(Triplet::new(ra[p], j, va[p] - vb[q])); p += 1; q += 1; }
        }
    }
    from_triplets(a.nrows(), a.ncols(), &t)
}

/// `double * SparseMatrix<complex>`: every stored value times the real scalar
fn sparse_scale(s: f64, a: &SprsC) -> SprsC {
    let t: Vec<TripC> = nonzeros(a.as_ref()).map(|(i, j, v)| Triplet::new(i, j, v * s)).collect();
    from_triplets(a.nrows(), a.ncols(), &t)
}

pub fn compute_combbed_field<'a>(field: &FaceRosyField<'a>, seam: &[bool]) -> FaceRosyField<'a> {
    let n = field.rosy_n;
    let mut fr = FaceRosyField::new_empty(field.mesh, field.rosy_n);
    // DEVIATION: Eigen's `resize` leaves the rows of faces the walk never reaches uninitialized (a mesh in several
    // pieces); here they are zero.
    fr.field = Mat::<Complex>::zeros(field.mesh.n_f as usize, n as usize);
    fr.connection = field.connection.clone();
    fr.compressed = field.compressed.clone();
    fr.singular = field.singular.clone();
    let mut turns = vec![0i32; fr.mesh.n_f as usize];
    let mut visit = vec![0i32; fr.mesh.n_f as usize];
    let mut matching_q: VecDeque<(i32, i32)> = VecDeque::new();
    matching_q.push_back((0, 0));

    loop {
        let (i, m) = matching_q.pop_front().unwrap();
        let r = n - m;
        if visit[i as usize] != 1 {
            visit[i as usize] = 1;
            turns[i as usize] = m;
            for k in 0..r { fr.field[(i as usize, k as usize)] = field.field[(i as usize, (m + k) as usize)]; }
            for k in 0..m { fr.field[(i as usize, (r + k) as usize)] = field.field[(i as usize, k as usize)]; }
            let face = fr.mesh.face(i);

            for h in face.adj_halfs() {
                let e = h.edge();
                let b = e.face0().id == i;
                let fn_ = if b { e.face1() } else { e.face0() };
                let mut nn = (if b { 1 } else { -1 }) * field.matching[e.id as usize];
                nn = (nn + m + 1000 * n) % n;
                if !e.is_boundary() && visit[fn_.id as usize] == 0 && !seam[e.id as usize] { matching_q.push_back((fn_.id, nn)); }
            }
        }
        if matching_q.is_empty() { break; }
    }

    fr.matching = vec![-1; fr.mesh.n_e as usize];
    for e in fr.mesh.edges().filter(|e| !e.is_boundary()) {
        fr.matching[e.id as usize] = (turns[e.face0().id as usize] - turns[e.face1().id as usize] + field.matching[e.id as usize] + 1000 * n) % n;
    }

    fr
}

/// `MatXd compute_extrinsic_field(field, rosyN)`: nF x 3 rosyN, the k-th direction of face f in columns 3k..3k+2
pub fn compute_extrinsic_field(field: &FaceRosyField<'_>, rosy_n: i32) -> Mat<f64> {
    metriko_check!(rosy_n as usize <= field.field.ncols(), "extrinsic field of {} directions from a {}-rosy field", rosy_n, field.field.ncols());
    let mut ext = Mat::<f64>::zeros(field.mesh.n_f as usize, 3 * rosy_n as usize);
    for f in field.mesh.faces() {
        for k in 0..rosy_n as usize {
            let c = field.field[(f.id as usize, k)];
            let d = normalized(c.re * f.basis_x() + c.im * f.basis_y());
            ext[(f.id as usize, 3 * k)] = d.x;
            ext[(f.id as usize, 3 * k + 1)] = d.y;
            ext[(f.id as usize, 3 * k + 2)] = d.z;
        }
    }
    ext
}

/// `cut_mesh_with_singularities(V, F, VF, VV, TT, TTi, singularities, cuts)`: returns `cuts`
pub fn cut_mesh_with_singularities_adj(
    v: &[Row3d],
    f: &[[i32; 3]],
    vf: &[Vec<i32>],
    vv: &[Vec<i32>],
    tt: &[[i32; 3]],
    tti: &[[i32; 3]],
    singularities: &[i32],
) -> Vec<[i32; 3]> {
    //first, get a spanning tree for the mesh (no missmatch needed)
    let zeros = vec![[0.; 3]; f.len()];
    let mut cuts = cut_mesh_from_singularities(v.len(), f, &zeros);

    let mut vertices_in_cut: BTreeSet<i32> = BTreeSet::new();
    for i in 0..cuts.len() {
        for j in 0..3 {
            if cuts[i][j] != 0 { vertices_in_cut.insert(f[i][j]); }
        }
    }

    //then, add all singularities one by using Dijkstra's algorithm
    for &sing in singularities {
        let mut path = vec![];
        let (vertex_found, _min_distance, previous) = dijkstra(sing, &vertices_in_cut, vv);
        if vertex_found == -1 {
            // this means that there are no cuts
            path.push(sing);
        } else {
            path = dijkstra_path(vertex_found, &previous);
        }

        vertices_in_cut.extend(path.iter().copied());

        //insert to cut
        for ii in 0..path.len() - 1 {
            let v0 = path[ii];
            let v1 = path[ii + 1];

            let mut vf0 = vf[v0 as usize].clone();
            vf0.sort_unstable();
            let mut vf1 = vf[v1 as usize].clone();
            vf1.sort_unstable();

            // rg::set_intersection of the sorted lists
            let mut common_face_v = vec![];
            let (mut p, mut q) = (0, 0);
            while p < vf0.len() && q < vf1.len() {
                if vf0[p] < vf1[q] { p += 1; }
                else if vf1[q] < vf0[p] { q += 1; }
                else { common_face_v.push(vf0[p]); p += 1; q += 1; }
            }
            metriko_check!(common_face_v.len() == 2, "cut edge {}-{} is shared by {} faces", v0, v1, common_face_v.len());

            let fi = common_face_v[0];
            let mut j = -1;
            for z in 0..3 {
                let a = f[fi as usize][z];
                let b = f[fi as usize][(z + 1) % 3];
                if (a == v0 && b == v1) || (a == v1 && b == v0) { j = z as i32; }
            }
            metriko_check!(j != -1, "cut edge {}-{} not found in face {}", v0, v1, fi);
            cuts[fi as usize][j as usize] = 1;
            cuts[tt[fi as usize][j as usize] as usize][tti[fi as usize][j as usize] as usize] = 1;
        }
    }
    cuts
}

//Wrapper of the above with only vertices and faces as mesh input
/// `cut_mesh_with_singularities(V, F, singularities, cuts)`: returns `cuts`
pub fn cut_mesh_with_singularities(v: &[Row3d], f: &[[i32; 3]], singularities: &[i32]) -> Vec<[i32; 3]> {
    let (vf, _vfi) = vertex_triangle_adjacency(v.len(), f);
    let vv = adjacency_list(f, false);
    let (tt, tti) = triangle_triangle_adjacency(f);
    cut_mesh_with_singularities_adj(v, f, &vf, &vv, &tt, &tti, singularities)
}

pub fn compute_seam(f: &FaceRosyField<'_>) -> Vec<bool> {
    let mut seam = vec![false; f.mesh.n_e as usize];
    let mut svids = vec![];
    for v in f.mesh.verts() {
        if f.singular[v.id as usize] != 0 { svids.push(v.id); }
    }
    let face2cut = cut_mesh_with_singularities(&f.mesh.pos, &f.mesh.idx, &svids);
    for i_f in 0..f.mesh.n_f as usize {
        for j in 0..3 {
            if face2cut[i_f][j] != 0 { seam[f.mesh.face2edge[i_f][j] as usize] = true; }
        }
    }
    seam
}
