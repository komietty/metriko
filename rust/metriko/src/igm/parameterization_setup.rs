//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/igm/parameterization_setup.ipp (see rust/PORTED_FROM).

use std::collections::BTreeSet;
use faer::sparse::Triplet;
use crate::igm::parameterization::*;
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets, nonzeros};

/// Eigen's `SparseMatrix::prune(reference, epsilon = dummy_precision)`: drops the entries with
/// `|v| <= |reference| * epsilon` (`dummy_precision` of double is 1e-12). A NaN is kept, as in Eigen.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn prune(a: &SprsD, reference: f64) -> SprsD {
    let eps = 1e-12;
    let t: Vec<TripD> = nonzeros(a.as_ref()).filter(|&(_, _, v)| !(v.abs() <= reference.abs() * eps)).map(|(i, j, v)| Triplet::new(i, j, v)).collect();
    from_triplets(a.nrows(), a.ncols(), &t)
}

impl<'a> RosyParameterization<'a> {
    pub(crate) fn compute_he2matching(&mut self) {
        let raw = self.raw;
        let nb = self.n_big;
        self.he2matching = vec![0; raw.n_h as usize];
        for h in raw.halfs() {
            let m = (if h.is_canonical() { -1 } else { 1 }) * self.matching[h.edge().id as usize]; // better not inversed
            self.he2matching[h.id as usize] = (m % nb + nb) % nb;
        }
    }

    pub(crate) fn compute_he2transidx(&mut self) {
        let raw = self.raw;
        self.he2transidx = vec![32767; raw.n_h as usize];
        let mut valence = vec![0f64; raw.n_v as usize];
        let mut claimed = vec![false; raw.n_e as usize];

        for e in raw.edges() {
            if self.seam[e.id as usize] {
                valence[e.vert0().id as usize] += 1.;
                valence[e.vert1().id as usize] += 1.;
            }
        }

        let mut tid = 1;

        for v in raw.verts() {
            if (valence[v.id as usize] == 2. && !self.is_inside_singular(v)) || valence[v.id as usize] == 0. { continue; }
            for c_h in v.adj_halfs_from(self.find_first_bndr_he(v), true) {
                let c_e = c_h.edge();
                if c_e.is_boundary() || !self.seam[c_e.id as usize] || claimed[c_e.id as usize] { continue; }
                self.he2transidx[c_h.id as usize] = tid;
                self.he2transidx[c_h.twin().id as usize] = -tid;
                claimed[c_e.id as usize] = true;

                // traces seam until reaching to singular or boundary vertex
                let mut jv = c_h.head();
                while valence[jv.id as usize] == 2. && !self.is_inside_singular(jv) && !jv.is_boundary() {
                    // the C++ default-constructed `Half nh` has id -1 when no unclaimed seam half is found
                    let mut nh = raw.half(-1);
                    for h in jv.adj_halfs() {
                        if self.seam[h.edge().id as usize] && !claimed[h.edge().id as usize] { nh = h; break; }
                    }

                    self.he2transidx[nh.id as usize] = tid;
                    self.he2transidx[nh.twin().id as usize] = -tid;
                    claimed[nh.edge().id as usize] = true;
                    jv = nh.head();
                }
                tid += 1;
            }
        }
        self.n_t = tid - 1;
    }

    pub fn setup(&mut self) {
        let raw = self.raw;
        let cut = self.cut;
        let nb = self.n_big as usize;
        let n = self.n as usize;

        // here we compute a permutation matrix
        let mut const_parm_mats: Vec<MatI> = vec![MatI::zero(nb, nb); nb];
        for k in 0..nb {
            for i in 0..nb { const_parm_mats[k].set((i + k) % nb, i, 1); }
        }

        let mut v_t: Vec<TripD> = vec![];
        let mut c_t: Vec<TripD> = vec![];
        // forming the constraints and the singularity positions
        let mut curr_constraint = 0i32;
        // this loop set up the transitions (vector field matching) across the cuts
        for v in raw.verts() {
            // 1: The initial corner gets the identity without any transition
            let mut perm_mats: Vec<MatI> = vec![];
            let mut perm_idcs: Vec<i32> = vec![];
            perm_mats.push(MatI::identity(nb, nb));
            perm_idcs.push(v.id);
            let mut i_vcut = -1;

            // remember uv = Rot * val + Transition
            let first = if v.is_boundary() { self.find_first_bndr_he(v) } else { self.find_first_seam_he(v) };
            for h in v.adj_halfs_from(first, true) {
                if h.is_boundary() { break; } // last boundary
                let j_vcut = cut.idx[h.face().id as usize][(h.next().crnr().id % 3) as usize];

                if j_vcut != i_vcut {
                    i_vcut = j_vcut;
                    for i in 0..perm_idcs.len() { assign_block(&mut v_t, &perm_mats[i], self.n_big * i_vcut, self.n_big * perm_idcs[i]); }
                }

                let hn = h.prev().twin();
                if !hn.is_boundary() && self.seam[hn.edge().id as usize] {
                    let p = &const_parm_mats[self.he2matching[hn.id as usize] as usize];
                    let t = self.he2transidx[hn.id as usize];
                    if t > 0 {
                        for m in perm_mats.iter_mut() { *m = p.mul(m); }
                        perm_mats.push(MatI::identity(nb, nb));
                        perm_idcs.push(raw.n_v + t - 1);
                    } else {
                        perm_mats.push(MatI::identity(nb, nb).neg());
                        perm_idcs.push(raw.n_v - t - 1);
                        for m in perm_mats.iter_mut() { *m = p.mul(m); }
                    }
                }
            }

            // 2: cleaning parmMats and permIdcs to see if there is a constraint or reveal singularity-from-transition
            if !v.is_boundary() {
                let temp: BTreeSet<i32> = perm_idcs.iter().copied().collect();
                let idcs: Vec<i32> = temp.into_iter().collect();
                let mut mats: Vec<MatI> = vec![MatI::zero(0, 0); idcs.len()];

                for j in 0..idcs.len() {
                    mats[j] = MatI::zero(nb, nb);
                    for k in 0..perm_idcs.len() {
                        if idcs[j] == perm_idcs[k] { mats[j].add_assign(&perm_mats[k]); }
                    }
                    if idcs[j] == v.id { mats[j].sub_assign(&MatI::identity(nb, nb)); }
                }

                if mats.iter().any(|m| m.abs_max() != 0) {
                    for j in 0..mats.len() {
                        assign_block(&mut c_t, &mats[j], self.n_big * curr_constraint, self.n_big * idcs[j]);
                    }
                    curr_constraint += 1;
                }
            }
        }

        self.vtrans2cut = from_triplets(nb * cut.n_v as usize, nb * self.n_r as usize, &v_t);
        self.vtrans2cut = prune(&self.vtrans2cut, 1e-3);

        self.constraint = from_triplets(nb * curr_constraint as usize, nb * self.n_r as usize, &c_t);
        self.constraint = prune(&self.constraint, 1e-3);

        // filtering out barycentric symmetry, including sign symmetry.
        // The parameterization should always only include n dof for the surface
        // Warning: this assumes n divides N!
        // integer variables are per single "d" packet, and the rounding is done for the N functions with projection over linRed
        let mut buff: Vec<TripD> = vec![];
        let mut i = 0;
        while i < nb as i32 * self.n_r {
            assign_block(&mut buff, &self.lreductor, i, i * self.n / self.n_big);
            i += self.n_big;
        }
        self.uncompress = from_triplets(nb * self.n_r as usize, n * self.n_r as usize, &buff);

        //----- indices to be integer (used for rounding on seams) -----
        let sv = raw.verts().find(|&v| self.is_inside_singular(v));
        let i_v = match sv { None => 0, Some(v) => v.id };
        self.fixed_idcs = lin_spaced(self.n, self.n * i_v, self.n * i_v + self.n - 1);
        self.integer_idcs = lin_spaced(self.n * self.n_t, self.n * raw.n_v, self.n * (raw.n_v + self.n_t) - 1);

        //----- indices of singular (used for rounding on singulars before rounding on seams) -----
        self.singular_idcs = vec![0; (self.n * self.n_s) as usize];
        let mut c = 0;
        let svs: Vec<i32> = raw.verts().filter(|&e| self.is_inside_singular(e)).map(|v| v.id).collect();
        for vid in svs {
            for j in 0..self.n { self.singular_idcs[c] = self.n * vid + j; c += 1; }
        }
    }
}

/// Eigen's `VecXi::LinSpaced(size, low, high)` for the cases used here, where `high - low + 1 == size` (consecutive
/// integers) or `size <= 1` (a single entry `high`).
fn lin_spaced(size: i32, low: i32, high: i32) -> Vec<i32> {
    if size <= 0 { return vec![]; }
    if size == 1 { return vec![high]; }
    assert!(high - low + 1 == size, "LinSpaced: only the consecutive case is ported");
    (low..=high).collect()
}
