//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tutte/tutte.h (see rust/PORTED_FROM).

/// A halfedge of the cut mesh lying on a tedge, with its parameter range along the thalf.
/// The C++ holds the `Half` itself; `compute_embedding_cut_hmesh` creates the mesh it refers to, so here it is the
/// half id in that mesh (`-1` for `Half()`).
#[derive(Clone, Copy, Debug)]
pub struct HalfData {
    pub half: i32,  //
    pub v0: f64,    // parameter value of tail
    pub v1: f64,    // parameter value of head
    pub thid: i32,  //
    pub tqid: i32,  //
    pub order: i32, // order inside thalf
}

impl Default for HalfData {
    fn default() -> Self { HalfData { half: -1, v0: 0., v1: 0., thid: -1, tqid: -1, order: -1 } }
}

impl HalfData {
    /// `operator<`: `std::tie(tqid, thid, order) < std::tie(rhs.tqid, rhs.thid, rhs.order)`
    pub fn lt(&self, rhs: &HalfData) -> bool { (self.tqid, self.thid, self.order) < (rhs.tqid, rhs.thid, rhs.order) }
}
