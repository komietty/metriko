//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Stage-by-stage comparison of the Rust port against the dumps of the C++ metriko (tools/dump/FORMAT.md).
//!
//! - [`npy`]: reader of the dumped arrays
//! - [`golden`]: the cases of testdata/cases.json and their dump directories
//! - [`compare`]: equality / tolerance checks with readable mismatch reports
//! - [`mesh`]: the input mesh as metriko_web's cli reads it (obj)
//! - [`summary`]: the summary lines and hash of metriko_web's cli

// mesh.rs mirrors the C++ obj loader of metriko_web statement by statement
#![allow(clippy::needless_lifetimes, clippy::type_complexity, clippy::needless_range_loop)]

pub mod compare;
pub mod golden;
pub mod mesh;
pub mod npy;
pub mod summary;
pub mod stage;

pub use compare::*;
pub use golden::{cases, golden_cases, Case, Golden};
pub use mesh::{grid_unit, load_mesh, InputMesh};
pub use npy::Npy;
pub use summary::Summary;
