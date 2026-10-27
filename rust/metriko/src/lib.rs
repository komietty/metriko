//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/lib.h (see rust/PORTED_FROM).

#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]
// The port mirrors the C++ expression by expression (rust/PORTING.md): negated comparisons keep the C++ NaN
// semantics, `(i % 4 + 4) % 4`, `x * -1`, explicit counters and loops are kept as written, and the module layout
// follows include/metriko file by file (hmesh/hmesh.rs, ...).
#![allow(
    clippy::neg_cmp_op_on_partial_ord,
    clippy::type_complexity,
    clippy::assign_op_pattern,
    clippy::explicit_counter_loop,
    clippy::module_inception,
    clippy::neg_multiply,
    clippy::manual_rem_euclid,
    clippy::manual_is_multiple_of,
    clippy::manual_find,
    clippy::while_let_loop,
    clippy::needless_lifetimes
)]

/// `METRIKO_FAIL`: the C++ throws `std::runtime_error`, caught only by `compute_quadrangulation`. Here it panics and
/// `compute_quadrangulation` turns the panic into `Err`.
#[macro_export]
macro_rules! metriko_fail {
    ($($arg:tt)*) => { panic!($($arg)*) };
}

/// `METRIKO_CHECK(cond, ...)`: fails with the message, or with the stringified condition when there is none.
#[macro_export]
macro_rules! metriko_check {
    ($cond:expr $(,)?) => { if !($cond) { $crate::metriko_fail!("{}", stringify!($cond)); } };
    ($cond:expr, $($arg:tt)+) => { if !($cond) { $crate::metriko_fail!($($arg)+); } };
}

pub mod common;
pub mod hmesh;
pub mod nvec;
pub mod igm;
pub mod quantization;
pub mod tmesh;
pub mod tutte;
pub mod qex;
pub mod patching;
pub mod solver;
pub mod igl;

pub use common::typedef::{Complex, Row3d};

use common::utilities::grid_unit;
use hmesh::hmesh::Hmesh;
use nvec::face_rosy_field::{compute_combbed_field, compute_extrinsic_field, compute_seam, FaceRosyField, FieldType};
use igm::parameterization::{compute_cut_mesh, RosyParameterization};
use quantization::quantization::{compute_quantization, validate_quantization};
use tmesh::emesh::Emesh;
use tmesh::emesh_validate::{repair_crossing_tedges, validate_collapse_done};
use tmesh::motorcycle::Mgrph;
use tutte::tutte_cutting::compute_embedding_cut_hmesh;
use tutte::tutte_params::compute_tutte_parameterization;
use tutte::tutte_slim::seamless_slim;
use qex::common::{Qedge, Qface, Qport};
use qex::gen_q_vert::generate_q_vert;
use qex::gen_q_port::{generate_eqvert_qport, generate_fqvert_qport, generate_vqvert_qport};
use qex::gen_q_edge::generate_q_edge;
use qex::gen_q_face::generate_q_faces;
use patching::{extract_quad_mesh, label_quad_patches};

pub const SLIM_MAX_ITER: i32 = 50;
pub const SLIM_REL_TOL: f64 = 1e-4;

/// The collapse loop of `compute_quadrangulation_impl`: zero length thalfs, then zero width tquads, until no zero
/// length thalf is left (at most 100 rounds). Returns the number of rounds run. Factored out of the impl so that the
/// parity tests run this very code.
pub fn collapse_zero_length(tm: &mut Emesh<'_>) -> i32 {
    use tmesh::emesh::Tqaux;
    let mut rounds = 0;
    for i in 0..100 {
        rounds = i + 1;
        // for (Ehalf th0 : tm->thalfs): a copy of each element, taken when it is reached
        let n = tm.thalfs.len();
        for k in 0..n {
            let th0 = tm.thalfs[k];
            if th0.id == -1 || th0.x != 0. { continue; }
            let th1 = *th0.twin(tm);
            if th0.tquad(tm).thids(th0.side(tm)).len() == 1 { continue; }
            if th1.tquad(tm).thids(th1.side(tm)).len() == 1 { continue; }
            tm.collapse_thalf(th0.id);
        }

        // for (const auto& [tqid, _] : tm->live_tquads())
        let n = tm.tquads.len();
        for k in 0..n {
            if tm.tquads[k].id == -1 { continue; }
            let tqid = tm.tquads[k].id;
            let mut a = Tqaux::default();
            if tm.collapse_tquad_prepare(tqid, &mut a) { tm.collapse_tquad_execute(&mut a); }
        }

        if !tm.thalfs.iter().any(|th| th.id != -1 && th.x == 0.) { break; }
    }
    rounds
}

/// The intermediate state of `compute_quadrangulation_impl` (the C++ `RemeshStage`). The C++ returns it owning
/// the meshes; here the T-mesh and the motorcycle graph borrow the mesh, so the stage is lent to a closure instead.
pub struct RemeshStage<'a> {
    pub hmesh: &'a Hmesh,
    pub mgrph: &'a Mgrph<'a>,
    pub emesh: &'a Emesh<'a>,
    /// the mesh with the t-mesh embedded, cut along the seam
    pub hmesh_cut: &'a Hmesh,
    /// per vertex of hmesh_cut: the slim result
    pub uv: &'a [[f64; 2]],
    pub q_ports: &'a [Qport],
    pub q_edges: &'a [Qedge],
    pub q_faces: &'a [Qface],
    pub q_pos: &'a [Row3d],
    pub q_idx: &'a [[i32; 4]],
    pub q_val: &'a [i32],
}

/// The quad mesh: vertex positions, quads (indices into `pos`) and the patch of the T-mesh each quad lies in.
#[derive(Clone, Debug, Default)]
pub struct RemeshResult {
    pub pos: Vec<Row3d>,
    pub idx: Vec<[i32; 4]>,
    pub val: Vec<i32>,
}

/// `compute_quadrangulation_impl(V, F, scale, calign = false, refine = true)`; `out` receives the final stage.
pub fn compute_quadrangulation_impl<R>(
    v: &[Row3d],
    f: &[[i32; 3]],
    scale: f64,
    calign: bool,
    refine: bool,
    out: impl FnOnce(&RemeshStage<'_>) -> R,
) -> R {
    let unit = grid_unit(v, scale);
    let scaled: Vec<Row3d> = v.iter().map(|p| *p / unit).collect();
    let mut hm = Hmesh::new(&scaled, f, false);

    // the field and the integer grid map borrow hm; they are dropped before hm->cfn is written
    let (seam, singular, matching, cfn) = {
        let rawf = FaceRosyField::new(&hm, 4, if calign { FieldType::CurvatureAligned } else { FieldType::Smoothest });
        let seam = compute_seam(&rawf);
        let cutm = compute_cut_mesh(&hm, &seam);
        let cmbf = compute_combbed_field(&rawf, &seam);
        let extf = compute_extrinsic_field(&cmbf, 4);

        let mut rp = RosyParameterization::new(&hm, &cutm, &extf, &cmbf.singular, &cmbf.matching, &seam, 4, scale);
        rp.local_injectivity = true;
        rp.setup();
        rp.integ();

        let mut cfn = vec![Complex::new(0., 0.); hm.n_f as usize * 3];
        for fc in hm.faces() {
            let i = fc.id as usize;
            cfn[i * 3] = Complex::new(rp.cfn[(i, 0)], rp.cfn[(i, 1)]);
            cfn[i * 3 + 1] = Complex::new(rp.cfn[(i, 4)], rp.cfn[(i, 5)]);
            cfn[i * 3 + 2] = Complex::new(rp.cfn[(i, 8)], rp.cfn[(i, 9)]);
        }
        (seam, cmbf.singular.clone(), cmbf.matching.clone(), cfn)
    };
    hm.cfn = cfn;
    let hm = hm;

    let mg = Mgrph::new(&hm, &matching, &singular);
    let mut tm = Emesh::new(&mg);
    let x = compute_quantization(&tm, &mg);
    validate_quantization(&tm, &x);
    tm.set_x(&x);

    collapse_zero_length(&mut tm);
    validate_collapse_done(&tm);

    tm.collapse_tedge_snap(false);
    tm.collapse_tedge_snap(true);
    // for (const auto& [teid, _] : tm->live_tedges())
    for k in 0..tm.tedges.len() {
        if tm.tedges[k].id == -1 { continue; }
        let teid = tm.tedges[k].id;
        tm.collapse_tedge_snap_dedup(teid);
    }

    repair_crossing_tedges(&mut tm);

    //--- cut the original mesh along the collapsed t-mesh ---//
    let (mut hm_emb, seam1, matching1, singular1, hdata) = compute_embedding_cut_hmesh(&hm, &tm, &seam, &matching, &singular);
    let hm_cut = compute_cut_mesh(&hm_emb, &seam1);

    let uv_ebd = compute_tutte_parameterization(&hm_emb, &tm, &seam1, &hdata);
    let uv_cut: Vec<[f64; 2]> = hm_cut.verts().map(|vt| uv_ebd[vt.half().next().crnr().id as usize]).collect();
    let uv = seamless_slim(&hm_emb, &hm_cut, &seam1, &singular1, &uv_cut, SLIM_MAX_ITER, SLIM_REL_TOL);

    // ------ qex on the slim result ------
    hm_emb.cfn = vec![Complex::new(0., 0.); hm_emb.n_f as usize * 3];
    for i in 0..hm_cut.n_f as usize {
        for j in 0..3 {
            let k = hm_cut.idx[i][j] as usize;
            hm_emb.cfn[i * 3 + j] = Complex::new(uv[k][0], uv[k][1]);
        }
    }

    qex::sanitization::sanitization(&mut hm_emb, &matching1, &singular1, 4);

    let mut q_ports: Vec<Qport> = vec![];
    let (vqvs, eqvs, fqvs) = generate_q_vert(&hm_emb);
    generate_vqvert_qport(&hm_emb, &vqvs, &mut q_ports);
    generate_eqvert_qport(&hm_emb, &eqvs, &mut q_ports);
    generate_fqvert_qport(&hm_emb, &fqvs, &mut q_ports);
    let q_edges = generate_q_edge(&hm_emb, &matching1, &mut q_ports);
    let q_faces = generate_q_faces(&q_ports, &q_edges);
    let (mut q_pos, q_idx) = extract_quad_mesh(&hm, &q_ports, &q_edges, &q_faces, refine);
    let q_val = label_quad_patches(&tm, &singular, &q_ports, &q_edges, &q_faces, &q_idx);
    for p in q_pos.iter_mut() { *p *= unit; }

    out(&RemeshStage {
        hmesh: &hm,
        mgrph: &mg,
        emesh: &tm,
        hmesh_cut: &hm_cut,
        uv: &uv,
        q_ports: &q_ports,
        q_edges: &q_edges,
        q_faces: &q_faces,
        q_pos: &q_pos,
        q_idx: &q_idx,
        q_val: &q_val,
    })
}

/// Turns a closed triangle mesh into a quad mesh whose quads have edge length `s` (a fraction of the bounding box
/// diagonal), every quad labelled with the T-mesh patch it lies in. Failures (the C++ exceptions) are returned as
/// `Err` with their message.
pub fn compute_quadrangulation(v: &[Row3d], f: &[[i32; 3]], s: f64) -> Result<RemeshResult, String> {
    let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_quadrangulation_impl(v, f, s, false, true, |q| RemeshResult {
            pos: q.q_pos.to_vec(),
            idx: q.q_idx.to_vec(),
            val: q.q_val.to_vec(),
        })
    }));
    match run {
        Ok(r) => Ok(r),
        Err(e) => {
            if let Some(m) = e.downcast_ref::<String>() { Err(m.clone()) }
            else if let Some(m) = e.downcast_ref::<&str>() { Err(m.to_string()) }
            else { Err("unknown error".to_string()) }
        }
    }
}
