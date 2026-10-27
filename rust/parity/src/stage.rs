//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

//! Rebuilding stage inputs from the dumps and comparing whole structures (Hmesh, Mgrph, Emesh) with them.

use metriko::hmesh::hmesh::Hmesh;
use metriko::hmesh::hmloc::HmLoc;
use metriko::tmesh::emesh::{Edata, Eedge, Ehalf, Emesh, Equad};
use metriko::tmesh::motorcycle::{JunctionType, Mgrph};
use metriko::Complex;
use crate::compare::*;
use crate::golden::Golden;

/// the `Hmesh` of a dumped stage, rebuilt from its `pos` / `idx` (topology only when `geometry` is false)
pub fn hmesh(g: &Golden, stage: &str, geometry: bool) -> Hmesh {
    let pos = g.load(&format!("{stage}.pos")).dvec3s();
    let idx = g.load(&format!("{stage}.idx")).rows_i32::<3>();
    Hmesh::new(&pos, &idx, !geometry)
}

/// every array of a dumped Hmesh: topology exactly, geometry bitwise
pub fn check_hmesh(case: &str, hm: &Hmesh, g: &Golden, stage: &str, geometry: bool) -> Result<(), String> {
    let n = |s: &str| format!("{case} {stage}.{s}");
    let l = |s: &str| g.load(&format!("{stage}.{s}"));
    check_bits_f64(&n("pos"), &flat3(&hm.pos), l("pos").f64s(), 3)?;
    check_eq_i32(&n("idx"), &flat(&hm.idx), l("idx").i32s(), 3)?;
    check_eq_i32(&n("edge2vert"), &flat(&hm.edge2vert), l("edge2vert").i32s(), 2)?;
    check_eq_i32(&n("face2edge"), &flat(&hm.face2edge), l("face2edge").i32s(), 3)?;
    check_eq_i32(&n("edge2face"), &flat(&hm.edge2face), l("edge2face").i32s(), 2)?;
    for (s, v) in [
        ("vert2half", &hm.vert2half), ("edge2half", &hm.edge2half), ("face2half", &hm.face2half),
        ("loop2half", &hm.loop2half), ("crnr2half", &hm.crnr2half), ("next", &hm.next), ("prev", &hm.prev),
        ("twin", &hm.twin), ("tail", &hm.tail), ("head", &hm.head), ("edge", &hm.edge), ("face", &hm.face),
        ("crnr", &hm.crnr),
    ] {
        check_eq_i32(&n(s), v, l(s).i32s(), 1)?;
    }
    check_eq_bool(&n("isBV"), &hm.is_bv, &l("isBV").bools())?;
    if !geometry { return Ok(()); }
    for (s, v) in [
        ("vertNormal", &hm.vert_normal), ("vertBasisX", &hm.vert_basis_x), ("vertBasisY", &hm.vert_basis_y),
        ("faceNormal", &hm.face_normal), ("faceBasisX", &hm.face_basis_x), ("faceBasisY", &hm.face_basis_y),
        ("baryCenter", &hm.bary_center),
    ] {
        check_bits_f64(&n(s), &flat3(v), l(s).f64s(), 3)?;
    }
    for (s, v) in [
        ("faceArea", &hm.face_area), ("baryDualArea", &hm.bary_dual_area), ("circDualArea", &hm.circ_dual_area),
        ("angleDefect", &hm.angle_defect), ("halfCotan", &hm.half_cotan), ("dihedralArg", &hm.dihedral_arg),
        ("heArgOnVert", &hm.he_arg_on_vert), ("heArgOnFace", &hm.he_arg_on_face), ("edgeCotan", &hm.edge_cotan),
    ] {
        check_bits_f64(&n(s), v, l(s).f64s(), 1)?;
    }
    Ok(())
}

/// the dumped `HmLoc` encoding: `(kind, id)` and `(a, b)` per FORMAT.md
pub fn encode_loc(l: &HmLoc) -> ([i32; 2], [f64; 2]) {
    match *l {
        HmLoc::OnV { id } => ([0, id], [0., 0.]),
        HmLoc::OnE { id, r } => ([1, id], [r, 0.]),
        HmLoc::OnF { id, xy } => ([2, id], [xy.re, xy.im]),
        HmLoc::OnH { id, r } => ([3, id], [r, 0.]),
        HmLoc::OnC { id } => ([4, id], [0., 0.]),
        HmLoc::OnP { id, uv } => ([5, id], [uv.re, uv.im]),
    }
}

pub fn decode_loc(ki: [i32; 2], val: [f64; 2]) -> HmLoc {
    match ki[0] {
        0 => HmLoc::OnV { id: ki[1] },
        1 => HmLoc::OnE { id: ki[1], r: val[0] },
        2 => HmLoc::OnF { id: ki[1], xy: Complex::new(val[0], val[1]) },
        3 => HmLoc::OnH { id: ki[1], r: val[0] },
        4 => HmLoc::OnC { id: ki[1] },
        5 => HmLoc::OnP { id: ki[1], uv: Complex::new(val[0], val[1]) },
        k => panic!("bad HmLoc kind {k}"),
    }
}

pub fn locs(g: &Golden, key: &str) -> Vec<HmLoc> {
    let ki = g.load(&format!("{key}_ki")).rows_i32::<2>();
    let val = g.load(&format!("{key}_val")).rows_f64::<2>();
    ki.iter().zip(&val).map(|(k, v)| decode_loc(*k, *v)).collect()
}

/// locations: kinds and ids exactly, values bitwise
pub fn check_locs(what: &str, rust: &[HmLoc], g: &Golden, key: &str) -> Result<(), String> {
    let (ki, val): (Vec<[i32; 2]>, Vec<[f64; 2]>) = rust.iter().map(encode_loc).unzip();
    check_eq_i32(&format!("{what} kind/id"), &flat(&ki), g.load(&format!("{key}_ki")).i32s(), 2)?;
    check_bits_f64(&format!("{what} value"), &flat(&val), g.load(&format!("{key}_val")).f64s(), 2)
}

fn ragged<T: Copy>(offsets: &[i32], rows: &[T], cols: usize) -> Vec<Vec<T>> {
    (0..offsets.len() - 1).map(|i| rows[offsets[i] as usize * cols..offsets[i + 1] as usize * cols].to_vec()).collect()
}

fn jt_code(jt: JunctionType) -> i32 {
    match jt { JunctionType::None => 0, JunctionType::F => 1, JunctionType::T => 2, JunctionType::C => 3 }
}

pub fn check_mgrph(case: &str, mg: &Mgrph<'_>, g: &Golden, stage: &str) -> Result<(), String> {
    let n = |s: &str| format!("{case} {stage}.{s}");
    let l = |s: &str| g.load(&format!("{stage}.{s}"));
    let dirs: Vec<Complex> = mg.mports.iter().map(|p| p.dir).collect();
    check_bits_c64(&n("mports_dir"), &dirs, l("mports_dir").c64s(), 1)?;
    let ids: Vec<i32> = mg.mports.iter().flat_map(|p| [p.crnr_id, p.this_id, p.next_id, p.prev_id]).collect();
    check_eq_i32(&n("mports_ids"), &ids, l("mports_ids").i32s(), 4)?;
    let jt: Vec<i32> = mg.mnodes.iter().map(|m| jt_code(m.jt)).collect();
    check_eq_i32(&n("mnodes_jt"), &jt, l("mnodes_jt").i32s(), 1)?;
    let nl: Vec<HmLoc> = mg.mnodes.iter().map(|m| m.loc).collect();
    check_locs(&n("mnodes_loc"), &nl, g, &format!("{stage}.mnodes_loc"))?;
    let adj_cpp = ragged(l("mnodes_adj_offsets").i32s(), l("mnodes_adj").i32s(), 2);
    let adj_rs: Vec<Vec<i32>> = mg.mnodes.iter().map(|m| flat(&m.adj)).collect();
    check_eq_i32(&n("mnodes_adj count"), &[adj_rs.len() as i32], &[adj_cpp.len() as i32], 1)?;
    for (i, (a, b)) in adj_rs.iter().zip(&adj_cpp).enumerate() { check_eq_i32(&n(&format!("mnodes[{i}].adj")), a, b, 2)?; }
    let cid: Vec<i32> = mg.mcurvs.iter().map(|c| c.id).collect();
    check_eq_i32(&n("mcurvs_id"), &cid, l("mcurvs_id").i32s(), 1)?;
    let sg_cpp = ragged(l("mcurvs_sgmt_offsets").i32s(), l("mcurvs_sgmts").i32s(), 5);
    for (i, (c, b)) in mg.mcurvs.iter().zip(&sg_cpp).enumerate() {
        let a: Vec<i32> = c.sgmts.iter().flat_map(|s| [s.fr_nid, s.to_nid, s.this_id, s.curv_id, s.face_id]).collect();
        check_eq_i32(&n(&format!("mcurvs[{i}].sgmts")), &a, b, 5)?;
    }
    let bl: Vec<HmLoc> = mg.mcurvs.iter().map(|c| c.buff.loc).collect();
    check_locs(&n("mcurvs_buff_loc"), &bl, g, &format!("{stage}.mcurvs_buff_loc"))?;
    let bd: Vec<Complex> = mg.mcurvs.iter().map(|c| c.buff.dir).collect();
    check_bits_c64(&n("mcurvs_buff_dir"), &bd, l("mcurvs_buff_dir").c64s(), 1)?;
    let be: Vec<bool> = mg.mcurvs.iter().map(|c| c.buff.end).collect();
    check_eq_bool(&n("mcurvs_buff_end"), &be, &l("mcurvs_buff_end").bools())
}

/// an Emesh: structure exactly; `x`, `r` and node positions bitwise when `bits`, else within `tol`
pub fn check_emesh(case: &str, tm: &Emesh<'_>, g: &Golden, stage: &str, tol: Option<Tol>) -> Result<(), String> {
    let n = |s: &str| format!("{case} {stage}.{s}");
    let l = |s: &str| g.load(&format!("{stage}.{s}"));
    let (ki, val): (Vec<[i32; 2]>, Vec<[f64; 2]>) = tm.tnodes.iter().map(encode_loc).unzip();
    check_eq_i32(&n("tnodes kind/id"), &flat(&ki), l("tnodes_ki").i32s(), 2)?;
    match tol {
        None => check_bits_f64(&n("tnodes value"), &flat(&val), l("tnodes_val").f64s(), 2)?,
        Some(t) => check_close_f64(&n("tnodes value"), &flat(&val), l("tnodes_val").f64s(), 2, t)?,
    }
    check_eq_bool(&n("fixed"), &tm.fixed, &l("fixed").bools())?;
    let te: Vec<i32> = tm.tedges.iter().map(|e| e.id).collect();
    check_eq_i32(&n("tedges_id"), &te, l("tedges_id").i32s(), 1)?;
    let nids_cpp = ragged(l("tedges_nids_offsets").i32s(), l("tedges_nids").i32s(), 1);
    for (i, (e, b)) in tm.tedges.iter().zip(&nids_cpp).enumerate() { check_eq_i32(&n(&format!("tedges[{i}].nids")), &e.nids, b, 1)?; }
    let th: Vec<i32> = tm.thalfs.iter().flat_map(|h| [h.id, h.twid, h.teid, h.tqid, h.cano as i32, h.bgn as i32, h.end as i32]).collect();
    check_eq_i32(&n("thalfs_ids"), &th, l("thalfs_ids").i32s(), 7)?;
    let xr: Vec<f64> = tm.thalfs.iter().flat_map(|h| [h.x, h.r]).collect();
    match tol {
        None => check_bits_f64(&n("thalfs_xr"), &xr, l("thalfs_xr").f64s(), 2)?,
        Some(t) => check_close_f64(&n("thalfs_xr"), &xr, l("thalfs_xr").f64s(), 2, t)?,
    }
    let tq: Vec<i32> = tm.tquads.iter().map(|q| q.id).collect();
    check_eq_i32(&n("tquads_id"), &tq, l("tquads_id").i32s(), 1)?;
    let data_cpp = ragged(l("tquads_data_offsets").i32s(), l("tquads_data").i32s(), 2);
    for (i, (q, b)) in tm.tquads.iter().zip(&data_cpp).enumerate() {
        let a: Vec<i32> = q.data.iter().flat_map(|d| [d.thid, d.side]).collect();
        check_eq_i32(&n(&format!("tquads[{i}].data")), &a, b, 2)?;
    }
    Ok(())
}

/// an Emesh holding the dumped arrays of `stage` (any stage after 11), on the mesh `hm` it was built on
pub fn emesh<'a>(hm: &'a metriko::hmesh::hmesh::Hmesh, g: &Golden, stage: &str) -> Emesh<'a> {
    let l = |s: &str| g.load(&format!("{stage}.{s}"));
    let mut tm = Emesh::empty(hm);
    tm.tnodes = locs(g, &format!("{stage}.tnodes"));
    tm.fixed = l("fixed").bools();
    let nids = ragged(l("tedges_nids_offsets").i32s(), l("tedges_nids").i32s(), 1);
    tm.tedges = l("tedges_id").i32s().iter().zip(nids).map(|(&id, nids)| Eedge { id, nids }).collect();
    let ids = l("thalfs_ids").rows_i32::<7>();
    let xr = l("thalfs_xr").rows_f64::<2>();
    tm.thalfs = ids.iter().zip(&xr).map(|(a, b)| Ehalf {
        id: a[0], twid: a[1], teid: a[2], tqid: a[3], cano: a[4] != 0, bgn: a[5] != 0, end: a[6] != 0, x: b[0], r: b[1],
    }).collect();
    let data = ragged(l("tquads_data_offsets").i32s(), l("tquads_data").i32s(), 2);
    tm.tquads = l("tquads_id").i32s().iter().zip(data).map(|(&id, d)| Equad {
        id, data: d.chunks(2).map(|p| Edata { thid: p[0], side: p[1] }).collect(),
    }).collect();
    tm
}
