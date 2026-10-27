//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/hmesh/hmloc.h (see rust/PORTED_FROM).

use crate::common::typedef::*;
use crate::common::utilities::find_strict_intersection;
use crate::common::predicates::is_inside_segment;
use crate::hmesh::hmesh::*;

/// `HmLoc = std::variant<HmLocOnV, HmLocOnE, HmLocOnF, HmLocOnH, HmLocOnC, HmLocOnP>`, alternatives in the same
/// order. Equality is the variant's: same alternative and equal members.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HmLoc {
    OnV { id: i32 },
    OnE { id: i32, r: f64 },
    OnF { id: i32, xy: Complex },
    OnH { id: i32, r: f64 },
    OnC { id: i32 },
    /// parameter space coord in face
    OnP { id: i32, uv: Complex },
}

impl Default for HmLoc {
    /// a value-initialized variant holds its first alternative
    fn default() -> Self { HmLoc::OnV { id: 0 } }
}

impl HmLoc {
    /// `std::variant::index()`
    pub fn index(&self) -> i32 {
        match self {
            HmLoc::OnV { .. } => 0,
            HmLoc::OnE { .. } => 1,
            HmLoc::OnF { .. } => 2,
            HmLoc::OnH { .. } => 3,
            HmLoc::OnC { .. } => 4,
            HmLoc::OnP { .. } => 5,
        }
    }
}

pub fn get_ptloc_pos(hm: &Hmesh, loc: &HmLoc) -> Row3d {
    match *loc {
        HmLoc::OnV { id } => hm.vert(id).pos(),
        HmLoc::OnC { id } => hm.crnr(id).vert().pos(),
        HmLoc::OnE { id, r } => hm.edge(id).lerp(r),
        HmLoc::OnH { id, r } => hm.half(id).lerp(r),
        HmLoc::OnF { id, xy } => hm.face(id).to_world(xy),
        HmLoc::OnP { .. } => metriko_fail!("no impl"),
    }
}

pub fn get_ptloc_dif(hm: &Hmesh, l0: &HmLoc, l1: &HmLoc) -> Row3d {
    get_ptloc_pos(hm, l1) - get_ptloc_pos(hm, l0)
}

pub fn get_ptloc_nml(hm: &Hmesh, loc: &HmLoc) -> Row3d {
    match *loc {
        HmLoc::OnV { id } => hm.vert(id).normal(),
        HmLoc::OnF { id, .. } => hm.face(id).normal(),
        HmLoc::OnE { id, .. } => hm.edge(id).nml(),
        HmLoc::OnH { id, .. } => hm.half(id).nml(),
        _ => metriko_fail!("no impl"),
    }
}

pub fn get_ptloc_verts(hm: &Hmesh, loc: &HmLoc) -> Vec<i32> {
    match *loc {
        HmLoc::OnV { id } => vec![id],
        HmLoc::OnE { id, .. } => { let e = hm.edge(id); vec![e.vert0().id, e.vert1().id] }
        HmLoc::OnH { id, .. } => { let h = hm.half(id); vec![h.tail().id, h.head().id] }
        HmLoc::OnF { id, .. } => { let [a, b, c] = hm.face(id).verts(); vec![a.id, b.id, c.id] }
        _ => vec![],
    }
}

pub fn get_ptloc_edges(hm: &Hmesh, loc: &HmLoc) -> Vec<i32> {
    let mut out = vec![];
    match *loc {
        HmLoc::OnV { id } => { for h in hm.vert(id).adj_halfs() { out.push(h.edge().id); } }
        HmLoc::OnE { id, .. } => { out.push(hm.edge(id).id); }
        HmLoc::OnH { id, .. } => { out.push(hm.half(id).edge().id); }
        _ => {}
    }
    out
}

pub fn get_ptloc_faces(hm: &Hmesh, loc: &HmLoc) -> Vec<i32> {
    let mut out = vec![];
    match *loc {
        HmLoc::OnV { id } => { for h in hm.vert(id).adj_halfs() { out.push(h.face().id); } }
        HmLoc::OnE { id, .. } => { out.push(hm.edge(id).face0().id); out.push(hm.edge(id).face1().id); }
        HmLoc::OnH { id, .. } => { let h = hm.half(id); out.push(h.face().id); out.push(h.twin().face().id); }
        HmLoc::OnF { id, .. } => { out.push(id); }
        _ => {}
    }
    out.retain(|&x| x != -1);
    out
}

pub fn try_get_face<'a>(hm: &'a Hmesh, a: &HmLoc, b: &HmLoc) -> Option<Face<'a>> {
    let fa = get_ptloc_faces(hm, a);
    let fb = get_ptloc_faces(hm, b);
    for &x in &fa {
        for &y in &fb {
            if x == y { return Some(hm.face(x)); }
        }
    }
    None
}

pub fn try_get_edge<'a>(hm: &'a Hmesh, a: &HmLoc, b: &HmLoc) -> Option<Edge<'a>> {
    let ea = get_ptloc_edges(hm, a);
    let eb = get_ptloc_edges(hm, b);
    for &x in &ea {
        for &y in &eb {
            if x == y { return Some(hm.edge(x)); }
        }
    }
    None
}

pub fn try_get_half_ratio<'a>(hm: &'a Hmesh, loc: &HmLoc) -> Option<(Half<'a>, f64)> {
    match *loc {
        HmLoc::OnE { id, r } => {
            let e = hm.edge(id);
            let h = e.half();
            Some((h, 1. - if h.is_canonical() { r } else { 1. - r }))
        }
        HmLoc::OnH { id, r } => Some((hm.half(id), 1. - r)),
        _ => None,
    }
}

pub fn try_get_edge_ratio<'a>(hm: &'a Hmesh, loc: &HmLoc) -> Option<(Edge<'a>, f64)> {
    match *loc {
        HmLoc::OnE { id, r } => Some((hm.edge(id), r)),
        HmLoc::OnH { id, r } => { let h = hm.half(id); Some((h.edge(), if h.is_canonical() { r } else { 1. - r })) }
        _ => None,
    }
}

pub fn is_in_face(f: Face<'_>, loc: &HmLoc) -> bool {
    match *loc {
        HmLoc::OnV { id } => { for h in f.adj_halfs() { if h.tail().id == id { return true; } } false }
        HmLoc::OnE { id, .. } => { for h in f.adj_halfs() { if h.edge().id == id { return true; } } false }
        HmLoc::OnH { id, .. } => { for h in f.adj_halfs() { if h.id == id { return true; } } false }
        HmLoc::OnF { id, .. } => id == f.id,
        _ => metriko_fail!("no impl"),
    }
}

/// The overload `find_strict_intersection(Face, HmLoc, HmLoc, HmLoc, HmLoc)`.
pub fn find_strict_intersection_in_face(f: Face<'_>, la: &HmLoc, lb: &HmLoc, lc: &HmLoc, ld: &HmLoc) -> bool {
    if !is_in_face(f, la) || !is_in_face(f, lb) || !is_in_face(f, lc) || !is_in_face(f, ld) { return false; }
    find_strict_intersection(
        f.to_local(get_ptloc_pos(f.m, la)),
        f.to_local(get_ptloc_pos(f.m, lb)),
        f.to_local(get_ptloc_pos(f.m, lc)),
        f.to_local(get_ptloc_pos(f.m, ld)),
        EPS,
    )
}

/// intersection except edge point and edge point
pub fn has_contact(f: Face<'_>, la: &HmLoc, lb: &HmLoc, lc: &HmLoc, ld: &HmLoc) -> bool {
    if !is_in_face(f, la) || !is_in_face(f, lb) || !is_in_face(f, lc) || !is_in_face(f, ld) { return false; }
    let a = f.to_local(get_ptloc_pos(f.m, la));
    let b = f.to_local(get_ptloc_pos(f.m, lb));
    let c = f.to_local(get_ptloc_pos(f.m, lc));
    let d = f.to_local(get_ptloc_pos(f.m, ld));
    find_strict_intersection(a, b, c, d, EPS)
        || is_inside_segment(a, b, c)
        || is_inside_segment(a, b, d)
        || is_inside_segment(c, d, a)
        || is_inside_segment(c, d, b)
}

pub fn is_in_star(v: Vert<'_>, loc: &HmLoc) -> bool {
    v.adj_halfs().any(|h| is_in_face(h.face(), loc))
}

pub fn loc_str(l: &HmLoc) -> String {
    match *l {
        HmLoc::OnV { id } => format!("V(id={})", id),
        HmLoc::OnC { id } => format!("C(id={})", id),
        HmLoc::OnE { id, r } => format!("E(id={}, r={})", id, r),
        HmLoc::OnH { id, r } => format!("H(id={}, r={})", id, r),
        HmLoc::OnF { id, xy } => format!("F(id={}, xy=({},{}))", id, xy.re, xy.im),
        HmLoc::OnP { id, uv } => format!("P(id={}, uv=({},{}))", id, uv.re, uv.im),
    }
}

pub fn get_chart_face<'a>(hm: &'a Hmesh, loc: &HmLoc) -> Face<'a> {
    match *loc {
        HmLoc::OnC { id } => hm.crnr(id).face(),
        HmLoc::OnH { id, .. } => hm.half(id).face(),
        HmLoc::OnP { id, .. } => hm.face(id),
        _ => metriko_fail!("{} has no chart", loc_str(loc)),
    }
}

pub fn to_chart_free(hm: &Hmesh, loc: &HmLoc) -> HmLoc {
    match *loc {
        HmLoc::OnC { id } => HmLoc::OnV { id: hm.crnr(id).vert().id },
        HmLoc::OnH { .. } => { let (e, r) = try_get_edge_ratio(hm, loc).expect("bad optional access"); HmLoc::OnE { id: e.id, r } }
        l => l,
    }
}

pub fn get_chart_uv(hm: &Hmesh, loc: &HmLoc) -> Complex {
    match *loc {
        HmLoc::OnC { id } => hm.crnr(id).uv(),
        HmLoc::OnH { id, r } => hm.half(id).uv(r),
        HmLoc::OnP { uv, .. } => uv,
        _ => metriko_fail!("{} has no chart", loc_str(loc)),
    }
}
