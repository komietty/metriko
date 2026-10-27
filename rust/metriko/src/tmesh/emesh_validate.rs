//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tmesh/emesh_validate.h (see rust/PORTED_FROM).

use crate::common::stdcxx::LibcxxIntMap;
use std::collections::BTreeSet;
use crate::common::stdcxx;
use crate::common::predicates::is_inside_segment;
use crate::common::typedef::*;
use crate::common::utilities::find_strict_intersection;
use crate::hmesh::hmloc::*;
use crate::tmesh::emesh::*;

// a contact between two tedges: either a strict segment crossing inside a face
// or a "touch" (a node of one tedge lying in the interior of the other's segment)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TeContact {
    pub fid: i32,
    pub te_seg: i32, // owner of the crossed / touched segment
    pub te_ndp: i32, // the other tedge (whose node or segment intrudes)
    pub cross: bool, // true: strict crossing, false: touch
}

#[derive(Clone, Copy, Debug)]
struct S { teid: i32, n0: i32, n1: i32 }

// tedges must not touch each other except at shared junction nodes
pub fn find_tedge_contacts(tm: &Emesh<'_>) -> Vec<TeContact> {
    let mut per_face: LibcxxIntMap<Vec<S>> = LibcxxIntMap::new();

    for e in tm.live_tedges() {
        let nids = &e.nids;
        let mut k = 0;
        while k + 1 < nids.len() {
            for i in get_ptloc_faces(tm.hm, &tm.tnodes[nids[k] as usize]) {
                if is_in_face(tm.hm.face(i), &tm.tnodes[nids[k + 1] as usize]) {
                    per_face.entry(i).push(S { teid: e.id, n0: nids[k], n1: nids[k + 1] });
                }
            }
            k += 1;
        }
    }

    let mut res: Vec<TeContact> = vec![];
    let mut seen: BTreeSet<(i32, i32)> = BTreeSet::new();
    let mut report = |fid: i32, te_seg: i32, te_ndp: i32, cross: bool| {
        if seen.insert(stdcxx::minmax(te_seg, te_ndp)) { res.push(TeContact { fid, te_seg, te_ndp, cross }); }
    };

    for (fid, segs) in per_face.iter() {
        let f = tm.hm.face(fid);
        for i in 0..segs.len() {
            for j in i + 1..segs.len() {
                let s = &segs[i];
                let t = &segs[j];
                if s.teid == t.teid { continue; }

                let a0 = f.to_local(get_ptloc_pos(tm.hm, &tm.tnodes[s.n0 as usize]));
                let a1 = f.to_local(get_ptloc_pos(tm.hm, &tm.tnodes[s.n1 as usize]));
                let b0 = f.to_local(get_ptloc_pos(tm.hm, &tm.tnodes[t.n0 as usize]));
                let b1 = f.to_local(get_ptloc_pos(tm.hm, &tm.tnodes[t.n1 as usize]));

                let shared = s.n0 == t.n0 || s.n0 == t.n1 || s.n1 == t.n0 || s.n1 == t.n1;
                if !shared && find_strict_intersection(a0, a1, b0, b1, EPS) { report(fid, s.teid, t.teid, true); continue; }
                // two distinct nodes at one point (a collapse can leave them): the tedges meet there without a shared node
                let coincide = |n: i32, m: i32, p: Complex, q: Complex| n != m && (p - q).norm() < EPS;
                if coincide(s.n0, t.n0, a0, b0) || coincide(s.n0, t.n1, a0, b1) || coincide(s.n1, t.n0, a1, b0) || coincide(s.n1, t.n1, a1, b1) { report(fid, s.teid, t.teid, false); continue; }
                if is_inside_segment(a0, a1, b0) || is_inside_segment(a0, a1, b1) { report(fid, s.teid, t.teid, false); continue; } // t's endpoint on s
                if is_inside_segment(b0, b1, a0) || is_inside_segment(b0, b1, a1) { report(fid, t.teid, s.teid, false); }           // s's endpoint on t
            }
        }
    }
    res
}

pub fn validate_no_crossing(tm: &Emesh<'_>, stage: &str) -> i32 {
    let cs = find_tedge_contacts(tm);
    for c in &cs {
        println!("[crossing] {}: {} at face {}, teid {} x teid {}", stage, if c.cross { "cross" } else { "touch" }, c.fid, c.te_seg, c.te_ndp);
    }
    if !cs.is_empty() { println!("[crossing] {}: {} total", stage, cs.len()); }
    cs.len() as i32
}

/// `repair_crossing_tedges(tm, max_iter = 10)`
pub fn repair_crossing_tedges(tm: &mut Emesh<'_>) { repair_crossing_tedges_n(tm, 10); }

// resolve contacts by re-tracing the intruding tedge (fall back to the other)
pub fn repair_crossing_tedges_n(tm: &mut Emesh<'_>, max_iter: i32) {
    for _ in 0..max_iter {
        let cs = find_tedge_contacts(tm);
        if cs.is_empty() { return; }
        for c in &cs {
            if tm.tedges[c.te_ndp as usize].id != -1 && tm.reroute_tedge(c.te_ndp) { println!("[crossing] rerouted teid {} (contact with teid {})", c.te_ndp, c.te_seg); continue; }
            if tm.tedges[c.te_seg as usize].id != -1 && tm.reroute_tedge(c.te_seg) { println!("[crossing] rerouted teid {} (contact with teid {})", c.te_seg, c.te_ndp); }
        }
    }
}

// live tquads with a side of zero total length: they map onto a segment or a point in tutte, so none may be
// left once the collapse is done
pub fn find_zero_width_tquads(tm: &Emesh<'_>) -> Vec<i32> {
    let mut res = vec![];
    for tq in tm.live_tquads() {
        for s in 0..4 {
            let mut x = 0.;
            for thid in tq.thids(s) { x += tm.thalfs[thid as usize].x; }
            if x == 0. { res.push(tq.id); break; }
        }
    }
    res
}

pub fn validate_collapse_done(tm: &Emesh<'_>) {
    let tqids = find_zero_width_tquads(tm);
    metriko_check!(tqids.is_empty(), "t-mesh collapse left {} zero-width tquads: {:?}", tqids.len(), tqids.iter().take(20).collect::<Vec<_>>());
}
