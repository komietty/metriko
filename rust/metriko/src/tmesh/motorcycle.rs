//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/tmesh/motorcycle.h (see rust/PORTED_FROM).

use std::collections::HashMap;
use crate::common::typedef::*;
use crate::common::utilities::*;
use crate::common::predicates::*;
use crate::common::stdcxx;
use crate::hmesh::hmesh::*;
use crate::hmesh::hmloc::*;
use crate::hmesh::hmray::*;
use crate::hmesh::utilities::get_face_uv;

pub const TOLERANCE_EDGE_AB: f64 = 1e-6; // tolerance on two curvs crash close to an edge
pub const TOLERANCE_EDGE_CD: f64 = 1e-6; // tolerance on two curvs crash close to an edge

#[derive(Clone, Copy, Debug, Default)]
pub struct Mport {
    pub dir: Complex,
    pub crnr_id: i32,
    pub this_id: i32,
    pub next_id: i32,
    pub prev_id: i32,
}

impl Mport {
    fn new(dir: Complex, crnr_id: i32) -> Self { Mport { dir, crnr_id, this_id: -1, next_id: -1, prev_id: -1 } }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum JunctionType { #[default] None, F, T, C }

#[derive(Clone, Debug, Default)]
pub struct Mnode {
    pub jt: JunctionType,
    pub loc: HmLoc,
    pub adj: Vec<[i32; 2]>, // [0]: curv_id, [1]: sgmt_id
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Msgmt {
    pub fr_nid: i32,
    pub to_nid: i32,
    pub this_id: i32,
    pub curv_id: i32,
    pub face_id: i32,
}

impl Default for Msgmt {
    fn default() -> Self { Msgmt { fr_nid: -1, to_nid: -1, this_id: -1, curv_id: -1, face_id: -1 } }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Mbuff {
    pub loc: HmLoc,   // curve front, chart-aware (OnC or OnH)
    pub dir: Complex, // direction in the chart of get_chart_face(loc)
    pub end: bool,
}

#[derive(Clone, Debug)]
pub struct Mcurv {
    pub id: i32,
    pub sgmts: Vec<Msgmt>,
    pub buff: Mbuff,
}

impl PartialEq for Mcurv {
    fn eq(&self, c: &Self) -> bool { self.id == c.id }
}

pub struct Mgrph<'a> {
    pub hm: &'a Hmesh,
    pub mports: Vec<Mport>,
    pub mnodes: Vec<Mnode>,
    pub mcurvs: Vec<Mcurv>,
}

impl<'a> Mgrph<'a> {
    pub fn new(hm: &'a Hmesh, matching: &[i32], singular: &[i32]) -> Mgrph<'a> {
        let mut mg = Mgrph { hm, mports: vec![], mnodes: vec![], mcurvs: vec![] };
        mg.gen_ports(singular);

        for v in hm.verts().filter(|v| singular[v.id as usize] != 0) {
            mg.mnodes.push(Mnode { jt: JunctionType::F, loc: HmLoc::OnV { id: v.id }, adj: vec![] });
        }

        // 1: Add the first segment for each curve
        mg.mcurvs.reserve(mg.mports.len());
        for i in 0..mg.mports.len() {
            let p = mg.mports[i];
            mg.mcurvs.push(Mcurv { id: p.this_id, sgmts: vec![], buff: Mbuff { loc: HmLoc::OnC { id: p.crnr_id }, dir: p.dir, end: false } });
            let k = mg.mcurvs.len() - 1;
            mg.add_segment(k);
        }

        // 2: Add further segments until every curve crash to another curve
        while mg.mcurvs.iter().any(|c| !c.buff.end) {
            for k in 0..mg.mcurvs.len() {
                if mg.mcurvs[k].buff.end { continue; }
                let (loc, dir) = cross_to_twin(hm, &mg.mcurvs[k].buff.loc, matching, mg.mcurvs[k].buff.dir);
                mg.mcurvs[k].buff.loc = loc;
                mg.mcurvs[k].buff.dir = dir;
                mg.add_segment(k);
            }
        }

        // 3: Collect node adjacency
        for mn in mg.mnodes.iter_mut() { mn.adj.clear(); }
        for mc in &mg.mcurvs {
            for i in 0..mc.sgmts.len() {
                let sg = &mc.sgmts[i];
                let ad = [mc.id, i as i32];
                if sg.fr_nid != -1 { mg.mnodes[sg.fr_nid as usize].adj.push(ad); }
                if sg.to_nid != -1 { mg.mnodes[sg.to_nid as usize].adj.push(ad); }
            }
        }

        mg.sort_node_adjacency();

        for c in mg.mcurvs.iter_mut() {
            for i in 0..c.sgmts.len() {
                c.sgmts[i].curv_id = c.id;
                c.sgmts[i].this_id = i as i32;
            }
        }
        mg
    }

    pub fn gen_ports(&mut self, singular: &[i32]) {
        let hm = self.hm;
        for v in hm.verts() {
            if singular[v.id as usize] == 0 { continue; }
            let mut buff0: Vec<Mport> = vec![]; // the outer scope buffer to assign next/prev
            let mut buff1: Vec<Mport> = vec![]; // the inner scope buffer

            for cr in v.adj_crnrs() {
                buff1.clear();
                let a = cr.uv();
                let b = cr.next().uv();
                let c = cr.prev().uv();
                let o = orientation(a, b, c);
                let ab = b - a;
                let ac = c - a;
                if o < 0. { panic!("Orientation should be ccw order"); }

                let mut r = 0;
                while r < 4 {
                    let d = get_quater_rot(r);
                    if !is_points_into(a, b, c, a + d, EPS) && r > 0 { break; }
                    r += 1;
                }

                for i in 0..4 {
                    let d = get_quater_rot(r - i);
                    if is_points_into(a, b, c, a + d, EPS) ||
                        d.arg() == ab.arg() ||
                        d.arg() == ac.arg()
                    { buff1.push(Mport::new(d, cr.id)); }
                }

                stdcxx::sort_by(&mut buff1, |p0, p1| dot(p0.dir, ab) > dot(p1.dir, ab));
                buff0.extend_from_slice(&buff1);
            }

            let base = self.mports.len() as i32;
            for i in 0..buff0.len() { buff0[i].this_id = base + i as i32; }
            for i in 0..buff0.len() {
                let s = buff0.len() as i32;
                let ii = i as i32;
                buff0[i].prev_id = buff0[((ii - 1 + s) % s) as usize].this_id;
                buff0[i].next_id = buff0[((ii + 1 + s) % s) as usize].this_id;
            }
            self.mports.extend_from_slice(&buff0);
        }
    }

    pub fn sort_node_adjacency(&mut self) {
        let hm = self.hm;
        for nid in 0..self.mnodes.len() {
            if self.mnodes[nid].adj.len() < 3 { continue; }
            let mut adj = std::mem::take(&mut self.mnodes[nid].adj);
            {
                let mnodes = &self.mnodes;
                let mcurvs = &self.mcurvs;
                let get_fid = |ad: &[i32; 2]| mcurvs[ad[0] as usize].sgmts[ad[1] as usize].face_id;
                let get_dir = |ad: &[i32; 2]| -> Complex {
                    let s = &mcurvs[ad[0] as usize].sgmts[ad[1] as usize];
                    let uv_a = get_face_uv(&mnodes[s.fr_nid as usize].loc, s.face_id, hm);
                    let uv_b = get_face_uv(&mnodes[s.to_nid as usize].loc, s.face_id, hm);
                    metriko_check!((uv_a - uv_b).norm() > 1e-8, "zero-length segment at node {}", nid);
                    if s.fr_nid == nid as i32 { uv_b - uv_a } else { uv_a - uv_b }
                };

                let sort_by_rank = |adj: &mut Vec<[i32; 2]>, rank: &dyn Fn(i32) -> i32| {
                    stdcxx::sort_by(adj, |a, b| {
                        let r_a = rank(get_fid(a));
                        let r_b = rank(get_fid(b));
                        if r_a != r_b { return r_a < r_b; }
                        cross(get_dir(a), get_dir(b)) > 0.
                    });
                };

                match mnodes[nid].loc {
                    HmLoc::OnV { id } => {
                        let mut r = 0;
                        let mut ccw_rank: HashMap<i32, i32> = HashMap::new();
                        for h in hm.vert(id).adj_halfs() { ccw_rank.insert(h.face().id, r); r += 1; }
                        sort_by_rank(&mut adj, &|fid| *ccw_rank.get(&fid).expect("unordered_map::at"));
                    }
                    HmLoc::OnE { id, .. } => {
                        let e = hm.edge(id);
                        sort_by_rank(&mut adj, &|fid| {
                            if fid == e.face0().id { return 0; }
                            if fid == e.face1().id { return 1; }
                            2
                        });
                    }
                    HmLoc::OnP { .. } => {
                        stdcxx::sort_by(&mut adj, |a, b| get_dir(a).arg() < get_dir(b).arg());
                    }
                    _ => metriko_fail!("unexpected location type"),
                }
            }
            self.mnodes[nid].adj = adj;
        }
    }

    /// `Mcurv::add_segment(hm)` of curve `k` (the C++ reaches the graph through `mg`)
    pub fn add_segment(&mut self, k: usize) {
        let hm = self.hm;
        let id = self.mcurvs[k].id;
        let fr = self.mcurvs[k].buff.loc;
        let bgn = self.mcurvs[k].sgmts.is_empty(); // add_segment always appends one, so only the first call sees it empty
        let uv0 = get_chart_uv(hm, &fr);
        let fid = get_chart_face(hm, &fr).id;
        self.mcurvs[k].buff.loc = find_ray_intersection(hm, &fr, self.mcurvs[k].buff.dir);
        let uv3 = get_chart_uv(hm, &self.mcurvs[k].buff.loc);

        let mut candidates: Vec<(f64, f64, Msgmt)> = vec![];

        {
            let my_vert = hm.crnr(self.mports[id as usize].crnr_id).vert();
            for c in &self.mcurvs {
                if c.id == id { continue; }
                if bgn && hm.crnr(self.mports[c.id as usize].crnr_id).vert() == my_vert { continue; }
                for s in &c.sgmts {
                    if s.face_id != fid { continue; }
                    let uv_a = get_face_uv(&self.mnodes[s.fr_nid as usize].loc, fid, hm);
                    let uv_b = get_face_uv(&self.mnodes[s.to_nid as usize].loc, fid, hm);
                    if let Some((ab, cd)) = find_strict_intersection_ratio(uv0, uv3, uv_a, uv_b, 0.) { candidates.push((ab, cd, *s)); } // closed segments
                }
            }
        }

        let fr_nid = if !self.mcurvs[k].sgmts.is_empty() {
            self.mcurvs[k].sgmts.last().unwrap().to_nid
        } else {
            let l = to_chart_free(hm, &fr);
            let i = self.mnodes.iter().position(|n| n.loc == l);
            metriko_check!(i.is_some(), "start node of the curve must exist");
            i.unwrap() as i32
        };

        let mut sg = Msgmt { fr_nid, curv_id: id, face_id: fid, ..Default::default() };

        // determine to_nid
        if candidates.is_empty() {
            let ml = to_chart_free(hm, &self.mcurvs[k].buff.loc); // OnC -> OnV, OnH -> OnE

            // if hit to other Mnode on a vert, return
            if matches!(ml, HmLoc::OnV { .. }) {
                for i in 0..self.mnodes.len() {
                    if self.mnodes[i].loc == ml {
                        self.mnodes[i].jt = JunctionType::T;
                        sg.to_nid = i as i32;
                        self.mcurvs[k].sgmts.push(sg);
                        self.mcurvs[k].buff.end = true;
                        return;
                    }
                }
            }

            // otherwise
            self.mnodes.push(Mnode { loc: ml, ..Default::default() });
            sg.to_nid = self.mnodes.len() as i32 - 1;
            self.mcurvs[k].sgmts.push(sg);
        }
        // intersection happens
        else {
            // rg::min by get<0>: the first minimum
            let mut best = 0;
            for i in 1..candidates.len() { if candidates[i].0 < candidates[best].0 { best = i; } }
            let (ab, cd, sg_cd) = candidates[best];
            metriko_check!(ab >= TOLERANCE_EDGE_AB, "a curve passes too close to a vertex (not handled yet)");
            let cd_bgn_snappable = cd < TOLERANCE_EDGE_CD;
            let cd_end_snappable = cd > 1. - TOLERANCE_EDGE_CD;

            let cv = sg_cd.curv_id as usize;
            let it = self.mcurvs[cv].sgmts.iter().position(|s| s.fr_nid == sg_cd.fr_nid && s.to_nid == sg_cd.to_nid).expect("find_if");
            self.mcurvs[k].buff.end = true;

            let snap = |mg: &mut Mgrph<'a>, sg: &mut Msgmt, nid: i32| -> bool {
                let jt = mg.mnodes[nid as usize].jt;
                if jt == JunctionType::F || jt == JunctionType::C { return false; }
                // if (jt == JunctionType::None) return false;

                for c in &mg.mcurvs {
                    for s in &c.sgmts {
                        if s.to_nid == nid && s.face_id == fid {
                            let o = get_face_uv(&mg.mnodes[nid as usize].loc, fid, hm);
                            let d0 = get_face_uv(&mg.mnodes[s.fr_nid as usize].loc, fid, hm) - o;
                            let d1 = get_face_uv(&mg.mnodes[sg.fr_nid as usize].loc, fid, hm) - o;
                            let l0 = d0.norm();
                            let l1 = d1.norm();
                            if l0 < EPS || l1 < EPS || dot(d0 / l0, d1 / l1) > EPS { return false; }
                        }
                    }
                }

                mg.mnodes[nid as usize].jt = JunctionType::T;
                sg.to_nid = nid;
                mg.mcurvs[k].sgmts.push(*sg);
                true
            };

            let it_fr = self.mcurvs[cv].sgmts[it].fr_nid;
            if cd_bgn_snappable && snap(self, &mut sg, it_fr) { println!("cd bgn snappable, fid: {}", fid); return; }
            let it_to = self.mcurvs[cv].sgmts[it].to_nid;
            if cd_end_snappable && snap(self, &mut sg, it_to) { println!("cd end snappable, fid: {}", fid); return; }

            { // hit in middle case, or close to singular node
                let mn = Mnode { jt: JunctionType::T, loc: HmLoc::OnP { id: fid, uv: lerp(uv0, uv3, ab) }, adj: vec![] };
                self.mnodes.push(mn);
                let nid = self.mnodes.len() as i32 - 1;

                sg.to_nid = nid;
                self.mcurvs[k].sgmts.push(sg);

                let mut s1 = self.mcurvs[cv].sgmts[it];
                s1.fr_nid = nid;
                self.mcurvs[cv].sgmts[it].to_nid = nid;

                self.mcurvs[cv].sgmts.insert(it + 1, s1);
            }
        }
    }
}
