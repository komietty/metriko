//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/cut_mesh_from_singularities.cpp (MPL-2.0).

use std::collections::VecDeque;
use crate::igl::edge_topology::edge_topology;
use crate::igl::triangle_triangle_adjacency::triangle_triangle_adjacency;

/// `igl::MeshCutter`
struct MeshCutter<'a> {
    n_v: usize,
    f: &'a [[i32; 3]],
    handle_mmatch: &'a [[f64; 3]],
    f_visited: Vec<i32>,
    tt: Vec<[i32; 3]>,
    tti: Vec<[i32; 3]>,
    e: Vec<[i32; 2]>,
    f2e: Vec<[i32; 3]>,
    e2f: Vec<[i32; 2]>,
}

impl<'a> MeshCutter<'a> {
    fn new(n_v: usize, f: &'a [[i32; 3]], handle_mmatch: &'a [[f64; 3]]) -> Self {
        let (tt, tti) = triangle_triangle_adjacency(f);
        let (e, f2e, e2f) = edge_topology(n_v, f);
        MeshCutter { n_v, f, handle_mmatch, f_visited: vec![], tt, tti, e, f2e, e2f }
    }

    fn is_rot_seam(&self, f0: i32, edge: i32) -> bool {
        // unsigned char MM = Handle_MMatch(f0,edge)
        let mm = self.handle_mmatch[f0 as usize][edge as usize] as u8;
        mm != 0
    }

    fn flood_fill(&mut self, start: i32, handle_seams: &mut [[i32; 3]]) {
        let mut d: VecDeque<i32> = VecDeque::new();
        // clean the visited flag
        self.f_visited[start as usize] = 1;
        d.push_back(start);

        while let Some(f) = d.pop_front() {
            for s in 0..3 {
                let mut g = self.tt[f as usize][s]; // f->FFp(s);
                let mut j = self.tti[f as usize][s]; // f->FFi(s);

                if j == -1 {
                    g = f;
                    j = s as i32;
                }

                if !self.is_rot_seam(f, s as i32) && !self.is_rot_seam(g, j) && self.f_visited[g as usize] == 0 {
                    handle_seams[f as usize][s] = 0;
                    handle_seams[g as usize][j as usize] = 0;
                    self.f_visited[g as usize] = 1;
                    d.push_back(g);
                }
            }
        }
    }

    fn retract(&self, handle_seams: &mut [[i32; 3]]) {
        let f = self.f;
        let mut e = vec![0i32; self.n_v]; // number of edges per vert
        for ei in 0..self.e.len() {
            //only need one face
            let mut f0 = self.e2f[ei][0];
            if f0 == -1 { f0 = self.e2f[ei][1]; }
            let mut k = 0;
            while k < 3 {
                if self.f2e[f0 as usize][k] == ei as i32 { break; }
                k += 1;
            }
            if handle_seams[f0 as usize][k] != 0 {
                e[f[f0 as usize][k] as usize] += 1;
                e[f[f0 as usize][(k + 1) % 3] as usize] += 1;
            }
        }

        let mut guard = 0;
        loop {
            let mut over = true;
            for fi in 0..f.len() {
                for s in 0..3 {
                    if handle_seams[fi][s] != 0 && !self.is_rot_seam(fi as i32, s as i32) { // never retract rot seams
                        if e[f[fi][s] as usize] == 1 {
                            // dissolve seam
                            handle_seams[fi][s] = 0;
                            if self.tt[fi][s] != -1 {
                                handle_seams[self.tt[fi][s] as usize][self.tti[fi][s] as usize] = 0;
                            }

                            e[f[fi][s] as usize] -= 1;
                            e[f[fi][(s + 1) % 3] as usize] -= 1;
                            over = false;
                        }
                    }
                }
            }

            let g = guard;
            guard += 1;
            if g > 10000 { over = true; }
            if over { break; }
        }
    }

    fn cut(&mut self) -> Vec<[i32; 3]> {
        let nf = self.f.len();
        self.f_visited = vec![0; nf];
        let mut handle_seams = vec![[1i32; 3]; nf];

        for fi in 0..nf {
            if self.f_visited[fi] == 0 {
                self.flood_fill(fi as i32, &mut handle_seams);
            }
        }

        self.retract(&mut handle_seams);

        for fi in 0..nf {
            for j in 0..3 {
                if self.is_rot_seam(fi as i32, j as i32) { handle_seams[fi][j] = 1; }
            }
        }
        handle_seams
    }
}

/// `igl::cut_mesh_from_singularities(V, F, Handle_MMatch, Handle_Seams)`: returns `Handle_Seams`, 1 on the face
/// edges of the cut graph (a spanning tree of the dual graph is left uncut, dangling seams are retracted, rotational
/// seams `Handle_MMatch(f, k) != 0` are always cut).
pub fn cut_mesh_from_singularities(n_v: usize, f: &[[i32; 3]], handle_mmatch: &[[f64; 3]]) -> Vec<[i32; 3]> {
    let mut mc = MeshCutter::new(n_v, f, handle_mmatch);
    mc.cut()
}
