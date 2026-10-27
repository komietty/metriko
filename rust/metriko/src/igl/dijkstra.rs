//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/dijkstra.cpp (MPL-2.0).

use std::collections::BTreeSet;
use crate::common::typedef::Row3d;
use crate::common::r#loop::DistId;

/// `igl::dijkstra(source, targets, VV, weights, min_distance, previous)`: Dijkstra from `source` on the vertex graph
/// `VV`, leaving vertex `u` at cost `weights[u]`, until the first vertex in `targets` is popped. Returns
/// `(that vertex or -1, min_distance, previous)`. The queue is a `std::set<std::pair<double, int>>`: it pops the
/// smallest `(distance, vertex)` pair.
pub fn dijkstra_weighted(source: i32, targets: &BTreeSet<i32>, vv: &[Vec<i32>], weights: &[f64]) -> (i32, Vec<f64>, Vec<i32>) {
    let num_v = vv.len();
    let mut min_distance = vec![f64::MAX; num_v];
    min_distance[source as usize] = 0.;
    let mut previous = vec![-1i32; num_v];
    let mut vertex_queue: BTreeSet<DistId> = BTreeSet::new();
    vertex_queue.insert(DistId(min_distance[source as usize], source));

    while let Some(DistId(dist, u)) = vertex_queue.pop_first() {
        if targets.contains(&u) { return (u, min_distance, previous); }

        // Visit each edge exiting u
        for &v in &vv[u as usize] {
            let distance_through_u = dist + weights[u as usize];
            if distance_through_u < min_distance[v as usize] {
                vertex_queue.remove(&DistId(min_distance[v as usize], v));

                min_distance[v as usize] = distance_through_u;
                previous[v as usize] = u;
                vertex_queue.insert(DistId(min_distance[v as usize], v));
            }
        }
    }
    //we should never get here
    (-1, min_distance, previous)
}

/// `igl::dijkstra(source, targets, VV, min_distance, previous)`: every vertex weighs 1.
pub fn dijkstra(source: i32, targets: &BTreeSet<i32>, vv: &[Vec<i32>]) -> (i32, Vec<f64>, Vec<i32>) {
    let weights = vec![1.0; vv.len()];
    dijkstra_weighted(source, targets, vv, &weights)
}

/// `igl::dijkstra(vertex, previous, path)`: the path from `vertex` back to the source along `previous`.
pub fn dijkstra_path(vertex: i32, previous: &[i32]) -> Vec<i32> {
    let mut path = vec![];
    let mut source = vertex;
    while source != -1 {
        path.push(source);
        source = previous[source as usize];
    }
    path
}

/// `igl::dijkstra(V, VV, source, targets, min_distance, previous)`: edges weigh their euclidean length.
pub fn dijkstra_euclidean(v: &[Row3d], vv: &[Vec<i32>], source: i32, targets: &BTreeSet<i32>) -> (i32, Vec<f64>, Vec<i32>) {
    let num_v = vv.len();
    let mut min_distance = vec![f64::INFINITY; num_v];
    min_distance[source as usize] = 0.;
    let mut previous = vec![-1i32; num_v];
    let mut vertex_queue: BTreeSet<DistId> = BTreeSet::new();
    vertex_queue.insert(DistId(min_distance[source as usize], source));

    while let Some(DistId(dist, u)) = vertex_queue.pop_first() {
        if targets.contains(&u) { return (u, min_distance, previous); }

        // Visit each edge exiting u
        for &w in &vv[u as usize] {
            let distance_through_u = dist + (v[u as usize] - v[w as usize]).length();
            if distance_through_u < min_distance[w as usize] {
                vertex_queue.remove(&DistId(min_distance[w as usize], w));

                min_distance[w as usize] = distance_through_u;
                previous[w as usize] = u;
                vertex_queue.insert(DistId(min_distance[w as usize], w));
            }
        }
    }
    (-1, min_distance, previous)
}
