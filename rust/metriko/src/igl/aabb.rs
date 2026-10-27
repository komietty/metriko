//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/AABB.cpp (MPL-2.0): `init(V, Ele)` and `squared_distance(V, Ele, P, sqrD, I, C)`
//--- for triangles in 3D, with the parts of barycenter.cpp, sort.cpp and point_simplex_squared_distance.cpp they call,
//--- and Eigen 5.0.1 AlignedBox (extend, contains, squaredExteriorDistance).

use crate::common::typedef::Row3d;
use crate::common::stdcxx;

/// `Eigen::AlignedBox<double, 3>`
#[derive(Clone, Copy, Debug)]
struct Box3 {
    min: [f64; 3],
    max: [f64; 3],
}

impl Box3 {
    /// `setEmpty()`: min = highest, max = lowest
    fn empty() -> Self { Box3 { min: [f64::MAX; 3], max: [f64::MIN; 3] } }

    /// `extend(p)`: `m_min = m_min.cwiseMin(p)`, `m_max = m_max.cwiseMax(p)` (std::min / std::max)
    fn extend(&mut self, p: Row3d) {
        let p = p.to_array();
        for k in 0..3 {
            if p[k] < self.min[k] { self.min[k] = p[k]; }
            if self.max[k] < p[k] { self.max[k] = p[k]; }
        }
    }

    fn contains(&self, p: Row3d) -> bool {
        let p = p.to_array();
        (0..3).all(|k| self.min[k] <= p[k]) && (0..3).all(|k| p[k] <= self.max[k])
    }

    fn squared_exterior_distance(&self, p: Row3d) -> f64 {
        let p = p.to_array();
        let mut dist2 = 0.;
        for k in 0..3 {
            if self.min[k] > p[k] {
                let aux = self.min[k] - p[k];
                dist2 += aux * aux;
            } else if p[k] > self.max[k] {
                let aux = p[k] - self.max[k];
                dist2 += aux * aux;
            }
        }
        dist2
    }

    /// the index of the first maximum of `diagonal() = max - min` (`maxCoeff(&max_d)`)
    fn longest_dim(&self) -> usize {
        let d = [self.max[0] - self.min[0], self.max[1] - self.min[1], self.max[2] - self.min[2]];
        let mut best = 0;
        for k in 1..3 { if d[k] > d[best] { best = k; } }
        best
    }
}

#[derive(Clone, Debug)]
struct Node {
    bbox: Box3,
    left: Option<usize>,
    right: Option<usize>,
    primitive: i32,
}

/// `igl::AABB<MatrixXd, 3>` over triangles. Nodes live in an arena; node 0 is the root.
#[derive(Clone, Debug)]
pub struct Aabb {
    nodes: Vec<Node>,
}

/// `igl::barycenter(V, F, BC)`: `((0 + v0) + v1) + v2`, divided by 3
fn barycenter(v: &[Row3d], ele: &[[i32; 3]]) -> Vec<[f64; 3]> {
    ele.iter().map(|t| {
        let mut bc = [0.; 3];
        for j in 0..3 {
            let p = v[t[j] as usize].to_array();
            for k in 0..3 { bc[k] += p[k]; }
        }
        for k in 0..3 { bc[k] /= 3.; }
        bc
    }).collect()
}

/// `igl::sort(X, 1, true, Y, IX)` on a (n, 3) matrix: `IX` column by column (`IX[d][i]`)
fn sort_columns(x: &[[f64; 3]]) -> [Vec<usize>; 3] {
    let n = x.len();
    let mut ix: [Vec<usize>; 3] = [vec![], vec![], vec![]];
    match n {
        2 => {
            // igl::sort2
            for d in 0..3 {
                let (mut a, mut b) = (x[0][d], x[1][d]);
                let (mut ai, mut bi) = (0, 1);
                if a > b { std::mem::swap(&mut a, &mut b); std::mem::swap(&mut ai, &mut bi); }
                ix[d] = vec![ai, bi];
            }
        }
        3 => {
            // igl::sort3
            for d in 0..3 {
                let (mut a, mut b, mut c) = (x[0][d], x[1][d], x[2][d]);
                let (mut ai, mut bi, mut ci) = (0, 1, 2);
                if a > b { std::mem::swap(&mut a, &mut b); std::mem::swap(&mut ai, &mut bi); }
                if b > c {
                    std::mem::swap(&mut b, &mut c);
                    std::mem::swap(&mut bi, &mut ci);
                    if a > b { std::mem::swap(&mut a, &mut b); std::mem::swap(&mut ai, &mut bi); }
                }
                let _ = (a, b, c);
                ix[d] = vec![ai, bi, ci];
            }
        }
        _ => {
            for d in 0..3 {
                // igl::sort(std::vector) : std::sort of the index map with IndexLessThan
                let data: Vec<f64> = x.iter().map(|r| r[d]).collect();
                let mut index_map: Vec<usize> = (0..n).collect();
                stdcxx::sort_by(&mut index_map, |&a, &b| data[a] < data[b]);
                ix[d] = index_map;
            }
        }
    }
    ix
}

/// `igl::point_simplex_squared_distance<3>(p, V, Ele, primitive, sqr_d, c)` for a triangle: `(sqr_d, c)`
pub fn point_triangle_squared_distance(p: Row3d, v: &[Row3d], ele: &[[i32; 3]], primitive: i32) -> (f64, Row3d) {
    let t = ele[primitive as usize];
    let a = v[t[0] as usize];
    let b = v[t[1] as usize];
    let c = v[t[2] as usize];
    let closest = (|| {
        let ab = b - a;
        let ac = c - a;
        let ap = p - a;
        let d1 = ab.dot(ap);
        let d2 = ac.dot(ap);
        if d1 <= 0.0 && d2 <= 0.0 { return a; }
        let bp = p - b;
        let d3 = ab.dot(bp);
        let d4 = ac.dot(bp);
        if d3 >= 0.0 && d4 <= d3 { return b; }
        let vc = d1 * d4 - d3 * d2;
        if a != b && vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            return a + ab * v;
        }
        let cp = p - c;
        let d5 = ab.dot(cp);
        let d6 = ac.dot(cp);
        if d6 >= 0.0 && d5 <= d6 { return c; }
        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            return a + ac * w;
        }
        let va = d3 * d6 - d5 * d4;
        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + (c - b) * w;
        }
        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        a + ab * v + ac * w // = u*a + v*b + w*c, u = va * denom = 1.0-v-w
    })();
    ((p - closest).length_squared(), closest)
}

impl Aabb {
    /// `tree.init(V, Ele)`
    pub fn init(v: &[Row3d], ele: &[[i32; 3]]) -> Aabb {
        let mut tree = Aabb { nodes: vec![] };
        let root = Node { bbox: Box3::empty(), left: None, right: None, primitive: -1 };
        tree.nodes.push(root);
        if v.is_empty() || ele.is_empty() { return tree; }
        let all_i: Vec<i32> = (0..ele.len() as i32).collect();
        let bc = barycenter(v, ele);
        // SI(i, d): which place i would be sorted into along d
        let is = sort_columns(&bc);
        let mut si = vec![[0i32; 3]; bc.len()];
        for d in 0..3 {
            for i in 0..is[d].len() { si[is[d][i]][d] = i as i32; }
        }
        tree.init_rec(0, v, ele, &si, &all_i);
        tree
    }

    /// `init(V, Ele, SI, I)`
    fn init_rec(&mut self, node: usize, v: &[Row3d], ele: &[[i32; 3]], si: &[[i32; 3]], ii: &[i32]) {
        if v.is_empty() || ele.is_empty() || ii.is_empty() { return; }
        let mut bbox = Box3::empty();
        // Compute bounding box
        for &i in ii {
            for c in 0..3 {
                bbox.extend(v[ele[i as usize][c] as usize]);
                bbox.extend(v[ele[i as usize][c] as usize]);
            }
        }
        self.nodes[node].bbox = bbox;
        if ii.len() == 1 {
            self.nodes[node].primitive = ii[0];
            return;
        }
        // Compute longest direction
        let max_d = bbox.longest_dim();
        // Can't use median on BC directly because many may have same value,
        // but can use median on sorted BC indices
        let sidi: Vec<i32> = ii.iter().map(|&i| si[i as usize][max_d]).collect();
        let med = {
            // std::nth_element on a copy; the ranks are distinct, so the n-th smallest value is unique
            let mut a = sidi.clone();
            let n = (a.len() - 1) / 2;
            *a.select_nth_unstable(n).1
        };
        let mut li = Vec::with_capacity(ii.len().div_ceil(2));
        let mut ri = Vec::with_capacity(ii.len() / 2);
        // Distribute left and right
        for k in 0..ii.len() {
            if sidi[k] <= med { li.push(ii[k]); } else { ri.push(ii[k]); }
        }
        if !li.is_empty() {
            let l = self.nodes.len();
            self.nodes.push(Node { bbox: Box3::empty(), left: None, right: None, primitive: -1 });
            self.nodes[node].left = Some(l);
            self.init_rec(l, v, ele, si, &li);
        }
        if !ri.is_empty() {
            let r = self.nodes.len();
            self.nodes.push(Node { bbox: Box3::empty(), left: None, right: None, primitive: -1 });
            self.nodes[node].right = Some(r);
            self.init_rec(r, v, ele, si, &ri);
        }
    }

    fn is_leaf(&self, node: usize) -> bool { self.nodes[node].left.is_none() && self.nodes[node].right.is_none() }

    /// `set_min`
    fn set_min(sqr_d_candidate: f64, i_candidate: i32, c_candidate: Row3d, sqr_d: &mut f64, i: &mut i32, c: &mut Row3d) {
        if sqr_d_candidate < *sqr_d {
            *i = i_candidate;
            *c = c_candidate;
            *sqr_d = sqr_d_candidate;
        }
    }

    /// `leaf_squared_distance(V, Ele, p, low_sqr_d, sqr_d, i, c)`
    fn leaf_squared_distance(&self, node: usize, v: &[Row3d], ele: &[[i32; 3]], p: Row3d, low_sqr_d: f64, sqr_d: &mut f64, i: &mut i32, c: &mut Row3d) {
        if low_sqr_d > *sqr_d {
            *sqr_d = low_sqr_d;
            return;
        }
        let prim = self.nodes[node].primitive;
        let (sqr_d_candidate, c_candidate) = point_triangle_squared_distance(p, v, ele, prim);
        Self::set_min(sqr_d_candidate, prim, c_candidate, sqr_d, i, c);
    }

    /// `squared_distance(V, Ele, p, low_sqr_d, up_sqr_d, i, c)`
    fn squared_distance_rec(&self, node: usize, v: &[Row3d], ele: &[[i32; 3]], p: Row3d, low_sqr_d: f64, up_sqr_d: f64, i: &mut i32, c: &mut Row3d) -> f64 {
        if low_sqr_d > up_sqr_d { return low_sqr_d; }
        let mut sqr_d = up_sqr_d;
        if self.is_leaf(node) {
            self.leaf_squared_distance(node, v, ele, p, low_sqr_d, &mut sqr_d, i, c);
        } else {
            let left = self.nodes[node].left.unwrap();
            let right = self.nodes[node].right.unwrap();
            let mut looked_left = false;
            let mut looked_right = false;
            let look_left = |sqr_d: &mut f64, i: &mut i32, c: &mut Row3d, looked_left: &mut bool| {
                let mut i_left = -1;
                let mut c_left = *c;
                let sqr_d_left = self.squared_distance_rec(left, v, ele, p, low_sqr_d, *sqr_d, &mut i_left, &mut c_left);
                Self::set_min(sqr_d_left, i_left, c_left, sqr_d, i, c);
                *looked_left = true;
            };
            let look_right = |sqr_d: &mut f64, i: &mut i32, c: &mut Row3d, looked_right: &mut bool| {
                let mut i_right = -1;
                let mut c_right = *c;
                let sqr_d_right = self.squared_distance_rec(right, v, ele, p, low_sqr_d, *sqr_d, &mut i_right, &mut c_right);
                Self::set_min(sqr_d_right, i_right, c_right, sqr_d, i, c);
                *looked_right = true;
            };

            // must look left or right if in box
            if self.nodes[left].bbox.contains(p) { look_left(&mut sqr_d, i, c, &mut looked_left); }
            if self.nodes[right].bbox.contains(p) { look_right(&mut sqr_d, i, c, &mut looked_right); }
            // if haven't looked left and could be less than current min, then look
            let left_up_sqr_d = self.nodes[left].bbox.squared_exterior_distance(p);
            let right_up_sqr_d = self.nodes[right].bbox.squared_exterior_distance(p);
            if left_up_sqr_d < right_up_sqr_d {
                if !looked_left && left_up_sqr_d < sqr_d { look_left(&mut sqr_d, i, c, &mut looked_left); }
                if !looked_right && right_up_sqr_d < sqr_d { look_right(&mut sqr_d, i, c, &mut looked_right); }
            } else {
                if !looked_right && right_up_sqr_d < sqr_d { look_right(&mut sqr_d, i, c, &mut looked_right); }
                if !looked_left && left_up_sqr_d < sqr_d { look_left(&mut sqr_d, i, c, &mut looked_left); }
            }
        }
        sqr_d
    }

    /// `squared_distance(V, Ele, p, i, c)` (up = infinity)
    pub fn squared_distance_point(&self, v: &[Row3d], ele: &[[i32; 3]], p: Row3d) -> (f64, i32, Row3d) {
        let mut i = -1;
        let mut c = Row3d::ZERO;
        let d = self.squared_distance_rec(0, v, ele, p, 0.0, f64::INFINITY, &mut i, &mut c);
        (d, i, c)
    }

    /// `squared_distance(V, Ele, P, sqrD, I, C)`: returns `(sqrD, I, C)`
    pub fn squared_distance(&self, v: &[Row3d], ele: &[[i32; 3]], p: &[Row3d]) -> (Vec<f64>, Vec<i32>, Vec<Row3d>) {
        let mut sqr_d = Vec::with_capacity(p.len());
        let mut ii = Vec::with_capacity(p.len());
        let mut cc = Vec::with_capacity(p.len());
        for &pp in p {
            let (d, i, c) = self.squared_distance_point(v, ele, pp);
            sqr_d.push(d);
            ii.push(i);
            cc.push(c);
        }
        (sqr_d, ii, cc)
    }
}
