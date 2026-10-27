//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/flip_avoiding_line_search.cpp and include/igl/line_search.cpp (MPL-2.0), for
//--- the 2D case metriko uses (uv with two columns; the volumetric branch, SolveP3 and get_min_pos_root_3D are not
//--- reachable with a 2 column uv and are not ported).

/// `igl::flip_avoiding::get_smallest_pos_quad_zero(a, b, c)`
fn get_smallest_pos_quad_zero(a: f64, b: f64, c: f64) -> f64 {
    let mut t1;
    let mut t2;
    if a.abs() > 1.0e-10 {
        let delta_in = b.powf(2.) - 4. * a * c;
        if delta_in <= 0. {
            return f64::INFINITY;
        }

        let delta = delta_in.sqrt(); // delta >= 0
        if b >= 0. {
            // avoid subtracting two similar numbers
            let bd = -b - delta;
            t1 = 2. * c / bd;
            t2 = bd / (2. * a);
        } else {
            let bd = -b + delta;
            t1 = bd / (2. * a);
            t2 = (2. * c) / bd;
        }

        if a < 0. { std::mem::swap(&mut t1, &mut t2); } // make t1 > t2
        // return the smaller positive root if it exists, otherwise return infinity
        if t1 > 0. {
            if t2 > 0. { t2 } else { t1 }
        } else {
            f64::INFINITY
        }
    } else {
        if b == 0. { return f64::INFINITY; } // just to avoid divide-by-zero
        t1 = -c / b;
        if t1 > 0. { t1 } else { f64::INFINITY }
    }
}

/// `igl::flip_avoiding::get_min_pos_root_2D(uv, F, d, f)`: the smallest positive step t at which face f degenerates
/// when its corners move from uv along d
fn get_min_pos_root_2d(uv: &[[f64; 2]], f: &[[i32; 3]], d: &[[f64; 2]], fi: usize) -> f64 {
    let v1 = f[fi][0] as usize;
    let v2 = f[fi][1] as usize;
    let v3 = f[fi][2] as usize;
    // get quadratic coefficients (ax^2 + b^x + c)
    let u11 = uv[v1][0];
    let u12 = uv[v1][1];
    let u21 = uv[v2][0];
    let u22 = uv[v2][1];
    let u31 = uv[v3][0];
    let u32 = uv[v3][1];

    let v11 = d[v1][0];
    let v12 = d[v1][1];
    let v21 = d[v2][0];
    let v22 = d[v2][1];
    let v31 = d[v3][0];
    let v32 = d[v3][1];

    let a = v11 * v22 - v12 * v21 - v11 * v32 + v12 * v31 + v21 * v32 - v22 * v31;
    let b = u11 * v22 - u12 * v21 - u21 * v12 + u22 * v11 - u11 * v32 + u12 * v31 + u31 * v12 - u32 * v11 + u21 * v32 - u22 * v31 - u31 * v22 + u32 * v21;
    let c = u11 * u22 - u12 * u21 - u11 * u32 + u12 * u31 + u21 * u32 - u22 * u31;

    get_smallest_pos_quad_zero(a, b, c)
}

/// `igl::flip_avoiding::compute_max_step_from_singularities(uv, F, d)` (2D branch)
fn compute_max_step_from_singularities(uv: &[[f64; 2]], f: &[[i32; 3]], d: &[[f64; 2]]) -> f64 {
    let mut max_step = f64::INFINITY;
    for fi in 0..f.len() {
        let min_positive_root = get_min_pos_root_2d(uv, f, d, fi);
        // std::min(max_step, min_positive_root)
        max_step = if min_positive_root < max_step { min_positive_root } else { max_step };
    }
    max_step
}

/// `igl::line_search(x, d, step_size, energy, cur_energy)`: halves the step until the energy decreases, at most 12
/// times. Updates x and returns the new energy.
pub fn line_search(
    x: &mut Vec<[f64; 2]>,
    d: &[[f64; 2]],
    mut step_size: f64,
    energy: &dyn Fn(&[[f64; 2]]) -> f64,
    cur_energy: f64,
) -> f64 {
    let old_energy = if cur_energy > 0. {
        cur_energy
    } else {
        energy(x) // no energy was given -> need to compute the current energy
    };
    let mut new_energy = old_energy;
    let mut cur_iter = 0;
    let max_step_size_iter = 12;

    while new_energy >= old_energy && cur_iter < max_step_size_iter {
        let new_x: Vec<[f64; 2]> = x.iter().zip(d.iter()).map(|(p, q)| [p[0] + step_size * q[0], p[1] + step_size * q[1]]).collect();

        let cur_e = energy(&new_x);
        if cur_e >= old_energy {
            step_size /= 2.;
        } else {
            *x = new_x;
            new_energy = cur_e;
        }
        cur_iter += 1;
    }
    new_energy
}

/// `igl::flip_avoiding_line_search(F, cur_v, dst_v, energy, cur_energy)`: a line search from cur_v toward dst_v
/// that stops short (0.8) of the first step at which a face would flip. Updates cur_v, returns the new energy.
pub fn flip_avoiding_line_search(
    f: &[[i32; 3]],
    cur_v: &mut Vec<[f64; 2]>,
    dst_v: &[[f64; 2]],
    energy: &dyn Fn(&[[f64; 2]]) -> f64,
    cur_energy: f64,
) -> f64 {
    let d: Vec<[f64; 2]> = dst_v.iter().zip(cur_v.iter()).map(|(p, q)| [p[0] - q[0], p[1] - q[1]]).collect();

    let min_step_to_singularity = compute_max_step_from_singularities(cur_v, f, &d);
    // std::min(1., min_step_to_singularity * 0.8)
    let s = min_step_to_singularity * 0.8;
    let max_step_size = if s < 1. { s } else { 1. };

    line_search(cur_v, &d, max_step_size, energy, cur_energy)
}
