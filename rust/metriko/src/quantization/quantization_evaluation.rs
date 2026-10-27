//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/quantization/quantization_evaluation.h (see rust/PORTED_FROM).

use faer::sparse::Triplet;
use crate::quantization::quantization_basisloop::*;
use crate::solver::matrix_ops::*;
use crate::tmesh::emesh::*;

pub fn construct_generating_vectors<F: Fn(&Comparator, &Comparator) -> bool + Sync>(
    tm: &Emesh<'_>,
    th2side: &[i32],
    r: &[f64],
    compare: F,
) -> SprsD {
    let canos: Vec<&Ehalf> = tm.thalfs.iter().filter(|th| th.cano).collect();

    // `#pragma omp parallel for`: every loop is independent of the others
    #[cfg(feature = "rayon")]
    let loops: Vec<Vec<i32>> = {
        use rayon::prelude::*;
        canos.par_iter().map(|th| gen_basis_loop(tm, th2side, r, th, &compare)).collect()
    };
    #[cfg(not(feature = "rayon"))]
    let loops: Vec<Vec<i32>> = canos.iter().map(|th| gen_basis_loop(tm, th2side, r, th, &compare)).collect();

    let mut t: Vec<TripD> = vec![];
    for i in 0..loops.len() {
        for &thid in &loops[i] {
            t.push(Triplet::new(i, tm.thalfs[thid as usize].teid as usize, 1.));
        }
    }

    let mut g = from_triplets(tm.tedges.len(), tm.tedges.len(), &t);

    for (rr, c, v) in nonzeros(g.as_ref()) {
        metriko_check!(v <= 2., "generating loop {} passes tedge {} more than twice", rr, c);
    }

    reduce_to_linearly_independent(&mut g);

    transpose(g.as_ref())
}
