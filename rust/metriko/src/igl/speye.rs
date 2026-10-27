//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libigl v2.6.4 include/igl/speye.cpp (MPL-2.0).

use faer::sparse::Triplet;
use crate::solver::matrix_ops::{SprsD, TripD, from_triplets};

/// `igl::speye(m, n, I)`: the m by n sparse identity (ones on the first min(m, n) diagonal entries)
pub fn speye_mn(m: usize, n: usize) -> SprsD {
    // size of diagonal
    let d = if m < n { m } else { n };
    let t: Vec<TripD> = (0..d).map(|i| Triplet::new(i, i, 1.0)).collect();
    from_triplets(m, n, &t)
}

/// `igl::speye(n, I)`
pub fn speye(n: usize) -> SprsD { speye_mn(n, n) }
