//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Ports of the Eigen 5.0.1 decompositions metriko relies on where faer has no equivalent with the same
//--- discrete behavior (rank revealing, pivot order, fill-reducing column order).

pub mod full_piv_lu;
pub mod colamd;
pub mod coletree;
pub mod sparse_qr;
