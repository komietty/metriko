//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of include/metriko/quantization/quantization_weight.h (see rust/PORTED_FROM).

/// `compute_weight(ideal_value, curr_value, num_tedges_in_tmesh, minus = false)`
pub fn compute_weight(ideal_value: f64, curr_value: f64, num_tedges_in_tmesh: i32, minus: bool) -> f64 {
    if minus && curr_value == 0. { return 1e9; }
    let d = (ideal_value - curr_value) * (if minus { -1. } else { 1. });
    let n = num_tedges_in_tmesh as f64;
    if d >= 1. { return 1. / (d + 1.); }
    if d < 0. { return n * n / (1. - d); }
    n / (d + 1.)
}
