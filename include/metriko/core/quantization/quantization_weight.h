//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_WEIGHT_H
#define METRIKO_QUANTIZATION_WEIGHT_H

namespace metriko {
    inline double compute_weight(
        const double ideal_value,
        const double curr_value,
        const int num_tedges_in_tmesh,
        const bool minus = false
    ) {
        if (minus && curr_value == 0) return 1e9;
        double d = (ideal_value - curr_value) * (minus ? -1. : 1.);
        double n = num_tedges_in_tmesh;
        if (d >= 1.) return 1. / (d + 1.);
        if (d < 0.)  return n * n / (1. - d);
        return n / (d + 1.);
    }
}

#endif
