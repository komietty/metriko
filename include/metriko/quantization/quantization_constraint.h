//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QUANTIZATION_CONSTRAINT_H
#define METRIKO_QUANTIZATION_CONSTRAINT_H
#include "metriko/tmesh/emesh.h"

namespace metriko {
    inline MatXd compute_constraint(const Emesh &tm) {
        MatXd M = MatXd::Zero(tm.tquads.size() * 2, tm.tedges.size());
        for (int iq = 0; iq < tm.tquads.size(); iq++) {
            const auto &tq = tm.tquads[iq];
            for (int ih = 0; ih < tq.data.size(); ih++) {
                const int side = tq.data[ih].side;
                const int teid = tm.thalfs[tq.data[ih].thid].teid;
                switch (side) {
                    case 0: { M(iq * 2 + 0, teid) += 1; break; }
                    case 2: { M(iq * 2 + 0, teid) -= 1; break; }
                    case 1: { M(iq * 2 + 1, teid) += 1; break; }
                    case 3: { M(iq * 2 + 1, teid) -= 1; break; }
                    default: throw std::invalid_argument("side has to be 0 ~ 3");
                }
            }
        }
        return M;
    }
}
#endif
