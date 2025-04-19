//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_CONSTRAINT_H
#define METRIKO_QUANTIZATION_CONSTRAINT_H
#include "metriko/core/tmesh/tmesh.h"

namespace metriko {
    inline MatXd compute_constraint(const Tmesh &tmesh) {
        MatXd M = MatXd::Zero(tmesh.nTQ * 2, tmesh.nTE);
        for (int iq = 0; iq < tmesh.nTQ; iq++) {
            const auto &tq = tmesh.tquads[iq];
            for (int ih = 0; ih < tq.thids.size(); ih++) {
                const auto &th = tmesh.thalfs[tq.thids[ih]];
                const int side = tq.sides[ih];
                const int teid = th.edge().id;
                if (side == 0) M(iq * 2 + 0, teid) = 1;
                else if (side == 2) M(iq * 2 + 0, teid) = -1;
                else if (side == 1) M(iq * 2 + 1, teid) = 1;
                else if (side == 3) M(iq * 2 + 1, teid) = -1;
                else throw std::invalid_argument("eval has to be 0 ~ 3");
            }
        }
        return M;
    }
}

#endif
