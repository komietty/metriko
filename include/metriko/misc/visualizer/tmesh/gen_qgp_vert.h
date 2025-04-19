//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_GEN_QGP_VERT_H
#define METRIKO_GEN_QGP_VERT_H
#include "metriko/core/tmesh/motorcycle.h"
#include "metriko/core/tmesh/tmesh.h"

namespace metriko::visualizer {

    struct SplitVert {
        Msgmt segment;
        complex uv;
    };

    struct SplitElem {
        Msgmt segment;
        complex uv;
        double distance;
    };

    inline std::vector<SplitVert> construct_verts_on_tedge(
        const Tmesh& tm,
        const VecXd& X,
        const VecXd& R,
        const int teid
    ) {
        //std::cout << X << std::endl;
        //assert((X.array() > 0).all());
        // Using original intrinsic length. It might be better using extrinsic length
        std::vector<SplitVert> verts;
        auto sum = R[teid];
        auto num = std::round(X[teid]);


        { // add first vert
            auto sg = tm.tedges[teid].seg_fr;
            auto uv = sg.fr.uv;
            verts.emplace_back(SplitVert{sg, uv});
        }

        if (num > 1) {
            Msgmt sg = tm.tedges[teid].seg_fr;
            complex uv = sg.fr.uv;
            double accum = 1e-3; // todo: fix
            double total = 0.;
            int count = 1;
            double div = sum / num;
            do {
                auto d = sg.to.uv - uv;
                auto l = abs(d);
                if (accum + l < div) {
                    accum += l;
                    if (sg.next_id == -1) break;
                    sg = sg.next();
                    uv = sg.fr.uv;
                } else {
                    complex _uv = uv + d / l * (div - accum);
                    verts.emplace_back(SplitVert{sg, _uv});
                    accum = 0;
                    total += div;
                    count++;
                    uv = _uv;
                }
            } while (total < R[teid] && count < num);
            assert(count == num);
        }
        if (num > 0) { // add last vert
            auto sg = tm.tedges[teid].seg_to;
            auto uv = sg.to.uv;
            verts.emplace_back(SplitVert{sg, uv});
        }

        return  verts;
    }
}

#endif
