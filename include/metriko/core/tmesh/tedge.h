//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_TEDGE_H
#define METRIKO_TEDGE_H
#include "motorcycle.h"

namespace metriko {
    class Tmesh;

    class Tedge {
    public:
        const int id;
        const Tmesh *tmesh;
        const Msgmt seg_fr;
        const Msgmt seg_to;

        Tedge(
            const Tmesh *tmesh,
            const int id,
            const Msgmt &sgfr,
            const Msgmt &sgto
        ) : tmesh(tmesh), id(id), seg_fr(sgfr), seg_to(sgto) { }

        complex uv_fr() const { return seg_fr.fr.uv; }
        complex uv_to() const { return seg_to.to.uv; }
        MvertType type_fr() const { return seg_fr.fr.type; }
        MvertType type_to() const { return seg_to.to.type; }

        auto segments() const {
            return seg_fr.curv->sgmts | vw::filter([&](auto &s) {
                return s.id >= seg_fr.id && s.id <= seg_to.id;
            });
        }

        //auto passing_halfedges() const {
        //    std::vector<Half> halfs;
        //    for (auto &s: segments()) {
        //        auto hfr = s.fr.half;
        //        auto hto = s.to.half;
        //        assert(hfr == std::nullopt);
        //        if (hto != std::nullopt) halfs.emplace_back(hto.value());
        //    }
        //    return halfs;
        //}
    };
}

#endif
