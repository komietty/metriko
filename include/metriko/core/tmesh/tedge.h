//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_TEDGE_H
#define METRIKO_TEDGE_H
#include "motorcycle.h"

namespace metriko {
    class Tmesh;

    class Tvert {
    public:
        const complex uv;
        const Mport *port = nullptr;

        explicit Tvert(
            const complex uv,
            const Mport *port
        ): uv(uv), port(port) {
        }

        bool is_sing() const { return port != nullptr; }
    };

    class Tsgmt {
    public:
        complex fr;
        complex to;
        Face face;
        bool sing; //is from vert singular or not
        CrashType side_fr; // temp
        CrashType side_to; // temp

        explicit Tsgmt(
            complex fr,
            complex to,
            Face face,
            bool sing,
            CrashType side_fr,
            CrashType side_to
        ): fr(fr), to(to), face(face), sing(sing), side_fr(side_fr), side_to(side_to) {
        }

        complex diff() const { return to - fr; }
    };

    class Tedge {
    public:
        const Tmesh *tmesh;
        const int id;
        std::vector<Tsgmt> sgmts;
        const Msgmt seg_fr;
        const Msgmt seg_to;

        Tedge(
            const Tmesh *tmesh,
            const int id,
            const Msgmt &seg_fr,
            const Msgmt &seg_to
        ) : tmesh(tmesh), id(id), seg_fr(seg_fr), seg_to(seg_to) {
            const Msgmt *curr = &seg_fr;
            while (true) {
                sgmts.emplace_back(Tsgmt(curr->fr.uv, curr->to.uv, curr->face, curr->id == 0, curr->fr.side,
                                         curr->to.side));
                if (curr->id == seg_to.id) break;
                curr = &curr->next();
            }

            assert(abs(sgmts.front().fr - seg_fr.fr.uv) < 1e-10);
            assert(abs(sgmts.back().to - seg_to.to.uv) < 1e-10);
        }


        complex uv_fr() const { return sgmts.front().fr; }
        complex uv_to() const { return sgmts.back().to; }
        CrashType side_fr() const { return sgmts.front().side_fr; }
        CrashType side_to() const { return sgmts.back().side_to; }
        //complex uv_fr() const { return seg_fr.fr.uv; }
        //complex uv_to() const { return seg_to.to.uv; }

        //auto segments() const {
        //    return seg_fr.curv->sgmts | vw::filter([&](auto &s) {
        //        return s.id >= seg_fr.id && s.id <= seg_to.id;
        //    });
        //}

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
