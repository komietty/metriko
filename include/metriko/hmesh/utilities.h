//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_HMESH_UTILITIES_H
#define METRIKO_HMESH_UTILITIES_H
#include "../common/utilities.h"
#include "../common/predicates.h"
#include "hmesh.h"
#include "hmloc.h"

namespace metriko {
inline Row3d conversion_2d_3d(const Face& f, const VecXc& cf, const complex uv) {
    auto [c1, c2, c3] = f.crnrs();
    return conversion_2d_3d(
        cf(c1.id),
        cf(c2.id),
        cf(c3.id),
        c1.vert().pos(),
        c2.vert().pos(),
        c3.vert().pos(),
        uv);
}

inline std::optional<Crnr> try_get_crnr(const Hmesh& hm, int vid, int fid) {
    Face f = hm.faces[fid];
    Vert v = hm.verts[vid];
    for (auto h: f.adjHalfs())
        if (h.crnr().vert() == v) return h.crnr();
    return std::nullopt;
}

inline std::optional<Half> try_get_half(const Hmesh& hm, int eid, int fid) {
    Face f = hm.faces[fid];
    Edge e = hm.edges[eid];
    for (auto h: f.adjHalfs())
        if (h.edge() == e) return h;
    return std::nullopt;
}

inline complex get_face_uv(const HmLoc& loc, int fid, const Hmesh& hm) {
    return std::visit(overloaded {
        [&](const HmLocOnP& f) -> complex { return f.uv; },
        [&](const HmLocOnV& v) -> complex { return hm.cfn[try_get_crnr(hm, v.id, fid).value().id]; },
        [&](const HmLocOnE& e) -> complex {
            auto h = try_get_half(hm, e.id, fid).value();
            return h.lerp_uv(h.isCanonical() ? e.r : 1 - e.r);
        },
        [&](const auto& _) -> complex { METRIKO_FAIL("no impl"); },
    }, loc);
}
}
#endif
