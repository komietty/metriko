#ifndef METRIKO_HMLOC_H
#define METRIKO_HMLOC_H
#include "metriko/core/common/utilities.h"
#include "metriko/core/hmesh/hmesh.h"

namespace metriko {

struct HmLocOnV { int id;             bool operator==(const HmLocOnV&) const = default; };
struct HmLocOnC { int id;             bool operator==(const HmLocOnC&) const = default; };
struct HmLocOnE { int id; double r;   bool operator==(const HmLocOnE&) const = default; };
struct HmLocOnH { int id; double r;   bool operator==(const HmLocOnH&) const = default; };
struct HmLocOnF { int id; complex xy; bool operator==(const HmLocOnF&) const = default; };
struct HmLocOnP { int id; complex uv; bool operator==(const HmLocOnP&) const = default; }; // parameter space coord in face

using HmLoc = std::variant<HmLocOnV, HmLocOnE, HmLocOnF, HmLocOnH, HmLocOnC, HmLocOnP>;

inline Row3d get_ptloc_pos(
    const Hmesh& hm,
    const HmLoc& hl
) {
    return std::visit(overloaded{
        [&](const HmLocOnV& l) -> Row3d { return hm.verts[l.id].pos(); },
        [&](const HmLocOnC& l) -> Row3d { return hm.crnrs[l.id].vert().pos(); },
        [&](const HmLocOnE& l) -> Row3d {
            Edge e = hm.edges[l.id];
            Row3d p0 = e.vert0().pos();
            Row3d p1 = e.vert1().pos();
            return p0 * (1 - l.r) + p1 * l.r;
        },
        [&](const HmLocOnH& l) -> Row3d {
            Half h = hm.halfs[l.id];
            Row3d p0 = h.tail().pos();
            Row3d p1 = h.head().pos();
            return p0 * (1 - l.r) + p1 * l.r;
        },
        [&](const HmLocOnF& l) -> Row3d {
            Face f = hm.faces[l.id];
            Row3d o = f.half().tail().pos();
            Row3d x = f.basisX();
            Row3d y = f.basisY();
            return o + x * l.xy.real() + y * l.xy.imag();
        },
        [&](const HmLocOnP& _) -> Row3d { throw std::runtime_error("no impl"); },
    }, hl);
}

inline std::optional<Crnr> try_get_crnr(const Hmesh& hm, int vid, int fid) {
    Face f = hm.faces[fid];
    Vert v = hm.verts[vid];
    for (auto h: f.adjHalfs()) if (h.crnr().vert() == v) return h.crnr();
    return std::nullopt;
}

inline std::optional<Half> try_get_half(const Hmesh& hm, int eid, int fid) {
    Face f = hm.faces[fid];
    Edge e = hm.edges[eid];
    for (auto h: f.adjHalfs()) if (h.edge() == e) return h;
    return std::nullopt;
}

}
#endif
