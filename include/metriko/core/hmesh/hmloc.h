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

inline Row3d get_ptloc_pos(const Hmesh& hm, const HmLoc& loc) {
    return std::visit(overloaded{
        [&](const HmLocOnV& l) -> Row3d { return hm.verts[l.id].pos(); },
        [&](const HmLocOnC& l) -> Row3d { return hm.crnrs[l.id].vert().pos(); },
        [&](const HmLocOnE& l) -> Row3d { return hm.edges[l.id].lerp(l.r); },
        [&](const HmLocOnH& l) -> Row3d { return hm.halfs[l.id].lerp(l.r); },
        [&](const HmLocOnF& l) -> Row3d {
            Face f = hm.faces[l.id];
            Row3d o = f.half().tail().pos();
            Row3d x = f.basisX();
            Row3d y = f.basisY();
            return o + x * l.xy.real() + y * l.xy.imag();
        },
        [&](const HmLocOnP& _) -> Row3d { throw std::runtime_error("no impl"); },
    }, loc);
}

inline Row3d get_ptloc_nml(const Hmesh& hm, const HmLoc& loc) {
    return std::visit(overloaded{
        [&](const HmLocOnV& l) -> Row3d { return hm.verts[l.id].normal(); },
        [&](const HmLocOnF& l) -> Row3d { return hm.faces[l.id].normal(); },
        [&](const HmLocOnE& l) -> Row3d {
            Edge  e = hm.edges[l.id];
            Row3d d = e.face0().normal() + e.face1().normal();
            return d.norm() > 0 ? d.normalized() : Row3d::Zero();
        },
        [&](const HmLocOnH& l) -> Row3d {
            Half  h = hm.halfs[l.id];
            Row3d d = h.face().normal() + h.twin().face().normal();
            return d.norm() > 0 ? d.normalized() : Row3d::Zero();
        },
        [&](const auto& _) -> Row3d { throw std::runtime_error("no impl"); },
    }, loc);
}

inline vec<int> get_ptloc_faces(const Hmesh& hm, const HmLoc& loc) {
    vec<int> out = {};
    std::visit(overloaded{
        [&](const HmLocOnV& l) { for (Half h : hm.verts[l.id].adjHalfs()) out.push_back(h.face().id); },
        [&](const HmLocOnE& l) { out.push_back(hm.edges[l.id].face0().id); out.push_back(hm.edges[l.id].face1().id); },
        [&](const HmLocOnH& l) { Half h = hm.halfs[l.id]; out.push_back(h.face().id); out.push_back(h.twin().face().id) ; },
        [&](const HmLocOnF& l) { out.push_back(l.id); },
        [&](const auto&) { },
    }, loc);
    std::erase(out, -1);
    return out;
}

inline vec<int> get_ptloc_edges(const Hmesh& hm, const HmLoc& loc) {
    vec<int> out = {};
    std::visit(overloaded{
        [&](const HmLocOnV& l) { for (Half h : hm.verts[l.id].adjHalfs()) out.push_back(h.edge().id); },
        [&](const HmLocOnE& l) { out.push_back(hm.edges[l.id].id); },
        [&](const HmLocOnH& l) { out.push_back(hm.halfs[l.id].edge().id); },
        [&](const auto&) { },
    }, loc);
    return out;
}

inline std::optional<Face> try_get_face(const Hmesh& hm, const HmLoc& a, const HmLoc& b) {
    auto fa = get_ptloc_faces(hm, a);
    auto fb = get_ptloc_faces(hm, b);
    for (int x : fa)
    for (int y : fb)
        if (x == y) return hm.faces[x];
    return std::nullopt;
}

inline std::optional<Edge> try_get_edge(const Hmesh& hm, const HmLoc& a, const HmLoc& b) {
    auto ea = get_ptloc_edges(hm, a);
    auto eb = get_ptloc_edges(hm, b);
    for (int x : ea)
    for (int y : eb)
        if (x == y) return hm.edges[x];
    return std::nullopt;
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

inline std::optional<std::pair<Half, double>> try_get_ratio(const Hmesh& hm, const HmLoc& l) {
    return std::visit(overloaded{
        [&](const HmLocOnE& l) -> std::optional<std::pair<Half, double>> {
            auto e = hm.edges[l.id];
            auto h = e.half();
            auto t = h.tail() == e.vert0() ? l.r : 1. - l.r;
            return std::pair(h, 1. - t);
        },
        [&](const HmLocOnH& l) -> std::optional<std::pair<Half, double>> {
            return std::pair(hm.halfs[l.id], 1. - l.r); // todo: seeems the ratio is flipped...
        },
        [&](const auto&) -> std::optional<std::pair<Half, double>> { return std::nullopt; },
    }, l);
}

inline bool is_in_face(Face face, const HmLoc& l) {
    return std::visit(overloaded{
        [&](const HmLocOnV& v) { for (Half h_: face.adjHalfs()) { if (h_.tail().id == v.id) return true; } return false; },
        [&](const HmLocOnE& e) { for (Half h_: face.adjHalfs()) { if (h_.edge().id == e.id) return true; } return false; },
        [&](const HmLocOnH& h) { for (Half h_: face.adjHalfs()) { if (h_.id == h.id)        return true; } return false; },
        [&](const HmLocOnF& f) { return f.id == face.id; },
        [&](const auto&) -> bool { throw std::runtime_error("not implemented"); },
    }, l);
};

inline auto loc_str(const HmLoc& l) {
    return std::visit(overloaded{
        [](const HmLocOnV& v){ return std::format("V(id={})", v.id); },
        [](const HmLocOnC& c){ return std::format("C(id={})", c.id); },
        [](const HmLocOnE& e){ return std::format("E(id={}, r={})", e.id, e.r); },
        [](const HmLocOnH& h){ return std::format("H(id={}, r={})", h.id, h.r); },
        [](const HmLocOnF& f){ return std::format("F(id={}, xy=({},{}))", f.id, f.xy.real(), f.xy.imag()); },
        [](const HmLocOnP& p){ return std::format("P(id={}, uv=({},{}))", p.id, p.uv.real(), p.uv.imag()); },
    }, l);
}
}
#endif
