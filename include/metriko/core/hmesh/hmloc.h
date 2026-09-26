//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
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
        [&](const HmLocOnF& l) -> Row3d { return hm.faces[l.id].to_world(l.xy); },
        [&](const HmLocOnP& _) -> Row3d { throw std::runtime_error("no impl"); },
    }, loc);
}

inline Row3d get_ptloc_nml(const Hmesh& hm, const HmLoc& loc) {
    return std::visit(overloaded{
        [&](const HmLocOnV& l) -> Row3d { return hm.verts[l.id].normal(); },
        [&](const HmLocOnF& l) -> Row3d { return hm.faces[l.id].normal(); },
        [&](const HmLocOnE& l) -> Row3d { return hm.edges[l.id].nml(); },
        [&](const HmLocOnH& l) -> Row3d { return hm.halfs[l.id].nml(); },
        [&](const auto& _) -> Row3d { throw std::runtime_error("no impl"); },
    }, loc);
}

inline vec<int> get_ptloc_verts(const Hmesh& hm, const HmLoc& loc) {
    return std::visit(overloaded{
        [&](const HmLocOnV& l) -> vec<int> { return {l.id}; },
        [&](const HmLocOnE& l) -> vec<int> { auto e = hm.edges[l.id]; return {e.vert0().id, e.vert1().id}; },
        [&](const HmLocOnH& l) -> vec<int> { auto h = hm.halfs[l.id]; return {h.tail().id, h.head().id}; },
        [&](const HmLocOnF& l) -> vec<int> { auto [a, b, c] = hm.faces[l.id].verts(); return {a.id, b.id, c.id}; },
        [&](const auto&)       -> vec<int> { return {}; },
    }, loc);
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

inline std::optional<std::pair<Half, double>> try_get_half_ratio(const Hmesh& hm, const HmLoc& loc) {
    using R = std::optional<std::pair<Half, double>>;
    return std::visit(overloaded{
        [&](const HmLocOnE& l) -> R { auto e = hm.edges[l.id]; auto h = e.half(); return std::pair(h, 1. - (h.isCanonical() ? l.r : 1. - l.r)); },
        [&](const HmLocOnH& l) -> R { return std::pair(hm.halfs[l.id], 1. - l.r); },
        [&](const auto&)       -> R { return std::nullopt; },
    }, loc);
}

inline std::optional<std::pair<Edge, double>> try_get_edge_ratio(const Hmesh& hm, const HmLoc& loc) {
    using R = std::optional<std::pair<Edge, double>>;
    return std::visit(overloaded{
        [&](const HmLocOnE& l) -> R { return std::pair{hm.edges[l.id], l.r}; },
        [&](const HmLocOnH& l) -> R { Half h = hm.halfs[l.id]; return std::pair{h.edge(), h.isCanonical() ? l.r : 1. - l.r}; },
        [&](const auto&)       -> R { return std::nullopt; },
    }, loc);
}

inline bool is_in_face(Face f, const HmLoc& loc) {
    return std::visit(overloaded{
        [&](const HmLocOnV& l) { for (Half h: f.adjHalfs()) { if (h.tail().id == l.id) return true; } return false; },
        [&](const HmLocOnE& l) { for (Half h: f.adjHalfs()) { if (h.edge().id == l.id) return true; } return false; },
        [&](const HmLocOnH& l) { for (Half h: f.adjHalfs()) { if (h.id == l.id)        return true; } return false; },
        [&](const HmLocOnF& l) { return l.id == f.id; },
        [&](const auto&) -> bool { throw std::runtime_error("not implemented"); },
    }, loc);
}

inline bool find_strict_intersection(
    const Face f,
    const HmLoc& la,
    const HmLoc& lb,
    const HmLoc& lc,
    const HmLoc& ld
) {
    if (!is_in_face(f, la) || !is_in_face(f, lb) || !is_in_face(f, lc) || !is_in_face(f, ld)) return false;
    return find_strict_intersection(
        f.to_local(get_ptloc_pos(*f.m, la)),
        f.to_local(get_ptloc_pos(*f.m, lb)),
        f.to_local(get_ptloc_pos(*f.m, lc)),
        f.to_local(get_ptloc_pos(*f.m, ld))
    );
}

inline bool is_in_star(Vert v, const HmLoc& loc) {
    return rg::any_of(v.adjHalfs(), [&](Half h) { return is_in_face(h.face(), loc); });
}

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
