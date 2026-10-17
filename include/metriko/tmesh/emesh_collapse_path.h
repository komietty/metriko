//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_EMESH_COMMON_H
#define METRIKO_EMESH_COMMON_H
#include "metriko/common/loop.h"
#include "metriko/hmesh/hmflat.h"
#include "emesh.h"

namespace metriko {
constexpr double UNFOLD_HAIR = 1e-4;

struct Unfolded {
    vec<int>    nids;  //
    vec<int>    heads; // the k-th thalf of the boundary starts at nids[heads[k]]
    vec<HmFlat> faces; //

    int size()      const { return nids.size(); }
    int next(int i) const { return (i + 1) % size(); }
    int prev(int i) const { return (i + size() - 1) % size(); }
    int at(int k, int nid) const { for (int i = heads[k]; ; i = next(i)) { if (nids[i] == nid) return i; } }
    complex uv(const Emesh& tm, int i) const { return flat_uv(faces[i], tm.tnodes[nids[i]]); }
};

inline Unfolded unfold(const Emesh& tm, const vec<int>& thids) {
    const Hmesh& hm = tm.hm;
    Unfolded c;

    for (int thid: thids) {
        auto& th = tm.thalfs[thid];
        auto  ns = th.tedge().nids;
        if (!th.cano) rg::reverse(ns);
        c.heads.push_back(c.nids.size());
        c.nids.insert(c.nids.end(), ns.begin(), ns.end() - 1);
    }

    int n = c.size();

    for (int k = 0; k < n; ++k) {
        auto& a = tm.tnodes[c.nids[k]];
        auto& b = tm.tnodes[c.nids[c.next(k)]];
        auto  f = try_get_face(hm, a, b);
        METRIKO_CHECK(f, "boundary segment {} -> {} lies in no face", loc_str(a), loc_str(b));
        c.faces.push_back(flat_own(*f));
    }

    for (int k = 1; k < n; ++k) {
        auto& q = c.faces[k - 1];
        auto& p = c.faces[k];
        auto& l = tm.tnodes[c.nids[k]];
        if (p.f == q.f)                               p = q;
        else if (auto* v = std::get_if<HmLocOnV>(&l)) p = flat_around(q, hm.verts[v->id], p.f, false).back();
        else                                          p = flat_across(q, try_get_edge(hm, l, l).value());
    }

    double gap = std::abs(flat_uv(c.faces[n - 1], tm.tnodes[c.nids[0]]) - c.uv(tm, 0));
    METRIKO_CHECK(gap < 1e-6, "the region does not unfold into one plane (closure error {})", gap);
    return c;
}

inline vec<HmLoc> retrace(const Emesh& tm, const Unfolded& c, int a, int b) {
    vec<complex> loop;
    for (int i = 0; i < c.size(); ++i) loop.push_back(c.uv(tm, i));
    auto pts = shortest_path_in_loop(loop, a, b, UNFOLD_HAIR);
    if (pts.empty()) return {};

    auto path = vec{tm.tnodes[c.nids[a]]};
    auto face = c.faces[a];
    for (int s = 0; s + 1 < pts.size(); ++s)
        for (auto& l: flat_trace(face, path.back(), pts[s + 1], s == 0 ? c.faces[c.prev(a)].f : face.f)) path.push_back(l);
    path.back() = tm.tnodes[c.nids[b]];
    return path;
}
}
#endif
