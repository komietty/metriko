//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_HMESH_CLEANUP_H
#define METRIKO_HMESH_CLEANUP_H
#include <set>
#include <igl/unique_edge_map.h>
#include <igl/edge_flaps.h>
#include <igl/edge_lengths.h>
#include <igl/flip_edge.h>
#include <igl/collapse_edge.h>
#include <igl/circulation.h>
#include <igl/remove_unreferenced.h>
#include <igl/per_vertex_normals.h>
#include "metriko/common/typedef.h"

// removes the slivers and the tiny edges of the input before the pipeline. the conditions (4) of bommes et al. 2013
// keep every triangle a margin away from degeneracy, which a zero-area or needle triangle makes meaningless. not part
// of the paper: the surface may move. each round
//   1. flips the long edge of every cap (a sliver with a corner above 120 deg)
//   2. collapses the short edge of every needle and every tiny edge
//   3. if neither changed anything, flips any edge of a sliver, and failing that splits the slivers
namespace metriko {
namespace cleanup_detail {
    constexpr double sliver_deg = 1.;     // a face whose smallest angle is below this is a sliver
    constexpr double tiny_ratio = 0.01;   // an edge shorter than this fraction of the median edge is tiny

    inline std::pair<int, int> key(const int a, const int b) { return a < b ? std::pair(a, b) : std::pair(b, a); }

    // corner angles in degrees. atan2 gives 0 at a zero-length edge where igl::internal_angles gives nan
    inline std::array<double, 3> angles(const Row3d& a, const Row3d& b, const Row3d& c) {
        auto at = [](const Row3d& o, const Row3d& p, const Row3d& q) {
            const Row3d u = p - o, w = q - o;
            return std::atan2(u.cross(w).norm(), u.dot(w)) * 180. / PI;
        };
        return { at(a, b, c), at(b, c, a), at(c, a, b) };
    }
    inline std::array<double, 3> angles(const MatXd& V, const MatXi& F, const int f) { return angles(V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2))); }
    inline double min_angle(const Row3d& a, const Row3d& b, const Row3d& c) { return rg::min(angles(a, b, c)); }
    inline double min_angle(const MatXd& V, const MatXi& F, const int f) { return rg::min(angles(V, F, f)); }

    inline double median_edge(const MatXd& V, const MatXi& F) {
        MatXi E, uE;
        VecXi EMAP;
        igl::unique_edge_map(F, E, uE, EMAP);
        vec<double> len(uE.rows());
        for (int e = 0; e < uE.rows(); ++e) len[e] = (V.row(uE(e, 0)) - V.row(uE(e, 1))).norm();
        rg::nth_element(len, len.begin() + len.size() / 2);
        return len[len.size() / 2];
    }

    inline bool clean(const MatXd& V, const MatXi& F, const double len_med) {
        MatXd L;
        igl::edge_lengths(V, F, L);
        for (int f = 0; f < F.rows(); ++f) if (min_angle(V, F, f) < sliver_deg) return false;
        return L.minCoeff() >= tiny_ratio * len_med;
    }

    // flips one edge of every sliver: the edge facing its largest angle when that is above 120 deg (a needle is left
    // to the collapses), or with any = true whichever edge is best. the smallest angle of the pair must grow and both
    // faces keep their orientation. a face is flipped at most once per call
    inline int flip_slivers(const MatXd& V, MatXi& F, const bool any) {
        MatXi E, uE;
        VecXi EMAP;
        vec<vec<int>> uE2E;
        igl::unique_edge_map(F, E, uE, EMAP, uE2E);
        MatXd N;
        igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, N);
        vec<int> val(V.rows(), 0);
        std::set<std::pair<int, int>> present;   // a flip must not create an edge that exists
        for (int e = 0; e < uE.rows(); ++e) { ++val[uE(e, 0)]; ++val[uE(e, 1)]; present.insert(key(uE(e, 0), uE(e, 1))); }

        const int m = F.rows();
        vec<bool> dirty(m, false);
        int n = 0;
        for (int f = 0; f < m; ++f) {
            if (dirty[f] || min_angle(V, F, f) >= sliver_deg) continue;
            const auto K = angles(V, F, f);
            const int cmax = rg::max_element(K) - K.begin();
            if (!any && K[cmax] <= 120.) continue;
            double best = -1;
            int best_ue = -1, best_f2 = -1;
            for (int c = 0; c < 3; ++c) {
                if (!any && c != cmax) continue;
                const int ue = EMAP(f + c * m);
                if (uE2E[ue].size() != 2) continue;
                const int h = uE2E[ue][0] % m == f ? uE2E[ue][1] : uE2E[ue][0];
                const int f2 = h % m;
                if (dirty[f2]) continue;
                // f = (v4, v1, v2) and f2 = (v3, v2, v1) become (v1, v3, v4) and (v2, v4, v3)
                const int v1 = F(f, (c + 1) % 3), v2 = F(f, (c + 2) % 3), v4 = F(f, c), v3 = F(f2, h / m);
                if (val[v1] <= 3 || val[v2] <= 3 || present.contains(key(v3, v4))) continue;
                const Row3d p1 = V.row(v1), p2 = V.row(v2), p3 = V.row(v3), p4 = V.row(v4);
                const double before = std::min(min_angle(V, F, f), min_angle(V, F, f2));
                const double after  = std::min(min_angle(p1, p3, p4), min_angle(p2, p4, p3));
                if (after <= before || after <= best) continue;
                const Row3d nn = N.row(v1) + N.row(v2) + N.row(v3) + N.row(v4);
                if (nn.norm() == 0) continue;
                if ((p3 - p1).cross(p4 - p1).dot(nn) <= 0 || (p4 - p2).cross(p3 - p2).dot(nn) <= 0) continue;
                best = after; best_ue = ue; best_f2 = f2;
            }
            if (best_ue < 0) continue;
            const int v1 = uE(best_ue, 0), v2 = uE(best_ue, 1);
            igl::flip_edge(F, E, uE, EMAP, uE2E, best_ue);
            const int v3 = uE(best_ue, 0), v4 = uE(best_ue, 1);
            present.erase(key(v1, v2));
            present.insert(key(v3, v4));
            --val[v1]; --val[v2]; ++val[v3]; ++val[v4];
            dirty[f] = dirty[best_f2] = true;
            ++n;
        }
        return n;
    }

    // collapses the shortest edge of every needle (a sliver whose shortest edge is under a quarter of the next one)
    // and every tiny edge, to whichever of its ends or its midpoint leaves the best smallest angle around it. no face
    // may flip; a needle must get better, and a tiny edge must not end below both the sliver bound and where it was
    inline int collapse_short(MatXd& V, MatXi& F, const double len_med) {
        MatXi uE, EF, EI;
        VecXi EMAP;
        igl::edge_flaps(F, uE, EMAP, EF, EI);
        MatXd N, L;
        igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, N);
        igl::edge_lengths(V, F, L);   // L(f, c) is the edge facing corner c, as EMAP(f + c * m)
        const int m = F.rows();
        const double tiny = tiny_ratio * len_med;

        vec<std::pair<double, int>> cand;   // (length, edge), shortest first
        for (int f = 0; f < m; ++f) {
            std::array<std::pair<double, int>, 3> l;
            for (int c = 0; c < 3; ++c) {
                l[c] = { L(f, c), EMAP(f + c * m) };
                if (l[c].first < tiny) cand.push_back(l[c]);
            }
            rg::sort(l);
            if (min_angle(V, F, f) < sliver_deg && l[0].first <= 0.25 * l[1].first) cand.push_back(l[0]);
        }
        rg::sort(cand);
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

        int n = 0;
        for (const auto& [len, e]: cand) {
            if (EF(e, 0) < 0 || EF(e, 1) < 0) continue;   // boundary: moving it would change the outline
            const int a = uE(e, 0), b = uE(e, 1);
            if (a == b) continue;                          // an earlier collapse already removed it
            vec<int> Nsv, Nsf, Ndv, Ndf;
            igl::circulation(e, true,  F, EMAP, EF, EI, Nsv, Nsf);
            igl::circulation(e, false, F, EMAP, EF, EI, Ndv, Ndf);
            vec<int> around = Nsf;
            around.insert(around.end(), Ndf.begin(), Ndf.end());
            const Row3d nref = (N.row(a) + N.row(b)).normalized();
            const bool is_tiny = (V.row(a) - V.row(b)).norm() < tiny;
            double before = 180.;
            for (int f: around) before = std::min(before, min_angle(V, F, f));

            // the smallest angle around the edge once a and b move to x, or -1 if a face would flip
            auto quality = [&](const Row3d& x) {
                double q = 180.;
                for (int f: around) {
                    if (f == EF(e, 0) || f == EF(e, 1)) continue;   // these two faces disappear
                    Row3d p[3], r[3];
                    for (int c = 0; c < 3; ++c) { const int v = F(f, c); p[c] = V.row(v); r[c] = v == a || v == b ? x : p[c]; }
                    const Row3d n0 = (p[1] - p[0]).cross(p[2] - p[0]), n1 = (r[1] - r[0]).cross(r[2] - r[0]);
                    // a sliver has no usable normal: the vertex normals stand in for it
                    if (n1.norm() == 0 || n1.dot(min_angle(p[0], p[1], p[2]) >= sliver_deg ? n0 : nref) <= 0) return -1.;
                    q = std::min(q, min_angle(r[0], r[1], r[2]));
                }
                return q;
            };

            double best_q = -1;
            Row3d best;
            for (const Row3d x: { Row3d(V.row(a)), Row3d(V.row(b)), Row3d(0.5 * (V.row(a) + V.row(b))) }) {
                const double q = quality(x);
                if (q < 0 || (is_tiny ? q < std::min(before, sliver_deg) : q <= before)) continue;
                if (q > best_q) { best_q = q; best = x; }
            }
            if (best_q < 0) continue;
            int e1, e2, f1, f2;   // collapse_edge also refuses the ones that break the link condition
            if (igl::collapse_edge(e, Eigen::RowVectorXd(best), Nsv, Nsf, Ndv, Ndf, V, F, uE, EMAP, EF, EI, e1, e2, f1, f2)) ++n;
        }
        if (n == 0) return 0;

        // collapse_edge blanks the faces it removes: drop them, then the vertices left orphaned
        MatXi G(F.rows(), 3);
        int g = 0;
        for (int f = 0; f < F.rows(); ++f) if (F(f, 0) != F(f, 1) && F(f, 1) != F(f, 2) && F(f, 2) != F(f, 0)) G.row(g++) = F.row(f);
        G.conservativeResize(g, 3);
        MatXd NV;
        MatXi NF;
        VecXi I;
        igl::remove_unreferenced(V, G, NV, NF, I);
        V = NV;
        F = NF;
        return n;
    }

    // splits the longest edge of every sliver, at the foot of the apex for a cap (a corner above 150 deg) and at the
    // midpoint otherwise. the new vertex lies on the edge, so the surface does not move. at most budget splits in all
    inline int split_slivers(MatXd& V, MatXi& F, const double min_len, int& budget) {
        MatXi E, uE;
        VecXi EMAP;
        vec<vec<int>> uE2E;
        igl::unique_edge_map(F, E, uE, EMAP, uE2E);
        MatXd L;
        igl::edge_lengths(V, F, L);
        const int m = F.rows();
        vec<bool> dirty(m, false);
        vec<Row3d> addV;
        vec<Row3i> addF;
        for (int f = 0; f < m && budget > 0; ++f) {
            if (dirty[f] || min_angle(V, F, f) >= sliver_deg) continue;
            int c = 0;   // the longest edge, facing corner c
            for (int i = 1; i < 3; ++i) if (L(f, i) > L(f, c)) c = i;
            if (L(f, c) < min_len) continue;
            const int ue = EMAP(f + c * m);
            if (uE2E[ue].size() > 2 || rg::any_of(uE2E[ue], [&](int h) { return dirty[h % m]; })) continue;
            const Row3d pa = V.row(F(f, (c + 1) % 3)), pb = V.row(F(f, (c + 2) % 3)), pc = V.row(F(f, c));
            const double t = rg::max(angles(V, F, f)) > 150. ? std::clamp((pc - pa).dot(pb - pa) / (pb - pa).squaredNorm(), 0.01, 0.99) : 0.5;
            const int nv = V.rows() + addV.size();
            addV.emplace_back((1 - t) * pa + t * pb);
            for (int h: uE2E[ue]) {
                const int g = h % m, k = h / m;   // face g, the edge facing its corner k
                const int w = F(g, k), x = F(g, (k + 1) % 3), y = F(g, (k + 2) % 3);
                F.row(g) << w, x, nv;
                addF.emplace_back(w, nv, y);
                dirty[g] = true;
            }
            --budget;
        }
        const int nV = V.rows(), nF = F.rows();
        V.conservativeResize(nV + addV.size(), 3);
        F.conservativeResize(nF + addF.size(), 3);
        for (size_t i = 0; i < addV.size(); ++i) V.row(nV + i) = addV[i];
        for (size_t i = 0; i < addF.size(); ++i) F.row(nF + i) = addF[i];
        return addV.size();
    }
}

inline void cleanup_mesh(MatXd& V, MatXi& F) {
    using namespace cleanup_detail;
    const double len_med = median_edge(V, F);
    int budget = F.rows();   // new vertices at most as many as the input faces
    // applies op until it changes nothing, at most 50 times, and returns the number of changes
    auto repeat = [](auto op) { int n = 0; for (int k = 0; k < 50; ++k) { const int c = op(); if (c == 0) break; n += c; } return n; };
    for (int it = 0; it < 200; ++it) {
        int n = repeat([&] { return flip_slivers(V, F, false); });
        n += repeat([&] { return collapse_short(V, F, len_med); });
        if (clean(V, F, len_med)) break;
        // stuck: flip any edge of the slivers left, and only then split
        if (n == 0 && flip_slivers(V, F, true) == 0 && split_slivers(V, F, 1e-3 * len_med, budget) == 0) break;
    }
}
}
#endif
