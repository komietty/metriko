#ifndef METRIKO_EXAMPLE_CLEANUP_H
#define METRIKO_EXAMPLE_CLEANUP_H
#include <print>
#include <set>
#include <igl/unique_edge_map.h>
#include <igl/flip_edge.h>
#include <igl/edge_flaps.h>
#include <igl/collapse_edge.h>
#include <igl/circulation.h>
#include <igl/remove_unreferenced.h>
#include <igl/per_vertex_normals.h>
#include <igl/qslim.h>
#include "metriko/core/common/utilities.h"

// Extrinsic cleanup of the input triangulation, to run before the quad pipeline.
//   delaunay_flips       removes caps    (a corner close to 180 deg) by flipping non-Delaunay edges
//   collapse_short_edges removes needles (one very short edge)       by collapsing them
// A flip cannot remove a needle and a collapse cannot remove a cap, so both passes are needed.
//
// Everything here is extrinsic: the flip replaces the diagonal of the quad formed by two
// triangles, which moves the surface unless the quad is flat. The dihedral guard is what keeps
// sharp features intact, so it must stay tight (20-30 deg). An intrinsic Delaunay triangulation
// would preserve the surface exactly, but its edges are geodesics and cannot be stored as (V, F).
namespace metriko::cleanup {

namespace detail {
    inline Row3d tri_normal(const MatXd& V, const int a, const int b, const int c) {
        Row3d p = V.row(a), q = V.row(b), r = V.row(c);
        return (q - p).cross(r - p);   // length is twice the area
    }

    // cotangent of the angle at corner c of face f, i.e. the angle opposite the edge (c+1, c+2).
    // an edge is locally Delaunay iff the two angles opposite it have cot sum >= 0
    inline double cot_opposite(const MatXd& V, const MatXi& F, const int f, const int c) {
        Row3d o = V.row(F(f, c)), p = V.row(F(f, (c + 1) % 3)), q = V.row(F(f, (c + 2) % 3));
        Row3d u = p - o, v = q - o;
        double s = u.cross(v).norm();
        return s > 0 ? u.dot(v) / s : 0.;
    }

    inline std::pair<int, int> key(const int a, const int b) { return a < b ? std::pair(a, b) : std::pair(b, a); }

    // shape quality: area over the square of the longest edge. an equilateral triangle sits at
    // 0.433 and a sliver tends to 0. below ~0.05 the face normal carries no usable direction
    inline double quality(const Row3d& a, const Row3d& b, const Row3d& c) {
        const double s = (b - a).cross(c - a).norm() / 2.;
        const double l = std::max({(b - a).norm(), (c - b).norm(), (a - c).norm()});
        return l > 0 ? s / (l * l) : 0.;
    }
    inline double quality(const MatXd& V, const MatXi& F, const int f) {
        return quality(Row3d(V.row(F(f, 0))), Row3d(V.row(F(f, 1))), Row3d(V.row(F(f, 2))));
    }
}

struct Stats {
    int    nV = 0, nF = 0, nE = 0;
    double len_min = 0, len_med = 0, len_max = 0;
    double ang_min = 0, ang_p05 = 0, ang_med = 0;   // smallest angle of a triangle, in degrees
    int    n_below_10 = 0;                          // triangles whose smallest angle is under 10 deg
    int    n_non_delaunay = 0;
    int    n_sharp = 0;                             // edges whose dihedral angle exceeds 30 deg
};

inline Stats measure(const MatXd& V, const MatXi& F) {
    Stats s;
    s.nV = (int)V.rows();
    s.nF = (int)F.rows();

    MatXi E, uE;
    VecXi EMAP;
    vec<vec<int>> uE2E;
    igl::unique_edge_map(F, E, uE, EMAP, uE2E);
    s.nE = (int)uE.rows();

    vec<double> len;
    len.reserve(s.nE);
    for (int e = 0; e < s.nE; ++e) len.push_back((V.row(uE(e, 0)) - V.row(uE(e, 1))).norm());
    rg::sort(len);
    if (!len.empty()) { s.len_min = len.front(); s.len_med = len[len.size() / 2]; s.len_max = len.back(); }

    vec<double> ang;
    ang.reserve(s.nF);
    for (int f = 0; f < s.nF; ++f) {
        double lo = 180.;
        for (int c = 0; c < 3; ++c) {
            Row3d o = V.row(F(f, c)), p = V.row(F(f, (c + 1) % 3)), q = V.row(F(f, (c + 2) % 3));
            Row3d u = p - o, v = q - o;
            lo = std::min(lo, std::atan2(u.cross(v).norm(), u.dot(v)) * 180. / PI);
        }
        ang.push_back(lo);
        if (lo < 10.) ++s.n_below_10;
    }
    rg::sort(ang);
    if (!ang.empty()) { s.ang_min = ang.front(); s.ang_p05 = ang[ang.size() / 20]; s.ang_med = ang[ang.size() / 2]; }

    constexpr double sharp_cos = 0.8660254037844387;   // cos(30 deg): the feature edges to preserve
    for (int ue = 0; ue < (int)uE2E.size(); ++ue) {
        if (uE2E[ue].size() != 2) continue;   // boundary edges are Delaunay by definition
        const int f1 = uE2E[ue][0] % s.nF, c1 = uE2E[ue][0] / s.nF;
        const int f2 = uE2E[ue][1] % s.nF, c2 = uE2E[ue][1] / s.nF;
        if (detail::cot_opposite(V, F, f1, c1) + detail::cot_opposite(V, F, f2, c2) < 0) ++s.n_non_delaunay;
        Row3d n1 = detail::tri_normal(V, F(f1, 0), F(f1, 1), F(f1, 2));
        Row3d n2 = detail::tri_normal(V, F(f2, 0), F(f2, 1), F(f2, 2));
        const double a1 = n1.norm(), a2 = n2.norm();
        if (a1 > 0 && a2 > 0 && n1.dot(n2) / (a1 * a2) < sharp_cos) ++s.n_sharp;
    }
    return s;
}

inline void print_stats(const char* tag, const Stats& s) {
    std::println("[cleanup] {:<6} nV {:>7} nF {:>7} | edge {:.5f} {:.5f} {:.5f} (min med max) | min angle {:.1f} {:.1f} {:.1f} deg (min p05 med), <10deg {:>5} | non-Delaunay {:>6} ({:.1f}%) | sharp {:>5}",
                 tag, s.nV, s.nF, s.len_min, s.len_med, s.len_max,
                 s.ang_min, s.ang_p05, s.ang_med, s.n_below_10,
                 s.n_non_delaunay, s.nE ? 100. * s.n_non_delaunay / s.nE : 0., s.n_sharp);
}

// flip every non-Delaunay edge whose two triangles are nearly coplanar. vertex positions never
// move. a pair containing a sliver is exempt from the coplanarity guard: a sliver has no usable
// normal, so the large dihedral it shows is not a feature but the defect itself. returns the
// number of flips applied
inline int delaunay_flips(const MatXd& V, MatXi& F, const double feature_angle_deg = 20.,
                          const double sliver_quality = 0.05, const int max_pass = 20) {
    const double cos_thr = std::cos(feature_angle_deg * PI / 180.);
    int total = 0;

    for (int pass = 0; pass < max_pass; ++pass) {
        MatXi E, uE;
        VecXi EMAP;
        vec<vec<int>> uE2E;
        igl::unique_edge_map(F, E, uE, EMAP, uE2E);

        // area-weighted vertex normals: the reference the orientation test uses, chosen because
        // slivers carry no area and so cannot corrupt it
        MatXd N;
        igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, N);

        // a flip whose new edge already exists makes the mesh non-manifold, and a flip that
        // takes the last-but-two edge off a vertex disconnects it: both need the counts below
        vec<int> val(V.rows(), 0);
        std::set<std::pair<int, int>> present;
        for (int e = 0; e < uE.rows(); ++e) {
            ++val[uE(e, 0)];
            ++val[uE(e, 1)];
            present.insert(detail::key(uE(e, 0), uE(e, 1)));
        }

        const int m = (int)F.rows();
        vec<bool> dirty(m, false);
        int n = 0;

        for (int ue = 0; ue < (int)uE2E.size(); ++ue) {
            if (uE2E[ue].size() != 2) continue;                      // boundary or non-manifold
            const int f1 = uE2E[ue][0] % m, c1 = uE2E[ue][0] / m;
            const int f2 = uE2E[ue][1] % m, c2 = uE2E[ue][1] / m;
            if (dirty[f1] || dirty[f2]) continue;                    // an earlier flip moved it this pass
            if (detail::cot_opposite(V, F, f1, c1) + detail::cot_opposite(V, F, f2, c2) >= 0) continue;

            const int v1 = F(f1, (c1 + 1) % 3), v2 = F(f1, (c1 + 2) % 3);   // the edge
            const int v4 = F(f1, c1),           v3 = F(f2, c2);             // the corners it faces
            if (val[v1] <= 3 || val[v2] <= 3) continue;
            if (present.contains(detail::key(v3, v4))) continue;

            Row3d n1 = detail::tri_normal(V, F(f1, 0), F(f1, 1), F(f1, 2));
            Row3d n2 = detail::tri_normal(V, F(f2, 0), F(f2, 1), F(f2, 2));
            const double a1 = n1.norm(), a2 = n2.norm();
            if (a1 == 0 || a2 == 0) continue;

            const double qmin = std::min(detail::quality(V, F, f1), detail::quality(V, F, f2));
            const bool reliable = qmin >= sliver_quality;
            if (reliable && n1.dot(n2) / (a1 * a2) < cos_thr) continue;   // a crease: flipping would move the surface

            // the two triangles the flip produces must be non-degenerate and keep the orientation
            Row3d p1 = V.row(v1), p2 = V.row(v2), p3 = V.row(v3), p4 = V.row(v4);
            Row3d nn = (N.row(v1) + N.row(v2) + N.row(v3) + N.row(v4));
            if (nn.norm() == 0) continue;
            nn.normalize();
            Row3d m1 = (p3 - p1).cross(p4 - p1);
            Row3d m2 = (p4 - p2).cross(p3 - p2);
            const double eps = 1e-10 * (a1 + a2);
            if (m1.norm() < eps || m2.norm() < eps) continue;
            if (m1.dot(nn) <= 0 || m2.dot(nn) <= 0) continue;

            // repairing a sliver is only worth it if the pair actually comes out better
            if (!reliable && std::min(detail::quality(p1, p3, p4), detail::quality(p2, p4, p3)) <= qmin) continue;

            igl::flip_edge(F, E, uE, EMAP, uE2E, ue);
            dirty[f1] = dirty[f2] = true;
            present.erase(detail::key(v1, v2));
            present.insert(detail::key(v3, v4));
            --val[v1]; --val[v2]; ++val[v3]; ++val[v4];
            ++n;
        }

        total += n;
        if (n == 0) break;
    }
    return total;
}

// collapse every edge shorter than min_len to its midpoint, unless the collapse would make a
// neighbouring triangle thinner than it already is. V and F are rebuilt, so any index into them
// is invalidated. boundary edges are left alone. returns the number of collapses applied
inline int collapse_short_edges(MatXd& V, MatXi& F, const double min_len, const double min_angle_deg = 10., const int max_pass = 10) {
    int total = 0;

    auto min_angle = [](const Row3d t[3]) {
        double lo = 180.;
        for (int c = 0; c < 3; ++c) {
            Row3d u = t[(c + 1) % 3] - t[c], w = t[(c + 2) % 3] - t[c];
            lo = std::min(lo, std::atan2(u.cross(w).norm(), u.dot(w)) * 180. / PI);
        }
        return lo;
    };

    for (int pass = 0; pass < max_pass; ++pass) {
        MatXi uE, EF, EI;
        VecXi EMAP;
        igl::edge_flaps(F, uE, EMAP, EF, EI);

        vec<std::pair<double, int>> cand;   // (length, edge), shortest first
        for (int e = 0; e < uE.rows(); ++e) {
            const double l = (V.row(uE(e, 0)) - V.row(uE(e, 1))).norm();
            if (l < min_len) cand.emplace_back(l, e);
        }
        rg::sort(cand);

        int n = 0;
        for (const auto& [l, e]: cand) {
            if (EF(e, 0) < 0 || EF(e, 1) < 0) continue;   // boundary: moving it would change the outline
            const int a = uE(e, 0), b = uE(e, 1);
            if (a == b) continue;                          // an earlier collapse already removed it
            Eigen::RowVectorXd p = 0.5 * (V.row(a) + V.row(b));

            vec<int> Nsv, Nsf, Ndv, Ndf;                   // the faces around each end of the edge
            igl::circulation(e, true,  F, EMAP, EF, EI, Nsv, Nsf);
            igl::circulation(e, false, F, EMAP, EF, EI, Ndv, Ndf);

            // every surviving neighbour must keep its orientation and must not get thinner than
            // it already is (a mesh that is bad here should not be made worse)
            bool ok = true;
            for (const vec<int>* nf: {&Nsf, &Ndf}) {
                for (int f: *nf) {
                    if (f == EF(e, 0) || f == EF(e, 1)) continue;   // these two faces disappear
                    Row3d before[3], after[3];
                    for (int c = 0; c < 3; ++c) {
                        const int v = F(f, c);
                        before[c] = V.row(v);
                        after[c]  = (v == a || v == b) ? Row3d(p) : before[c];
                    }
                    Row3d n_old = (before[1] - before[0]).cross(before[2] - before[0]);
                    Row3d n_new = (after[1]  - after[0]).cross(after[2]  - after[0]);
                    if (n_new.norm() == 0 || n_new.dot(n_old) <= 0) { ok = false; break; }
                    if (min_angle(after) < std::min(min_angle_deg, min_angle(before))) { ok = false; break; }
                }
                if (!ok) break;
            }
            if (!ok) continue;

            int e1, e2, f1, f2;   // collapse_edge also refuses the ones that break the link condition
            if (igl::collapse_edge(e, p, Nsv, Nsf, Ndv, Ndf, V, F, uE, EMAP, EF, EI, e1, e2, f1, f2)) ++n;
        }

        total += n;
        if (n == 0) break;

        // collapse_edge blanks the faces it removes: drop them, then the vertices left orphaned
        MatXi G(F.rows(), 3);
        int g = 0;
        for (int f = 0; f < F.rows(); ++f)
            if (F(f, 0) != F(f, 1) && F(f, 1) != F(f, 2) && F(f, 2) != F(f, 0)) G.row(g++) = F.row(f);
        G.conservativeResize(g, 3);

        MatXd NV;
        MatXi NF;
        VecXi I;
        igl::remove_unreferenced(V, G, NV, NF, I);
        V = NV;
        F = NF;
    }
    return total;
}

// relax the vertices sitting on a fold, an edge whose two triangles face away from each other.
// a fold between two well shaped triangles is not a feature and not something a flip can repair:
// flipping it would invert a face. moving the vertex is the only way, so this is the one pass
// here that changes the surface. the default 150 deg means "the surface turns back on itself":
// at 90 deg this would round off the creases of a CAD model, where a right angle is a feature.
// returns the folds left
inline int smooth_folds(MatXd& V, const MatXi& F, const double fold_angle_deg = 150.,
                        const int iters = 20, const double step = 0.5) {
    const double cos_fold = std::cos(fold_angle_deg * PI / 180.);
    const int m = (int)F.rows();

    MatXi E, uE;
    VecXi EMAP;
    vec<vec<int>> uE2E;
    igl::unique_edge_map(F, E, uE, EMAP, uE2E);

    vec<vec<int>> ring(V.rows()), faces(V.rows());
    for (int e = 0; e < uE.rows(); ++e) {
        ring[uE(e, 0)].push_back(uE(e, 1));
        ring[uE(e, 1)].push_back(uE(e, 0));
    }
    for (int f = 0; f < m; ++f) for (int c = 0; c < 3; ++c) faces[F(f, c)].push_back(f);

    auto folded = [&](const int ue) {
        if (uE2E[ue].size() != 2) return false;
        const int f1 = uE2E[ue][0] % m, f2 = uE2E[ue][1] % m;
        Row3d n1 = detail::tri_normal(V, F(f1, 0), F(f1, 1), F(f1, 2));
        Row3d n2 = detail::tri_normal(V, F(f2, 0), F(f2, 1), F(f2, 2));
        const double a = n1.norm(), b = n2.norm();
        return a > 0 && b > 0 && n1.dot(n2) / (a * b) < cos_fold;
    };

    int left = 0;
    for (int it = 0; it < iters; ++it) {
        std::set<int> bad;
        left = 0;
        for (int ue = 0; ue < (int)uE2E.size(); ++ue)
            if (folded(ue)) { ++left; bad.insert(uE(ue, 0)); bad.insert(uE(ue, 1)); }
        if (bad.empty()) break;

        int moved = 0;
        for (const int v: bad) {
            if (ring[v].empty()) continue;
            Row3d c = Row3d::Zero();
            for (const int u: ring[v]) c += V.row(u);
            c /= (double)ring[v].size();

            const Row3d old = V.row(v);
            const Row3d now = old + step * (c - old);

            // keep the move only if no face of the ring turns over or gets thinner
            double q_old = 1., q_new = 1.;
            bool flip = false;
            for (const int f: faces[v]) {
                Row3d a[3], b[3];
                for (int k = 0; k < 3; ++k) {
                    a[k] = V.row(F(f, k));
                    b[k] = F(f, k) == v ? now : a[k];
                }
                Row3d na = (a[1] - a[0]).cross(a[2] - a[0]);
                Row3d nb = (b[1] - b[0]).cross(b[2] - b[0]);
                if (nb.norm() == 0 || nb.dot(na) <= 0) { flip = true; break; }
                q_old = std::min(q_old, detail::quality(a[0], a[1], a[2]));
                q_new = std::min(q_new, detail::quality(b[0], b[1], b[2]));
            }
            // no face of the ring may turn over or get thinner. relaxing this to a floor lets
            // the vertex travel through the surface and come out inverted on the far side
            if (flip || q_new < q_old) continue;
            V.row(v) = now;
            ++moved;
        }
        if (moved == 0) break;
    }
    return left;
}

// the whole cleanup: flip, collapse the shortest edges, flip again. smooth_folds is deliberately
// left out: it is the only pass that moves a vertex across a fold, and on nefertiti it puts four
// faces on the far side of the input surface while removing only 4 of the 16 worst folds.
//
// short_ratio is a fraction of the median edge length. measured across nefertiti, fandisk and
// spot, 0.2 keeps the features (fandisk 722 -> 722 sharp edges) and adds no inverted face; 0.5
// costs nefertiti 44% of its sharp edges and inverts faces on two of the three meshes, which is
// decimation rather than cleanup. pass 0 to skip the collapses entirely
inline void cleanup_mesh(MatXd& V, MatXi& F, const double feature_angle_deg = 20., const double short_ratio = 0.2) {
    delaunay_flips(V, F, feature_angle_deg);
    if (short_ratio > 0) collapse_short_edges(V, F, short_ratio * measure(V, F).len_med);
    delaunay_flips(V, F, feature_angle_deg);
}

// reduce the face count with quadric error decimation, cleaning before and after. qslim minimises
// the geometric error and ignores triangle shape, so on its own it produces slivers and faces
// lying on the wrong side of the surface: on fandisk at half the faces it leaves a 0.0 deg
// smallest angle and 57 inverted faces. the flips around it recover almost all of that (57 -> 0),
// which is why this wrapper exists instead of a bare qslim call.
//
// note that fewer triangles means fewer of them per quad, so gridscale has to grow with the
// reduction. on a model whose features live on creases (fandisk) do not go below half: the
// vertices crowd onto the creases and the triangles there cannot stay well shaped.
// returns false if qslim failed, leaving V and F cleaned but not decimated
inline bool decimate_and_clean(MatXd& V, MatXi& F, const size_t target_faces, const double feature_angle_deg = 20.) {
    cleanup_mesh(V, F, feature_angle_deg);
    if (target_faces >= (size_t)F.rows()) return true;

    MatXd U;
    MatXi G;
    VecXi J, I;
    if (!igl::qslim(V, F, target_faces, U, G, J, I)) return false;
    V = U;
    F = G;

    delaunay_flips(V, F, feature_angle_deg);
    return true;
}
}
#endif
