#ifndef METRIKO_TUTTE_CUTTING_H
#define METRIKO_TUTTE_CUTTING_H
#include <ranges>
#include <set>
#include <unordered_map>
#include "./tutte.h"
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/tmesh/tmesh_mut.h"

namespace metriko {

// split points on an original halfedge: (ratio, vid)
using EdgeSplits = std::set<std::pair<double, int>>;

struct SegRec {
    int i0;
    int i1;
    HalfData d0;
    HalfData d1;
};

struct IntermidiateData {

};

inline bool is_on_side(const vec<vec<int>>& sides, int vid, int sid) {
    return vid < sides.size() && rg::contains(sides[vid], sid);
}

inline int common_side(const vec<vec<int>>& sides, int vid0, int vid1) {
    if (vid0 >= sides.size()) return -1;
    for (int l: sides[vid0])
        if (is_on_side(sides, vid1, l)) return l;
    return -1;
}

inline int common_side(const vec<vec<int>>& sides, int vid0, int vid1, int vid2) {
    if (vid0 >= sides.size()) return -1;
    for (int l: sides[vid0])
        if (is_on_side(sides, vid1, l) &&
            is_on_side(sides, vid2, l)) return l;
    return -1;
}

// vertex chain along an original halfedge, tail -> head, split points included
inline vec<int> chain_of(const vec<EdgeSplits>& splits, Half h) {
    vec chain = { h.tail().id };
    for (const auto& [r, vid]: vw::reverse(splits[h.id])) chain.push_back(vid);  // descending r = tail -> head
    chain.push_back(h.head().id);
    return chain;
}

// triangulate one face against the segments crossing it. the resulting pieces may
// be CONCAVE (traced tedges bend inside a face), so the trace allows reflex turns
// and the pieces are ear-clipped rather than fanned.
inline void face_cutting(
    const Face& f,                 // face to be cut
    const vec<Row2i>& sgs,         // segment endpoints (cut-mesh vertex indices)
    const vec<EdgeSplits>& splits, // split points per original halfedge
    const vec<vec<int>>& sides,    // vertex -> tquad * side
    const vec<Row3d>& vpos,        // vertex position
    vec<Row3i>& tris               // output triangles (appended, CCW wrt f)
) {
    /// 1: collect the directed edges bounding the pieces:
    ///    segments (both directions) + subdivided boundary halfedges (one direction)
    vec<Row2i> halfs;
    for (auto s: sgs) {
        halfs.emplace_back(s.x(), s.y());
        halfs.emplace_back(s.y(), s.x());
    }
    for (Half h: f.adjHalfs()) {
        const auto& sp = splits[h.id];
        if (sp.empty()) { halfs.emplace_back(h.tail().id, h.head().id); continue; }
        halfs.emplace_back(sp.begin()->second, h.head().id);
        auto prev = sp.begin();
        auto curr = std::next(sp.begin());
        for (; curr != sp.end(); ++curr, ++prev) { halfs.emplace_back(curr->second, prev->second); }
        halfs.emplace_back(h.tail().id, sp.rbegin()->second);
    }

    const Row3d org = f.half().tail().pos();
    const Row3d bx  = f.basisX();
    const Row3d by  = f.basisY();
    auto pos2 = [&](int vid) { Row3d d = vpos[vid] - org; return complex(d.dot(bx), d.dot(by)); };
    auto cr   = [](complex u, complex v) { return (std::conj(u) * v).imag(); };

    constexpr double penalty = 100.; // dominates min_angle in (-pi, pi]

    /// 2: trace the pieces by the leftmost-turn rule (max signed CCW angle).
    ///    reflex turns are legal, and every chain must close.
    while (!halfs.empty()) {
        auto h = halfs.back();
        halfs.pop_back();

        vec poly = {h};
        int sta = h.x();
        int cur = h.y();

        while (cur != sta) {
            int prev_v = poly.back().x();
            int curr_v = poly.back().y();
            const Row3d& p_prev = vpos[prev_v];
            const Row3d& p_curr = vpos[curr_v];
            const Row3d  d_prev = (p_curr - p_prev).normalized();

            auto   best_it = halfs.end();
            double best    = -10;   // signed turn angle in (-pi, pi]
            for (auto it = halfs.begin(); it != halfs.end(); ++it) {
                if (it->x() != curr_v) continue;
                if (it->y() == prev_v) continue;   // no immediate u-turn
                Row3d d_cand = (vpos[it->y()] - p_curr).normalized();
                double ang = std::atan2(f.normal().dot(d_prev.cross(d_cand)), d_prev.dot(d_cand));
                if (ang > best) { best = ang; best_it = it; }
            }

            if (best_it == halfs.end()) throw std::runtime_error(std::format("[cut]: open chain (face {})", f.id));

            poly.push_back(*best_it);
            cur = best_it->y();
            halfs.erase(best_it);
        }

        /// 3: ear-clip the piece in the face plane (concave-capable, keeps CCW
        ///    orientation so no flipped faces can be produced)
        auto cyc = poly | vw::transform([](const Row2i& p) { return p.x(); }) | rg::to<vec<int>>();

        // interior angles of a CCW triangle in the face plane (all positive);
        // negative for a CW triangle, so it also ranks invalid choices last
        auto min_angle = [&](int va, int vb, int vc) {
            complex a = pos2(va), b = pos2(vb), c = pos2(vc);
            return std::min({std::arg((c - a) / (b - a)),
                             std::arg((a - b) / (c - b)),
                             std::arg((b - c) / (a - c))});
        };

        assert(cyc.size() >= 3);

        while (cyc.size() > 3) {
            const size_t n = cyc.size();
            size_t best_k = n;
            double best_q = std::numeric_limits<double>::lowest();
            for (size_t k = 0; k < n; ++k) {
                int va = cyc[(k + n - 1) % n];
                int vb = cyc[k];
                int vc = cyc[(k + 1) % n];
                complex a = pos2(va), b = pos2(vb), c = pos2(vc);
                if (cr(b - a, c - b) <= EPS) continue;   // reflex or flat corner: not an ear

                bool empty = true;                        // no other cycle vertex inside the ear
                for (size_t m = 0; m < n && empty; ++m) {
                    if (m == k || m == (k + 1) % n || m == (k + n - 1) % n) continue;
                    complex q = pos2(cyc[m]);
                    empty = cr(b - a, q - a) < -EPS || cr(c - b, q - b) < -EPS || cr(a - c, q - c) < -EPS;
                }
                if (!empty) continue;

                // quality of this clip: the worst interior angle it commits to.
                // a quad also fixes the remaining triangle, so score the pair —
                // this picks the better diagonal and avoids slivers whenever the
                // other choice is available
                double q = min_angle(va, vb, vc);
                if (n == 4) q = std::min(q, min_angle(vc, cyc[(k + 2) % n], va));
                if (common_side(sides, va, vb, vc) >= 0)                         q -= penalty;
                if (n == 4 && common_side(sides, vc, cyc[(k + 2) % n], va) >= 0) q -= penalty;
                if (q > best_q) { best_q = q; best_k = k; }
            }
            if (best_k == n) throw std::runtime_error(std::format("[cut]: ear clipping failed (face {}, size {})", f.id, cyc.size()));
            if (best_q < -4) std::println("[cut] uv-degenerate triangle unavoidable (face {})", f.id);
            tris.emplace_back(cyc[(best_k + n - 1) % n], cyc[best_k], cyc[(best_k + 1) % n]);
            cyc.erase(cyc.begin() + (long)best_k);
        }

        if (common_side(sides, cyc[0], cyc[1], cyc[2]) >= 0) std::println("[cut] uv-degenerate triangle unavoidable (face {})", f.id);

        tris.emplace_back(cyc[0], cyc[1], cyc[2]);
    }
}

inline void split_face(
    const Hmesh& hm,
    const vec<bool>& seam0,
    const vec<SegRec>& sgs,
    const vec<EdgeSplits>& h_aux,
    const vec<vec<int>>& sides,
    vec<Row3d>& vpos,
    vec<Row3i>& tris
) {
    std::set<std::pair<int, int>> sg_set;
    std::set<std::pair<int, int>> sm_sub;

    for (auto& s: sgs) sg_set.insert(std::minmax(s.i0, s.i1));

    for (Edge e: hm.edges) {
        if (!seam0[e.id]) continue;
        auto ch = chain_of(h_aux, e.half());
        for (size_t k = 0; k + 1 < ch.size(); ++k)
            sm_sub.insert(std::minmax(ch[k], ch[k + 1]));
    }

    for (bool again = true; std::exchange(again, false);) {
        std::map<std::pair<int, int>, int> tri_of;

        for (int i = 0; i < tris.size(); ++i)
        for (int j = 0; j < 3; ++j)
            tri_of[{tris[i][j], tris[i][(j + 1) % 3]}] = i;

        for (int i = 0; i < tris.size() && !again; ++i) {
            int l = common_side(sides, tris[i].x(), tris[i].y(), tris[i].z());
            if (l < 0) continue;
            for (int j = 0; j < 3; ++j) {
                int a = tris[i][j];
                int b = tris[i][(j + 1) % 3];
                int c = tris[i][(j + 2) % 3];
                if (sg_set.contains(std::minmax(a, b))) continue;
                if (sm_sub.contains(std::minmax(a, b))) { std::println("[cut] [warn]: cut on seam happens"); continue; } // todo: might cause edge case
                int k = tri_of.at({b, a});
                int d = tris[k][0] + tris[k][1] + tris[k][2] - a - b;
                int m = vpos.size();
                if (is_on_side(sides, d, l)) continue;
                Row3d mid = (vpos[a] + vpos[b]) / 2;
                tris[i] = {a, m, c};
                tris[k] = {b, m, d};
                tris.emplace_back(m, b, c);
                tris.emplace_back(m, a, d);
                vpos.emplace_back(mid);
                again = true;
                break;
            }
        }

        for (int i = 0; i < tris.size() && !again; ++i) {
        for (int j = 0; j < 3; ++j) {
            int a = tris[i][j];
            int b = tris[i][(j + 1) % 3];
            int c = tris[i][(j + 2) % 3];
            if (int s = common_side(sides, a, b); s >= 0 && !is_on_side(sides, c, s) && !sg_set.contains(std::minmax(a, b))) {
                if (sm_sub.contains(std::minmax(a, b))) { std::println("[cut] [warn]: cut on seam happens"); continue; } // todo: might cause edge case
                int k = tri_of.at({b, a});
                int d = tris[k][0] + tris[k][1] + tris[k][2] - a - b;
                int m = vpos.size();
                if (is_on_side(sides, d, s)) continue;
                Row3d mid = (vpos[a] + vpos[b]) / 2;
                tris[i] = {a, m, c};
                tris[k] = {b, m, d};
                tris.emplace_back(m, b, c);
                tris.emplace_back(m, a, d);
                vpos.emplace_back(mid);
                again = true;
                break;
            }
        }
        }
    }
}

inline std::unique_ptr<Hmesh> compute_embedding_cut_hmesh(
    const Hmesh& hm,
    const TmeshMut& tm,
    const vec<bool>& seam0,
    const VecXi& matching0,
    const VecXi& singular0,
          vec<bool>& seam1,
          VecXi& matching1,
          VecXi& singular1,
    vec<HalfData>& data
) {
    std::map<int, vec<Row2i>> cuts; // face id -> segment endpoints

    vec<SegRec> sgms; // flat annotation list; matched to cut halfedges at the end

    vec auxs(hm.nH, EdgeSplits{});

    vec<Row3d> vpos;
    for (auto p: hm.pos.rowwise()) { vpos.emplace_back(p); }

    // nid -> cut-mesh vertex index: tnode identity replaces epsilon-based position
    // matching. creates the vertex and registers the edge split exactly once.
    umap<int, int> vid_of_nid;
    vec<vec<int>> sides; // cut vertex -> tquad * side

    auto vid_of = [&](int nid) -> int {
        if (!vid_of_nid.contains(nid)) {
            vid_of_nid[nid] = std::visit(overloaded{
                [&](const HmLocOnV& v) { return v.id; },
                [&](const auto& l) {
                    vpos.emplace_back(get_ptloc_pos(hm, l));
                    int i = vpos.size() - 1;
                    if (auto hr = try_get_ratio(hm, l)) {
                        auto [h0, r] = hr.value();
                        auxs[h0.id].emplace(r, i);
                        auxs[h0.twin().id].emplace(1 - r, i);
                    }
                    return i;
                },
            }, tm.tnodes[nid]);
        }
        return vid_of_nid[nid];
    };

    for (const auto& th: tm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids = tm.tedges[th.teid].nids;
        const int   nsgs = nids.size() - 1;

        const auto& tw = tm.thalfs[th.twid];
        const int   sa = th.tqid * 4 + tm.tquads[th.tqid].side_of(th);
        const int   sb = tw.tqid * 4 + tm.tquads[tw.tqid].side_of(tw);

        double total = 0;
        vec<double> len(nsgs);
        for (int k = 0; k < nsgs; ++k) {
            len[k] = (get_ptloc_pos(hm, tm.tnodes[nids[k + 1]]) - get_ptloc_pos(hm, tm.tnodes[nids[k]])).norm();
            total += len[k];
        }
        if (total < EPS) throw std::runtime_error("[cut]: zero-length tedge");

        double sum = 0;
        for (int k = 0; k < nsgs; ++k) {
            auto& fr = tm.tnodes[nids[k]];
            auto& to = tm.tnodes[nids[k + 1]];
            auto  v0 = sum / total; sum += len[k];
            auto  v1 = sum / total;
            auto  i0 = vid_of(nids[k]);
            auto  i1 = vid_of(nids[k + 1]);
            for (int i: {i0, i1}) {
                if (i >= sides.size()) sides.resize(vpos.size());
                auto& ss = sides[i];
                if (!rg::contains(ss, sa)) ss.push_back(sa);
                if (!rg::contains(ss, sb)) ss.push_back(sb);
            }
            auto eo = try_get_edge(hm, fr, to);
            auto fo = try_get_face(hm, fr, to);
            auto t  = nsgs - k - 1;
            if (!eo.has_value()) cuts[fo.value().id].emplace_back(i0, i1);
            sgms.push_back({
                .i0 = i0, .i1 = i1,
                .d0 = HalfData{.v0 = v0,     .v1 = v1,     .thid = th.id,   .tqid = th.tqid,                 .order = k},
                .d1 = HalfData{.v0 = 1 - v1, .v1 = 1 - v0, .thid = th.twid, .tqid = tm.thalfs[th.twid].tqid, .order = t}
            });
        }
    }

    vec<Row3i> tris;
    for (Face f: hm.faces) {
        auto h0 = f.half();
        auto h1 = f.half().next();
        auto h2 = f.half().prev();
        if (cuts[f.id].empty()  &&
            auxs[h0.id].empty() &&
            auxs[h1.id].empty() &&
            auxs[h2.id].empty()
        ) { tris.emplace_back(h0.tail().id, h1.tail().id, h2.tail().id); }
        else { face_cutting(f, cuts[f.id], auxs, sides, vpos, tris); }
    }

    split_face(hm, seam0, sgms, auxs, sides, vpos, tris);

    MatXi face_info(tris.size(), 3);
    MatXd vert_info(vpos.size(), 3);
    for (int i = 0; i < vpos.size(); i++) { vert_info.row(i) = vpos[i]; }
    for (int i = 0; i < tris.size(); i++) { face_info.row(i) << tris[i][0], tris[i][1], tris[i][2]; }

    auto hm_cut = std::make_unique<Hmesh>(vert_info, face_info);
    seam1     = vec(hm_cut->nE, false);
    matching1 = VecXi::Zero(hm_cut->nE);
    singular1 = VecXi::Zero(hm_cut->nV);
    singular1.head(hm.nV) = singular0;

    std::map<std::pair<int, int>, Half> half_by_verts;
    for (Half h: hm_cut->halfs)
        half_by_verts.insert({{h.tail().id, h.head().id}, h});

    for (Edge e: hm.edges) {
        if (!seam0[e.id]) continue;
        auto ch = chain_of(auxs, e.half());
        for (int k = 0; k + 1 < ch.size(); ++k) {
            Half h = half_by_verts.at({ch[k], ch[k + 1]});   // same direction as h
            seam1[h.edge().id] = true;
            matching1(h.edge().id) = h.isCanonical() ? matching0(e.id) : -matching0(e.id);
        }
    }

    // annotate the cut halfedges: both directions pushed adjacently, twin = index ^ 1
    data.clear();
    for (auto& [i0, i1, d0, d1]: sgms) {
        auto it = half_by_verts.find({i0, i1});
        if (it == half_by_verts.end()) throw std::runtime_error("No half_by_verts");

        d0.half = it->second;
        d1.half = it->second.twin();
        d0.twin = data.size() + 1;
        d1.twin = data.size();
        data.push_back(d0);
        data.push_back(d1);
    }

    std::sort(data.begin(), data.end());
    return hm_cut;
}
}
#endif
