#ifndef METRIKO_QEX_GEN_Q_EDGE_H
#define METRIKO_QEX_GEN_Q_EDGE_H
#include "common.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::qex {
    // sign of orientation with the global tolerance
    inline int orient_sign(const complex a, const complex b, const complex c) {
        double o = orientation(a, b, c);
        return o > EPS ? 1 : o < -EPS ? -1 : 0;
    }

    // closed point-in-triangle test (points on an edge or a vertex count as inside)
    inline bool is_inside_face_closed(const Face f, const VecXc &cf, const complex p) {
        auto a = cf(f.id * 3), b = cf(f.id * 3 + 1), c = cf(f.id * 3 + 2);
        auto s = orient_sign(a, b, c);
        if (s == 0) return false;
        return orient_sign(a, b, p) * s >= 0 &&
               orient_sign(b, c, p) * s >= 0 &&
               orient_sign(c, a, p) * s >= 0;
    }

    // does the closed segment [p,q] meet the half-open segment (a,b]?
    inline bool meets_half_open(const complex a, const complex b, const complex p, const complex q) {
        const int op = orient_sign(a, b, p), oq = orient_sign(a, b, q);
        if (op == 0 && oq == 0) {// collinear: overlap along the line
            double L2 = std::norm(b - a);
            auto t = [&](const complex x) { return ((x - a) * std::conj(b - a)).real() / L2; };
            double t0 = std::min(t(p), t(q)), t1 = std::max(t(p), t(q));
            return t1 > EPS && t0 <= 1 + EPS;
        }
        if (op * oq > 0) return false; // p and q strictly on one side of the ray line
        int oa = orient_sign(p, q, a);
        int ob = orient_sign(p, q, b);
        if (oa * ob > 0)        return false; // a and b strictly on one side of the edge line
        if (oa == 0 && ob != 0) return false; // the only contact is the excluded point a
        return true;
    }

    // p's outgoing direction is the reverse of the ray (a, d) traced in face f, compared in 3d so that
    // ports living in a neighbouring chart can be matched without accumulating transitions
    inline bool is_reverse_port(const Hmesh &hm, const VecXc &cf, const complex a, const complex d, const Face f, const Qport &p) {
        Row3d pb1 = conversion_2d_3d(f, cf, a);
        Row3d pb2 = conversion_2d_3d(f, cf, a + d);
        Row3d pa1 = conversion_2d_3d(hm.faces[p.fid], cf, p.uv);
        Row3d pa2 = conversion_2d_3d(hm.faces[p.fid], cf, p.uv + p.dir);
        Row3d da = (pa2 - pa1).normalized();
        Row3d db = (pb2 - pb1).normalized();
        return std::abs(1. + da.dot(db)) < 1e-7;
    }

    // PICK_NEXT_EDGE of Ebke et al. 2013 (Alg. 5): the edge of f other than the one entered through that
    // meets (a,b]. if two edges meet it (the ray passes through a vertex or runs along an edge), take the one
    // with fewer endpoints on the ray line, which steps around the vertex until the ray leaves properly
    inline std::optional<Half> pick_next_half(
        const VecXc &cf,
        const complex a,
        const complex b,
        const Face f,
        const int skip_hid
    ) {
        std::optional<Half> pick;
        int best = 3;
        for (Half h: f.adjHalfs()) {
            if (h.id == skip_hid) continue;
            const complex p = cf(h.next().crnr().id);
            const complex q = cf(h.prev().crnr().id);
            if (!meets_half_open(a, b, p, q)) continue;
            const int on_line = (orient_sign(a, b, p) == 0) + (orient_sign(a, b, q) == 0);
            if (on_line < best) { best = on_line; pick = h; }
        }
        return pick;
    }

    inline std::vector<Qedge> generate_q_edge(
        const Hmesh &hm,
        const VecXc &cfn,
        const VecXi &matching,
        std::vector<Qport> &qports
    ) {
        std::vector<Qedge> qedges;
        VecXc heR;
        VecXc heT;
        compute_trs_matrix(hm, cfn, matching, 4, heR, heT);

        // ports grouped by carrier for the arrival lookup
        std::map<int, vec<Qport*>> byV, byE, byF;
        for (Qport &p: qports) {
            if      (p.vid >= 0) byV[p.vid].push_back(&p);
            else if (p.eid >= 0) byE[p.eid].push_back(&p);
            else                 byF[p.fid].push_back(&p);
        }

        for (Qport &pfr: qports) {
            if (pfr.isConnected) continue;
            complex a = pfr.uv;
            complex d = pfr.dir;
            complex b = a + d;
            int fid  = pfr.fid;
            int e_in = -1;
            if (pfr.eid >= 0) { // an eqvert port starts on its own edge: never cross back over it
                Half h0 = hm.edges[pfr.eid].half();
                e_in = (h0.face().id == fid ? h0 : h0.twin()).id;
            }

            for (int step = 0; step < 1000; ++step) {
                Face f = hm.faces[fid];

                if (is_inside_face_closed(f, cfn, b)) {
                    // the target grid point lies in the closure of f: find its q-vertex (vertex, edge, then face)
                    Qport* hit = nullptr;
                    for (Half h: f.adjHalfs()) {
                        if (std::abs(cfn(h.next().crnr().id) - b) >= EPS || !byV.contains(h.tail().id)) continue;
                        for (Qport* p: byV[h.tail().id])
                            if (!p->isConnected && p != &pfr && is_reverse_port(hm, cfn, a, d, f, *p)) { hit = p; break; }
                        if (hit) break;
                    }
                    if (!hit) {
                        Row3d pb = conversion_2d_3d(f, cfn, b);
                        for (Half h: f.adjHalfs()) {
                            if (orient_sign(cfn(h.next().crnr().id), cfn(h.prev().crnr().id), b) != 0 || !byE.contains(h.edge().id)) continue;
                            for (Qport* p: byE[h.edge().id])
                                if (!p->isConnected && p != &pfr && (p->pos - pb).norm() < 1e-6 && is_reverse_port(hm, cfn, a, d, f, *p)) { hit = p; break; }
                            if (hit) break;
                        }
                    }
                    if (!hit && byF.contains(fid)) {
                        for (Qport* p: byF[fid])
                            if (!p->isConnected && p != &pfr && std::abs(p->uv - b) < EPS && std::abs(p->dir + d) < EPS) { hit = p; break; }
                    }
                    if (hit) {
                        pfr.isConnected = true;
                        hit->isConnected = true;
                        qedges.emplace_back(pfr, *hit);
                    }
                    break; // reached b: paired, or left dangling
                }

                auto nh = pick_next_half(cfn, a, b, f, e_in);
                if (!nh) break; // numerically inconsistent chart: leave dangling
                if (nh->twin().isBoundary()) throw std::runtime_error("not implemented yet");
                auto r = heR(nh->id);
                auto t = heT(nh->id);
                a = r * a + t;
                b = r * b + t;
                d = r * d;
                fid  = nh->twin().face().id;
                e_in = nh->twin().id;
            }
        }

        for (const Qport& p: qports)
            if (!p.isConnected) std::println("[qedge] unpaired port {} (vid {}, eid {}, fid {})", p.idx, p.vid, p.eid, p.fid);

        return qedges;
    }
}

#endif
