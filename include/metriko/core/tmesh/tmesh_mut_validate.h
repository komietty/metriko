#ifndef METRIKO_TMESH_MUT_VALIDATE_H
#define METRIKO_TMESH_MUT_VALIDATE_H
#include <set>
#include "tmesh_mut.h"

namespace metriko {

// a contact between two tedges: either a strict segment crossing inside a face
// or a "touch" (a node of one tedge lying in the interior of the other's segment)
struct TeContact {
    int  fid;
    int  te_seg;  // owner of the crossed / touched segment
    int  te_ndp;  // the other tedge (whose node or segment intrudes)
    bool cross;   // true: strict crossing, false: touch
};

// tedges must not touch each other except at shared junction nodes
inline vec<TeContact> find_tedge_contacts(const TmeshMut& tm) {
    const Hmesh& hm = tm.hm;

    auto faces_of = [&](const HmLoc& l, vec<int>& out) {
        std::visit(overloaded{
            [&](const HmLocOnV& v) { for (Half h: hm.verts[v.id].adjHalfs()) out.push_back(h.face().id); },
            [&](const HmLocOnE& e) { out.push_back(hm.edges[e.id].face0().id); out.push_back(hm.edges[e.id].face1().id); },
            [&](const HmLocOnH& h) { Half hh = hm.halfs[h.id]; out.push_back(hh.face().id); out.push_back(hh.twin().face().id); },
            [&](const HmLocOnF& f) { out.push_back(f.id); },
            [&](const auto&)       {},
        }, l);
    };

    struct Seg { int teid; int n0; int n1; };
    umap<int, vec<Seg>> per_face;
    for (const auto& [teid, nids]: tm.live_tedges()) {
        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            vec<int> fids;
            faces_of(tm.tnodes[nids[k]], fids);
            for (int fid: fids)
                if (fid != -1 && is_in_face(hm.faces[fid], tm.tnodes[nids[k + 1]]))
                    per_face[fid].push_back({teid, nids[k], nids[k + 1]});
        }
    }

    vec<TeContact> res;
    std::set<std::pair<int, int>> seen;
    auto report = [&](int fid, int te_seg, int te_ndp, bool cross) {
        if (seen.insert(std::minmax(te_seg, te_ndp)).second)
            res.push_back({fid, te_seg, te_ndp, cross});
    };

    for (auto& [fid, segs]: per_face) {
        Face f = hm.faces[fid];
        const Row3d org = f.half().tail().pos();
        const Row3d bx  = f.basisX();
        const Row3d by  = f.basisY();
        auto pos2 = [&](int nid) { Row3d d = get_ptloc_pos(hm, tm.tnodes[nid]) - org; return complex(d.dot(bx), d.dot(by)); };

        // p strictly inside segment a-b (2d)
        auto on_seg = [&](const complex p, const complex a, const complex b) {
            complex ab = b - a, ap = p - a;
            if (std::norm(ab) < EPS * EPS) return false;
            double t = (std::conj(ab) * ap).real() / std::norm(ab);
            return std::abs((std::conj(ab) * ap).imag()) / std::abs(ab) < EPS && t > EPS && t < 1 - EPS;
        };

        for (size_t i = 0; i < segs.size(); ++i)
        for (size_t j = i + 1; j < segs.size(); ++j) {
            auto& s = segs[i];
            auto& t = segs[j];
            if (s.teid == t.teid) continue;
            bool shared = s.n0 == t.n0 || s.n0 == t.n1 || s.n1 == t.n0 || s.n1 == t.n1;

            double rab, rcd;
            if (!shared && find_strict_intersection(pos2(s.n0), pos2(s.n1), pos2(t.n0), pos2(t.n1), rab, rcd)) {
                report(fid, s.teid, t.teid, true);
                continue;
            }
            if ((t.n0 != s.n0 && t.n0 != s.n1 && on_seg(pos2(t.n0), pos2(s.n0), pos2(s.n1))) ||
                (t.n1 != s.n0 && t.n1 != s.n1 && on_seg(pos2(t.n1), pos2(s.n0), pos2(s.n1)))) { report(fid, s.teid, t.teid, false); continue; }
            if ((s.n0 != t.n0 && s.n0 != t.n1 && on_seg(pos2(s.n0), pos2(t.n0), pos2(t.n1))) ||
                (s.n1 != t.n0 && s.n1 != t.n1 && on_seg(pos2(s.n1), pos2(t.n0), pos2(t.n1)))) { report(fid, t.teid, s.teid, false); }
        }
    }
    return res;
}

inline int validate_no_crossing(const TmeshMut& tm, const char* stage) {
    auto cs = find_tedge_contacts(tm);
    for (auto& [fid, te_seg, te_ndp, cross]: cs)
        std::println("[crossing] {}: {} at face {}, teid {} x teid {}", stage, cross ? "cross" : "touch", fid, te_seg, te_ndp);
    std::println("[crossing] {}: {} total", stage, cs.size());
    return (int)cs.size();
}

// resolve contacts by re-tracing the intruding tedge (fall back to the other)
inline void repair_crossing_tedges(TmeshMut& tm, const int max_iter = 10) {
    for (int it = 0; it < max_iter; ++it) {
        auto cs = find_tedge_contacts(tm);
        if (cs.empty()) return;
        for (auto& c: cs) {
            if (tm.tedges[c.te_ndp].id != -1 && tm.reroute_tedge(c.te_ndp)) { std::println("[crossing] rerouted teid {} (contact with teid {})", c.te_ndp, c.te_seg); continue; }
            if (tm.tedges[c.te_seg].id != -1 && tm.reroute_tedge(c.te_seg)) { std::println("[crossing] rerouted teid {} (contact with teid {})", c.te_seg, c.te_ndp); }
        }
    }
}
}
#endif
