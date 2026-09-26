#ifndef METRIKO_EMESH_VALIDATE_H
#define METRIKO_EMESH_VALIDATE_H
#include <set>
#include "emesh.h"

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
inline vec<TeContact> find_tedge_contacts(const Emesh& tm) {
    struct S { int teid; int n0; int n1; };
    umap<int, vec<S>> per_face;

    for (const auto& [teid, nids]: tm.live_tedges()) {
    for (int k = 0; k + 1 < nids.size(); ++k) {
    for (int i: get_ptloc_faces(tm.hm, tm.tnodes[nids[k]])) {
        if (is_in_face(tm.hm.faces[i], tm.tnodes[nids[k + 1]]))
            per_face[i].push_back({.teid = teid, .n0 = nids[k], .n1 = nids[k + 1]});
    }}}

    vec<TeContact> res;
    std::set<std::pair<int, int>> seen;
    auto report = [&](int fid, int te_seg, int te_ndp, bool cross) {
        if (seen.insert(std::minmax(te_seg, te_ndp)).second)
            res.push_back({.fid=fid, .te_seg=te_seg, .te_ndp=te_ndp, .cross=cross});
    };

    for (auto& [fid, segs]: per_face) {
        auto f = tm.hm.faces[fid];
        for (size_t i = 0; i < segs.size(); ++i)
        for (size_t j = i + 1; j < segs.size(); ++j) {
            auto& s = segs[i];
            auto& t = segs[j];
            if (s.teid == t.teid) continue;

            auto a0 = f.to_local(get_ptloc_pos(tm.hm, tm.tnodes[s.n0])),
                 a1 = f.to_local(get_ptloc_pos(tm.hm, tm.tnodes[s.n1])),
                 b0 = f.to_local(get_ptloc_pos(tm.hm, tm.tnodes[t.n0])),
                 b1 = f.to_local(get_ptloc_pos(tm.hm, tm.tnodes[t.n1]));

            bool shared = s.n0 == t.n0 || s.n0 == t.n1 || s.n1 == t.n0 || s.n1 == t.n1;
            if (!shared && find_strict_intersection(a0, a1, b0, b1)) { report(fid, s.teid, t.teid, true); continue; }
            if (is_inside_segment(a0, a1, b0) || is_inside_segment(a0, a1, b1)) { report(fid, s.teid, t.teid, false); continue; } // t's endpoint on s
            if (is_inside_segment(b0, b1, a0) || is_inside_segment(b0, b1, a1)) { report(fid, t.teid, s.teid, false); }           // s's endpoint on t
        }
    }
    return res;
}

inline int validate_no_crossing(const Emesh& tm, const char* stage) {
    auto cs = find_tedge_contacts(tm);
    for (auto& [fid, te_seg, te_ndp, cross]: cs)
        std::println("[crossing] {}: {} at face {}, teid {} x teid {}", stage, cross ? "cross" : "touch", fid, te_seg, te_ndp);
    if (!cs.empty()) std::println("[crossing] {}: {} total", stage, cs.size());
    return (int)cs.size();
}

// resolve contacts by re-tracing the intruding tedge (fall back to the other)
inline void repair_crossing_tedges(Emesh& tm, const int max_iter = 10) {
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
