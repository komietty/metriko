#ifndef METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#define METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#include <set>
#include <queue>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "metriko/core/tmesh/emesh.h"
#include "metriko/core/qex/common.h"

namespace metriko::visualizer {

// replay the t-mesh on the quad graph and colour the quads by tquad.
//
// singular vertices are q-vertices, so their quad vertex is exact. around a node the tedges are chained CCW
// combinatorially: consecutive tedges share a tquad (left of one = right of the next), and that tquad spans one
// quad-edge slot when the node is its corner, two when the node lies inside one of its sides. at a singular only
// the rotation of this chain against the quad edges is unknown; it is chosen by a joint vote of all its tedges
// (3d direction of the first tedge segment vs the quad edge). every tedge then walks exactly x quad edges, which
// fixes its far node, and the tedges there are chained the same way. no rim pinning and no nearest-vertex search.
inline void visualize_quad_patch(
    const Hmesh& hm,
    const Emesh& tm,
    const VecXi& singular,
    const vec<qex::Qface>& qfaces,
    const MatXd& qv,    // displayed quad mesh vertices (from visualize_qfaces)
    const MatXi& qidx   // displayed quad mesh faces, corner (i, j) <-> qfaces[i].qhalfs[j].port1()
) {
    const int l = (int)qfaces.size();

    // directed edge -> (quad, corner); vertex -> outgoing heads; embedding vertex -> quad vertex
    std::map<std::pair<int, int>, std::pair<int, int>> de;
    umap<int, vec<int>> out;
    umap<int, int> qv_of_vid;
    for (int i = 0; i < l; ++i)
    for (int j = 0; j < 4; ++j) {
        int a = qidx(i, j), b = qidx(i, (j + 1) % 4);
        de[{a, b}] = {i, j};
        out[a].push_back(b);
        if (int v = qfaces[i].qhalfs[j].port1().vid; v >= 0) qv_of_vid[v] = a;
    }
    auto rot      = [&](int v, int w) { auto [q, j] = de.at({v, w}); return qidx(q, (j + 3) % 4); };   // next head CCW around v
    auto straight = [&](int t, int v) { return rot(v, rot(v, t)); };                                 // arrived t -> v: continue

    umap<int, int> xs, te2th;                 // teid -> steps, cano thalf
    umap<int, std::pair<int, int>> sides;     // teid -> (left, right) tqid along nids
    for (auto& th: tm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        xs[th.teid]    = (int)std::round(th.x);
        te2th[th.teid] = th.id;
        sides[th.teid] = {th.tqid, tm.thalfs[th.twid].tqid};
    }
    umap<int, vec<std::pair<int, bool>>> node_te;   // node -> (teid, leaves the node forward)
    for (auto& [teid, nids]: tm.live_tedges()) {
        node_te[nids.front()].emplace_back(teid, true);
        node_te[nids.back()].emplace_back(teid, false);
    }
    auto dir_at = [&](int teid, int nid) -> Row3d {   // 3d direction of the first segment leaving nid
        auto& nids = tm.tedges[teid].nids;
        int nb = nids.front() == nid ? nids[1] : nids[nids.size() - 2];
        return (get_ptloc_pos(hm, tm.tnodes[nb]) - get_ptloc_pos(hm, tm.tnodes[nid])).normalized();
    };
    auto lr_out  = [&](int te, bool fwd) { auto [lq, rq] = sides.at(te); return fwd ? std::pair(lq, rq) : std::pair(rq, lq); };
    auto side_in = [&](int te, int P) {   // side of tquad P bounded by tedge te
        int c = te2th.at(te);
        int h = tm.thalfs[c].tqid == P ? c : tm.thalfs[c].twid;
        return tm.tquads[P].side_of(tm.thalfs[h]);
    };
    // CCW chain of the tedges at a node starting from (te, fwd): (teid, fwd, slot offset)
    auto chain_at = [&](int nid, int te, bool fwd) {
        vec<std::tuple<int, bool, int>> ch = {{te, fwd, 0}};
        const auto& ports = node_te[nid];
        int off = 0;
        for (size_t it = 1; it < ports.size(); ++it) {
            int P = lr_out(te, fwd).first;   // the tquad CCW after te
            auto nx = rg::find_if(ports, [&](auto& pr) { return pr.first != te && lr_out(pr.first, pr.second).second == P; });
            if (nx == ports.end()) break;
            off += side_in(te, P) == side_in(nx->first, P) ? 2 : 1;
            ch.emplace_back(nx->first, nx->second, off);
            te  = nx->first;
            fwd = nx->second;
        }
        return ch;
    };

    vec<double> tqid_of_quad(l, -1);
    std::set<std::pair<int, int>> track;
    vec<bool> done(tm.tedges.size(), false);
    umap<int, int> node_qv;   // node -> quad vertex
    struct Job { int teid; bool fwd; int s; int b; };
    std::queue<Job> jobs;

    // seeds: singular nodes
    int nseed = 0, bad_chain = 0;
    for (auto& [nid, ports]: node_te) {
        auto* lv = std::get_if<HmLocOnV>(&tm.tnodes[nid]);
        if (!lv || !singular(lv->id) || !qv_of_vid.contains(lv->id)) continue;
        int s = qv_of_vid.at(lv->id);
        node_qv[nid] = s;
        ++nseed;
        vec<int> ring = {out[s].front()};   // quad edges around s, CCW
        while ((int)ring.size() < (int)out[s].size()) ring.push_back(rot(s, ring.back()));
        auto ch = chain_at(nid, ports[0].first, ports[0].second);
        if (ch.size() != ports.size()) ++bad_chain;
        const int n = (int)ring.size();
        int best_r = 0;
        double best = -std::numeric_limits<double>::infinity();
        for (int r = 0; r < n; ++r) {
            double sc = 0;
            for (auto [te, fwd, off]: ch) sc += (Row3d(qv.row(ring[(off + r) % n])) - Row3d(qv.row(s))).normalized().dot(dir_at(te, nid));
            if (sc > best) { best = sc; best_r = r; }
        }
        for (auto [te, fwd, off]: ch) jobs.push({te, fwd, s, ring[(off + best_r) % n]});
    }

    // walks
    int failed = 0, sing_wrong = 0;
    while (!jobs.empty()) {
        Job jb = jobs.front(); jobs.pop();
        if (done[jb.teid]) continue;
        done[jb.teid] = true;
        auto [lq, rq] = sides.at(jb.teid);
        if (!jb.fwd) std::swap(lq, rq);
        int a = jb.s, b = jb.b;
        bool ok = true;
        for (int k = 0, steps = xs.at(jb.teid); k < steps; ++k) {
            tqid_of_quad[de.at({a, b}).first] = lq;
            tqid_of_quad[de.at({b, a}).first] = rq;
            track.insert(std::minmax(a, b));
            if (k + 1 == steps) break;
            if ((int)out[b].size() != 4) { ok = false; break; }   // hit an irregular vertex early
            int c = straight(a, b);
            a = b;
            b = c;
        }
        if (!ok) { ++failed; continue; }
        const auto& nids = tm.tedges[jb.teid].nids;
        int nid = jb.fwd ? nids.back() : nids.front();
        if (auto* lv = std::get_if<HmLocOnV>(&tm.tnodes[nid]); lv && singular(lv->id)) {   // singular far end: must agree
            if (!qv_of_vid.contains(lv->id) || qv_of_vid.at(lv->id) != b) ++sing_wrong;
            continue;
        }
        if (node_qv.contains(nid)) continue;
        node_qv[nid] = b;
        if ((int)out[b].size() != 4) { ++failed; continue; }
        int ring[4] = {a};   // slots CCW from the reverse of arrival
        for (int k = 1; k < 4; ++k) ring[k] = rot(b, ring[k - 1]);
        for (auto [te, fwd, off]: chain_at(nid, jb.teid, !jb.fwd)) if (!done[te]) jobs.push({te, fwd, b, ring[off % 4]});
    }
    int unreached = 0;
    for (auto& [teid, nids]: tm.live_tedges()) if (xs[teid] > 0 && !done[teid]) ++unreached;

    // flood the patch interiors without crossing a track
    std::queue<int> que;
    for (int i = 0; i < l; ++i) if (tqid_of_quad[i] >= 0) que.push(i);
    while (!que.empty()) {
        int q = que.front(); que.pop();
        for (int j = 0; j < 4; ++j) {
            int a = qidx(q, j), b = qidx(q, (j + 1) % 4);
            if (track.contains(std::minmax(a, b))) continue;
            int nb = de.at({b, a}).first;
            if (tqid_of_quad[nb] < 0) { tqid_of_quad[nb] = tqid_of_quad[q]; que.push(nb); }
        }
    }
    const int unlabeled = (int)rg::count(tqid_of_quad, -1.);

    // show: nodes (singular / junction), tracks, patches
    std::vector<glm::vec3> ps, pe;
    for (auto& [nid, v]: node_qv) {
        auto* lv = std::get_if<HmLocOnV>(&tm.tnodes[nid]);
        (lv && singular(lv->id) ? ps : pe).emplace_back(qv(v, 0), qv(v, 1), qv(v, 2));
    }
    auto* pc_s = polyscope::registerPointCloud("tedge start", ps);
    pc_s->setPointColor({0.1, 0.8, 0.1});
    pc_s->setPointRadius(0.004);
    pc_s->resetTransform();
    auto* pc_e = polyscope::registerPointCloud("tedge end", pe);
    pc_e->setPointColor({0.9, 0.1, 0.1});
    pc_e->setPointRadius(0.003);
    pc_e->resetTransform();

    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    for (auto [a, b]: track) {
        size_t base = ns.size();
        ns.emplace_back(qv(a, 0), qv(a, 1), qv(a, 2));
        ns.emplace_back(qv(b, 0), qv(b, 1), qv(b, 2));
        es.push_back({base, base + 1});
    }
    auto* cn = polyscope::registerCurveNetwork("tedge tracks on the quad mesh", ns, es);
    cn->setRadius(0.0015);
    cn->resetTransform();

    auto* patch = polyscope::registerSurfaceMesh("quad patch", qv, qidx);
    patch->setShadeStyle(polyscope::MeshShadeStyle::Flat);
    patch->setEdgeWidth(1.);
    patch->addFaceScalarQuantity("tqid", tqid_of_quad)->setEnabled(true);

    std::println("[quad patch] singular seeds {} | inconsistent chains {} | walks failed {} | unreached tedges {} | wrong landing on a singular {} | junction vertices {} | track edges {} | unlabeled quads {}",
                 nseed, bad_chain, failed, unreached, sing_wrong, (int)pe.size(), (int)track.size(), unlabeled);
}

}
#endif
