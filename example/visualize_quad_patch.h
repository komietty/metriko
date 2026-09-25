#ifndef METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#define METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#include <set>
#include <queue>
#include <cmath>
#include <algorithm>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "metriko/core/tmesh/emesh.h"
#include "metriko/core/qex/common.h"

namespace metriko::visualizer {

// replay the t-mesh on the quad graph and colour the quads by tquad, with no geometry.
//
// singular vertices are q-vertices, so their quad vertex is exact. around a node the tedges are chained CCW by
// the tquad boundary order: rotating CCW from an outgoing tedge sweeps through the tquad on its left and reaches
// that tquad's previous boundary thalf; the slot advances by one quad edge, or two when the node lies inside a
// side. every tedge walks exactly x quad edges, and its landing vertex IS its far node. the only unknown is the
// rotation of the chain at a singular: the first singular is anchored through a tedge to another singular (only
// one rotation lands there), every other node is anchored by the first track arriving at it.
inline void visualize_quad_patch(
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
    // CCW chain of the tedges at a node starting from (te, fwd): (teid, fwd, slot offset). rotating CCW around
    // the node from an outgoing tedge sweeps through the tquad P on its left and reaches P's previous boundary
    // thalf, which ends at the node. this follows P's boundary order, so a tquad touching itself at the node
    // (two of its corners on one vertex) is handled; the slot advances by 2 when the node lies inside a side of P
    auto chain_at = [&](int nid, int te, bool fwd) {
        vec<std::tuple<int, bool, int>> ch = {{te, fwd, 0}};
        const int n = (int)node_te[nid].size();
        int off = 0;
        for (int it = 1; it < n; ++it) {
            const int c = te2th.at(te);
            const auto& th = tm.thalfs[fwd ? c : tm.thalfs[c].twid];   // P on its left, pointing away from the node
            const auto& tq = tm.tquads[th.tqid];
            auto cur = rg::find(tq.data, th.id, &Edata::thid);
            auto prv = circular_prev(tq.data, cur);
            const auto& th2 = tm.thalfs[prv->thid];              // ends at the node
            off += prv->side == cur->side ? 2 : 1;
            te  = th2.teid;
            fwd = !th2.cano;                                     // outgoing from the node: backwards along nids if cano
            ch.emplace_back(te, fwd, off);
        }
        return ch;
    };

    vec<double> tqid_of_quad(l, -1);
    std::set<std::pair<int, int>> track;
    vec<bool> done(tm.tedges.size(), false);
    umap<int, int> node_qv;   // node -> quad vertex
    struct Job { int teid; bool fwd; int s; int b; };
    std::queue<Job> jobs;

    auto sing_qv = [&](int nid) {   // quad vertex of a singular node, -1 otherwise
        auto* lv = std::get_if<HmLocOnV>(&tm.tnodes[nid]);
        return lv && singular(lv->id) && qv_of_vid.contains(lv->id) ? qv_of_vid.at(lv->id) : -1;
    };
    auto land = [&](int s, int b, int steps) {   // straight walk of `steps` edges from s through b: landing vertex or -1
        int a = s;
        for (int k = 0; k + 1 < steps; ++k) {
            if ((int)out[b].size() != 4) return -1;
            int c = straight(a, b);
            a = b;
            b = c;
        }
        return b;
    };
    auto ring_at = [&](int s, int first) {   // quad edge heads around s, CCW from `first`
        vec<int> ring = {first};
        while ((int)ring.size() < (int)out[s].size()) ring.push_back(rot(s, ring.back()));
        return ring;
    };

    int by_pair = 0, by_arrival = 0, failed = 0, sing_wrong = 0;
    // anchor a singular through a tedge to another singular: the rotation whose walk lands on it
    auto seed_by_pair = [&](int nid) {
        const int s = sing_qv(nid);
        auto ring = ring_at(s, out[s].front());
        const int n = (int)ring.size();
        auto ch = chain_at(nid, node_te[nid][0].first, node_te[nid][0].second);
        for (auto [te, fwd, off]: ch) {
            const auto& nids = tm.tedges[te].nids;
            int far = fwd ? nids.back() : nids.front();
            int t = sing_qv(far);
            if (t < 0 || far == nid) continue;
            vec<int> good;
            for (int r = 0; r < n; ++r) if (land(s, ring[(off + r) % n], xs.at(te)) == t) good.push_back(r);
            if (good.size() != 1) continue;
            node_qv[nid] = s;
            ++by_pair;
            for (auto [te2, fwd2, off2]: ch) jobs.push({te2, fwd2, s, ring[(off2 + good[0]) % n]});
            return true;
        }
        return false;
    };
    auto run = [&] {
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
            if (int t = sing_qv(nid); t >= 0 && t != b) { ++sing_wrong; continue; }   // landed on a singular: must agree
            if (node_qv.contains(nid)) continue;
            node_qv[nid] = b;
            if (sing_qv(nid) >= 0) ++by_arrival;
            auto ring = ring_at(b, a);   // slots CCW from the reverse of arrival
            for (auto [te, fwd, off]: chain_at(nid, jb.teid, !jb.fwd)) if (!done[te]) jobs.push({te, fwd, b, ring[off % ring.size()]});
        }
    };
    for (auto& [nid, ports]: node_te) if (sing_qv(nid) >= 0 && !node_qv.contains(nid)) { seed_by_pair(nid); run(); }
    int unanchored = 0;
    for (auto& [nid, ports]: node_te) if (sing_qv(nid) >= 0 && !node_qv.contains(nid)) ++unanchored;
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
    vec<double> s_vid, s_tedges, s_tracks, s_valence;   // per singular: vertex id, tedges at it, quad edges on a track, quad valence
    for (auto& [nid, v]: node_qv) {
        auto* lv = std::get_if<HmLocOnV>(&tm.tnodes[nid]);
        if (!(lv && singular(lv->id))) { pe.emplace_back(qv(v, 0), qv(v, 1), qv(v, 2)); continue; }
        ps.emplace_back(qv(v, 0), qv(v, 1), qv(v, 2));
        int ntrack = 0;
        for (int w: out[v]) if (track.contains(std::minmax(v, w))) ++ntrack;
        s_vid.push_back(lv->id);
        s_tedges.push_back((double)node_te[nid].size());
        s_tracks.push_back(ntrack);
        s_valence.push_back((double)out[v].size());
    }
    auto* pc_s = polyscope::registerPointCloud("tedge start", ps);
    pc_s->setPointColor({0.1, 0.8, 0.1});
    pc_s->setPointRadius(0.004);
    pc_s->addScalarQuantity("vertex id", s_vid);
    pc_s->addScalarQuantity("tedges", s_tedges);
    pc_s->addScalarQuantity("tracks", s_tracks);
    pc_s->addScalarQuantity("quad valence", s_valence);
    pc_s->resetTransform();
    for (size_t k = 0; k < s_vid.size(); ++k)
        if (s_tracks[k] != s_tedges[k]) std::println("[quad patch] singular vertex {}: {} tedges but {} quad edges on a track (quad valence {})", (int)s_vid[k], (int)s_tedges[k], (int)s_tracks[k], (int)s_valence[k]);
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

    // colours: greedy graph colouring of the patch adjacency, so neighbouring patches never share a colour index
    std::map<int, std::set<int>> adj;
    for (int i = 0; i < l; ++i)
    for (int j = 0; j < 4; ++j) {
        int a = qidx(i, j), b = qidx(i, (j + 1) % 4);
        int nb = de.at({b, a}).first;
        int ta = (int)tqid_of_quad[i], tb = (int)tqid_of_quad[nb];
        if (ta >= 0 && tb >= 0 && ta != tb) { adj[ta].insert(tb); adj[tb].insert(ta); }
    }
    vec<std::pair<int, int>> order;   // (degree, tqid), most constrained first
    for (auto& [t, ns]: adj) order.emplace_back((int)ns.size(), t);
    rg::sort(order, std::greater<>{});
    umap<int, int> colour_of;
    for (auto [deg, t]: order) {
        std::set<int> used;
        for (int u: adj[t]) if (colour_of.contains(u)) used.insert(colour_of.at(u));
        int c = 0;
        while (used.contains(c)) ++c;
        colour_of[t] = c;
    }
    // colour index as a face scalar; consecutive indices are spread by the golden ratio so that neighbours map far apart
    vec<double> colour(l, -1);
    for (int i = 0; i < l; ++i) {
        int t = (int)tqid_of_quad[i];
        if (t >= 0 && colour_of.contains(t)) colour[i] = std::fmod(colour_of.at(t) * 0.618033988749895, 1.);
    }

    auto* patch = polyscope::registerSurfaceMesh("quad patch", qv, qidx);
    patch->setShadeStyle(polyscope::MeshShadeStyle::Flat);
    patch->setEdgeWidth(1.);
    patch->addFaceScalarQuantity("tqid", tqid_of_quad);
    patch->addFaceScalarQuantity("patch colour", colour)->setEnabled(true);

    std::println("[quad patch] singulars anchored by a pair {} | by arrival {} | unanchored {} | walks failed {} | unreached tedges {} | wrong landing on a singular {} | junction vertices {} | track edges {} | unlabeled quads {}",
                 by_pair, by_arrival, unanchored, failed, unreached, sing_wrong, (int)pe.size(), (int)track.size(), unlabeled);
}

}
#endif
