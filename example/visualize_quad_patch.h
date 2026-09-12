#ifndef METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#define METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#include <queue>
#include <igl/remove_duplicate_vertices.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "metriko/core/tmesh/emesh.h"
#include "metriko/core/qex/common.h"

namespace metriko::visualizer {

// assign each extracted quad the id of the t-mesh patch (tquad) it belongs to,
// and show the result as a face scalar on its own surface mesh.
//
// the assignment is combinatorial: every tedge is a straight path of exactly x
// (its quantized length) edges on the quad graph, so the t-mesh is REPLAYED on
// the quad connectivity — no re-race, no geometric point location. tedges
// anchored at singularities are walked first (singular vertices sit exactly on
// grid corners); each finished walk fixes its far junction, which anchors the
// tedges meeting there. geometry is used only once per seed, to pick the
// initial direction at a singularity.
inline void visualize_quad_patch(
    const Hmesh& hm,
    const Emesh& tm,
    const VecXi& singular,
    const vec<qex::Qface>& qfaces
) {
    const int l = (int)qfaces.size();
    MatXd pos(l * 4, 3);
    MatXi idx(l, 4);
    for (int i = 0; i < l; i++)
    for (int j = 0; j < 4; j++) {
        pos.row(i * 4 + j) = qfaces[i].qhalfs[j].port1().pos;
        idx(i, j) = i * 4 + j;
    }

    // weld the quad soup (combinatorics only, positions unsmoothed)
    MatXd qv;
    VecXi SVI, SVJ;
    const double weps = 1e-7 * (pos.colwise().maxCoeff() - pos.colwise().minCoeff()).norm();
    igl::remove_duplicate_vertices(pos, weps, qv, SVI, SVJ);
    MatXi qidx(l, 4);
    for (int i = 0; i < l; ++i) for (int j = 0; j < 4; ++j) qidx(i, j) = SVJ(idx(i, j));

    // directed edge -> (quad, corner); vertex -> outgoing edge heads
    std::map<std::pair<int, int>, std::pair<int, int>> de;
    umap<int, vec<int>> out;
    for (int i = 0; i < l; ++i)
    for (int j = 0; j < 4; ++j) {
        int a = qidx(i, j), b = qidx(i, (j + 1) % 4);
        de[{a, b}] = {i, j};
        out[a].push_back(b);
    }
    // next outgoing edge head, rotating CCW around v from (v -> w)
    auto rot = [&](int v, int w) {
        auto [q, j] = de.at({v, w});
        return qidx(q, (j + 3) % 4);   // prev corner of the quad
    };
    // arrived via (t -> v): the straight continuation is two rotations
    auto straight = [&](int t, int v) { return rot(v, rot(v, t)); };

    auto find_qv = [&](const Row3d& p) {
        for (int i = 0; i < (int)qv.rows(); ++i)
            if ((Row3d(qv.row(i)) - p).norm() < weps * 10) return i;
        return -1;
    };

    umap<int, int> xs;                       // teid -> steps
    umap<int, int> te2th;                    // teid -> cano thalf id
    umap<int, std::pair<int, int>> sides;    // teid -> (left, right) tqid
    for (auto& th: tm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        xs[th.teid]    = (int)std::round(th.x);
        te2th[th.teid] = th.id;
        sides[th.teid] = {th.tqid, tm.thalfs[th.twid].tqid};
    }

    // node -> incident tedges (with which endpoint touches the node)
    umap<int, vec<std::pair<int, bool>>> node_te;
    for (auto& [teid, nids]: tm.live_tedges()) {
        node_te[nids.front()].emplace_back(teid, true);
        node_te[nids.back()].emplace_back(teid, false);
    }
    // outgoing 3d direction of a tedge at one of its endpoint nodes
    auto dir_at = [&](int teid, int nid) -> Row3d {
        auto& nids = tm.tedges[teid].nids;
        int nb = nids.front() == nid ? nids[1] : nids[nids.size() - 2];
        return (get_ptloc_pos(hm, tm.tnodes[nb]) - get_ptloc_pos(hm, tm.tnodes[nid])).normalized();
    };

    std::vector<double> quad_tqid(l, -1);
    std::set<std::pair<int, int>> bnd;       // traced (undirected) quad edges
    vec done(tm.tedges.size(), false);

    struct Job { int teid; bool fwd; int s; int b; };   // walk from quad vertex s, first head b
    std::queue<Job> jobs;

    // walk x quad edges; the cano thalf's tquad lies on the LEFT of the nids
    // direction (same convention as the cut annotation). returns the far vertex
    // and the tail of the arriving edge, or {-1, -1}
    auto walk = [&](const Job& jb) -> std::pair<int, int> {
        auto [lq, rq] = sides.at(jb.teid);
        if (!jb.fwd) std::swap(lq, rq);
        int steps = xs.at(jb.teid);
        int a = jb.s, b = jb.b;
        for (int k = 0; k < steps; ++k) {
            quad_tqid[de.at({a, b}).first] = lq;
            quad_tqid[de.at({b, a}).first] = rq;
            bnd.insert(std::minmax(a, b));
            if (k + 1 < steps) {
                if ((int)out[b].size() != 4) return {-1, -1};   // early singular: mismatch
                int c = straight(a, b);
                a = b;
                b = c;
            }
        }
        return {b, a};
    };

    auto lr_out = [&](int te, bool afn) {     // left/right of the OUTGOING direction
        auto [lq, rq] = sides.at(te);
        return afn ? std::pair(lq, rq) : std::pair(rq, lq);
    };

    // seeds: at each singular vertex, match its tedges to its quad edges JOINTLY.
    // independent per-tedge angle matching can send two tedges down the same
    // separatrix under distortion; instead anchor ONE pair by angle and assign
    // the rest by CCW chaining (every wedge at a singularity spans one slot)
    umap<int, vec<std::pair<int, bool>>> sing_te;
    for (auto& [teid, nids]: tm.live_tedges()) {
        if (xs.at(teid) <= 0 || nids.size() < 2) { done[teid] = true; continue; }
        for (bool fwd: {true, false}) {
            int nid = fwd ? nids.front() : nids.back();
            auto* lv = std::get_if<HmLocOnV>(&tm.tnodes[nid]);
            if (lv && singular(lv->id) != 0) sing_te[nid].emplace_back(teid, fwd);
        }
    }
    for (auto& [nid, ports]: sing_te) {
        int vid = std::get<HmLocOnV>(tm.tnodes[nid]).id;
        int s = find_qv(hm.verts[vid].pos());
        if (s < 0) continue;

        vec<int> ring = {out[s].front()};   // quad edges around s, CCW
        while ((int)ring.size() < (int)out[s].size()) ring.push_back(rot(s, ring.back()));
        if (ports.size() != ring.size())
            std::println("[tqid] singular {}: {} tedges vs {} quad edges", vid, ports.size(), ring.size());

        // anchor: the (tedge, slot) pair with the best angular agreement
        int a_p = 0, a_s = 0;
        double bd = -2;
        for (int p = 0; p < (int)ports.size(); ++p) {
            Row3d d3 = dir_at(ports[p].first, nid);
            for (int k = 0; k < (int)ring.size(); ++k) {
                double d = (Row3d(qv.row(ring[k])) - Row3d(qv.row(s))).normalized().dot(d3);
                if (d > bd) { bd = d; a_p = p; a_s = k; }
            }
        }

        // the rest by CCW chaining from the anchor
        int  cur  = ports[a_p].first;
        bool afn  = ports[a_p].second;
        int  slot = a_s;
        jobs.push({cur, afn, s, ring[slot]});
        for (size_t it = 1; it < ports.size(); ++it) {
            int P = lr_out(cur, afn).first;
            auto nx = rg::find_if(ports, [&](auto& pr) {
                return pr.first != cur && lr_out(pr.first, pr.second).second == P;
            });
            if (nx == ports.end()) break;
            slot = (slot + 1) % (int)ring.size();
            jobs.push({nx->first, nx->second, s, ring[slot]});
            cur = nx->first;
            afn = nx->second;
        }
    }

    int failed = 0;
    while (!jobs.empty()) {
        Job jb = jobs.front();
        jobs.pop();
        if (done[jb.teid]) continue;
        done[jb.teid] = true;
        auto [vj, tail] = walk(jb);
        if (vj < 0) { ++failed; continue; }

        // the walk fixed the far junction: anchor the tedges meeting there.
        // slots are assigned combinatorially: going CCW around the node,
        // consecutive tedges share a patch (left of the current = right of the
        // next), and that patch spans one slot when the node is its corner,
        // two when its boundary runs straight through
        auto& nids = tm.tedges[jb.teid].nids;
        int nid = jb.fwd ? nids.back() : nids.front();
        if ((int)out[vj].size() != 4) continue;   // singular end: those seed themselves

        int ring[4];   // quad-edge slots, CCW from the reverse of arrival
        ring[0] = tail;
        for (int k = 1; k < 4; ++k) ring[k] = rot(vj, ring[k - 1]);

        auto side_in = [&](int te, int P) {       // side index of te's thalf bounding P
            int c = te2th.at(te);
            int h = tm.thalfs[c].tqid == P ? c : tm.thalfs[c].twid;
            return tm.tquads[P].side_of(tm.thalfs[h]);
        };

        auto ports = node_te[nid];
        int  cur  = jb.teid;
        bool afn  = !jb.fwd;   // the arrival tedge, oriented outgoing at this node
        int  slot = 0;
        for (size_t it = 1; it < ports.size(); ++it) {
            int P = lr_out(cur, afn).first;       // the wedge CCW after cur
            auto nx = rg::find_if(ports, [&](auto& pr) {
                return pr.first != cur && lr_out(pr.first, pr.second).second == P;
            });
            if (nx == ports.end()) break;         // inconsistent incidence data
            slot += side_in(cur, P) == side_in(nx->first, P) ? 2 : 1;
            if (!done[nx->first]) jobs.push({nx->first, nx->second, vj, ring[slot % 4]});
            cur = nx->first;
            afn = nx->second;
        }
    }
    int unreached = 0;
    for (auto& [teid, nids]: tm.live_tedges()) if (!done[teid]) ++unreached;
    if (failed || unreached) std::println("[tqid] replay: {} failed, {} unreached", failed, unreached);

    // flood fill the patch interiors, never crossing a traced edge
    std::queue<int> que;
    for (int i = 0; i < l; ++i) if (quad_tqid[i] >= 0) que.push(i);
    while (!que.empty()) {
        int q = que.front(); que.pop();
        for (int j = 0; j < 4; ++j) {
            int a = qidx(q, j), b = qidx(q, (j + 1) % 4);
            if (bnd.contains(std::minmax(a, b))) continue;
            int nb = de.at({b, a}).first;
            if (quad_tqid[nb] < 0) { quad_tqid[nb] = quad_tqid[q]; que.push(nb); }
        }
    }
    if (int u = (int)rg::count(quad_tqid, -1.); u > 0) std::println("[tqid] {} quads unlabeled", u);

    // greedy-color the patch adjacency graph so neighboring patches never share
    // a color (raw tqids give near-identical colormap values to neighbors)
    umap<int, int> color;
    for (auto& tq: tm.live_tquads()) {
        std::set<int> used;
        for (auto& d: tq.data) {
            int nb = tm.thalfs[tm.thalfs[d.thid].twid].tqid;
            if (auto it = color.find(nb); it != color.end()) used.insert(it->second);
        }
        int c = 0;
        while (used.contains(c)) ++c;
        color[tq.id] = c;
    }
    std::vector<double> quad_color(l, -1);
    for (int i = 0; i < l; ++i)
        if (quad_tqid[i] >= 0) quad_color[i] = color[(int)quad_tqid[i]];

    auto* surf = polyscope::registerSurfaceMesh("quad patch", qv, qidx);
    surf->setShadeStyle(polyscope::MeshShadeStyle::Flat);
    surf->setEdgeWidth(1.);
    surf->addFaceScalarQuantity("tqid", quad_tqid);
    surf->addFaceScalarQuantity("patch color", quad_color)->setEnabled(true);

    // patch boundaries: the traced (replayed) quad edges
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t c = 0;
        for (auto& [a, b]: bnd) {
            ns.emplace_back(qv(a, 0), qv(a, 1), qv(a, 2));
            ns.emplace_back(qv(b, 0), qv(b, 1), qv(b, 2));
            es.push_back({c, c + 1});
            c += 2;
        }
        auto* cn = polyscope::registerCurveNetwork("quad patch boundary", ns, es);
        cn->setMaterial("flat");
        cn->setColor({0.05, 0.05, 0.05});
        cn->setRadius(0.0012);
        cn->resetTransform();
    }
}
}

#endif
