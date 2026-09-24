#ifndef METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#define METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#include <set>
#include <queue>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "metriko/core/tmesh/emesh.h"
#include "metriko/core/tutte/tutte.h"
#include "metriko/core/qex/common.h"

namespace metriko::visualizer {

// 1: show the tedge end nodes at their quad vertices: singular nodes as starts, T-junctions as ends.
//    node -> embedding-mesh vertex comes from the cut annotation (hdata: order 0 tail / last order head of the
//    cano thalf), embedding vertex -> quad vertex from the vertex-q-vertex port carrying that vid
// 2: walk every tedge along the quad graph from its end nodes, 3: flood fill the tquad patches between the tracks
inline void visualize_quad_patch(
    const Emesh& tm,
    const VecXi& singular,
    const vec<HalfData>& hdata,
    const VecXc& cfn,   // the uv qex was run on
    const vec<qex::Qface>& qfaces,
    const MatXd& qv,    // displayed quad mesh vertices (from visualize_qfaces)
    const MatXi& qidx   // displayed quad mesh faces, corner (i, j) <-> qfaces[i].qhalfs[j].port1()
) {
    const int l = (int)qfaces.size();

    // vertex -> outgoing edge heads; embedding vertex -> quad vertex (through the vertex-q-vertex ports)
    umap<int, int> qv_of_vid;
    std::map<std::pair<int, int>, vec<std::pair<int, int>>> corners_of;   // (vid, fid) -> quad corners (i, j) whose port sits there
    std::map<std::pair<int, int>, std::pair<int, int>> de;                // directed quad edge (a, b) -> (quad, corner); the quad lies on its left
    umap<int, int> valence;
    for (int i = 0; i < l; ++i)
    for (int j = 0; j < 4; ++j) {
        const auto& p = qfaces[i].qhalfs[j].port1();
        de[{qidx(i, j), qidx(i, (j + 1) % 4)}] = {i, j};
        valence[qidx(i, j)]++;
        if (p.vid < 0) continue;
        qv_of_vid[p.vid] = qidx(i, j);
        corners_of[{p.vid, p.fid}].emplace_back(i, j);
    }
    auto rot      = [&](int v, int w) { auto [q, j] = de.at({v, w}); return qidx(q, (j + 3) % 4); };   // next head CCW around v
    auto straight = [&](int t, int v) { return rot(v, rot(v, t)); };                                 // arrived t -> v: continue

    // every tedge end node -> embedding vertex, from the cano thalf's cut halfedges (order 0 tail, last order head).
    // a node is a start if it sits on a singular vertex, otherwise an end (T-junction)
    std::set<int> starts, ends;
    int missing = 0, nsing = 0;
    for (const auto& d: hdata) {
        const auto& th = tm.thalfs[d.thid];
        if (!th.cano) continue;
        const auto& nids = tm.tedges[th.teid].nids;
        for (auto [k, nid, v]: {std::tuple{0, nids.front(), d.half.tail().id}, std::tuple{(int)nids.size() - 2, nids.back(), d.half.head().id}}) {
            if (d.order != k) continue;
            auto* lv = std::get_if<HmLocOnV>(&tm.tnodes[nid]);
            bool sing = lv && singular(lv->id) != 0;
            if (!qv_of_vid.contains(v)) { ++missing; continue; }
            (sing ? starts : ends).insert(qv_of_vid.at(v));
        }
    }
    for (int i = 0; i < singular.size(); ++i) if (singular(i) != 0) ++nsing;

    // 1: start (singular) and end (T-junction) quad vertices
    auto to_glm = [&](const std::set<int>& vs) {
        std::vector<glm::vec3> ps;
        for (int v: vs) ps.emplace_back(qv(v, 0), qv(v, 1), qv(v, 2));
        return ps;
    };
    auto* pc_s = polyscope::registerPointCloud("tedge start", to_glm(starts));
    pc_s->setPointColor({0.1, 0.8, 0.1});
    pc_s->setPointRadius(0.004);
    pc_s->resetTransform();
    auto* pc_e = polyscope::registerPointCloud("tedge end", to_glm(ends));
    pc_e->setPointColor({0.9, 0.1, 0.1});
    pc_e->setPointRadius(0.003);
    pc_e->resetTransform();

    // 2: walk every tedge on the quad graph. it leaves its end node along the cut halfedge h (tail = node); the port
    //    of the iso-line along h lives in h.face() (see generate_vqvert_qport) and points along h in that chart, and
    //    its q-edge is the first quad edge. from there the walk goes straight until it lands on another node vertex.
    //    the quad left of the walk belongs to the thalf's tquad, the quad right of it to the twin's
    std::set<int> node_qvs = starts;
    node_qvs.insert(ends.begin(), ends.end());
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::set<std::pair<int, int>> track;
    vec<double> tqid_of_quad(l, -1);
    int no_port = 0, stuck = 0;
    auto walk_tedge = [&](Half h, int tq_left, int tq_right) {
        const int v = h.tail().id;
        const complex e = cfn(h.prev().crnr().id) - cfn(h.next().crnr().id);   // head - tail in h.face()'s chart
        int a = -1, b = -1;
        for (auto [i, j]: corners_of[{v, h.face().id}]) {
            const complex d = qfaces[i].qhalfs[j].port1().dir;
            if (std::abs((d * std::conj(e)).imag()) > EPS * std::abs(e) || (d * std::conj(e)).real() <= 0) continue;
            a = qidx(i, j);
            b = qidx(i, (j + 1) % 4);
            break;
        }
        if (a < 0) { ++no_port; return; }
        while (true) {
            tqid_of_quad[de.at({a, b}).first] = tq_left;
            tqid_of_quad[de.at({b, a}).first] = tq_right;
            if (track.insert(std::minmax(a, b)).second) {
                size_t base = ns.size();
                ns.emplace_back(qv(a, 0), qv(a, 1), qv(a, 2));
                ns.emplace_back(qv(b, 0), qv(b, 1), qv(b, 2));
                es.push_back({base, base + 1});
            }
            if (node_qvs.contains(b)) return;
            if (valence.at(b) != 4) { ++stuck; return; }   // irregular vertex that is not a node: inconsistent
            const int c = straight(a, b);
            a = b;
            b = c;
        }
    };
    for (const auto& d: hdata) {
        const auto& th = tm.thalfs[d.thid];
        if (!th.cano) continue;
        const int last = (int)tm.tedges[th.teid].nids.size() - 2;
        const int tq_cano = th.tqid, tq_twin = tm.thalfs[th.twid].tqid;
        if (d.order == 0)    walk_tedge(d.half,        tq_cano, tq_twin);   // forward along nids
        if (d.order == last) walk_tedge(d.half.twin(), tq_twin, tq_cano);   // backward along nids
    }

    // 3: flood the patch interiors without crossing a track
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

    auto* patch = polyscope::registerSurfaceMesh("quad patch", qv, qidx);
    patch->setShadeStyle(polyscope::MeshShadeStyle::Flat);
    patch->setEdgeWidth(1.);
    patch->addFaceScalarQuantity("tqid", tqid_of_quad)->setEnabled(true);
    auto* cn = polyscope::registerCurveNetwork("tedge tracks on the quad mesh", ns, es);
    cn->setRadius(0.0015);
    cn->resetTransform();
    std::println("[quad patch] singular vertices {} | start vertices {} | end vertices {} | end nodes without a vertex-q-vertex {} | track edges {} | tedge ends without a matching port {} | walks stopped at an irregular non-node vertex {} | unlabeled quads {}", nsing, starts.size(), ends.size(), missing, es.size(), no_port, stuck, unlabeled);
}

}
#endif
