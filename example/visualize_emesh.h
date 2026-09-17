#ifndef METRIKO_VISUALIZE_EMESH_H
#define METRIKO_VISUALIZE_EMESH_H
#include "./visualize_common.h"
#include "metriko/core/tmesh/emesh.h"

namespace metriko::visualizer {
// allowed corridor of a shortest-path query: the admissible sub-segment [r0, r1] of every edge
inline void visualize_allowed_range(
    const Hmesh& hm,
    const vec<std::tuple<int, double, double>>& allowed,
    const std::string& name = "allowed range",
    const bool show = true,
    const double scale = 0.0015
) {
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    vec<double> ids, spans;
    for (auto& [eid, r0, r1]: allowed) {
        Edge  e = hm.edges[eid];
        Row3d a = e.lerp(r0);
        Row3d b = e.lerp(r1);
        es.push_back({ns.size(), ns.size() + 1});
        ns.emplace_back(a.x(), a.y(), a.z());
        ns.emplace_back(b.x(), b.y(), b.z());
        ids.push_back(eid);
        spans.push_back(r1 - r0);
    }
    auto* cn = polyscope::registerCurveNetwork(name, ns, es);
    cn->setMaterial("flat");
    cn->setEnabled(show);
    cn->setRadius(scale);
    cn->addEdgeScalarQuantity("eid", ids)->setEnabled(true);
    cn->addEdgeScalarQuantity("span", spans);
    cn->resetTransform();
}

inline void visualize_tedge_mut_collapsed(
    const Hmesh& hm,
    const Emesh& tm,
    const bool show = true,
    const double scale = 0.001
) {
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    vec<double> ids;
    size_t c = 0;
    for (const auto& [id, nids]: tm.live_tedges()) {
        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            Row3d a = get_ptloc_pos(hm, tm.tnodes[nids[k]]);
            Row3d b = get_ptloc_pos(hm, tm.tnodes[nids[k + 1]]);
            ns.emplace_back(a.x(), a.y(), a.z());
            ns.emplace_back(b.x(), b.y(), b.z());
            es.push_back({c, c + 1});
            c += 2;
            ids.push_back(id);
        }
    }
    auto* cn = polyscope::registerCurveNetwork("tedges_collapsed", ns, es);
    cn->setMaterial("flat");
    cn->setEnabled(show);
    cn->addEdgeScalarQuantity("teid", ids)->setEnabled(true);
    cn->setRadius(scale);
    cn->resetTransform();
}

inline void visualize_tedge_mut_snapped(
    const Hmesh& hm,
    const Emesh& tm,
    const bool show = true,
    const double scale = 0.001
) {
    // quantized length per tedge (both thalfs of a tedge carry the same x)
    umap<int, double> x_of;
    for (const auto& th: tm.thalfs) if (th.id != -1) x_of[th.teid] = th.x;

    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::vector<double> ids, xs;
    size_t c = 0;
    for (const auto& [id, nids]: tm.live_tedges()) {
        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            Row3d a = get_ptloc_pos(hm, tm.tnodes[nids[k]]);
            Row3d b = get_ptloc_pos(hm, tm.tnodes[nids[k + 1]]);
            ns.emplace_back(a.x(), a.y(), a.z());
            ns.emplace_back(b.x(), b.y(), b.z());
            es.push_back({c, c + 1});
            c += 2;
            ids.push_back(id);
            xs.push_back(x_of.contains(id) ? x_of.at(id) : -1);
        }
    }
    auto* cn = polyscope::registerCurveNetwork("tedges_snapped", ns, es);
    cn->setMaterial("flat");
    cn->setEnabled(show);
    cn->addEdgeScalarQuantity("teid", ids)->setEnabled(true);
    cn->addEdgeScalarQuantity("x", xs);
    cn->setRadius(scale);
    cn->resetTransform();
}

inline void visualize_tquad_mut_collapsed(
    const Hmesh& hm,
    const Emesh& tm,
    const bool show = true,
    const double scale = 0.001
) {

    for (const auto& [id, data] : tm.live_tquads()) {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> eside, ex, ey, er, ethid;   // per-edge (thalf) params
        size_t c = 0;
        for (const Edata& d : data) {
            const Ehalf& th = tm.thalfs[d.thid];
            const Eedge& te = tm.tedges[th.teid];
            for (size_t i = 0; i + 1 < te.nids.size(); ++i) {
                Row3d a = get_ptloc_pos(hm, tm.tnodes[te.nids[i]]);
                Row3d b = get_ptloc_pos(hm, tm.tnodes[te.nids[i + 1]]);
                ns.emplace_back(a.x(), a.y(), a.z());
                ns.emplace_back(b.x(), b.y(), b.z());
                es.push_back({c, c + 1}); c += 2;
                eside.push_back(d.side);
                ex.push_back(th.x);
                ey.push_back(th.x > 0 ? 1 : 0);
                er.push_back(th.r);
                ethid.push_back(d.thid);
            }
        }
        if (ns.empty()) continue;
        auto* cn = polyscope::registerCurveNetwork(std::format("tq {:03}", id), ns, es);
        cn->addEdgeScalarQuantity("side", eside);
        auto cx = cn->addEdgeScalarQuantity("x", ex);
        auto cy = cn->addEdgeScalarQuantity("y", ey);
        cn->setEnabled(show);
        cn->addEdgeScalarQuantity("r", er);
        cn->addEdgeScalarQuantity("thid", ethid);
        cn->setMaterial("flat");
        cn->setRadius(scale);
        cn->resetTransform();
    }
}

inline void visualize_non_snapped_tnodes(
    const Hmesh& hm,
    const Emesh& tm,
    const bool show = true,
    const double scale = 0.003
) {
    std::vector<glm::vec3> ps;
    std::vector<double> type, ids, nids_;
    std::set<int> seen;   // shared nodes (junctions/crossings) appear in several chains
    for (const auto& [teid, nids]: tm.tedges) {
        if (teid == -1) continue;
        for (int nid: nids) {
            if (!seen.insert(nid).second) continue;
            double t = -1, id = -1;
            std::visit(overloaded{
                [&](const HmLocOnE& e) { t = 0; id = e.id; },
                [&](const HmLocOnH& h) { t = 1; id = h.id; },
                [&](const HmLocOnF& f) { t = 2; id = f.id; },
                [&](const auto&)       {},
            }, tm.tnodes[nid]);
            if (t < 0) continue;   // OnV: snapped, skip
            Row3d p = get_ptloc_pos(hm, tm.tnodes[nid]);
            ps.emplace_back(p.x(), p.y(), p.z());
            type.push_back(t);
            ids.push_back(id);
            nids_.push_back(nid);
        }
    }
    auto* pc = polyscope::registerPointCloud("unsnapped tnodes", ps);
    pc->setEnabled(show);
    pc->addScalarQuantity("type (0:E 1:H 2:F)", type)->setEnabled(true);
    pc->addScalarQuantity("elem id", ids);
    pc->addScalarQuantity("node id", nids_);
    pc->setPointRadius(scale);
}

inline void visualize_face_collinear_error(
    const Hmesh& hm,
    const Emesh& tm,
    const bool show = true
) {
    std::set<int> bad;   // hm face ids violating the snap_1 criterion
    for (auto& tq: tm.live_tquads()) {
        for (int side = 0; side < 4; side++) {
            std::set<int>  nids;
            umap<int, int> count;
            for (int thid: tq.thids(side))
                for (int nid: tm.tedges[tm.thalfs[thid].teid].nids) nids.insert(nid);
            for (int nid: nids) {
                if (auto* l = std::get_if<HmLocOnV>(&tm.tnodes[nid]))
                    for (Face f: hm.verts[l->id].adjHalfs() | vw::transform(&Half::face)) count[f.id]++;
            }
            for (auto& [fid, c]: count) if (c >= 3) bad.insert(fid);
        }}
    std::println("[collinear] {} faces violate snap_1", bad.size());
    // TEMP debug: which tedge / position each vertex of a violating face belongs to
    for (int fid: bad) {
        std::print("[collinear] face {}:", fid);
        for (Vert v: hm.faces[fid].verts())
            for (auto& [teid, nids]: tm.live_tedges())
                for (int k = 0; k < nids.size(); ++k)
                    if (auto* l = std::get_if<HmLocOnV>(&tm.tnodes[nids[k]]); l && l->id == v.id) std::print(" V{}=te{}[{}/{}]", v.id, teid, k, nids.size() - 1);
        std::println("");
    }

    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    size_t c = 0;
    for (int fid: bad) {
        for (Half h: hm.faces[fid].adjHalfs()) {
            Row3d a = h.tail().pos();
            Row3d b = h.head().pos();
            ns.emplace_back(a.x(), a.y(), a.z());
            ns.emplace_back(b.x(), b.y(), b.z());
            es.push_back({c, c + 1});
            c += 2;
        }
    }
    if (!ns.empty()) {
        auto* cn = polyscope::registerCurveNetwork("collinear faces", ns, es);
        cn->setMaterial("flat");
        cn->setEnabled(show);
        cn->setColor({1., 0.2, 0.1});
        cn->setRadius(0.0015);
        cn->resetTransform();
    }

    // TEMP debug: vertices of the violating faces, marked whether they are junctions
    {
        std::vector<glm::vec3> ps; std::vector<double> vids, junction;
        for (int fid: bad)
        for (Vert v: hm.faces[fid].verts()) {
            int hits = 0;
            for (auto& [teid, nids]: tm.live_tedges())
                for (int k = 0; k < nids.size(); ++k)
                    if (auto* l = std::get_if<HmLocOnV>(&tm.tnodes[nids[k]]); l && l->id == v.id && (k == 0 || k + 1 == nids.size())) ++hits;
            ps.emplace_back(v.pos().x(), v.pos().y(), v.pos().z());
            vids.push_back(v.id);
            junction.push_back(hits >= 2);   // endpoint of two or more tedges
        }
        if (!ps.empty()) {
            auto* pc = polyscope::registerPointCloud("collinear face vertices", ps);
            pc->addScalarQuantity("vid", vids);
            pc->addScalarQuantity("junction", junction)->setEnabled(true);
            pc->setPointRadius(0.002);
            pc->setEnabled(show);
        }
    }
}
}

#endif
