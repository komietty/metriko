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

// a chain of locations on the mesh, as a polyline
inline void visualize_path(
    const Hmesh& hm,
    const vec<HmLoc>& path,
    const std::string& name,
    const bool show = true,
    const double scale = 0.0025
) {
    vec<glm::vec3> ps;
    vec<std::array<size_t, 2>> es;
    for (size_t k = 0; k < path.size(); ++k) {
        Row3d p = get_ptloc_pos(hm, path[k]);
        ps.emplace_back(p.x(), p.y(), p.z());
        if (k + 1 < path.size()) es.push_back({k, k + 1});
    }
    auto* cn = polyscope::registerCurveNetwork(name, ps, es);
    cn->setMaterial("flat");
    cn->setEnabled(show);
    cn->setRadius(scale);
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
    vec<double> teids;
    vec<double> tqids0;
    vec<double> tqids1;
    size_t c = 0;
    for (const auto& [id, nids]: tm.live_tedges()) {
        int tqid0 = -1;
        int tqid1 = -1;
        for (const Ehalf& th : tm.thalfs) {
            if (th.id != -1 && th.teid == id) {
                tqid0 = th.tqid;
                tqid1 = tm.thalfs[th.twid].tqid;
                break;
            }
        }

        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            Row3d a = get_ptloc_pos(hm, tm.tnodes[nids[k]]);
            Row3d b = get_ptloc_pos(hm, tm.tnodes[nids[k + 1]]);
            ns.emplace_back(a.x(), a.y(), a.z());
            ns.emplace_back(b.x(), b.y(), b.z());
            es.push_back({c, c + 1});
            c += 2;
            teids.push_back(id);
            tqids0.push_back(tqid0);
            tqids1.push_back(tqid1);
        }
    }
    auto* cn = polyscope::registerCurveNetwork("tedges_collapsed", ns, es);
    cn->setMaterial("flat");
    cn->setEnabled(show);
    cn->addEdgeScalarQuantity("teid", teids)->setEnabled(true);
    cn->addEdgeScalarQuantity("tqid0", tqids0)->setEnabled(false);
    cn->addEdgeScalarQuantity("tqid1", tqids1)->setEnabled(false);
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
    const double scale = 0.001,
    const std::string& append = "",
    const vec<int>& only = {}   // when given, draw just these tquads
) {

    for (const auto& [id, data] : tm.live_tquads()) {
        if (!only.empty() && !rg::contains(only, id)) continue;
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
        auto* cn = polyscope::registerCurveNetwork(std::format("tq {} {:03}", append, id), ns, es);
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
}
#endif
