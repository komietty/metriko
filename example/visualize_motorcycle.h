#ifndef METRIKO_VISUALIZE_MOTORCYCLE_H
#define METRIKO_VISUALIZE_MOTORCYCLE_H
#include "./visualize_common.h"
#include "metriko/core/hmesh/hmloc.h"
#include "metriko/core/hmesh/utilities.h"
#include "metriko/core/tmesh/motorcycle.h"

namespace metriko::visualizer {
inline void visualize_motorcycle_graph(
    const Mgrph& mg,
    const VecXc& uv,
    bool show = true
) {
    vec<glm::vec3> l0;
    vec<glm::vec3> mnodes;
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    vec<double> mcid;
    vec<double> mnid;
    vec<double> modes_val;

    vec<bool> mnodes_reserved = vec(mg.mnodes.size(), false);
    size_t counter = 0;

    auto get_loc_type_value = [](const HmLoc& loc) -> double {
        if (std::holds_alternative<HmLocOnV>(loc)) return 1.;
        if (std::holds_alternative<HmLocOnE>(loc)) return 2.;
        if (std::holds_alternative<HmLocOnP>(loc)) return 3.;
        return 0.;
    };

    for (const auto& c: mg.mcurvs) {
    for (const auto& s: c.sgmts) {
        auto iFr = s.fr_nid;
        auto iTo = s.to_nid;
        auto nFr = mg.mnodes[iFr];
        auto nTo = mg.mnodes[iTo];
        auto uv1 = get_face_uv(nFr.loc, s.face_id, mg.hm, mg.cf);
        auto uv2 = get_face_uv(nTo.loc, s.face_id, mg.hm, mg.cf);
        Row3d p1 = conversion_2d_3d(mg.hm.faces[s.face_id], uv, uv1);
        Row3d p2 = conversion_2d_3d(mg.hm.faces[s.face_id], uv, uv2);
        if (abs(uv1 - uv2) < EPS) { l0.emplace_back(p1.x(), p1.y(), p1.z()); }

        if (!mnodes_reserved[iFr]) {
            mnodes.emplace_back(p1.x(), p1.y(), p1.z());
            modes_val.emplace_back(get_loc_type_value(nFr.loc));
            mnid.emplace_back(iFr);
            mnodes_reserved[iFr] = true;
        }

        if (!mnodes_reserved[iTo]) {
            mnodes.emplace_back(p2.x(), p2.y(), p2.z());
            modes_val.emplace_back(get_loc_type_value(nTo.loc));
            mnid.emplace_back(iTo);
            mnodes_reserved[iTo] = true;
        }

        ns.emplace_back(p1.x(), p1.y(), p1.z());
        ns.emplace_back(p2.x(), p2.y(), p2.z());
        es.emplace_back(std::array{counter, counter + 1});
        mcid.emplace_back(c.id);
        counter += 2;
    }}

    {
        auto p = polyscope::registerPointCloud("zero len edge", l0);
        p->setMaterial("flat");
        p->setPointRadius(0.003);
    }

    {
        auto p = polyscope::registerPointCloud("mnodes", mnodes);
        p->addScalarQuantity("type", modes_val);
        p->addScalarQuantity("mnid", mnid);
        p->setMaterial("flat");
        p->setPointRadius(0.003);
    }

    auto c = polyscope::registerCurveNetwork("motorcycle graph", ns, es);
    c->setMaterial("flat");
    c->setColor(glm::vec4(.0, .0, .0, 1.));
    auto v_mcid = c->addEdgeScalarQuantity("mcid", mcid);
    v_mcid->setEnabled(true);
    v_mcid->setColorMap("magma");
    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.001);
    c->setMaterial("flat");
}

inline void visualize_node_adjacency(const Mgrph& mg, const VecXc& uv, bool show = true) {
    std::vector<glm::vec3> pts;
    std::vector<double> adj_order;
    std::vector<double> nid_list;

    for (int nid = 0; nid < mg.mnodes.size(); ++nid) {
        const auto& mn = mg.mnodes[nid];
        if (mn.adj.size() < 3) continue;
        std::vector<glm::vec3> local_pts;

        for (int i = 0; i < mn.adj.size(); ++i) {
            const auto& as = mn.adj[i];
            const auto& sg = mg.mcurvs[as.x()].sgmts[as.y()];

            bool is_outgoing = sg.fr_nid == nid;
            int fid = sg.face_id;

            Row3d pA = conversion_2d_3d(mg.hm.faces[fid], uv, get_face_uv(mg.mnodes[sg.fr_nid].loc, fid, mg.hm, mg.cf));
            Row3d pB = conversion_2d_3d(mg.hm.faces[fid], uv, get_face_uv(mg.mnodes[sg.to_nid].loc, fid, mg.hm, mg.cf));
            Row3d p0 = is_outgoing ? pA : pB;
            Row3d p1 = is_outgoing ? pB : pA;
            Row3d pt = p0 * 0.85 + p1 * 0.15;
            glm::vec3 gpt(pt.x(), pt.y(), pt.z());

            pts.push_back(gpt);
            local_pts.push_back(gpt);

            adj_order.push_back(i);
            nid_list.push_back(nid);
        }
    }

    auto pc = polyscope::registerPointCloud("CCW Adjacency Points", pts);
    pc->setEnabled(show);
    pc->setPointRadius(0.002);

    auto q_order = pc->addScalarQuantity("adj_index", adj_order);
    q_order->setEnabled(true);
    q_order->setColorMap("turbo");
    pc->addScalarQuantity("node_id", nid_list);

}
}

#endif
