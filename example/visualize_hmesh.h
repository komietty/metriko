#ifndef TMESH_MUT_COLLAPSE_TQUAD_CPP_VISUALIZE_HMESH_H
#define TMESH_MUT_COLLAPSE_TQUAD_CPP_VISUALIZE_HMESH_H
#include "./visualize_common.h"

namespace metriko::visualizer {
template <typename PType, typename IType>
polyscope::SurfaceMesh* visualize_mesh(
    const PType& pos,
    const IType& idx,
    const bool show = true,
    const std::string& name = "mesh"
) {
    auto* surf = polyscope::registerSurfaceMesh(name, pos, idx);
    surf->setEdgeWidth(0.7);
    surf->setEnabled(show);
    surf->setMaterial("flat");
    surf->setSurfaceColor(glm::vec3(0.3, 0.3, 0.3));
    return surf;
}

template <typename PType, typename IType, typename UType>
polyscope::SurfaceMesh* visualize_mesh_with_uv(
    const PType& pos,
    const IType& idx,
    const UType& uv,
    const bool show_hm = true,
    const bool show_uv = true,
    const std::string& name = "mesh"
) {
    auto* surf = polyscope::registerSurfaceMesh(name, pos, idx);
    auto* prms = surf->addParameterizationQuantity("uv_param", uv);
    prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
    prms->setCheckerSize(1.);
    prms->setEnabled(show_uv);
    surf->setEdgeWidth(0.7);
    surf->setEnabled(show_hm);
    surf->setMaterial("flat");
    surf->setSurfaceColor(glm::vec3(0.3, 0.3, 0.3));
    return surf;
}

inline void visualize_seam(
    const Hmesh& hmesh,
    const vec<bool>& seam,
    const VecXi& matching = VecXi(),
    const std::string& name = "seam",
    const bool show = true
) {
    bool use_matching = matching.rows() > 0;
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    vec<double> ms;
    size_t counter = 0;
    for (auto e: hmesh.edges) {
        if (seam[e.id]) {
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{counter, counter + 1});
            if (use_matching) ms.emplace_back(matching[e.id]);
            counter += 2;
        }
    }
    auto c = polyscope::registerCurveNetwork(name, ns, es);
    if(use_matching) c->addEdgeScalarQuantity("matching", ms);
    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.001);
}
}

#endif
