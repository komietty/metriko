#ifndef METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#define METRIKO_EXAMPLE_VISUALIZE_QUAD_PATCH_H
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "metriko/patching.h"

namespace metriko::visualizer {

// show a labeled quad mesh: singular / junction nodes, tedge tracks, and the patches coloured by tquad
inline void visualize_quad_patch(
    const MatXd& qv,             // quad mesh vertices
    const MatXi& qidx,           // quad mesh faces
    const QuadPatch& patch  // labels of qidx, from label_quad_patches
) {
    auto at = [&](int v) { return glm::vec3(qv(v, 0), qv(v, 1), qv(v, 2)); };

    std::vector<glm::vec3> ps, pe;
    vec<double> s_vid, s_tedges, s_tracks, s_valence;
    for (auto& s: patch.singulars) {
        ps.push_back(at(s.qv));
        s_vid.push_back(s.vid);
        s_tedges.push_back(s.tedges);
        s_tracks.push_back(s.tracks);
        s_valence.push_back(s.valence);
        if (s.tracks != s.tedges) std::println("[quad patch] singular vertex {}: {} tedges but {} quad edges on a track (quad valence {})", s.vid, s.tedges, s.tracks, s.valence);
    }
    for (int v: patch.junctions) pe.push_back(at(v));
    auto* pc_s = polyscope::registerPointCloud("tedge start", ps);
    pc_s->setPointColor({0.1, 0.8, 0.1});
    pc_s->setPointRadius(0.002);
    pc_s->addScalarQuantity("vertex id", s_vid);
    pc_s->addScalarQuantity("tedges", s_tedges);
    pc_s->addScalarQuantity("tracks", s_tracks);
    pc_s->addScalarQuantity("quad valence", s_valence);
    pc_s->setEnabled(false);
    pc_s->resetTransform();
    auto* pc_e = polyscope::registerPointCloud("tedge end", pe);
    pc_e->setPointColor({0.9, 0.1, 0.1});
    pc_e->setPointRadius(0.0015);
    pc_e->resetTransform();
    pc_e->setEnabled(false);

    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    for (auto [a, b]: patch.track) {
        size_t base = ns.size();
        ns.push_back(at(a));
        ns.push_back(at(b));
        es.push_back({base, base + 1});
    }
    auto* cn = polyscope::registerCurveNetwork("tedge tracks on the quad mesh", ns, es);
    cn->setMaterial("flat");
    cn->setColor({0., 0., 0.});
    cn->setRadius(0.001);
    cn->resetTransform();

    auto* surf = polyscope::registerSurfaceMesh("quad patch", qv, qidx);
    surf->setShadeStyle(polyscope::MeshShadeStyle::Smooth);
    surf->setEdgeWidth(1.);
    surf->addFaceScalarQuantity("tqid", patch.tqid_of_quad);
    auto* pcol = surf->addFaceScalarQuantity("patch colour", patch.colour);
    pcol->setColorMap("coolwarm");
    pcol->setEnabled(true);

    if (!patch.ok())
        std::println("[quad patch] anchors {} | singulars with no consistent rotation {} | with several consistent rotations {} | unreached tedges {} | junction vertices {} | track edges {} | unlabeled quads {}",
                     patch.anchors, patch.no_rotation, patch.several, patch.unreached, (int)patch.junctions.size(), (int)patch.track.size(), patch.unlabeled);
}

}
#endif
