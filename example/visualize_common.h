#ifndef TMESH_MUT_COLLAPSE_TQUAD_CPP_VISUALIZE_COMMON_H
#define TMESH_MUT_COLLAPSE_TQUAD_CPP_VISUALIZE_COMMON_H
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/hmesh/hmesh.h"

namespace metriko::visualizer {

inline void visualize_init() {
    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
}

}

#endif
