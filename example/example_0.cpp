#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <igl/upsample.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/hmesh/hpath.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "metriko/core/tutte/tutte.h"
#include "common.h"
#include "visualize_hmesh.h"
#include "visualize_vectorfield.h"

using namespace metriko;
static int N = 4;
static VecXc uv2;
static MatXd V;
static MatXi F;

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);

    FaceRosyField rawf(hm, N, FieldType::Smoothest);
    rawf.computeMatching(MatchingType::Principal);
    auto seam = compute_seam(rawf);
    auto cutm = compute_cut_mesh(hm, seam);
    auto cmbf = compute_combbed_field(rawf, seam);

    MatXd ext(hm.nF, 3 * N);
    for (Face f: hm.faces) {
    for (int k = 0; k < N; ++k) {
        complex c = cmbf->field(f.id, k);
        ext.block(f.id, 3 * k, 1, 3) = (c.real() * f.basisX() + c.imag() * f.basisY()).normalized();
    }}

    RosyParameterization rp(hm, *cutm, ext, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();

    uv2.resize(hm.nF * 3);
    for (const Face f: hm.faces) {
        uv2(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        uv2(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        uv2(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }

    const std::string cache = std::format("{}.{}.cache", argv[1], argv[2]);
    save_cache(cache, uv2, cmbf->matching, cmbf->singular, seam);
    std::println("saved cache: {}", cache);

    visualizer::visualize_init();
    auto surf = visualizer::visualize_mesh_with_uv(hm.pos, hm.idx, uv2);
    visualizer::visualize_frosy_field(surf, hm, rawf, *cmbf);
    visualizer::visualize_seam(hm, seam);

    polyscope::show();
    return 0;
}
