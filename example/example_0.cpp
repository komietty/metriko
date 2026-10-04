#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <igl/upsample.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/nvec/face_rosy_field.h"
#include "metriko/igm/parameterization.h"
#include "metriko/quantization/quantization.h"
#include "metriko/hmesh/hpath.h"
#include "metriko/tmesh/emesh.h"
#include "metriko/tutte/tutte.h"
#include "common_io.h"
#include "common_visualizer.h"
#include "cleanup.h"

using namespace metriko;
static int N = 4;
static MatXd V;
static MatXi F;

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    //cleanup::cleanup_mesh(V, F);
    cleanup::decimate_and_clean(V, F, 1000000);
    V /= grid_unit(V, std::stod(argv[2]));   // grid units: one quad edge = 1
    Hmesh hm(V, F);

    FaceRosyField rawf(hm, N, FieldType::Smoothest);
    auto seam = compute_seam(rawf);
    auto cutm = compute_cut_mesh(hm, seam);
    auto cmbf = compute_combbed_field(rawf, seam);

    MatXd ext(hm.nF, 3 * N);
    for (Face f: hm.faces) {
    for (int k = 0; k < N; ++k) {
        complex c = cmbf->field(f.id, k);
        ext.block(f.id, 3 * k, 1, 3) = (c.real() * f.basisX() + c.imag() * f.basisY()).normalized();
    }}

    std::cout << "start params" << std::endl;

    RosyParameterization rp(hm, *cutm, ext, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();

    hm.cfn.resize(hm.nF * 3);
    for (const Face f: hm.faces) {
        hm.cfn(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        hm.cfn(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        hm.cfn(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }

    const std::string cache = std::format("{}.{}.cache", argv[1], argv[2]);
    save_cache(cache, hm.cfn, cmbf->matching, cmbf->singular, seam);
    std::println("saved cache: {}", cache);
    /* */

    visualizer::visualize_init();
    auto surf = visualizer::visualize_mesh(hm);
    visualizer::visualize_frosy_field(surf, rawf, *cmbf);
    visualizer::visualize_seam(hm, seam);
    visualizer::visualize_wrong_cones(surf, hm, cmbf->singular);

    polyscope::show();
    return 0;
}
