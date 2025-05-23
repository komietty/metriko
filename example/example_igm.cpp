#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include <igl/readOBJ.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>

using namespace metriko;
int N = 4;
MatXd V;
MatXi F;

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    auto mesh = std::make_unique<Hmesh>(V, F);
    auto rawf = std::make_unique<FaceRosyField>(*mesh, N, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Curl);
    auto seam = compute_seam(*rawf);
    auto cutm = compute_cut_mesh(*mesh, seam);
    auto cmbf = compute_combbed_field(*rawf, seam);
    MatXd extf = compute_extrinsic_field(*mesh, *cmbf, N);

    RosyParameterization rp(*mesh, *cutm, extf, cmbf->singular, cmbf->matching, seam, N);
    rp.seamless = true;
    rp.localInjectivity = false;
    rp.roundSeams = false;
    rp.verbose = false;

    MatXd uv(mesh->nF * 3, 2);
    rp.setup();
    rp.integ();
    for (const Face f: mesh->faces) {
        uv.row(f.id * 3 + 0) << rp.cfn(f.id, 0), rp.cfn(f.id, 1);
        uv.row(f.id * 3 + 1) << rp.cfn(f.id, 4), rp.cfn(f.id, 5);
        uv.row(f.id * 3 + 2) << rp.cfn(f.id, 8), rp.cfn(f.id, 9);
    }

    /// ---- visualize mesh ---- ///
    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    const auto surf = polyscope::registerSurfaceMesh("mesh", V, F);
    const auto prms = surf->addParameterizationQuantity("params", uv);
    prms->setEnabled(true);
    prms->setStyle(polyscope::ParamVizStyle::GRID);
    prms->setCheckerSize(1);

    ///--- visuailize singularities ---///
    std::vector<glm::vec3> singPos;
    std::vector<double> singVal;
    for (Vert v: mesh->verts) {
        if (int s = cmbf->singular[v.id]; s != 0) {
            Row3d p = v.pos();
            singPos.emplace_back(p.x(), p.y(), p.z());
            singVal.emplace_back(s);
        }
    }
    auto pc = polyscope::registerPointCloud("singularity", singPos);
    pc->addScalarQuantity("singular nums", singVal);
    pc->setEnabled(false);
    pc->setPointRadius(0.01);
    pc->resetTransform();

    ///--- visuailize seam ---///
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    size_t counterS = 0;
    for (auto e: mesh->edges) {
        if (seam[e.id]) {
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{counterS, counterS + 1});
            counterS += 2;
        }
    }
    auto c = polyscope::registerCurveNetwork("seam", ns, es);
    c->resetTransform();
    c->setRadius(0.02);
    c->setEnabled(false);

    polyscope::show();
    return 0;
}
