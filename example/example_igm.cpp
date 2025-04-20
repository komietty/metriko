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
    auto rawf = std::make_unique<FaceRosyField>(*mesh, N, FieldType::CurvatureAligned);
    rawf->computeMatching(MatchingType::Curl);
    auto seam = compute_seam(*rawf);
    auto cutm = compute_cut_mesh(*mesh, seam);
    auto cmbf = compute_combbed_field(*rawf, seam);
    MatXd cmbExtRosy(mesh->nF, 3 * N);
    MatXd cmbExtZero(mesh->nF, 3);
    for (Face f: mesh->faces) {
        complex c0 = cmbf->field(f.id, 0);
        complex c1 = cmbf->field(f.id, 1);
        complex c2 = cmbf->field(f.id, 2);
        complex c3 = cmbf->field(f.id, 3);
        cmbExtZero.row(f.id) = (c0.real() * f.basisX() + c0.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 0, 1, 3) = (c0.real() * f.basisX() + c0.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 3, 1, 3) = (c1.real() * f.basisX() + c1.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 6, 1, 3) = (c2.real() * f.basisX() + c2.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 9, 1, 3) = (c3.real() * f.basisX() + c3.imag() * f.basisY()).normalized();
    }

    RosyParameterization rp(*mesh, *cutm, cmbExtRosy, cmbf->singular, cmbf->matching, seam, N);
    rp.seamless = true;
    rp.localInjectivity = false;
    //data.roundSeams = true;
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
    surf->addFaceVectorQuantity("cmb field", cmbExtZero);
    prms->setEnabled(false);
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
    auto pc = polyscope::registerPointCloud("face field singulars", singPos);
    pc->addScalarQuantity("face field singular nums", singVal);
    pc->setEnabled(true);
    pc->setPointRadius(0.002);
    pc->resetTransform();

    ///--- visuailize seam ---///
    std::vector<glm::vec3> nodesS;
    std::vector<std::array<size_t, 2>> edgesS;
    size_t counterS = 0;
    for (auto e: mesh->edges) {
        if (seam[e.id]) {
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            nodesS.emplace_back(p1.x(), p1.y(), p1.z());
            nodesS.emplace_back(p2.x(), p2.y(), p2.z());
            edgesS.emplace_back(std::array{counterS, counterS + 1});
            counterS += 2;
        }
    }
    auto curvS = polyscope::registerCurveNetwork("seam", nodesS, edgesS);
    curvS->resetTransform();
    curvS->setRadius(0.003);

    polyscope::show();
    return 0;
}
