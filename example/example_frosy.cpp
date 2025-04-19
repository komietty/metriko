#include <igl/readOBJ.h>
#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include "metriko/core/vectorfield/face_rosy_field.h"

using namespace metriko;

int main(int argc, char *argv[]) {
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    MatXd V;
    MatXi F;
    igl::readOBJ("/Users/saki/dev/models/spot.obj", V, F);
    auto mesh = std::make_unique<Hmesh>(V, F);

    polyscope::init();
    int rosyN = 4;
    auto surf = polyscope::registerSurfaceMesh("mesh", mesh->pos, mesh->idx);
    surf->setSurfaceColor({0, 10./ 255., 27./ 255.});
    surf->setSmoothShade(true);

    auto rawf = std::make_unique<FaceRosyField>(*mesh, rosyN, FieldType::CurvatureAligned);
    rawf->computeMatching(MatchingType::Principal);
    auto seam = compute_seam(*rawf);
    auto cmbf = compute_combbed_field(*rawf, seam);
    MatXd rawInt(mesh->nF, 2);
    MatXd cmbInt(mesh->nF, 2);
    MatXd rawExt(mesh->nF, 3 * rosyN);
    MatXd cmbExt(mesh->nF, 3 * rosyN);
    for (Face f: mesh->faces) {
        complex rc0 = rawf->field(f.id, 0);
        complex rc1 = rawf->field(f.id, 1);
        complex rc2 = rawf->field(f.id, 2);
        complex rc3 = rawf->field(f.id, 3);
        complex cc0 = cmbf->field(f.id, 0);
        complex cc1 = cmbf->field(f.id, 1);
        complex cc2 = cmbf->field(f.id, 2);
        complex cc3 = cmbf->field(f.id, 3);
        rawInt.row(f.id) = Row2d(rc0.real(), rc0.imag()).normalized();
        cmbInt.row(f.id) = Row2d(cc0.real(), cc0.imag()).normalized();

        rawExt.block(f.id, 0, 1, 3) = (rc0.real() * f.basisX() + rc0.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 3, 1, 3) = (rc1.real() * f.basisX() + rc1.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 6, 1, 3) = (rc2.real() * f.basisX() + rc2.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 9, 1, 3) = (rc3.real() * f.basisX() + rc3.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 0, 1, 3) = (cc0.real() * f.basisX() + cc0.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 3, 1, 3) = (cc1.real() * f.basisX() + cc1.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 6, 1, 3) = (cc2.real() * f.basisX() + cc2.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 9, 1, 3) = (cc3.real() * f.basisX() + cc3.imag() * f.basisY()).normalized();
    }
    auto rawFQ = surf->addFaceVectorQuantity("raw ext", rawExt.block(0, 0, rawExt.rows(), 3));
    auto cmbFQ = surf->addFaceVectorQuantity("cmb ext", cmbExt.block(0, 0, cmbExt.rows(), 3));
    rawFQ->setEnabled(false);
    cmbFQ->setEnabled(true);
    rawFQ->setVectorLengthScale(0.004);
    cmbFQ->setVectorLengthScale(0.004);

    std::vector<glm::vec3> singPos;
    std::vector<double> singVal;
    for (Vert v: mesh->verts) {
        if (int s = rawf->singular[v.id]; s != 0) {
            Row3d p = v.pos();
            singPos.emplace_back(p.x(), p.y(), p.z());
            singVal.emplace_back(s);
        }
    }
    auto pc = polyscope::registerPointCloud("face field singulars", singPos);
    pc->addScalarQuantity("face field singular nums", singVal);
    pc->setEnabled(true);
    pc->setPointRadius(0.005);
    pc->resetTransform();

    /*--- visuailize seam ---*/
    std::vector<glm::vec3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    size_t counter = 0;
    for (auto e: mesh->edges) {
        if(seam[e.id]) {
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            nodes.emplace_back(p1.x(), p1.y(), p1.z());
            nodes.emplace_back(p2.x(), p2.y(), p2.z());
            edges.emplace_back(std::array{ counter, counter + 1 });
            counter += 2;
        }
    }
    auto pn = polyscope::registerCurveNetwork("seam", nodes, edges);
    pn->resetTransform();
    pn->setRadius(0.001);

    polyscope::show();
}
