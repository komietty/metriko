#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "igl/upsample.h"
#include "metriko/core/subdivide.h"
#include "metriko/core/subdivide_with_emesh.h"
#include "metriko/core/tutte/convex_conbinatin_map.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/emesh_collapse_ehalf.h"
#include "metriko/core/tutte/emesh_postprocess.h"
#include "metriko/core/tutte/tutte_params.h"
#include "common.h"

using namespace metriko;
int N = 4;
MatXd flatV;
MatXi flatF;
MatXd uv1;                          // real number uv
VecXc uv2;                          // complex number uv
std::vector<bool> seam;
std::unique_ptr<Hmesh> hm;
std::unique_ptr<FaceRosyField> rawf;
std::unique_ptr<FaceRosyField> cmbf;

MatXd V;
MatXi F;

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    hm = std::make_unique<Hmesh>(V, F);
    rawf = std::make_unique<FaceRosyField>(*hm, N, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Principal);
    auto seam = compute_seam(*rawf);
    auto cutm = compute_cut_mesh(*hm  , seam);
    cmbf = compute_combbed_field(*rawf, seam);
    MatXd cmbExtRosy(hm->nF, 3 * N);
    MatXd cmbExtZero(hm->nF, 3);
    for (Face f: hm->faces) {
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

    double t_integ0 = omp_get_wtime();
    RosyParameterization rp(*hm, *cutm, cmbExtRosy, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();
    double t_integ1 = omp_get_wtime();
    std::cout << "[time] compute_integration: " << (t_integ1 - t_integ0) << " s" << std::endl;

    uv1.resize(hm->nF * 3, 2);
    uv2.resize(hm->nF * 3);
    for (const Face f: hm->faces) {
        uv1.row(f.id * 3 + 0) << rp.cfn(f.id, 0), rp.cfn(f.id, 1);
        uv1.row(f.id * 3 + 1) << rp.cfn(f.id, 4), rp.cfn(f.id, 5);
        uv1.row(f.id * 3 + 2) << rp.cfn(f.id, 8), rp.cfn(f.id, 9);
    }
    for (const Face f: hm->faces) {
        uv2(f.id * 3 + 0) = complex{uv1(f.id * 3 + 0, 0), uv1(f.id * 3 + 0, 1)};
        uv2(f.id * 3 + 1) = complex{uv1(f.id * 3 + 1, 0), uv1(f.id * 3 + 1, 1)};
        uv2(f.id * 3 + 2) = complex{uv1(f.id * 3 + 2, 0), uv1(f.id * 3 + 2, 1)};
    }

    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;


    { /// ---- visualize mesh ---- ///
        const auto surf = polyscope::registerSurfaceMesh("mesh", hm->pos, hm->idx);
        const auto prms = surf->addParameterizationQuantity("params", uv1);
        surf->setEdgeWidth(0.7);
        surf->setEnabled(false);
        prms->setStyle(polyscope::ParamVizStyle::GRID);
        prms->setCheckerSize(1);

        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> ms;
        size_t counter = 0;
        for (auto e: hm->edges) {
            if (seam[e.id]) {
                Row3d p1 = e.half().tail().pos();
                Row3d p2 = e.half().head().pos();
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                int m = cmbf->matching[e.id];
                es.emplace_back(std::array{counter, counter + 1});
                ms.emplace_back(m);
                counter += 2;
            }
        }
        auto c = polyscope::registerCurveNetwork("seam", ns, es);
        c->addEdgeScalarQuantity("matching", ms);
        //c->setEnabled(false);
        c->resetTransform();
        c->setRadius(0.001);
    }

    ///--- gen mport, medge ---///
    auto graph = mc::MotorcycleGraph(*hm, uv2, cmbf->matching, cmbf->singular);
    auto tmesh = metriko::tm::Tmesh(graph.mcurvs);

    VecXd R(tmesh.nTE);
    for (int i = 0; i < tmesh.nTE; i++) R[i] = tmesh.tedges[i].len;

    VecXd X = compute_quantization(tmesh, R);
    validate_quantization(tmesh, X);
    visualizer::visualize_tedge(tmesh, uv2, &X, &R);

    std::vector<tutte::HalfData> half_data;
    std::set<tutte::HalfData> half_set;
    std::vector<bool> seam_cut;
    auto hm1 = tutte::compute_embedding_cut_hmesh(*hm, tmesh, uv2, R, seam, seam_cut, half_set, half_data);
    auto em1 = tutte::Emesh(*hm1, tmesh, half_set, X);

    { /// ---- visualize mesh ---- ///
        const auto surf = polyscope::registerSurfaceMesh("cut mesh", hm1->pos, hm1->idx);
        surf->setEnabled(false);
        surf->setEdgeWidth(1);
        surf->setEnabled(true);
    }

    polyscope::show();
    return 0;
}