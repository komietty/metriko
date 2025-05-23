#include <igl/readOBJ.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/misc/visualizer/tmesh/visualize_qgp_edge.h"
#include "metriko/misc/visualizer/tmesh/visualize_qgp_mesh.h"
#include "metriko/misc/visualizer/tmesh/visualize_tedge.h"

using namespace metriko;
int main(int argc, char **argv) {
    MatXd V;
    MatXi F;
    int N = 4;
    igl::readOBJ(argv[1], V, F);
    auto mesh = std::make_unique<Hmesh>(V, F);
    auto rawf = std::make_unique<FaceRosyField>(*mesh, N, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Principal);
    auto seam = compute_seam(*rawf);
    auto cutm = compute_cut_mesh(*mesh, seam);
    auto cmbf = compute_combbed_field(*rawf, seam);
    MatXd extf = compute_extrinsic_field(*mesh, *cmbf, N);

    RosyParameterization rp(*mesh, *cutm, extf, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();
    MatXd uv1(mesh->nF * 3, 2);
    VecXc uv2(mesh->nF * 3);

    for (const Face f: mesh->faces) {
        uv1.row(f.id * 3 + 0) << rp.cfn(f.id, 0), rp.cfn(f.id, 1);
        uv1.row(f.id * 3 + 1) << rp.cfn(f.id, 4), rp.cfn(f.id, 5);
        uv1.row(f.id * 3 + 2) << rp.cfn(f.id, 8), rp.cfn(f.id, 9);
    }
    for (const Face f: mesh->faces) {
        uv2(f.id * 3 + 0) = complex{uv1(f.id * 3 + 0, 0), uv1(f.id * 3 + 0, 1)};
        uv2(f.id * 3 + 1) = complex{uv1(f.id * 3 + 1, 0), uv1(f.id * 3 + 1, 1)};
        uv2(f.id * 3 + 2) = complex{uv1(f.id * 3 + 2, 0), uv1(f.id * 3 + 2, 1)};
    }

    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

    /// ---- visualize mesh ---- ///
    {
        const auto surf = polyscope::registerSurfaceMesh("mesh", V, F);
        const auto prms = surf->addParameterizationQuantity("params", uv1);
        surf->setEnabled(false);
        surf->addFaceVectorQuantity("cmb field", extf);
        prms->setStyle(polyscope::ParamVizStyle::GRID);
        prms->setCheckerSize(1);
    }

    ///--- visuailize seam ---///
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2> > es;
        size_t counter = 0;
        for (auto e: mesh->edges) {
            if (seam[e.id]) {
                Row3d p1 = e.half().tail().pos();
                Row3d p2 = e.half().head().pos();
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array{counter, counter + 1});
                counter += 2;
            }
        }
        auto c = polyscope::registerCurveNetwork("seam", ns, es);
        c->setEnabled(false);
        c->resetTransform();
        c->setRadius(0.003);
    }

    ///--- gen mport, medge ---///
    auto graph = MotorcycleGraph(*mesh, uv2, cmbf->matching, cmbf->singular);
    auto tmesh = Tmesh(graph.mcurvs);
    VecXd R = VecXd::Zero(tmesh.nTE);
    for (int i = 0; i < tmesh.nTE; i++) {
        bool bgn = false;
        const auto &te = tmesh.tedges[i];
        for (const Msgmt &seg: te.seg_fr.curv->sgmts) {
            if (seg == te.seg_fr) bgn = true;
            if (bgn) {
                R[i] += std::abs(seg.diff());
                if (seg == te.seg_to) break;
            }
        }
    }

    VecXd X = compute_quantization(tmesh, R);
    validate_quantization(tmesh, X);
    visualizer::visualize_tedge(tmesh, uv2, &X, &R);

    std::vector<std::vector<visualizer::SplitVert> > split_verts(tmesh.nTE);
    for (int i = 0; i < tmesh.nTE; i++)
        split_verts[i] = visualizer::construct_verts_on_tedge(tmesh, X, R, i);

    MatXd VQ;
    MatXi FQ;
    std::vector<glm::vec3> CQ;

    for (int ii = 0; ii < tmesh.tquads.size(); ii++) {
        std::vector<visualizer::SplitArc> arcs1(0);
        std::vector<visualizer::SplitArc> arcs2(0);
        const Tquad &tq_draw = tmesh.tquads[ii];
        std::vector<std::vector<visualizer::SplitElem> > sv_quad;
        gen_split_table(tmesh, tq_draw, split_verts, sv_quad);
        visualizer::visualize_split_edge(tmesh, R, uv2, cmbf->matching, sv_quad, ii, arcs1, arcs2);
        visualizer::visualize_split_mesh(arcs1, arcs2, tmesh, uv2, ii, VQ, FQ, CQ);
    } {
        auto surf = polyscope::registerSurfaceMesh("quad mesh", VQ, FQ);
        surf->addFaceColorQuantity("color", CQ)->setEnabled(true);
        surf->setShadeStyle(polyscope::MeshShadeStyle::Flat);
        surf->setEdgeWidth(1.);
    }
    polyscope::show();
    return 0;
}
