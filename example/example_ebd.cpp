#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tutte/embedding.h"
#include "metriko/misc/visualizer/tmesh/visualize_tedge.h"
#include "../include/metriko/core/tutte/visualize_tedge_ebd.h"
#include "igl/upsample.h"
#include "metriko/core/tutte/convex_conbinatin_map.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/tutte_params.h"

// 1: computing speed up
// 2: tmesh collapsing
// 3: update vector field

using namespace metriko;
int N = 4;
MatXd flatV;
MatXi flatF;
MatXd uv1;                          // real number uv
VecXc uv2;                          // complex number uv
std::vector<bool> seam;
std::unique_ptr<Hmesh> mesh;
std::unique_ptr<FaceRosyField> rawf;
std::unique_ptr<FaceRosyField> cmbf;

MatXd V;
MatXi F;

int main(int argc, char** argv) {

    {
        std::cout << "OMP num threads = " << omp_get_max_threads() << std::endl;
        const int N_bench = 10'000'000;
        std::vector<double> tmp(N_bench);

        // シリアル版
        double t0 = omp_get_wtime();
        for (int i = 0; i < N_bench; ++i) {
            tmp[i] = i;
        }
        double t1 = omp_get_wtime();

        // 並列版
        double t2 = omp_get_wtime();
        #pragma omp parallel for
        for (int i = 0; i < N_bench; ++i) {
            tmp[i] = i;
        }
        double t3 = omp_get_wtime();

        std::cout << "[bench] serial:  " << (t1 - t0) << " s\n";
        std::cout << "[bench] parallel:" << (t3 - t2) << " s\n";
        std::cout << "[bench] speedup:  " << (t1 - t0) / (t3 - t2) << " x\n";
    }

    std::cout << "Available :SIMD Instructions: "<< Eigen::SimdInstructionSetsInUse() << std::endl;
    igl::readOBJ(argv[1], V, F);

    //igl::upsample(V, F, 1);

    mesh = std::make_unique<Hmesh>(V, F);
    rawf = std::make_unique<FaceRosyField>(*mesh, N, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Principal);
    auto seam = compute_seam(*rawf);
    auto cutm = compute_cut_mesh(*mesh, seam);
    cmbf = compute_combbed_field(*rawf, seam);
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

    RosyParameterization rp(*mesh, *cutm, cmbExtRosy, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();

    uv1.resize(mesh->nF * 3, 2);
    uv2.resize(mesh->nF * 3);
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
    const auto surf = polyscope::registerSurfaceMesh("mesh", V, F);
    //const auto surf = polyscope::registerSurfaceMesh("cut_mesh", cutm->pos, cutm->idx);
    const auto prms = surf->addParameterizationQuantity("params", uv1);
    surf->setEdgeWidth(1);
    surf->setEnabled(false);
    //surf->addFaceVectorQuantity("cmb field", cmbExtZero);
    prms->setStyle(polyscope::ParamVizStyle::GRID);
    prms->setCheckerSize(1);

    ///--- visuailize seam ---///
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
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
        c->setRadius(0.001);
    }

    ///--- gen mport, medge ---///
    auto graph = MotorcycleGraph(*mesh, uv2, cmbf->matching, cmbf->singular);
    auto tmesh = Tmesh(graph.mcurvs);
    VecXd R = VecXd::Zero(tmesh.nTE);
    for (int i = 0; i < tmesh.nTE; i++) {
        bool bgn = false;
        const auto& te = tmesh.tedges[i];
        for (const Msgmt& seg: te.seg_fr.curv->sgmts) {
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




    std::set<tutte::HalfData> half_data;
    std::vector<bool> seam_cut;

    double t_cut0 = omp_get_wtime();
    Hmesh hm_cut = tutte::compute_embedding_cut_hmesh(*mesh, tmesh, uv2, R, seam, seam_cut, half_data);
    double t_cut1 = omp_get_wtime();
    std::cout << "[time] compute_embedding_cut_hmesh: " << (t_cut1 - t_cut0) << " s" << std::endl;

    double t_tutte0 = omp_get_wtime();
    MatXd uv = tutte::compute_tutte_parameterization(hm_cut, tmesh, seam_cut, half_data, X);
    double t_tutte1 = omp_get_wtime();
    std::cout << "[time] compute_tutte_parameterization: " << (t_tutte1 - t_tutte0) << " s" << std::endl;

    ///--- visualize cut mesh ---///
    const auto surf_cut = polyscope::registerSurfaceMesh("cut_1", hm_cut.pos, hm_cut.idx);
    surf_cut->setEdgeWidth(1);

    {
        auto prms1 = surf_cut->addParameterizationQuantity("params_1", uv);
        prms1->setEnabled(true);
        prms1->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
        prms1->setCheckerSize(1);

    }

    /*
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> val;
        size_t count = 0;

        for (const auto& hd: half_data) {
            if (hd.tqid != 0) continue;
            Row3d p1 = hd.half.tail().pos();
            Row3d p2 = hd.half.head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            val.emplace_back(hd.v0);
            val.emplace_back(hd.v1);
            es.emplace_back(std::array{count, count + 1});
            count += 2;
        }
        auto c = polyscope::registerCurveNetwork("cut half data", ns, es);
        auto v = c->addNodeScalarQuantity("val", val);
        c->setEnabled(true);
        v->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.002);
    }
    */

    ///--- visuailize seam of cut mesh ---///
    ///auto hm_cut_cut = compute_cut_mesh(hm_cut, seam_cut);
    /*
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t counter = 0;
        for (auto e: hm_cut.edges) {
            if (seam_cut[e.id]) {
                Row3d p1 = e.half().tail().pos();
                Row3d p2 = e.half().head().pos();
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array{counter, counter + 1});
                counter += 2;
            }
        }
        auto c = polyscope::registerCurveNetwork("seam_cut", ns, es);
        c->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.001);
    }
    */

    /*
    std::vector<std::vector<visualizer::SplitVert> > split_verts(tmesh.nTE);
    for (int i = 0; i < tmesh.nTE; i++)
        split_verts[i] = visualizer::construct_verts_on_tedge(tmesh, X, R, i);


    std::vector<EmbeddedTEdge> etes;
    std::vector<EmbeddedTHalf> eths;
    std::vector passthrough(mesh->nE, false);
    for (int i = 0; i < tmesh.nTE; i++) {
        auto res = gen_embedded_tedge_easy(*mesh, uv2, tmesh.tedges[i], passthrough);
        if (res.has_value()) {
            etes.emplace_back(res.value());
            eths.emplace_back(res.value(), true);
            eths.emplace_back(res.value(), false);
        } else throw std::runtime_error("failed to generate embedded tedge");
    }

    reassign_quantization_values(*mesh, X, etes);
    visualizer::visualize_embedding(*mesh, etes, X);

    MatXd uv = compute_tutte_parameterization(*mesh, tmesh, etes, eths, seam, X);
    auto prms1 = surf->addParameterizationQuantity("params_1", uv);
    prms1->setEnabled(true);
    prms1->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
    prms1->setCheckerSize(1);
    */

    /*
    std::vector<int> b_;
    std::vector<Row2d> bc_;

    for (auto v: hm_cut_cut->verts) {
        if (v.isBoundary()) {
            Row2d val = uv.row(v.half().next().crnr().id);
            b_.push_back(v.id);
            bc_.push_back(val);
        }
    }

    Eigen::VectorXi b;
    Eigen::MatrixXd bc;

    b = Eigen::Map<VecXi>(b_.data(), b_.size());
    bc.resize(bc_.size(), 2);
    for (int i = 0; i < static_cast<int>(bc_.size()); ++i) { bc.row(i) = bc_[i]; }
    double soft_const_p = 1e10;
    Eigen::MatrixXd uv_init(hm_cut_cut->nV, 2);
    for (auto v: hm_cut_cut->verts) {
        uv_init.row(v.id) = uv.row(v.half().next().crnr().id);
    }
    igl::SLIMData sData;
    sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;

    slim_precompute(hm_cut_cut->pos, hm_cut_cut->idx, uv_init, sData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, soft_const_p);
    slim_solve(sData, 50);

    const auto surf_cut_cut = polyscope::registerSurfaceMesh("cut_2", hm_cut_cut->pos, hm_cut_cut->idx);
    auto prms2 = surf_cut_cut->addVertexParameterizationQuantity("params_2", sData.V_o);
    prms2->setEnabled(true);
    prms2->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
    prms2->setCheckerSize(1);
    */

    polyscope::show();
    return 0;
}
