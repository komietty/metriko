#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "igl/false_barycentric_subdivision.h"
#include "igl/upsample.h"
#include "metriko/core/subdivide.h"
#include "metriko/core/subdivide_with_emesh.h"
#include "metriko/core/tutte/convex_conbinatin_map.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/emesh_collapse_ehalf.h"
#include "metriko/core/tutte/emesh_collapse_equad.h"
#include "metriko/core/tutte/emesh_postprocess.h"
#include "metriko/core/tutte/tutte_params.h"
#include "common.h"

using namespace metriko;
int N = 4;
MatXd flatV;
MatXi flatF;
MatXd uv1;                            // real number uv
VecXc uv2;                            // complex number uv
std::vector<bool> seam;
std::unique_ptr<Hmesh> hm0;
std::unique_ptr<FaceRosyField> rawf;
std::unique_ptr<FaceRosyField> cmbf;

MatXd V;
MatXi F;

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    hm0  = std::make_unique<Hmesh>(V, F);
    rawf = std::make_unique<FaceRosyField>(*hm0, N, FieldType::CurvatureAligned);
    //rawf = std::make_unique<FaceRosyField>(*hm0, N, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Principal);
    auto seam = compute_seam(*rawf);
    auto cutm = compute_cut_mesh(*hm0 , seam);
    cmbf = compute_combbed_field(*rawf, seam);
    MatXd cmbExtRosy(hm0->nF, 3 * N);
    for (Face f: hm0->faces) {
        complex c0 = cmbf->field(f.id, 0);
        complex c1 = cmbf->field(f.id, 1);
        complex c2 = cmbf->field(f.id, 2);
        complex c3 = cmbf->field(f.id, 3);
        cmbExtRosy.block(f.id, 0, 1, 3) = (c0.real() * f.basisX() + c0.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 3, 1, 3) = (c1.real() * f.basisX() + c1.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 6, 1, 3) = (c2.real() * f.basisX() + c2.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 9, 1, 3) = (c3.real() * f.basisX() + c3.imag() * f.basisY()).normalized();
    }

    double t_integ0 = omp_get_wtime();
    RosyParameterization rp(*hm0, *cutm, cmbExtRosy, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();
    double t_integ1 = omp_get_wtime();
    std::cout << "[time] compute_integration: " << (t_integ1 - t_integ0) << " s" << std::endl;

    uv1.resize(hm0->nF * 3, 2);
    uv2.resize(hm0->nF * 3);
    for (Face f: hm0->faces) {
        uv1.row(f.id * 3 + 0) << rp.cfn(f.id, 0), rp.cfn(f.id, 1);
        uv1.row(f.id * 3 + 1) << rp.cfn(f.id, 4), rp.cfn(f.id, 5);
        uv1.row(f.id * 3 + 2) << rp.cfn(f.id, 8), rp.cfn(f.id, 9);
    }
    for (Face f: hm0->faces) {
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
        surf->setEdgeWidth(0.7);
        surf->setEnabled(false);
        prms->setStyle(polyscope::ParamVizStyle::GRID);
        prms->setCheckerSize(1);

        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t counter = 0;
        for (auto e: hm0->edges) {
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
        c->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.0005);
    }

    ///--- gen mport, medge ---///
    auto graph = mc::MotorcycleGraph(*hm0, uv2, cmbf->matching, cmbf->singular);
    auto tm    = metriko::tm::Tmesh(graph.mcurvs);

    VecXd R(tm.nTE);
    for (int i = 0; i < tm.nTE; i++) R[i] = tm.tedges[i].len;

    double t_quantize0 = omp_get_wtime();

    VecXd X = compute_quantization(tm, R);
    validate_quantization(tm, X);

    double t_quantize1 = omp_get_wtime();
    std::cout << "[time] compute_quantization: " << (t_quantize1 - t_quantize0) << " s" << std::endl;

    //visualizer::visualize_tedge(tm, uv2, &X, &R);


    std::vector<tutte::HalfData> half_data;
    std::set<tutte::HalfData> half_set;
    std::vector<bool> seam1;

    auto hm1 = tutte::compute_embedding_cut_hmesh(*hm0, tm, uv2, R, seam, seam1, half_set, half_data);

    ///--- visuailize seam of cut mesh ---
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t counter = 0;
        for (auto e: hm1->edges) {
            if (seam1[e.id]) {
                Row3d p1 = e.half().tail().pos();
                Row3d p2 = e.half().head().pos();
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array{counter, counter + 1});
                counter += 2;
            }
        }
        auto c = polyscope::registerCurveNetwork("seam of cut mesh", ns, es);
        c->setEnabled(false);
        c->resetTransform();
        c->setRadius(0.002);
    }

    auto em0 = std::make_unique<tutte::Emesh>(*hm1, tm, half_set, X);

    std::vector<bool> seam2;
    auto hm2 = twelve_subdivide_3(*hm1, seam1, seam2);
    auto em2 = tutte::upgrade_emesh(*em0, *hm2);

    hm1 = std::move(hm2);
    em0 = std::move(em2);
    //seam2 = seam1;

    auto success_collapse_quads = true;
    auto success_collapse_halfs = true;

    for (auto& eq: em0->equads) {
        if (!em0->collapse_quad(eq.id)) {
            success_collapse_quads = false;
            break;
        }
    }

    if (!success_collapse_quads) {
        auto surf = polyscope::registerSurfaceMesh("cut_1", hm1->pos, hm1->idx);
        surf->setEdgeWidth(0.7);
        polyscope::show();
        return 0;
    }

    for (auto& eq: em0->equads) {
        if (eq.id == -1) { continue; }

        std::vector<int> ehids_to_collapse;
        for (auto [ehid, side]: eq.edata) {
            auto& eh = em0->ehalfs[ehid];
            if (eh.x == 0) { ehids_to_collapse.push_back(ehid); }
        }

        for (int ehid : ehids_to_collapse) {
            auto& eh = em0->ehalfs[ehid];
            if (eh.eqid == -1) { continue; }
            std::cout << "----------- ehid to collapse: " << ehid << std::endl;

            if (eq.id == 15) {
                for (auto [ehid_, side_]: eq.edata)
                    std::cout << "eqid: " << eq.id << ", ehid: " << ehid_ << ", side: " << side_ << ", x: " << em0->ehalfs[ehid_].x << std::endl;
            }

            if (!em0->collapse_half(ehid)) {
                success_collapse_halfs = false;
                std::cout << "collapse failed at ehid: " << ehid << std::endl;
                goto finished_collapse_ehalfs;
            }
        }
    }

    finished_collapse_ehalfs:

    if (!success_collapse_halfs) {
        auto surf = polyscope::registerSurfaceMesh("cut_1", hm1->pos, hm1->idx);
        surf->setEdgeWidth(0.7);
        polyscope::show();
        return 0;
    }

    bool top_check = em0->check_topology();
    std::cout << "top check: " << top_check << std::endl;

    auto half_data_em = tutte::compute_half_data(*em0);

    double t_tutte0 = omp_get_wtime();
    MatXd uv;
    bool tutte_success = tutte::compute_tutte_parameterization(*hm1, *em0, seam2, half_data_em, uv);
    double t_tutte1 = omp_get_wtime();
    std::cout << "[time] compute_tutte_parameterization: " << (t_tutte1 - t_tutte0) << " s" << std::endl;

    ///--- visualize cut mesh ---///
    {
        for (auto& eq: em0->equads) { eq.debug_draw(); }
        auto surf = polyscope::registerSurfaceMesh("cut_1", hm1->pos, hm1->idx);
        auto prms = surf->addParameterizationQuantity("uv", uv);
        surf->setEdgeWidth(0.7);
        surf->setEnabled(false);
        prms->setEnabled(true);
        prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
        prms->setCheckerSize(1);
    }


    /*
    { // cut half data 1
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> val1; //
        std::vector<double> val2; // first, or crash information
        std::vector<double> val3; // tqid
        size_t count = 0;

        for (const auto& hd: half_data) {
            if (!tm.thalfs[hd.thid].cano) continue;
            //if (hd.tqid != 3) continue;
            Row3d p1 = hd.half.tail().pos();
            Row3d p2 = hd.half.head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            val1.emplace_back(hd.v0);
            val1.emplace_back(hd.v1);
            val2.emplace_back(hd.first ? 1 : 0);
            val2.emplace_back(hd.crash ? 1 : 0);
            val3.emplace_back(hd.tqid);
            es.emplace_back(std::array{count, count + 1});
            count += 2;
        }
        auto c = polyscope::registerCurveNetwork("cut half data 1", ns, es);
        auto v1 = c->addNodeScalarQuantity("val1", val1);
        auto v2 = c->addNodeScalarQuantity("val2", val2);
        c->addEdgeScalarQuantity("tqid", val3);
        c->setEnabled(false);
        v2->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.002);
    }
    { // cut half data 2
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> val1; //
        std::vector<double> val2; // first, or crash information
        std::vector<double> val3; // tqid
        size_t count = 0;

        for (const auto& hd: half_data) {
            if (tm.thalfs[hd.thid].cano) continue;
            //if (hd.tqid != 3) continue;
            Row3d p1 = hd.half.tail().pos();
            Row3d p2 = hd.half.head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            val1.emplace_back(hd.v0);
            val1.emplace_back(hd.v1);
            val2.emplace_back(hd.first ? 1 : 0);
            val2.emplace_back(hd.crash ? 1 : 0);
            val3.emplace_back(hd.tqid);
            es.emplace_back(std::array{count, count + 1});
            count += 2;
        }
        auto c = polyscope::registerCurveNetwork("cut half data 2", ns, es);
        auto v1 = c->addNodeScalarQuantity("val1", val1);
        auto v2 = c->addNodeScalarQuantity("val2", val2);
        c->addEdgeScalarQuantity("tqid", val3);
        c->setEnabled(false);
        v2->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.002);
    }
    */

    if (tutte_success) {
        auto hm3 = compute_cut_mesh(*hm1, seam2);

        std::vector<int> b_;
        std::vector<Row2d> bc_;
        VecXi b;
        MatXd bc;

        for (auto v: hm3->verts) {
            if (v.isBoundary()) {
                Row2d val = uv.row(v.half().next().crnr().id);
                b_.push_back(v.id);
                bc_.push_back(val);
            }
        }

        b = Eigen::Map<VecXi>(b_.data(), b_.size());
        bc.resize(bc_.size(), 2);


        for (int i = 0; i < (int)bc_.size(); ++i) { bc.row(i) = bc_[i]; }
        double soft_const_p = 1e5;
        MatXd uv_init(hm3->nV, 2);
        for (auto v: hm3->verts) {
            uv_init.row(v.id) = uv.row(v.half().next().crnr().id);
        }

        int flipped_triangles = 0;
        int degenerate_triangles = 0;

        for (int i = 0; i < hm3->nF; ++i) {
            auto f = hm3->faces[i];
            auto h0 = f.half();
            int v0 = h0.tail().id;
            int v1 = h0.next().tail().id;
            int v2 = h0.prev().tail().id;

            Vec2d p0 = uv_init.row(v0);
            Vec2d p1 = uv_init.row(v1);
            Vec2d p2 = uv_init.row(v2);

            // 2Dの符号付き面積（外積のZ成分）
            double area = (p1.x() - p0.x()) * (p2.y() - p0.y()) - (p1.y() - p0.y()) * (p2.x() - p0.x());

            // 判定を細かく分ける
            if (area < -1e-10) {
                flipped_triangles++;

                std::vector<glm::vec3> ns;
                std::vector<std::array<size_t, 2>> es;
                size_t counter = 0;
                for (Half h: f.adjHalfs()) {
                    Row3d p1 = h.tail().pos();
                    Row3d p2 = h.head().pos();
                    ns.emplace_back(p1.x(), p1.y(), p1.z());
                    ns.emplace_back(p2.x(), p2.y(), p2.z());
                    es.emplace_back(std::array{counter, counter + 1});
                    counter += 2;
                }
                auto c = polyscope::registerCurveNetwork("flipped uv face " + std::to_string(f.id), ns, es);
                c->setEnabled(true);
                c->resetTransform();
                c->setRadius(0.001);


                std::cout << "[Bad Triangle] FID: " << i << ", Area (Negative): " << area << std::endl;
            } else if (area <= 1e-10) {
                // 非常に薄い、または面積ゼロ
                degenerate_triangles++;
                std::cout << "[Bad Triangle] FID: " << i << ", Area (Zero): " << area << std::endl;
            }
        }
        std::cout << "[SLIM Check] Flipped: " << flipped_triangles << ", Degenerate (Zero): " << degenerate_triangles << std::endl;

        igl::SLIMData sData;
        sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;

        slim_precompute(hm3->pos, hm3->idx, uv_init, sData, sData.slim_energy, b, bc, soft_const_p);
        slim_solve(sData, 50);

        std::cout << "compute success: slim result: " << (sData.V_o - uv_init).norm() << std::endl;

        auto surf = polyscope::registerSurfaceMesh("cut_2", hm3->pos, hm3->idx);
        auto prms = surf->addVertexParameterizationQuantity("uv", sData.V_o);
        surf->setEdgeWidth(0.7);
        prms->setEnabled(true);
        prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
        prms->setCheckerSize(1);
    }

    polyscope::show();
    return 0;
}
