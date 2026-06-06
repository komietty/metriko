#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "common.h"
#include "igl/upsample.h"
#include "metriko/core/tmesh/emesh_subdivide.h"
#include "metriko/core/tmesh/emesh_collapse_util.h"
#include "metriko/core/tmesh/emesh_collapse_ehalf.h"
#include "metriko/core/tmesh/emesh_collapse_equad.h"
#include "metriko/core/tmesh/emesh_postprocess.h"
#include "metriko/core/tmesh/emesh_tutte_params.h"

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
    for (Face f: hm->faces) {
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


    /// ---- visualize mesh ---- ///
    {
        auto s = visualizer::visualize_mesh_with_uv(hm->pos, hm->idx, uv1, "base_mesh", false);
        //visualizer::visualize_seam(*hm, seam, "base_seam", cmbf->matching);
        //visualizer::visualize_frosy_field(s, *hm, *rawf, *cmbf, N);
    }

    ///--- gen mport, medge ---///
    auto mg = mc::Mgrph(*hm, uv2, cmbf->matching, cmbf->singular);
    visualizer::visualize_motorcycle_graph(mg, uv2);
    //visualizer::visualize_node_adjacency(mg, uv2);

    auto tm = Tmesh(mg);
    VecXd X = compute_quantization(tm, mg);
    validate_quantization(tm, X);
    visualizer::visualize_tedge(tm, mg, uv2, &X);
    //visualizer::visualize_tedge(tm, mg, uv2);
    //visualizer::debug_tquad_sides(tm, mg, uv2);

    // todo: early return!
    polyscope::show(); return 0;


    auto sdiv_data = generate_tracked_data(*hm, uv2);
    //compute_midpoint_subdivision(sdiv_data, 1);
    //compute_barycentric_subdivision(sdiv_data, 1);
    compute_midpoint_subdivision(sdiv_data, 2);
    //compute_barycentric_subdivision(sdiv_data, 1);
    visualizer::visualize_tracked_mesh(sdiv_data, *hm, uv2, "sdiv_data");


    // =======================================================================
    // 2. 3D座標の再構築と、高解像度 Hmesh のインスタンス化
    // =======================================================================
    vec dense_pos_vec(sdiv_data.num_verts, Row3d());
    vec visited(sdiv_data.num_verts, false);

    for (size_t i = 0; i < sdiv_data.polygons.size(); ++i) {
        int parent_fid = sdiv_data.face2parent[i];
        for (size_t j = 0; j < 3; ++j) {
            int vid = sdiv_data.polygons[i][j];
            if (!visited[vid]) {
                dense_pos_vec[vid] = conversion_2d_3d(hm->faces[parent_fid], uv2, sdiv_data.uvs[i][j]);
                visited[vid] = true;
            }
        }
    }

    MatXd dV(sdiv_data.num_verts, 3);
    MatXi dF(sdiv_data.polygons.size(), 3);
    for (int i = 0; i < sdiv_data.num_verts; ++i)       { dV.row(i) = dense_pos_vec[i]; }
    for (int i = 0; i < sdiv_data.polygons.size(); ++i) { dF.row(i) << sdiv_data.polygons[i][0], sdiv_data.polygons[i][1], sdiv_data.polygons[i][2]; }
    auto dense_hm = std::make_unique<Hmesh>(dV, dF);
    auto dense_sm = compute_dense_seam(seam, *dense_hm, sdiv_data);
    //visualizer::visualize_seam(*dense_hm, dense_sm, "dense_seam", VecXi(), false);

    // 4. Emesh
    std::cout << "[Info] Constructing Emesh (Running Dijkstra)..." << std::endl;
    auto m2dv = mnode2dense_v(tm, mg, *hm, uv2, sdiv_data);
    Emesh em(tm, mg, *dense_hm, sdiv_data, m2dv, X);
    //em.collapse_equad(9);
    //em.collapse_ehalf(85);
    visualizer::visualize_mapped_mnodes(m2dv, *dense_hm);
    visualizer::visualize_eedge(em, X, "emesh_");

    // todo: early return!
    polyscope::show(); return 0;



    auto half_data = tutte::compute_half_data(em);

    int debug_tqid = 7;
    visualizer::visualize_half_data(half_data, debug_tqid, "HalfData (tqid=" + std::to_string(debug_tqid) + ")");

    MatXd dense_uv;
    tutte::compute_tutte_parameterization(*dense_hm, em, dense_sm, half_data, dense_uv);
    visualizer::visualize_mesh_with_uv(dense_hm->pos, dense_hm->idx, dense_uv, "result_mesh");

     {
        auto hm3 = compute_cut_mesh(*dense_hm, dense_sm);
        std::vector<int> b_;
        std::vector<Row2d> bc_;
        VecXi b;
        MatXd bc;

        for (auto v: hm3->verts) {
            if (v.isBoundary()) {
                Row2d val = dense_uv.row(v.half().next().crnr().id);
                b_.push_back(v.id);
                bc_.push_back(val);
            }
        }

        b = Eigen::Map<VecXi>(b_.data(), b_.size());
        bc.resize(bc_.size(), 2);


        for (int i = 0; i < bc_.size(); ++i) { bc.row(i) = bc_[i]; }
        double soft_const_p = 1e5;
        MatXd uv_init(hm3->nV, 2);
        for (auto v: hm3->verts) {
            uv_init.row(v.id) = dense_uv.row(v.half().next().crnr().id);
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
            } else if (area < 1e-12) {
                degenerate_triangles++;

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

                auto c = polyscope::registerCurveNetwork("degenerated uv face " + std::to_string(f.id), ns, es);
                c->setEnabled(true);
                c->resetTransform();
                c->setRadius(0.001);

                std::cout << "[Bad Triangle] FID: " << i << ", Area (Zero): " << area << std::endl;
            }
        }
        std::cout << "[SLIM Check] Flipped: " << flipped_triangles << ", Degenerate (Zero): " << degenerate_triangles << std::endl;

        if (true) {
            igl::SLIMData sData;
            sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;

            slim_precompute(hm3->pos, hm3->idx, uv_init, sData, sData.slim_energy, b, bc, soft_const_p);
            slim_solve(sData, 10);

            std::cout << "compute success: slim result: " << (sData.V_o - uv_init).norm() << std::endl;

            auto surf = polyscope::registerSurfaceMesh("slim result", hm3->pos, hm3->idx);
            auto prms = surf->addVertexParameterizationQuantity("uv", sData.V_o);
            surf->setEdgeWidth(0.7);
            prms->setEnabled(true);
            prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
            prms->setCheckerSize(1);
        }
    }

    polyscope::show();
    return 0;
}