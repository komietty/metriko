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
#include "metriko/core/subdivide.h"
#include "metriko/core/tmesh/emesh_collapse_util.h"
#include "metriko/core/tmesh/emesh_collapse_ehalf.h"
#include "metriko/core/tmesh/emesh_collapse_equad.h"
#include "metriko/core/tmesh/emesh_subdivide.h"

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
        //visualizer::visualize_frosy_field(surf, *hm, *rawf, *cmbf, N);
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
        c->setEnabled(false);
        c->resetTransform();
        c->setRadius(0.001);
    }
    //{
    //    std::vector<bool> seam2;
    //    std::vector<bool> seam3;
    //    std::vector<bool> seam4;
    //    auto hm2 = four_subdivide(*hm, seam, seam2);
    //    auto hm3 = four_subdivide(*hm2, seam2, seam3);
    //    auto hm4 = four_subdivide(*hm3, seam3, seam4);
    //    const auto surf = polyscope::registerSurfaceMesh("mesh_subdiv", hm4->pos, hm4->idx);
    //    surf->setEdgeWidth(0.7);
    //}

    ///--- gen mport, medge ---///
    auto mg = mc::Mgrph(*hm, uv2, cmbf->matching, cmbf->singular);
    //visualizer::visualize_motorcycle_graph(mg, uv2);
    //visualizer::visualize_node_adjacency(mg, uv2);
    auto tm = Tmesh(mg);
    VecXd X = compute_quantization(tm, mg);
    //validate_quantization(tmesh, X);
    visualizer::visualize_tedge(tm, mg, uv2, &X);
    //visualizer::visualize_tedge(tm, mg, uv2);
    //visualizer::debug_tquad_sides(tm, mg, uv2);

    auto sdiv_data = compute_midpoint_subdivision(*hm, uv2, 2);
    visualizer::visualize_tracked_mesh(sdiv_data, *hm, uv2, "sdiv_data (Level 2)");



    // =======================================================================
    // 2. 3D座標の再構築と、高解像度 Hmesh のインスタンス化
    // =======================================================================
    std::vector<Row3d> dense_pos_vec(sdiv_data.num_verts);
    std::vector<bool> visited(sdiv_data.num_verts, false);

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

    MatXd dense_V(sdiv_data.num_verts, 3);
    MatXi dense_F(sdiv_data.polygons.size(), 3);

    for (int i = 0; i < sdiv_data.num_verts; ++i)       { dense_V.row(i) = dense_pos_vec[i]; }
    for (int i = 0; i < sdiv_data.polygons.size(); ++i) { dense_F.row(i) << sdiv_data.polygons[i][0], sdiv_data.polygons[i][1], sdiv_data.polygons[i][2]; }
    auto dense_hm = std::make_unique<Hmesh>(dense_V, dense_F);

    // =======================================================================
    // 3. Mnode -> dense_vid のマッピング辞書の自動構築 (安全版)
    // =======================================================================
    std::map<int, int> mnode2dense_v;
    for (const auto& te : tm.tedges) {
        auto map_node = [&](int nid, int fid) {
            if (mnode2dense_v.contains(nid)) return;
            complex target_uv = mc::get_face_uv(mg.mnodes[nid], fid, *hm, uv2);

            int best_vid = -1;
            double min_dist = 1e9;
            for (size_t i = 0; i < sdiv_data.polygons.size(); ++i) {
                if (sdiv_data.face2parent[i] != fid) continue;
                for (int j = 0; j < 3; ++j) {
                    double d = std::abs(sdiv_data.uvs[i][j] - target_uv);
                    if (d < min_dist) { min_dist = d; best_vid = sdiv_data.polygons[i][j]; }
                }
            }
            assert(best_vid != -1 && "Corresponding dense vertex not found!");
            mnode2dense_v[nid] = best_vid;
        };
        map_node(te.fr_nid, te.segs.front().face_id);
        map_node(te.to_nid, te.segs.back().face_id);
    }

    // =======================================================================
    // 4. Emesh の構築と可視化
    // =======================================================================
    std::cout << "[Info] Constructing Emesh (Running Dijkstra)..." << std::endl;
    Emesh em(tm, mg, *dense_hm, sdiv_data, mnode2dense_v, X);
    verts_inside(em.equads[9], em, *dense_hm);
    //em.collapse_ehalf(85);
    em.collapse_equad(9);
    //verts_inside(em.equads[0], em, *dense_hm);
    //verts_inside(em.equads[16], em, *dense_hm);

    //visualizer::visualize_eedge(em, X, "emesh_");
    for (const auto& eq : em.equads) {
        if (eq.id == -1) continue;
        visualizer::visualize_equad(em, eq);
    }


    polyscope::show();
    return 0;
}