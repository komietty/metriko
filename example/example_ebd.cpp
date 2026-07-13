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
#include "metriko/core/hmesh/hpath.h"
#include "metriko/core/tmesh/tmesh_mut.h"

using namespace metriko;
int N = 4;
MatXd flatV;
MatXi flatF;
MatXd uv1;
VecXc uv2;
std::vector<bool> seam;
std::unique_ptr<Hmesh> hm;
std::unique_ptr<FaceRosyField> rawf;
std::unique_ptr<FaceRosyField> cmbf;

MatXd V;
MatXi F;

void validate_tmeshmut(const TmeshMut& tm) {
    for (const TquadMut& tq : tm.tquads) {
        if (tq.data.empty()) continue;
        double s[4] = {0, 0, 0, 0};
        for (const auto& [thid, side] : tq.data) {
            assert(thid != -1);
            auto& th_cur = tm.thalfs[thid];        if (th_cur.twid == -1 || th_cur.teid == -1 || th_cur.tqid == -1) throw std::runtime_error("th_cur");
            auto& th_twn = tm.thalfs[th_cur.twid]; if (th_twn.twid == -1 || th_twn.teid == -1 || th_twn.tqid == -1) throw std::runtime_error("th_twn");
            int x = th_cur.x;
            //if (x == 0) throw std::runtime_error("x == 0");
            s[side] += x;
        }
        if (std::abs(s[0] - s[2]) > 1e-6 || std::abs(s[1] - s[3]) > 1e-6) {
            std::println("s0: {}, s1: {}, s2: {}, s2: {}, tqid: {}", s[0], s[1], s[2], s[3], tq.id);
            throw new std::runtime_error("s[0] != s[2] || s[1] != s[3]");
        }
    }
}

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

    RosyParameterization rp(*hm, *cutm, cmbExtRosy, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();

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
    visualizer::visualize_mesh_with_uv(hm->pos, hm->idx, uv1, "base_mesh", true);
    visualizer::visualize_seam(*hm, seam, "seam", {}, false);

    ///--- gen mport, medge ---///
    auto mg = mc::Mgrph(*hm, uv2, cmbf->matching, cmbf->singular);
    auto tm = Tmesh(mg);
    VecXd X = compute_quantization(tm, mg);
    validate_quantization(tm, X);
    //assert(tm.check_non_zero_tquad(X));
    visualizer::visualize_tedge(tm, mg, uv2, &X);

    TmeshMut tmm(mg, tm, X);

    for (int i = 0; i < tmm.tquads.size(); i++) {
        int tqid = std::min(i, (int)tmm.tquads.size() - 1);
        auto allowed = tmm.allowed_range(tqid);
        std::vector<glm::vec3> rns;
        std::vector<std::array<size_t, 2>> res;
        size_t rc = 0;
        for (auto& [eid, r0, r1] : allowed) {
            Half h = hm->edges[eid].half();
            Row3d p0 = get_ptloc_pos(*hm, HmLoc(HmLocOnH{h.id, r0}));
            Row3d p1 = get_ptloc_pos(*hm, HmLoc(HmLocOnH{h.id, r1}));
            rns.emplace_back(p0.x(), p0.y(), p0.z());
            rns.emplace_back(p1.x(), p1.y(), p1.z());
            res.push_back({rc, rc + 1});
            rc += 2;
        }
        auto* reg = polyscope::registerCurveNetwork(std::format("tquad {:03} allowed region", tqid), rns, res);
        reg->setColor({0.2, 0.6, 1.0});
        reg->setRadius(0.0005);
        reg->resetTransform();
        reg->setEnabled(false);
    }

    for (int i = 0; i < 5; ++i) {
        // collapse thalf
        for (ThalfMut th0 : tmm.thalfs) {
            if (th0.id == -1) continue;
            auto& th1 = tmm.thalfs[th0.twid];
            auto& tq0 = tmm.tquads[th0.tqid];
            auto& tq1 = tmm.tquads[th1.tqid];
            if (th1.id == -1) continue;
            if (th0.x != 0)   continue;
            if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
            if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
            std::cout << "th collapse: " << th0.id << std::endl;
            tmm.collapse_thalf(th0.id);
            validate_tmeshmut(tmm);
        }

        // collapse tquad
        for (const TquadMut& tq : tmm.tquads) {
            if (tq.id == -1) continue;
            Tqchain chain;
            if (tmm.collapse_tquad_chain_prepare(tq.id, chain)) {
                std::cout << "tq collapse: " << tq.id << std::endl;
                tmm.collapse_tquad_chain_execute(chain);
                validate_tmeshmut(tmm);
            }

            //Tqaux tqaux;
            //if (tmm.collapse_tquad_prepare(tq.id, tqaux)) {
            //    std::cout << "tq collapse: " << tq.id << std::endl;
            //    tmm.collapse_tquad_execute(tq.id, tqaux);
            //    validate_tmeshmut(tmm);
            //}
        }
    }

    // collapse tquad simple chain
    /* */
    for (auto& tq_: tmm.tquads) {
        if (tq_.id == -1) continue;
        continue;
        //if (tq_.id > 69) continue;

        Tqchain chain;
        if (tmm.collapse_tquad_chain_prepare(tq_.id, chain)) {
            std::cout << "tq collapse: " << tq_.id << std::endl;
            //if (tq_.id != 69)
                tmm.collapse_tquad_chain_execute(chain);
            validate_tmeshmut(tmm);
            if (tq_.id != 69) continue;
            std::vector<glm::vec3> pcs;
            std::vector<double> vals, adjcs;
            for (const auto& p : chain.pts) {
                Row3d q = get_ptloc_pos(*hm, p.loc);
                pcs.emplace_back(q.x(), q.y(), q.z());
                vals.push_back(p.val);
                adjcs.push_back(p.adj);
                std::println("vals: {}, ahjs: {}", p.val, p.adj);
            }
            auto* pc = polyscope::registerPointCloud(std::format("chain checkpoints: {}", tq_.id), pcs);
            pc->addScalarQuantity("val", vals)->setEnabled(true);
            pc->addScalarQuantity("adj", adjcs);
            pc->setPointRadius(0.004);
        }

        // visualize a thalf sequence colored by its index (= order in the list)
        //auto viz_order = [&](const vec<int>& thids, const std::string& name) {
        //    if (thids.empty()) return;
        //    std::vector<glm::vec3> ns;
        //    std::vector<std::array<size_t, 2>> es;
        //    std::vector<double> order;
        //    size_t c = 0;
        //    for (size_t i = 0; i < thids.size(); ++i) {
        //        const TedgeMut& te = tmm.tedges[tmm.thalfs[thids[i]].teid];
        //        for (size_t k = 0; k + 1 < te.nids.size(); ++k) {
        //            Row3d a = get_ptloc_pos(*hm, tmm.tnodes[te.nids[k]]);
        //            Row3d b = get_ptloc_pos(*hm, tmm.tnodes[te.nids[k + 1]]);
        //            ns.emplace_back(a.x(), a.y(), a.z());
        //            ns.emplace_back(b.x(), b.y(), b.z());
        //            es.push_back({c, c + 1}); c += 2;
        //            order.push_back((double)i);
        //        }
        //    }
        //    auto* cn = polyscope::registerCurveNetwork(name, ns, es);
        //    cn->addEdgeScalarQuantity(std::format("order: {}", tq_.id), order)->setEnabled(true);
        //    cn->setRadius(0.002);
        //};
        //viz_order(chain.thids_z, "chain thids_z");
        //viz_order(chain.thids_t, "chain thids_t");
        //viz_order(chain.thids_b, "chain thids_b");
    }
    //for (auto tqid_ : {103}) {
    //    Tqchain chain;
    //    if (tmm.collapse_tquad_chain_prepare(tqid_, chain)) {
    //        tmm.collapse_tquad_chain_execute(chain);
    //        if (tqid_ == 103) {
    //            std::vector<glm::vec3> pcs;
    //            std::vector<double> vals, adjcs;
    //            for (const auto& p : chain.pts) {
    //                Row3d q = get_ptloc_pos(*hm, p.loc);
    //                pcs.emplace_back(q.x(), q.y(), q.z());
    //                vals.push_back(p.val);
    //                adjcs.push_back(p.adj);
    //                std::println("vals: {}, ahjs: {}", p.val, p.adj);
    //            }
    //            auto* pc = polyscope::registerPointCloud(std::format("chain checkpoints: {}", tqid_), pcs);
    //            pc->addScalarQuantity("val", vals)->setEnabled(true);
    //            pc->addScalarQuantity("adj", adjcs);
    //            pc->setPointRadius(0.004);
    //        }
    //    }
    //}

    /*
    // second loop
    for (int i = 0; i < 20; ++i) {
        // collapse thalf
        for (ThalfMut th0 : tmm.thalfs) {
            auto& th1 = tmm.thalfs[th0.twid];
            auto& tq0 = tmm.tquads[th0.tqid];
            auto& tq1 = tmm.tquads[th1.tqid];
            if (th0.id == -1) continue;
            if (th1.id == -1) continue;
            if (th0.x != 0)   continue;
            if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
            if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
            std::cout << "th collapse: " << th0.id << std::endl;
            tmm.collapse_thalf(th0.id);
            validate_tmeshmut(tmm);
        }

        for (const TquadMut& tq : tmm.tquads) {
            if (tq.id == -1) continue;
            std::cout << "tq collapse: " << tq.id << std::endl;
            Tqaux tqaux;
            if (tmm.collapse_tquad_prepare(tq.id, tqaux)) {
                tmm.collapse_tquad_execute(tq.id, tqaux);
                validate_tmeshmut(tmm);
            }
        }
    }
    */

    // debug view
    for (const TquadMut& tq : tmm.tquads) {
        if (tq.id == -1) continue;
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> eside, ex, ey, er, ethid;   // per-edge (thalf) params
        size_t c = 0;
        for (const TdataMut& d : tq.data) {
            const ThalfMut& th = tmm.thalfs[d.thid];
            const TedgeMut& te = tmm.tedges[th.teid];
            for (size_t i = 0; i + 1 < te.nids.size(); ++i) {
                Row3d a = get_ptloc_pos(*hm, tmm.tnodes[te.nids[i]]);
                Row3d b = get_ptloc_pos(*hm, tmm.tnodes[te.nids[i + 1]]);
                ns.emplace_back(a.x(), a.y(), a.z());
                ns.emplace_back(b.x(), b.y(), b.z());
                es.push_back({c, c + 1}); c += 2;
                eside.push_back(d.side);
                ex.push_back(th.x);
                ey.push_back(th.x > 0 ? 1 : 0);
                er.push_back(th.r);
                ethid.push_back(d.thid);
            }
        }
        if (ns.empty()) continue;
        auto* cn = polyscope::registerCurveNetwork(std::format("tq {:03}", tq.id), ns, es);
        cn->addEdgeScalarQuantity("side", eside);
        auto cx = cn->addEdgeScalarQuantity("x", ex);
        auto cy = cn->addEdgeScalarQuantity("y", ey);
        cy->setEnabled(true);
        cn->addEdgeScalarQuantity("r", er);
        cn->addEdgeScalarQuantity("thid", ethid);
        cn->setMaterial("flat");
        cn->setRadius(0.0007); cn->resetTransform();
    }

    polyscope::show(); return 0; // early return!

    /*
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
    */
}