//
// Created by saki on 2026/07/14.
//
// demo: collapse the TmeshMut, then cut the original hmesh along the collapsed
// t-mesh (emesh_cutting) and display the resulting cut mesh with its patch
// boundaries.
#include <fstream>
#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/tutte_params.h"
#include "metriko/core/qex/sanitization.h"
#include "metriko/core/qex/gen_q_vert.h"
#include "metriko/core/qex/gen_q_port.h"
#include "metriko/core/qex/gen_q_edge.h"
#include "metriko/core/qex/gen_q_face.h"
#include "metriko/core/qex/refinement.h"
#include "common.h"
#include "visualize_hmesh.h"
#include "visualize_qex.h"
#include "visualize_tmesh_mut.h"

using namespace metriko;
static MatXd V;
static MatXi F;
static VecXc uv2;
static VecXi matching;
static VecXi singular;
static vec<bool> seam;

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);
    if (!load_cache(std::format("{}.{}.cache", argv[1], argv[2]), uv2, matching, singular, seam)) throw std::runtime_error("the cache does not exist");

    ///--- gen mport, medge ---///
    auto mg = mc::Mgrph(hm, uv2, matching, singular);
    auto tm = Tmesh(mg);
    auto X  = compute_quantization(tm, mg);
    TmeshMut tmm(mg, tm, X);

    for (int i = 0; i < 10; ++i) {
        for (const ThalfMut& th0: tmm.thalfs) {
            if (th0.id == -1) continue;
            auto& th1 = tmm.thalfs[th0.twid];
            if (th1.id == -1) continue;
            if (th0.x != 0)   continue;
            auto& tq0 = tmm.tquads[th0.tqid];
            auto& tq1 = tmm.tquads[th1.tqid];
            if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
            if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
            std::cout << "th collapse: " << th0.id << std::endl;
            tmm.collapse_thalf(th0.id);
        }
        for (const TquadMut& tq: tmm.live_tquads()) {
            Tqchain chain;
            if (tmm.collapse_tquad_chain_prepare(tq.id, chain)) {
                std::cout << "tq collapse: " << tq.id << std::endl;
                tmm.collapse_tquad_chain_execute(chain);
            }
        }
    }

    tmm.collapse_tedge_snap(false);
    tmm.collapse_tedge_snap(true);
    for (const auto& [teid, _] : tmm.live_tedges()) { tmm.collapse_tedge_snap_dedup(teid); }

    ///--- validate tedge nid chains: duplicated / backtracking nodes break the cut ---///
    for (const auto& th: tmm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids = tmm.tedges[th.teid].nids;
        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            Row3d a = get_ptloc_pos(hm, tmm.tnodes[nids[k]]);
            Row3d b = get_ptloc_pos(hm, tmm.tnodes[nids[k + 1]]);
            if (nids[k] == nids[k + 1] || (a - b).norm() < 1e-12) std::println("[warn] tedge {} (thid {}, x {}): duplicate node at {} (nid {} / {}, locs {} / {})", th.teid, th.id, th.x, k, nids[k], nids[k + 1], loc_str(tmm.tnodes[nids[k]]), loc_str(tmm.tnodes[nids[k + 1]]));
        }
        for (size_t k = 0; k + 2 < nids.size(); ++k) {
            if (nids[k] == nids[k + 2]) std::println("[warn] tedge {} (thid {}, x {}): backtrack at {} (nid {})", th.teid, th.id, th.x, k, nids[k]);
        }
    }

    ///--- cut the original mesh along the collapsed t-mesh ---///
    vec<bool> seam1;
    VecXi matching1;
    VecXi singular1;
    vec<HalfData> hdata;
    auto hm_emb = compute_embedding_cut_hmesh(hm, tmm, seam, matching, singular, seam1, matching1, singular1, hdata);
    auto hm_cut = compute_cut_mesh(*hm_emb, seam1);


    // --- validate: every halfedge must have its opposite pair. ---
    {
        std::map<std::pair<int, int>, int> cnt;
        for (int i = 0; i < hm_emb->idx.rows(); ++i)
        for (int j = 0; j < 3; ++j) {
            int a = hm_emb->idx(i, j);
            int b = hm_emb->idx(i, (j + 1) % 3);
            cnt[{a, b}]++;
        }
        int unpaired = 0, duplicated = 0;
        for (const auto& [k, c]: cnt) {
            if (c > 1)                              { if (duplicated++ < 10) std::println("[error] duplicated halfedge ({} -> {}) x{}", k.first, k.second, c); }
            if (!cnt.contains({k.second, k.first})) { if (unpaired++   < 10) std::println("[error] unpaired halfedge ({} -> {})", k.first, k.second); }
        }

        if (unpaired || duplicated) throw std::runtime_error(std::format("cut mesh validation failed: {} unpaired, {} duplicated halfedges", unpaired, duplicated));
        std::println("halfedge pairing: all {} halfedges paired", cnt.size());
    }

    visualizer::visualize_init();
    auto* base = visualizer::visualize_mesh(hm.pos, hm.idx, false, "base mesh");
    auto* embd = visualizer::visualize_mesh(hm_emb->pos, hm_emb->idx, false, "embd mesh");
    visualizer::visualize_tedge(tm, mg, uv2, &X, {}, "tedge", false);
    visualizer::visualize_seam(*hm_emb, seam1, VecXi(), "cut seam", false);
    visualizer::visualize_non_snapped_tnodes(hm, tmm, false);
    visualizer::visualize_tedge_mut_collapsed(hm, tmm, false);
    visualizer::visualize_face_collinear_error(hm, tmm, true);

    ///--- tutte parameterization (pre-SLIM initial uv) ---///
    MatXd uv;
    if (compute_tutte_parameterization(*hm_emb, tmm, seam1, hdata, uv)) {
        embd->addParameterizationQuantity("tutte uv", uv);
        igl::SLIMData sData;

        {
            MatXd uv_init(hm_cut->nV, 2);
            for (auto v: hm_cut->verts) uv_init.row(v.id) = uv.row(v.half().next().crnr().id);

            // pin the seam (boundary) vertices softly to the tutte uv
            std::vector<int>   b_;
            std::vector<Row2d> bc_;
            for (auto v: hm_cut->verts) {
                if (!v.isBoundary()) continue;
                b_.push_back(v.id);
                bc_.emplace_back(uv_init.row(v.id));
            }
            VecXi b = Eigen::Map<VecXi>(b_.data(), b_.size());
            MatXd bc(bc_.size(), 2);
            for (int i = 0; i < bc_.size(); ++i) bc.row(i) = bc_[i];

            sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;
            slim_precompute(hm_cut->pos, hm_cut->idx, uv_init, sData, sData.slim_energy, b, bc, 1e5);
            slim_solve(sData, 10);

            std::println("[slim] displacement: {}", (sData.V_o - uv_init).norm());
            auto* surf = polyscope::registerSurfaceMesh("slim result", hm_cut->pos, hm_cut->idx);
            auto* prms = surf->addVertexParameterizationQuantity("uv", sData.V_o);
            surf->setEnabled(false);
            surf->setEdgeWidth(0.7);
            prms->setEnabled(true);
            prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
            prms->setCheckerSize(1);
        }

        // ------ qex on the slim result ------
        {
            // per-corner uv from the per-vertex slim result: hm_cut and hm2
            // share the face matrix, so corner (i, j) <-> vertex hm2->idx(i, j)
            VecXc cfn(hm_emb->nF * 3);
            for (int i = 0; i < hm_cut->nF; ++i) {
            for (int j = 0; j < 3; ++j) {
                int k = hm_cut->idx(i, j);
                cfn(i * 3 + j) = complex(sData.V_o(k, 0), sData.V_o(k, 1));
            }}

            qex::sanitization(*hm_emb, matching1, singular1, 4, cfn);

            vec<qex::Qport> q_ports;
            vec<qex::Qvert> vqvs, eqvs, fqvs;
            qex::generate_q_vert(*hm_emb, cfn, vqvs, eqvs, fqvs);

            visualizer::visualize_qverts(vqvs, eqvs, fqvs, 0.001, false);

            qex::generate_vqvert_qport(*hm_emb, cfn, vqvs, q_ports);
            qex::generate_eqvert_qport(*hm_emb, cfn, eqvs, q_ports);
            qex::generate_fqvert_qport(*hm_emb, fqvs, q_ports);

            visualizer::visualize_qports(*hm_emb, cfn, q_ports, 0.001, false);
            auto qedges = qex::generate_q_edge(*hm_emb, cfn, matching1, q_ports);
            auto qfaces = qex::generate_q_faces(q_ports, qedges);

            visualizer::visualize_qedges(qedges);
            visualizer::visualize_qfaces(hm, qfaces, true);
        }
    } else std::println("[tutte] compute_tutte_parameterization failed");

    polyscope::show(); return 0;
}
