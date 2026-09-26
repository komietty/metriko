//
// Created by saki on 2026/07/14.
//
// demo: collapse the Emesh, then cut the original hmesh along the collapsed
// t-mesh (emesh_cutting) and display the resulting cut mesh with its patch
// boundaries.
#include <fstream>
#include <limits>
#include <chrono>
#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/emesh.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/tutte_params.h"
#include "metriko/core/qex/sanitization.h"
#include "metriko/core/qex/gen_q_vert.h"
#include "metriko/core/qex/gen_q_port.h"
#include "metriko/core/qex/gen_q_edge.h"
#include "metriko/core/qex/gen_q_face.h"
#include "common_io.h"
#include "common_visualizer.h"
#include "cleanup.h"

using namespace metriko;
static MatXd V;
static MatXi F;
static VecXc uv2;
static VecXi matching;
static VecXi singular;
static vec<bool> seam;

int main(int argc, char** argv) {
    const auto t_start = std::chrono::steady_clock::now();
    auto t_lap = t_start;
    auto lap = [&](const char* name) { auto t = std::chrono::steady_clock::now(); std::println("[time] {:<22} {:8.3f} s", name, std::chrono::duration<double>(t - t_lap).count()); t_lap = t; };
    igl::readOBJ(argv[1], V, F);
    //cleanup::decimate_and_clean(V, F,  100000);
    Hmesh hm(V, F);
    if (!load_cache(std::format("{}.{}.cache", argv[1], argv[2]), uv2, matching, singular, seam)) throw std::runtime_error("the cache does not exist");

    ///--- load the collapsed t-mesh produced by example_1 ---///
    Emesh em(hm);
    if (!load_emesh(std::format("{}.{}.em", argv[1], argv[2]), em))
        throw std::runtime_error("the em cache does not exist. run example_1 first");
    lap("load");

    ///--- validate tedge nid chains: duplicated / backtracking nodes break the cut ---///
    for (const auto& th: em.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids = em.tedges[th.teid].nids;
        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            Row3d a = get_ptloc_pos(hm, em.tnodes[nids[k]]);
            Row3d b = get_ptloc_pos(hm, em.tnodes[nids[k + 1]]);
            if (nids[k] == nids[k + 1] || (a - b).norm() < 1e-12) std::println("[warn] tedge {} (thid {}, x {}): duplicate node at {} (nid {} / {}, locs {} / {})", th.teid, th.id, th.x, k, nids[k], nids[k + 1], loc_str(em.tnodes[nids[k]]), loc_str(em.tnodes[nids[k + 1]]));
        }
        for (size_t k = 0; k + 2 < nids.size(); ++k) {
            if (nids[k] == nids[k + 2]) std::println("[warn] tedge {} (thid {}, x {}): backtrack at {} (nid {})", th.teid, th.id, th.x, k, nids[k]);
        }
    }

    { // cross-tedge duplicate nodes: two distinct tnodes at the same location
        std::map<std::string, std::pair<int, int>> seen; // loc string -> (nid, teid)
        for (const auto& th: em.thalfs) {
            if (th.id == -1 || !th.cano) continue;
            for (int nid: em.tedges[th.teid].nids) {
                auto key = loc_str(em.tnodes[nid]);
                if (auto [it, ins] = seen.try_emplace(key, nid, th.teid); !ins && it->second.first != nid) std::println("[warn] duplicate tnode: {} (nid {} in teid {} / nid {} in teid {})", key, it->second.first, it->second.second, nid, th.teid);
            }
        }
    }

    { // boundary of every live tquad must be a closed loop in order
        for (const auto& tq: em.live_tquads()) {
            auto& d = tq.data;
            for (size_t k = 0; k < d.size(); ++k) {
                const auto& a = em.thalfs[d[k].thid];
                const auto& b = em.thalfs[d[(k + 1) % d.size()].thid];
                if (a.loc_to() != b.loc_fr()) std::println("[warn] tquad {}: boundary broken between thid {} and thid {} ({} vs {})", tq.id, a.id, b.id, loc_str(a.loc_to()), loc_str(b.loc_fr()));
            }
        }
    }

    ///--- cut the original mesh along the collapsed t-mesh ---///
    vec<bool> seam1;
    VecXi matching1;
    VecXi singular1;
    vec<HalfData> hdata;
    auto hm_emb = compute_embedding_cut_hmesh(hm, em, seam, matching, singular, seam1, matching1, singular1, hdata);
    auto hm_cut = compute_cut_mesh(*hm_emb, seam1);
    lap("cut");


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
    visualizer::visualize_seam(*hm_emb, seam1, VecXi(), "cut seam", false);
    visualizer::visualize_unsnapped_tnodes(hm, em, false);
    visualizer::visualize_tedges(hm, em, "tedges collapsed", false);

    ///--- tutte parameterization (pre-SLIM initial uv) ---///
    MatXd uv =compute_tutte_parameterization(*hm_emb, em, seam1, hdata);
    lap("tutte");
    {
        embd->addParameterizationQuantity("tutte uv", uv);
        igl::SLIMData sData;

        auto uv_at = [&](Crnr c) { return complex(uv(c.id, 0), uv(c.id, 1)); };
        vec<double> flag(hm_emb->nF, 0.);    // 0 ok, 1 inverted, 2 exactly degenerate
        vec<double> sarea(hm_emb->nF, 0.);   // signed uv area: how close a face is to folding
        int bad = 0;
        for (Face f: hm_emb->faces) {
            auto [a, b, c] = f.crnrs();
            double o = orientation(uv_at(a), uv_at(b), uv_at(c));
            sarea[f.id] = o / 2.;
            if (o > 0) continue;
            flag[f.id] = o < 0 ? 1. : 2.;
            bad++;
        }
        std::println("[tutte] {} degenerate/inverted faces", bad);

        // where they sit on the surface
        embd->addFaceScalarQuantity("tutte flip (1:inv 2:degen)", flag)->setEnabled(true);
        embd->addFaceScalarQuantity("tutte signed uv area", sarea);

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
            constexpr int    slim_max_iter = 50;
            constexpr double slim_rel_tol  = 1e-4;
            double prev = std::numeric_limits<double>::infinity();
            for (int i = 0; i < slim_max_iter; ++i) {
                slim_solve(sData, 1);
                if (std::abs(prev - sData.energy) < slim_rel_tol * std::abs(sData.energy)) break;
                std::println("[slim] iter {} energy {}", i, sData.energy);
                prev = sData.energy;
            }

            std::println("[slim] displacement: {}", (sData.V_o - uv_init).norm());
            lap("slim");
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
            lap("qex sanitization");

            vec<qex::Qport> q_ports;
            vec<qex::Qvert> vqvs, eqvs, fqvs;
            qex::generate_q_vert(*hm_emb, cfn, vqvs, eqvs, fqvs);
            lap("qex q_vert");

            //visualizer::visualize_qverts(vqvs, eqvs, fqvs, 0.001, false);

            qex::generate_vqvert_qport(*hm_emb, cfn, vqvs, q_ports);
            qex::generate_eqvert_qport(*hm_emb, cfn, eqvs, q_ports);
            qex::generate_fqvert_qport(*hm_emb, fqvs, q_ports);
            lap("qex q_port");

            //visualizer::visualize_qports(*hm_emb, cfn, q_ports, 0.001, false);

            auto qedges = qex::generate_q_edge(*hm_emb, cfn, matching1, q_ports);
            lap("qex q_edge");
            auto qfaces = qex::generate_q_faces(q_ports, qedges);
            lap("qex q_face");
            std::println("[time] {:<22} {:8.3f} s", "total (before visualize)", std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count());

            //visualizer::visualize_qedges(qedges);
            auto [qv, qidx] = extract_quad_mesh(hm, qfaces, true);
            lap("quad mesh refinement");
            visualizer::visualize_quad_patch(qv, qidx, label_quad_patches(em, singular, qfaces, qidx));
            lap("quad patch");
        }
    }

    polyscope::show(); return 0;
}
