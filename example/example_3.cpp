//
// example_3: the whole pipeline in one run, for checking a mesh end to end.
//   stage 0 (example_0): vectorfield + rosy parameterization
//   stage 1 (example_1): quantization + t-mesh collapse / snap
//   stage 2 (example_2): cutting + tutte + slim + qex
// usage: example_3 <mesh.obj> <scale>
#include <limits>
#include <chrono>
#include <igl/readOBJ.h>
#include <igl/slim.h>

#include "cleanup.h"
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/emesh.h"
#include "metriko/core/tmesh/emesh_validate.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/tutte_params.h"
#include "metriko/core/qex/sanitization.h"
#include "metriko/core/qex/gen_q_vert.h"
#include "metriko/core/qex/gen_q_port.h"
#include "metriko/core/qex/gen_q_edge.h"
#include "metriko/core/qex/gen_q_face.h"
#include "common_io.h"
#include "common_visualizer.h"

using namespace metriko;
static constexpr int N = 4;

int main(int argc, char** argv) {
    const auto t_start = std::chrono::steady_clock::now();
    auto t_lap = t_start;
    auto print_time = [](const char* name, auto dur) { std::println("[time] {:<26} {:8.3f} s", name, std::chrono::duration<double>(dur).count()); };
    auto lap = [&](const char* name) { auto t = std::chrono::steady_clock::now(); print_time(name, t - t_lap); t_lap = t; };

    MatXd V;
    MatXi F;
    igl::readOBJ(argv[1], V, F);
    cleanup::cleanup_mesh(V, F);
    Hmesh hm(V, F);
    lap("load");

    ///--- stage 0: vectorfield + rosy parameterization ---///
    FaceRosyField rawf(hm, N, FieldType::Smoothest);
    rawf.computeMatching(MatchingType::Principal);
    auto seam = compute_seam(rawf);
    auto cutm = compute_cut_mesh(hm, seam);
    auto cmbf = compute_combbed_field(rawf, seam);
    auto extf = compute_extrinsic_field(*cmbf, N);
    lap("field");

    RosyParameterization rp(hm, *cutm, extf, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();
    lap("parameterization");

    VecXc uv2(hm.nF * 3);
    for (const Face f: hm.faces) {
        uv2(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        uv2(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        uv2(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }
    const VecXi& matching = cmbf->matching;
    const VecXi& singular = cmbf->singular;
    save_cache(std::format("{}.{}.cache", argv[1], argv[2]), uv2, matching, singular, seam);

    ///--- stage 1: quantization + t-mesh collapse / snap ---///
    auto mg = Mgrph(hm, uv2, matching, singular);
    auto tm = Tmesh(mg);
    auto X  = compute_quantization(tm, mg);
    Emesh em(mg, tm, X);
    lap("motorcycle + quantization");

    for (int i = 0; i < 100; ++i) {
        for (Ehalf th0 : em.thalfs) {
            if (th0.id == -1 || th0.x != 0) continue;
            auto& th1 = em.thalfs[th0.twid];
            auto& tq0 = em.tquads[th0.tqid];
            auto& tq1 = em.tquads[th1.tqid];
            if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
            if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
            em.collapse_thalf(th0.id);
        }
        for (const auto& [tqid, _] : em.live_tquads())
            if (Tqchain c; em.collapse_tquad_chain_prepare(tqid, c)) em.collapse_tquad_chain_execute(c);
        if (rg::none_of(em.thalfs, [](const Ehalf& th) { return th.id != -1 && th.x == 0; })) break;
    }
    lap("collapse");

    em.collapse_tedge_snap(false);
    em.collapse_tedge_snap(true);
    for (const auto& [teid, _] : em.live_tedges()) em.collapse_tedge_snap_dedup(teid);
    repair_crossing_tedges(em);
    if (validate_no_crossing(em, "after repair") > 0) throw std::runtime_error("tedge contacts remain");
    save_emesh(std::format("{}.{}.em", argv[1], argv[2]), em);
    lap("snap + repair");

    ///--- stage 2: cut, tutte, slim, qex ---///
    vec<bool> seam1;
    VecXi matching1, singular1;
    vec<HalfData> hdata;
    auto hm_emb = compute_embedding_cut_hmesh(hm, em, seam, matching, singular, seam1, matching1, singular1, hdata);
    auto hm_cut = compute_cut_mesh(*hm_emb, seam1);
    lap("cut");

    MatXd uv = compute_tutte_parameterization(*hm_emb, em, seam1, hdata);
    lap("tutte");

    igl::SLIMData sData;
    MatXd uv_init(hm_cut->nV, 2);
    for (auto v: hm_cut->verts) uv_init.row(v.id) = uv.row(v.half().next().crnr().id);
    std::vector<int>   b_;   // the seam (boundary) vertices are pinned softly to the tutte uv
    std::vector<Row2d> bc_;
    for (auto v: hm_cut->verts) if (v.isBoundary()) { b_.push_back(v.id); bc_.emplace_back(uv_init.row(v.id)); }
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
        prev = sData.energy;
    }
    lap("slim");

    // per-corner uv from the per-vertex slim result: hm_cut and hm_emb share the face matrix
    VecXc cfn(hm_emb->nF * 3);
    for (int i = 0; i < hm_cut->nF; ++i)
    for (int j = 0; j < 3; ++j) cfn(i * 3 + j) = complex(sData.V_o(hm_cut->idx(i, j), 0), sData.V_o(hm_cut->idx(i, j), 1));
    qex::sanitization(*hm_emb, matching1, singular1, N, cfn);
    lap("qex sanitization");

    vec<qex::Qport> q_ports;
    vec<qex::Qvert> vqvs, eqvs, fqvs;
    qex::generate_q_vert(*hm_emb, cfn, vqvs, eqvs, fqvs);
    lap("qex q_vert");
    qex::generate_vqvert_qport(*hm_emb, cfn, vqvs, q_ports);
    qex::generate_eqvert_qport(*hm_emb, cfn, eqvs, q_ports);
    qex::generate_fqvert_qport(*hm_emb, fqvs, q_ports);
    lap("qex q_port");
    auto qedges = qex::generate_q_edge(*hm_emb, cfn, matching1, q_ports);
    lap("qex q_edge");
    auto qfaces = qex::generate_q_faces(q_ports, qedges);
    lap("qex q_face");
    print_time("total", std::chrono::steady_clock::now() - t_start);

    ///--- visualize ---///
    visualizer::visualize_init();
    visualizer::visualize_mesh(hm.pos, hm.idx, false, "base mesh");
    visualizer::visualize_tedges(hm, em, "tedges snapped", false);
    auto [qv, qidx] = visualizer::visualize_qfaces(hm, qfaces, true, false);
    visualizer::visualize_quad_patch(qv, qidx, label_quad_patches(em, singular, qfaces, qidx));
    polyscope::show();
    return 0;
}
