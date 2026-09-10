#ifndef HMLOC_H_LIB_H
#define HMLOC_H_LIB_H
#include <igl/slim.h>
#include "core/hmesh/hmesh.h"
#include "core/vectorfield/face_rosy_field.h"
#include "core/igm/parameterization.h"
#include "core/quantization/quantization.h"
#include "core/tmesh/tmesh.h"
#include "core/tmesh/tmesh_mut.h"
#include "core/tmesh/tmesh_mut_validate.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/tutte_params.h"
#include "metriko/core/qex/sanitization.h"
#include "metriko/core/qex/gen_q_vert.h"
#include "metriko/core/qex/gen_q_port.h"
#include "metriko/core/qex/gen_q_edge.h"
#include "metriko/core/qex/gen_q_face.h"

namespace metriko {

struct RemeshResult {
    std::unique_ptr<Hmesh>         hmesh;
    std::unique_ptr<Tmesh>         tmesh;
    std::unique_ptr<TmeshMut>      emesh;
    std::unique_ptr<FaceRosyField> cmbf;
    std::unique_ptr<mc::Mgrph>     mgrph;
    vec<bool>                      seam;
    MatXd                          cfn_d;
    VecXc                          cfn_c;
    MatXd                          qnt_x;
    vec<qex::Qport>                q_ports;
    vec<qex::Qedge>                q_edges;
    vec<qex::Qface>                q_faces;
};

inline RemeshResult compute_remesh(
    const MatXd& V,
    const MatXi& F,
    const double scale,
    const bool c_aligned = false
) {
    auto res  = RemeshResult();
    auto hm   = std::make_unique<Hmesh>(V, F);
    auto rawf = FaceRosyField(*hm, 4, c_aligned? FieldType::CurvatureAligned : FieldType::Smoothest); rawf.computeMatching(MatchingType::Principal);
    auto seam = compute_seam(rawf);
    auto cutm = compute_cut_mesh(*hm, seam);
    auto cmbf = compute_combbed_field(rawf, seam);
    auto extf = compute_extrinsic_field(*cmbf, 4);

    auto rp = RosyParameterization(*hm, *cutm, extf, cmbf->singular, cmbf->matching, seam, 4, scale);
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();

    VecXc cfn_c(hm->nF * 3);
    for (const Face f: hm->faces) {
        cfn_c(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        cfn_c(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        cfn_c(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }

    auto mg = std::make_unique<mc::Mgrph>(*hm, cfn_c, cmbf->matching, cmbf->singular);
    auto tm = std::make_unique<Tmesh>(*mg);
    auto X  = compute_quantization(*tm, *mg);
    validate_quantization(*tm, X);


    auto em = std::make_unique<TmeshMut>(*mg, *tm, X);

    for (int i = 0; i < 10; ++i) {
        // collapse thalf
        for (ThalfMut th0 : em->thalfs) {
            if (th0.id == -1 || th0.x != 0) continue;
            auto& th1 = em->thalfs[th0.twid];
            auto& tq0 = em->tquads[th0.tqid];
            auto& tq1 = em->tquads[th1.tqid];
            if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
            if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
            em->collapse_thalf(th0.id);
        }

        // collapse tquad
        for (const auto& [tqid, data] : em->live_tquads()) {
            Tqchain chain;
            if (em->collapse_tquad_chain_prepare(tqid, chain)) em->collapse_tquad_chain_execute(chain);
        }
    }

    em->collapse_tedge_snap(false);
    em->collapse_tedge_snap(true);
    for (const auto& [teid, _] : em->live_tedges()) { em->collapse_tedge_snap_dedup(teid); }

    repair_crossing_tedges(*em);


    ///--- cut the original mesh along the collapsed t-mesh ---///
    vec<bool> seam1;
    VecXi matching1;
    VecXi singular1;
    vec<HalfData> hdata;
    auto hm_emb = compute_embedding_cut_hmesh(*hm, *em, seam, cmbf->matching, cmbf->singular, seam1, matching1, singular1, hdata);
    auto hm_cut = compute_cut_mesh(*hm_emb, seam1);

    MatXd uv;
    bool flag = compute_tutte_parameterization(*hm_emb, *em, seam1, hdata, uv);
    if (!flag) throw std::runtime_error("compute_tutte parameterization failed");

    igl::SLIMData sData;
    MatXd uv_init(hm_cut->nV, 2);
    for (auto v: hm_cut->verts) uv_init.row(v.id) = uv.row(v.half().next().crnr().id);

    vec<int>   b_;
    vec<Row2d> bc_;
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
    slim_solve(sData, 50);

    // ------ qex on the slim result ------
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
    qex::generate_vqvert_qport(*hm_emb, cfn, vqvs, q_ports);
    qex::generate_eqvert_qport(*hm_emb, cfn, eqvs, q_ports);
    qex::generate_fqvert_qport(*hm_emb, fqvs, q_ports);
    auto q_edges = qex::generate_q_edge(*hm_emb, cfn, matching1, q_ports);
    auto q_faces = qex::generate_q_faces(q_ports, q_edges);

    res.hmesh = std::move(hm);
    res.cmbf  = std::move(cmbf);
    res.mgrph = std::move(mg);
    res.tmesh = std::move(tm);
    res.emesh = std::move(em);
    res.seam  = seam;
    res.cfn_d = std::move(rp.cfn);
    res.cfn_c = std::move(cfn_c);
    res.qnt_x = std::move(X);
    res.q_ports = std::move(q_ports);
    res.q_edges = std::move(q_edges);
    res.q_faces = std::move(q_faces);
    return res;
}
}
#endif
