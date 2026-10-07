#ifndef METRIKO_LIB_H
#define METRIKO_LIB_H
#include <igl/slim.h>
#include "hmesh/hmesh.h"
#include "nvec/face_rosy_field.h"
#include "igm/parameterization.h"
#include "quantization/quantization.h"
#include "tmesh/emesh.h"
#include "tmesh/emesh_validate.h"
#include "metriko/tutte/tutte_cutting.h"
#include "metriko/tutte/tutte_params.h"
#include "metriko/qex/sanitization.h"
#include "metriko/qex/gen_q_vert.h"
#include "metriko/qex/gen_q_port.h"
#include "metriko/qex/gen_q_edge.h"
#include "metriko/qex/gen_q_face.h"
#include "patching.h"

namespace metriko {

constexpr int    slim_max_iter = 50;
constexpr double slim_rel_tol  = 1e-4;

struct RemeshStage {
    std::unique_ptr<Hmesh>         hmesh;
    std::unique_ptr<Mgrph>         mgrph;
    std::unique_ptr<Emesh>         emesh;
    std::unique_ptr<FaceRosyField> cmbf;
    vec<bool>                      seam;
    MatXd                          qnt_x;
    vec<qex::Qport>                q_ports;
    vec<qex::Qedge>                q_edges;
    vec<qex::Qface>                q_faces;
    MatXd                          q_pos;
    MatXi                          q_idx;
    MatXi                          q_val;
};

struct RemeshResult {
    MatXd pos;
    MatXi idx;
    MatXi val;
};

inline RemeshStage compute_quadrangulation_impl(
    const MatXd& V,
    const MatXi& F,
    const double scale,
    const bool calign = false,
    const bool refine = true
) {
    auto unit = grid_unit(V, scale);
    auto hm   = std::make_unique<Hmesh>(MatXd(V / unit), F);
    auto rawf = FaceRosyField(*hm, 4, calign? FieldType::CurvatureAligned : FieldType::Smoothest);
    auto seam = compute_seam(rawf);
    auto cutm = compute_cut_mesh(*hm, seam);
    auto cmbf = compute_combbed_field(rawf, seam);
    auto extf = compute_extrinsic_field(*cmbf, 4);

    auto rp = RosyParameterization(*hm, *cutm, extf, cmbf->singular, cmbf->matching, seam, 4, scale);
    rp.localInjectivity = true;
    rp.setup();
    rp.integ();

    hm->cfn = VecXc(hm->nF * 3);
    for (const Face f: hm->faces) {
        hm->cfn(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        hm->cfn(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        hm->cfn(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }

    auto mg = std::make_unique<Mgrph>(*hm, cmbf->matching, cmbf->singular);
    auto tm = std::make_unique<Emesh>(*mg);
    auto X  = compute_quantization(*tm, *mg);
    validate_quantization(*tm, X);
    tm->set_x(X);

    for (int i = 0; i < 100; ++i) {
        for (Ehalf th0 : tm->thalfs) {
            if (th0.id == -1 || th0.x != 0) continue;
            const auto& th1 = th0.twin();
            if (th0.tquad().thids(th0.side()).size() == 1) continue;
            if (th1.tquad().thids(th1.side()).size() == 1) continue;
            tm->collapse_thalf(th0.id);
        }

        for (const auto& [tqid, _] : tm->live_tquads())
            if (Tqaux a; tm->collapse_tquad_prepare(tqid, a)) tm->collapse_tquad_execute(a);

        if (rg::none_of(tm->thalfs, [](const Ehalf& th) { return th.id != -1 && th.x == 0; })) break;
    }
    validate_collapse_done(*tm);

    tm->collapse_tedge_snap(false);
    tm->collapse_tedge_snap(true);
    for (const auto& [teid, _] : tm->live_tedges()) { tm->collapse_tedge_snap_dedup(teid); }

    repair_crossing_tedges(*tm);

    ///--- cut the original mesh along the collapsed t-mesh ---///
    vec<bool> seam1;
    VecXi matching1;
    VecXi singular1;
    vec<HalfData> hdata;
    auto hm_emb = compute_embedding_cut_hmesh(*hm, *tm, seam, cmbf->matching, cmbf->singular, seam1, matching1, singular1, hdata);
    auto hm_cut = compute_cut_mesh(*hm_emb, seam1);

    MatXd uv = compute_tutte_parameterization(*hm_emb, *tm, seam1, hdata);

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
    slim_precompute(hm_cut->pos, hm_cut->idx, uv_init, sData, sData.slim_energy, b, bc, 1e9);

    double prev = std::numeric_limits<double>::infinity();
    for (int i = 0; i < slim_max_iter; ++i) {
        slim_solve(sData, 1);
        if (std::abs(prev - sData.energy) < slim_rel_tol * std::abs(sData.energy)) break;
        prev = sData.energy;
    }

    // ------ qex on the slim result ------
    hm_emb->cfn = VecXc(hm_emb->nF * 3);
    for (int i = 0; i < hm_cut->nF; ++i) {
    for (int j = 0; j < 3; ++j) {
        int k = hm_cut->idx(i, j);
        hm_emb->cfn(i * 3 + j) = complex(sData.V_o(k, 0), sData.V_o(k, 1));
    }}

    qex::sanitization(*hm_emb, matching1, singular1, 4);

    vec<qex::Qport> q_ports;
    vec<qex::Qvert> vqvs, eqvs, fqvs;
    qex::generate_q_vert(*hm_emb, vqvs, eqvs, fqvs);
    qex::generate_vqvert_qport(*hm_emb, vqvs, q_ports);
    qex::generate_eqvert_qport(*hm_emb, eqvs, q_ports);
    qex::generate_fqvert_qport(*hm_emb, fqvs, q_ports);
    auto q_edges = qex::generate_q_edge(*hm_emb, matching1, q_ports);
    auto q_faces = qex::generate_q_faces(q_ports, q_edges);
    auto [q_pos, q_idx] = extract_quad_mesh(*hm, q_faces, refine);
    auto q_val = label_quad_patches(*tm, cmbf->singular, q_faces, q_idx);
    q_pos *= unit;

    return {
        std::move(hm),
        std::move(mg),
        std::move(tm),
        std::move(cmbf),
        seam,
        X,
        std::move(q_ports),
        std::move(q_edges),
        std::move(q_faces),
        std::move(q_pos),
        std::move(q_idx),
        std::move(q_val)
    };
}

inline RemeshResult compute_quadrangulation(
    const MatXd& V,
    const MatXi& F,
    const double scale
) {
    auto res = compute_quadrangulation_impl(V, F, scale);
    return {
        std::move(res.q_pos),
        std::move(res.q_idx),
        std::move(res.q_val)
    };
}
}
#endif
