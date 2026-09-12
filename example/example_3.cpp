//
// example_3: full pipeline in one run.
//   stage 0 (example_0): vectorfield + rosy parameterization
//   stage 1 (example_1): quantization + t-mesh collapse / snap
//   stage 2 (example_2): cutting + tutte + slim + qex
// usage: example_3 <mesh.obj> <scale>
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
#include "visualize_vectorfield.h"
#include "visualize_qex.h"
#include "visualize_tmesh_mut.h"
#include "visualize_quad_patch.h"

using namespace metriko;
static int N = 4;
static MatXd V;
static MatXi F;
static VecXc uv2;

static void validate_tmeshmut(const TmeshMut& tm) {
    for (const TquadMut& tq : tm.tquads) {
        if (tq.data.empty()) continue;
        double s[4] = {0, 0, 0, 0};
        for (const auto& [thid, side] : tq.data) {
            assert(thid != -1);
            auto& th_cur = tm.thalfs[thid];        if (th_cur.twid == -1 || th_cur.teid == -1 || th_cur.tqid == -1) throw std::runtime_error("th_cur");
            auto& th_twn = tm.thalfs[th_cur.twid]; if (th_twn.twid == -1 || th_twn.teid == -1 || th_twn.tqid == -1) throw std::runtime_error("th_twn");
            s[side] += th_cur.x;
        }
        if (std::abs(s[0] - s[2]) > 1e-6 || std::abs(s[1] - s[3]) > 1e-6) {
            std::println("s0: {}, s1: {}, s2: {}, s2: {}, tqid: {}", s[0], s[1], s[2], s[3], tq.id);
            throw std::runtime_error("s[0] != s[2] || s[1] != s[3]");
        }
    }
}

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);

    ///--- stage 0: vectorfield + rosy parameterization ---///
    FaceRosyField rawf(hm, N, FieldType::Smoothest);
    rawf.computeMatching(MatchingType::Principal);
    auto seam = compute_seam(rawf);
    auto cutm = compute_cut_mesh(hm, seam);
    auto cmbf = compute_combbed_field(rawf, seam);

    MatXd ext(hm.nF, 3 * N);
    for (Face f: hm.faces) {
    for (int k = 0; k < N; ++k) {
        complex c = cmbf->field(f.id, k);
        ext.block(f.id, 3 * k, 1, 3) = (c.real() * f.basisX() + c.imag() * f.basisY()).normalized();
    }}

    RosyParameterization rp(hm, *cutm, ext, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = false;
    rp.localInjectivity = true;
    rp.verbose = false;
    rp.setup();
    rp.integ();

    uv2.resize(hm.nF * 3);
    for (const Face f: hm.faces) {
        uv2(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        uv2(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        uv2(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }
    const VecXi& matching = cmbf->matching;
    const VecXi& singular = cmbf->singular;

    save_cache(std::format("{}.{}.cache", argv[1], argv[2]), uv2, cmbf->matching, cmbf->singular, seam);
    std::println("saved cache");

    ///--- stage 1: quantization + t-mesh collapse / snap ---///
    auto mg = mc::Mgrph(hm, uv2, matching, singular);
    auto tm = Tmesh(mg);
    auto X  = compute_quantization(tm, mg);
    validate_quantization(tm, X);
    assert(tm.check_non_zero_tquad(X));

    TmeshMut tmm(mg, tm, X);

    for (int i = 0; i < 10; ++i) {
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
        for (const auto& [tqid, data] : tmm.tquads) {
            Tqchain chain;
            if (tmm.collapse_tquad_chain_prepare(tqid, chain)) {
                std::cout << "tq collapse: " << tqid << std::endl;
                tmm.collapse_tquad_chain_execute(chain);
                validate_tmeshmut(tmm);
            }
        }
    }

    tmm.collapse_tedge_snap(false);
    tmm.collapse_tedge_snap(true);
    for (const auto& [teid, _] : tmm.live_tedges()) { tmm.collapse_tedge_snap_dedup(teid); }

    save_tmm(std::format("{}.{}.tmm", argv[1], argv[2]), tmm);
    std::println("saved tmm cache");

    ///--- stage 2: cut along the collapsed t-mesh ---///
    // validate tedge nid chains: duplicated / backtracking nodes break the cut
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

    ///--- visualize stages 0 and 1 ---///
    visualizer::visualize_init();
    auto surf0 = visualizer::visualize_mesh_with_uv(hm.pos, hm.idx, uv2);
    visualizer::visualize_frosy_field(surf0, hm, rawf, *cmbf);
    visualizer::visualize_seam(hm, seam);
    visualizer::visualize_tedge(tm, mg, uv2, &X);
    visualizer::visualize_tedge_mut_snapped(hm, tmm, true);
    visualizer::visualize_tedge_mut_collapsed(hm, tmm, false);

    ///--- stage 2 continued: tutte parameterization (pre-SLIM initial uv) ---///
    MatXd uv;
    try {
        uv = compute_tutte_parameterization(*hm_emb, tmm, seam1, hdata);
    } catch (std::exception& e) {
        std::println(stderr, "[tutte] {}", e.what());
        polyscope::show();
        return 1;
    }

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
        auto* tutt = surf->addVertexParameterizationQuantity("tutte uv", uv_init);
        surf->setEnabled(false);
        surf->setEdgeWidth(0.7);
        prms->setEnabled(true);
        prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
        prms->setCheckerSize(1);
        tutt->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
        tutt->setCheckerSize(1);
    }

    // ------ qex on the slim result ------
    {
        // per-corner uv from the per-vertex slim result: hm_cut and hm_emb
        // share the face matrix, so corner (i, j) <-> vertex hm_cut->idx(i, j)
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

        { // TODO TEMP: validate port cycles per qvert
            int i = 0;
            while (i < (int)q_ports.size()) {
                int j = i;
                while (j < (int)q_ports.size() && (q_ports[j].pos - q_ports[i].pos).norm() < 1e-12) ++j;
                const int s = j - i;

                int cur = q_ports[i].idx, cnt = 0;
                do { cur = q_ports[cur].next_id; ++cnt; } while (cur != q_ports[i].idx && cnt <= s);
                if (cnt != s) {
                    std::println("[qport] broken cycle: group at port {} (size {}, vid {}, eid {}, fid {})",
                                 q_ports[i].idx, s, q_ports[i].vid, q_ports[i].eid, q_ports[i].fid);
                    visualizer::visualize_qport_group(*hm_emb, cfn, q_ports, q_ports[i].idx);
                }

                // angular order: the cycle must be a rotation of the
                // angle-sorted order. counting descents tolerates a gap
                // wider than pi, which is legitimate for small groups
                // (e.g. 3 ports at a valence-3 singularity)
                vec<Row3d> dirs;
                for (int k = i; k < j; ++k) {
                    auto& p = q_ports[k];
                    Row3d d = (conversion_2d_3d(hm_emb->faces[p.fid], cfn, p.uv + p.dir)
                             - conversion_2d_3d(hm_emb->faces[p.fid], cfn, p.uv)).normalized();
                    dirs.push_back(d);
                }
                Row3d n = Row3d::Zero();
                for (int k = 0; k < s; ++k) n += dirs[k].cross(dirs[(k + 1) % s]);
                n.normalize();
                Row3d bx = (dirs[0] - n * n.dot(dirs[0])).normalized();
                Row3d by = n.cross(bx);
                int descents = 0;
                for (int k = 0; k < s; ++k) {
                    double t0 = std::atan2(dirs[k].dot(by), dirs[k].dot(bx));
                    double t1 = std::atan2(dirs[(k + 1) % s].dot(by), dirs[(k + 1) % s].dot(bx));
                    if (t1 < t0) ++descents;
                }
                if (descents != 1) {
                    std::println("[qport] non-CCW cycle: group at port {} (vid {}, eid {}, fid {}, descents {})",
                                 q_ports[i].idx, q_ports[i].vid, q_ports[i].eid, q_ports[i].fid, descents);
                    visualizer::visualize_qport_group(*hm_emb, cfn, q_ports, q_ports[i].idx);
                }
                i = j;
            }
        }

        auto qedges = qex::generate_q_edge(*hm_emb, cfn, matching1, q_ports);
        auto qfaces = qex::generate_q_faces(q_ports, qedges);

        visualizer::visualize_qedges(qedges);
        visualizer::visualize_quad_patch(hm, tmm, singular, qfaces);
    }

    polyscope::show(); return 0;
}
