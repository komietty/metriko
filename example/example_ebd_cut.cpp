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
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/tutte_params.h"

using namespace metriko;

// binary cache for the expensive field -> parameterization stage (debug convenience).
// keyed by mesh path + gridscale; delete the .cache file when N or solver flags change.
static bool load_cache(const std::string& p, VecXc& uv2, VecXi& matching, VecXi& singular, std::vector<bool>& seam) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return false;
    int64_t nu, nm, ns, ne;
    f.read((char*)&nu, 8); f.read((char*)&nm, 8); f.read((char*)&ns, 8); f.read((char*)&ne, 8);
    uv2.resize(nu); matching.resize(nm); singular.resize(ns);
    f.read((char*)uv2.data(),      nu * (int64_t)sizeof(complex));
    f.read((char*)matching.data(), nm * (int64_t)sizeof(int));
    f.read((char*)singular.data(), ns * (int64_t)sizeof(int));
    std::vector<char> sb(ne);
    f.read(sb.data(), ne);
    seam.assign(sb.begin(), sb.end());
    return (bool)f;
}

static void save_cache(const std::string& p, const VecXc& uv2, const VecXi& matching, const VecXi& singular, const std::vector<bool>& seam) {
    std::ofstream f(p, std::ios::binary);
    int64_t nu = uv2.size(), nm = matching.size(), ns = singular.size(), ne = (int64_t)seam.size();
    f.write((char*)&nu, 8); f.write((char*)&nm, 8); f.write((char*)&ns, 8); f.write((char*)&ne, 8);
    f.write((char*)uv2.data(),      nu * (int64_t)sizeof(complex));
    f.write((char*)matching.data(), nm * (int64_t)sizeof(int));
    f.write((char*)singular.data(), ns * (int64_t)sizeof(int));
    std::vector<char> sb(seam.begin(), seam.end());
    f.write(sb.data(), ne);
}

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: example_ebd_cut <mesh.obj> <gridscale>\n"; return 1; }
    constexpr int N = 4;

    ///--- field -> parameterization (same pipeline as example_ebd; cached) ---///
    MatXd V; MatXi F;
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);

    VecXc uv2; VecXi matching, singular; std::vector<bool> seam;
    const std::string cache = std::format("{}.{}.cache", argv[1], argv[2]);

    if (load_cache(cache, uv2, matching, singular, seam)) {
        std::println("loaded cache: {}", cache);
    } else {
        FaceRosyField rawf(hm, N, FieldType::Smoothest);
        rawf.computeMatching(MatchingType::Principal);
        seam = compute_seam(rawf);
        auto cutm = compute_cut_mesh(hm, seam);
        auto cmbf = compute_combbed_field(rawf, seam);

        MatXd ext(hm.nF, 3 * N);
        for (Face f: hm.faces)
            for (int k = 0; k < N; ++k) {
                complex c = cmbf->field(f.id, k);
                ext.block(f.id, 3 * k, 1, 3) = (c.real() * f.basisX() + c.imag() * f.basisY()).normalized();
            }

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

        matching = cmbf->matching;
        singular = cmbf->singular;
        save_cache(cache, uv2, matching, singular, seam);
        std::println("saved cache: {}", cache);
    }

    ///--- motorcycle graph -> tmesh -> quantization ---///
    auto mg = mc::Mgrph(hm, uv2, matching, singular);
    auto tm = Tmesh(mg);
    VecXd X = compute_quantization(tm, mg);
    TmeshMut tmm(mg, tm, X);

    ///--- collapse: rounds of (thalf pass + chain collapse); later rounds pick up
    ///    the zero edges created by earlier collapses ---///
    for (int i = 0; i < 5; ++i) {
        for (const ThalfMut& th0: tmm.thalfs) {
            if (th0.id == -1) continue;
            auto& th1 = tmm.thalfs[th0.twid];
            if (th1.id == -1) continue;
            if (th0.x != 0)   continue;
            auto& tq0 = tmm.tquads[th0.tqid];
            auto& tq1 = tmm.tquads[th1.tqid];
            if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
            if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
            tmm.collapse_thalf(th0.id);
        }
        for (const TquadMut& tq: tmm.tquads) {
            if (tq.id == -1) continue;
            Tqchain chain;
            if (tmm.collapse_tquad_chain_prepare(tq.id, chain)) tmm.collapse_tquad_chain_execute(chain);
        }
    }

    for (const auto& [teid, nids] : tmm.tedges) { if (!nids.empty()) tmm.collapse_tedge_short_segment(teid); }
    for (const auto& [teid, nids] : tmm.tedges) { if (!nids.empty()) tmm.collapse_tedge_vert_snapping(teid); }
    for (const auto& [teid, nids] : tmm.tedges) { if (!nids.empty()) tmm.collapse_tedge_short_segment(teid); }

    ///--- validate tedge nid chains: duplicated / backtracking nodes break the cut ---///
    for (const auto& th: tmm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids = tmm.tedges[th.teid].nids;
        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            Row3d a = get_ptloc_pos(hm, tmm.tnodes[nids[k]]);
            Row3d b = get_ptloc_pos(hm, tmm.tnodes[nids[k + 1]]);
            if (nids[k] == nids[k + 1] || (a - b).norm() < 1e-12)
                std::println("[warn] tedge {} (thid {}, x {}): duplicate node at {} (nid {} / {}, locs {} / {})",
                             th.teid, th.id, th.x, k, nids[k], nids[k + 1],
                             loc_str(tmm.tnodes[nids[k]]), loc_str(tmm.tnodes[nids[k + 1]]));
        }
        for (size_t k = 0; k + 2 < nids.size(); ++k) {
            if (nids[k] == nids[k + 2])
                std::println("[warn] tedge {} (thid {}, x {}): backtrack at {} (nid {})", th.teid, th.id, th.x, k, nids[k]);
        }
    }

    ///--- cut the original mesh along the collapsed t-mesh ---///
    vec<bool> seam1;
    vec<HalfData> hdata;
    auto hm_cut = compute_embedding_cut_hmesh(hm, tmm, seam, seam1, hdata);
    std::println("cut mesh: V {} F {} (original: V {} F {})", hm_cut->pos.rows(), hm_cut->idx.rows(), hm.pos.rows(), hm.idx.rows());
    std::println("patch boundary halfedges: {}", hdata.size());
    std::println("euler characteristic: {} (closed manifold iff 2E == 3F: {})",
                 hm_cut->pos.rows() - (long)hm_cut->nE + hm_cut->idx.rows(),
                 2 * hm_cut->nE == 3 * hm_cut->idx.rows());

    ///--- validate: every halfedge must have its opposite pair.
    ///    unpaired -> hole, duplicated -> flipped face / non-manifold ---///
    {
        std::map<std::pair<int, int>, int> cnt;
        for (int i = 0; i < hm_cut->idx.rows(); ++i)
            for (int j = 0; j < 3; ++j) {
                int a = hm_cut->idx(i, j);
                int b = hm_cut->idx(i, (j + 1) % 3);
                cnt[{a, b}]++;
            }
        int unpaired = 0, duplicated = 0;
        for (const auto& [k, c]: cnt) {
            if (c > 1)                            { if (duplicated++ < 10) std::println("[error] duplicated halfedge ({} -> {}) x{}", k.first, k.second, c); }
            if (!cnt.contains({k.second, k.first})) { if (unpaired++   < 10) std::println("[error] unpaired halfedge ({} -> {})", k.first, k.second); }
        }
        if (unpaired || duplicated)
            throw std::runtime_error(std::format("cut mesh validation failed: {} unpaired, {} duplicated halfedges", unpaired, duplicated));
        std::println("halfedge pairing: all {} halfedges paired", cnt.size());
    }

    ///--- visualize ---///
    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

    auto* base = polyscope::registerSurfaceMesh("base mesh", hm.pos, hm.idx);
    base->setEnabled(false);

    //{ // TEMP: face 5018 segments
    //    std::vector<glm::vec3> ns; std::vector<std::array<size_t, 2>> es; size_t c = 0;
    //    for (const auto& th: tmm.thalfs) {
    //        if (th.id == -1 || !th.cano) continue;
    //        const auto& nids = tmm.tedges[th.teid].nids;
    //        for (size_t k = 0; k + 1 < nids.size(); ++k) {
    //            if (emesh::common_face(hm, tmm.tnodes[nids[k]], tmm.tnodes[nids[k + 1]]) != 5018) continue;
    //            Row3d a = get_ptloc_pos(hm, tmm.tnodes[nids[k]]);
    //            Row3d b = get_ptloc_pos(hm, tmm.tnodes[nids[k + 1]]);
    //            ns.emplace_back(a.x(), a.y(), a.z()); ns.emplace_back(b.x(), b.y(), b.z());
    //            es.push_back({c, c + 1}); c += 2;
    //        }
    //    }
    //    polyscope::registerCurveNetwork("face 5018 segments", ns, es);
    //}

    auto* cut = polyscope::registerSurfaceMesh("cut mesh", hm_cut->pos, hm_cut->idx);
    cut->setEdgeWidth(1.0);

    ///--- tutte parameterization (pre-SLIM initial uv) ---///
    std::sort(hdata.begin(), hdata.end());   // compute_tutte_parameterization groups halfedges via equal_range on tqid
    MatXd uv;
    if (compute_tutte_parameterization(*hm_cut, tmm, seam1, hdata, uv)) {
        cut->addParameterizationQuantity("tutte uv", uv);

        ///--- SLIM (refine the tutte uv) ---///
        {
            // cut along the seams so the uv becomes per-vertex
            auto hm2 = compute_cut_mesh(*hm_cut, seam1);

            MatXd uv_init(hm2->nV, 2);
            for (auto v: hm2->verts) uv_init.row(v.id) = uv.row(v.half().next().crnr().id);

            // pin the seam (boundary) vertices softly to the tutte uv
            std::vector<int> b_;
            std::vector<Row2d> bc_;
            for (auto v: hm2->verts) {
                if (!v.isBoundary()) continue;
                b_.push_back(v.id);
                bc_.push_back(uv_init.row(v.id));
            }
            VecXi b = Eigen::Map<VecXi>(b_.data(), (int)b_.size());
            MatXd bc((int)bc_.size(), 2);
            for (int i = 0; i < (int)bc_.size(); ++i) bc.row(i) = bc_[i];

            igl::SLIMData sData;
            sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;
            slim_precompute(hm2->pos, hm2->idx, uv_init, sData, sData.slim_energy, b, bc, 1e5);
            slim_solve(sData, 100);
            std::println("[slim] displacement: {}", (sData.V_o - uv_init).norm());

            auto* surf = polyscope::registerSurfaceMesh("slim result", hm2->pos, hm2->idx);
            auto* prms = surf->addVertexParameterizationQuantity("uv", sData.V_o);
            surf->setEdgeWidth(0.7);
            prms->setEnabled(true);
            prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
            prms->setCheckerSize(1);
        }
    } else {
        std::println("[tutte] compute_tutte_parameterization failed");
    }

    { // patch boundaries (t-mesh edges on the cut mesh), colored by tquad id
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> tqids, thids, v0s;
        size_t c = 0;
        for (const auto& d: hdata) {
            Row3d a = d.half.tail().pos();
            Row3d b = d.half.head().pos();
            ns.emplace_back(a.x(), a.y(), a.z());
            ns.emplace_back(b.x(), b.y(), b.z());
            es.push_back({c, c + 1}); c += 2;
            tqids.push_back(d.tqid);
            thids.push_back(d.thid);
            v0s.push_back(d.v0);
        }
        auto* cn = polyscope::registerCurveNetwork("patch boundaries", ns, es);
        cn->addEdgeScalarQuantity("tqid", tqids)->setEnabled(true);
        cn->addEdgeScalarQuantity("thid", thids);
        cn->addEdgeScalarQuantity("v0", v0s);
        cn->setRadius(0.0015);
    }

    { // seam edges propagated onto the cut mesh
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t c = 0;
        for (Edge e: hm_cut->edges) {
            if (!seam1[e.id]) continue;
            Row3d a = e.vert0().pos();
            Row3d b = e.vert1().pos();
            ns.emplace_back(a.x(), a.y(), a.z());
            ns.emplace_back(b.x(), b.y(), b.z());
            es.push_back({c, c + 1}); c += 2;
        }
        if (!ns.empty()) {
            auto* cn = polyscope::registerCurveNetwork("seam (cut)", ns, es);
            cn->setColor({0.9, 0.3, 0.2});
            cn->setRadius(0.001);
            cn->setEnabled(false);
        }
    }

    polyscope::show();
    return 0;
}
