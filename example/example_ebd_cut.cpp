//
// Created by saki on 2026/07/14.
//
// demo: collapse the TmeshMut, then cut the original hmesh along the collapsed
// t-mesh (emesh_cutting) and display the resulting cut mesh with its patch
// boundaries.
#include <igl/readOBJ.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "metriko/core/emesh/emesh_cutting.h"

using namespace metriko;

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: example_ebd_cut <mesh.obj> <gridscale>\n"; return 1; }
    constexpr int N = 4;

    ///--- field -> parameterization (same pipeline as example_ebd) ---///
    MatXd V; MatXi F;
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);
    FaceRosyField rawf(hm, N, FieldType::Smoothest);
    rawf.computeMatching(MatchingType::Principal);
    auto seam = compute_seam(rawf);
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

    VecXc uv2(hm.nF * 3);
    for (const Face f: hm.faces) {
        uv2(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        uv2(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        uv2(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }

    ///--- motorcycle graph -> tmesh -> quantization ---///
    auto mg = mc::Mgrph(hm, uv2, cmbf->matching, cmbf->singular);
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
    auto hm_cut = emesh::compute_embedding_cut_hmesh(hm, tmm, seam, seam1, hdata);
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

    auto* cut = polyscope::registerSurfaceMesh("cut mesh", hm_cut->pos, hm_cut->idx);
    cut->setEdgeWidth(1.0);

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
