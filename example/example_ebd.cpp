#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <igl/upsample.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/hmesh/hpath.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "metriko/core/tutte/tutte.h"
#include "common.h"

using namespace metriko;
static int N = 4;
static VecXc uv2;
static MatXd V;
static MatXi F;

static void validate_tmeshmut(const TmeshMut& tm) {
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
    Hmesh hm(V, F);

    VecXi matching;
    VecXi singular;
    std::vector<bool> seam;

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

    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

    ///--- gen mport, medge ---///
    auto mg = mc::Mgrph(hm, uv2, matching, singular);
    auto tm = Tmesh(mg);
    VecXd X = compute_quantization(tm, mg);
    validate_quantization(tm, X);
    assert(tm.check_non_zero_tquad(X));
    visualizer::visualize_tedge(tm, mg, uv2, &X);

    TmeshMut tmm(mg, tm, X);

    std::println("---- collapse bgn");

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

    std::println("---- collapse end");

    // debug view
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> ids;
        size_t c = 0;
        for (const auto& [id, nids]: tmm.live_tedges()) {
            for (size_t k = 0; k + 1 < nids.size(); ++k) {
                Row3d a = get_ptloc_pos(hm, tmm.tnodes[nids[k]]);
                Row3d b = get_ptloc_pos(hm, tmm.tnodes[nids[k + 1]]);
                ns.emplace_back(a.x(), a.y(), a.z());
                ns.emplace_back(b.x(), b.y(), b.z());
                es.push_back({c, c + 1});
                c += 2;
                ids.push_back(id);
            }
        }
        auto* cn = polyscope::registerCurveNetwork("tedges_collapsed", ns, es);
        cn->addEdgeScalarQuantity("teid", ids)->setEnabled(true);
        cn->setRadius(0.0015);
        cn->resetTransform();
    }

    auto t0 = std::chrono::steady_clock::now();
    tmm.collapse_tedge_snap(false);
    tmm.collapse_tedge_snap(true);
    for (const auto& [teid, _] : tmm.live_tedges()) { tmm.collapse_tedge_snap_dedup(teid); }
    auto t1 = std::chrono::steady_clock::now();
    std::println("[time] snap: {:.3f}s", std::chrono::duration<double>(t1 - t0).count());

    vec<bool> seam1;
    vec<HalfData> hdata;
    //auto hm_emb = compute_embedding_cut_hmesh(hm, tmm, seam, seam1, hdata);
    //auto hm_cut = compute_cut_mesh(*hm_emb, seam1);

    auto* base = polyscope::registerSurfaceMesh("base mesh", hm.pos, hm.idx);           base->setEnabled(false);
    //auto* embd = polyscope::registerSurfaceMesh("embd mesh", hm_emb->pos, hm_emb->idx); embd->setEdgeWidth(1.);

    { // tnodes not snapped to a vertex, by carrier type
        std::vector<glm::vec3> ps;
        std::vector<double> type, ids, nids_;
        std::set<int> seen;   // shared nodes (junctions/crossings) appear in several chains
        for (const auto& [teid, nids]: tmm.tedges) {
            if (teid == -1) continue;
            for (int nid: nids) {
                if (!seen.insert(nid).second) continue;
                double t = -1, id = -1;
                std::visit(overloaded{
                    [&](const HmLocOnE& e) { t = 0; id = e.id; },
                    [&](const HmLocOnH& h) { t = 1; id = h.id; },
                    [&](const HmLocOnF& f) { t = 2; id = f.id; },
                    [&](const auto&)       {},
                }, tmm.tnodes[nid]);
                if (t < 0) continue;   // OnV: snapped, skip
                Row3d p = get_ptloc_pos(hm, tmm.tnodes[nid]);
                ps.emplace_back(p.x(), p.y(), p.z());
                type.push_back(t);
                ids.push_back(id);
                nids_.push_back(nid);
            }
        }
        auto* pc = polyscope::registerPointCloud("unsnapped tnodes", ps);
        pc->addScalarQuantity("type (0:E 1:H 2:F)", type)->setEnabled(true);
        pc->addScalarQuantity("elem id", ids);
        pc->addScalarQuantity("node id", nids_);
        pc->setPointRadius(0.002);
    }

    // debug view
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> ids;
        size_t c = 0;
        for (const auto& [id, nids]: tmm.live_tedges()) {
            for (size_t k = 0; k + 1 < nids.size(); ++k) {
                Row3d a = get_ptloc_pos(hm, tmm.tnodes[nids[k]]);
                Row3d b = get_ptloc_pos(hm, tmm.tnodes[nids[k + 1]]);
                ns.emplace_back(a.x(), a.y(), a.z());
                ns.emplace_back(b.x(), b.y(), b.z());
                es.push_back({c, c + 1});
                c += 2;
                ids.push_back(id);
            }
        }
        auto* cn = polyscope::registerCurveNetwork("tedges_snapped", ns, es);
        cn->addEdgeScalarQuantity("teid", ids)->setEnabled(true);
        cn->setRadius(0.0015);
        cn->resetTransform();
    }
    /*
    for (const auto& [id, data] : tmm.live_tquads()) {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> eside, ex, ey, er, ethid;   // per-edge (thalf) params
        size_t c = 0;
        for (const TdataMut& d : data) {
            const ThalfMut& th = tmm.thalfs[d.thid];
            const TedgeMut& te = tmm.tedges[th.teid];
            for (size_t i = 0; i + 1 < te.nids.size(); ++i) {
                Row3d a = get_ptloc_pos(hm, tmm.tnodes[te.nids[i]]);
                Row3d b = get_ptloc_pos(hm, tmm.tnodes[te.nids[i + 1]]);
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
        auto* cn = polyscope::registerCurveNetwork(std::format("tq {:03}", id), ns, es);
        cn->addEdgeScalarQuantity("side", eside);
        auto cx = cn->addEdgeScalarQuantity("x", ex);
        auto cy = cn->addEdgeScalarQuantity("y", ey);
        cy->setEnabled(true);
        cn->addEdgeScalarQuantity("r", er);
        cn->addEdgeScalarQuantity("thid", ethid);
        cn->setMaterial("flat");
        cn->setRadius(0.001); cn->resetTransform();
    }
    */

    polyscope::show(); return 0;
}