#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <igl/upsample.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "common.h"
#include "visualize_hmesh.h"
#include "visualize_tmesh_mut.h"

using namespace metriko;
static VecXc uv2;
static MatXd V;
static MatXi F;
static VecXi matching;
static VecXi singular;
static vec<bool> seam;

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
            throw std::runtime_error("s[0] != s[2] || s[1] != s[3]");
        }
    }
}

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);
    if (!load_cache(std::format("{}.{}.cache", argv[1], argv[2]), uv2, matching, singular, seam)) throw std::runtime_error("the cache does not exist");

    ///--- gen mport, medge ---///
    auto mg = mc::Mgrph(hm, uv2, matching, singular);
    auto tm = Tmesh(mg);
    auto X  = compute_quantization(tm, mg);
    validate_quantization(tm, X);
    assert(tm.check_non_zero_tquad(X));

    visualizer::visualize_init();
    visualizer::visualize_tedge(tm, mg, uv2, &X);

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

    visualizer::visualize_mesh(hm.pos, hm.idx);
    visualizer::visualize_non_snapped_tnodes(hm, tmm, false);
    visualizer::visualize_tedge_mut_snapped(hm, tmm, true);
    visualizer::visualize_tedge_mut_collapsed(hm, tmm, false);
    visualizer::visualize_tquad_mut_collapsed(hm, tmm, false);
    visualizer::visualize_face_collinear_error(hm, tmm, true);

    polyscope::show(); return 0;
}