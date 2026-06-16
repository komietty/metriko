#include <string>
#include "pipeline.h"
#include "check.h"

using namespace metriko;

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: test_tmesh <gridscale> <mesh.obj> [more.obj ...]\n"; return 2; }
    const int    N     = 4;
    const double scale = std::stod(argv[1]);

    for (int a = 2; a < argc; ++a) {
        const char* mesh = argv[a];
        TmeshPipeline P(mesh, scale, N);
        CHECK(P.ok);
        const Hmesh&     hm  = *P.hm;
        const mc::Mgrph& mg  = *P.mg;
        const Tmesh&     tm  = *P.tm;
        TmeshMut&        tmm = *P.tmm;

        CHECK(tm.nTH == 2 * tm.nTE);
        CHECK(tmm.thalfs.size() == 2 * tmm.tedges.size());
        CHECK(tmm.tnodes.size() == mg.mnodes.size());
        CHECK(!tmm.tquads.empty());
        CHECK(tmm.tquads.size() == tm.tquads.size());

        // collapse thalf
        for (ThalfMut th0 : tmm.thalfs) {
            auto& th1 = tmm.thalfs[th0.twid];
            auto& tq0 = tmm.tquads[th0.tqid];
            auto& tq1 = tmm.tquads[th1.tqid];
            if (th0.id == -1) continue;
            if (th1.id == -1) continue;
            if (th0.x != 0) continue;
            if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
            if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
            tmm.collapse_thalf(th0.id);
        }

        // collapse tquad
        for (const TquadMut& tq: tmm.tquads) {
            Tqaux tqaux;
            if (tmm.collapse_tquad_prepare(tq.id, tqaux))
                tmm.collapse_tquad_execute(tq.id, tqaux);
        }

        auto opp_balanced = [&](const TmeshMut& m) -> bool {
            for (const TquadMut& q : m.tquads) {
                if (q.data.empty()) continue;
                double s[4] = {0, 0, 0, 0};
                for (const auto& [thid, side] : q.data) {
                    int x = m.thalfs[thid].x;
                    if (x == 0) return false;
                    s[side] += x;
                }
                if (std::abs(s[0] - s[2]) > 1e-6 || std::abs(s[1] - s[3]) > 1e-6) { return false; }
            }
            return true;
        };
        CHECK(opp_balanced(tmm));

        std::cout << "[test_tmesh] OK  " << mesh << std::endl;
    }
    return 0;
}