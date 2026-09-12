#include <string>
#include "pipeline.h"
#include "check.h"

using namespace metriko;

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: test_tmesh <gridscale> <vectorfield> <mesh.obj> [more.obj ...]\n"; return 2; }
    const int       N     = 4;
    const double    scale = std::stod(argv[1]);
    const FieldType ft    = parse_field_type(argv[2]);

    for (int a = 3; a < argc; ++a) {
        const char* mesh = argv[a];
        TmeshPipeline P(mesh, scale, ft);
        CHECK(P.ok);
        const Hmesh&     hm  = *P.hm;
        const Mgrph& mg  = *P.mg;
        const Tmesh&     tm  = *P.tm;
        Emesh&        em = *P.em;

        CHECK(tm.nTH == 2 * tm.nTE);
        CHECK(em.thalfs.size() == 2 * em.tedges.size());
        CHECK(em.tnodes.size() == mg.mnodes.size());
        CHECK(!em.tquads.empty());
        CHECK(em.tquads.size() == tm.tquads.size());

        if (ft == FieldType::Smoothest) {
            for (int i = 0; i < 5; ++i) {
                for (const Ehalf& th0: em.thalfs) {
                    if (th0.id == -1) continue;             // guards before any indexing:
                    auto& th1 = em.thalfs[th0.twid];       // cleared thalfs have twid/tqid == -1
                    if (th1.id == -1) continue;
                    if (th0.x != 0) continue;
                    auto& tq0 = em.tquads[th0.tqid];
                    auto& tq1 = em.tquads[th1.tqid];
                    if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
                    if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
                    em.collapse_thalf(th0.id);
                }

                for (const Equad& tq: em.tquads) {
                    if (tq.id == -1) continue;
                    Tqchain chain;
                    if (em.collapse_tquad_chain_prepare(tq.id, chain)) em.collapse_tquad_chain_execute(chain);
                }
            }
        } else {
            for (int i = 0; i < 20; ++i) {
                for (const Ehalf& th0 : em.thalfs) {
                    if (th0.id == -1) continue;
                    auto& th1 = em.thalfs[th0.twid];
                    if (th1.id == -1) continue;
                    if (th0.x != 0)   continue;
                    auto& tq0 = em.tquads[th0.tqid];
                    auto& tq1 = em.tquads[th1.tqid];
                    if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
                    if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
                    em.collapse_thalf(th0.id);
                }

                for (const Equad& tq: em.tquads) {
                    if (tq.id == -1) continue;
                    Tqchain chain;
                    if (em.collapse_tquad_chain_prepare(tq.id, chain)) em.collapse_tquad_chain_execute(chain);
                }
            }
        }


        auto opp_balanced = [&](const Emesh& m) -> bool {
            for (const Equad& q : m.tquads) {
                if (q.data.empty()) continue;
                if (q.id == -1) continue;
                double s[4] = {0, 0, 0, 0};
                for (const auto& [thid, side] : q.data) {
                    int x = m.thalfs[thid].x;
                    if (x == 0) {
                        std::println("thid: {}, tqid: {}", thid, m.thalfs[thid].tqid);
                        return false;
                    }
                    s[side] += x;
                }
                std::cout << "s[0]: " << s[0] << std::endl;
                std::cout << "s[1]: " << s[1] << std::endl;
                std::cout << "s[2]: " << s[2] << std::endl;
                std::cout << "s[3]: " << s[3] << std::endl;
                if (std::abs(s[0] - s[2]) > 1e-6 || std::abs(s[1] - s[3]) > 1e-6) { return false; }
            }
            return true;
        };
        CHECK(opp_balanced(em));

        std::cout << "[test_tmesh] OK  " << mesh << std::endl;
    }
    return 0;
}