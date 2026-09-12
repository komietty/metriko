// Compute region(allowed_range) for every tquad of the given mesh and save it to JSON.
// Each entry keeps the (eid, r0, r1) triple:
//   { "<tqid>": [[eid, r0, r1], [eid, r0, r1], ...], ... }
//   - array elements are sorted ascending by (eid, r0, r1).
//
// usage: dump_region <mesh.obj> <gridscale> <out.json>
//
// build & run (from repository root):
//   cmake --build test/cmake-build-debug --target dump_region -j8
//   ./test/cmake-build-debug/dump_region /Users/saki/dev/models/spot.obj 0.05 test/fixtures/tquad_region/spot.json
//   ./test/cmake-build-debug/dump_region /Users/saki/dev/models/torus.obj 0.05 test/fixtures/tquad_region/torus.json
//   note: only right after editing CMakeLists, reconfigure first: cmake -S test -B test/cmake-build-debug
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "pipeline.h"

using namespace metriko;

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: dump_region <mesh.obj> <gridscale> <vectorfield> <out.json>\n"; return 2; }

    TmeshPipeline P(argv[1], std::stod(argv[2]), parse_field_type(argv[3]));
    if (!P.ok) { std::cerr << "failed to load mesh: " << argv[1] << "\n"; return 1; }
    auto& em = *P.em;

    std::ofstream out(argv[4]);
    if (!out.good()) { std::cerr << "cannot open output: " << argv[4] << "\n"; return 1; }
    out << std::setprecision(17);

    const int n = (int)em.tquads.size();
    out << "{\n";
    for (int tqid = 0; tqid < n; ++tqid) {
        auto rng = em.allowed_range_tquads({tqid});     // vec<tuple<eid, r0, r1>>
        std::sort(rng.begin(), rng.end());      // ascending by (eid, r0, r1)
        out << "  \"" << tqid << "\": [";
        for (size_t i = 0; i < rng.size(); ++i) {
            auto& [eid, r0, r1] = rng[i];
            out << (i ? "," : "") << "[" << eid << "," << r0 << "," << r1 << "]";
        }
        out << "]" << (tqid + 1 < n ? "," : "") << "\n";
    }
    out << "}\n";

    std::cout << "[dump_region] " << argv[1] << " tquads=" << n << " -> " << argv[4] << "\n";
    return 0;
}
