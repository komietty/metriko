// end-to-end smoke test: lib.h must run from an OBJ to quad faces without
// throwing, and the t-mesh handed to the cut must be free of tedge contacts.
//
// usage: test_remesh <gridscale> <vectorfield> <mesh.obj> [more.obj ...]
#include <igl/readOBJ.h>
#include <string>
#include "metriko/lib.h"
#include "pipeline.h"   // parse_field_type
#include "check.h"

using namespace metriko;

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: test_remesh <gridscale> <vectorfield> <mesh.obj> [more.obj ...]\n"; return 2; }
    const double    scale = std::stod(argv[1]);
    const FieldType ft    = parse_field_type(argv[2]);

    for (int a = 3; a < argc; ++a) {
        const char* mesh = argv[a];
        MatXd V; MatXi F;
        igl::readOBJ(mesh, V, F);
        CHECK(V.rows() > 0 && F.rows() > 0);

        try {
            auto res = compute_remesh(V, F, scale, ft == FieldType::CurvatureAligned);

            CHECK(res.hmesh && res.mgrph && res.tmesh && res.emesh);
            CHECK(validate_no_crossing(*res.emesh, "test") == 0);                              // no tedge contacts
            CHECK(!res.q_faces.empty());                                                        // reached quad extraction
            CHECK(rg::all_of(res.q_ports, [](const qex::Qport& p) { return p.isConnected; }));  // every port paired
            for (const auto& qf: res.q_faces) CHECK(qf.qhalfs.size() == 4);

            std::cout << "[test_remesh] OK  " << mesh << "  quads=" << res.q_faces.size() << std::endl;
        } catch (const std::exception& ex) {
            std::cerr << "FAIL: compute_remesh threw: " << ex.what() << "  (" << mesh << ")\n";
            return 1;
        }
    }
    return 0;
}
