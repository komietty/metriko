// end-to-end smoke test: lib.h must run from an OBJ to quad faces without
// throwing, and the t-mesh handed to the cut must be free of tedge contacts.
//
// usage: test_remesh <gridscale> <nvec> <mesh.obj> [more.obj ...]
#include <igl/readOBJ.h>
#include <string>
#include "metriko/lib.h"
#include "pipeline.h"   // parse_field_type
#include "check.h"

using namespace metriko;

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: test_remesh <gridscale> <nvec> <mesh.obj> [more.obj ...]\n"; return 2; }
    const double scale = std::stod(argv[1]);
    const FieldType ft = parse_field_type(argv[2]);

    for (int a = 3; a < argc; ++a) {
        const char* mesh = argv[a];
        MatXd V; MatXi F;
        igl::readOBJ(mesh, V, F);
        CHECK(V.rows() > 0 && F.rows() > 0);

        try {
            const auto im = compute_quadrangulation_impl(V, F, scale, ft == FieldType::CurvatureAligned);
            CHECK(im.hmesh && im.mgrph && im.emesh);
            CHECK(validate_no_crossing(*im.emesh, "test") == 0);
            CHECK(!im.q_faces.empty());
            CHECK(rg::all_of(im.q_ports, [](const qex::Qport& p) { return p.isConnected; }));
            for (const auto& [qhalfs]: im.q_faces) CHECK(qhalfs.size() == 4);
            CHECK(im.q_idx.rows() == (int)im.q_faces.size() && im.q_val.rows() == im.q_idx.rows());
            CHECK((im.q_val.col(0).array() >= 0).all());
            std::cout << "[test_remesh] OK  " << mesh << "  quads=" << im.q_idx.rows() << std::endl;
        } catch (const std::exception& ex) {
            std::cerr << "FAIL: compute_remesh threw: " << ex.what() << "  (" << mesh << ")\n";
            return 1;
        }
    }
    return 0;
}
