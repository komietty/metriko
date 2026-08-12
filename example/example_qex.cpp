#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include <igl/readOBJ.h>
#include "visualize_qex.h"
#include "visualize_hmesh.h"
#include "metriko/core/qex/sanitization.h"
#include "metriko/core/qex/gen_q_vert.h"
#include "metriko/core/qex/gen_q_port.h"
#include "metriko/core/qex/gen_q_edge.h"
#include "metriko/core/qex/gen_q_face.h"
using namespace metriko;

// try testing horsers.obj with a grid scale of 0.005
int main(int argc, char **argv) {
    MatXd V;
    MatXi F;
    int N = 4;
    igl::readOBJ(argv[1], V, F);
    auto mesh = std::make_unique<Hmesh>(V, F);
    auto rawf = std::make_unique<FaceRosyField>(*mesh, N, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Curl);
    auto seam = compute_seam(*rawf);
    auto cutm = compute_cut_mesh(*mesh, seam);
    auto cmbf = compute_combbed_field(*rawf, seam);
    auto extf = compute_extrinsic_field(*mesh, *cmbf, N);

    RosyParameterization rp(*mesh, *cutm, extf, cmbf->singular, cmbf->matching, seam, N, std::stod(argv[2]));
    rp.seamless = true;
    rp.localInjectivity = true;
    rp.roundSeams = false;
    rp.verbose = false;
    rp.setup();
    rp.integ();
    MatXd uv1(mesh->nF * 3, 2);
    VecXc uv2(mesh->nF * 3);
    for (const Face f: mesh->faces) {
        uv1.row(f.id * 3 + 0) << rp.cfn(f.id, 0), rp.cfn(f.id, 1);
        uv1.row(f.id * 3 + 1) << rp.cfn(f.id, 4), rp.cfn(f.id, 5);
        uv1.row(f.id * 3 + 2) << rp.cfn(f.id, 8), rp.cfn(f.id, 9);
    }

    for (const Face f: mesh->faces) {
        complex a = complex{uv1(f.id * 3 + 0, 0), uv1(f.id * 3 + 0, 1)};
        complex b = complex{uv1(f.id * 3 + 1, 0), uv1(f.id * 3 + 1, 1)};
        complex c = complex{uv1(f.id * 3 + 2, 0), uv1(f.id * 3 + 2, 1)};
        assert(orientation(a, b, c) > 0);
        uv2(f.id * 3 + 0) = a;
        uv2(f.id * 3 + 1) = b;
        uv2(f.id * 3 + 2) = c;
    }

    qex::sanitization(*mesh, cmbf->matching, cmbf->singular, 4, uv2);

    for (const Face f: mesh->faces) {
        assert(orientation(
            uv2(f.id * 3 + 0),
            uv2(f.id * 3 + 1),
            uv2(f.id * 3 + 2)) > 0);
    }

    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

    visualizer::visualize_mesh_with_uv(mesh->pos, mesh->idx, uv1, false);
    visualizer::visualize_seam(*mesh, seam);

    std::vector<qex::Qvert> vqvs;
    std::vector<qex::Qvert> eqvs;
    std::vector<qex::Qvert> fqvs;
    qex::generate_q_vert(*mesh, uv2, vqvs, eqvs, fqvs);

    visualizer::visualize_qverts(vqvs, eqvs, fqvs);

    std::vector<qex::Qport> q_ports;
    qex::generate_vqvert_qport(*mesh, uv2, vqvs, q_ports);
    qex::generate_eqvert_qport(*mesh, uv2, eqvs, q_ports);
    qex::generate_fqvert_qport(*mesh, fqvs, q_ports);

    visualizer::visualize_qports(*mesh, uv2, q_ports);

    auto qedges = qex::generate_q_edge(*mesh, uv2, cmbf->matching, q_ports);
    auto qfaces = qex::generate_q_faces(q_ports, qedges);

    visualizer::visualize_qedges(qedges);
    visualizer::visualize_qfaces(*mesh, qfaces);

    polyscope::show();
    return 0;
}
