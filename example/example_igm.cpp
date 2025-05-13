#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include <igl/readOBJ.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/qex/sanitization.h"
#include "metriko/core/qex/gen_q_vert.h"
#include "metriko/core/qex/gen_q_port.h"
#include "metriko/core/qex/gen_q_edge.h"
#include "metriko/core/qex/gen_q_face.h"

using namespace metriko;
int N = 4;
MatXd V;
MatXi F;

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    auto mesh = std::make_unique<Hmesh>(V, F);
    auto rawf = std::make_unique<FaceRosyField>(*mesh, N, FieldType::Smoothest);
    rawf->computeMatching(MatchingType::Curl);
    auto seam = compute_seam(*rawf);
    auto cutm = compute_cut_mesh(*mesh, seam);
    auto cmbf = compute_combbed_field(*rawf, seam);
    MatXd cmbExtRosy(mesh->nF, 3 * N);
    MatXd cmbExtZero(mesh->nF, 3);
    for (Face f: mesh->faces) {
        complex c0 = cmbf->field(f.id, 0);
        complex c1 = cmbf->field(f.id, 1);
        complex c2 = cmbf->field(f.id, 2);
        complex c3 = cmbf->field(f.id, 3);
        cmbExtZero.row(f.id) = (c0.real() * f.basisX() + c0.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 0, 1, 3) = (c0.real() * f.basisX() + c0.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 3, 1, 3) = (c1.real() * f.basisX() + c1.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 6, 1, 3) = (c2.real() * f.basisX() + c2.imag() * f.basisY()).normalized();
        cmbExtRosy.block(f.id, 9, 1, 3) = (c3.real() * f.basisX() + c3.imag() * f.basisY()).normalized();
    }

    RosyParameterization rp(*mesh, *cutm, cmbExtRosy, cmbf->singular, cmbf->matching, seam, N);
    rp.seamless = true;
    rp.localInjectivity = true;
    rp.roundSeams = false;
    rp.verbose = false;

    MatXd uv(mesh->nF * 3, 2);
    rp.setup();
    rp.integ();
    for (const Face f: mesh->faces) {
        uv.row(f.id * 3 + 0) << rp.cfn(f.id, 0), rp.cfn(f.id, 1);
        uv.row(f.id * 3 + 1) << rp.cfn(f.id, 4), rp.cfn(f.id, 5);
        uv.row(f.id * 3 + 2) << rp.cfn(f.id, 8), rp.cfn(f.id, 9);
    }

    /// ---- visualize mesh ---- ///
    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    const auto surf = polyscope::registerSurfaceMesh("mesh", V, F);
    const auto prms = surf->addParameterizationQuantity("params", uv);
    surf->addFaceVectorQuantity("cmb field", cmbExtZero);
    prms->setEnabled(true);
    prms->setStyle(polyscope::ParamVizStyle::GRID);
    prms->setCheckerSize(1);

    ///--- visuailize singularities ---///
    std::vector<glm::vec3> singPos;
    std::vector<double> singVal;
    for (Vert v: mesh->verts) {
        if (int s = cmbf->singular[v.id]; s != 0) {
            Row3d p = v.pos();
            singPos.emplace_back(p.x(), p.y(), p.z());
            singVal.emplace_back(s);
        }
    }
    auto pc = polyscope::registerPointCloud("face field singulars", singPos);
    pc->addScalarQuantity("face field singular nums", singVal);
    pc->setEnabled(true);
    pc->setPointRadius(0.002);
    pc->resetTransform();

    ///--- visuailize seam ---///
    std::vector<glm::vec3> nodesS;
    std::vector<std::array<size_t, 2>> edgesS;
    size_t counterS = 0;
    for (auto e: mesh->edges) {
        if (seam[e.id]) {
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            nodesS.emplace_back(p1.x(), p1.y(), p1.z());
            nodesS.emplace_back(p2.x(), p2.y(), p2.z());
            edgesS.emplace_back(std::array{counterS, counterS + 1});
            counterS += 2;
        }
    }
    auto curvS = polyscope::registerCurveNetwork("seam", nodesS, edgesS);
    curvS->resetTransform();
    curvS->setRadius(0.003);

    // qex
    qex::sanitization(*mesh, cmbf->matching, cmbf->singular, 4, uv);
    std::vector<qex::Qvert> vqverts;
    std::vector<qex::Qvert> eqverts;
    std::vector<qex::Qvert> fqverts;
    std::vector<qex::Qport> q_ports;
    std::vector<qex::Qedge> q_edges;

    qex::generate_q_vert(*mesh, uv, vqverts, eqverts, fqverts);
    qex::generate_vqvert_qport(*mesh, uv, vqverts, q_ports);
    qex::generate_eqvert_qport(*mesh, uv, eqverts, q_ports);
    qex::generate_fqvert_qport(*mesh, uv, fqverts, q_ports);

    //--- display q_vert ---//
    std::vector<glm::vec3> VQV;
    std::vector<glm::vec3> EQV;
    std::vector<glm::vec3> FQV;
    for (qex::Qvert q: vqverts) { VQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z()); }
    for (qex::Qvert q: eqverts) { EQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z()); }
    for (qex::Qvert q: fqverts) { FQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z()); }
    auto vq = polyscope::registerPointCloud("VQV", VQV);
    auto eq = polyscope::registerPointCloud("EQV", EQV);
    auto fq = polyscope::registerPointCloud("FQV", FQV);
    vq->setEnabled(true);
    eq->setEnabled(true);
    fq->setEnabled(true);
    vq->setPointRadius(0.003);
    eq->setPointRadius(0.003);
    fq->setPointRadius(0.003);
    vq->resetTransform();
    eq->resetTransform();
    fq->resetTransform();

    //--- display vert_q_port ---//
    std::vector<glm::vec3> QP0, QP1, QP2, QP3;
    std::vector<int> QP0_idx, QP1_idx, QP2_idx, QP3_idx;
    std::vector<int> QP0_next, QP1_next, QP2_next, QP3_next;
    std::vector<int> QP0_prev, QP1_prev, QP2_prev, QP3_prev;
    std::vector<int> QP0_fid, QP1_fid, QP2_fid, QP3_fid;
    std::vector<double> QP0_u, QP1_u, QP2_u, QP3_u;
    std::vector<double> QP0_v, QP1_v, QP2_v, QP3_v;

    double len = 0.1;
    for (qex::Qport& q: q_ports) {
        Face f = mesh->faces[q.fid];
        Row3d p = qex::uv2pos(f, uv, q.uvw + q.dir * len);
        if      (q.dir == Row2d(1, 0))  { QP0.emplace_back(p.x(), p.y(), p.z()); QP0_idx.emplace_back(q.idx); QP0_fid.emplace_back(f.id); QP0_u.emplace_back(q.uvw.x()); QP0_v.emplace_back(q.uvw.y()); }
        else if (q.dir == Row2d(0, 1))  { QP1.emplace_back(p.x(), p.y(), p.z()); QP1_idx.emplace_back(q.idx); QP1_fid.emplace_back(f.id); QP1_u.emplace_back(q.uvw.x()); QP1_v.emplace_back(q.uvw.y()); }
        else if (q.dir == Row2d(-1, 0)) { QP2.emplace_back(p.x(), p.y(), p.z()); QP2_idx.emplace_back(q.idx); QP2_fid.emplace_back(f.id); QP2_u.emplace_back(q.uvw.x()); QP2_v.emplace_back(q.uvw.y()); }
        else if (q.dir == Row2d(0, -1)) { QP3.emplace_back(p.x(), p.y(), p.z()); QP3_idx.emplace_back(q.idx); QP3_fid.emplace_back(f.id); QP3_u.emplace_back(q.uvw.x()); QP3_v.emplace_back(q.uvw.y()); }
        else { std::cout << q.dir << std::endl; }

        if      (q.dir == Row2d(1, 0))  { QP0_next.emplace_back(q.next_id); QP0_prev.emplace_back(q.prev_id); }
        else if (q.dir == Row2d(0, 1))  { QP1_next.emplace_back(q.next_id); QP1_prev.emplace_back(q.prev_id); }
        else if (q.dir == Row2d(-1, 0)) { QP2_next.emplace_back(q.next_id); QP2_prev.emplace_back(q.prev_id); }
        else if (q.dir == Row2d(0, -1)) { QP3_next.emplace_back(q.next_id); QP3_prev.emplace_back(q.prev_id); }
    }

    auto qp0 = polyscope::registerPointCloud("QP0", QP0);
    auto qp1 = polyscope::registerPointCloud("QP1", QP1);
    auto qp2 = polyscope::registerPointCloud("QP2", QP2);
    auto qp3 = polyscope::registerPointCloud("QP3", QP3);
    qp0->resetTransform(); qp0->setPointRadius(0.003);
    qp1->resetTransform(); qp1->setPointRadius(0.003);
    qp2->resetTransform(); qp2->setPointRadius(0.003);
    qp3->resetTransform(); qp3->setPointRadius(0.003);
    qp0->addScalarQuantity("QP0_idx", QP0_idx); qp0->addScalarQuantity("QP0_fid", QP0_fid); qp0->addScalarQuantity("QP0_u", QP0_u); qp0->addScalarQuantity("QP0_v", QP0_v);
    qp1->addScalarQuantity("QP1_idx", QP1_idx); qp1->addScalarQuantity("QP1_fid", QP1_fid); qp1->addScalarQuantity("QP1_u", QP1_u); qp1->addScalarQuantity("QP1_v", QP1_v);
    qp2->addScalarQuantity("QP2_idx", QP2_idx); qp2->addScalarQuantity("QP2_fid", QP2_fid); qp2->addScalarQuantity("QP2_u", QP2_u); qp2->addScalarQuantity("QP2_v", QP2_v);
    qp3->addScalarQuantity("QP3_idx", QP3_idx); qp3->addScalarQuantity("QP3_fid", QP3_fid); qp3->addScalarQuantity("QP3_u", QP3_u); qp3->addScalarQuantity("QP3_v", QP3_v);

    qp0->addScalarQuantity("QP0_next", QP0_next); qp0->addScalarQuantity("QP0_prev", QP0_prev);
    qp1->addScalarQuantity("QP1_next", QP1_next); qp1->addScalarQuantity("QP1_prev", QP1_prev);
    qp2->addScalarQuantity("QP2_next", QP2_next); qp2->addScalarQuantity("QP2_prev", QP2_prev);
    qp3->addScalarQuantity("QP3_next", QP3_next); qp3->addScalarQuantity("QP3_prev", QP3_prev);

    //*
    qex::generate_q_edge(*mesh, uv, cmbf->matching, q_ports, q_edges);
    qex::generate_q_faces(*mesh, q_ports, q_edges);

    std::vector<std::array<size_t, 2>> QE;
    std::vector<glm::vec3> QN;
    std::vector<double> p1;
    std::vector<double> p2;
    size_t counter = 0;
    for (qex::Qedge& q: q_edges) {
        QN.emplace_back(q.port1.pos.x(), q.port1.pos.y(), q.port1.pos.z());
        QN.emplace_back(q.port2.pos.x(), q.port2.pos.y(), q.port2.pos.z());
        QE.emplace_back(std::array{counter, counter + 1});
        p1.emplace_back(q.port1.idx);
        p2.emplace_back(q.port2.idx);
        counter += 2;
    }
    auto q_edge_curv = polyscope::registerCurveNetwork("q edge", QN, QE);
    q_edge_curv->resetTransform();
    q_edge_curv->setRadius(0.002);
    q_edge_curv->setEnabled(true);
    q_edge_curv->addEdgeScalarQuantity("p1 idx", p1);
    q_edge_curv->addEdgeScalarQuantity("p2 idx", p2);
    //*/

    polyscope::show();
    return 0;
}
