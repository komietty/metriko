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

/*
using json = nlohmann::json;

void toJson(const std::string& path, const VecXc& cfn) {
    json j;
    std::vector<double> data;
    for (complex c: cfn) {
        data.push_back(c.real());
        data.push_back(c.imag());
    }
    j["cfn"] = data;
    std::ofstream o(path.c_str());
    o << j.dump() << std::endl;
    o.close();
}

void fromJson(const std::string& path, MatXd& uv1, VecXc& uv2) {
    if (std::ifstream ifs(path); ifs.good()) {
        json j;
        ifs >> j;
        std::vector<double> data = j["cfn"];
        uv1.resize(data.size() / 2, 2);
        uv2.resize(data.size() / 2);
        for (int i = 0; i < data.size(); i += 2) {
            uv1.row(i / 2) << data[i], data[i + 1];
            uv2(i / 2) = complex(data[i], data[i + 1]);
        }
    }
}
*/

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

    //*
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
        assert(qex::orientation(a, b, c) > 0);
        uv2(f.id * 3 + 0) = a;
        uv2(f.id * 3 + 1) = b;
        uv2(f.id * 3 + 2) = c;
    }

    qex::sanitization(*mesh, cmbf->matching, cmbf->singular, 4, uv2);

    for (const Face f: mesh->faces) {
        assert(qex::orientation(
            uv2(f.id * 3 + 0),
            uv2(f.id * 3 + 1),
            uv2(f.id * 3 + 2)) > 0);
    }

    //*/

    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;

    // visualize mesh
    {
        const auto surf = polyscope::registerSurfaceMesh("mesh", V, F);
        const auto prms = surf->addParameterizationQuantity("params", uv1);
        surf->addFaceVectorQuantity("cmb field", extf);
        surf->setEnabled(false);
        prms->setEnabled(true);
        prms->setStyle(polyscope::ParamVizStyle::GRID);
        prms->setGridColors(std::pair(glm::vec3(1, 1, 1), glm::vec3(0.2, 0.2, 0.2)));
        prms->setCheckerSize(1);
    }

    // visuailize seam
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2> > es;
        size_t i = 0;
        for (auto e: mesh->edges) {
            if (!seam[e.id]) continue;
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{i, i + 1});
            i += 2;
        }
        auto c = polyscope::registerCurveNetwork("seam", ns, es);
        c->setEnabled(false);
        c->resetTransform();
        c->setRadius(0.003);
    }

    // qex
    std::vector<qex::Qvert> vqvs;
    std::vector<qex::Qvert> eqvs;
    std::vector<qex::Qvert> fqvs;
    qex::generate_q_vert(*mesh, uv2, vqvs, eqvs, fqvs);

    // display qvert
    std::vector<glm::vec3> VQV;
    std::vector<glm::vec3> EQV;
    std::vector<glm::vec3> FQV;
    for (const auto &q: vqvs) { VQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z()); }
    for (const auto &q: eqvs) { EQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z()); }
    for (const auto &q: fqvs) { FQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z()); }
    auto vq = polyscope::registerPointCloud("VQV", VQV);
    auto eq = polyscope::registerPointCloud("EQV", EQV);
    auto fq = polyscope::registerPointCloud("FQV", FQV);
    vq->setEnabled(true);
    eq->setEnabled(false);
    fq->setEnabled(false);
    vq->setPointRadius(0.0005);
    eq->setPointRadius(0.0005);
    fq->setPointRadius(0.0005);
    vq->resetTransform();
    eq->resetTransform();
    fq->resetTransform();

    std::vector<qex::Qport> q_ports;
    qex::generate_vqvert_qport(*mesh, uv2, vqvs, q_ports);
    qex::generate_eqvert_qport(*mesh, uv2, eqvs, q_ports);
    qex::generate_fqvert_qport(*mesh, fqvs, q_ports);

    // display qport
    /*
    std::vector<glm::vec3> QP;
    std::vector<int> QP_idx, QP_fid, QP_dir, QP_n, QP_p;
    std::vector<int> QP_sid_vert, QP_sid_edge, QP_sid_face;
    std::vector<double> QP_u, QP_v;
    std::vector<double> QP_flag(q_ports.size(), 0);
    for (const auto& q: q_ports) {
        Face f = mesh->faces[q.fid];
        Row3d p = conversion_2d_3d(f, uv2, q.uv + q.dir * 0.15);
        QP.emplace_back(p.x(), p.y(), p.z());
        QP_idx.emplace_back(q.idx);
        QP_fid.emplace_back(f.id);
        QP_sid_vert.emplace_back(q.vid);
        QP_sid_edge.emplace_back(q.eid);
        QP_sid_face.emplace_back(q.fid);
        QP_u.emplace_back(q.uv.real());
        QP_v.emplace_back(q.uv.imag());
        QP_n.emplace_back(q.next_id);
        QP_p.emplace_back(q.prev_id);
        if      (equal(q.dir, complex(1, 0)))  QP_dir.emplace_back(0);
        else if (equal(q.dir, complex(0, 1)))  QP_dir.emplace_back(1);
        else if (equal(q.dir, complex(-1, 0))) QP_dir.emplace_back(2);
        else if (equal(q.dir, complex(0, -1))) QP_dir.emplace_back(3);
    }

    auto qp = polyscope::registerPointCloud("QP", QP);
    qp->setEnabled(false);
    qp->resetTransform();
    qp->setPointRadius(0.003);
    qp->addScalarQuantity("QP_idx", QP_idx);
    qp->addScalarQuantity("QP_fid", QP_fid);
    qp->addScalarQuantity("QP_dir", QP_dir);
    qp->addScalarQuantity("QP_sid_vert", QP_sid_vert);
    qp->addScalarQuantity("QP_sid_edge", QP_sid_edge);
    qp->addScalarQuantity("QP_sid_face", QP_sid_face);
    qp->addScalarQuantity("QP_u", QP_u);
    qp->addScalarQuantity("QP_v", QP_v);
    qp->addScalarQuantity("QP0_next", QP_n);
    qp->addScalarQuantity("QP0_prev", QP_p);
    qp->addScalarQuantity("QP_flag", QP_flag);
    //*/

    auto qedges = qex::generate_q_edge(*mesh, uv2, cmbf->matching, q_ports);
    auto qfaces = qex::generate_q_faces(q_ports, qedges);

    /*
    std::vector<std::array<size_t, 2>> QE;
    std::vector<glm::vec3> QN;
    std::vector<double> p1;
    std::vector<double> p2;
    size_t counter = 0;
    for (qex::Qedge& q: qedges) {
        QN.emplace_back(q.port1.pos.x(), q.port1.pos.y(), q.port1.pos.z());
        QN.emplace_back(q.port2.pos.x(), q.port2.pos.y(), q.port2.pos.z());
        QE.emplace_back(std::array{counter, counter + 1});
        p1.emplace_back(q.port1.idx);
        p2.emplace_back(q.port2.idx);
        counter += 2;
    }
    auto q_edge_curv = polyscope::registerCurveNetwork("q edge", QN, QE);
    q_edge_curv->resetTransform();
    q_edge_curv->setRadius(0.001);
    q_edge_curv->setEnabled(false);
    q_edge_curv->addEdgeScalarQuantity("p1 idx", p1);
    q_edge_curv->addEdgeScalarQuantity("p2 idx", p2);
    //*/

    //*
    {
        std::vector<std::array<size_t, 4> > QF;
        int l = (int) qfaces.size();
        MatXd pos(l * 4, 3);
        MatXi idx(l, 4);
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < 4; j++) {
                pos.row(i * 4 + j) = qfaces[i].qhalfs[j].port1().pos;
                idx(i, j) = i * 4 + j;
            }
        }
        auto surf = polyscope::registerSurfaceMesh("quad mesh!", pos, idx);
        surf->setShadeStyle(polyscope::MeshShadeStyle::Flat);
        surf->setEdgeWidth(1.);
    }
    //*/

    polyscope::show();
    return 0;
}
