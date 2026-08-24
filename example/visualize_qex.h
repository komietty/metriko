#ifndef TMESH_MUT_COLLAPSE_TQUAD_CPP_VISUALIZE_QEX_H
#define TMESH_MUT_COLLAPSE_TQUAD_CPP_VISUALIZE_QEX_H
#include "./visualize_common.h"
#include "metriko/core/qex/sanitization.h"
#include "metriko/core/qex/gen_q_vert.h"
#include "metriko/core/qex/gen_q_port.h"
#include "metriko/core/qex/gen_q_edge.h"
#include "metriko/core/qex/gen_q_face.h"
#include "metriko/core/qex/refinement.h"

namespace metriko::visualizer {

inline void visualize_qverts(
    const vec<qex::Qvert>& vqvs,
    const vec<qex::Qvert>& eqvs,
    const vec<qex::Qvert>& fqvs,
    const double scale = 0.001,
    const bool show = true
) {
    vec<glm::vec3> VQV;
    vec<glm::vec3> EQV;
    vec<glm::vec3> FQV;
    for (const auto &q: vqvs) VQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z());
    for (const auto &q: eqvs) EQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z());
    for (const auto &q: fqvs) FQV.emplace_back(q.pos.x(), q.pos.y(), q.pos.z());
    auto vq = polyscope::registerPointCloud("VQV", VQV);
    auto eq = polyscope::registerPointCloud("EQV", EQV);
    auto fq = polyscope::registerPointCloud("FQV", FQV);
    vq->setEnabled(show);
    eq->setEnabled(show);
    fq->setEnabled(show);
    vq->setPointRadius(scale);
    eq->setPointRadius(scale);
    fq->setPointRadius(scale);
    vq->resetTransform();
    eq->resetTransform();
    fq->resetTransform();
}

inline void visualize_qports(
    const Hmesh& mesh,
    const VecXc& uv,
    const vec<qex::Qport>& q_ports,
    const double scale = 0.001,
    const bool show = true
) {
    std::vector<glm::vec3> QP;
    std::vector<int> QP_idx, QP_fid, QP_dir, QP_n, QP_p;
    std::vector<int> QP_sid_vert, QP_sid_edge, QP_sid_face;
    std::vector<double> QP_u, QP_v;
    std::vector<double> QP_flag(q_ports.size(), 0);
    for (const auto& q: q_ports) {
        auto f = mesh.faces[q.fid];
        auto p = conversion_2d_3d(f, uv, q.uv + q.dir * 0.15);
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
    qp->setEnabled(show);
    qp->resetTransform();
    qp->setPointRadius(scale);
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
}

// show one qvert's port group: coincident points carrying the port directions,
// the cycle order, and outlines of the faces the ports live in. every structure
// is named with ids so offenders reported in the log can be inspected one by one
inline void visualize_qport_group(
    const Hmesh& hm,
    const VecXc& cfn,
    const vec<qex::Qport>& qports,
    int pid   // any port of the group (the group = contiguous ports at one position)
) {
    int first = pid, last = pid;
    auto same = [&](int a, int b) { return (qports[a].pos - qports[b].pos).norm() < 1e-12; };
    while (first > 0 && same(first - 1, pid)) --first;
    while (last + 1 < (int)qports.size() && same(last + 1, pid)) ++last;

    std::vector<glm::vec3> ps, ds;
    std::vector<double> order, idxs, fids;
    for (int k = first; k <= last; ++k) {
        const auto& p = qports[k];
        ps.emplace_back(p.pos.x(), p.pos.y(), p.pos.z());
        Row3d d = (conversion_2d_3d(hm.faces[p.fid], cfn, p.uv + p.dir)
                 - conversion_2d_3d(hm.faces[p.fid], cfn, p.uv)).normalized();
        ds.emplace_back(d.x(), d.y(), d.z());
        order.push_back(k - first);
        idxs.push_back(p.idx);
        fids.push_back(p.fid);
    }
    const auto& p0 = qports[first];
    auto* pc = polyscope::registerPointCloud(
        std::format("bad qvert p{} (vid {} eid {})", p0.idx, p0.vid, p0.eid), ps);
    pc->addVectorQuantity("dir", ds)->setEnabled(true);
    pc->addScalarQuantity("cycle order", order)->setEnabled(true);
    pc->addScalarQuantity("port idx", idxs);
    pc->addScalarQuantity("fid", fids);
    pc->setPointRadius(0.002);
    pc->resetTransform();

    std::set<double> fset(fids.begin(), fids.end());
    for (double fd: fset) {
        int fid = (int)fd;
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t c = 0;
        for (Half h: hm.faces[fid].adjHalfs()) {
            Row3d a = h.tail().pos();
            Row3d b = h.head().pos();
            ns.emplace_back(a.x(), a.y(), a.z());
            ns.emplace_back(b.x(), b.y(), b.z());
            es.push_back({c, c + 1});
            c += 2;
        }
        auto* cn = polyscope::registerCurveNetwork(std::format("bad face {}", fid), ns, es);
        cn->setMaterial("flat");
        cn->setColor({1., 0.15, 0.1});
        cn->setRadius(0.0012);
        cn->resetTransform();
    }
}

inline void visualize_qedges(
    const vec<qex::Qedge>& qedges,
    const double scale = 0.001
) {
    std::vector<std::array<size_t, 2>> QE;
    std::vector<glm::vec3> QN;
    std::vector<double> p1;
    std::vector<double> p2;
    size_t counter = 0;
    for (auto& q: qedges) {
        QN.emplace_back(q.port1.pos.x(), q.port1.pos.y(), q.port1.pos.z());
        QN.emplace_back(q.port2.pos.x(), q.port2.pos.y(), q.port2.pos.z());
        QE.emplace_back(std::array{counter, counter + 1});
        p1.emplace_back(q.port1.idx);
        p2.emplace_back(q.port2.idx);
        counter += 2;
    }
    auto q_edge_curv = polyscope::registerCurveNetwork("q_edges", QN, QE);
    q_edge_curv->setMaterial("flat");
    q_edge_curv->resetTransform();
    q_edge_curv->setRadius(scale);
    q_edge_curv->setEnabled(false);
    q_edge_curv->addEdgeScalarQuantity("p1 idx", p1);
    q_edge_curv->addEdgeScalarQuantity("p2 idx", p2);
}

inline void visualize_qfaces(
    const Hmesh& hm,
    const vec<qex::Qface>& qfaces,
    const bool refine = true
) {
    std::vector<std::array<size_t, 4> > QF;
    int l = qfaces.size();
    MatXd pos(l * 4, 3);
    MatXi idx(l, 4);
    for (int i = 0; i < l; i++) {
    for (int j = 0; j < 4; j++) {
        pos.row(i * 4 + j) = qfaces[i].qhalfs[j].port1().pos;
        idx(i, j) = i * 4 + j;
    }}

    MatXd pos_refined;
    MatXi idx_refined;
    if (refine) {
        qex::refinement_hmesh(pos, idx, hm.pos, hm.idx, pos_refined, idx_refined);
    } else {
        pos_refined = pos;
        idx_refined = idx;
    }
    auto* surf = polyscope::registerSurfaceMesh("quad mesh", pos_refined, idx_refined);
    surf->setShadeStyle(polyscope::MeshShadeStyle::Flat);
    surf->setEdgeWidth(1.);
}

}

#endif
