#ifndef GEN_Q_FACE_H
#define GEN_Q_FACE_H
#include "common.h"

namespace metriko::qex {

    inline Qface generate_q_face(
        const std::vector<Qport> &q_ports,
        std::vector<Qhalf> &q_halfs,
        const Qhalf& qh
    ) {
        Qface qf;
        //std::cout << "====== new qf " <<std::endl;

        qf.qhalfs.emplace_back(qh);
        //std::cout << "qh ports: " << qh.port1().idx << ", " << qh.port2().idx << std::endl;
        Qport qp_curr = qh.port2();
        //std::cout << "qp_curr.idx: " << qp_curr.idx <<std::endl;

        for (int i = 0; i < 3; i++) {
            auto it = std::ranges::find_if(q_halfs, [&](const Qhalf& qh_) {
                return qh_.port1().idx == q_ports[qp_curr.prev_id].idx;
            });
            if (it == q_halfs.end()) continue;
            assert(it != q_halfs.end());
            qf.qhalfs.emplace_back(*it);
            //std::cout << "qh ports: " << it->port1().idx << ", " << it->port2().idx << std::endl;
            qp_curr = it->port2();
            //std::cout << "qp_curr.idx: " << qp_curr.idx <<std::endl;
        }
        return qf;
    }

    inline std::vector<Qface> generate_q_faces(
        const Hmesh &mesh,
        const std::vector<Qport> &q_ports,
        const std::vector<Qedge> &q_edges
    ) {
        std::vector<Qface> q_faces;
        std::vector<Qhalf> qhalfs;
        for (int i = 0; i < q_edges.size(); i++) {
            const Qedge& qe = q_edges[i];
            qhalfs.emplace_back(qe, i * 2 + 0, true);
            qhalfs.emplace_back(qe, i * 2 + 1, false);
        }


        std::vector visit(qhalfs.size(), false);
        int counter1 = 0;
        while (true) {
            auto it = std::ranges::find_if(qhalfs, [&](const Qhalf& qh){ return !visit[qh.idx];});
            if (it == qhalfs.end()) break;
            const Qface qf = generate_q_face(q_ports, qhalfs, *it);
            q_faces.emplace_back(qf);
            for (const Qhalf& qh_: qf.qhalfs) { visit[qh_.idx] = true; }
            counter1++;
        }
        int l = q_faces.size();

        std::vector<std::array<size_t, 2>> QE;
        std::vector<glm::vec3> QN;
        size_t counter = 0;
        MatXd pos(l * 4, 3);
        MatXi idx(l, 4);

        //for (const auto&[qhalfs]: q_faces) {
        for (int i = 0; i < l; i++) {
            const Qface& qf = q_faces[i];
            for (const Qhalf& qh: qf.qhalfs) {
                qex::Qport port1 = qh.port1();
                qex::Qport port2 = qh.port2();
                QN.emplace_back(port1.pos.x(), port1.pos.y(), port1.pos.z());
                QN.emplace_back(port2.pos.x(), port2.pos.y(), port2.pos.z());
                QE.emplace_back(std::array{counter, counter + 1});
                counter += 2;
            }
        }

        for (int i = 0; i < l; i++) {
            for (int j = 0; j < 4; j++) {
                pos.row(i * 4 + j) = q_faces[i].qhalfs[j].port1().pos;
                idx(i, j) = i * 4 + j;
            }
        }
        auto surf = polyscope::registerSurfaceMesh("quad mesh!", pos, idx);
        surf->setShadeStyle(polyscope::MeshShadeStyle::Flat);
        surf->setEdgeWidth(1.);

        auto q_edge_curv = polyscope::registerCurveNetwork("q edge test", QN, QE);
        q_edge_curv->resetTransform();
        q_edge_curv->setRadius(0.005);
        q_edge_curv->setEnabled(true);

        return q_faces;
    }
}

#endif
