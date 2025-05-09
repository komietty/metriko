//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_VISUALIZE_TEDGE_EBD_H
#define METRIKO_VISUALIZE_TEDGE_EBD_H

namespace metriko::visualizer {

    inline void visualize_embedding(
        const Hmesh&  hmesh,
        const std::vector<EmbeddedTEdge>& embedded_t_edges,
        const VecXd& X
    ) {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> x;
        std::vector<double> x_sum;
        std::vector<double> teids;

        for (auto& ete : embedded_t_edges) {
            for(int i = 0; i < ete.vids.size() - 1; i++) {
                int vid0 = ete.vids[i];
                int vid1 = ete.vids[i + 1];
                Row3d p0 = hmesh.verts[vid0].pos();
                Row3d p1 = hmesh.verts[vid1].pos();
                ns.emplace_back(p0.x(), p0.y(), p0.z());
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                es.emplace_back(std::array{ns.size() - 2, ns.size() - 1});
                x.emplace_back(ete.vals[i + 1]);
                x_sum.emplace_back(X[ete.teid]);
                teids.emplace_back(ete.teid);
            }
            //for(Half h: ete.halfs) {
            //    Row3d p0 = h.tail().pos();
            //    Row3d p1 = h.head().pos();
            //    ns.emplace_back(p0.x(), p0.y(), p0.z());
            //    ns.emplace_back(p1.x(), p1.y(), p1.z());
            //    es.emplace_back(std::array{ns.size() - 2, ns.size() - 1});
            //    //x.emplace_back(ete.vals[i + 1]);
            //    x_sum.emplace_back(X[ete.teid]);
            //    teids.emplace_back(ete.teid);
            //}
        }

        auto c = polyscope::registerCurveNetwork("embeddings", ns, es);
        c->addEdgeScalarQuantity("x", x)->setEnabled(true);
        c->addEdgeScalarQuantity("x_sum", x_sum)->setEnabled(false);
        c->addEdgeScalarQuantity("teid", teids)->setEnabled(false);
        c->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.0005);
        c->setMaterial("flat");
    }
}

#endif
