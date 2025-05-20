#ifndef GEN_Q_FACE_H
#define GEN_Q_FACE_H
#include "common.h"

namespace metriko::qex {

    inline std::vector<Qface> generate_q_faces(
        const std::vector<Qport> &qps,
        const std::vector<Qedge> &qes
    ) {
        std::vector<Qface> qfs;
        std::vector<Qhalf> qhs;
        for (int i = 0; i < qes.size(); i++) {
            auto &qe = qes[i];
            qhs.emplace_back(qe, i * 2 + 0, true);
            qhs.emplace_back(qe, i * 2 + 1, false);
        }

        std::vector visit(qhs.size(), false);
        while (true) {
            auto it1 = rg::find_if(qhs, [&](auto &qh) { return !visit[qh.idx]; });
            if (it1 == qhs.end()) break;
            std::vector qfhs{*it1};
            visit[it1->idx] = true;
            Qport qp = it1->port2();
            for (int i = 0; i < 3; i++) {
                auto it2 = rg::find_if(qhs, [&](auto &qh) {
                    return qh.port1().idx == qps[qp.prev_id].idx;
                });
                //assert(it2 != qhs.end());
                if (it2 != qhs.end()) {
                    qfhs.emplace_back(*it2);
                    visit[it2->idx] = true;
                    qp = it2->port2();
                }
            }
            if (qfhs.size() == 4) qfs.emplace_back(qfhs);
        }
        return qfs;
    }
}

#endif
