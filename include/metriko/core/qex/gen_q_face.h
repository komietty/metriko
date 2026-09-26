#ifndef METRIKO_QEX_GEN_Q_FACE_H
#define METRIKO_QEX_GEN_Q_FACE_H
#include "common.h"

namespace metriko::qex {

    inline vec<Qface> generate_q_faces(
        const vec<Qport> &qps,
        const vec<Qedge> &qes
    ) {
        vec<Qface> qfs;
        vec<Qhalf> qhs;
        for (int i = 0; i < qes.size(); i++) {
            auto &qe = qes[i];
            qhs.emplace_back(qe, i * 2 + 0, true);
            qhs.emplace_back(qe, i * 2 + 1, false);
        }

        vec po2qh(qps.size(), -1);
        vec visit(qhs.size(), false);

        for (const auto &qh: qhs) po2qh[qh.port1().idx] = qh.idx;

        for (size_t i = 0; i < qhs.size(); ++i) {
            if (visit[i]) continue;
            vec qfhs{qhs[i]};
            visit[i] = true;
            auto *qp = &qhs[i].port2();
            for (int i = 0; i < 3; i++) {
                int h = po2qh[qps[qp->prev_id].idx]; if (h < 0) break;
                qfhs.emplace_back(qhs[h]);
                visit[h] = true;
                qp = &qhs[h].port2();
            }
            if (qfhs.size() == 4) {
                if (qps[qp->prev_id].idx != qhs[i].port1().idx) throw std::runtime_error("not closed ring");
                qfs.emplace_back(qfhs);
            }
        }
        return qfs;
    }
}

#endif
