#ifndef GEN_Q_FACE_H
#define GEN_Q_FACE_H
#include "common.h"

namespace metriko::qex {
    inline Qface generate_q_face(
        const std::vector<Qport> &qps,
        const std::vector<Qhalf> &qhs,
        const Qhalf &qh
    ) {
        Qface qf;
        qf.qhalfs.emplace_back(qh);
        Qport curr = qh.port2();

        for (int i = 0; i < 3; i++) {
            auto it = rg::find_if(qhs, [&](const Qhalf &qh_) { return qh_.port1().idx == qps[curr.prev_id].idx; });
            assert(it != qhs.end());
            qf.qhalfs.emplace_back(*it);
            curr = it->port2();
        }
        return qf;
    }

    inline std::vector<Qface> generate_q_faces(
        const std::vector<Qport> &qps,
        const std::vector<Qedge> &qes
    ) {
        std::vector<Qface> qfs;
        std::vector<Qhalf> qhs;
        for (int i = 0; i < qes.size(); i++) {
            const Qedge &qe = qes[i];
            qhs.emplace_back(qe, i * 2 + 0, true);
            qhs.emplace_back(qe, i * 2 + 1, false);
        }

        std::vector visit(qhs.size(), false);
        while (true) {
            auto it = rg::find_if(qhs, [&](const Qhalf &qh) { return !visit[qh.idx]; });
            if (it == qhs.end()) break;
            Qface qf = generate_q_face(qps, qhs, *it);
            qfs.emplace_back(qf);
            for (const Qhalf &qh_: qf.qhalfs) visit[qh_.idx] = true;
        }
        return qfs;
    }
}

#endif
