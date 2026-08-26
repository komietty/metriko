#ifndef METRIKO_QEX_GEN_Q_FACE_H
#define METRIKO_QEX_GEN_Q_FACE_H
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
                if (it2 == qhs.end()) {
                    // TODO TEMP: report and skip instead of asserting; the
                    // partial face is dropped by the size check below
                    const auto& p = qps[qp.prev_id];
                    std::println("[qface] no qhalf starts at port {} (vid {}, eid {}, fid {}, pos {} {} {})",
                                 p.idx, p.vid, p.eid, p.fid, p.pos.x(), p.pos.y(), p.pos.z());
                    break;
                }
                qfhs.emplace_back(*it2);
                visit[it2->idx] = true;
                qp = it2->port2();
            }
            if (qfhs.size() == 4) {
                // the walk must return to the start port; otherwise some qvert's
                // port cycle (next/prev) is inconsistent
                if (qps[qp.prev_id].idx != it1->port1().idx)
                    std::println("[qface] non-closing walk at port {} (vid {}, eid {}, fid {}, pos {} {} {})", qp.idx, qp.vid, qp.eid, qp.fid, qp.pos.x(), qp.pos.y(), qp.pos.z());
                qfs.emplace_back(qfhs);
            }
        }
        return qfs;
    }
}

#endif
