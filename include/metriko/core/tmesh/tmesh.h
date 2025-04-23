//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_TMESH_H
#define METRIKO_TMESH_H
#include "tquad.h"

namespace metriko {
    class Tmesh {
    public:
        std::vector<Tquad> tquads;
        std::vector<Thalf> thalfs;
        std::vector<Tedge> tedges;
        std::vector<int> stemming_thalfs; // todo: remove...
        VecXi th2sing; // -1 if thalf is not from singular, otherwise vertex id
        VecXi th2quad;
        VecXi th2side;
        VecXi th2iter;
        size_t nTQ;
        size_t nTE;
        size_t nTH;

        explicit Tmesh() = default;
        explicit Tmesh(const std::vector<Mcurv>& splines) {
            for (auto &me: splines) {
                Msgmt start = me.sgmts.front();
                for (auto it = me.sgmts.begin(); it != me.sgmts.end(); ++it) {
                    if (it->to.side != None) {
                        int s = tedges.size();
                        tedges.emplace_back(this, s, start, *it);
                        thalfs.emplace_back(this, s, s * 2 + 0, s * 2 + 1, true);
                        thalfs.emplace_back(this, s, s * 2 + 1, s * 2 + 0, false);
                        if (it->next_id != -1) start = me.sgmts[it->next_id];
                    }
                }
            }

            std::vector visit(thalfs.size(), false);

            while (rg::any_of(visit, [](const bool f) { return !f; })) {
                auto it = rg::find(visit, false);
                auto id = std::distance(visit.begin(), it);
                auto tq = Tquad(splines, thalfs, thalfs[id]);
                tquads.emplace_back(tq);
                for (int thid: tq.thids) visit[thid] = true;
            }

            nTE = tedges.size();
            nTH = thalfs.size();
            nTQ = tquads.size();

            th2side.resize(nTH);
            th2quad.resize(nTH);
            th2iter.resize(nTH);
            th2sing.resize(nTH);

            for (int i = 0; i < nTQ; i++) {
                const Tquad &tq = tquads[i];
                for (int j = 0; j < tq.thids.size(); j++) {
                    th2quad[tq.thids[j]] = i;
                    th2side[tq.thids[j]] = tq.sides[j];
                    th2iter[tq.thids[j]] = j;
                }
            }

            th2sing.setConstant(-1);
            for (auto &me: splines) {
                auto it = rg::find_if(thalfs, [&](auto &th) { return th.edge().seg_fr == me.sgmts.front(); });
                assert(it != thalfs.end());
                stemming_thalfs.emplace_back(it->id);
                th2sing[it->id] = me.port.vert.id;
            }
        }

        // should be in tquad?
        int next_thid(const int curr) const {
            int iQ = th2quad[curr];
            int iT = th2iter[curr];
            auto& tq = tquads[iQ];
            return tq.thids[(iT + 1) % tq.thids.size()];
        }

        // should be in tquad?
        int prev_thid(const int curr) const {
            int iQ = th2quad[curr];
            int iT = th2iter[curr];
            auto& tq = tquads[iQ];
            return tq.thids[(iT - 1 + tq.thids.size()) % tq.thids.size()];
        }
    };
}

#include "tmesh.ipp"
#endif
