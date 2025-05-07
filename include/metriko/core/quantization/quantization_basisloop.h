//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_BASISLOOP_H
#define METRIKO_QUANTIZATION_BASISLOOP_H
#include "metriko/core/tmesh/tmesh.h"

namespace metriko {
    class Comparator {
    public:
        int length;
        double weight;
        int thid;

        bool operator==(const Comparator& rhs) const { return   length == rhs.length && weight == rhs.weight;  }
        bool operator!=(const Comparator& rhs) const { return !(length == rhs.length && weight == rhs.weight); }
        bool operator< (const Comparator& rhs) const { return length < rhs.length || (length == rhs.length && weight > rhs.weight); } // flip weight comparator
        bool operator> (const Comparator& rhs) const { return length > rhs.length || (length == rhs.length && weight < rhs.weight); } // same
        bool operator<=(const Comparator& rhs) const { return !(*this > rhs); }
        bool operator>=(const Comparator& rhs) const { return !(*this < rhs); }
    };

    template<typename Func>
    std::vector<int> gen_basis_loop(
        const std::vector<Tquad>& tquads,
        const std::vector<Thalf>& thalfs,
        const VecXi& thalf_tquad_table,
        const VecXi& thalf_sides_table,
        const VecXd& R,
        const Thalf& bgn,
        Func compare
    ) {
        std::priority_queue<Comparator, std::vector<Comparator>, decltype(compare)> q(compare);
        std::unordered_map<int, int> m;
        std::vector<int> visited;
        std::vector<std::vector<int>> result;

        int counter = 0;
        q.emplace(Comparator{0, 0, bgn.id});

        do {
            int len = q.top().length;
            Thalf prev = thalfs[q.top().thid];
            q.pop();

            int tquad_idx = thalf_tquad_table[prev.id];
            int sides_idx = thalf_sides_table[prev.id];
            auto& tq = tquads[tquad_idx];

            std::vector<int> pair_idcs;
            for(int i = 0; i < tq.thids.size(); i++) {
                int thid = tq.thids[i];
                int side = tq.sides[i];
                if (side == (sides_idx + 2) % 4) {
                    pair_idcs.emplace_back(thid);
                }
            }

            for (int pair_idx: pair_idcs) {
                const auto &pair = thalfs[pair_idx];
                const auto &curr = pair.twin();

                if (curr.id == bgn.id && counter > 0) {
                    int idx = prev.id;
                    std::vector<int> res;
                    res.emplace_back(idx);
                    do {
                        idx = m[idx];
                        res.emplace_back(idx);
                    } while (idx != bgn.id);
                    return res;
                }

                bool visited_ = false;
                int tmp_idx = prev.id;
                while (tmp_idx != bgn.id) {
                    if (curr.id == tmp_idx) { visited_ = true; break; }
                    tmp_idx = m[tmp_idx];
                }
                if (!visited_ && rg::none_of(visited, [&](int id) { return id == curr.id; })) {
                    q.emplace(Comparator{len - 1, R[curr.edge().id], curr.id});
                    m.emplace(curr.id, prev.id);
                    visited.emplace_back(curr.id);
                }
            }
            counter++;
        } while (!q.empty());
    }
}

#endif
