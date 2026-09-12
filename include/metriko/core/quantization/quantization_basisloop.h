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
        std::priority_queue<Comparator, vec<Comparator>, decltype(compare)> q(compare);
        std::unordered_map<int, int> m;
        vec<int> visited;
        vec<vec<int>> result;

        int counter = 0;
        q.emplace(Comparator{.length = 0, .weight = 0, .thid = bgn.id});

        do {
            int len = q.top().length;
            Thalf prev = thalfs[q.top().thid];
            q.pop();

            int tqid = thalf_tquad_table[prev.id];
            int side = thalf_sides_table[prev.id];
            for (int thid: tquads[tqid].thids_by_side((side + 2) % 4)) {
                const auto &pair = thalfs[thid];
                const auto &curr = pair.twin();

                if (curr.id == bgn.id && counter > 0) {
                    int idx = prev.id;
                    vec<int> res;
                    res.emplace_back(idx);
                    do { idx = m[idx]; res.emplace_back(idx); }
                    while (idx != bgn.id);
                    return res;
                }

                bool visited_ = false;
                int tmp_idx = prev.id;
                while (tmp_idx != bgn.id) {
                    if (curr.id == tmp_idx) { visited_ = true; break; }
                    tmp_idx = m[tmp_idx];
                }
                if (!visited_ && rg::none_of(visited, [&](int id) { return id == curr.id; })) {
                    q.emplace(Comparator{.length = len - 1, .weight = R[curr.edge().id], .thid = curr.id});
                    m.emplace(curr.id, prev.id);
                    visited.emplace_back(curr.id);
                }
            }
            counter++;
        } while (!q.empty());
        throw std::runtime_error("cannot find a loop");
    }
}

#endif
