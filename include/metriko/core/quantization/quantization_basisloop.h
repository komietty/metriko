//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_BASISLOOP_H
#define METRIKO_QUANTIZATION_BASISLOOP_H
#include <queue>
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
    vec<int> gen_basis_loop(
        const vec<Tquad>& tquads,
        const vec<Thalf>& thalfs,
        const VecXi& thalf_tquad_table,
        const VecXi& thalf_sides_table,
        const VecXd& R,
        const Thalf& bgn,
        Func compare
    ) {
        std::priority_queue<Comparator, vec<Comparator>, decltype(compare)> q(compare);
        vec parents(thalfs.size(), -1);    // parent thalf in the search tree
        vec visited(thalfs.size(), false); // pushed once already

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
                    do { idx = parents[idx]; res.emplace_back(idx); }
                    while (idx != bgn.id);
                    return res;
                }

                // every ancestor of prev except bgn has been pushed, and bgn is
                // handled above, so the pushed flag covers the ancestor test too
                if (!visited[curr.id]) {
                    q.emplace(Comparator{.length = len - 1, .weight = R[curr.edge().id], .thid = curr.id});
                    parents[curr.id] = prev.id;
                    visited[curr.id] = true;
                }
            }
            counter++;
        } while (!q.empty());
        throw std::runtime_error("cannot find a loop");
    }
}
#endif
