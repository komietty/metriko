//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QUANTIZATION_BASISLOOP_H
#define METRIKO_QUANTIZATION_BASISLOOP_H
#include <queue>
#include "metriko/tmesh/emesh.h"

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

    // side of every thalf within its tquad. emesh keeps it only inside the tquad data, so the quantization,
    // which looks it up for every step of the search, builds the table once
    inline VecXi thalf_sides(const Emesh& tm) {
        VecXi sides = VecXi::Constant(tm.thalfs.size(), -1);
        for (const Equad& tq: tm.live_tquads())
            for (const Edata& d: tq.data) sides[d.thid] = d.side;
        return sides;
    }

    template<typename Func>
    vec<int> gen_basis_loop(
        const Emesh& tm,
        const VecXi& th2side,
        const VecXd& R,
        const Ehalf& bgn,
        Func compare
    ) {
        const auto& thalfs = tm.thalfs;
        std::priority_queue<Comparator, vec<Comparator>, decltype(compare)> q(compare);
        vec parents(thalfs.size(), -1);    // parent thalf in the search tree
        vec visited(thalfs.size(), false); // pushed once already

        int counter = 0;
        q.emplace(Comparator{.length = 0, .weight = 0, .thid = bgn.id});

        do {
            int len = q.top().length;
            const Ehalf& prev = thalfs[q.top().thid];
            q.pop();

            for (int thid: tm.tquads[prev.tqid].thids((th2side[prev.id] + 2) % 4)) {
                const auto &curr = thalfs[thalfs[thid].twid];

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
                    q.emplace(Comparator{.length = len - 1, .weight = R[curr.teid], .thid = curr.id});
                    parents[curr.id] = prev.id;
                    visited[curr.id] = true;
                }
            }
            counter++;
        } while (!q.empty());
        METRIKO_FAIL("no loop through thalf {}", bgn.id);
    }
}
#endif
