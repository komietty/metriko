//
// Created by saki on 2025/12/01.
//

#ifndef TMESH_H_TUTTE_CUTTING_UPSAMPLE_H
#define TMESH_H_TUTTE_CUTTING_UPSAMPLE_H
#include "tutte_cutting.h"

namespace metriko::tutte {
inline void compute_embedded_halfs(
    const Hmesh& hm,      // input hmesh
    const Tmesh& tm,      // input tmesh
    const VecXc& cf,      // input corner function of naive parameterization
    set<HalfData>& h_data //
) {
    vec<bool> visit = vec(hm.nV, false);

    int temp = 0;
    for (const Thalf& th: tm.thalfs) {
        if (!th.cano) continue;
        if (++temp > 15) {}
        const Tedge& te = th.edge();
        vec<int> vids0;
        vec<int> vids1;
        vec<Half> hs;
        for (const auto& s: te.segments()) {
            Face f = s.face;
            complex fr = s.fr.uv;
            complex to = s.to.uv;
            complex uv0 = cf[f.half().crnr().id];
            complex uv1 = cf[f.half().next().crnr().id];
            complex uv2 = cf[f.half().prev().crnr().id];
            Vert vt0 = f.half().next().head();
            Vert vt1 = f.half().tail();
            Vert vt2 = f.half().head();

            auto pairs = std::vector{
                std::pair(uv0, vt0.id),
                std::pair(uv1, vt1.id),
                std::pair(uv2, vt2.id)
            };

            rg::sort(pairs, [&fr](auto &a, auto &b) { return norm(fr - a.first) < norm(fr - b.first); });
            for (int i = 0; i < 3; i++) {
                auto vid = pairs[i].second;
                if (!visit[vid]) { vids0.emplace_back(vid); break; }
            }

            rg::sort(pairs, [&to](auto &a, auto &b) { return norm(to - a.first) < norm(to - b.first); });
            for (int i = 0; i < 3; i++) {
                auto vid = pairs[i].second;
                if (!visit[vid]) { vids0.emplace_back(vid); break; }
            }
        }

        for (int i = 0; i < vids0.size(); ++i) {
            if (i == 0 || vids0[i] != vids0[i - 1]) vids1.emplace_back(vids0[i]);
        }

        for (int i = 1; i < vids1.size() - 1; ++i) { visit[vids1[i]] = true; }

        for (int i = 0; i < vids1.size() - 1; ++i) {
            int i0 = vids1[i];
            int i1 = vids1[i + 1];
            auto it = rg::find_if(hm.halfs, [&](const Half& h) { return h.tail().id == i0 && h.head().id == i1; });
            //assert(it != hm.halfs.end());
            if (it == hm.halfs.end()) { std::cout << "ERROR: " << i0 << " " << i1 << std::endl;}
            else hs.emplace_back(*it);
        }

        double sum = 0;
        double cur = 0;
        for (const auto& h: hs) { sum += h.len(); }
        for (int i = 0; i < hs.size(); ++i) {
            const auto& h = hs[i];
            double v0 = cur;
            cur += h.len() / sum;
            double v1 = cur;
            std::cout << th.id << std::endl;
            h_data.emplace(HalfData{h, v0, v1, th.id, tm.th2quad(th.id), i});
        }
    }
}
}

#endif
