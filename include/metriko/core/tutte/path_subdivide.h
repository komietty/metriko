#ifndef TMESH_H_PATH_SUBDIVIDE_H
#define TMESH_H_PATH_SUBDIVIDE_H
#include <cstddef>
#include <ctime>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "emesh.h"
#include "tutte_cutting.h"
#include "glm/vec4.hpp"
#include "glm/ext/matrix_common.hpp"
#include "metriko/core/common/typedef.h"
#include "metriko/core/hmesh/hmesh.h"

namespace metriko::tutte {

inline std::pair<Hmesh, Tmesh> subdivide_hmesh_by_tmesh(
    const Hmesh& hm,     // input hmesh
    const tm::Tmesh& tm  // input tmesh
) {
    for (const auto& th: tm.thalfs) {
        if (!th.cano) continue;
        const auto& te = th.edge();
        const double r = R[te.id];
        double sum = 0;
        int order = 0;
        for (const auto& s: te.segments()) {
            auto len = abs(s.diff());
            auto v0 = sum / r; sum += len;
            auto v1 = sum / r;
            auto aux = AuxSgmt{s, th, order, te.n_segments() - order, {v0, v1}, {1 - v1, 1 - v0}};
            order++;
            cuts[s.face.id].emplace_back(aux);
        }
    }

    return hm_cut;
}
}

#endif
