//
// Created by saki on 2026/07/14.
//

#ifndef TMESH_MUT_COLLAPSE_TQUAD_CPP_EMESH_H
#define TMESH_MUT_COLLAPSE_TQUAD_CPP_EMESH_H
#include <tuple>
#include "metriko/core/hmesh/hmesh.h"

namespace metriko {
struct HalfData {
    Half half;
    double v0;
    double v1;
    int thid;
    int tqid;
    int twin;
    int order; // order inside thalf

    bool operator<(const HalfData& rhs) const noexcept {
        return std::tie(tqid, thid, order) < std::tie(rhs.tqid, rhs.thid, rhs.order);
    }
};
}

#endif
