//
// Created by saki on 2026/07/14.
//

#ifndef TMESH_MUT_COLLAPSE_TQUAD_CPP_EMESH_H
#define TMESH_MUT_COLLAPSE_TQUAD_CPP_EMESH_H
#include <tuple>
#include "metriko/core/hmesh/hmesh.h"

namespace metriko {
struct HalfData {
    Half half = Half(); //
    double v0;          // parameter value of tail
    double v1;          // parameter value of head
    int thid = -1;      //
    int tqid = -1;      //
    int twin = -1;      //
    int order = -1;     // order inside thalf

    bool operator<(const HalfData& rhs) const noexcept { return std::tie(tqid, thid, order) < std::tie(rhs.tqid, rhs.thid, rhs.order); }
};
}

#endif
