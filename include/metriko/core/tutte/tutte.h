//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_TUTTE_H
#define METRIKO_TUTTE_H
#include <tuple>
#include "metriko/core/hmesh/hmesh.h"

namespace metriko {
struct HalfData {
    Half half = Half(); //
    double v0 {};       // parameter value of tail
    double v1 {};       // parameter value of head
    int thid  = -1;     //
    int tqid  = -1;     //
    int order = -1;     // order inside thalf

    bool operator<(const HalfData& rhs) const noexcept { return std::tie(tqid, thid, order) < std::tie(rhs.tqid, rhs.thid, rhs.order); }
};
}
#endif
