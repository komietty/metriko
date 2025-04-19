//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_BASE_FIELD_H
#define METRIKO_BASE_FIELD_H
#include "metriko/core/hmesh/hmesh.h"

namespace metriko {
    class BaseVectorField {
    public:
        const Hmesh& mesh;
        const int rosyN;
        MatXc field;
        VecXc connection;
        VecXc compressed;
        VecXi matching;
        VecXi singular;
        BaseVectorField(const Hmesh& m, const int n) : mesh(m), rosyN(n) { }
    };

    enum class FieldType { UnSpecified, Smoothest, CurvatureAligned };
    enum class MatchingType { Principal, Curl };
}

#endif
