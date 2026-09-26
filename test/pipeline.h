#ifndef METRIKO_TEST_PIPELINE_H
#define METRIKO_TEST_PIPELINE_H
#include <string>
#include "metriko/core/vectorfield/base_field.h"

namespace metriko {
inline FieldType parse_field_type(const std::string& s) {
    return s == "curvature_aligned" ? FieldType::CurvatureAligned : FieldType::Smoothest;
}
}
#endif
