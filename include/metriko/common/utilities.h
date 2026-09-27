//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_COMMON_UTILITIES_H
#define METRIKO_COMMON_UTILITIES_H
#include <iostream>
#include "typedef.h"

namespace metriko {
template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };

constexpr auto circular_prev = [](auto& c, auto it) { return it == c.begin() ? std::prev(c.end()) : std::prev(it); };
constexpr auto circular_next = [](auto& c, auto it) { auto n = std::next(it); return n == c.end() ? c.begin() : n; };

inline complex get_quater_rot(int i) {
    switch ((i % 4 + 4) % 4) {
        case 0: return {1, 0};
        case 1: return {0, 1};
        case 2: return {-1, 0};
        case 3: return {0, -1};
        default: throw std::invalid_argument("Invalid index");
    }
}

inline bool equal(
    const complex a,
    const complex b,
    const double tolerance = 1e-10
) {
    return abs(a - b) < tolerance;
}

inline complex lerp(
    const complex a,
    const complex b,
    const double t
) {
    return a * (1. - t) + b * t;
}

inline complex normalize(
    const complex a
) {
    double l = abs(a);
    METRIKO_CHECK(l > 0, "normalize of a zero vector");
    return a / l;
}

inline double dot(
    const complex a,
    const complex b
) {
    return a.real() * b.real() + a.imag() * b.imag();
}

inline double cross(
    const complex a,
    const complex b
) {
    return a.real() * b.imag() - a.imag() * b.real();
}

struct minmax_int {
    int min_x;
    int min_y;
    int max_x;
    int max_y;
};

inline minmax_int get_minmax_int(std::vector<complex> uvs) {
    auto xs = vw::transform(uvs, [](complex uv) { return uv.real(); });
    auto ys = vw::transform(uvs, [](complex uv) { return uv.imag(); });
    return minmax_int{
        .min_x = static_cast<int>(std::floor(rg::min(xs))),
        .min_y = static_cast<int>(std::floor(rg::min(ys))),
        .max_x = static_cast<int>(std::ceil(rg::max(xs))),
        .max_y = static_cast<int>(std::ceil(rg::max(ys)))
    };
};

inline complex calc_coefficient(
    const complex ori,
    const complex pt1,
    const complex pt2,
    const complex tgt
) {
    auto v1  = pt1 - ori;
    auto v2  = pt2 - ori;
    auto dif = tgt - ori;
    auto det = cross(v1, v2);
    METRIKO_CHECK(det != 0, "degenerate uv triangle");
    return {cross(dif, v2) / det, cross(v1, dif) / det};
}

// Find the intersection of a line passing through a and b and another line passing through c and d.
// Return true if segments are not parallel. Check 0 <= ratio <= 1 if you want a segment-segment intersection.
inline bool find_extended_intersection(
    const complex &a,
    const complex &b,
    const complex &c,
    const complex &d,
    double &ratio_a2b,
    double &ratio_c2d
) {
    auto e1  = d - c, e2 = b - a, w = a - c;
    auto det = cross(e1, e2);
    if (det == 0) return false;
    ratio_c2d = cross(w, e2) / det;
    ratio_a2b = cross(w, e1) / det;
    return true;
}

inline bool find_strict_intersection(
    const complex &a,  // segment AB
    const complex &b,  // segment AB
    const complex &c,  // segment CD
    const complex &d,  // segment CD
    const double eps = EPS
) {
    double ratio_a2b; // ratio in AB
    double ratio_c2d;// ratio in CD
    return find_extended_intersection(a, b, c, d, ratio_a2b, ratio_c2d) &&
           ratio_a2b >= 0 + eps &&
           ratio_a2b <= 1 - eps &&
           ratio_c2d >= 0 + eps &&
           ratio_c2d <= 1 - eps;
}

inline bool find_strict_intersection(
    const complex &a,  // segment AB
    const complex &b,  // segment AB
    const complex &c,  // segment CD
    const complex &d,  // segment CD
    double &ratio_a2b, // ratio in AB
    double &ratio_c2d, // ratio in CD
    const double eps = EPS
) {
    return find_extended_intersection(a, b, c, d, ratio_a2b, ratio_c2d) &&
           ratio_a2b >= 0 + eps &&
           ratio_a2b <= 1 - eps &&
           ratio_c2d >= 0 + eps &&
           ratio_c2d <= 1 - eps;
}

inline Row3d conversion_2d_3d(
    const complex o2, // origin of 2d
    const complex a2, //
    const complex b2, //
    const Row3d &o3,  // origin of 3d
    const Row3d &a3,  //
    const Row3d &b3,  //
    const complex uv  // target uv value
) {
    const complex c = calc_coefficient(o2, a2, b2, uv);
    return o3 + (a3 - o3) * c.real() + (b3 - o3) * c.imag();
}
}
#endif
