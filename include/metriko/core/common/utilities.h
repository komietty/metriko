//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_COMMON_UTILITIES_H
#define METRIKO_COMMON_UTILITIES_H
#include <iostream>
#include "typedef.h"

namespace metriko {
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
        assert(l > 0);
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

    inline bool isccw(
        const complex a,
        const complex b,
        const complex c
    ) {
        const double v = cross(b - a, c - a);
        if (v == 0) std::cout << "Edge case!: cross == 0" << std::endl;
        return v > 0;
    }

    inline double orientation(
        const complex a,
        const complex b,
        const complex c
    ) {
        return cross(b - a, c - a);
    }

    inline bool points_into(
        const complex dir,
        const complex a,
        const complex b,
        const complex c
    ) {
        return isccw(a, b, a + dir) && isccw(a, a + dir, c);
    }

    inline int inside(
        const complex a,
        const complex b,
        const complex c,
        const complex tgt,
        const double precision = 1e-6
    ) {
        if (abs(tgt - a) < precision || abs(tgt - b) < precision || abs(tgt - c) < precision) return 2;
        return 0;
    }

    inline complex calc_coefficient(
        const complex origin,
        const complex point1,
        const complex point2,
        const complex target
    ) {
        const complex v1 = point1 - origin;
        const complex v2 = point2 - origin;
        const complex d = target - origin;
        Mat2d m;
        m << v1.real(), v2.real(), v1.imag(), v2.imag();
        assert(m.determinant() != 0);
        const Row2d x = m.inverse() * Row2d(d.real(), d.imag()).transpose();
        return {x[0], x[1]};
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
        Mat2d m;
        m << d.real() - c.real(), -(b - a).real(), d.imag() - c.imag(), -(b - a).imag();
        if (m.determinant() == 0) return false;
        Vec2d v = m.inverse() * Vec2d(-c.real() + a.real(), -c.imag() + a.imag());
        ratio_a2b = v.y();
        ratio_c2d = v.x();
        return true;
    }

    inline bool find_strict_intersection(
        const complex &a,
        const complex &b,
        const complex &c,
        const complex &d,
        double &ratio_a2b,
        double &ratio_c2d
    ) {
        return find_extended_intersection(a, b, c, d, ratio_a2b, ratio_c2d) &&
               ratio_a2b > 0 &&
               ratio_a2b < 1 &&
               ratio_c2d > 0 &&
               ratio_c2d < 1;
    }
}

#endif
