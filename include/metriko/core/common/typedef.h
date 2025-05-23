//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_TYPEDEF_H
#define METRIKO_TYPEDEF_H
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <ranges>
#include <__ranges/views.h>

namespace metriko {
    constexpr double PI = M_PI;
    constexpr double TwoPI = 2 * M_PI;
    constexpr double EPS = 1e-9;
    using complex = std::complex<double>;
    using Row2i = Eigen::RowVector2i;
    using Row2d = Eigen::RowVector2d;
    using Row3i = Eigen::RowVector3i;
    using Row3d = Eigen::RowVector3d;
    using Mat2i = Eigen::Matrix2i;
    using Mat2d = Eigen::Matrix2d;
    using Vec2d = Eigen::Matrix<double, 2, 1>;
    using Vec3d = Eigen::Matrix<double, 3, 1>;
    using VecXi = Eigen::Matrix<int, Eigen::Dynamic, 1>;
    using VecXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
    using VecXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using VecXc = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
    using SprsI = Eigen::SparseMatrix<int>;
    using SprsD = Eigen::SparseMatrix<double>;
    using SprsC = Eigen::SparseMatrix <std::complex<double>>;
    using MatXi = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
    using MatXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using MatXc = Eigen::Matrix <std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
    using TripI = Eigen::Triplet<int>;
    using TripD = Eigen::Triplet<double>;
    using TripC = Eigen::Triplet <std::complex<double>>;
    namespace rg = std::ranges;
    namespace vw = std::views;
}

#endif
