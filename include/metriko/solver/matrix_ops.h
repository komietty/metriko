//
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_MATRIX_OPS_H
#define METRIKO_MATRIX_OPS_H
#include "metriko/common/typedef.h"
#include "metriko/solver/sparse_solver.h"

namespace metriko {
template<class Scalar> struct SparseEntry { int row; int col; Scalar value; };
using SprsEntry = SparseEntry<double>;

// the nonzeros of a column-major sparse matrix as a range of {row, col, value}, column by column, without copying.
// the matrix must outlive the range and keep its sparsity pattern while it is walked
template<class Scalar>
auto nonzeros(const Eigen::SparseMatrix<Scalar>& S) {
    return vw::iota(0, (int)S.outerSize()) | vw::transform([&S](int k) {
        const int b = S.outerIndexPtr()[k];
        const int e = S.isCompressed() ? S.outerIndexPtr()[k + 1] : b + S.innerNonZeroPtr()[k];
        return vw::iota(b, e) | vw::transform([&S, k](int p) { return SparseEntry<Scalar>{ (int)S.innerIndexPtr()[p], k, S.valuePtr()[p] }; });
    }) | vw::join;
}

// the solutions of C x = b for a sparse C whose nonzeros sit in a few columns, as x = x0 + Z y. the columns C does not
// touch stay free: identity columns of Z, in their order. the touched ones are confined to the kernel of C restricted to
// them, computed densely by FullPivLU, whose vectors follow in Z. x0 solves C x0 = b and is zero off the touched columns
inline SprsD sparse_null_space(const SprsD& C, const VecXd& b, VecXd& x0) {
    const int nx = C.cols();
    vec<int> cols;   // the touched columns

    for (const auto& [r, c, v]: nonzeros(C)) cols.push_back(c);

    rg::sort(cols);
    cols.erase(rg::unique(cols).begin(), cols.end());
    VecXi loc = VecXi::Constant(nx, -1);
    for (int j = 0; j < cols.size(); ++j) loc(cols[j]) = j;

    vec<TripD> T;
    int ny = 0;
    for (int i = 0; i < nx; ++i) if (loc(i) == -1) T.emplace_back(i, ny++, 1.);
    x0 = VecXd::Zero(nx);

    if (!cols.empty()) {
        MatXd Cd = MatXd::Zero(C.rows(), cols.size());
        for (const auto& [r, c, v]: nonzeros(C)) Cd(r, loc(c)) = v;

        Eigen::FullPivLU<MatXd> lu(Cd);
        const VecXd xc = lu.solve(b);
        for (int j = 0; j < cols.size(); ++j) x0(cols[j]) = xc(j);

        if (lu.dimensionOfKernel() > 0) {
            const MatXd K = lu.kernel();
            for (int k = 0; k < K.cols(); ++k, ++ny)
            for (int j = 0; j < K.rows(); ++j)
                if (K(j, k) != 0) T.emplace_back(cols[j], ny, K(j, k));
        }
    }

    SprsD Z(nx, ny);
    Z.setFromTriplets(T.begin(), T.end());
    return Z;
}

// the homogeneous case, C x = 0: x = Z y
inline SprsD sparse_null_space(const SprsD& C) {
    VecXd x0;
    return sparse_null_space(C, VecXd::Zero(C.rows()), x0);
}

inline void reduce_to_linearly_independent(SprsD& mat) {
    if (mat.rows() == 0) return;
    SparseRankQR qr;
    qr.compute(mat.transpose());
    int rank = qr.rank();
    const VecXi &idcs = qr.colsPermutation().indices(); // the remaining row idcs of original mtx

    VecXi rid_map = VecXi::Constant(mat.rows(), -1);
    for (int j = 0; j < rank; ++j) rid_map(idcs(j)) = j; // original mat row idx -> QR decomped mat row idx

    vec<TripD> T;
    T.reserve(mat.nonZeros());
    for (const auto& [r, c, v]: nonzeros(mat)) if (int rid = rid_map(r); rid != -1) T.emplace_back(rid, c, v);

    mat.resize(rank, mat.cols());
    mat.setFromTriplets(T.begin(), T.end());
}

template<typename Scalar>
void sparse_block(
    const MatXi &idcs,
    const vec<Eigen::SparseMatrix<Scalar> *> &mats,
    Eigen::SparseMatrix<Scalar> &result
) {
    //assessing dimensions
    int row_oft = idcs.rows();
    int col_oft = idcs.cols();
    VecXi row_offsets = VecXi::Zero(row_oft);
    VecXi col_offsets = VecXi::Zero(col_oft);
    for (int i = 1; i < row_oft; i++) row_offsets(i) = row_offsets(i - 1) + mats[idcs(i - 1, 0)]->rows();
    for (int i = 1; i < col_oft; i++) col_offsets(i) = col_offsets(i - 1) + mats[idcs(0, i - 1)]->cols();

    result.conservativeResize(
        row_offsets(row_oft - 1) + mats[idcs(row_oft - 1, 0)]->rows(),
        col_offsets(col_oft - 1) + mats[idcs(0, col_oft - 1)]->cols()
    );

    vec<Eigen::Triplet<Scalar>> T;
    for (int i = 0; i < row_offsets.size(); i++)
    for (int j = 0; j < col_offsets.size(); j++)
    for (const auto& [r, c, v]: nonzeros(*mats[i]))
        T.emplace_back(row_offsets(i) + r, col_offsets(j) + c, v);

    result.setFromTriplets(T.begin(), T.end());
}

// for each row, the columns holding a nonzero on it, in column order
inline vec<vec<int>> cols_by_row(const SprsD& S) {
    vec<vec<int>> res(S.rows());
    for (const auto& [r, c, _]: nonzeros(S)) res[r].push_back(c);
    return res;
}
}
#endif
