//
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_SOLVER_COMMON_H
#define METRIKO_SOLVER_COMMON_H
#include "metriko/core/common/typedef.h"

namespace metriko {
inline void reduce_to_linearly_independent(SprsD& mat) {
    if (mat.rows() == 0) return;
    Eigen::SparseQR<SprsD, Eigen::COLAMDOrdering<int>> qr;
    qr.compute(mat.transpose());
    int rank = qr.rank();
    const VecXi &idcs = qr.colsPermutation().indices(); // the remaining row idcs of original mtx

    VecXi rid_map = VecXi::Constant(mat.rows(), -1);
    for (int j = 0; j < rank; ++j) rid_map(idcs(j)) = j; // original mat row idx -> QR decomped mat row idx

    std::vector<TripD> T;
    T.reserve(mat.nonZeros());
    for (int k = 0; k < mat.outerSize(); ++k) {
    for (SprsD::InnerIterator it(mat, k); it; ++it) {
        if (int rid = rid_map(it.row()); rid != -1) T.emplace_back(rid, it.col(), it.value());
    }}

    mat.resize(rank, mat.cols());
    mat.setFromTriplets(T.begin(), T.end());
}

template<typename Scalar>
void sparse_block(
    const MatXi &idcs,
    const std::vector<Eigen::SparseMatrix<Scalar> *> &mats,
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

    std::vector<Eigen::Triplet<Scalar>> T;
    for (int i = 0; i < row_offsets.size(); i++)
    for (int j = 0; j < col_offsets.size(); j++)
    for (int k = 0; k < mats[i]->outerSize(); ++k)
    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(*(mats[i]), k); it; ++it)
        T.push_back(Eigen::Triplet<Scalar>(row_offsets(i) + it.row(), col_offsets(j) + it.col(), it.value()));

    result.setFromTriplets(T.begin(), T.end());
}

struct SprsEntry { int row; int col; double value; };

// all nonzeros, column by column
inline vec<SprsEntry> nonzeros(const SprsD& S) {
    vec<SprsEntry> res;
    res.reserve(S.nonZeros());
    for (int k = 0; k < S.outerSize(); ++k)
    for (SprsD::InnerIterator it(S, k); it; ++it)
        res.push_back({.row=(int)it.row(), .col=(int)it.col(), .value=it.value()});
    return res;
}

// for each row, the columns holding a nonzero on it, in column order
inline vec<vec<int>> cols_by_row(const SprsD& S) {
    vec<vec<int>> res(S.rows());
    for (const auto& [r, c, _]: nonzeros(S)) res[r].push_back(c);
    return res;
}
}

#endif
