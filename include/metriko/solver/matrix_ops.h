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

// a column-major sparse matrix written column by column in parallel: column c gets count(c) entries,
// written by fill(c, rows, values) in ascending row order
template<class Count, class Fill>
void sparse_by_columns(int nr, int nc, SprsD& res, Count&& count, Fill&& fill) {
    res.resize(nr, nc);
    int* op = res.outerIndexPtr();
    #pragma omp parallel for
    for (int c = 0; c < nc; ++c) op[c + 1] = count(c);
    for (int c = 0; c < nc; ++c) op[c + 1] += op[c];
    res.resizeNonZeros(op[nc]);
    #pragma omp parallel for
    for (int c = 0; c < nc; ++c) fill(c, res.innerIndexPtr() + op[c], res.valuePtr() + op[c]);
}

// the matrices, all with the same columns, stacked on top of each other
inline void sparse_vstack(const vec<const SprsD*>& mats, SprsD& res) {
    int nr = 0;
    for (auto m: mats) nr += m->rows();
    sparse_by_columns(nr, mats[0]->cols(), res,
        [&](int c) { int n = 0; for (auto m: mats) n += m->col(c).nonZeros(); return n; },
        [&](int c, int* ri, double* rv) {
            int oft = 0;
            for (auto m: mats) {
                for (SprsD::InnerIterator it(*m, c); it; ++it) { *ri++ = oft + it.row(); *rv++ = it.value(); }
                oft += m->rows();
            }
        });
}

// for each row, the columns holding a nonzero on it, in column order
inline vec<vec<int>> cols_by_row(const SprsD& S) {
    vec<vec<int>> res(S.rows());
    for (const auto& [r, c, _]: nonzeros(S)) res[r].push_back(c);
    return res;
}
}
#endif
