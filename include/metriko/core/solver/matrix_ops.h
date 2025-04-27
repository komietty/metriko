//
//--- Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_SOLVER_COMMON_H
#define METRIKO_SOLVER_COMMON_H
#include "metriko/core/common/typedef.h"

namespace metriko {

    inline void reduce_to_linearly_independent(SprsD& mat) {
        if (mat.rows() != 0) {
            Eigen::SparseQR<SprsD, Eigen::COLAMDOrdering<int>> qr;
            qr.compute(mat.transpose());
            int rank = qr.rank();
            VecXi idcs = qr.colsPermutation().indices(); //creating a sliced perm matrix

            std::vector<TripD> T;
            for (int k = 0; k < mat.outerSize(); ++k) {
            for (SprsD::InnerIterator it(mat, k); it; ++it) {
            for (int j = 0; j < rank; j++) {
                if (it.row() == idcs(j)) T.emplace_back(j, it.col(), it.value());
            }}}

            mat.resize(rank, mat.cols());
            mat.setFromTriplets(T.begin(), T.end());
        }
    }

    template<typename Scalar>
    void sparse_block (
            const Eigen::MatrixXi &blockIndices,
            const std::vector<Eigen::SparseMatrix<Scalar>*> &blockMats,
            Eigen::SparseMatrix<Scalar> &result
    ) {
        //assessing dimensions
        VecXi blockRowOffsets = VecXi::Zero(blockIndices.rows());
        VecXi blockColOffsets = VecXi::Zero(blockIndices.cols());
        for (int i = 1; i < blockIndices.rows(); i++)
            blockRowOffsets(i) = blockRowOffsets(i - 1) + blockMats[blockIndices(i - 1, 0)]->rows();

        for (int i = 1; i < blockIndices.cols(); i++)
            blockColOffsets(i) = blockColOffsets(i - 1) + blockMats[blockIndices(0, i - 1)]->cols();

        int rowSize = blockRowOffsets(blockIndices.rows() - 1) + blockMats[blockIndices(blockIndices.rows() - 1, 0)]->rows();
        int colSize = blockColOffsets(blockIndices.cols() - 1) + blockMats[blockIndices(0, blockIndices.cols() - 1)]->cols();

        result.conservativeResize(rowSize, colSize);
        std::vector<Eigen::Triplet<Scalar>> resultTriplets;
        for (int i = 0; i < blockRowOffsets.size(); i++)
        for (int j = 0; j < blockColOffsets.size(); j++)
        for (int k = 0; k < blockMats[i]->outerSize(); ++k)
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(*(blockMats[i]), k); it; ++it)
            resultTriplets.push_back(Eigen::Triplet<Scalar>(blockRowOffsets(i) + it.row(), blockColOffsets(j) + it.col(), it.value()));

        result.setFromTriplets(resultTriplets.begin(), resultTriplets.end());
    }

}

#endif
