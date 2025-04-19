//
//--- Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_SOLVER_COMMON_H
#define METRIKO_SOLVER_COMMON_H
#include "metriko/core/common/typedef.h"

namespace metriko {

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
