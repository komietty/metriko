#ifndef SUITESPARSETEST_INCLUDE_H
#define SUITESPARSETEST_INCLUDE_H
#include "Eigen/Sparse"
#include <SuiteSparseQR.hpp>
#include <cholmod.h>
#include <umfpack.h>


namespace metriko::suitesparse {

template<typename K> using VecX = Eigen::Matrix<K, Eigen::Dynamic, 1>;
template<typename K> using Sprs = Eigen::SparseMatrix<K>;

enum class SType { UNSYMMETRIC = 0, SYMMETRIC };

template <typename T>
class LinearSolver {
public:
    LinearSolver(const Sprs<T>& mat) : nRows(mat.rows()), nCols(mat.cols()) {}
    virtual ~LinearSolver() {}
    virtual void solve(VecX<T>& x, const VecX<T>& rhs) = 0;
protected:
    size_t nRows, nCols;
};

class CholmodContext {
public:
    CholmodContext() { cholmod_l_start(&context); }
    ~CholmodContext() { cholmod_l_finish(&context); }
    void setSimplicial() { context.supernodal = CHOLMOD_SIMPLICIAL; }
    void setSupernodal() { context.supernodal = CHOLMOD_SIMPLICIAL; }
    void setLL() { context.final_ll = true; }
    void setLDL() {context.final_ll = false; }
    operator cholmod_common*() { return &context; }
    cholmod_common context;
};

namespace {
int flagForStype(SType s) {
    switch (s) {
        case SType::UNSYMMETRIC: { return 0; }
        case SType::SYMMETRIC:   { return 1; }
    }
    return -1;
}
}

// Type helper. This type is 'double' if T == 'float', and T otherwise
//template<typename T>
//struct SOLVER_ENTRYTYPE {
//    typedef typename std::conditional<std::is_same<T, float>::value, double, T>::type type;
//};
// === Conversion functions
template <typename T>
cholmod_sparse* toCholmod(Eigen::SparseMatrix<T, Eigen::ColMajor>& A, CholmodContext& context, SType stype = SType::UNSYMMETRIC);

// Convert a vector
template <typename T>
cholmod_dense* toCholmod(const VecX<T>& v, CholmodContext& context);

template <>
cholmod_dense* toCholmod(const VecX<double>& v, CholmodContext& context) {
    size_t N = v.rows();
    cholmod_dense* cVec = cholmod_l_allocate_dense(N, 1, N, CHOLMOD_REAL, context);
    auto* cVecD = (double*)cVec->x;
    for (size_t i = 0; i < N; i++) cVecD[i] = v(i);
    return cVec;
}

template <>
cholmod_dense* toCholmod(const VecX<std::complex<double>>& v, CholmodContext& context) {
    size_t N = v.rows();
    cholmod_dense* cVec = cholmod_l_allocate_dense(N, 1, N, CHOLMOD_COMPLEX, context);
    auto* cVecC = (std::complex<double>*)cVec->x;
    for (size_t i = 0; i < N; i++) cVecC[i] = v(i);
    return cVec;
}

template <>
cholmod_sparse* toCholmod(Sprs<double>& A, CholmodContext& context, SType stype) {
    A.makeCompressed();
    // Allocate spase
    size_t Nentries = A.nonZeros();
    size_t Ncols = A.cols();
    size_t Nrows = A.rows();

    cholmod_sparse* cMat = cholmod_l_allocate_sparse(Nrows, Ncols, Nentries, true, true, flagForStype(stype), CHOLMOD_REAL, context);

    // Pull out useful pointers
    auto* values = (double*)cMat->x;
    auto* rowIndices = (SuiteSparse_long*)cMat->i;
    auto* colStart = (SuiteSparse_long*)cMat->p;

    // Copy
    for (size_t iEntry = 0; iEntry < Nentries; iEntry++) {
        values[iEntry] = A.valuePtr()[iEntry];
        rowIndices[iEntry] = A.innerIndexPtr()[iEntry];
    }
    for (size_t iCol = 0; iCol < Ncols; iCol++) {
        colStart[iCol] = A.outerIndexPtr()[iCol];
    }
    colStart[Ncols] = Nentries;
    return cMat;
}

template <>
cholmod_sparse* toCholmod(Sprs<std::complex<double>>& A, CholmodContext& context, SType stype) {
    A.makeCompressed();
    // Allocate spase
    size_t Nentries = A.nonZeros();
    size_t Ncols = A.cols();
    size_t Nrows = A.rows();

    cholmod_sparse* cMat = cholmod_l_allocate_sparse(Nrows, Ncols, Nentries, true, true, flagForStype(stype), CHOLMOD_COMPLEX, context);

    // Pull out useful pointers
    auto* values = (std::complex<double>*)cMat->x;
    auto* rowIndices = (SuiteSparse_long*)cMat->i;
    auto* colStart = (SuiteSparse_long*)cMat->p;

    // Copy
    for (size_t iEntry = 0; iEntry < Nentries; iEntry++) {
        values[iEntry] = A.valuePtr()[iEntry];
        rowIndices[iEntry] = A.innerIndexPtr()[iEntry];
    }
    for (size_t iCol = 0; iCol < Ncols; iCol++) {
        colStart[iCol] = A.outerIndexPtr()[iCol];
    }
    colStart[Ncols] = Nentries;

    return cMat;
}

template <typename T>
void toEigen(cholmod_dense* cVec, VecX<T>& xOut) {
    size_t N = cVec->nrow;
    xOut = VecX<T>(N);
    // Type wizardry. This type is 'double' if T == 'float', and T otherwise
    // Needed because cholmod always uses double precision
    typedef typename std::conditional<std::is_same<T, float>::value, double, T>::type SCALAR_TYPE;
    auto* cVecS = (SCALAR_TYPE*)cVec->x;
    for (size_t i = 0; i < N; i++) xOut(i) = cVecS[i];
}
template void toEigen(cholmod_dense* cVec, Eigen::Matrix<double, Eigen::Dynamic, 1>& xOut);
template void toEigen(cholmod_dense* cVec, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& xOut);
}

#endif
