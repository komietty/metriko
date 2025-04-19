#ifndef SUITESPARSETEST_SUITESPARSE_SQUARE_H
#define SUITESPARSETEST_SUITESPARSE_SQUARE_H
#include "suitesparse_utilities.h"
#include "suitesparse_umfpack.h"

namespace metriko::suitesparse {

template <typename T>
struct SquareSolverInternals {
    CholmodContext context;
    cholmod_sparse* cMat = nullptr;
    void* symbolicFactorization = nullptr;
    void* numericFactorization = nullptr;
};

template <typename T>
class SquareSolver : public LinearSolver<T> {
public:
    SquareSolver(Eigen::SparseMatrix<T>& mat);
    ~SquareSolver();
    void solve(VecX<T>& x, const VecX<T>& rhs) override;

protected:
    std::unique_ptr<SquareSolverInternals<T>> internals;
};

template <typename T>
SquareSolver<T>::SquareSolver(Sprs<T>& mat) : LinearSolver<T>(mat), internals(new SquareSolverInternals<T>()) {
    mat.makeCompressed();
    if (internals->cMat != nullptr) cholmod_l_free_sparse(&internals->cMat, internals->context);
    internals->cMat = toCholmod(mat, internals->context);
    umfFactor<T>(this->nRows, internals->cMat, internals->symbolicFactorization, internals->numericFactorization);
}

template <typename T>
SquareSolver<T>::~SquareSolver() {
    if (internals->cMat != nullptr) {
        cholmod_l_free_sparse(&internals->cMat, internals->context);
        internals->cMat = nullptr;
    }
    if (internals->symbolicFactorization != nullptr) {
        umfpack_dl_free_symbolic(&internals->symbolicFactorization);
    }
    if (internals->numericFactorization != nullptr) {
        umfpack_dl_free_numeric(&internals->numericFactorization);
    }
}


template <typename T>
void SquareSolver<T>::solve(VecX<T>& x, const VecX<T>& rhs) {
    size_t N = this->nRows;
    umfSolve<T>(N, internals->cMat, internals->numericFactorization, x, rhs);
}
}
#endif
