#ifndef SUITESPARSETEST_SUITESPARSE_PDEFINITE_H
#define SUITESPARSETEST_SUITESPARSE_PDEFINITE_H
#include "suitesparse_utilities.h"

namespace metriko::suitesparse {

template <typename T>
struct PSDSolverInternals {
    CholmodContext context;
    cholmod_sparse* cMat = nullptr;
    cholmod_factor* factorization = nullptr;
};

template <typename T>
class PositiveDefiniteSolver : public LinearSolver<T> {
public:
    PositiveDefiniteSolver(Sprs<T>& mat);
    ~PositiveDefiniteSolver();
    void solve(VecX<T>& x, const VecX<T>& rhs) override;

protected:
    std::unique_ptr<PSDSolverInternals<T>> internals;
};

template <typename T>
PositiveDefiniteSolver<T>::PositiveDefiniteSolver(Sprs<T>& mat): LinearSolver<T>(mat),  internals(new PSDSolverInternals<T>()) {
    if (internals->cMat != nullptr) cholmod_l_free_sparse(&internals->cMat, internals->context);
    internals->cMat = toCholmod(mat, internals->context, SType::SYMMETRIC);

    internals->context.setSimplicial(); // must use simplicial for LDLt
    internals->context.setLDL();        // ensure we get an LDLt internals->factorization
    internals->factorization = cholmod_l_analyze(internals->cMat, internals->context);
    bool success = (bool)cholmod_l_factorize(internals->cMat, internals->factorization, internals->context);

    if(!success) {
        throw std::runtime_error("failure in cholmod_l_factorize");
    }
    if(internals->context.context.status == CHOLMOD_NOT_POSDEF) {
        throw std::runtime_error("matrix is not positive definite");
    }
}

template <typename T>
PositiveDefiniteSolver<T>::~PositiveDefiniteSolver() {
    if (internals->cMat != nullptr) {
        cholmod_l_free_sparse(&internals->cMat, internals->context);
        internals->cMat = nullptr;
    }
    if (internals->factorization != nullptr) {
        cholmod_l_free_factor(&internals->factorization, internals->context);
    }
}

template <typename T>
void PositiveDefiniteSolver<T>::solve(VecX<T>& x, const VecX<T>& rhs) {
    cholmod_dense* inVec = toCholmod(rhs, internals->context);
    cholmod_dense* outVec = cholmod_l_solve(CHOLMOD_A, internals->factorization, inVec, internals->context);
    toEigen(outVec, x);
    cholmod_l_free_dense(&outVec, internals->context);
    cholmod_l_free_dense(&inVec, internals->context);
}
}
#endif
