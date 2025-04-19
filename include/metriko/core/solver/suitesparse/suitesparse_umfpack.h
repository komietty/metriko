#ifndef SUITESPARSETEST_SUITESPARSE_UMFPACK_H
#define SUITESPARSETEST_SUITESPARSE_UMFPACK_H
#include "suitesparse_utilities.h"

namespace metriko::suitesparse {

template <typename T>
void umfFactor(size_t N, cholmod_sparse* mat, void*& symbolicFac, void*& numericFac);

template <>
void umfFactor<double>(size_t N, cholmod_sparse* mat, void*& symbolicFac, void*& numericFac) {
    auto* cMat_p = (SuiteSparse_long*)mat->p;
    auto* cMat_i = (SuiteSparse_long*)mat->i;
    auto* cMat_x = (double*)mat->x;
    umfpack_dl_symbolic(N, N, cMat_p, cMat_i, cMat_x, &symbolicFac, nullptr, nullptr);
    umfpack_dl_numeric(cMat_p, cMat_i, cMat_x, symbolicFac, &numericFac, nullptr, nullptr);
}
template <>
void umfFactor<float>(size_t N, cholmod_sparse* mat, void*& symbolicFac, void*& numericFac) {
    auto* cMat_p = (SuiteSparse_long*)mat->p;
    auto* cMat_i = (SuiteSparse_long*)mat->i;
    auto* cMat_x = (double*)mat->x;
    umfpack_dl_symbolic(N, N, cMat_p, cMat_i, cMat_x, &symbolicFac, nullptr, nullptr);
    umfpack_dl_numeric(cMat_p, cMat_i, cMat_x, symbolicFac, &numericFac, nullptr, nullptr);
}
template <>
void umfFactor<std::complex<double>>(size_t N, cholmod_sparse* mat, void*& symbolicFac, void*& numericFac) {
auto* cMat_p = (SuiteSparse_long*)mat->p;
auto* cMat_i = (SuiteSparse_long*)mat->i;
auto* cMat_x = (double*)mat->x;
umfpack_zl_symbolic(N, N, cMat_p, cMat_i, cMat_x, nullptr, &symbolicFac, nullptr, nullptr);
umfpack_zl_numeric(cMat_p, cMat_i, cMat_x, nullptr, symbolicFac, &numericFac, nullptr, nullptr);
}

template <typename T>
void umfSolve(size_t N, cholmod_sparse* mat, void* numericFac, VecX<T>& x, const VecX<T>& rhs);

template <>
void umfSolve<double>(size_t N, cholmod_sparse* mat, void* numericFac, VecX<double>& x, const VecX<double>& rhs) {
    x = VecX<double>(N);
    auto* cMat_p = (SuiteSparse_long*)mat->p;
    auto* cMat_i = (SuiteSparse_long*)mat->i;
    auto* cMat_x = (double*)mat->x;
    umfpack_dl_solve(UMFPACK_A, cMat_p, cMat_i, cMat_x, &(x[0]), &(rhs[0]), numericFac, nullptr, nullptr);
}
template <>
void umfSolve<std::complex<double>>(size_t N, cholmod_sparse* mat, void* numericFac, VecX<std::complex<double>>& x, const VecX<std::complex<double>>& rhs) {
// Note: the ordering of std::complex is specified by the standard, so this certainly works
x = VecX<std::complex<double>>(N);
auto* cMat_p = (SuiteSparse_long*)mat->p;
auto* cMat_i = (SuiteSparse_long*)mat->i;
auto* cMat_x = (double*)mat->x;
umfpack_zl_solve(UMFPACK_A, cMat_p, cMat_i, cMat_x, nullptr, (double*)&(x[0]), nullptr, (double*)&(rhs[0]), nullptr, numericFac, nullptr, nullptr);
}
}
#endif
