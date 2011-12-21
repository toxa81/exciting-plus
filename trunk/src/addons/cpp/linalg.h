#ifndef __LINALG_H__
#define __LINALG_H__

#include <stdint.h>
#include <complex>
#include "config.h"

extern "C" void FORTFUNC(zgemm)(
    const char *transa, 
    const char *transb, 
    int32_t *m, 
    int32_t *n, 
    int32_t *k, 
    std::complex<double> *alpha, 
    std::complex<double> *a, 
    int32_t *lda, 
    std::complex<double> *b, 
    int32_t *ldb, 
    std::complex<double> *beta, 
    std::complex<double> *c, 
    int32_t *ldc, 
    int32_t transalen,
    int32_t transblen);

extern "C" void zgemm_cublas(
    int transa,
    int transb, 
    int32_t *m, 
    int32_t *n, 
    int32_t *k,
    std::complex<double> *alpha, 
    std::complex<double> *a, 
    int32_t *lda, 
    std::complex<double> *b, 
    int32_t *ldb, 
    std::complex<double> *beta, 
    std::complex<double> *c, 
    int32_t *ldc);
    
template<int N> void zgemm(
    int transa, 
    int transb, 
    int32_t *m, 
    int32_t *n, 
    int32_t *k,
    std::complex<double> *alpha, 
    std::complex<double> *a, 
    int32_t *lda, 
    std::complex<double> *b, 
    int32_t *ldb, 
    std::complex<double> *beta, 
    std::complex<double> *c, 
    int32_t *ldc) {
    
    if (N == blas_worker)
    {
        const char *trans[] = {"N", "T", "C"};
        FORTFUNC(zgemm)(trans[transa], trans[transb], m, n, k, alpha, 
            a, lda, b, ldb, beta, c, ldc, (int32_t)1, (int32_t)1);
    }
    
    if (N == cublas_worker)
    {    
        zgemm_cublas(transa, transb, m, n, k, alpha, a, lda, b, ldb, 
            beta, c, ldc);
    }
}

#endif // __LINALG_H__

