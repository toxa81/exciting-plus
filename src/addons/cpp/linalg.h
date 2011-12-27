#ifndef __LINALG_H__
#define __LINALG_H__

#include <stdint.h>
#include <string.h>
#include <iostream>
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
    int32_t m, 
    int32_t n, 
    int32_t k,
    std::complex<double> alpha, 
    std::complex<double> *a, 
    int32_t lda, 
    std::complex<double> *b, 
    int32_t ldb, 
    std::complex<double> beta, 
    std::complex<double> *c, 
    int32_t ldc);
    
template<int N> void zgemm(
    int transa, 
    int transb, 
    int32_t m, 
    int32_t n, 
    int32_t k,
    std::complex<double> alpha, 
    std::complex<double> *a, 
    int32_t lda, 
    std::complex<double> *b, 
    int32_t ldb, 
    std::complex<double> beta, 
    std::complex<double> *c, 
    int32_t ldc) {
    
    if (N == blas_worker)
    {
        const char *trans[] = {"N", "T", "C"};
        FORTFUNC(zgemm)(trans[transa], trans[transb], &m, &n, &k, &alpha, 
            a, &lda, b, &ldb, &beta, c, &ldc, (int32_t)1, (int32_t)1);
    }
    
    if (N == cublas_worker)
    {    
        zgemm_cublas(transa, transb, m, n, k, alpha, a, lda, b, ldb, 
            beta, c, ldc);
    }
}

extern "C" void FORTFUNC(zhemm)(
    const char *side,
    const char *uplo,
    int32_t *m,
    int32_t *n,
    std::complex<double> *alpha,
    std::complex<double> *a,
    int32_t *lda,
    std::complex<double> *b,
    int32_t *ldb,
    std::complex<double> *beta,
    std::complex<double> *c,
    int32_t *ldc,
    int32_t sidelen,
    int32_t uplolen);
    
template <int N> void zhemm(
    int side,
    int uplo,
    int32_t m,
    int32_t n,
    std::complex<double> alpha,
    std::complex<double> *a,
    int32_t lda,
    std::complex<double> *b,
    int32_t ldb,
    std::complex<double> beta,
    std::complex<double> *c,
    int32_t ldc)
{
    if (N == blas_worker) 
    {
        const char *sidestr[] = {"L", "R"};
        const char *uplostr[] = {"U", "L"};
        FORTFUNC(zhemm)(sidestr[side], uplostr[uplo], &m, &n, &alpha, a, &lda, 
            b, &ldb, &beta, c, &ldc, (int32_t)1, (int32_t)1);
    }
}    
    

extern "C" void FORTFUNC(zhegvx)(
    int32_t *itype,
    const char *jobz,
    const char *range,
    const char *uplo,
    int32_t *n,
    std::complex<double> *a,
    int32_t *lda,
    std::complex<double> *b,
    int32_t *ldb,
    double *vl,
    double *vu,
    int32_t *il,
    int32_t *lu,
    double *abstol,
    int32_t *m,
    double *w,
    std::complex<double> *z,
    int32_t *ldz,
    std::complex<double> *work,
    int32_t *lwork,
    double *rwork,
    int32_t *iwork,
    int32_t *ifail,
    int32_t *info,
    int32_t jobzlen,
    int32_t rangelen,
    int32_t uplolen);

extern "C" int32_t FORTFUNC(ilaenv)(
    int32_t *ispec,
    const char *name,
    const char *opts,
    int32_t *n1,
    int32_t *n2,
    int32_t *n3,
    int32_t *n4,
    int32_t namelen,
    int32_t optslen);

template<int N> void zhegv(
  int32_t n,
  int32_t nv,
  double abstol,
  std::complex<double> *a,
  std::complex<double> *b,
  double *eval,
  std::complex<double> *evec,
  int ld)
{
    if (N == lapack_worker) 
    {
        int ispec = 1;
        int n1 = -1;
        int nb = FORTFUNC(ilaenv)(&ispec, "ZHETRD", "U",  &n, &n1, &n1, &n1, (int32_t)6, (int32_t)1);
        int lwork = (nb + 1) * n;
        int *iwork = new int[5 * n];
        int *ifail = new int[n];
        double *w = new double[n];
        double *rwork = new double[7 * n];
        std::complex<double> *work = new std::complex<double>[lwork];
        n1 = 1;
        double vl = 0.0;
        double vu = 0.0;
        int m;
        int info;
        FORTFUNC(zhegvx)(&n1, "V", "I", "U", &n, a, &n, b, &n, &vl, &vu, &n1, 
            &nv, &abstol, &m, w, evec, &ld, work, &lwork, rwork, iwork, ifail,
            &info, (int32_t)1, (int32_t)1, (int32_t)1);
        if (info)
        {
           std::cout << "lapack diagonalisation failed" << std::endl;
           std::cout << "info = " << info << std::endl;
           std::cout << "matrix size = " << n << std::endl;
           exit(0);
        }
        memcpy(eval, w, nv * sizeof(double));
        delete iwork;
        delete ifail;
        delete w;
        delete rwork;
        delete work;
    }
}

extern "C" void FORTFUNC(zgemv)(
    const char *transa,
    int32_t *m,
    int32_t *n,
    std::complex<double> *alpha,
    std::complex<double> *a,
    int32_t *lda,
    std::complex<double> *x,
    int32_t *incx,
    std::complex<double> *beta,
    std::complex<double> *y,
    int32_t *incy,
    int32_t translen);

template<int N> void zgemv(
    int transa,
    int32_t m,
    int32_t n,
    std::complex<double> alpha,
    std::complex<double> *a,
    int32_t lda,
    std::complex<double> *x,
    int32_t incx,
    std::complex<double> beta,
    std::complex<double> *y,
    int32_t incy)
{
    if (N == blas_worker) 
    {
        const char *trans[] = {"N", "T", "C"};
        FORTFUNC(zgemv)(trans[transa], &m, &n, &alpha, a, &lda, x, &incx, 
            &beta, y, &incy, (int32_t)1);
    }
}

#endif // __LINALG_H__

