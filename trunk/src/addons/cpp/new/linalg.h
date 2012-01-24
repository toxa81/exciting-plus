#ifndef __LINALG_H__
#define __LINALG_H__

#include <stdint.h>
#include <string.h>
#include <iostream>
#include <complex>
#include <vector>
#include <cstdlib>
#include "config.h"

typedef std::complex<double> complex16;

extern "C" void FORTRAN(zgemm)(const char *transa, const char *transb, int32_t *m, int32_t *n, 
                               int32_t *k, complex16 *alpha, complex16 *a, int32_t *lda, 
                               complex16 *b, int32_t *ldb, complex16 *beta, complex16 *c, 
                               int32_t *ldc, int32_t transalen, int32_t transblen);

extern "C" void zgemm_cublas(int transa, int transb, int32_t m, int32_t n, int32_t k, 
                             complex16 alpha, complex16 *a, int32_t lda, complex16 *b, 
                             int32_t ldb, complex16 beta, complex16 *c, int32_t ldc);
    
template<int N> void zgemm(int transa, int transb, int32_t m, int32_t n, int32_t k, complex16 alpha, 
                           complex16 *a, int32_t lda, complex16 *b, int32_t ldb, complex16 beta, 
                           complex16 *c, int32_t ldc) 
{
    
    if (N == blas_worker)
    {
        const char *trans[] = {"N", "T", "C"};
        FORTRAN(zgemm)(trans[transa], trans[transb], &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, 
            c, &ldc, (int32_t)1, (int32_t)1);
    }
    
    if (N == cublas_worker)
    {    
        zgemm_cublas(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

extern "C" void FORTRAN(zhemm)(const char *side, const char *uplo, int32_t *m, int32_t *n, 
                               complex16 *alpha, complex16 *a, int32_t *lda, complex16 *b,
                               int32_t *ldb, complex16 *beta, complex16 *c, int32_t *ldc,
                               int32_t sidelen, int32_t uplolen);
    
template <int N> void zhemm(int side, int uplo, int32_t m, int32_t n, complex16 alpha,
                            complex16 *a, int32_t lda, complex16 *b, int32_t ldb, complex16 beta,
                            complex16 *c, int32_t ldc)
{
    if (N == blas_worker) 
    {
        const char *sidestr[] = {"L", "R"};
        const char *uplostr[] = {"U", "L"};
        FORTRAN(zhemm)(sidestr[side], uplostr[uplo], &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, 
            &ldc, (int32_t)1, (int32_t)1);
    }
}    
    

extern "C" void FORTRAN(zhegvx)(int32_t *itype, const char *jobz, const char *range, 
                                const char *uplo, int32_t *n, complex16 *a, int32_t *lda,
                                complex16 *b, int32_t *ldb, double *vl, double *vu, int32_t *il,
                                int32_t *iu, double *abstol, int32_t *m, double *w, complex16 *z,
                                int32_t *ldz, complex16 *work, int32_t *lwork, double *rwork,
                                int32_t *iwork, int32_t *ifail, int32_t *info, int32_t jobzlen,
                                int32_t rangelen, int32_t uplolen);

extern "C" int32_t FORTRAN(ilaenv)(int32_t *ispec, const char *name, const char *opts, int32_t *n1,
                                   int32_t *n2, int32_t *n3, int32_t *n4, int32_t namelen, 
                                   int32_t optslen);

template<int N> void zhegv(int32_t n, int32_t nv, double abstol, complex16 *a, complex16 *b,
                           double *eval, complex16 *z, int32_t ldz)
{
    if (N == lapack_worker) 
    {
        int ispec = 1;
        int n1 = -1;
        int nb = FORTRAN(ilaenv)(&ispec, "ZHETRD", "U",  &n, &n1, &n1, &n1, (int32_t)6, (int32_t)1);
        int lwork = (nb + 1) * n;
        std::vector<int> iwork(5 * n);
        std::vector<int> ifail(n);
        std::vector<double> w(n);
        std::vector<double> rwork(7 * n);
        std::vector< std::complex<double> > work(lwork);
        n1 = 1;
        double vl = 0.0;
        double vu = 0.0;
        int m;
        int info;
        FORTRAN(zhegvx)(&n1, "V", "I", "U", &n, a, &n, b, &n, &vl, &vu, &n1, 
            &nv, &abstol, &m, &w[0], z, &ldz, &work[0], &lwork, &rwork[0], 
            &iwork[0], &ifail[0], &info, (int32_t)1, (int32_t)1, (int32_t)1);
        if (info)
        {
           std::cout << "zhegvx diagonalization failed" << std::endl
                     << "info = " << info << std::endl
                     << "matrix size = " << n << std::endl;
           exit(0);
        }
        memcpy(eval, &w[0], nv * sizeof(double));
    }
}

extern "C" void FORTRAN(zheev)(const char *jobz, const char *uplo, int32_t *n, complex16 *a,
                               int32_t *lda, double *w, double *work, int32_t *lwork, double *rwork,
                               int32_t *info, int32_t jobzlen, int32_t uplolen);

template<int N> void zheev(int32_t n, complex16 *a, int32_t lda, double *eval)
{
    if (N == lapack_worker)
    {
        int ispec = 1;
        int n1 = -1;
        int nb = FORTRAN(ilaenv)(&ispec, "ZHETRD", "U",  &n, &n1, &n1, &n1, (int32_t)6, (int32_t)1);
        int lwork = (nb + 1) * n;
        std::vector<double> work(lwork * 2);
        std::vector<double> rwork(3* n + 2);
        int info;
        FORTRAN(zheev)("V", "U", &n, a, &lda, eval, &work[0], &lwork, &rwork[0], &info, (int32_t)1, (int32_t)1);
        if (info)
        {
           std::cout << "zheev diagonalization failed" << std::endl
                     << "info = " << info << std::endl
                     << "matrix size = " << n << std::endl;
           exit(0);
        }
    }
}

extern "C" void FORTRAN(zcopy)(int32_t *n, complex16 *zx, int32_t *incx, complex16 *zy, int32_t *incy);

void zcopy(int32_t n, complex16 *zx, int32_t incx, complex16 *zy, int32_t incy);

/*
extern "C" void FORTRAN(zgemv)(const char *transa,
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

template<int N> void zgemv(int transa,
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
*/
#endif // __LINALG_H__

