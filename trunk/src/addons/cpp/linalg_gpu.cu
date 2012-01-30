#include "gpu_interface.h"
#include "linalg_gpu.h"

extern "C" void gpu_zgemms(int transa, int transb, int32_t m, int32_t n, int32_t k, 
                           complex16 alpha, complex16 *a, int32_t lda, complex16 *b, 
                           int32_t ldb, complex16 beta, complex16 *c, int32_t ldc)
{
    assert(sizeof(cuDoubleComplex) == sizeof(complex16));
    
    const cublasOperation_t trans[] = {CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C};

    if (cublasZgemm(cublas_handle(), trans[transa], trans[transb], m, n, k, (cuDoubleComplex *)(&alpha), (cuDoubleComplex *)a, lda, 
                    (cuDoubleComplex *)b, ldb, (cuDoubleComplex *)(&beta), (cuDoubleComplex *)c, ldc) != CUBLAS_STATUS_SUCCESS)
    {
        printf("failed to execute cublasZgemm() \n");
        exit(0);
    }
}

extern "C" void gpu_zhegvx(int32_t n, int32_t nv, double abstol, void *a, void *b,
                           double *eval, void *z, int32_t ldz)
{
    magma_int_t m1, info;

    magma_int_t nb = magma_get_zhetrd_nb(n);
    magma_int_t lwork = 2 * n * (nb + 1);
    magma_int_t lrwork = 7 * n;
    magma_int_t liwork = 6 * n;
    
    cuDoubleComplex *h_work;
    double *rwork, *w1;
    magma_int_t *iwork, *ifail;
    
    w1 = (double *)malloc(n * sizeof(double));
    h_work = (cuDoubleComplex *)malloc(lwork * sizeof(cuDoubleComplex));
    rwork = (double *)malloc(lrwork * sizeof(double));
    iwork = (magma_int_t *)malloc(liwork * sizeof(magma_int_t));
    ifail = iwork + 5 * n;

    magma_zhegvx(1, 'V', 'I', 'U', n, (cuDoubleComplex *)a, n, (cuDoubleComplex *)b, n, 0.0, 0.0, 1, nv, abstol, 
                 &m1, w1, (cuDoubleComplex *)z, ldz, h_work, lwork, rwork, iwork, ifail, &info);

    memcpy(eval, &w1[0], nv * sizeof(double)); 
    
    free(iwork);
    free(rwork);
    free(w1);
    free(h_work);
}
 
