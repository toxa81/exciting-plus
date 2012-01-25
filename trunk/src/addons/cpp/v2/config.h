#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "constants.h"
#include "typedefs.h"

#define FORTRAN(x) x##_

//const int gemm_worker = blas_worker;

const bool check_evecfv = false;

const bool check_scalar_wf = false;

//const bool use_gpu = (gemm_worker == cublas_worker) || (lapack_worker == magma_worker);

const implementation lapw_impl = cpu;

#endif // __CONFIG_H__
