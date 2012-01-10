#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "constants.h"

#define FORTRAN(x) x##_

const int gemm_worker = blas_worker;

const bool check_evecfv = true;

const bool check_scalar_wf = true;

#endif // __CONFIG_H__
