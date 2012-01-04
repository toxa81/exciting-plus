#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "constants.h"

#define FORTFUNC(x) x##_

const int gemm_worker = blas_worker;

const bool check_evecfv = false; 

#endif // __CONFIG_H__
