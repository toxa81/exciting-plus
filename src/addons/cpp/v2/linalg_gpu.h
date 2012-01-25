#ifndef __LINALG_GPU_H__
#define __LINALG_GPU_H__

#include <stdint.h>
#include <string.h>
#include <iostream>
#include <complex>
#include <vector>
#include <cstdlib>
#include "config.h"

extern "C" void zgemm_cublas(int transa, int transb, int32_t m, int32_t n, int32_t k, 
                             complex16 alpha, complex16 *a, int32_t lda, complex16 *b, 
                             int32_t ldb, complex16 beta, complex16 *c, int32_t ldc);
    
#endif // __LINALG_GPU_H__

