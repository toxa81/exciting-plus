#ifndef _MDARRAY_H_
#define _MDARRAY_H_

#include "constants.h"
#include "mdarray_base.h"
#ifdef _GPU_
#include "mdarray_base_gpu.h"
#else
#include "mdarray_base_cpu.h"
#endif

template <typename T, int ND> class mdarray : public mdarray_base_impl<T,ND> 
{
};

// 1d specialization
template <typename T> class mdarray<T,1> : public mdarray_base_impl<T,1> 
{
    public:
  
        mdarray() 
        {
        }

        mdarray(T *data_ptr, const dimension& j0)
        {
            std::vector<dimension> vd;
            vd.push_back(j0);
            this->init_dimensions(vd);
            this->set_ptr(data_ptr);
        }
    
        inline T& operator()(const int i0) 
        {
            assert(this->mdarray_ptr);
            assert(i0 >= this->d[0].start() && i0 <= this->d[0].end());
            
            int i = this->offset[0] + i0;
            return this->mdarray_ptr[i];
        }
};

// 2d specialization
template <typename T> class mdarray<T,2> : public mdarray_base_impl<T,2> 
{
    public:
  
        mdarray() 
        {
        }

        mdarray(T *data_ptr, const dimension& j0, const dimension& j1)
        {
            std::vector<dimension> vd;
            vd.push_back(j0);
            vd.push_back(j1);
            this->init_dimensions(vd);
            this->set_ptr(data_ptr);
        }
    
        inline T& operator()(const int i0, const int i1) 
        {
            assert(this->mdarray_ptr);
            assert(i0 >= this->d[0].start() && i0 <= this->d[0].end());
            assert(i1 >= this->d[1].start() && i1 <= this->d[1].end());
            
            int i = this->offset[0] + i0 + i1 * this->offset[1];
            return this->mdarray_ptr[i];
        }
};

// 3d specialization
template <typename T> class mdarray<T,3> : public mdarray_base_impl<T,3> 
{
    public:
  
        mdarray() 
        {
        }

        mdarray(T *data_ptr, const dimension& j0, const dimension& j1, const dimension& j2)
        {
            std::vector<dimension> vd;
            vd.push_back(j0);
            vd.push_back(j1);
            vd.push_back(j2);
            this->init_dimensions(vd);
            this->set_ptr(data_ptr);
        }
    
        inline T& operator()(const int i0, const int i1, const int i2) 
        {
            assert(this->mdarray_ptr);
            assert(i0 >= this->d[0].start() && i0 <= this->d[0].end());
            assert(i1 >= this->d[1].start() && i1 <= this->d[1].end());
            assert(i2 >= this->d[2].start() && i2 <= this->d[2].end());
            
            int i = this->offset[0] + i0 + i1 * this->offset[1] + i2 * this->offset[2];
            return this->mdarray_ptr[i];
        }
};

#endif // _MDARRAY_H_

