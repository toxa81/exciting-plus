#ifndef _MDARRAY_BASE_GPU_H_
#define _MDARRAY_BASE_GPU_H_

#include <iostream>

template <typename T, int ND> class mdarray_base_impl : public mdarray_base<T,ND> 
{
    public:
        
        void copy_to_device() 
        {
            std::cout << "GPU copy to device" << std::endl;
        }
};

#endif // _MDARRAY_BASE_GPU_H_


