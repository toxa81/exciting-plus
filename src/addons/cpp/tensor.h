#ifndef __TENSOR_H__
#define __TENSOR_H__

#include <vector>
#include <cstring>
#include <string>
#include <iostream>
#include <algorithm>
#include <assert.h>

class t_index 
{
    public:
  
        t_index() : start(0), end(-1), size(0) 
        {
        }
    
        t_index(int start, int end) : start(start), end(end) 
        {
            assert(end >= start);
            size = end - start + 1;
        };
    
        int start;
        int end;
        int size;
};
    
template <typename T, int D> class tensor_base 
{
    public:
  
        tensor_base() : data(0), size_(0), allocated(false)
        {
        }
    
        tensor_base(const tensor_base<T,D>& src) 
        {
            copy(src);
            data = new T[size_];
            memcpy(data, src.data, size_ * sizeof(T));
            allocated = true;
        }
    
        ~tensor_base() 
        {
            if (allocated) delete[] data;
        }
    
        void copy(const tensor_base<T,D>& src) 
        {
            size_ = src.size_;
            for (int i = 0; i < D; i++) 
            { 
                offset[i] = src.offset[i];
                index[i] = src.index[i];
            }
        }
    
        void tensor_init(const std::vector<t_index>& vidx) 
        {
            assert(vidx.size() == D);
            
            size_ = 1;
            for (int i = 0; i < D; i++) 
            {
                index[i] = vidx[i];
                size_ *= index[i].size;
            }
            
            offset[0] = -index[0].start;
            int n = 1;
            for (int i = 1; i < D; i++) 
            {
                n *= index[i-1].size;
                offset[i] = n;
                offset[0] -= offset[i] * index[i].start;
            }
        }
    
        tensor_base<T,D>& operator=(tensor_base<T,D> rhs) 
        {
            copy(rhs);
            std::swap(this->data, rhs.data);
            std::swap(this->allocated, rhs.allocated);
            return *this;
        }
    
        int size() 
        {
            return size_;
        }
        
        inline int size(int i)
        {
            return index[i].size;
        } 
    
    protected:
  
        T* data;
        t_index index[D];
        int size_;
        int offset[D];
        bool allocated;
};

template <typename T, int D> class tensor : public tensor_base<T,D> 
{
};

// 1d specialization
template <typename T> class tensor<T,1> : public tensor_base<T,1> 
{
    public:
  
        tensor() 
        {
        }
    
        tensor(T *data_, const int n0) 
        {
            init(data_, t_index(0, n0 - 1));
        }
    
        tensor(T *data_, const t_index& j0) 
        {
            init(data_, j0);
        }
    
        void init(T *data_, const t_index& j0) 
        {
            std::vector<t_index> vidx;
            vidx.push_back(j0);
            this->tensor_init(vidx);
            this->data = data_;
        }
    
        inline T& operator()(const int i0) 
        {
            assert(i0 >= this->index[0].start && i0 <= this->index[0].end);
            
            int i = this->offset[0] + i0;
            return this->data[i];
        } 
};

// 2d specialization
template <typename T> class tensor<T,2> : public tensor_base<T,2> 
{
    public:
        
        tensor() 
        {
        }
        
        tensor(const int n0, const int n1)
        {
            init(0, t_index(0, n0 - 1), t_index(0, n1 - 1));
        }
  
        tensor(T *data_, const int n0, const int n1) 
        {
            init(data_, t_index(0, n0 - 1), t_index(0, n1 - 1));
        }
        
        tensor(const t_index& j0, const t_index& j1) 
        {
            init(0, j0, j1);
        }
    
        tensor(T *data_, const t_index& j0, const t_index& j1) 
        {
            init(data_, j0, j1);
        }
    
        void init(T *data_, const t_index& j0, const t_index& j1) 
        {
            std::vector<t_index> vidx;
            vidx.push_back(j0);
            vidx.push_back(j1);
            this->tensor_init(vidx);
            
            if (data_) 
                this->data = data_;
            else 
            {
                this->data = new T[this->size_];
                memset(this->data, 0, this->size_ * sizeof(T));
                this->allocated = true;
            }
        }
    
        inline T& operator()(const int i0, const int i1) 
        {
            assert(i0 >= this->index[0].start && i0 <= this->index[0].end);
            assert(i1 >= this->index[1].start && i1 <= this->index[1].end);
      
            int i = this->offset[0] + i0 + this->offset[1] * i1;
            return this->data[i];
        } 
};

// 3d specialization
template <typename T> class tensor<T,3> : public tensor_base<T,3> 
{
    public:
    
        tensor() 
        {
        }
  
        tensor(T *data_, const int n0, const int n1, const int n2) 
        {
            init(data_, t_index(0, n0 - 1), t_index(0, n1 - 1), t_index(0, n2 - 1));
        }
    
        tensor(T *data_, const t_index& j0, const t_index& j1, const t_index& j2) 
        {
            init(data_, j0, j1, j2);
        }
    
        void init(T *data_, const t_index& j0, const t_index& j1, const t_index& j2) 
        {
            std::vector<t_index> vidx;
            vidx.push_back(j0);
            vidx.push_back(j1);
            vidx.push_back(j2);
            this->tensor_init(vidx);
            this->data = data_;
        }
    
        inline T& operator()(const int i0, const int i1, const int i2) 
        {
            assert(i0 >= this->index[0].start && i0 <= this->index[0].end);
            assert(i1 >= this->index[1].start && i1 <= this->index[1].end);
            assert(i2 >= this->index[2].start && i2 <= this->index[2].end);
      
            int i = this->offset[0] + i0 + this->offset[1] * i1 + this->offset[2] * i2;
            return this->data[i];
        } 
};

// 4d specialization
template <typename T> class tensor<T,4> : public tensor_base<T,4> 
{
    public:
    
        tensor() 
        {
        }
  
        tensor(T *data_, const int n0, const int n1, const int n2, const int n3) 
        {
            init(data_, t_index(0, n0 - 1), t_index(0, n1 - 1), t_index(0, n2 - 1), t_index(0, n3 - 1));
        }
    
        tensor(T *data_, const t_index& j0, const t_index& j1, const t_index& j2, const t_index& j3) 
        {
            init(data_, j0, j1, j2, j3);
        }
    
        void init(T *data_, const t_index& j0, const t_index& j1, const t_index& j2, const t_index& j3) 
        {
            std::vector<t_index> vidx;
            vidx.push_back(j0);
            vidx.push_back(j1);
            vidx.push_back(j2);
            vidx.push_back(j3);
            this->tensor_init(vidx);
            this->data = data_;
        }
    
        inline T& operator()(const int i0, const int i1, const int i2, const int i3) 
        {
            assert(i0 >= this->index[0].start && i0 <= this->index[0].end);
            assert(i1 >= this->index[1].start && i1 <= this->index[1].end);
            assert(i2 >= this->index[2].start && i2 <= this->index[2].end);
            assert(i3 >= this->index[3].start && i3 <= this->index[3].end);            
      
            int i = this->offset[0] + i0 + this->offset[1] * i1 + this->offset[2] * i2 + this->offset[3] * i3;
            return this->data[i];
        } 
};

// 5d specialization
template <typename T> class tensor<T,5> : public tensor_base<T,5> 
{
    public:
    
        tensor() 
        {
        }
  
        tensor(T *data_, const int n0, const int n1, const int n2, const int n3, const int n4) 
        {
            init(data_, t_index(0, n0 - 1), t_index(0, n1 - 1), t_index(0, n2 - 1), t_index(0, n3 - 1), t_index(0, n4 - 1));
        }
    
        tensor(T *data_, const t_index& j0, const t_index& j1, const t_index& j2, const t_index& j3, const t_index& j4) 
        {
            init(data_, j0, j1, j2, j3, j4);
        }
    
        void init(T *data_, const t_index& j0, const t_index& j1, const t_index& j2, const t_index& j3, const t_index& j4) 
        {
            std::vector<t_index> vidx;
            vidx.push_back(j0);
            vidx.push_back(j1);
            vidx.push_back(j2);
            vidx.push_back(j3);
            vidx.push_back(j4);
            this->tensor_init(vidx);
            this->data = data_;
        }
    
        inline T& operator()(const int i0, const int i1, const int i2, const int i3, const int i4) 
        {
            assert(i0 >= this->index[0].start && i0 <= this->index[0].end);
            assert(i1 >= this->index[1].start && i1 <= this->index[1].end);
            assert(i2 >= this->index[2].start && i2 <= this->index[2].end);
            assert(i3 >= this->index[3].start && i3 <= this->index[3].end);            
            assert(i4 >= this->index[4].start && i4 <= this->index[4].end);            
      
            int i = this->offset[0] + i0 + this->offset[1] * i1 + this->offset[2] * i2 + this->offset[3] * i3 + this->offset[4] * i4;
            return this->data[i];
        } 
};

// 6d specialization
template <typename T> class tensor<T,6> : public tensor_base<T,6> 
{
    public:
    
        tensor() 
        {
        }
  
        tensor(T *data_, const int n0, const int n1, const int n2, const int n3, const int n4, const int n5) 
        {
            init(data_, t_index(0, n0 - 1), t_index(0, n1 - 1), t_index(0, n2 - 1), t_index(0, n3 - 1), t_index(0, n4 - 1), t_index(0, n5 - 1));
        }
    
        tensor(T *data_, const t_index& j0, const t_index& j1, const t_index& j2, const t_index& j3, const t_index& j4, const t_index& j5) 
        {
            init(data_, j0, j1, j2, j3, j4, j5);
        }
    
        void init(T *data_, const t_index& j0, const t_index& j1, const t_index& j2, const t_index& j3, const t_index& j4, const t_index& j5) 
        {
            std::vector<t_index> vidx;
            vidx.push_back(j0);
            vidx.push_back(j1);
            vidx.push_back(j2);
            vidx.push_back(j3);
            vidx.push_back(j4);
            vidx.push_back(j5);
            this->tensor_init(vidx);
            this->data = data_;
        }
    
        inline T& operator()(const int i0, const int i1, const int i2, const int i3, const int i4, const int i5) 
        {
            assert(i0 >= this->index[0].start && i0 <= this->index[0].end);
            assert(i1 >= this->index[1].start && i1 <= this->index[1].end);
            assert(i2 >= this->index[2].start && i2 <= this->index[2].end);
            assert(i3 >= this->index[3].start && i3 <= this->index[3].end);            
            assert(i4 >= this->index[4].start && i4 <= this->index[4].end);            
            assert(i5 >= this->index[5].start && i5 <= this->index[5].end);            
      
            int i = this->offset[0] + i0 + this->offset[1] * i1 + this->offset[2] * i2 + this->offset[3] * i3 + this->offset[4] * i4 + this->offset[5] * i5;
            return this->data[i];
        } 
};



#endif // __TENSOR_H__
