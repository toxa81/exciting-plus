#include <iostream>
#include <vector>
#include "tensor.h"
#include "timer.h"

int main(int argc, char **argv) 
{
    //tensor<std::vector<int>,2> t2;
    //t2=tensor<std::vector<int>,2>(3,4);
    //tensor<int,2> ti2(10,20);
    //for (int i = 1; i > 0; i--)
    //    std::cout << i << std::endl;
    
    timer *t = new timer("t1");
    double d = 0;
    for (int i = 0; i< 100; i++)
      for (int j=0; j<100; j++)
        for (int k=0; k<100000; k++)
          d += 1.0*i*j/double(k+1);
    
    std::cout << d << std::endl;
    
    delete t;
    
    timer::print("t1");
       
    return 0;
}
