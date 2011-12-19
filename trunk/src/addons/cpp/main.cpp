#include <iostream>
#include "tensor.h"

int main(int argc, char** argv) {
  tensor<int,2> lmidx(50,100);
  
  
  for (int i=0; i<50; i++) {
   for (int j=0; j<100; j++) {
     std::cout << lmidx(i,j) << std::endl;
   }
  }
  
  //int intgv_[] = { -14,         -14,         -26,          15,         15,          27};
  //tensor<int,1> ias2is(intgv_, 6);
  //tensor<int,1> ias2is1 = tensor<int,1>(intgv_, 6);
  //tensor<int,1> ias2is2;
  //ias2is2 = ias2is;
  //ias2is1 = ias2is2;
  //tensor<int,1> ias2is3 = tensor<int,1>(intgv_, 6);
 // std::cout << ias2is1(0) << std::endl;
  //ias2is1(0) -= 1;
  //std::cout << ias2is1(0) << std::endl;
  
  

  //std::cout << "size = " << ias2is1.size() << std::endl;
  //std::cout << "size = " << ias2is2.size() << std::endl;
  //std::cout << "element = " << ias2is(5) << std::endl;
  //std::cout << "element = " << ias2is1(5) << std::endl;
  //ias2is2(5) -= 1;
  //std::cout << "element = " << ias2is2(5) << std::endl;

  //std::cout << ias2is2.size() << std::endl;
  //tensor<int,1> intgv(intgv_,3);
  //std::cout << intgv.size() << std::endl;
  //ias2is = intgv;
  //std::cout << intgv(0) << " " << intgv(1) << std::endl;
  //std::cout << intgv.offset[0] << " " << intgv.offset[1] << std::endl;
  //return 0;
}
