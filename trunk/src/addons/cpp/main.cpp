#include <iostream>
#include <vector>
#include "tensor.h"

int main(int argc, char **argv) 
{
    tensor<std::vector<int>,2> t2;
    t2=tensor<std::vector<int>,2>(3,4);
    tensor<int,2> ti2(10,20);

    return 0;
}
