#include <iostream>

#include "parma_util.hh"


int main()
{   

    auto generator = PARMA::ParmaGen();
    for (int i=0; i<100; i++)
    {
        auto particle = generator.Generate();
        std::cout<<particle.x<<" "<<particle.y<<" "<<particle.x<<" "<<std::endl;
    }
    return 0;
}