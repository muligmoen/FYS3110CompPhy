#include <iostream>
#include <cstring>
#include <cstdlib>
#include <chrono>

#include "lattice.hpp"
#include "ising.hpp"

int main(int argc, char **argv)
{
  
  if ((argc > 4) && (!std::strcmp(argv[1],"-P"))) { // Python called functions
    const int N = std::atoi(argv[3]);
    const int L = std::atoi(argv[4]);
    const double beta = std::atof(argv[5]);
    const long int seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    Ising model(L, seed, beta);
    model.init_rand();
    
    if (!std::strcmp(argv[2],"M")) { // Magnetisation
      for (int ii=0; ii<N; ii++){
        std::cout << model.get_magnetisation() << " ";
        model.try_flip();
      }
    }
    if (!std::strcmp(argv[2],"E")){ //Energy
      for (int ii=0; ii<N; ii++){
        std::cout << model.get_energy() << " ";
        model.try_flip();
      }
    }
    
  }
}