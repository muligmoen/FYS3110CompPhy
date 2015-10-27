#include <iostream>
#include <ctime>

#include "lattice.hpp"
//#include "metropolis.hpp"
#include "ising.hpp"

int main()
{

  const int L = 20;
  
  
  const int seed = std::clock();
  const double beta = 1;
  Ising model(L, seed, beta);
  model.set_print_format(print_t::arrows);
  
  
  std::cout << model.get_magnetisation() << " " << model.recompute_magnetisation() << std::endl;
  model.try_flip(1000);
  std::cout << model.get_magnetisation() << " " << model.recompute_magnetisation() << std::endl;
}