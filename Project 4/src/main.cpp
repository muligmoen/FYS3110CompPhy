#include <iostream>

#include "lattice.hpp"
#include "metropolis.hpp"

int main()
{
  const int lx = 20;
  const int ly = 10;
  Lattice lat(lx, ly, init::random);
  lat.set_print_format(print_t::arrows);
  
  
  const double beta = 1;
  
  double *exp_beta = new double[2];
  pre_compute_exp(exp_beta, beta);
  
  for (int ii = 0; ii<10000000; ii++){
    int x,y;
    get_indexes(x, y, lx, ly);
    try_flip(lat, x, y, exp_beta);
    
    if (ii%1000000==0) {
      std::cout << lat;
      std::cout << "Sum of the spins = " << lat.sum_spins() << std::endl;
      std::cout << "Energy of the lattice = " << lat.energy() << std::endl;
    } 
  }
}