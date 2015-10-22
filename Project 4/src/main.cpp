#include <iostream>

#include "lattice.hpp"

int main()
{
  Lattice lat(20, 10, init::random);
  lat.set_print_format(print_t::arrows);
  
  
  std::cout << lat;
  
  std::cout << "Sum of the spins = " << lat.sum_spins() << std::endl;
  std::cout << "Energy of the lattice = " << lat.energy(1) << std::endl;
}