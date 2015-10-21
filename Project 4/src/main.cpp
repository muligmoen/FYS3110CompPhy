#include <iostream>

#include "lattice.hpp"

int main()
{
  Lattice lat(20, 10, init::random);
  lat.set_print_format(print_t::arrows);
  
  
  std::cout << lat;
}