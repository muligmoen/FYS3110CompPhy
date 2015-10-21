#include <iostream>

#include "lattice.hpp"

int main()
{
  Lattice lat(20,10, init::random);
  
  std::cout << lat;
}