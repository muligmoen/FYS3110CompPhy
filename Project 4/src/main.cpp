#include <iostream>

#include "lattice.hpp"
#include "metropolis.hpp"

int main()
{
  const int Nsteps = 10000;
  const int plot_every = Nsteps/30;
  const int lx = 20;
  const int ly = 20;
  Lattice lat(lx, ly, init::random);
  lat.set_print_format(print_t::numbers);
  
  const double beta = 1;
  
  std::cout << "Parameters begin\n";
  std::cout << lx << " Lx\n";
  std::cout << ly << " Ly\n";
  std::cout << beta << " beta\n";
  std::cout << "Initial state = \n";
  std::cout << lat;
  std::cout << "Parameters end\n";
  
  double exp_beta[2];
  
  
  int easy_flips = 0;
  int hard_flips = 0;
  int not_flips = 0;
  
  for (int ii = 0; ii<Nsteps; ii++){
    int x,y;
    get_indexes(x, y, lx, ly);
    const int code = try_flip(lat, x, y, exp_beta);
    
    if (code==FlipCodes::FLIPPED) easy_flips++;
    if (code==FlipCodes::RAND_FLIPPED) hard_flips++;
    if (code==FlipCodes::NOT_FLIPPED) not_flips++;
    
    if (ii%plot_every==0) {
      std::cout << lat << "\n";
      //std::cout << "Sum of the spins = " << lat.sum_spins() << std::endl;
      //std::cout << "Energy of the lattice = " << lat.energy() << std::endl;
    } 
  }
  std::cout << "Simulation end";
  
  //std::cout << "Less energy flips = " << easy_flips << std::endl;
  //std::cout << "Temperature flips = " << hard_flips << std::endl;
  //std::cout << "Not flipped = " << not_flips << std::endl;
  
}