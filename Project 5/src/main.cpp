#include <iostream>
#include "vector.hpp"
#include "diffusion.hpp"

#include <cstring>
#include <chrono>
#include <cstdlib>
#include <cmath>

using namespace diffusion;

int main(int argc, char **argv)
{ 
  if (argc<4){
    std::cout << "Usage: " << argv[0] << " Delta_t Delta_x T" << std::endl;
    std::exit(1);
  }
  const double Delta_t = std::strtod(argv[1], nullptr);
  const double Delta_x = std::strtod(argv[2], nullptr);
  const double T = std::strtod(argv[3], nullptr);

  
  const int N = 1/Delta_x; // Resolution in x-domain
  const int Ntimesteps = T/Delta_t;
  
  
  // Matrix things
  const double alpha = Delta_t/(Delta_x*Delta_x);
  auto f_steady = [](double x){return 1.0 - x;};
  const auto u_steady = init_vector(0, 1, N, f_steady);
  
  // MC thingies
  const int Nsteps = 2*N*N*T; // MC timesteps
  const int Nparticles = 1000;
  const int MC_repetitions = 10;
  
  const auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  
  //Analytical solution
  auto v_ANA = Analytical(T, N);
  auto u_ANA = u_steady + v_ANA;
  std::cout <<  u_ANA << "\n";
  
  
  // Forward Euler
  auto v_EF = forward_euler(-u_steady, alpha, Ntimesteps);
  auto u_EF = u_steady + v_EF;
  std::cout << u_EF << " " << Error(u_EF, u_ANA) << "\n";
  
  // Backward Euler
  auto v_EB = backward_euler(-u_steady, alpha, Ntimesteps);
  std::cout << u_steady + v_EB << "\n";
  
  // Crank Nicolson
  auto v_CN = Crank_Nicolson(-u_steady, alpha, Ntimesteps);
  std::cout << u_steady + v_CN << "\n";
  
  
  
  //Monte Carlo equal steplength
  Vector<int> sum_MC(N, [](){return 0;});
  Vector<int> quadrate_sum_MC(N, [](){return 0;});
    
  for (int ii=0; ii<MC_repetitions; ii++){
    const auto vec = Monte_Carlo(N, Nsteps, Nparticles, seed+ii);
    sum_MC += vec;
    quadrate_sum_MC += vec*vec;
  }
  
  std::cout << normalise(sum_MC, sum_MC[0]) << "\n";
  
  Vector<double> variance(N);
  for (int ii=0; ii<N; ii++){
    variance[ii] = std::sqrt(sum_MC[ii]*sum_MC[ii] - quadrate_sum_MC[ii]);
  }
    
  //std::cout << variance << "\n";
  
  /*
  //
    const double L0 = Delta_x;
    const int bins = N;
    //auto sim = Monte_Carlo_gaussian(Nsteps, bins, Nparticles, L0, seed);
    std::cout << normalise(Monte_Carlo_gaussian(Nsteps, bins, 100*Nparticles, L0, seed), 0.0)
              << std::endl;
  }
  */
  
}
