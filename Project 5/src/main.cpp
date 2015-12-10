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
  std::cout << u_ANA << "\n";
  
  
  // Forward Euler
  auto v_EF = forward_euler(-u_steady, alpha, Ntimesteps);
  auto u_EF = u_steady + v_EF;
  std::cout << u_EF << "\n" << Error(u_EF, u_ANA) << "\n";
  
  // Backward Euler
  auto v_EB = backward_euler(-u_steady, alpha, Ntimesteps);
  auto u_EB = u_steady + v_EB;
  std::cout << u_EB << "\n" << Error(u_EB, u_ANA) << "\n";
  
  // Crank Nicolson
  auto v_CN = Crank_Nicolson(-u_steady, alpha, Ntimesteps);
  auto u_CN = u_steady + v_CN;
  std::cout << u_CN << "\n" << Error(u_CN, u_ANA) << "\n";
  
  
  //Monte Carlo equal steplength
  Vector<int> sum_MC(N, [](){return 0;});
  Vector<int> quadrate_sum_MC(N, [](){return 0;});
    
  for (int ii=0; ii<MC_repetitions; ii++){
    const auto vec = Monte_Carlo(N, Nsteps, Nparticles, seed+ii);
    sum_MC += vec;
    quadrate_sum_MC += vec*vec;
  }
  
  auto u_MC = normalise(sum_MC, Nparticles*MC_repetitions);
  std::cout << u_MC << "\n" << Error(u_MC, u_ANA) << "\n";
  
  Vector<double> variance(N);
  auto u_MCsquare = normalise(quadrate_sum_MC, Nparticles*Nparticles*MC_repetitions);
  
  for (int ii=0; ii<N; ii++){
    variance[ii] = std::sqrt(u_MCsquare[ii] - u_MC[ii]*u_MC[ii]);
  }
    
  std::cout << variance << "\n";
  
  // Monte Carlo gaussian
  const double L0_factor = 2; // how much Delta_x is shrinked (changes number of timesteps)
                              // this should be a number > 1 to ensure no bins are skipped
  const double L0 = Delta_x/L0_factor;
  
  Vector<int> sum_MCG(N, [](){return 0;});
  Vector<int> quadrate_sum_MCG(N, [](){return 0;});
    
  for (int ii=0; ii<MC_repetitions; ii++){
    const auto vec = Monte_Carlo_gaussian(Nsteps*L0_factor*L0_factor, N, Nparticles, L0, seed+ii);
    sum_MCG += vec;
    quadrate_sum_MCG += vec*vec;
  }
  
  const double norm_factor = (double)sum_MCG[0]/(double)MC_repetitions;
  auto u_MCG = normalise(sum_MCG, norm_factor*MC_repetitions);
  std::cout << u_MCG << "\n" << Error(u_MCG, u_ANA) << "\n";
  
  auto u_MCGsquare = normalise(quadrate_sum_MCG, norm_factor*norm_factor*MC_repetitions);
  
  Vector<double> varianceG(N);
  for (int ii=0; ii<N; ii++){
    varianceG[ii] = std::sqrt(u_MCGsquare[ii] - u_MCG[ii]*u_MCG[ii]);
  }
  
  std::cout << varianceG << "\n";
  
}
