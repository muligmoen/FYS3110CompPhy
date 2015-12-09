#include <iostream>
#include "vector.hpp"
#include "diffusion.hpp"

#include <cstring>
#include <chrono>

using namespace diffusion;

int main(int argc, char **argv)
{
  const int N = 10; // Resolution in x-domain
  const double t = 0.1; // Final time
  
  // Matrix things
  const double alpha = 0.5; // stability criteria for forward Euler
  const int Ntimesteps = t*N*N/alpha; // timesteps for matrix thingies
  
  // MC thingies
  const int Nsteps = 2*N*N*t; // MC timesteps
  const int Nparticles = 1000;
  const int MC_repetitions = 10;
  
  const auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  
  
  // Crank Nicolson continous plotting
  if (argc>1 && !std::strcmp(argv[1], "plot")){
    
    auto f_steady = [](double x){return 1.0 - x;};
    
    const auto u_steady = init_vector(0, 1, N, f_steady);
    
    auto vec = not_ends(-u_steady);
    
    std::cout << N << "\n";
    std::cout << Ntimesteps << "\n";
    
    
    const double s = 0.1;
    for (int ii=0; ii<Ntimesteps; ii++){
      vec = Crank_Nicolson(vec, s, 10);
      std::cout << u_steady + add_ends(vec, 0.0, 0.0) << std::endl;
    }
  }
  
  // Crank Nicolson
  if (argc>1 && !std::strcmp(argv[1], "CN")){
    auto f_steady = [](double x){return 1.0 - x;};
    
    const auto u_steady = init_vector(0, 1, N, f_steady);

    auto vec = not_ends(-u_steady); // not including the ends (which are zero)
    
    const double s = 1;
    //const int steps = 100;
    auto v = Crank_Nicolson(vec, s, Ntimesteps);
    
    
    std::cout << u_steady + add_ends(v, 0.0, 0.0) << std::endl;
    
  }
  
  //Analytical solution
  if (argc>1 && !std::strcmp(argv[1], "ANA")){
    
    std::cout << N << std::endl;
    std::cout << Ntimesteps << std::endl;
    
    auto f_steady = [](double x){return 1.0 - x;};
    
    const auto u_steady = init_vector(0, 1, N, f_steady);
    
    const double dt = 0.01;
    for (int ii=0; ii<Ntimesteps; ii++){
      const double t = dt*ii; 
      std::cout << Analytical(t, N) + u_steady << std::endl;
    }
  }
    
  //Monte Carlo equal steplength
  if (argc>1 && !std::strcmp(argv[1], "plotMC")){
    
    std::cout << N << std::endl;
    std::cout << Ntimesteps << std::endl;
  
    const int Nparticles = 1000;
    
    Vector<int> vec(N, [](){return 0;});
    vec[0] = 100;
  

    for (int ii=0; ii<Ntimesteps; ii++){
      vec = Monte_Carlo(vec, 1, Nparticles, seed+ii);
      std::cout << vec << std::endl;
    }
  }
  
  //Monte Carlo equal steplength
  if (argc>1 && !std::strcmp(argv[1], "MC")){

    Vector<int> sum_MC(N);
    //Vector<int> quadrate_sum_MC(N);
    
    for (int ii=0; ii<MC_repetitions; ii++){
      const auto vec = Monte_Carlo(N, Nsteps, Nparticles, seed+ii);
      sum_MC += vec;
      //quadrate_sum_MC += vec*vec;
    }
    
    
    std::cout << normalise(sum_MC) << "\n";
    
    //std::cout << normalise(Monte_Carlo(N, Nsteps, Nparticles, seed)) << std::endl;
  }
  
  //Monte Carlo unequal steplength
  if (argc>1 && !std::strcmp(argv[1], "MCG")){

    const double L0 = 1.0/(double)N;
    const int bins = N;
    //auto sim = Monte_Carlo_gaussian(Nsteps, bins, Nparticles, L0, seed);
    std::cout << normalise(Monte_Carlo_gaussian(Nsteps, bins, 100*Nparticles, L0, seed))
              << std::endl;
  }
  
}
