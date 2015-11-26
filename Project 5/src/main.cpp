#include <iostream>
#include "vector.hpp"
#include "diffusion.hpp"

#include <cstring>

using namespace diffusion;

int main(int argc, char **argv)
{
  if (argc>1 && !std::strcmp(argv[1], "CN")){
    const int N = 8;
    const int Nt = 50;
    
    auto f_steady = [](double x){return 1.0 - x;};
    
    const auto u_steady = init_vector(0, 1, N, f_steady);
    
    const auto init_vec = -u_steady;

    auto vec = init_vec;
    
    std::cout << N << std::endl;
    std::cout << Nt << std::endl;
    
    const double s = 1;
    for (int ii=0; ii<Nt; ii++){
      vec = Crank_Nicolson(vec, s, 1);
      std::cout << u_steady + vec << std::endl;
    }
  }
  
  if (argc>1 && !std::strcmp(argv[1], "ANA")){
    const int N = 100;
    const int Nt = 30;
    
    std::cout << N << std::endl;
    std::cout << Nt << std::endl;
    
    auto f_steady = [](double x){return 1.0 - x;};
    
    const auto u_steady = init_vector(0, 1, N, f_steady);
    
    const double dt = 0.01;
    for (int ii=0; ii<Nt; ii++){
      const double t = dt*ii; 
      std::cout << Analytical(t, N) + u_steady << std::endl;
    }
  }
    
    
  
  //std::cout << N << std::endl;
  //std::cout << Nt << std::endl;
  
  //Vector<int> vec(N, [](){return 0;});
  //vec[0] = 100;
  
  //std::cout << Monte_Carlo_gaussian(100000, 10, 1000, 0.01) << std::endl;

  //for (int ii=0; ii<Nt; ii++){
  //  vec = Monte_Carlo(vec, 1, 100);
  //  std::cout << vec << std::endl;
  //}
}
