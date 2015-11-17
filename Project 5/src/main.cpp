#include <iostream>
#include "vector.hpp"
#include "diffusion.hpp"

using namespace diffusion;

int main()
{
  const int N = 8;
  
  auto f_steady = [](double x){return 1 - x;};
  
  auto u_steady = init_vector(0, 1, N, f_steady);
  
  auto init_vec = -u_steady;

  auto vec = init_vec;
  
  const double s = 0.3;
  for (int ii=0; ii<100; ii++){
    vec = Crank_Nicolson(vec, s, 1);
  
    std::cout << vec << std::endl;
  }
}
