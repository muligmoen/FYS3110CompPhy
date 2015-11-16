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
  
  init_vec[0] = 0;
  init_vec[N-1] = 0;
  
  const double s = 4;
  auto vec = Crank_Nicolson(init_vec, s, 100);
  
  std::cout << vec << std::endl;
}
