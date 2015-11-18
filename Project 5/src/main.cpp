#include <iostream>
#include "vector.hpp"
#include "diffusion.hpp"

using namespace diffusion;

int main()
{
  
  //const int N = 800;
  //const int Nt = 20000;
  /*
  auto f_steady = [](double x){return 1 - x;};
  
  const auto u_steady = init_vector(0, 1, N, f_steady);
  
  const auto init_vec = -u_steady;

  auto vec = init_vec;
  
  std::cout << N << std::endl;
  std::cout << Nt << std::endl;
  
  const double s = 3;
  for (int ii=0; ii<Nt; ii++){
    vec = Crank_Nicolson(vec, s, 1);
    std::cout << u_steady + vec << std::endl;
  }
  */
  //std::cout << N << std::endl;
  //std::cout << Nt << std::endl;
  
  //Vector<int> vec(N, [](){return 0;});
  //vec[0] = 100;
  
  std::cout << Monte_Carlo_gaussian(100000, 10, 1000, 0.01) << std::endl;

  //for (int ii=0; ii<Nt; ii++){
  //  vec = Monte_Carlo(vec, 1, 100);
  //  std::cout << vec << std::endl;
  //}
}
