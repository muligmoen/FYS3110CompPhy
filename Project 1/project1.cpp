#include <iostream>
#include <armadillo>
#include <cmath>
#include <cstdlib>

#include "tridiag_solver.h"


int main( int argc, char *argv[] )
{
  if (argc < 3 || atoi(argv[1]) < 2)
  {
    std::cout << "Usage: " << argv[0] << " N filename" << std::endl;
    std::cout << "With N > 1" << std::endl;
    exit(1);
  }
  
  const int N = atoi(argv[1]);
  
  arma::Col<double> u(N);
  u_theory(u);
  

  
  arma::Col<double> v = matrix_alg(0, 1, N);
  
  std::cout << std::log10(max_relative_error(v, u)) << std::endl;
  
  Writer fileprinter(argv[2], N);
  fileprinter.print(v);
  
  fileprinter.print(u);
  
  return 0;
}


