#include <iostream>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "tridiag_solver.h"


int main( int argc, char *argv[] )
{
  if (argc < 4 || atoi(argv[1]) < 2)
  {
    std::cout << "Usage: " << argv[0] << " N filename <type>" << std::endl;
    std::cout << "With N > 1, and <type> as M for matrix, S for sparse matrix\n" <<
     "and T for Thomas algorithm" << std::endl;
    exit(1);
  }
  
  const int N = atoi(argv[1]);
  
  arma::Col<double> u(N);
  u_theory(u);
  

  std::clock_t t0, t1;
  
  
  arma::Col<double> v;
  
  t0 = std::clock();
  const char alg_type = argv[3][0];
  if (alg_type == 'M')
  {
    v = matrix_alg(0, 1, N, false);
  } else if (alg_type == 'S') {
    v = matrix_alg(0, 1, N);
  } else if (alg_type == 'T') {
    v = thomas_alg(0, 1, N);
  } else {
    std::cout << "Please use one of the supported algorithms" << std::endl;
    exit(1);
  }
  t1 = std::clock();
  
  
  auto time_ticks = t1-t0;
  //std::cout << "Number of clock cycles used: " << time_ticks << std::endl;
  double time_used = time_ticks/(double)CLOCKS_PER_SEC;
  std::cout << "Time used: " << time_used << " seconds." << std::endl;
  
  double err = std::log10(max_relative_error(v, u));
  std::cout << "Log10 of relative error:" << err << std::endl;
  
  Writer fileprinter(argv[2]);
  
  fileprinter.print("N", N);
  fileprinter.print("Time used", time_used);
  fileprinter.print("max(log10(rel_error))", err);

  fileprinter.print(v);
  
  fileprinter.print(u);
  
  return 0;
}


