#include <iostream>
#include <cstdlib>

#include <armadillo>
#include "Jacobi_rotation.hpp"
#include "filewriter.hpp"
#include "helper_files.hpp"

#include "unittest++/UnitTest++.h"

inline double hydrogen(double r)
{
  return -1/r;  
}

int main(int argc, char *argv[])
{
  UnitTest::RunAllTests();
  std::cout << std::endl;
  int N;
  double rho_0, rho_inf, omega_r;
  bool two_electrons;
  
  rho_0 = 0;
  check_args(argc, argv, N, rho_0, rho_inf, two_electrons, omega_r);
  
  
  //arma::Mat<double> A = ham_matrix(N, rho_0, rho_inf, two_electrons, omega_r);
  
  arma::Mat<double> A = ham_matrix(N, rho_0, rho_inf, hydrogen);
  
  arma::Mat<double> S = identity(N);
  
  rotate_to_diag_with_eigvec(A, S, 1e-8);
  

  auto E = min_three_diag(A);
  
  
  FileWriter file("test.txt");
  file.print(N, rho_0, rho_inf, two_electrons, omega_r);
  file.print(E);
  
  for (int iii=0; iii<3; iii++)
  {
    auto Evec = get_eigv(S, E.indexes[iii]);
    file.print(Evec);
  }
  
  return 0;
}


