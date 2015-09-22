#include <iostream>
#include <cstdlib>

#include <armadillo>
#include "Jacobi_rotation.hpp"

#include "unittest++/UnitTest++.h"



arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf, double (*V)(double));

inline double potential(double rho);


int main(int argc, char *argv[])
{
  UnitTest::RunAllTests();
  
  int N = -1;
  int rho_inf;
  
  if (argc < 3)
  {
    N = 100;
    rho_inf = 10;
    
    std::cout << "Running standard settings" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "rho_inf = " << rho_inf << std::endl;
  }
  else 
  {
    N = std::atoi(argv[1]);
    rho_inf = std::atof(argv[2]);
    if (N<1 || rho_inf<=0)
    {
      std::cerr << "Invalid arguments, exiting" << std::endl;
      std::exit(1);
    }
  }

  double rho_0 = 0.0;
  
  arma::Mat<double> A = ham_matrix(N, rho_0, rho_inf, potential);
  
  double max_err = 100;
  
  while (max_err > 1e-10)
  {
    int k, l;
    max_err_offdiag(A, k, l, max_err);
    
    double cos, sin;
    find_cos_sin(A(k,k), A(l,l), A(k,l), cos, sin);
    
    rotate(A, cos, sin, k, l);
  }
  
  arma::vec eigs = A.diag();
  
  arma::vec eig = sort(eigs);
  
  for (int iii=0; iii<std::min(10, N); iii++)
  {
    std::cout << eig(iii) << std::endl;
  }
  
  return 0;
}



arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf, double (*V)(double))
{
  arma::Mat<double> A(N, N, arma::fill::zeros);
  
  double h = (rho_inf-rho_0)/(double)N;
  double inv_h_square = 1/(h*h);
  A.diag(1) -= inv_h_square;
  A.diag(-1) -= inv_h_square;
  
  for (int iii=0; iii<N; iii++)
  {
    double rho = rho_0 + iii*h;
    A(iii,iii) = 2*inv_h_square + V(rho);
  }
  
  return A;
}

inline double potential(double rho)
{
  return rho*rho;
}
