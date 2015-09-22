#include <iostream>
#include <cstdlib>

#include <armadillo>
#include "Jacobi_rotation.hpp"

#include "unittest++/UnitTest++.h"



arma::Mat<double> get_matrix(int N, double rho_0, double rho_inf);


int main(int argc, char *argv[])
{
  UnitTest::RunAllTests();
  
  int N = -1;
  int MAX_iter = -1;
  
  if (argc < 3)
  {
    N = 100;
    MAX_iter = 10000;
    std::cout << "Running standard settings" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "Max iterations = " << MAX_iter << std::endl;
  }
  else 
  {
    N = std::atoi(argv[1]);
    MAX_iter = std::atof(argv[2]);
    if (N<1 || MAX_iter<1)
    {
      std::cerr << "Invalid arguments, exiting" << std::endl;
      std::exit(1);
    }
  }

  double rho_0 = 0.0;
  double rho_inf = 10;
  
  
  arma::Mat<double> A = get_matrix(N, rho_0, rho_inf);
  

  auto eigs =  eig_sym(A);

 
  for (int iii=0; iii<MAX_iter; iii++)
  {
    int k, l;
    double max_err;
    max_err_offdiag(A, k, l, max_err);
    
    
    if (max_err < 1e-10)
    {
      break;
    }
    
    double cos, sin;
    find_cos_sin(A(k,k), A(l,l), A(k,l), cos, sin);
    
    rotate(A, cos, sin, k, l);

  }
  
  arma::Col<double> eigs2(N);
 
  for (int iii=0; iii<N; iii++)
  {
    eigs2(iii) = A(iii,iii);
  }
  eigs2 = sort(eigs2);
  
  
  for (int iii=0; iii<std::min(10, N); iii++)
  {
    std::cout << eigs(iii) << " " << eigs2(iii) << "\n";
  }
  return 0;
}



arma::Mat<double> get_matrix(int N, double rho_0, double rho_inf)
{
  arma::Mat<double> A(N, N, arma::fill::zeros);
  
  double h = (rho_inf-rho_0)/(double)N;
  double inv_h_square = 1/(h*h);
  A.diag(1) -= inv_h_square;
  A.diag(-1) -= inv_h_square;
  
  for (int iii=0; iii<N; iii++)
  {
    double rho = rho_0 + iii*h;
    A(iii,iii) = 2*inv_h_square + rho*rho;
  }
  
  return A;
}
