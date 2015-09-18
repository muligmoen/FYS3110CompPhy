#include <iostream>

#include <armadillo>
#include "Jacobi_rotation.h"

#include "unittest++/UnitTest++.h"

void do_nothing(...) { }

double potential(double r);

arma::Mat<double> get_matrix(int N, double rho_0, double rho_inf, double (*V)(double));

TEST(Potential)
{
  double tolerance = 1e-10;
  double k = 1.0;
  
  CHECK_CLOSE(potential(5), k/(5.0*5.0), tolerance);
}

int main()
{
  UnitTest::RunAllTests();
  /*
  const int N = 5;
  
  double rho_0 = 0.01;
  double rho_inf = 10;
  
  
  arma::Mat<double> A = get_matrix(N, rho_0, rho_inf, potential);
  
  auto eigs =  eig_sym(A);
  
  int MAX_iter = 10;
  for (int iii=0; iii<MAX_iter; iii++)
  {
    int k, l;
    double max_err;
    max_err_offdiag(A, k, l, max_err);
    
    
    if (max_err < 1e-5)
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
  eigs2= sort(eigs2);
  
  for (int iii=0; iii<10; iii++)
  {
    std::cout << eigs(iii) << " " << eigs2(iii) << "\n";
  }*/
}


double potential(double r)
{
  double k = 1.0;
  return k/(r*r);
}

arma::Mat<double> get_matrix(int N, double rho_0, double rho_inf, double (*V)(double))
{
  arma::Mat<double> A(N, N, arma::fill::zeros);
  
  double h = (rho_inf-rho_0)/(double)N;
  double inv_h_square = 1/(h*h);
  A.diag(1) += inv_h_square;
  A.diag(-1) += inv_h_square;
  
  for (int iii=0; iii<N; iii++)
  {
    double rho = rho_0 + iii*h;
    A(iii,iii) = 2*inv_h_square + V(rho);
  }
  
  return A;
}
