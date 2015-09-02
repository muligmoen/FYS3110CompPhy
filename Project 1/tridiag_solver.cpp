//#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "tridiag_solver.h"


double u_theory(double x)
{
  return 1.0 - (1.0 - std::exp(-10.0))*x - std::exp(-10.0*x);
}

arma::Col<double> u_theory(const int N, double x0, double x1)
{
  arma::Col<double> u(N);
  double h = 1.0/(double)N;
  for (int iii=0; iii<N; iii++)
  {
    u[iii] = u_theory(iii*h);
  }
  return u;
}

double f(double x)
{
  return 100.0*std::exp(-10.0*x);
}

arma::Col<double> f_column(const double x0, const double x1, const int N)
{
  arma::Col<double> b(N);
  
  double h = (x1-x0)/(double)N;
  
  for (int iii=0; iii<N; iii++)
  {
    double x = x0 + h*iii;
    b[iii] = f(x)*h*h;
  }
  return b;
}

arma::Mat<double> second_deriv_matr(const int&N)
{
  arma::Mat<double> A(N, N);
  A.diag(0) += 2.0;
  A.diag(1) += -1;
  A.diag(-1) += -1;
  return A;
}

arma::SpMat<double> spsecond_deriv_matr(const int& N)
{
  arma::SpMat<double> A(N, N);
  for (int diag=0; diag < N; diag++)
  {
    A(diag, diag) = 2; //diagonal
    if (diag < N - 1) //upper band
    {
      A(diag, diag+1) = -1;
    }
    if (diag > 0) //lower band
    {
      A(diag, diag-1) = -1;
    }
  }
  return A;
}


arma::Col<double> thomas_alg(const double x0, const double x1, const int N)
{
  
  double *cprime = new double[N];
  double *dprime = new double[N];
  
  arma::Col<double> v(N);
  
  const arma::Col<double> b_tilde = f_column(x0, x1, N);
  
  
  const double a = -1;
  const double b = 2;
  const double c = -1;
  
  
  cprime[0] = c/b;
  dprime[0] = b_tilde[0]/b;
  
  for (int iii=1; iii<N; iii++)
  {
    cprime[iii] = c/(b-a*cprime[iii-1]);
    dprime[iii] = (b_tilde[iii] - a*dprime[iii-1])/(b-a*cprime[iii-1]);
  }
  
  v[N] = dprime[N];
  for (int iii=N-1; iii>0; iii--)
  {
    v[iii] = dprime[iii] - cprime[iii]*v[iii+1];
  }
  v[0] = 0;
  
  delete[] cprime;
  delete[] dprime;
  return v;
}


arma::Col<double> matrix_alg(const double x0, const double x1, const int N, const bool SPARSE)
{
  arma::Col<double> f = f_column(0, 1, N);
  
  if (SPARSE)
  {
    arma::SpMat<double> A = spsecond_deriv_matr(N);
    arma::Col<double> v = spsolve(A, f);
    return v;
  } else {
    arma::Mat<double> A = second_deriv_matr(N); 
    arma::Col<double> v = solve(A, f);
    return v;
  }
}


double max_relative_error(const arma::Col<double> &v, const arma::Col<double> &u)
{
  const int N = v.n_elem;
  double max_error = 0;
  for (int iii=0; iii<N; iii++)
  {
    double err = std::abs((v[iii] - u[iii])/u[iii]);
    if (err > max_error && u[iii]) // Checking not null element
    {
      max_error = err;
    }
  }
  return max_error;
}


Writer::Writer(const char* name) : outf(name) 
{
  outf.precision(15);
};

void Writer::print(const arma::Col<double> &vec)
{
  outf << "\n";
  vec.raw_print(outf);
}

void Writer::print(const char *text, const double value)
{
  outf << value << "\t" << text << std::endl;
}

void Writer::print(const char *text, const int value)
{
  outf << value << "\t" << text << std::endl;
}

