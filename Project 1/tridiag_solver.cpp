#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "tridiag_solver.h"




arma::Col<double> f_column(double x0, double x1, int N, double (*f)(double), bool squared_h)
{
  arma::Col<double> b(N);
  
  double h = (x1-x0)/(double)N;
  
  if (squared_h)
  {
    for (int iii=0; iii<N; iii++)
    {
      double x = x0 + h*iii;
      b[iii] = f(x)*h*h;
    }
    return b;
  } else {
    for (int iii=0; iii<N; iii++)
    {
      double x = x0 + h*iii;
      b[iii] = f(x);
    }
    return b;
  }
}

arma::Mat<double> second_deriv_matr(int N)
{
  arma::Mat<double> A(N, N);
  A.diag(0) += -2.0;
  A.diag(1) += 1;
  A.diag(-1) += 1;
  return A;
}

arma::SpMat<double> spsecond_deriv_matr(int N)
{
  arma::SpMat<double> A(N, N);
  for (int diag=0; diag < N; diag++)
  {
    A(diag, diag) = -2; //diagonal
    if (diag < N - 1) //upper band
    {
      A(diag, diag+1) = 1;
    }
    if (diag > 0) //lower band
    {
      A(diag, diag-1) = 1;
    }
  }
  return A;
}


arma::Col<double> matrix_alg(double x0, double x1, int N, double (*f)(double), bool SPARSE)
{
  arma::Col<double> f_col = -f_column(0, 1, N, f, true);
  
  if (SPARSE)
  {
    arma::SpMat<double> A = spsecond_deriv_matr(N);
    arma::Col<double> v = spsolve(A, f_col);
    return v;
  } else {
    arma::Mat<double> A = second_deriv_matr(N); 
    arma::Col<double> v = solve(A, f_col);
    return v;
  }
}

void thomas_alg(double *v, double x0, double x1, int N, double (*f)(double))
{
  
  double cprime[N];
  double dprime[N];
  
  arma::Col<double> b_tilde = f_column(x0, x1, N, f, true);
  
  // the diagonal band does not change
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

double max_relative_error(const double* u, const arma::Col<double>& v)
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

void Writer::print(const char *text, double value)
{
  outf << value << "\t" << text << std::endl;
}

void Writer::print(const char *text, int value)
{
  outf << value << "\t" << text << std::endl;
}

void Writer::print(const arma::Col<double> &vec)
{
  outf << "\n";
  vec.raw_print(outf);
}

void Writer::print(int N, const double *vec)
{
  outf << "\n";
  for (int iii=0; iii<N; iii++)
  {
    outf << vec[iii] << "\n";
  }
}

void Writer::newline()
{
  outf << "\n";
}
