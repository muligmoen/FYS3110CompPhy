#ifndef functions_h
#define functions_h

#include <cmath>

const double pi = 4*std::atan(1);

inline double square_diff(const double x, const double y)
{
  return (x-y)*(x-y); 
}

inline double square_sum(const double x, const double y, const double z)
{
  return x*x + y*y + z*z;
}

void gauss_laguerre(double *x_return, double *w_return,
                    const int n, const double alf);

void gauss_legendre(double *x, double *w, const int n, 
                    const double x1, const double x2);


double loop_6dim(const int N, const double* x, const double* w, 
                 const double alpha);

double loop_6dim(const int Nr, const int Ntheta, const int Nphi,
                 const double *r, const double *theta, const double *phi,
                 const double *wr, const double *wtheta, const double *wphi,
                 const double alpha);

#endif