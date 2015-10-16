#ifndef functions_h
#define functions_h

#include <cmath>
#include <functional>
#include <random>
#include <chrono>



const double pi = 4*std::atan(1);
const int MAXIT = 10;
const double EPS = 3e-14;
const double ZERO = 1.0E-10;
const double tolerance = 1e-9;




inline double square_sum(const double x, const double y, const double z)
{
  return x*x + y*y + z*z;
}

void gauss_laguerre(double *x_return, double *w_return,
                    const int n, const double alf);

void gauss_legendre(double *x, double *w, const int n, 
                    const double x1, const double x2);


double cartesian_loop(const int N, const double* x, const double* w, 
                 const double alpha);

double polar_loop(const int Nr, const int Ntheta, const int Nphi,
                 const double *r, const double *theta, const double *phi,
                 const double *wr, const double *wtheta, const double *wphi,
                 const double alpha);

// uniform number generator, could be replaced by a simple function or lambda
// function to give a generator with the ranges specified
std::_Bind<std::uniform_real_distribution<double>(std::default_random_engine)> 
            uniform_distribution(const double lower, const double upper);
            

#endif