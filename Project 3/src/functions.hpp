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

inline double cos_beta(const double theta1, const double theta2,
                       const double phi1, const double phi2)
{
   return std::cos(theta1)*std::cos(theta2) +
            std::sin(theta1)*std::sin(theta2)*std::cos(phi1 - phi2);
}


inline double square_sum(const double x, const double y, const double z)
{
  return x*x + y*y + z*z;
}

void gauss_laguerre(double *x_return, double *w_return,
                    const int n, const double alf);

void gauss_legendre(double *x, double *w, const int n, 
                    const double x1, const double x2);


double sum_elements_6dim_cartesian(const int N, const double* x, const double* w, 
                 const double alpha);



double sum_elements_6dim_polar(const int Nr, const int Ntheta, const int Nphi,
                 const double *r, const double *theta, const double *phi,
                 const double *wr, const double *wtheta, const double *wphi);

// uniform number generator, could be replaced by a simple function or lambda
// function to give a generator with the ranges specified
std::_Bind<std::uniform_real_distribution<double>(std::mt19937_64)> 
            uniform_distribution(const double lower, const double upper);
            
std::_Bind<std::exponential_distribution<double>(std::mt19937_64)> 
            exponential_distribution(const double lambda);
            

#endif