#include <iostream>
#include <cmath>
#include <cstdlib>

#include "functions.hpp"

double loop_6dim(int N, const double *x, const double *w, 
                 double (*f)(double, double))
{
  double sum = 0;
  for (int x1i = 0; x1i<N; x1i++){
  for (int y1i = 0; y1i<N; y1i++){
  for (int z1i = 0; z1i<N; z1i++){
  for (int x2i = 0; x2i<N; x2i++){
  for (int y2i = 0; y2i<N; y2i++){
  for (int z2i = 0; z2i<N; z2i++){
  const double rdiff_square = square_diff((x[x1i]-x[x2i]), 
           (x[y1i] - x[y2i]), (x[z1i] - x[z2i]));
  if (rdiff_square > 1e-9){ 
    const double r1 = x[x1i]*x[x1i] + x[y1i]*x[y1i] + x[z1i]*x[z1i];
    const double r2 = x[x2i]*x[x2i] + x[y2i]*x[y2i] + x[z2i]*x[z2i];
    sum += f(r1,r2)/std::sqrt(rdiff_square)*
                w[x1i]*w[y1i]*w[z1i]*w[x2i]*w[y2i]*w[z2i];
  }}}}}}}
  return sum;
}

double int_func(double x)
{
  return std::exp(-x);
}

int main()
{
  const int N = 15;
  const double limit = 5;
  
  if (false)
  {
    const double pi = 4*atan(1);
    double analytical = 5*pi*pi/double(16*16);
    std::cout << "I = " << analytical << " Analytical" << std::endl;
  }
  if (true)//brute force
  {
    double *x = new double[N];
    const double delta = 2*limit/(N-1);
    
    const double N_inv = 1/(double)N;
    double *w = new double[N];
    
    for (int ii=0; ii<N; ii++)
    {
      x[ii] = -limit + delta*ii;
      w[ii] = N_inv;
    }
  
    auto brute = [](double r1, double r2){return std::exp(-2*(r1+r2));};
    std::cout << "I = " << loop_6dim(N, x, w, brute) <<
               " Brute force" << std::endl;
    delete[] x;
    delete[] w;
  }
  if (false)
  {
    double *x = new double[N];
    double *w = new double[N];
    
    gauss_laguerre(x, w, N, 2);
    
    auto gauss_lag = [](double r1, double r2){return r1+r2;};
    std::cout << "I = " << loop_6dim(N, x, w, gauss_lag) <<
                " Gauss-Laguerre" << std::endl;
    
    delete[] x;
    delete[] w;
  }
}