#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <ctime>
#include <functional>

#include "functions.hpp"




int main(int argc, char **argv)
{
  int N;
  double limit;
  if (argc < 3) {
    N = 15;
    limit = 5;
  } else {
    N = std::atoi(argv[1]);
    limit = std::atof(argv[2]);
  }

  
  if (true) { // analytical
    double analytical = 5*pi*pi/double(16*16);
    std::cout << "I = " << analytical << "\tAnalytical" << std::endl;
  }
  if (false) { // brute brute force
    double *x = new double[N];
    const double delta = 2*limit/(N-1);
    
    const double N_inv = 1/(double)N;
    double *w = new double[N];
    
    for (int ii=0; ii<N; ii++)
    {
      x[ii] = -limit + delta*ii;
      w[ii] = N_inv;
    }
  
    
    std::cout << "I = " << loop_6dim(N, x, w, 2) <<
               "\tBrute force" << std::endl;
    delete[] x;
    delete[] w;
  }
  if (false) { // gauss legendre
    double *x = new double[N];
    double *w = new double[N];
    
    gauss_legendre(x, w, N, -limit, limit);
    

    std::cout << "I = " << loop_6dim(N, x, w, 2)
              << "\tLegendre" << std::endl;
    
    delete[] x;
    delete[] w;
  }
  if (false) { // brute monte carlo
    
    const int dim = 6;
    auto seed = std::time(nullptr);
    std::default_random_engine gen(seed);
    std::uniform_real_distribution<double> RNG(-limit, limit);
    auto get_num = std::bind(RNG, gen);
    
    double x[dim];
    double result = 0;
    
    for (int jj=0; jj<N; jj++) {
      for (int ii=0; ii<dim; ii++) {
        x[ii] = get_num();
      }
      
      const double rdiff_sum_square = square_sum( x[0] - x[3], 
                                    x[1] - x[4], x[2] - x[5] );
      if (rdiff_sum_square > 1e-9){ 
        const double r1 = std::sqrt(square_sum(x[0], x[1], x[2]));
        const double r2 = std::sqrt(square_sum(x[3], x[4], x[5]));
      
        result += std::exp(-4*(r1 + r2))/std::sqrt(rdiff_sum_square);
      }
    }
    
    result *= 64*limit*limit*limit*limit*limit*limit; // normation from 
                                                      // change of limits
    result /= (double)N;
    
    std::cout << "I = " << result
              << "\tBrute Monte Carlo" << std::endl;

  }
  if (true) { // laguerre
    double *r = new double[N];
    double *wr = new double[N];
    gauss_laguerre(r, wr, N, 2);
    
    double *theta = new double[N];
    double *wtheta = new double[N];
    gauss_legendre(theta, wtheta, N, 0, pi);
    
    double *phi = new double[N];
    double *wphi = new double[N];
    gauss_legendre(phi, wphi, N, 0, 2*pi);
    
    double result = loop_6dim(N, N, N, r, theta, phi, wr, wtheta, wphi, 2);
    
    std::cout << "I = " << result << "\tGauss-Laguerre and Gauss-Legendre"
              << std::endl;
  }
    
    
  
}