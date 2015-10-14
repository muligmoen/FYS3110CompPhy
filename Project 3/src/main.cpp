#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <ctime>
#include <functional>

#include "functions.hpp"




int main(int argc, char **argv)
{
  auto seed = std::time(nullptr);
  std::default_random_engine gen(seed);
    
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
  if (false) { // gauss legendre
    double *x = new double[N];
    double *w = new double[N];
    
    gauss_legendre(x, w, N, -limit, limit);
    

    std::cout << "I = " << loop_6dim(N, x, w, 2)
              << "\tLegendre" << std::endl;
    
    delete[] x;
    delete[] w;
  }
  if (false) { // brute monte carlo cartesian
    
    const int dim = 6;
    //auto seed = std::time(nullptr);
    //std::default_random_engine gen(seed);
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
  if (false) { // laguerre
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
              
    delete[] r;
    delete[] theta;
    delete[] phi;
    delete[] wr;
    delete[] wtheta;
    delete[] wphi;
  }
  if (false) { // brute monte carlo radial
    std::uniform_real_distribution<double> r_rng(0, limit);
    std::uniform_real_distribution<double> theta_rng(0, pi);
    std::uniform_real_distribution<double> phi_rng(0, 2*pi);
    
    auto r_rand = std::bind(r_rng, gen);
    auto theta_rand = std::bind(theta_rng, gen);
    auto phi_rand = std::bind(phi_rng, gen);
    
    using std::cos;
    using std::sin;
    
    const int dim = 2;
    double r[dim];
    double theta[dim];
    double phi[dim];
    
    double sum = 0;
    
    for (int jj=0; jj<N; jj++) {
      for (int ii = 0; ii< dim; ii++) {
        r[ii] = r_rand();
        theta[ii] = theta_rand();
        phi[ii] = phi_rand();
      }
      
      const double cos_beta = cos(theta[0])*cos(theta[1]) +
                             sin(theta[0])*sin(theta[1])*cos(phi[0]-phi[1]);
                             
      const double r12_square = r[0]*r[0] + r[1]*r[1]
                                  - 2*r[0]*r[1]*cos_beta;
      if (r12_square > tolerance) {
        sum += r[0]*r[0]*r[1]*r[1]*sin(theta[0])*sin(theta[1])/std::sqrt(r12_square);
      }
    }
    
    //sum *= (limit*limit)*(pi*pi)*(2*pi*2*pi);
    sum /= (double)N;
    
    std::cout << "I = " << sum << "\tMonte Carlo radial" << std::endl;
    
    
  }
  if (true) { // monte carlo importance sampling
    
    std::uniform_real_distribution<double> r_rng(0, 1);
    std::uniform_real_distribution<double> theta_rng(-1, 1);
    std::uniform_real_distribution<double> phi_rng(0, 2*pi);
    
    
    auto cos_theta = std::bind(theta_rng, gen);
    auto phi = std::bind(phi_rng, gen);
    auto r_rnd = std::bind(r_rng, gen);
    
    auto cos_beta = [&cos_theta, &phi](){
      double thet1 = cos_theta();
      double thet2 = cos_theta();
      double phi_1 = phi();
      double phi_2 = phi();
      return thet1*thet2+ std::sqrt(1-thet1*thet1)*std::sqrt(1-thet2*thet2)*std::cos(phi_1 - phi_2);
    };
    
    auto r = [&r_rnd](double lambda){
      double x = r_rnd();
      return -lambda*(1.0-x);
    };
    
    const double alpha = 2;
    const double lambda = 1/(double)(2*alpha);
    const double norm_factor = (2*2)*(2*pi*2*pi);
    
    double result = 0;
    for (int ii=0; ii<N; ii++)
    {
      const double r1 = r(lambda);
      const double r2 = r(lambda);
      const double cos_b = cos_beta();
      result += r1*r1*r2*r2/std::sqrt(r1*r1+r2*r2 - 2*r1*r2*cos_b);
    }
    result *= norm_factor;
    result /= N*(2*alpha*2*alpha);
    
    
    std::cout << "I = " << result << "\tMonte Carlo importance" << std::endl;
  }
    
  
}