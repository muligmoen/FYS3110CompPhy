#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <omp.h>

#include "functions.hpp"

            

void process_args(const int argc, char **argv, int &N, double &limit, 
                  bool *methods);


int main(int argc, char **argv)
{
  int N;
  double limit;
  bool method[6] = {}; // all init to false
  process_args(argc, argv, N, limit, method);
  
  const double alpha = 2;
  
  

  
  if (method[0]) { // analytical
    double analytical = 5*pi*pi/double(16*16);
    std::cout << "I = " << analytical << "\tAnalytical" << std::endl;
  }
  if (method[1]) { // Gauss Quadrature Legendre
    double *x = new double[N];
    double *w = new double[N];
    
    gauss_legendre(x, w, N, -limit, limit);
    

    std::cout << "I = " << sum_elements_6dim_cartesian(N, x, w, 2)
              << "\tLegendre" << std::endl;
    
    delete[] x;
    delete[] w;
  }
  if (method[2]) { // Gauss Quadrature Laguerre
    double *r = new double[N];
    double *wr = new double[N];
    gauss_laguerre(r, wr, N, 2);
    
    double *theta = new double[N];
    double *wtheta = new double[N];
    gauss_legendre(theta, wtheta, N, 0, pi);
    
    double *phi = new double[N];
    double *wphi = new double[N];
    gauss_legendre(phi, wphi, N, 0, 2*pi);
    
    double result = sum_elements_6dim_polar(N, N, N, r, theta, phi, wr, wtheta, wphi);
    
    result /=  std::pow(2*alpha, 5);
    
    std::cout << "I = " << result << "\tGauss-Laguerre and Gauss-Legendre"
              << std::endl;
              
    delete[] r;
    delete[] theta;
    delete[] phi;
    delete[] wr;
    delete[] wtheta;
    delete[] wphi;
  }
  if (method[3]) { // brute monte carlo cartesian
    
    const int dim = 6;
    
    auto get_num = uniform_distribution(-limit, limit);
    
    double x[dim];
    double sum = 0;
    double sum_squares = 0;
    
    #pragma omp parallel for private(x) reduction(+:sum,sum_squares)
    for (int jj=0; jj<N; jj++) {
      for (int ii=0; ii<dim; ii++) {
        x[ii] = get_num();
      }
      
      const double rdiff_sum_square = square_sum( x[0] - x[3], 
                                    x[1] - x[4], x[2] - x[5] );
      if (rdiff_sum_square > tolerance){ 
        const double r1 = std::sqrt(square_sum(x[0], x[1], x[2]));
        const double r2 = std::sqrt(square_sum(x[3], x[4], x[5]));
      
        const double part_sum = std::exp(-2*alpha*(r1 + r2))/std::sqrt(rdiff_sum_square);
        sum += part_sum;
        sum_squares += part_sum*part_sum;
      }
    }
    
    const double norm_factor = std::pow(2*limit, 6); // normation from change of limits
    sum *= norm_factor;
    sum_squares *= norm_factor*norm_factor;
    
    
    const double mean = sum/(double)N;
    const double variance = sum_squares/(double)N - mean*mean;
    
    
    std::cout << "I = " << mean << "\tsigma = " << std::sqrt(variance/(double)N)
              << "\tMonte Carlo cartesian" << std::endl;

  }
  if (method[4]) { // brute monte carlo radial
    
    auto r_rand = uniform_distribution(0, limit);
    auto theta_rand = uniform_distribution(0, pi);
    auto phi_rand = uniform_distribution(0, 2*pi);
    
    
    const int dim = 2;
    double r[dim];
    double theta[dim];
    double phi[dim];
    
    double sum = 0;
    double sum_squares = 0;
    
    #pragma omp parallel for private(r, theta, phi) reduction(+:sum,sum_squares)
    for (int jj=0; jj<N; jj++) {
      for (int ii = 0; ii< dim; ii++) {
        r[ii] = r_rand();
        theta[ii] = theta_rand();
        phi[ii] = phi_rand();
      }
      
      const double cosbeta = cos_beta(theta[0], theta[1], phi[0], phi[1]);
                             
      const double r12_square = r[0]*r[0] + r[1]*r[1]
                                  - 2*r[0]*r[1]*cosbeta;
      if (r12_square > tolerance) {
        const double dV = r[0]*r[0]*r[1]*r[1]*std::sin(theta[0])*std::sin(theta[1]);
        const double part_sum = std::exp(-2*alpha*(r[0]+r[1]))/std::sqrt(r12_square)*dV;
        sum += part_sum;
        sum_squares += part_sum*part_sum;
      }
    }
    const double norm_factor = (pi*pi)*(2*pi*2*pi)*(limit*limit);
    sum *= norm_factor;
    sum_squares *= norm_factor*norm_factor;
    
    const double mean = sum/(double)N;
    const double variance = sum_squares/(double)N  - mean*mean;
    
    
    std::cout << "I = " << mean << "\tsigma = " << std::sqrt(variance/(double)N)
              << "\tMonte Carlo polar" << std::endl;
    
    
  }
  if (method[5]) { // monte carlo importance sampling
    
    auto cos_theta = uniform_distribution(-1, 1);
    auto phi = uniform_distribution(0, 2*pi);
    
    // mapping cos(beta) from uniform distributions
    auto cos_beta = [&cos_theta, &phi](){
      double x1 = cos_theta();
      double x2 = cos_theta();
      double phi_1 = phi();
      double phi_2 = phi();
      return x1*x2 + std::sqrt(1-x1*x1)*std::sqrt(1-x2*x2)*std::cos(phi_1 - phi_2);
    };
    
    
    const double lambda = 1.0/(2*alpha);
    
    // Mapping to exponential distribution with a lambda func
    //auto x = uniform_distribution(0, 1);
    //auto r = [&x, lambda](){return -lambda*std::log(1-x());};
    auto r = exponential_distribution(lambda);
    
    
    double sum = 0;
    double sum_squares = 0;
    
    #pragma omp parallel for reduction(+:sum,sum_squares)
    for (int ii=0; ii<N; ii++)
    {
      const double r1 = r();
      const double r2 = r();
      const double cos_b = cos_beta();
      const double dV = r1*r1*r2*r2;
      const double r12_square = r1*r1 + r2*r2 - 2*r1*r2*cos_b;
      if (r12_square > tolerance) {
        const double part_sum = dV/std::sqrt(r12_square);
        sum += part_sum;
        sum_squares += part_sum*part_sum;
      }
    }
    
    const double norm_factor = (2*2)*(2*pi*2*pi)*(lambda*lambda);
    sum *= norm_factor;
    sum_squares *= norm_factor*norm_factor;
    
    
    const double mean = sum/(double)N;
    const double variance = sum_squares/(double)N - mean*mean;
    
    
    std::cout << "I = " << mean <<  "\tsigma = " << std::sqrt(variance/(double)N)
              << "\tMonte Carlo importance" << std::endl;
  }
    
  
}


void process_args(const int argc, char **argv, int &N, double &limit, 
                  bool *methods)
{
  if (argc < 2) {
    std::cerr << "Usage : " << argv[0] << " N <-l lim> <methods>\n" <<
                 "Where the limit can be specified with the -l flag\n\n" << 
                 "Methods available are \n" << 
                 "ANA : analytical solution\n" <<
                 "GLE : gauss legendre quadrature\n" <<
                 "GLA : gauss laguerre quadrature\n" <<
                 "MCC : Monte Carlo with cartesian coordinates\n" << 
                 "MCP : Monte Carlo with polar coordinates\n" << 
                 "MCI : Monte Carlo with importance sampling" << std::endl;
    std::exit(1);
  } else {
    N = std::atof(argv[1]);
  }
  
  limit = 5;
  
  for (int SearchI=0; SearchI < argc; SearchI++)
  {
    if (!std::strcmp(argv[SearchI], "-l")) limit = std::atof(argv[++SearchI]);
    if (!std::strcmp(argv[SearchI], "ANA")) methods[0] = true;
    if (!std::strcmp(argv[SearchI], "GLE")) methods[1] = true;
    if (!std::strcmp(argv[SearchI], "GLA")) methods[2] = true;
    if (!std::strcmp(argv[SearchI], "MCC")) methods[3] = true;
    if (!std::strcmp(argv[SearchI], "MCP")) methods[4] = true;
    if (!std::strcmp(argv[SearchI], "MCI")) methods[5] = true;
  }
  
}


