#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <omp.h>

#include "functions.hpp"


// uniform number generator
auto uniform_distribution(const double lower, const double upper)
{
  static const auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  static std::default_random_engine gen(seed);
  
  std::uniform_real_distribution<double> dist(lower, upper);
  auto uniform = std::bind(dist, gen);
  return uniform;
}

void process_args(int argc, char **argv, int &N, double &limit, 
                  bool *met);



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
  if (method[1]) { // gauss legendre
    double *x = new double[N];
    double *w = new double[N];
    
    gauss_legendre(x, w, N, -limit, limit);
    

    std::cout << "I = " << cartesian_loop(N, x, w, 2)
              << "\tLegendre" << std::endl;
    
    delete[] x;
    delete[] w;
  }
  if (method[2]) { // laguerre
    double *r = new double[N];
    double *wr = new double[N];
    gauss_laguerre(r, wr, N, 2);
    
    double *theta = new double[N];
    double *wtheta = new double[N];
    gauss_legendre(theta, wtheta, N, 0, pi);
    
    double *phi = new double[N];
    double *wphi = new double[N];
    gauss_legendre(phi, wphi, N, 0, 2*pi);
    
    double result = polar_loop(N, N, N, r, theta, phi, wr, wtheta, wphi, 2);
    
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
    double result = 0;
    
    #pragma omp parallel for private(x) reduction(+:result)
    for (int jj=0; jj<N; jj++) {
      for (int ii=0; ii<dim; ii++) {
        x[ii] = get_num();
      }
      
      const double rdiff_sum_square = square_sum( x[0] - x[3], 
                                    x[1] - x[4], x[2] - x[5] );
      if (rdiff_sum_square > 1e-9){ 
        const double r1 = std::sqrt(square_sum(x[0], x[1], x[2]));
        const double r2 = std::sqrt(square_sum(x[3], x[4], x[5]));
      
        result += std::exp(-2*alpha*(r1 + r2))/std::sqrt(rdiff_sum_square);
      }
    }
    
    result *= 64*limit*limit*limit*limit*limit*limit; // normation from 
                                                      // change of limits
    result /= (double)N;
    
    std::cout << "I = " << result
              << "\tBrute Monte Carlo" << std::endl;

  }
  if (method[4]) { // brute monte carlo radial
    
    auto r_rand = uniform_distribution(0, limit);
    auto theta_rand = uniform_distribution(0, pi);
    auto phi_rand = uniform_distribution(0, 2*pi);
    
    using std::cos;
    using std::sin;
    
    const int dim = 2;
    double r[dim];
    double theta[dim];
    double phi[dim];
    
    double sum = 0;
    
    #pragma omp parallel for private(r, theta, phi) reduction(+:sum)
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
        const double dV = r[0]*r[0]*r[1]*r[1]*sin(theta[0])*sin(theta[1]);
        sum += std::exp(-2*alpha*(r[0]+r[1]))/std::sqrt(r12_square)*dV;
      }
    }
    
    sum *= (pi*pi)*(2*pi*2*pi)*(limit*limit);
    sum /= (double)N;
    
    std::cout << "I = " << sum << "\tMonte Carlo radial" << std::endl;
    
    
  }
  if (method[5]) { // monte carlo importance sampling
    
    auto cos_theta = uniform_distribution(-1, 1);
    auto phi = uniform_distribution(0, 2*pi);
    
    // mapping cos(beta) from uniform distributions
    auto cos_beta = [&cos_theta, &phi](){
      double thet1 = cos_theta();
      double thet2 = cos_theta();
      double phi_1 = phi();
      double phi_2 = phi();
      return thet1*thet2 + std::sqrt(1-thet1*thet1)*std::sqrt(1-thet2*thet2)*std::cos(phi_1 - phi_2);
    };
    
    
    const double alpha = 2;
    const double lambda = 1/(double)(2*alpha);
    
    // Mapping to exponential distribution with a lambda func
    auto x = uniform_distribution(0, 1);
    auto r = [&x, lambda](){return -lambda*std::log(1-x());};
    
    const double norm_factor = (2*2)*(2*pi*2*pi)*(lambda*lambda);
    
    double result = 0;
    
    #pragma omp parallel for reduction(+:result)
    for (int ii=0; ii<N; ii++)
    {
      const double r1 = r();
      const double r2 = r();
      const double cos_b = cos_beta();
      result += r1*r1*r2*r2/std::sqrt(r1*r1 + r2*r2 - 2*r1*r2*cos_b);
    }
    result *= norm_factor;
    
    
    result /= N;
    
    
    std::cout << "I = " << result << "\tMonte Carlo importance" << std::endl;
  }
    
  
}


void process_args(int argc, char **argv, int &N, double &limit, 
                  bool *met)
{
  if (argc < 2) {
    std::cerr << "Usage : " << argv[0] << " N <-l lim> <methods>\n" <<
                 "Where the limit can be specified with the -l flag\n\n" << 
                 "Methods available are \n" << 
                 "ANA -> analytical solution\n" <<
                 "GLE -> gauss legendre quadrature\n" <<
                 "GLA -> gauss laguerre quadrature\n" <<
                 "MCC -> Monte Carlo with cartesian coordinates\n" << 
                 "MCR -> Monte Carlo with spherical coordinates\n" << 
                 "MCI -> Monte Carlo with importance sampling" << std::endl;
    std::exit(1);
  } else {
    N = std::atof(argv[1]);
  }
  
  limit = 5;
  
  for (int SearchI=0; SearchI < argc; SearchI++)
  {
    if (!std::strcmp(argv[SearchI], "-l")) limit = std::atof(argv[SearchI+1]);
    if (!std::strcmp(argv[SearchI], "ANA")) met[0] = true;
    if (!std::strcmp(argv[SearchI], "GLE")) met[1] = true;
    if (!std::strcmp(argv[SearchI], "GLA")) met[2] = true;
    if (!std::strcmp(argv[SearchI], "MCC")) met[3] = true;
    if (!std::strcmp(argv[SearchI], "MCR")) met[4] = true;
    if (!std::strcmp(argv[SearchI], "MCI")) met[5] = true;
  }
  
}