#include <cmath>

#include <omp.h>

#include "functions.hpp"




double gammln(double);

void gauss_laguerre(double *x_return, double *w_return,
                    const int n, const double alf)
{
  double *x = new double[n+1];
  double *w = new double[n+1];
  double p2, p3, pp, z;

  for (int i=1;i<=n;i++) {
    if (i == 1) {
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {
      double ai=i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
              (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
    }
    for (int its=0;its<MAXIT;its++) {
      double p1=1.0;
      //double p2=0.0;
      p2 = 0.0;
      for (int j=1;j<=n;j++) {
        p3=p2;
        p2=p1;
        p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      pp=(n*p1-(n+alf)*p2)/z;
      double z1 = z;
      z = z1 - p1/pp;
      if (std::abs(z-z1) <= EPS) break;
    }
    x[i] = z;
    w[i] = -std::exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
  }
  for (int ii=0; ii<n; ii++)
  {
    x_return[ii] = x[ii+1];
    w_return[ii] = w[ii+1];
  }
  delete[] x;
  delete[] w;
}


double gammln( double xx)
{
  const double cof[6]={76.18009172947146,    -86.50532032941677,
                       24.01409824083091,    -1.231739572450155,
                       0.1208650973866179e-2,-0.5395239384953e-5};

  double y = xx;
  double x = xx;
  double tmp = x + 5.5;
  tmp -= (x+0.5)*std::log(tmp);
  double ser=1.000000000190015;
  for (int j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+std::log(2.5066282746310005*ser/x);
}



void gauss_legendre(double *x, double *w, const int n,
                    const double x1, const double x2)
{
  for(int ii = 1; ii <= (n+1)/2; ii++) 
  { 
    double z = std::cos(pi * (ii - 0.25)/(n + 0.5));
    double z1;
    double pp;
    do {
      double p1 = 1.0;
      double p2 = 0.0;


      for(int jj = 1; jj <= n; jj++) 
      {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * jj - 1.0) * z * p2 - (jj - 1.0) * p3)/jj;
      }

      pp = n * (z * p1 - p2)/(z * z - 1.0);
      z1 = z;
      z  = z1 - p1/pp;                   // Newton's method
    } while(std::abs(z - z1) > ZERO);
    
    double xm = 0.5 * (x2 + x1); // change of variables
    double xl = 0.5 * (x2 - x1);
  
    x[ii - 1] = xm - xl *z;
    x[n - ii] = xm + xl*z;
      
    w[ii - 1] = 2.0 * xl/((1.0 - z * z) * pp * pp);
    w[n - ii] = w[ii-1];
  }
}




double sum_elements_6dim_cartesian(const int N, const double *x, const double *w, 
                 const double alpha)
{
  double sum = 0;
  for (int ii = 0; ii<N; ii++){
  for (int jj = 0; jj<N; jj++){
  for (int kk = 0; kk<N; kk++){
  for (int ll = 0; ll<N; ll++){
  for (int mm = 0; mm<N; mm++){
  for (int nn = 0; nn<N; nn++){
    const double rdiff_sum_square = square_sum((x[ii]-x[ll]), 
                            (x[jj] - x[mm]), (x[kk] - x[nn]));
    if (rdiff_sum_square > tolerance){ 
      const double r1 = std::sqrt(square_sum(x[ii], x[jj], x[kk]));
      const double r2 = std::sqrt(square_sum(x[ll], x[mm], x[nn]));
      const double weigths = w[ii]*w[jj]*w[kk]*w[ll]*w[mm]*w[nn];
      sum += std::exp(-2*alpha*(r1 + r2))/std::sqrt(rdiff_sum_square)*weigths;
  }}}}}}}
  return sum;
}



double sum_elements_6dim_polar(const int Nx, const int Ntheta, const int Nphi,
                 const double *x, const double *theta, const double *phi,
                 const double *wr, const double *wtheta, const double *wphi)
{

  double sum = 0;
  for (int ii = 0; ii<Nx; ii++){
  for (int jj = 0; jj<Nx; jj++){
  for (int kk = 0; kk<Ntheta; kk++){
  for (int ll = 0; ll<Ntheta; ll++){
  for (int mm = 0; mm<Nphi; mm++){
  for (int nn = 0; nn<Nphi; nn++){
    const double cosbeta = cos_beta(theta[kk], theta[ll],
                                       phi[mm], phi[nn]);
    
    const double x12_square = x[ii]*x[ii] + x[jj]*x[jj]
                                     - 2*x[ii]*x[jj]*cosbeta;
    if (x12_square > tolerance) {
      const double weights = wr[ii]*wr[jj]*wtheta[kk]*wtheta[ll]*wphi[mm]*wphi[nn];
      sum += weights*std::sin(theta[kk])*std::sin(theta[ll])/std::sqrt(x12_square);
    }
  }}}}}}
  //sum *= 2; // magic factor
  
  return sum;
}


std::_Bind<std::uniform_real_distribution<double>(std::default_random_engine)> 
            uniform_distribution(const double lower, const double upper)
/*
 * Binds a generator and a uniform distribution 
 * to a new function which gives the number in the uniform 
 * range from lower to upper
 */
{
  static const auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  static std::default_random_engine gen(seed);
  
  std::uniform_real_distribution<double> dist(lower, upper);
  auto uniform = std::bind(dist, gen);
  return uniform;
}
