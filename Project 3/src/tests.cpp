
#include "catch.hpp"

#include "functions.hpp"




TEST_CASE( "SANITY" )
{
  REQUIRE( 1==1 );
  REQUIRE( 2 !=  2+1 );
}


TEST_CASE( "SQUARE_SUM" , "[functions]")
{
  CHECK(square_sum(1, -6, 9) == Approx ( 1 + 6*6 + 9*9 ) );
}

TEST_CASE( "COS(BETA)" , "[functions]")
{
  const double theta1 = 0.4;
  const double theta2 = -0.2;
  const double phi1 = -1;
  const double phi2 = 0.15;
  
  const double cosB = cos_beta(theta1, theta2, phi1, phi2);
  const double expect = 0.8710982688;
  
  CHECK( cosB == Approx(expect) );
}

TEST_CASE( "GaussLaguerreIntegrateExponential" , "[functions] [GQ]")
{
  //Testing the function xe(-x), -> alpha=1, -> g(x) = 1
  const int N = 2;
  double *x = new double[N];
  double *w = new double[N];
  
  gauss_laguerre(x, w, N, 1);
  
  double sum = 0;
  for (int ii=0; ii<N; ii++)
  {
    sum += w[ii];
  }
  delete[] x;
  delete[] w;
  CHECK( sum == Approx(1) );
}

TEST_CASE( "GaussLaguerreIntegrateExponential2" , "[functions] [GQ]")
{
  //Testing the function 5x^2 e^(-x), -> alpha = 2, -> g(x) = 5
  const int N = 5;
  double *x = new double[N];
  double *w = new double[N];
  gauss_laguerre(x, w, N, 2);
  
  double sum = 0;
  for (int ii=0; ii<N; ii++)
  {
    sum += 5*w[ii];
  }
  delete[] x;
  delete[] w;
  CHECK( sum == Approx(10) );
}

TEST_CASE( "GaussLegendreWeightsMeshpoints", "[functions] [GQ]")
{
  const int N = 2;
  auto x = new double[N];
  auto w = new double[N];
  
  gauss_legendre(x, w, N, -1, 1);
  
  CHECK( x[0] == Approx(-0.57735027).epsilon(1e-5));
  CHECK( x[1] == Approx(0.57735027).epsilon(1e-5));
  
  CHECK( w[0] == Approx(1.0).epsilon(1e-5) );
  CHECK( w[1] == Approx(1.0).epsilon(1e-5) );
  
  delete[] x;
  delete[] w;
}

TEST_CASE( "GaussLegendreIntegratePolynomial", "[functions] [GQ]")
{
  const int N = 3;
  // func to evaluate g(x) = 9*x^5 - 2*x^4 + 10*x^2 - 19.4 from [-1,5]
  auto f = [](double x){return 9*x*x*x*x*x - 2*x*x*x*x + 10*x*x - 19.4;};
  double *x = new double[N];
  double *w = new double[N];
  gauss_legendre(x, w, N, -1, 5);
  
  double sum = 0;
  for (int ii=0; ii<N; ii++)
  {
    sum += w[ii]*f(x[ii]);
  }
  
  CHECK( sum == Approx(22489.2).epsilon(0.1));
}

TEST_CASE( "GaussLegendreLoop", "[functions] [GQ]")
{
  const int dim = 2;
  const double alpha = 2;
  
  const double x[dim] = {-0.57735027, 0.57735027};
  const double w[dim] = {1, 0.5};
  double sum = 0;
  
  for (int ii=0; ii<dim; ii++) {
  for (int jj=0; jj<dim; jj++) {
  for (int kk=0; kk<dim; kk++) {
  for (int ll=0; ll<dim; ll++) {
  for (int mm=0; mm<dim; mm++) {
  for (int nn=0; nn<dim; nn++) {
    const double x12_square = square_sum(x[ii]-x[ll], 
                                           x[jj]-x[mm], x[kk]-x[nn]);
    if (x12_square > tolerance) {
      const double weights = w[ii]*w[jj]*w[kk]*w[ll]*w[mm]*w[nn];
      const double x1 = std::sqrt(square_sum(x[ii], x[jj], x[kk]));
      const double x2 = std::sqrt(square_sum(x[ll], x[mm], x[nn]));
      sum += weights*std::exp(-2*alpha*(x1 + x2))/std::sqrt(x12_square);
    }
  }}}}}}
  
  const double comp = sum_elements_6dim_cartesian(dim, x, w, alpha);
  CHECK( sum == Approx(comp).epsilon(1e-5) );
  
}

TEST_CASE( "GaussLaguerreLoop", "[functions] [GQ]")
{
  const int dim = 2;
  const double r[dim] = {3, 8};
  const double theta[dim] = {-0.5, 0.5};
  const double w[dim] = {1, 0.5};
  
  
  double sum = 0;
  for (int ii=0; ii<dim; ii++) {
  for (int jj=0; jj<dim; jj++) {
  for (int kk=0; kk<dim; kk++) {
  for (int ll=0; ll<dim; ll++) {
  for (int mm=0; mm<dim; mm++) {
  for (int nn=0; nn<dim; nn++) {
    const double cosB = cos_beta(theta[kk], theta[ll], 
                                    theta[mm], theta[nn]);
    const double r_12_square = r[ii]*r[ii] + r[jj]*r[jj] - 2*r[ii]*r[jj]*cosB;
    if (r_12_square > tolerance) {
      const double weights = w[ii]*w[jj]*w[kk]*w[ll]*w[mm]*w[nn];
      sum += std::sin(theta[kk])*std::sin(theta[ll])*weights/std::sqrt(r_12_square);
    }
  }}}}}}
  
  
  
  const double expect = sum_elements_6dim_polar(dim, dim, dim, r, theta, theta, w, w, w);
  
  CHECK( sum == Approx(expect) );
  
  
  
}