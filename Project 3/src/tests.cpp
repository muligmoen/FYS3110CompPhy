#include <iostream>

#include "catch.hpp"
#include "functions.hpp"



//const double tolerance = 1e-10;

TEST_CASE( "SANITY" )
{
  REQUIRE( 1==1 );
  REQUIRE( 2 !=  2+1 );
}

TEST_CASE( "SQUARE_DIFF" , "[functions]")
{
  CHECK( square_diff(5, 2) == Approx( 9.0 ) );
}

TEST_CASE( "SQUARE_SUM" , "[functions]")
{
  CHECK(square_sum(1, -6, 9) == Approx ( 1 + 6*6 + 9*9 ) );
}


TEST_CASE( "GaussLaguerre" , "[functions] [Gauss]")
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

TEST_CASE( "GaussLaguerre2" , "[functions] [Gauss]")
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

TEST_CASE( "GaussLegendre", "[functions] [Gauss]")
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

TEST_CASE( "GaussLegendreIntegrate", "[functions] [Gauss]")
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
