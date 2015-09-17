 
#include <armadillo>

#include "Jacobi_rotation.h"
#include "unittest++/UnitTest++.h"

SUITE(Jacobi)
{
  double tolerance = 1e-10;
  TEST(Rotate)
  {
    double off_number = 1/std::sqrt(2);
    arma::Mat<double> A;
    A << 1 << off_number << arma::endr
      << off_number << 3;
      
    double tau = (3-1)/((int)2*off_number);
    double t = -tau - std::sqrt(1+ tau*tau);
    double c = 1/std::sqrt(1+t*t);
    double s = t*c;
    
    
    arma::Mat<double> S;
    S << c << s << arma::endr
      << -s << c;
    
    arma::Mat<double> B = S.t()*A*S;
    rotate(A, c, s, 0, 1);
    
    
    for (int iii=0; iii<2; iii++)
    {
      for (int jjj=0; jjj<2; jjj++)
      {
        CHECK_CLOSE(A(iii,jjj), B(iii,jjj), tolerance);
      }
    }
  }
  TEST(max_index)
  {
    arma::Mat<double> A;
    A << 1 << 0 << 6 << arma::endr
      << 0 << 3 << 9;
    int k = 0;
    int l = 0;
    double err = 0;
    
    max_err_offdiag(A, k, l, err);
    
    CHECK_EQUAL(k, 1);
    CHECK_EQUAL(l, 2);
    CHECK_EQUAL(err, 9);
  }
  TEST(sin_cos)
  {
    
   /* A = [1        off_number]
    *     [off_number        3]
    */
   
    double off_number = 1/std::sqrt(2);
    double tau = (3-1)/((int)2*off_number);
    double t = -tau - std::sqrt(1+ tau*tau);
    double c = 1/std::sqrt(1+t*t);
    double s = t*c;

    
    double x,y;
    find_cos_sin(1, 3, off_number, x, y);
    
    CHECK_CLOSE(x, c, tolerance);
    CHECK_CLOSE(y, s, tolerance);
  }
}