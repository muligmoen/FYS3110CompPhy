#include <iostream>
#include <cmath>
#include <cassert>

#include <armadillo>



void rotate(arma::Mat<double> &A, double c, double s, int k, int l)
{
  for (int iii=0; iii<(int)A.n_rows; iii++)
  {
    if (iii != k && iii != l)
    {
      A(iii,k) = A(iii, k)*c - A(iii, l)*s;
      A(iii, l) = A(iii, l)*c + A(iii, l)*s;
    }
  }
  A(k, k) = A(k, k)*c*c - 2*A(k,l)*c*s + A(l, l)*s*s;
  A(l, l) = A(l, l)*c*c + 2*A(l,k)*c*s + A(k, k)*s*s;
  A(k,l) = 0;
}

void test_rotate()
{
  arma::Mat<double> A;
  A << 1 << 5 << arma::endr
    << 0 << 3;
    

    
  double tau = 1/5.;
  double t = -tau - std::sqrt(1+ 1/(tau*tau));
  double c = 1/std::sqrt(1+t*t);
  double s = t*c;
  std::cout << c << " " << s << std::endl;
  
  rotate(A, c, s, 0, 1);
  
  arma::Mat<double> B;
  B << c*c-5*s*c+3*s*s << 5*c*c-2*s*c << arma::endr
    << -5*s*s-2*s*c << s*s+5*s*c+3*c*c;
    
  std::cout << A << std::endl;
  std::cout << B << std::endl;
}

void max_offdiag_indexes(const arma::Mat<double> &A, int &k, int &l)
{
  double max_error = 0;
  
  for (int iii=0; iii<(int)A.n_rows; iii++)
  {
    for (int jjj=0; jjj<(int)A.n_cols; jjj++)
    {
      if ((jjj != iii) && (std::abs(A(iii,jjj)) > max_error))
      {
	max_error = std::abs(A(iii,jjj));
	k = iii;
	l = jjj;
      }
    }
  }
}

void test_max_indexes()
{
  arma::Mat<double> A;
  A << 1 << 0 << 6 << arma::endr
    << 0 << 3 << 9;
  int k = 0;
  int l = 0;
  max_offdiag_indexes(A, k, l);
  assert(k==1 && l==2);
}

int main()
{
  //test_rotate();
  /*arma::Mat<double> A;
  A << 1 << 5 << 6 << arma::endr
    << 0 << 3 << 7;
  */
  test_max_indexes();
  //int k,l;
  //max_offdiag_indexes(A, k, l);
  //std::cout << k << " " << l << std::endl;
}

