#include <iostream>
#include <cmath>
#include <cassert>

#include <armadillo>

void find_cos_sin(double a_kk, double a_ll, double a_kl, double &c, double &s)
{
  double tau = (a_ll - a_kk)/(2*a_kl);
  double t = tau > 0 ? -tau -std::sqrt(1+tau*tau) : tau -std::sqrt(1+tau*tau);
  
  c = 1/std::sqrt(1+t*t);
  s = t*c;  
}

void rotate(arma::Mat<double> &A, double c, double s, int k, int l)
{
  double a_kk = A(k, k);
  double a_ll = A(l, l);
  double a_lk = A(l, k);
  double a_kl = A(k, l);
  
  for (int iii=0; iii<(int)A.n_rows; iii++)
  {
    A(iii, k) = A(iii, k)*c - A(iii, l)*s;
    A(iii, l) = A(iii, l)*c + A(iii, l)*s;
  }
  
  A(k, k) = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s;
  A(l, l) = a_ll*c*c + 2*a_lk*c*s + a_kk*s*s;
  A(k,l) = (a_kk-a_ll)*c*s + a_kl*(c*c-s*s); // DEBUG
  A(l,k) = (a_kk-a_ll)*c*s + a_kl*(c*c-s*s); // DEBUG
}

void test_rotate()
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
  
  std::cout << A-B << std::endl;
}

void test_find_cos_sin()
{
  double off_number = 1/std::sqrt(2);
  /*
   * A = [1        off_number]
   *     [off_number        3]
   */
  
  double tau = (3-1)/((int)2*off_number);
  double t = -tau - std::sqrt(1+ tau*tau);
  double c = 1/std::sqrt(1+t*t);
  double s = t*c;

  
  double x,y;
  find_cos_sin(1, 3, off_number, x, y);
  
  std::cout << x << " " << c << "\n"
     << y << " " << s << std::endl;
  
}

void max_offdiag_indexes(const arma::Mat<double> &A, int &k, int &l, double &err)
{
  err = 0;
  
  for (int iii=0; iii<(int)A.n_rows; iii++)
  {
    for (int jjj=0; jjj<(int)A.n_cols; jjj++)
    {
      if ((jjj != iii) && (std::abs(A(iii,jjj)) > err))
      {
	err = std::abs(A(iii,jjj));
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
  double err = 0;
  
  max_offdiag_indexes(A, k, l, err);
  assert(k==1 && l==2);
}

int main()
{
  test_rotate();
  test_max_indexes();
  test_find_cos_sin();
}

