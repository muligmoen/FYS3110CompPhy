#include <iostream>
#include <cmath>
#include <cassert>

#include <armadillo>

#include "Jacobi_rotation.h"
#include "unittest++/UnitTest++.h"

double sum_offdiag(const arma::mat &A) // Can be used for tolerance check
{
double sum = 0;
  for (int ii=0; ii < A.n_cols; ii++)
  {
     for (int jj=0; jj < A.n_cols; jj++)
        {
          if(ii != jj){
            sum = sum + std::abs(A(ii,jj));
            }
        }
  }
  return sum;
}

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

void max_err_offdiag(const arma::Mat<double> &A, int &k, int &l, double &err)
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

