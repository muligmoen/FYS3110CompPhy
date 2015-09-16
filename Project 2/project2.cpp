#include <iostream>
#include <armadillo>

template <typename T>
inline T abs(T num)
{
  return num ? num > 0 : -num;
}

void largest_value(const arma::mat &A, int &k, int &l)
{
double val = 0;
  for (int ii=0; ii < A.n_cols; ii++)
  {
     for (int jj=0; jj < A.n_cols; jj++)
        {
          if((ii != jj) && (std::abs(A(ii,jj))) > val) {
            val = std::abs(A(ii,jj));
            k=ii;
            l=jj;
            }
         
        }
  }
}

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

int main()
{
  arma::Mat<double> A;
  A << 1 << 5 << arma::endr
    << 0 << 3;
   
  std::cout << A << std::endl;
}

