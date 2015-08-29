#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <cstdlib>


double f(double x);

double u_theory(double x);

void u_theory(arma::Col<double> &u);

arma::SpMat<double> second_deriv_matr(const int &N);

arma::Mat<double> second_deriv_matr(const int&N, bool sparse);

arma::Col<double> f_column(const double x0, const double x1, const int N);

class Writer
{
private:
  std::ofstream outf;
public:
  Writer(const char* name, const int matrix_size);
  void print(arma::Col<double> &vec);
};


int main( int argc, char *argv[] )
{
  if (argc < 3 || atoi(argv[1]) < 2)
  {
    std::cout << "Usage: " << argv[0] << " N filename" << std::endl;
    std::cout << "With N > 1" << std::endl;
    exit(1);
  }
  
  const int matrix_size = atoi(argv[1]);
  
  
  arma::Col<double> f = f_column(0, 1, matrix_size);
  
  
  // non-sparse methods (arma version < 5)
  //arma::Mat<double> A = second_deriv_matr(matrix_size, false); 
  //arma::Col<double> v = solve(A, f);
  
  arma::SpMat<double> A = second_deriv_matr(matrix_size);
  arma::Col<double> v = spsolve(A, f);
  
  
  arma::Col<double> u(matrix_size);
  
  u_theory(u);
  
  //std::cout << "Standard deviation of (u_{theory}-u_{computed}) = " << stddev(u-v) << std::endl;
  

  Writer fileprinter(argv[2], matrix_size);
  fileprinter.print(v);
  
  fileprinter.print(u);
  
  return 0;
}


double f(double x)
{
  return 100.0*std::exp(-10.0*x);
}

double u_theory(double x)
{
  return 1.0 - (1.0 - std::exp(-10.0))*x - std::exp(-10.0*x);
}

void u_theory(arma::Col<double> &u)
{
  const int N = u.n_elem;
  double h = 1.0/(double)N;
  for (int iii=0; iii<N; iii++)
  {
    u[iii] = u_theory(iii*h);
  }
}

arma::Col<double> f_column(const double x0, const double x1, const int N)
{
  arma::Col<double> b(N);
  
  double h = (x1-x0)/(double)N;
  
  for (int iii=0; iii<N; iii++)
  {
    double x = x0 + h*iii;
    b[iii] = f(x)*h*h;
  }
  return b;
}

arma::Mat<double> second_deriv_matr(const int&N, bool sparse)
{
  if (sparse)
  {
    std::cout << "Please use second_deriv_matr(N) instead" << std::endl;
  }
  arma::Mat<double> A(N, N);
  A.diag(0) += 2.0;
  A.diag(1) += -1;
  A.diag(-1) += -1;
  return A;
}

arma::SpMat<double> second_deriv_matr(const int& N)
{
  arma::SpMat<double> A(N, N);
  for (int diag=0; diag < N; diag++)
  {
    A(diag, diag) = 2; //diagonal
    if (diag < N - 1) //upper band
    {
      A(diag, diag+1) = -1;
    }
    if (diag > 0) //lower band
    {
      A(diag, diag-1) = -1;
    }
  }
  return A;
}


Writer::Writer(const char* name, const int matrix_size) : outf(name)
{
  outf << matrix_size << "\n\n";
};


void Writer::print(arma::Col<double> &vec)
{
  outf.precision(10);
  vec.raw_print(outf);
  outf << "\n";
}

