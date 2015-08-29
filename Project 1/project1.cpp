#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <cstdlib>



double f(double x);

double u_theory(double x);

void fill_matrix(arma::Mat<double> &A);

void fill_matrix(arma::SpMat<double> &A);

void f_column(arma::Col<double> &b);

void u_theory(arma::Col<double> &u);


class Writer
{
private:
  std::ofstream outf;
public:
  Writer(const char* name, const int matrix_size);
  ~Writer();
  void print(arma::Col<double> &vec);
};

int main( int argc, char *argv[] )
{
  if (argc < 3 || atoi(argv[1]) < 2)
  {
    std::cout << "Usage: " << argv[0] << " N filename" << std::endl;
    std::cout << "With N > 1" << std::endl;
    exit(0);
  }
  
  const int matrix_size = atoi(argv[1]);
  
  //Use the top constructor if armadillo version < 5
  //arma::Mat<double> A(matrix_size, matrix_size); 
  arma::SpMat<double> A(matrix_size, matrix_size);
  
  
  fill_matrix(A);
  
  arma::Col<double> b(matrix_size);
  f_column(b);
  
  
  //Top for armadillo < 5
  //arma::Col<double> v = solve(A, b);
  arma::Col<double> v = spsolve(A, b);
   
  
  
  arma::Col<double> u(matrix_size);
  
  u_theory(u);
  
  //std::cout << "Standard deviation of (u_{theory}-u_{computed}) = " << stddev(u-v) << std::endl;
  

  Writer fileprinter(argv[2], matrix_size);
  fileprinter.print(b);
  
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

void f_column(arma::Col<double> &b)
{
  const int N = b.n_elem;
  double h = 1.0/(double)N;
  for (int iii=0; iii<N; iii++)
  {
    b[iii] = f(iii*h)*h*h;
  }
}

void fill_matrix(arma::Mat<double> &A)
{
  A.diag(0) += 2.0;
  A.diag(-1) += -1.0;
  A.diag(1) += -1.0;
}


void fill_matrix(arma::SpMat<double> &A)
{
  const int sizeA = A.n_rows;
  for (int diag=0; diag<sizeA; diag++)
  {
    A(diag, diag) = 2;
  }
  for (int diag=0; diag < sizeA; diag++)
  {
    if (diag < sizeA - 1)
    {
      A(diag, diag+1) = -1;
    }
  }
  for (int diag=0; diag < sizeA; diag++)
  {
    if (diag > 0)
    {
      A(diag, diag-1) = -1;
    }
  }
}

Writer::Writer(const char *name, const int array_length) : outf(name)
{
  outf << array_length << "\n\n";
};

Writer::~Writer() { };

void Writer::print(arma::Col<double> &vec)
{
  outf.precision(10);
  vec.raw_print(outf);
  outf << "\n";
}

