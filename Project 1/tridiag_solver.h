#ifndef tridiag_H
#define tridiag_H

#include <armadillo>
#include <fstream>

double f(double x);

double u_theory(double x);

void u_theory(arma::Col<double> &u);

arma::SpMat<double> spsecond_deriv_matr(const int &N);

arma::Mat<double> second_deriv_matr(const int&N);

arma::Col<double> f_column(const double x0, const double x1, const int N);

arma::Col<double> thomas_alg(const double x0, const double x1, const int N);

arma::Col<double> matrix_alg(const double x0, const double x1, const int N, const bool SPARSE=true);

double max_relative_error(const arma::Col<double> &v, const arma::Col<double> &u);

class Writer
{
private:
  std::ofstream outf;
public:
  Writer(const char* name);
  void print(const arma::Col<double> &vec);
  void print(const char *text, const double value);
  void print(const char *text, const int value);
};

#endif
