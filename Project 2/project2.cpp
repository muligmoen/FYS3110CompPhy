#include <iostream>
#include <cstdlib>

#include <armadillo>
#include "Jacobi_rotation.hpp"
#include "filewriter.hpp"

#include "unittest++/UnitTest++.h"



arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf,
			     bool two_electrons=false, double omega_r=1);

// struct Energies; //defined in filewriter.hpp


Energies min_three_diag(const arma::Mat<double> &A);

arma::vec get_eigv(const arma::Mat<double> &S, int index);

void check_args(int argc, char *argv[], int &N, double &rho_0, double &rho_inf,
		bool &two_electron, double &omega_r)
{
  if (argc>2)
  {
    N = std::atoi(argv[1]);
    rho_inf = std::atof(argv[2]);
    rho_0 = 0;
  } else if (argc>3)
  {
    two_electron = true;
    omega_r = std::atof(argv[3]);
  } else
  {
    std::cerr << "Usage: " << argv[0] << " <N> <rho_inf> <optional: omega_r> \n\n";
    std::cerr << "<N> an integer which gives the size and approximation\n";
    std::cerr << "<rho_inf> max limit of rho to simulate infinity\n";
    std::cerr << "<omega_r> given if it is a two_electron problem" << std::endl;
    std::exit(1);
  }
}

int main(int argc, char *argv[])
{
  UnitTest::RunAllTests();
  
  int N;
  double rho_0, rho_inf, omega_r;
  bool two_electrons;
  
  check_args(argc, argv, N, rho_0, rho_inf, two_electrons, omega_r);
  
  
  arma::Mat<double> A = ham_matrix(N, rho_0, rho_inf);
  
  arma::Mat<double> S(N, N, arma::fill::eye);
  
  rotate_to_diag_with_eigvec(A, S, 1e-10);

  auto E = min_three_diag(A);
  
  
  FileWriter file("test.txt");
  file.print(N, rho_0, rho_inf);
  file.print(E);
  
  for (int iii=0; iii<3; iii++)
  {
    auto Evec = get_eigv(S, E.indexes[iii]);
    file.print(Evec);
  }
  
  return 0;
}



arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf,
			     bool two_electrons, double omega_r)
{
  arma::Mat<double> A(N, N, arma::fill::zeros);
  
  double h = (rho_inf-rho_0)/(double)N;
  double inv_h_square = 1/(h*h);
  A.diag(1) -= inv_h_square;
  A.diag(-1) -= inv_h_square;
  
  
  if (two_electrons)
  {
    double omega_square = omega_r*omega_r;
    for (int iii=0; iii<N; iii++)
    {
      double rho = rho_0 + iii*h;
      A(iii,iii) = 2*inv_h_square + omega_square*rho*rho + 1.0/rho;
    }
  } else
  {
    for (int iii=0; iii<N; iii++)
    {
      double rho = rho_0 + iii*h;
      A(iii,iii) = 2*inv_h_square + rho*rho;
    }
  }
  
  return A;
}



Energies min_three_diag(const arma::Mat<double> &A)
{
  Energies E;
  
  double *min = E.Energy;
  int *ind = E.indexes;
  
  for (int iii=0; iii<3; iii++)
  {
    min[iii] = arma::datum::inf;
    ind[iii] = -1;
  }
  
  for (int iii=0; iii<(int)A.n_rows; iii++)
  {
    double num = A(iii,iii);
    if (num<min[2])
    {
      min[2] = num;
      ind[2] = iii;
      if (num < min[1])
      {
	min[2] = min[1];
	ind[2] = ind[1];
	min[1] = num;
	ind[1] = iii;
	if (num < min[0])
	{  
	  min[1] = min[0];
	  ind[1] = ind[0];
	  min[0] = num;
	  ind[0] = iii;
	}
      }
    }
  }
  
  return E;
}


arma::vec get_eigv(const arma::Mat<double> &S, int index)
{
  int N = S.n_cols;
  arma::vec x(N);
  
  double square_sum = 0;
  
  for (int jjj=0; jjj<N; jjj++)
  {
    double S_ij = S(index, jjj);
    x(jjj) = S_ij;
    square_sum += S_ij*S_ij;
  }
  x /= square_sum;
  
  return x;
}
