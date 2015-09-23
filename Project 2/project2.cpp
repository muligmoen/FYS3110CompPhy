#include <iostream>
#include <cstdlib>

#include <armadillo>
#include "Jacobi_rotation.hpp"
#include "filewriter.hpp"

#include "unittest++/UnitTest++.h"



arma::Mat<double> ham_matrix(int N, double rho_0, double rho_inf,
			     bool two_electrons=false, double omega_r=1);

struct Energies
{
  double Energy[3];
  int indexes[3];
};

//int main(int argc, char *argv[])
int main()
{
  UnitTest::RunAllTests();
  /*
  int N = -1;
  int rho_inf;
  
  if (argc < 3)
  {
    N = 100;
    rho_inf = 10;
    
    std::cout << "Running standard settings" << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "rho_inf = " << rho_inf << std::endl;
  }
  else 
  {
    N = std::atoi(argv[1]);
    rho_inf = std::atof(argv[2]);
    if (N<1 || rho_inf<=0)
    {
      std::cerr << "Invalid arguments, exiting" << std::endl;
      std::exit(1);
    }
  }*/
/*
  double rho_0 = 0.0;
  arma::Mat<double> A = ham_matrix(N, rho_0, rho_inf);
  
  arma::Mat<double> S(N, N, arma::fill::eye);
  
  rotate_to_diag_with_eigvec(A, S, 1e-10);
*/

  FileWriter file("test.txt");
  file.print(10, 0, 10, true, 20);
  
  
  Energies E;
  E.Energy[0] = 3;
  
  std::cout << E.Energy[0] << std::endl;
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

