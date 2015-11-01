#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <chrono>

#include <omp.h>

#include "lattice.hpp"
#include "ising.hpp"

//! Makes uniformly distributed points from T0 to T1 inclusive
void linspace(const double T0, const double T1, const int N, double* T);

void print_to_file(const double *T, const double *E, const double* Cv, const double* M, const double* chi, const double* ar,
                   const int N, const char* filename = "test.txt");



int main(int argc, char **argv)
{
  //Problem e) stud of critical temperature
  if (argc>2 && !std::strcmp(argv[1], "e")){
    const int dim = std::atoi(argv[2]);
    
    const int Ntemps = 15;
    const int N = 100000000;
  

  
    double *T, *E, *Cv, *M, *chi, *ar;
    T = new double[Ntemps];
    E = new double[Ntemps];
    Cv = new double[Ntemps];
    M = new double[Ntemps];
    chi = new double[Ntemps];
    ar = new double[Ntemps];
    
    linspace(2.0, 2.4, Ntemps, T);

  
    const auto global_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

    #pragma omp parallel for default(shared)
    for (int ii = 0; ii<Ntemps; ii++){
      const auto seed = global_seed + ii;
      find_statistics(N, dim, T[ii], E[ii], Cv[ii], M[ii], chi[ii], seed, ar[ii]);
    }
    
    print_to_file(T, E, Cv, M, chi, ar, Ntemps);
    
    delete[] T;
    delete[] E;
    delete[] Cv;
    delete[] M;
    delete[] chi;
    delete[] ar;
  }
  
}





void linspace(const double T0, const double T1, const int N, double* T)
{
  const double deltaT = (T1 - T0)/(N-1);
  for (int ii=0; ii<N; ii++){
    T[ii] = T0 + deltaT*ii;
  }
}

void print_to_file(const double *T, const double* E, const double* Cv, const double* M, const double* chi, const double* ar,
                   const int N, const char* filename)
{
  std::ofstream outf(filename);
  outf << N << " N" << "\n\n\n\n\n\n\n\n\n";
  outf << "T E M Cv chi acceptance_rate\n";
  for (int ii=0; ii<N; ii++){
    outf << T[ii] << " " << E[ii] << " " << M[ii] << " " << Cv[ii] << " " << chi[ii] << " " << ar[ii] << "\n";
  }
}

