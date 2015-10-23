#ifndef metropolis_hpp
#define metropolis_hpp

#include "lattice.hpp"


enum FlipCodes: int
{
  NOT_FLIPPED,
  FLIPPED,
  RAND_FLIPPED,
};


void get_indexes(int &x, int &y, const int Lx, const int Ly);

double rand_uniform();

int try_flip(Lattice& lat, const int x, const int y, const double* exp_beta);

void pre_compute_exp(double* exp, const double beta);
// Changes exp to: exp = [exp(-beta*6), exp(-beta*12)]


#endif