
#include "lattice.hpp"
#include "metropolis.hpp"

#include <chrono>
#include <random>
#include <cmath>

auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 generator(seed);

double rand_uniform()
{
  std::uniform_real_distribution<double> distribution(0,1);
  return distribution(generator);
}

void get_indexes(int &x, int &y, const int Lx, const int Ly)
{
  std::uniform_int_distribution<int> Lx_distribution(0, Lx);
  x = Lx_distribution(generator);
  
  std::uniform_int_distribution<int> Ly_distribution(0, Ly);
  y = Ly_distribution(generator);
}

// exp_bet = [exp(-beta*6), exp(-beta*12)]
int try_flip(Lattice& lat, const int x, const int y, const double* exp_beta)
{
  const int energy_diff = lat.dE(x, y);
  
  if (energy_diff <= 0) {
    lat(x,y) *= -1;
    return FlipCodes::FLIPPED;
  } else {
    
    const double comparator = rand_uniform();
    
    double exp_bet;
    if (energy_diff == 6) {
      exp_bet = exp_beta[0];
    } else if (energy_diff == 12) {
      exp_bet = exp_beta[1];
    } else {
      throw "Energy error";
    }
    
    if (comparator < exp_bet) {
      lat(x,y) *= -1;
      return FlipCodes::RAND_FLIPPED;
    } else {
      return FlipCodes::NOT_FLIPPED;
    }
  }
}


void pre_compute_exp(double* exp, const double beta)
{
  exp[0] = std::exp(-beta*6);
  exp[1] = std::exp(-beta*12);
}


