
#include "ising.hpp"
#include <cmath>
#include <iostream>
#include <functional>

Ising::Ising(const int L, const long int seed, const double Jbeta) : L(L), init_seed(seed),
                                            lat(L, L), generator(init_seed)
{ 
  set_beta(Jbeta);
  Magnetisation = lat.sum_spins();
  Energy = lat.energy();
}

Ising::~Ising() { }

void Ising::set_beta(const double Jbeta)
{
  exp_Jbeta[0] = std::exp(-Jbeta*4);
  exp_Jbeta[1] = std::exp(-Jbeta*8);
}

void Ising::set_print_format(std::string (*print_func)(lat_t))
{
  lat.set_print_format(print_func);
}


int Ising::try_flip()
{
  int x, y;
  rand_pos(x,y);
  
  const int energy_diff = lat.dE(x, y);
  int S0 = lat(x,y);
  
  if (energy_diff <= 0) {
    flip(x,y, S0, energy_diff);
    return FlipCodes::FLIPPED;
    
  } else {
    
    const double comparator = rand_uniform();
    const double exp_bet = (energy_diff==4) ? exp_Jbeta[0] : exp_Jbeta[1];
    
    if (comparator < exp_bet) {
      flip(x,y,S0, energy_diff);
      return FlipCodes::RAND_FLIPPED;
    } else {
      return FlipCodes::NOT_FLIPPED;
    }
  }
}

void Ising::try_flip(const int N)
{
  for (int ii = 0; ii<N; ii++) {
    this->try_flip();
  }
}

void Ising::flip(const int x, const int y, const int S, const int dE)
{
  lat(x,y) *= -1;
  const int dS = -2*S;
  Magnetisation += dS;
  Energy += dE;
}

void Ising::flip(const int x, const int y)
{
  const double dE = lat.dE(x,y);
  const int S = lat(x,y);
  flip(x,y, S, dE);
}



double Ising::rand_uniform()
{
  static std::uniform_real_distribution<double> dist(0, 1);
  return dist(generator);
}

void Ising::init_rand()
{
  static std::uniform_int_distribution<lat_t> dist(0, 1);
  auto z = std::bind(dist, generator);
  
  for (int yy=0;yy<L; yy++){
    for (int xx=0; xx<L; xx++){
      lat(xx, yy) = 2*z()-1;
    }
  }
  Energy = lat.energy();
  Magnetisation = lat.sum_spins();
}

void Ising::init_up()
{
  for (int yy=0;yy<L; yy++){
    for (int xx=0; xx<L; xx++){
      lat(xx, yy) = 1;
    }
  }
  Energy = lat.energy();
  Magnetisation = lat.sum_spins();
}



void Ising::rand_pos(int& x, int& y)
{
  static std::uniform_int_distribution<int> L_dist(0, L-1);
  x = L_dist(generator);
  y = L_dist(generator);
}

long int Ising::get_init_seed() const
{
  return init_seed;
}

int Ising::get_energy() const
{
  return Energy;
}

int Ising::get_magnetisation() const
{
  return Magnetisation;
}


std::ostream& operator<< (std::ostream &out, const Ising &ising)
{
  out << ising.lat;
  return out;
}

int Ising::recompute_energy() const
{
  return lat.energy();
}

int Ising::recompute_magnetisation() const
{
  return lat.sum_spins();
}

lat_t* Ising::buffer()
{
  return lat.buffer();
}
