#ifndef ising_hpp
#define ising_hpp

#include <iostream>
#include <cmath>

#include "lattice.hpp"

#include <random>

enum FlipCodes: int
{
  NOT_FLIPPED,
  FLIPPED,
  RAND_FLIPPED,
};

class Ising
{
private:
  const int L;
  const long int init_seed;
  const double exp_Jbeta[2];
  Lattice lat;
  std::mt19937 generator;
  
  int Energy;
  int Magnetisation;
  
public:
  Ising(const int L, const long int seed, const double Jbeta);
  ~Ising();
  
  void set_print_format(std::string (*print_func)(lat_t));
  
  long int get_init_seed() const;
  
  int get_energy() const;
  int recompute_energy() const;
  
  int get_magnetisation() const;
  int recompute_magnetisation() const;
  
  int rand_init(); // gives -1, 1
  void rand_pos(int &x, int &y); // gives 0,1,2,...,L-1
  double rand_uniform(); // gives 0,...,1
  
  void flip(const int x, const int y, const int dS, const int dE);
  
  int try_flip();
  
  void try_flip(const int N);
  
  lat_t* buffer();
  
  friend std::ostream& operator<< (std::ostream &out, const Ising &ising);
  
};












#endif