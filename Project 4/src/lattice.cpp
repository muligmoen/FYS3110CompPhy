#include <iostream>
#include <random>
#include <chrono>
#include <functional>

#include "lattice.hpp"


Lattice::Lattice(const int lx, const int ly) : Lx(lx), Ly(ly)
{
  lattice = new lat_t*[Ly];
  for (int ii=0; ii<Ly; ii++){
    lattice[ii] = new lat_t[Lx];
  }
}

Lattice::Lattice(const int lx, const int ly, lat_t (*init)(int, int)) : Lx(lx), Ly(ly)
{
  lattice = new lat_t*[Ly];
  for (int ii=0; ii<Ly; ii++){
    lattice[ii] = new lat_t[Lx];
  }
  
  for (int yy=0; yy<Ly; yy++){
    for (int xx=0; xx<Lx; xx++){
      lattice[yy][xx] = init(xx, yy);
    }
  }
}

Lattice::~Lattice()
{
  for (int ii=0; ii<Ly; ii++){
    delete[] lattice[ii];
  }
  delete[] lattice;
}




lat_t Lattice::operator()(const int x, const int y) const
{
  if (x < 0) {
    return Lattice::operator()(Lx + x, y);
  } else if (y < 0) {
    return Lattice::operator()(x, Ly + y);
  } else {
    return lattice[y%Ly][x%Lx];
  }
}

lat_t& Lattice::operator()(const int x, const int y)
{
  if (x < 0) {
    return Lattice::operator()(Lx + x, y);
  } else if (y < 0) {
    return Lattice::operator()(x, Ly + y);
  } else {
    return lattice[y%Ly][x%Lx];
  }
}



std::ostream& operator<<(std::ostream& out, const Lattice &lat)
{
  for (int ii=0; ii<lat.Ly; ii++){
    for (int jj=0; jj<lat.Lx; jj++){
      if (lat(jj, ii) > 0) {
        out << "8";
      } else {
        out << "-";
      }
    }
    out << "\n";
  }
  return out;
}






lat_t init::zeros(int, int)
{
  return 0;
}

lat_t init::ones(int, int)
{
  return 1;
}

namespace init { namespace rand {
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  static std::mt19937_64 generator(seed);
  const static std::uniform_int_distribution<lat_t> distribution(0,1);
  auto rand = std::bind(distribution, generator);
}}

lat_t init::random(int, int)
{  
  if (init::rand::rand() == 0) {
    return -1;
  } else {
    return 1;
  }
}