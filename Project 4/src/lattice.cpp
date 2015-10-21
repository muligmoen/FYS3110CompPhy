#include <iostream>
#include <sstream>
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
  print_format = print_t::numbers;
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
  print_format = print_t::numbers;
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
      out << lat.print_format(lat(jj,ii));
    }
    out << "\n";
  }
  return out;
}


void Lattice::set_print_format(std::string (*print_func)(lat_t))
{
  print_format = print_func;
}

std::string print_t::numbers(lat_t num)
{
  std::ostringstream stringy;
  stringy << (int)num;
  return stringy.str();
}

std::string print_t::arrows(lat_t num)
{
  std::ostringstream stringy;
  if (num > 0) {
    stringy << "▲";
  } else if (num < 0) {
    stringy << "▼";
  } else {
    stringy << "-";
  }
  return stringy.str();
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
  std::mt19937_64 generator(seed);
  const std::uniform_int_distribution<lat_t> distribution(0,1);
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