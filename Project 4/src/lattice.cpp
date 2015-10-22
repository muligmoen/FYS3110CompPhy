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
  int x_index = x%Lx;
  if (x_index < 0) x_index += Lx;
  
  int y_index = y%Ly;
  if (y_index <0) y_index += Ly;
  
  return lattice[y_index][x_index];
}

lat_t& Lattice::operator()(const int x, const int y)
{
  int x_index = x%Lx;
  if (x_index < 0) x_index += Lx;
  
  int y_index = y%Ly;
  if (y_index <0) y_index += Ly;
  
  return lattice[y_index][x_index];
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
  stringy << (int)num << " ";
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

std::string print_t::crazy(lat_t num)
{
  std::ostringstream stringy;
  if (num>0) {
    stringy << "  ໒( • ͜ʖ • )७  ";
  } else {
    stringy << " ლ( ◕ 益 ◕ ) ლ ";
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


double Lattice::sum_spins() const
{
  double sum = 0;
  for (int ii=0; ii<Ly; ii++) {
    for (int jj=0; jj<Lx; jj++) {
      sum += this->operator()(jj, ii);
    }
  }
  return sum;
}

double Lattice::energy(const double J) const
{
  double sum_E = 0;
  for (int y = 0; y<Ly; y++) {
    for (int x = 0; x<Lx; x++) {
      double E = this->operator()(x, y);
      E *= this->operator()(x+1,y) + this->operator()(x-1,y) // horisontal neighbours
          + this->operator()(x,y-1) + this->operator()(x,y+2); // vertical neighbours
      sum_E += E;
    }
  }
  return -J*sum_E;
}