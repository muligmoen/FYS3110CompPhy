#ifndef lattice_hpp
#define lattice_hpp

#include <iostream>
#include <cstdint>

typedef int8_t lat_t; // lattice type


class Lattice
{
private:
  lat_t **lattice;
  const int Lx;
  const int Ly;
public:
  Lattice(const int lx, const int ly);
  Lattice(const int lx, const int ly, lat_t (*init)(int, int));
  ~Lattice();
  
  lat_t &operator() (const int x, const int y);
  lat_t operator() (const int x, const int y) const;
  
  friend std::ostream& operator<< (std::ostream &out, const Lattice &lat);
};


namespace init {
  lat_t zeros(int, int);
  lat_t ones(int, int);
  lat_t random(int, int);
}

#endif