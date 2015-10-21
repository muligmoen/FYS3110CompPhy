#ifndef lattice_hpp
#define lattice_hpp

#include <iostream>
#include <cstdint>

typedef int8_t lat_t; // lattice type
typedef uint8_t lat_it; // lattice index type


class Lattice
{
private:
  lat_t **lattice;
  const lat_it Lx;
  const lat_it Ly;
public:
  Lattice(const lat_it lx, const lat_it ly);
  Lattice(const lat_it lx, const lat_it ly, lat_t (*init)(lat_it, lat_it));
  ~Lattice();
  
  lat_t &operator() (const int x, const int y);
  lat_t operator() (const int x, const int y) const;
  
  lat_t &operator() (const lat_it x, const lat_it y);
  
  friend std::ostream& operator<< (std::ostream &out, const Lattice &lat);
  
};

namespace init {
  lat_t zeros(lat_it, lat_it);
  
  lat_t ones(lat_it, lat_it);
  
  lat_t random(lat_it, lat_it);
}

#endif