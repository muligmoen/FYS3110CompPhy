#ifndef lattice_hpp
#define lattice_hpp

#include <iostream>
#include <cstdint>

typedef std::int_fast8_t lat_t; // lattice type


class Lattice
{
private:
  lat_t **lattice;
  const int Lx;
  const int Ly;
  std::string (*print_format)(lat_t);
  
public:
  Lattice(const int lx, const int ly);
  Lattice(const int lx, const int ly, lat_t (*init)(int, int));
  ~Lattice();
  
  void get_size(int &Lx, int&Ly) const;
  
  lat_t &operator() (const int x, const int y);
  lat_t operator() (const int x, const int y) const;
  
  void set_print_format(std::string (*print_func)(lat_t));
  friend std::ostream& operator<< (std::ostream &out, const Lattice &lat);
  
  
  int sum_spins() const;
  int energy() const;
  int energy(const int x, const int y) const;
  int dE(const int x, const int y) const;
};

namespace print_t {
  std::string arrows(lat_t num);
  std::string numbers(lat_t num);
  std::string crazy(lat_t num);
}

namespace init {
  lat_t zeros(int, int);
  lat_t ones(int, int);
  lat_t random(int, int);
}

#endif