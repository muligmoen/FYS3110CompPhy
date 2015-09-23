#ifndef filewriter_h
#define filewriter_h

#include <armadillo>
#include <fstream>


class FileWriter
{
private:
  std::ofstream outf;
public:
  FileWriter(const char* name);
  void print(int N, double rho_0, double rho_inf, 
	     bool two_electrons=false, double omega_r=1.0);
};









#endif
