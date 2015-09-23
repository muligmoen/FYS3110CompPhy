#include <fstream>

#include "filewriter.hpp"



FileWriter::FileWriter(const char* name) : outf(name) 
{
  outf << std::boolalpha;
}

void FileWriter::print(int N, double rho_0, double rho_inf, 
		       bool two_electrons, double omega_r)
{
  outf << N << " : size of system\n";
  outf << rho_0 << " " << rho_inf << " : limits of rho\n";
  outf << two_electrons << " : two electrons\n";
  outf << omega_r << " : omega_r\n" << std::endl;
}
