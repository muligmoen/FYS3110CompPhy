#include "integrate.hpp"



Vector<double> forward_euler(const Vector<double> &init_vec,
                             const double s, const int steps)
{
  auto new_vec = init_vec;
  const int N = init_vec.size();
  
  const double a = 1-2*s;
  const double b = s;
  for (int ii=0; ii<steps; ii++){
    multiply_inplace(new_vec, a, b);
    new_vec[0] = 0;
    new_vec[N-1] = 0;
  }
  return new_vec;
}

Vector<double> backward_euler(const Vector<double>& init_vec,
                                const double s, const int steps)
{
  auto vec = init_vec;
  const int N = init_vec.size();
  
  const double a = 1 + 2*s;
  const double b = -s;
  for (int ii=0; ii<steps; ii++){
    solve_inplace(vec, a, b);
    vec[0] = 0;
    vec[N-1] = 0;
  }
  
  return vec;
}


Vector<double> Crank_Nicolson(const Vector<double>& init_vec,
                              const double s, const int steps)
{
  auto vec = init_vec;
  const int N = init_vec.size();
  
  const double a1 = 1 + 2*s;
  const double b1 = -s;
  
  const double a2 = 1 - 2*s;
  const double b2 = s;
  
  for (int ii=0; ii<steps; ii++){
    multiply_inplace(vec, a2, b2);
    solve_inplace(vec, a1, b1);
    vec[0] = 0;
    vec[N-1] = 0;
  }
  
  return vec;
}
