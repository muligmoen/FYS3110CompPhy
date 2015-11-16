#include "vector.hpp"
#include "diffusion.hpp"



Vector<double> diffusion::forward_euler(const Vector<double> &init_vec,
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

Vector<double> diffusion::backward_euler(const Vector<double>& init_vec,
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


Vector<double> diffusion::Crank_Nicolson(const Vector<double>& init_vec,
                              const double s, const int steps)
{
  auto vec = init_vec;
  const int N = init_vec.size();
  
  const double a_forward = 1 + s;
  const double b_forward = -s/2;
  
  const double a_backward = 1 - s;
  const double b_backward = s/2;
  
  for (int ii=0; ii<steps; ii++){
    multiply_inplace(vec, a_backward, b_backward);
    solve_inplace(vec, a_forward, b_forward);
    vec[0] = 0;
    vec[N-1] = 0;
  }
  
  return vec;
}
