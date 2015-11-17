#include "vector.hpp"
#include "diffusion.hpp"



Vector<double> diffusion::forward_euler(const Vector<double> &init_vec,
                             const double alpha, const int steps)
{
  auto new_vec = init_vec;
  
  const double a = 1-2*alpha;
  const double b = alpha;
  for (int ii=0; ii<steps; ii++){
    multiply_inplace(new_vec, a, b);
  }
  return new_vec;
}

Vector<double> diffusion::backward_euler(const Vector<double>& init_vec,
                                const double alpha, const int steps)
{
  auto vec = init_vec;
  
  const double a = 1 + 2*alpha;
  const double b = -alpha;
  for (int ii=0; ii<steps; ii++){
    solve_inplace(vec, a, b);
  }
  
  return vec;
}


Vector<double> diffusion::Crank_Nicolson(const Vector<double>& init_vec,
                              const double alpha, const int steps)
{
  auto vec = init_vec;
  
  const double a_forward = 1 + alpha;
  const double b_forward = -alpha/2;
  
  const double a_backward = 1 - alpha;
  const double b_backward = alpha/2;
  
  for (int ii=0; ii<steps; ii++){
    multiply_inplace(vec, a_backward, b_backward);
    solve_inplace(vec, a_forward, b_forward);
  }
  
  return vec;
}
