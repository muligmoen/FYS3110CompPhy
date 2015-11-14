#include <iostream>
#include "vector.hpp"

Vector<double> forward_euler(const Vector<double> &init_vec,
                             const double s, const int steps);

Vector<double> backward_euler(const Vector<double> &init_vec,
                              const double s, const int steps);

Vector<double> init_vector(const double x0, const double x1, const int N,
                 double (*f)(double));



int main()
{
  const int N = 25;
  //Vector<double> vec(N);
  
  auto f_steady = [](double x){return 1 - x;};
  
  auto u_steady = init_vector(0, 1, N, f_steady);
  
  auto init_vec = -u_steady;
  
  init_vec[0] = 0;
  init_vec[N-1] = 0;
  
  const double s = 0.2;
  auto vec = forward_euler(init_vec, s, 100);
  
  std::cout << vec + u_steady << std::endl;
}


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

Vector<double> init_vector(const double x0, const double x1, const int N,
                           double (*f)(double))
{
  Vector<double> vec(N);
  const double dx = (x1 - x0)/(N-1);
  
  for (int ii = 0; ii<N; ii++){
    const double x = x0 + dx*ii;
    vec[ii] = f(x);
  }
  return vec;
}
