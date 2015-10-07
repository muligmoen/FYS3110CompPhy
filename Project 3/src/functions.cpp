#include <cmath>

inline double hydrogen_1s(const double r, const double alpha)
{
  return std::exp(-alpha*r);
}

inline double square_diff(const double x, const double y)
{
  return (x-y)*(x-y); 
}

inline double inverse_diff_distance(const double x1, const double y1, const double z1,
                                    const double x2, const double y2, const double z2)
{
  const double sum_squared_diff = square_diff(x1, x2) +
                                    square_diff(y1, y2) + square_diff(z1,z2);
  return 1/std::sqrt(sum_squared_diff);
}