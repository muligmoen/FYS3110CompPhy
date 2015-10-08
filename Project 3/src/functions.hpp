

#ifndef functions_h
#define functions_h

inline double square_diff(const double x, const double y)
{
  return (x-y)*(x-y); 
}

inline double square_diff(const double x, const double y, const double z)
{
  return x*x + y*y + z*z;
}

void gauss_laguerre(double *x_return, double *w_return,
                    const int n, const double alf);

#endif