#ifndef vector_hpp
#define vector_hpp

#include <ostream>

/*!
 *\file vector.hpp 
 * 
 * This file contains the Vector class and some functions for 
 * tridiagonal matrices.
 * 
 */

//! This class holds an array of type T
/*! 
 * This class supports copy, and accessing with operator[].
 * 
 * Some functions are defined as friend functions, and solves sparse multiplication
 * for a special kind of matrice.
 * 
 * There are no checks of any type or size, use vectors with equal size!
 */
template <typename T> 
class Vector;



//! Solves the equation d = Au
/*!
 * A is on the form
 * \f$
 * \begin{pmatrix}
 * a & b & 0 & 0\\
 * b & a & b & 0\\
 * 0 & b & a & b\\
 * \cdots & & 0 & 1
 * \end{pmatrix}\f$
 * and u is the returned vector
 */
template <typename T>
Vector<T> multiply(const Vector<T>& d, const T a, const T b);

//! Same as multiply() but works inplace
template <typename T>
Vector<T>& multiply_inplace(Vector<T> &u, const T a, const T b);

//! Solves d = Au with d and A known
/*!
 * A is a supersparse matrix, with diagonal element a and 
 * first offdiagonal element as b, and the rest as zero.
 */
template <typename T>
Vector<T> solve(const Vector<T>& d, const T a, const T b);

//! Same as solve() but works inplace
template <typename T>
Vector<T>& solve_inplace(Vector<T> &d, const T a, const T b);

//! Prints the vector elements in the format [ v0 v1 v2 ... vN ]
template <typename T>
std::ostream& operator<<(std::ostream& out, const Vector<T> &vector);



template <typename T>
class Vector
{
private:
  //! The size of the vector
  int N;
  
  //! The array elements
  T* vec;
  
public:
  //! Creates the vector
  Vector(const int size) : N(size)
  {
    vec = new T[N];
  }
  
  //! Copy constructor (deep copy)
  Vector(const Vector<T> &other) : N(other.N)
  {
    vec = new T[N];
    for (int ii=0; ii<N; ii++){
      vec[ii] = other[ii];
    }
  }
  
  //! Destructs the vector
  ~Vector()
  {
    delete[] vec;
  }
  
  //! Copies one vector to the other
  Vector<T>& operator=(const Vector<T> &other)
  {
    for (int ii=0;ii<N; ii++){
      vec[ii] = other[ii];
    }
    return *this;
  }
  
  //! Accesing of the elements
  T& operator[](const int N)
  {
    return vec[N];
  }

  //! Constant version of accessor
  T operator[](const int N) const
  {
    return vec[N];
  }
  
  
  //! Addition of two vectors
  Vector<T> operator+(const Vector<T> other)
  {
    Vector<T> new_vec(N);
    
    for (int ii=0; ii<N; ii++){
      new_vec[ii] = vec[ii] + other[ii];
    }
    
    return new_vec;
  }
  
  //! Unary negative operator
  Vector<T> operator-()
  {
    Vector<T> new_vec(N);
    for (int ii=0; ii<N; ii++){
      new_vec[ii] = -vec[ii];
    }
    return new_vec;
  }
  
  //! Gets the size of the vector
  int size() const
  {
    return N;
  }
  
  friend Vector<T> multiply<>(const Vector<T>& d, const T a, const T b);
  
  friend Vector<T>& multiply_inplace<>(Vector<T> &u, const T a, const T b);

  friend Vector<T> solve<>(const Vector<T>& d, const T a, const T b);

  friend Vector<T>& solve_inplace<>(Vector<T> &d, const T a, const T b);
  
  friend std::ostream& operator<< <>(std::ostream& out, const Vector<T> &vector);
  
};

template <typename T>
Vector<T> multiply(const Vector<T>& d, const T a, const T b)
{
  const int N = d.size();
  Vector<T> u(N);
  
  for (int ii=1; ii<N-1; ii++){
    u[ii] = b*d[ii-1] + a*d[ii] + b*d[ii+1];
  }
  
  u[0] = a*d[0] + b*d[1];
  u[N-1] = b*d[N-2] + a*d[N-1];
  
  return u;
}

template <typename T>
Vector<T>& multiply_inplace(Vector<T> &u, const T a, const T b)
{
  const int N = u.size();
  
  T u_previous = u[0];
  u[0] = a*u[0] + b*u[1];
  
  for (int ii=1; ii<N-1; ii++){
    const T u_temp = u[ii];
    u[ii] = b*u_previous + a*u[ii] + b*u[ii+1];
    u_previous = u_temp;
  }
  
  u[N-1] = b*u_previous + a*u[N-1];
  return u;
}


template <typename T>
Vector<T> solve(const Vector<T>& d, const T a, const T b)
{
  const int N = d.size();
  Vector<T> u(N);
   
  Vector<T> beta(N);
  Vector<T> gamma(N);
   
  beta[0] = d[0]/a;
  gamma[0] = b/a;
    
  for (int ii=1; ii<N; ii++){ // forward substition/sweep
    beta[ii] = (d[ii] - b*beta[ii-1])/(a - b*gamma[ii-1]);
    gamma[ii] = b/(a - b*gamma[ii-1]);
  }
    
  u[N-1] = beta[N-1];
    
  for (int ii=N-1; ii>0; ii--){ // backwards substitution
    u[ii-1] = beta[ii-1] - gamma[ii-1]*u[ii];
  }
  return u;
}




template <typename T>
Vector<T>& solve_inplace(Vector<T> &d, const T a, const T b)
{
  const int N = d.size();
    
  Vector<T> beta(N);
  Vector<T> gamma(N);
   
  beta[0] = d[0]/a;
  gamma[0] = b/a;
   
  for (int ii=1; ii<N; ii++){ // forward substition/sweep
    beta[ii] = (d[ii] - b*beta[ii-1])/(a - b*gamma[ii-1]);
    gamma[ii] = b/(a - b*gamma[ii-1]);
  }
    
  d[N-1] = beta[N-1];
  
  for (int ii=N-1; ii>0; ii--){ // backwards substitution
    d[ii-1] = beta[ii-1] - gamma[ii-1]*d[ii];
  }
  return d;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Vector<T> &vector)
{
  out << "[ ";
  for (int ii=0; ii<vector.size(); ii++){
    out << vector[ii] << " ";
  }
  out << "]";
  return out;
}

//! Sets up f(x) in the ranges specified
/*
* x0, x1 are the range of the range, linearly spaced, and f(x) are
* computed for all the points
* 
* This function is inlined to avoid the one-definition rule
*/
inline Vector<double> init_vector(const double x0, const double x1, const int N,
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
#endif