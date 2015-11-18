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


/*! 
 * 
 * @brief This class holds a vector of type T
 * 
 * This class supports copy, and accessing with operator[].
 * 
 * Some functions are defined as friend functions, and solves sparse multiplication
 * for a special kind of matrice.
 * 
 * There are no checks of any type or size, use vectors with equal size!
 * 
 * This should be used as a replacement for pointers to arrays, as 
 * it is destructed automatically and has some types overloaded
 */
template <typename T> 
class Vector;



 
/*!
 * 
 * @brief Solves the equation d = Au
 * 
 * @param d The input vector
 * @param a Diagonal element of A
 * @param b Non-diagonal element of A
 * 
 * A is on the form
 * \f$
 * \begin{pmatrix}
 * a & b & 0 & \cdots\\
 * b & a & b & 0\\
 * 0 & b & a & b\\
 *  \cdots  &   &  \ddots\\
 * & & b & a
 * \end{pmatrix}\f$
 */
template <typename T>
Vector<T> multiply(const Vector<T>& d, const T a, const T b);

/*!
 * @brief Same as multiply() but works inplace
 * 
 * @see multiply()
 */
template <typename T>
Vector<T>& multiply_inplace(Vector<T> &u, const T a, const T b);


/*!
 * @brief Solves d = Au with d and A known
 * 
 * @param d The input vector
 * @param a The diagonal element of A
 * @param b The off-diagonal element of A
 * 
 * @see multiply() For a description of A
 */
template <typename T>
Vector<T> solve(const Vector<T>& d, const T a, const T b);

/*!
 * @brief Solves d = Au
 * 
 * Inplace version of solve()
 * 
 * @see solve()
 */
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
  
  //! The vector elements
  T* vec;
  
public:
  /*! 
   * @brief Creates the vector
   * 
   * @param size Size of vector
   * 
   * @return A vector with uninitialized values
   */
  Vector(const int size) : N(size)
  {
    vec = new T[N];
  }
  
  /*! 
   * 
   * @brief Function constructor
   * 
   * @param size Size of vector
   * @param f Function which decides the elements in the vector
   * 
   * @return A constructed Vector with the elements given by f
   * 
   * By calling a function this vector can be initialised with 
   * whatever the caller wants. This is mostly used for zeroing
   * the vector. If all elements are to be set to zero, one 
   * way is by calling with a lambda function,
   * 
   * Vector< int > vector(5, [](){return 0;})
   * 
   */
  Vector(const int size, T (*f)() ) : N(size)
  {
    vec = new T[N];
    for (int ii=0; ii<N; ii++){
      vec[ii] = f();
    }
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
  Vector<T> operator+(const Vector<T> other) const
  {
    Vector<T> new_vec(N);
    
    for (int ii=0; ii<N; ii++){
      new_vec[ii] = vec[ii] + other[ii];
    }
    
    return new_vec;
  }
  
  //! Unary negative operator
  Vector<T> operator-() const
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